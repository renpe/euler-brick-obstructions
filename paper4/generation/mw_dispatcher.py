"""
Mordell-Weil dispatcher for the pub schema.

For each fibre $(m, n)$ in pub.master_hits with hit count in
[min_hits, max_hits], spawn a Sage subprocess (mw_fibre_worker.sage)
with a hard timeout. Insert resulting bricks via pub_db.insert_master_hit().

Usage:
    python3 mw_dispatcher.py [workers=4] [min_hits=4] [max_hits=7]

Skip rules:
  - Skip fibres already covered by Mordell-Weil provenance.
  - Skip fibres tagged with classical families that we don't want to
    duplicate (Saunderson / Lenhart / Himane-T1/T2/T3).
"""
import os
import sys
import json
import subprocess
from multiprocessing import Pool
from time import time

# Common helpers
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db

SAGE = pub_db.SAGE_BIN
WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mw_fibre_worker.sage")
TIMEOUT_SEC = 60
MAX_SCALAR = 3
DESCENT_LIMIT = 15

# Skip families: any fibre that already has a hit tagged with one of these
# is excluded. This matches §5.1 of the paper, where MW generation
# focuses on fibres outside the classical parametrisations.
SKIP_FAMILIES = {"Saunderson", "Lenhart", "Himane-T1", "Himane-T2", "Himane-T3"}


def run_one(args):
    """Worker: spawn Sage subprocess for one fibre."""
    m, n = args
    t0 = time()
    try:
        proc = subprocess.run(
            [SAGE, WORKER, str(m), str(n), str(MAX_SCALAR), str(DESCENT_LIMIT)],
            capture_output=True, text=True, timeout=TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired:
        return (m, n, "timeout", [], time() - t0)
    if proc.returncode != 0:
        return (m, n, f"exit {proc.returncode}", [], time() - t0)
    candidates = []
    status = "fail"
    for line in proc.stdout.strip().split("\n"):
        line = line.strip()
        if not line.startswith("{"): continue
        try:
            j = json.loads(line)
            status = j.get("status", "fail")
            candidates = j.get("candidates", [])
        except Exception:
            continue
    return (m, n, status, candidates, time() - t0)


def select_target_fibres(cur, min_hits: int, max_hits: int) -> list[tuple[int, int]]:
    """Return list of (m, n) fibres meeting the hit-count and skip criteria."""
    cur.execute("""
        WITH fiber_counts AS (
            SELECT m, n, count(*) AS hits
            FROM pub.master_hits
            GROUP BY m, n
            HAVING count(*) >= %s AND count(*) <= %s
        ),
        mw_done AS (
            SELECT DISTINCT mh.m, mh.n
            FROM pub.master_hits mh
            JOIN pub.hit_groups hg ON hg.hit_id = mh.id
            JOIN pub.generator_groups g ON g.id = hg.group_id
            WHERE g.short_name LIKE 'MW-%%' AND hg.is_primary
        ),
        family_fibres AS (
            SELECT DISTINCT mh.m, mh.n
            FROM pub.master_hits mh
            JOIN pub.hit_groups hg ON hg.hit_id = mh.id
            JOIN pub.generator_groups g ON g.id = hg.group_id
            WHERE g.short_name = ANY(%s)
        )
        SELECT m, n FROM fiber_counts
        WHERE (m, n) NOT IN (SELECT m, n FROM mw_done)
          AND (m, n) NOT IN (SELECT m, n FROM family_fibres)
        ORDER BY m + n
    """, (min_hits, max_hits, list(SKIP_FAMILIES)))
    return [(int(m), int(n)) for m, n in cur.fetchall()]


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    min_hits = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    max_hits = int(sys.argv[3]) if len(sys.argv) > 3 else 7

    print(f"workers={workers}, hits ∈ [{min_hits}, {max_hits}]", flush=True)

    conn = pub_db.connect(autocommit=False)
    cur = conn.cursor()

    fibres = select_target_fibres(cur, min_hits, max_hits)
    print(f"  {len(fibres)} target fibres", flush=True)
    if not fibres:
        return

    # Existing hits for dedup (so we don't try to re-insert)
    cur.execute("SELECT a, b, m, n FROM pub.master_hits")
    seen = set((int(a), int(b), int(m), int(n)) for a, b, m, n in cur.fetchall())

    n_done = 0
    n_inserted_total = 0
    n_timeout = 0
    n_nogens = 0
    t0 = time()

    with Pool(processes=workers) as pool:
        for m, n, status, candidates, dt in pool.imap_unordered(run_one, fibres, chunksize=1):
            n_done += 1
            inserted = 0
            if status == "timeout":
                n_timeout += 1
            elif status == "nogens":
                n_nogens += 1
            elif status == "ok":
                # Ensure the MW-{m}-{n} provenance group exists
                mw_gid = pub_db.ensure_mw_group(cur, m, n)
                conn.commit()
                for ab in candidates:
                    a, b = int(ab[0]), int(ab[1])
                    if (a, b, m, n) in seen: continue
                    try:
                        hit_id = pub_db.insert_master_hit(cur, a, b, m, n, mw_gid)
                        if hit_id is not None:
                            inserted += 1
                            seen.add((a, b, m, n))
                    except Exception:
                        conn.rollback()
                        cur = conn.cursor()
                conn.commit()
                n_inserted_total += inserted
            print(f"  [{n_done}/{len(fibres)}] ({m},{n}) {status} "
                  f"in {dt:.1f}s, +{inserted} (gesamt: {n_inserted_total}, "
                  f"timeout={n_timeout}, nogens={n_nogens})", flush=True)

    print(f"\nFertig in {(time()-t0)/60:.1f} min: {n_inserted_total} neu, "
          f"timeout={n_timeout}, nogens={n_nogens}")
    conn.close()


if __name__ == "__main__":
    main()
