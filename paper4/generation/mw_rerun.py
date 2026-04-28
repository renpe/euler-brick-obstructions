"""
Mordell-Weil rerun for under-covered high-rank fibres (pub schema).

Same architecture as mw_dispatcher.py, but with a longer subprocess
timeout, higher descent_second_limit and a higher MAX_SCALAR. Used after
the systematic rank-scan to harvest hits from newly-found
analytic_rank ≥ 2 fibres.

Usage:
    python3 mw_rerun.py                                 # default RERUN_FIBERS
    python3 mw_rerun.py '[[41,32],[97,40],[125,36]]'    # custom list
    python3 mw_rerun.py FIBERS_JSON TIMEOUT DESCENT MAX_SCALAR
"""
import os
import sys
import json
import subprocess
from multiprocessing import Pool
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db

SAGE = pub_db.SAGE_BIN
WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mw_fibre_worker.sage")

# Defaults (can be overridden via argv)
TIMEOUT_SEC = 600
DESCENT_LIMIT = 20
MAX_SCALAR = 5

# Default fibres (under-covered high-rank, from §3 of the paper):
RERUN_FIBERS = [
    (77, 48),
    (592, 539),
    (88, 7),
]


def run_one(args):
    m, n = args
    t0 = time()
    try:
        proc = subprocess.run(
            [SAGE, WORKER, str(m), str(n), str(MAX_SCALAR), str(DESCENT_LIMIT)],
            capture_output=True, text=True, timeout=TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired:
        return (m, n, "timeout", 0, [], time() - t0)
    if proc.returncode != 0:
        return (m, n, f"exit {proc.returncode}", 0, [], time() - t0)
    candidates = []
    status = "fail"
    gens_count = 0
    for line in proc.stdout.strip().split("\n"):
        line = line.strip()
        if not line.startswith("{"): continue
        try:
            j = json.loads(line)
            status = j.get("status", "fail")
            candidates = j.get("candidates", [])
            gens_count = j.get("gens_count", 0)
        except Exception:
            continue
    return (m, n, status, gens_count, candidates, time() - t0)


def main():
    global TIMEOUT_SEC, DESCENT_LIMIT, MAX_SCALAR
    if len(sys.argv) > 1:
        fibres = [tuple(p) for p in json.loads(sys.argv[1])]
    else:
        fibres = RERUN_FIBERS
    if len(sys.argv) > 2: TIMEOUT_SEC = int(sys.argv[2])
    if len(sys.argv) > 3: DESCENT_LIMIT = int(sys.argv[3])
    if len(sys.argv) > 4: MAX_SCALAR = int(sys.argv[4])

    print(f"Rerun {len(fibres)} fibres, timeout={TIMEOUT_SEC}s, "
          f"descent_limit={DESCENT_LIMIT}, MAX_SCALAR={MAX_SCALAR}", flush=True)
    print(f"Fibres: {fibres}", flush=True)

    conn = pub_db.connect(autocommit=False)
    cur = conn.cursor()
    cur.execute("SELECT a, b, m, n FROM pub.master_hits")
    seen = set((int(a), int(b), int(m), int(n)) for a, b, m, n in cur.fetchall())
    print(f"  {len(seen)} existing hits loaded\n", flush=True)

    n_inserted_total = 0
    t0 = time()

    workers = min(20, len(fibres))
    with Pool(processes=workers) as pool:
        for m, n, status, gens_count, candidates, dt in \
                pool.imap_unordered(run_one, fibres, chunksize=1):
            inserted = 0
            if status == "ok":
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
            print(f"  ({m},{n}) {status} gens={gens_count} cand={len(candidates)} "
                  f"+{inserted} (gesamt: {n_inserted_total}) in {dt/60:.1f}min",
                  flush=True)

    print(f"\nFertig in {(time()-t0)/60:.1f} min: {n_inserted_total} new bricks")
    conn.close()


if __name__ == "__main__":
    main()
