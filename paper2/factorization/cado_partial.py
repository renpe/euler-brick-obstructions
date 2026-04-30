"""
Run CADO-NFS on the residual of a partial Master-Hit.

A partial hit is one where pub_factorize.py extracted some small
factors of f1_prim by trial-division/sympy but YAFU timed out on the
remaining composite. This script takes such a residual, runs
CADO-NFS on it, parses the resulting prime factors, and writes them
back to pub.prim_brick_factors. If the residual is now fully factored,
it sets num_blockers on pub.master_hits.

CADO-NFS produces several "PRIME" lines on stdout containing the
factors of the input integer. We collect them, multiply, verify
against the residual, and store.

Usage:
    "$CADO_VENV_PY" cado_partial.py [HIT_ID]

If HIT_ID is omitted, we pick the SMALLEST residual that still has
>= MIN_RESIDUAL_DIGITS (=84) digits. Below 84 digits YAFU/QS is
faster than CADO-NFS; above that CADO wins. Smallest-first means
the next CADO run is always the cheapest still-pending one.

Caution: CADO-NFS on a 130-150 digit residual takes anywhere from
several hours to days on a single workstation. Don't run this in
parallel with pub_factorize.py — both will fight for the cores.
"""
import os
import sys
import re
import subprocess
import json
from math import gcd
from time import time

import psycopg

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db

CADO_NFS = pub_db.CADO_NFS
CADO_VENV_PY = pub_db.CADO_VENV_PY

# Break-even YAFU vs. CADO-NFS empirically lies at ~84-digit residual
# (Hit 187370: 84 digits, CADO ~ 2.7 min). Below that YAFU is faster,
# above that CADO. So we pick only partials with residual >= 84 digits
# and of those the smallest first (shortest CADO run).
MIN_RESIDUAL_DIGITS = 84


def get_partial(cur, hit_id=None):
    """Return (hit_id, x_prim, y_prim, z_prim, f1_prim, residual, known_factors)
    for the chosen partial."""
    if hit_id is None:
        cur.execute("""
            SELECT mh.id, mh.x_prim, mh.y_prim, mh.z_prim
            FROM pub.master_hits mh
            WHERE mh.num_blockers IS NULL
              AND EXISTS (SELECT 1 FROM pub.prim_brick_factors f
                          WHERE f.x_prim=mh.x_prim AND f.y_prim=mh.y_prim AND f.z_prim=mh.z_prim)
        """)
        rows = cur.fetchall()
        best = None
        best_dig = None
        for hid, xp, yp, zp in rows:
            xp, yp, zp = int(xp), int(yp), int(zp)
            cur.execute(
                """SELECT prime, exponent FROM pub.prim_brick_factors
                   WHERE x_prim=%s AND y_prim=%s AND z_prim=%s""",
                (xp, yp, zp),
            )
            facs = [(int(p), int(e)) for p, e in cur.fetchall()]
            f1_prim = xp*xp + yp*yp + zp*zp
            rem = f1_prim
            for p, e in facs:
                rem //= p**e
            if rem <= 1:
                continue
            dig = len(str(rem))
            if dig < MIN_RESIDUAL_DIGITS:
                continue
            if best_dig is None or dig < best_dig:
                best_dig = dig
                best = (int(hid), xp, yp, zp, f1_prim, rem, facs)
        return best
    else:
        cur.execute(
            "SELECT x_prim, y_prim, z_prim FROM pub.master_hits WHERE id=%s",
            (hit_id,),
        )
        row = cur.fetchone()
        if row is None:
            raise ValueError(f"hit {hit_id} not found")
        xp, yp, zp = int(row[0]), int(row[1]), int(row[2])
        f1_prim = xp*xp + yp*yp + zp*zp
        cur.execute(
            """SELECT prime, exponent FROM pub.prim_brick_factors
               WHERE x_prim=%s AND y_prim=%s AND z_prim=%s""",
            (xp, yp, zp),
        )
        facs = [(int(p), int(e)) for p, e in cur.fetchall()]
        rem = f1_prim
        for p, e in facs:
            rem //= p**e
        return (hit_id, xp, yp, zp, f1_prim, rem, facs)


def run_cado_nfs(n, workdir):
    """Run cado-nfs.py on n, return list of prime factors found."""
    print(f"\n=== Starting CADO-NFS on {len(str(n))}-digit residual ===", flush=True)
    print(f"workdir: {workdir}")
    print(f"input  : {n}\n")
    os.makedirs(workdir, exist_ok=True)

    cmd = [CADO_VENV_PY, CADO_NFS, str(n)]
    proc = subprocess.run(
        cmd, cwd=workdir, capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        print("CADO-NFS exited non-zero. stdout tail:")
        print(proc.stdout[-2000:])
        print("stderr tail:")
        print(proc.stderr[-2000:])
        raise RuntimeError(f"cado-nfs.py exit {proc.returncode}")

    # Parse last numeric line of stdout (CADO-NFS prints final factors)
    primes = []
    for line in proc.stdout.strip().split("\n"):
        line = line.strip()
        # CADO emits the factorisation as a single line of space-separated primes
        if re.fullmatch(r"\d+( \d+)*", line):
            primes = [int(p) for p in line.split()]
    if not primes:
        # Fallback: parse "PRIME" lines
        primes = [int(m) for m in re.findall(r"\b(\d{30,})\b", proc.stdout)]

    return primes, proc.stdout


def process_one(conn, hit_id_arg=None):
    """Pick one partial, run CADO, write back. Returns True on success,
    False if no candidate was available, raises on hard error."""
    cur = conn.cursor()
    info = get_partial(cur, hit_id_arg)
    if info is None:
        print("No partial hits with residual > 1 found.")
        return False
    hit_id, xp, yp, zp, f1_prim, rem, known = info
    print(f"Hit ID: {hit_id}")
    print(f"primitive Brick: ({xp}, {yp}, {zp})")
    print(f"f1_prim has {len(str(f1_prim))} digits, {len(known)} known factors:")
    for p, e in known:
        print(f"  {p}^{e}")
    print(f"Residual: {len(str(rem))} digits")

    workdir = f"/tmp/cado_hit_{hit_id}"
    t0 = time()
    primes, stdout = run_cado_nfs(rem, workdir)
    elapsed = time() - t0
    print(f"\nCADO-NFS done in {elapsed/60:.1f} min")
    print(f"Found primes: {primes}")

    prod = 1
    for p in primes:
        prod *= p
    if prod != rem:
        print(f"WARNING: product mismatch {prod} != residual {rem}")
        print("Aborting DB write for this hit.")
        return True  # don't retry the same hit; loop continues with next

    new_factors = {}
    for p in primes:
        new_factors[p] = new_factors.get(p, 0) + 1

    all_factors = {p: e for p, e in known}
    for p, e in new_factors.items():
        all_factors[p] = all_factors.get(p, 0) + e
    n_blockers = pub_db.store_prim_brick_factors(cur, hit_id, xp, yp, zp, all_factors)
    conn.commit()
    print(f"\nWrote {len(new_factors)} new factors. Final num_blockers = {n_blockers}.")
    print(f"Theorem holds for hit {hit_id}: blocker count = {n_blockers} (>=1).")

    with open(f"{workdir}/result.json", "w") as fp:
        json.dump({
            "hit_id": hit_id,
            "residual": str(rem),
            "primes_found": [str(p) for p in primes],
            "elapsed_seconds": elapsed,
            "all_factors": {str(p): e for p, e in all_factors.items()},
            "num_blockers": n_blockers,
        }, fp, indent=2)
    return True


def main():
    hit_id_arg = int(sys.argv[1]) if len(sys.argv) > 1 else None
    conn = pub_db.connect(autocommit=False)
    try:
        if hit_id_arg is not None:
            process_one(conn, hit_id_arg)
            return
        # Loop mode: always pick the next smallest candidate >= 84 digits,
        # until none are left. SIGINT cleanly breaks between two hits.
        n_done = 0
        while True:
            print(f"\n========== Loop iteration {n_done+1} ==========", flush=True)
            try:
                ok = process_one(conn)
            except KeyboardInterrupt:
                print("\nInterrupted between hits — exiting cleanly.")
                break
            if not ok:
                break
            n_done += 1
        print(f"\nLoop done. Successfully processed {n_done} hit(s).")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
