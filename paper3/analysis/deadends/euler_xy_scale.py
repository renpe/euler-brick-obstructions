"""
Massive brute-force intersection test: for each (x_prim, y_prim) from pub.master_hits
enumerate alpha = num/den with |num|, |den| <= BOUND and check:
  (a) P(num,den) = num^4 - 2*(x^2+y^2)*num^2*den^2 + (x^2-y^2)^2*den^4 is a square
  (b) If yes, does the corresponding z value also lie in cuboid position
      (x^2+y^2+z^2 = square)?

Only (b) interests us: each hit would be a perfect cuboid candidate.

Technical observation: bricks in pub often have alpha >> BOUND (e.g. (44,117,240)
has alpha=511). Brute force therefore does not find all existing bricks but
only those whose elliptic curve has "small" rational points. That is
nevertheless a valuable spot-check - because a perfect cuboid with an alpha
position at small height would be an "easy" one that we have not yet
discovered.

Usage:
    python3 euler_xy_scale.py [N_PAIRS=10000] [BOUND=200] [WORKERS=8]
"""
from __future__ import annotations

import os
import sys
from math import gcd, isqrt
from multiprocessing import Pool
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


N_PAIRS = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
BOUND = int(sys.argv[2]) if len(sys.argv) > 2 else 200
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 8


def is_perfect_square(n: int) -> int | None:
    """Returns sqrt(n) if n is a perfect square, else None."""
    if n < 0:
        return None
    r = isqrt(n)
    if r*r == n:
        return r
    return None


def process_pair(args):
    """Brute-force over alpha = num/den. Returns
       (n_alpha_pts, n_cuboid_candidates, list_of_cuboids)
    list_of_cuboids: [(x, y, z_num, z_den, w, num, den)] for each cuboid candidate.
    """
    x, y, bound = args
    A = x*x + y*y                  # 2*(coefficient of alpha^2)/2
    B = (x*x - y*y)**2             # constant term
    n_alpha = 0
    cuboids = []
    # 2A*num^2*den^2 in the inner loop
    # P(num, den) * den^4/den^4 = num^4 - 2A*num^2*den^2 + B*den^4
    den2_table = [d*d for d in range(bound + 1)]
    den4_table = [d*d for d in den2_table]
    for num in range(1, bound + 1):
        num2 = num*num
        num4 = num2*num2
        for den in range(1, bound + 1):
            if gcd(num, den) != 1:
                continue
            den2 = den2_table[den]
            den4 = den4_table[den]
            P_int = num4 - 2*A*num2*den2 + B*den4
            if P_int <= 0:
                continue
            w = is_perfect_square(P_int)
            if w is None:
                continue
            n_alpha += 1
            # z = w / (2*num*den) (positive root)
            # f1 = x^2 + y^2 + z^2 = (4*num^2*den^2*A + w^2) / (4*num^2*den^2)
            T = 4*num2*den2*A + P_int
            s = is_perfect_square(T)
            if s is not None:
                # *** cuboid candidate
                z_num = w
                z_den = 2*num*den
                cuboids.append((x, y, z_num, z_den, w, num, den))
    return (n_alpha, len(cuboids), cuboids)


def main():
    print(f"N_PAIRS={N_PAIRS}, BOUND={BOUND}, WORKERS={WORKERS}\n")
    conn = pub_db.connect()
    cur = conn.cursor()
    print(f"Loading smallest {N_PAIRS} (x_prim, y_prim) from DB...")
    cur.execute("""
        SELECT x_prim, y_prim, GREATEST(x_prim, y_prim) AS mxy
        FROM pub.master_hits
        GROUP BY x_prim, y_prim
        ORDER BY mxy ASC
        LIMIT %s
    """, (N_PAIRS,))
    pairs = [(int(x), int(y), BOUND) for x, y, _ in cur.fetchall()]
    conn.close()
    print(f"Loaded: {len(pairs)} pairs\n")

    t0 = time()
    n_alpha_total = 0
    n_cuboid_total = 0
    all_cuboids = []
    last_report = t0

    with Pool(processes=WORKERS) as pool:
        for i, (n_alpha, n_cub, cubs) in enumerate(
                pool.imap_unordered(process_pair, pairs, chunksize=4)):
            n_alpha_total += n_alpha
            n_cuboid_total += n_cub
            if cubs:
                all_cuboids.extend(cubs)
                for cub in cubs:
                    print(f"  *** CUBOID CANDIDATE: {cub}", flush=True)
            now = time()
            if now - last_report > 10:
                rate = (i + 1) / (now - t0)
                eta = (len(pairs) - (i + 1)) / rate if rate > 0 else 0
                print(f"  {i+1}/{len(pairs)}  alpha hits={n_alpha_total}  "
                      f"cuboid={n_cuboid_total}  {rate:.1f} pairs/s  "
                      f"ETA {eta/60:.1f} min", flush=True)
                last_report = now

    print(f"\n========== Summary ==========")
    print(f"Checked: {len(pairs)} (x,y) pairs, BOUND={BOUND}")
    print(f"Total alpha hits (Euler bricks proportional to (x,y)): {n_alpha_total}")
    print(f"Cuboid candidates: {n_cuboid_total}")
    print(f"Runtime: {(time()-t0)/60:.1f} min")
    if all_cuboids:
        print(f"\n*** ALARM: {len(all_cuboids)} cuboid candidates ***")
        for x, y, zn, zd, w, n, d in all_cuboids:
            print(f"  (x={x}, y={y}, z={zn}/{zd}, alpha={n}/{d})")


if __name__ == "__main__":
    main()
