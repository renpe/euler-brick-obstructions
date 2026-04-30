"""
Exploration: how do hits from a (g_+, g_-) fiber distribute over (m, n)?

We pick e.g. (g_+ = 9, g_- = 1) - the largest Saund+ fiber with 14086 hits -
and look whether the points are concentrated on few (m,n) sub-fibers or
broadly scattered. If concentrated: we can set up an elliptic (sub)curve
per (g_+, g_-, m, n) and compute Mordell-Weil.

Usage:
    python3 gpm_subfiber_explore.py [g_plus] [g_minus]
    Default: 9 1
"""
from __future__ import annotations

import os
import sys
from collections import Counter
from math import gcd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def main():
    target_gp = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    target_gm = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    conn = pub_db.connect()
    cur = conn.cursor()
    cur.execute("""
        SELECT id, a, b, m, n, num_blockers
        FROM pub.master_hits
    """)
    rows = cur.fetchall()

    in_fiber = []
    for hid, a, b, m, n, nb in rows:
        a, b, m, n = int(a), int(b), int(m), int(n)
        A = a*m + b*n; B = a*n + b*m
        C = a*m - b*n; D = a*n - b*m
        gp = gcd(A, B); gm = gcd(abs(C), abs(D))
        if gp == target_gp and gm == target_gm:
            in_fiber.append((int(hid), a, b, m, n, nb))

    print(f"(g_+={target_gp}, g_-={target_gm}): {len(in_fiber)} hits.\n")

    # --- Distribution over (m, n) -----------------------------------
    mn_counter = Counter((m, n) for _, _, _, m, n, _ in in_fiber)
    print(f"Number of distinct (m, n) sub-fibers: {len(mn_counter)}")
    print(f"Maximum hits per (m, n):              {max(mn_counter.values())}")
    print(f"Median hits per (m, n):               "
          f"{sorted(mn_counter.values())[len(mn_counter)//2]}")
    print(f"Sub-fibers with only 1 hit:           "
          f"{sum(1 for v in mn_counter.values() if v == 1)}")
    print(f"Sub-fibers with >=10 hits:            "
          f"{sum(1 for v in mn_counter.values() if v >= 10)}")
    print()

    print("=== Top 15 (m, n) sub-fibers in this fiber ===")
    print(f"{'m':>5} {'n':>5} {'#hits':>6} {'min_b':>6} {'max_b':>6}")
    for (m, n), c in mn_counter.most_common(15):
        bs = [int(nb) for hid, a, b, mm, nn, nb in in_fiber
              if mm == m and nn == n and nb is not None]
        if bs:
            print(f"{m:>5} {n:>5} {c:>6} {min(bs):>6} {max(bs):>6}")
        else:
            print(f"{m:>5} {n:>5} {c:>6} {'-':>6} {'-':>6}")

    # --- Distribution over (a, b) -----------------------------------
    print()
    ab_counter = Counter((a, b) for _, a, b, _, _, _ in in_fiber)
    print(f"Number of distinct (a, b) values:   {len(ab_counter)}")
    print(f"Top (a, b) with most hits (= different m,n):")
    for (a, b), c in ab_counter.most_common(10):
        print(f"  (a,b)=({a},{b}): {c} different (m,n)")

    # --- Single blockers within this fiber --------------------------
    sb = [(hid, a, b, m, n) for hid, a, b, m, n, nb in in_fiber if nb == 1]
    print(f"\n=== Single blockers in (g_+={target_gp}, g_-={target_gm}) ===")
    print(f"Count: {len(sb)}")
    for hid, a, b, m, n in sb[:20]:
        print(f"  hit_id={hid}: (a,b,m,n)=({a},{b},{m},{n})")

    conn.close()


if __name__ == "__main__":
    main()
