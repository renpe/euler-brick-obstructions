"""
Verify Conjecture C (single-blocker hits are semi-scaled) on the public DB.

For each Master-Hit (a,b,m,n) with num_blockers = 1, compute
    g_+ = gcd(am+bn, an+bm)
    g_- = gcd(|am-bn|, |an-bm|)
and partition by (g_+, g_-) into the four strata
    (k>1, 1)  (1, k>1)  (1, 1)  (both > 1).

Conjecture C asserts: every single-blocker hit lies in (k,1) U (1,k).

Reverse spot-check: on a random sample of M_BLOCK >= 2 hits, every
non-single-blocker hit should have min(g_+, g_-) >= 2.
"""
from __future__ import annotations

import os
import sys
from math import gcd

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(THIS), "_common"))
import pub_db  # type: ignore


def gpm(a: int, b: int, m: int, n: int) -> tuple[int, int]:
    A, B = a * m + b * n, a * n + b * m
    C, D = abs(a * m - b * n), abs(a * n - b * m)
    return gcd(A, B), gcd(C, D)


def main() -> int:
    conn = pub_db.connect()
    cur = conn.cursor()

    cur.execute(
        "SELECT a, b, m, n FROM pub.master_hits WHERE num_blockers = 1"
    )
    rows = cur.fetchall()
    print(f"single-blocker hits: {len(rows)}")

    buckets = {"(k,1)": 0, "(1,k)": 0, "(1,1)": 0, "(both>1)": 0}
    for a, b, m, n in rows:
        gp, gm = gpm(int(a), int(b), int(m), int(n))
        if gp == 1 and gm == 1:
            buckets["(1,1)"] += 1
        elif gp > 1 and gm == 1:
            buckets["(k,1)"] += 1
        elif gp == 1 and gm > 1:
            buckets["(1,k)"] += 1
        else:
            buckets["(both>1)"] += 1

    print("(g_+, g_-) distribution of single-blocker hits:")
    for k, v in buckets.items():
        print(f"  {k:>10s}: {v}")

    semi_scaled = buckets["(k,1)"] + buckets["(1,k)"]
    total = len(rows)
    if semi_scaled == total:
        print(f"OK: all {total} single-blocker hits are semi-scaled.")
        return 0
    else:
        print(
            f"FAIL: {total - semi_scaled} single-blocker hits "
            f"are NOT semi-scaled."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
