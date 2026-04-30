"""
Sweep pub.master_hits for perfect cuboids: x² + y² + z² = □.

A perfect cuboid would falsify the conjecture. The check runs in
seconds even on hundreds of thousands of hits, since we only need
isqrt of x²+y²+z² and a comparison.

Usage:
    python3 perfect_check.py
"""
import os
import sys
from math import isqrt
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def main():
    conn = pub_db.connect()
    cur = conn.cursor()
    cur.execute("SELECT id, x, y, z FROM pub.master_hits")
    n = 0
    perfect = []
    t0 = time()
    for hit_id, x, y, z in cur:
        x, y, z = int(x), int(y), int(z)
        d2 = x*x + y*y + z*z
        s = isqrt(d2)
        if s * s == d2:
            perfect.append((hit_id, x, y, z, s))
        n += 1
    print(f"Checked {n} bricks in {time()-t0:.1f}s.")
    print(f"Perfect cuboids found: {len(perfect)}")
    for p in perfect[:10]:
        print(f"  {p}")
    conn.close()
    return 0 if not perfect else 1


if __name__ == "__main__":
    sys.exit(main())
