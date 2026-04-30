"""
Empirical multi-angle sweep over the 1.28M primitive bricks (x_prim, y_prim,
z_prim). Looks for structures that may have been overlooked in the literature.

Angle A - f1_prim multiplicity: how often do bricks share the same
            space diagonal squared?
Angle B - (x, y) -> z multiplicity: at fixed (x_prim, y_prim), how many
            z's occur?
Angle C - proximity to "perfect": for each brick check the nearest integer to
            sqrt(x^2+y^2) resp. sqrt(x^2+z^2) resp. sqrt(y^2+z^2).
Angle D - Projective distribution in P^2: clustering on rational curves?
Angle E - Heron triangle consistency: for each brick check the three
            face diagonal squares (x^2+y^2, x^2+z^2, y^2+z^2).
"""
from __future__ import annotations

import os
import sys
from collections import Counter, defaultdict
from math import gcd, isqrt
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def is_square(n):
    if n < 0:
        return False
    r = isqrt(n)
    return r * r == n


def main():
    conn = pub_db.connect()
    cur = conn.cursor()
    print("Loading bricks...")
    cur.execute("SELECT DISTINCT x_prim, y_prim, z_prim FROM pub.master_hits")
    bricks = [(int(x), int(y), int(z)) for x, y, z in cur.fetchall()]
    print(f"Distinct primitive bricks: {len(bricks)}\n")

    t0 = time()

    # =========== Angle A: f1_prim multiplicity =====================
    print("=== A. f1_prim multiplicity ===")
    f1_count = Counter()
    for x, y, z in bricks:
        f1 = x*x + y*y + z*z
        f1_count[f1] += 1
    multi = [(f1, c) for f1, c in f1_count.items() if c > 1]
    print(f"f1_prim occurring >1 times: {len(multi)}")
    if multi:
        multi.sort(key=lambda kv: -kv[1])
        print(f"Top-10 f1_prim with most frequent brick count:")
        for f1, c in multi[:10]:
            print(f"  f1={f1} ({len(str(f1))} digits): {c} bricks")
    print(f"Maximum multiplicity: {max(f1_count.values())}\n")

    # =========== Angle B: (x,y) -> z multiplicity ===================
    print("=== B. (x_prim, y_prim) -> z_prim multiplicity ===")
    xy_to_z = defaultdict(set)
    for x, y, z in bricks:
        # Sorted pairs, so that (x,y) ~ (y,x)
        a, b = min(x, y), max(x, y)
        xy_to_z[(a, b)].add(z)
    multi_xy = [(xy, zs) for xy, zs in xy_to_z.items() if len(zs) > 1]
    print(f"(x,y) pairs with >1 z: {len(multi_xy)} of {len(xy_to_z)}")
    if multi_xy:
        multi_xy.sort(key=lambda kv: -len(kv[1]))
        print(f"Top-10 (x,y) with most z's:")
        for (a, b), zs in multi_xy[:10]:
            print(f"  ({a},{b}): z in {sorted(zs)[:5]}{'...' if len(zs) > 5 else ''} "
                  f"(total {len(zs)})")
    print(f"Median z count per (x,y): "
          f"{sorted(len(zs) for zs in xy_to_z.values())[len(xy_to_z)//2]}\n")

    # =========== Angle C: face diagonal squares =====================
    print("=== C. Face diagonal squares (Pythagoras test) ===")
    face_xy_sq = 0
    face_xz_sq = 0
    face_yz_sq = 0
    two_faces_sq = 0
    three_faces_sq = 0
    sample_two = []
    for x, y, z in bricks[:200000]:  # Sample 200K for speed
        a = is_square(x*x + y*y)
        b = is_square(x*x + z*z)
        c = is_square(y*y + z*z)
        if a:
            face_xy_sq += 1
        if b:
            face_xz_sq += 1
        if c:
            face_yz_sq += 1
        n = sum([a, b, c])
        if n >= 2:
            two_faces_sq += 1
            if len(sample_two) < 5:
                sample_two.append((x, y, z, a, b, c))
        if n == 3:
            three_faces_sq += 1
    print(f"(Sample 200K)")
    print(f"  x^2+y^2 a square:  {face_xy_sq:>8}")
    print(f"  x^2+z^2 a square:  {face_xz_sq:>8}")
    print(f"  y^2+z^2 a square:  {face_yz_sq:>8}")
    print(f"  >=2 squares (Heron trio): {two_faces_sq}")
    print(f"  All 3 squares (Euler brick): {three_faces_sq}")
    if sample_two:
        print(f"Example bricks with >=2 square faces:")
        for x, y, z, a, b, c in sample_two:
            print(f"  ({x},{y},{z}) - xy^2={a}, xz^2={b}, yz^2={c}")
    print()

    # =========== Angle D: gcd structure =============================
    print("=== D. gcd structures ===")
    gcd_xy = Counter()
    gcd_xz = Counter()
    gcd_yz = Counter()
    gcd_all3 = Counter()
    for x, y, z in bricks:
        gcd_xy[gcd(x, y)] += 1
        gcd_xz[gcd(x, z)] += 1
        gcd_yz[gcd(y, z)] += 1
        gcd_all3[gcd(gcd(x, y), z)] += 1
    print(f"gcd(x,y)=1: {gcd_xy[1]}; !=1: {sum(c for g,c in gcd_xy.items() if g!=1)}")
    print(f"gcd(x,z)=1: {gcd_xz[1]}; !=1: {sum(c for g,c in gcd_xz.items() if g!=1)}")
    print(f"gcd(y,z)=1: {gcd_yz[1]}; !=1: {sum(c for g,c in gcd_yz.items() if g!=1)}")
    print(f"gcd(x,y,z): only g=1 should appear (primitive!)")
    for g, c in sorted(gcd_all3.items())[:5]:
        print(f"  gcd(x,y,z)={g}: {c}")
    print()

    # =========== Angle E: log scaling ===============================
    print("=== E. Scaling log(f1_prim) vs log(xyz) ===")
    import math
    ratios = []
    for x, y, z in bricks[:50000]:
        f1 = x*x + y*y + z*z
        ratio = math.log(f1) / math.log(x*y*z)
        ratios.append(ratio)
    ratios.sort()
    print(f"  log(f1)/log(xyz): min={ratios[0]:.3f}, median={ratios[len(ratios)//2]:.3f}, "
          f"max={ratios[-1]:.3f}")
    # Expected: ~2/3 if (x,y,z) of equal size
    print(f"  (theoretical value for x~y~z: 2/3 = 0.667; greater discrepancy "
          f"= one-sided dominated)")

    print(f"\nDone in {time()-t0:.1f}s.")
    conn.close()


if __name__ == "__main__":
    main()
