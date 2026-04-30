"""
Empirischer Multi-Angle-Sweep über die 1.28M primitiven Bricks (x_prim, y_prim,
z_prim). Sucht nach Strukturen, die in der Literatur evtl. übersehen wurden.

Winkel A — f1_prim-Multiplizität: wie oft teilen sich Bricks dieselbe
            Raumdiagonale²?
Winkel B — (x, y) → z Multiplizität: bei festem (x_prim, y_prim), wie viele
            z's kommen vor?
Winkel C — Nähe zu „perfekt": für jeden Brick die nächste ganze Zahl zu
            √(x²+y²) bzw. √(x²+z²) bzw. √(y²+z²) prüfen.
Winkel D — Projektive Verteilung in ℙ²: Clustern auf rationalen Kurven?
Winkel E — Heron-Triangel-Konsistenz: für jeden Brick die drei
            Face-Diagonal-Quadrate (x²+y², x²+z², y²+z²) prüfen.
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
    print("Lade Bricks...")
    cur.execute("SELECT DISTINCT x_prim, y_prim, z_prim FROM pub.master_hits")
    bricks = [(int(x), int(y), int(z)) for x, y, z in cur.fetchall()]
    print(f"Distinkte primitive Bricks: {len(bricks)}\n")

    t0 = time()

    # =========== Winkel A: f1_prim-Multiplizität ====================
    print("=== A. f1_prim-Multiplizität ===")
    f1_count = Counter()
    for x, y, z in bricks:
        f1 = x*x + y*y + z*z
        f1_count[f1] += 1
    multi = [(f1, c) for f1, c in f1_count.items() if c > 1]
    print(f"f1_prim die >1 Mal vorkommen: {len(multi)}")
    if multi:
        multi.sort(key=lambda kv: -kv[1])
        print(f"Top-10 f1_prim mit häufigster Brick-Anzahl:")
        for f1, c in multi[:10]:
            print(f"  f1={f1} ({len(str(f1))} Stellen): {c} Bricks")
    print(f"Maximum-Multiplizität: {max(f1_count.values())}\n")

    # =========== Winkel B: (x,y) → z Multiplizität ==================
    print("=== B. (x_prim, y_prim) → z_prim Multiplizität ===")
    xy_to_z = defaultdict(set)
    for x, y, z in bricks:
        # Sortierte Paare, damit (x,y)~(y,x)
        a, b = min(x, y), max(x, y)
        xy_to_z[(a, b)].add(z)
    multi_xy = [(xy, zs) for xy, zs in xy_to_z.items() if len(zs) > 1]
    print(f"(x,y)-Paare mit >1 z: {len(multi_xy)} von {len(xy_to_z)}")
    if multi_xy:
        multi_xy.sort(key=lambda kv: -len(kv[1]))
        print(f"Top-10 (x,y) mit den meisten z's:")
        for (a, b), zs in multi_xy[:10]:
            print(f"  ({a},{b}): z ∈ {sorted(zs)[:5]}{'...' if len(zs) > 5 else ''} "
                  f"(insgesamt {len(zs)})")
    print(f"Median z-Anzahl pro (x,y): "
          f"{sorted(len(zs) for zs in xy_to_z.values())[len(xy_to_z)//2]}\n")

    # =========== Winkel C: Face-Diagonal-Quadrate ===================
    print("=== C. Face-Diagonal-Quadrate (Pythagoras-Test) ===")
    face_xy_sq = 0
    face_xz_sq = 0
    face_yz_sq = 0
    two_faces_sq = 0
    three_faces_sq = 0
    sample_two = []
    for x, y, z in bricks[:200000]:  # Sample 200K für Speed
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
    print(f"  x²+y² ein Quadrat:  {face_xy_sq:>8}")
    print(f"  x²+z² ein Quadrat:  {face_xz_sq:>8}")
    print(f"  y²+z² ein Quadrat:  {face_yz_sq:>8}")
    print(f"  ≥2 Quadrate (Heron-Trio): {two_faces_sq}")
    print(f"  Alle 3 Quadrate (Euler-Brick): {three_faces_sq}")
    if sample_two:
        print(f"Beispiel-Bricks mit ≥2 quadrischen Faces:")
        for x, y, z, a, b, c in sample_two:
            print(f"  ({x},{y},{z}) — xy²={a}, xz²={b}, yz²={c}")
    print()

    # =========== Winkel D: ggT-Struktur =============================
    print("=== D. ggT-Strukturen ===")
    gcd_xy = Counter()
    gcd_xz = Counter()
    gcd_yz = Counter()
    gcd_all3 = Counter()
    for x, y, z in bricks:
        gcd_xy[gcd(x, y)] += 1
        gcd_xz[gcd(x, z)] += 1
        gcd_yz[gcd(y, z)] += 1
        gcd_all3[gcd(gcd(x, y), z)] += 1
    print(f"gcd(x,y)=1: {gcd_xy[1]}; ≠1: {sum(c for g,c in gcd_xy.items() if g!=1)}")
    print(f"gcd(x,z)=1: {gcd_xz[1]}; ≠1: {sum(c for g,c in gcd_xz.items() if g!=1)}")
    print(f"gcd(y,z)=1: {gcd_yz[1]}; ≠1: {sum(c for g,c in gcd_yz.items() if g!=1)}")
    print(f"gcd(x,y,z): only g=1 should appear (primitive!)")
    for g, c in sorted(gcd_all3.items())[:5]:
        print(f"  gcd(x,y,z)={g}: {c}")
    print()

    # =========== Winkel E: log-Skalierung ===========================
    print("=== E. Skalierung log(f1_prim) vs log(xyz) ===")
    import math
    ratios = []
    for x, y, z in bricks[:50000]:
        f1 = x*x + y*y + z*z
        ratio = math.log(f1) / math.log(x*y*z)
        ratios.append(ratio)
    ratios.sort()
    print(f"  log(f1)/log(xyz): min={ratios[0]:.3f}, median={ratios[len(ratios)//2]:.3f}, "
          f"max={ratios[-1]:.3f}")
    # Erwartet: ~2/3 wenn (x,y,z) gleich groß
    print(f"  (theoretischer Wert für x≈y≈z: 2/3 = 0.667; größere Diskrepanz "
          f"= einseitig dominiert)")

    print(f"\nFertig in {time()-t0:.1f}s.")
    conn.close()


if __name__ == "__main__":
    main()
