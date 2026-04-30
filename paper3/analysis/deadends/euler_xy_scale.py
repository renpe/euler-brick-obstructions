"""
Massiver brute-force-Schnitt-Test: für jedes (x_prim, y_prim) aus pub.master_hits
enumeriere α = num/den mit |num|, |den| ≤ BOUND und prüfe:
  (a) P(num,den) = num⁴ − 2·(x²+y²)·num²·den² + (x²−y²)²·den⁴ ist Quadrat
  (b) Falls ja, liegt der zugehörige z-Wert auch in Cuboid-Position
      (x²+y²+z² = Quadrat)?

Nur (b) interessiert uns: jeder Treffer wäre ein perfekter Cuboid-Kandidat.

Technische Beobachtung: Bricks in pub haben oft α >> BOUND (z.B. (44,117,240)
hat α=511). Brute-force findet also nicht alle existierenden Bricks, sondern
nur diejenigen, deren elliptische Kurve „kleine" rationale Punkte hat. Das
ist trotzdem ein wertvoller Stichprobentest — weil ein perfekter Cuboid mit
einer α-Position bei kleiner Höhe ein „leichter" wäre, den wir bisher nicht
entdeckt hätten.

Aufruf:
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
    """Brute-force über α = num/den. Liefert
       (n_alpha_pts, n_cuboid_candidates, list_of_cuboids)
    list_of_cuboids: [(x, y, z_num, z_den, w, num, den)] für jeden Cuboid-Kandidat.
    """
    x, y, bound = args
    A = x*x + y*y                  # 2·(Koeffizient von α²)/2
    B = (x*x - y*y)**2             # konstanter Term
    n_alpha = 0
    cuboids = []
    # 2A·num²·den² in der inneren Schleife
    # P(num, den) · den⁴/den⁴ = num⁴ − 2A·num²·den² + B·den⁴
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
            # z = w / (2·num·den) (positive root)
            # f1 = x² + y² + z² = (4·num²·den²·A + w²) / (4·num²·den²)
            T = 4*num2*den2*A + P_int
            s = is_perfect_square(T)
            if s is not None:
                # ★★★ Cuboid-Kandidat
                z_num = w
                z_den = 2*num*den
                cuboids.append((x, y, z_num, z_den, w, num, den))
    return (n_alpha, len(cuboids), cuboids)


def main():
    print(f"N_PAIRS={N_PAIRS}, BOUND={BOUND}, WORKERS={WORKERS}\n")
    conn = pub_db.connect()
    cur = conn.cursor()
    print(f"Lade kleinste {N_PAIRS} (x_prim, y_prim) aus DB...")
    cur.execute("""
        SELECT x_prim, y_prim, GREATEST(x_prim, y_prim) AS mxy
        FROM pub.master_hits
        GROUP BY x_prim, y_prim
        ORDER BY mxy ASC
        LIMIT %s
    """, (N_PAIRS,))
    pairs = [(int(x), int(y), BOUND) for x, y, _ in cur.fetchall()]
    conn.close()
    print(f"Geladen: {len(pairs)} Paare\n")

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
                    print(f"  ★★★ CUBOID-KANDIDAT: {cub}", flush=True)
            now = time()
            if now - last_report > 10:
                rate = (i + 1) / (now - t0)
                eta = (len(pairs) - (i + 1)) / rate if rate > 0 else 0
                print(f"  {i+1}/{len(pairs)}  α-Treffer={n_alpha_total}  "
                      f"cuboid={n_cuboid_total}  {rate:.1f} pairs/s  "
                      f"ETA {eta/60:.1f} min", flush=True)
                last_report = now

    print(f"\n========== Zusammenfassung ==========")
    print(f"Geprüft: {len(pairs)} (x,y)-Paare, BOUND={BOUND}")
    print(f"Gesamt α-Treffer (Euler-Bricks proportional zu (x,y)): {n_alpha_total}")
    print(f"Cuboid-Kandidaten: {n_cuboid_total}")
    print(f"Laufzeit: {(time()-t0)/60:.1f} min")
    if all_cuboids:
        print(f"\n★★★ ALARM: {len(all_cuboids)} Cuboid-Kandidaten ★★★")
        for x, y, zn, zd, w, n, d in all_cuboids:
            print(f"  (x={x}, y={y}, z={zn}/{zd}, α={n}/{d})")


if __name__ == "__main__":
    main()
