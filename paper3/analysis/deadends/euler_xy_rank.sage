"""
Konstruiert für Stichprobe von (x_prim, y_prim) aus pub.master_hits die
elliptische Kurve E_{x,y}, deren rationale Punkte den Euler-Bricks (x,y,z)
mit gegebenem (x,y) entsprechen.

Herleitung: Aus u² = x²+z² und v² = y²+z² folgt u²-v² = x²-y², also
(u+v)(u-v) = (x-y)(x+y). Setze α = u+v, β = u-v, β = (x²-y²)/α. Dann
   z² · 4α² = α⁴ - 2(x²+y²)α² + (x²-y²)²
Mit w = 2αz:
   E_{x,y}:  w² = α⁴ − 2·(x²+y²)·α² + (x²−y²)²

Genus 1 (smooth quartic) — eine elliptische Kurve, deren rationale Punkte
genau die Euler-z's zu (x,y) parametrisieren.

Stichprobe:
  - Kleinste 50 (x,y) aus DB (für Speed)
  - 38 Twin-Paare (rank ≥ 1 erwartet)
  - 50 zufällige mittlere Größe
Output: Rangstatistik.

Aufruf:  sage euler_xy_rank.sage
"""
from sage.all import *
import psycopg
import sys
import os
import signal
from collections import Counter, defaultdict
from math import gcd as pygcd
from time import time

pari.allocatemem(int(2e9))

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"
TIMEOUT = int(30)


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(int(0))


def build_E(x, y):
    """Liefert E_{x,y} oder None bei Fehler."""
    A = x*x + y*y
    B = (x*x - y*y)**2
    R = PolynomialRing(QQ, ['s', 'w'])
    s, w = R.gens()
    eqn = w**2 - (s**4 - 2*A*s**2 + B)
    coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
    return EllipticCurve(coeffs)


def rank_info(E):
    try:
        result = with_timeout(TIMEOUT, lambda: E.pari_curve().ellrank())
        return int(result[0]), int(result[1])
    except Exception:
        return None, None


def main():
    conn = psycopg.connect(DB)
    cur = conn.cursor()
    print("Lade Stichprobe...")

    # 1. Kleinste 50 (x_prim, y_prim) (sortiert nach max(x,y))
    cur.execute("""
        SELECT x_prim, y_prim, GREATEST(x_prim, y_prim) AS mxy
        FROM pub.master_hits
        GROUP BY x_prim, y_prim
        ORDER BY mxy ASC
        LIMIT 50
    """)
    smallest = [(int(x), int(y)) for x, y, _ in cur.fetchall()]

    # 2. Twin-Paare: alle (x,y) mit ≥2 verschiedenen z's
    cur.execute("""
        SELECT x_prim, y_prim FROM pub.master_hits
        GROUP BY x_prim, y_prim HAVING COUNT(DISTINCT z_prim) > 1
    """)
    twins_xy = [(int(x), int(y)) for x, y in cur.fetchall()]

    # 3. Mittlere Größe — 30 Stichproben
    cur.execute("""
        SELECT x_prim, y_prim FROM pub.master_hits
        WHERE x_prim BETWEEN 100 AND 10000
        GROUP BY x_prim, y_prim
        ORDER BY RANDOM()
        LIMIT 30
    """)
    medium = [(int(x), int(y)) for x, y in cur.fetchall()]

    conn.close()

    print(f"Stichprobe: {len(smallest)} smallest, {len(twins_xy)} twins, "
          f"{len(medium)} medium")

    # Vereinheitliche x ≤ y und dedup
    def norm(xy):
        x, y = xy
        return (min(x, y), max(x, y))

    samples = []
    seen = set()
    for label, lst in [("small", smallest), ("twin", twins_xy), ("medium", medium)]:
        for xy in lst:
            xy = norm(xy)
            key = (label, xy)
            if xy in seen:
                continue
            seen.add(xy)
            samples.append((label, xy))

    print(f"Total dedupe: {len(samples)}\n")

    # Berechne ranks
    rank_dist = Counter()
    twin_results = []
    print(f"{'label':<8} {'x':>10} {'y':>10} {'rk_lo':>5} {'rk_hi':>5} "
          f"{'time':>6}")
    for label, (x, y) in samples:
        t0 = time()
        try:
            E = build_E(x, y)
        except Exception as ex:
            print(f"{label:<8} {x:>10} {y:>10} curve_fail: {ex}")
            continue
        rl, ru = rank_info(E)
        elapsed = time() - t0
        if rl is None:
            r_str = "?-?"
            key = (label, "?")
        else:
            r_str = f"{rl}-{ru}"
            key = (label, rl)
        rank_dist[key] += 1
        if label == "twin":
            twin_results.append((x, y, rl, ru))
        print(f"{label:<8} {x:>10} {y:>10} "
              f"{rl if rl is not None else '?':>5} "
              f"{ru if ru is not None else '?':>5} "
              f"{elapsed:>6.1f}")
        sys.stdout.flush()

    print("\n=== Rang-Verteilung ===")
    for (label, r), c in sorted(rank_dist.items()):
        print(f"  {label} rank={r}: {c}")

    print("\n=== Twin-Detail ===")
    for x, y, rl, ru in twin_results:
        print(f"  ({x},{y}): rk = {rl}-{ru}")


if __name__ == "__main__":
    main()
