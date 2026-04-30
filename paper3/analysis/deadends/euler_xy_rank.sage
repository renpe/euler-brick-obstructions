"""
For a sample of (x_prim, y_prim) from pub.master_hits, construct the
elliptic curve E_{x,y} whose rational points correspond to the Euler
bricks (x,y,z) with given (x,y).

Derivation: from u^2 = x^2+z^2 and v^2 = y^2+z^2 we get u^2-v^2 = x^2-y^2,
so (u+v)(u-v) = (x-y)(x+y). Set alpha = u+v, beta = u-v, beta = (x^2-y^2)/alpha. Then
   z^2 * 4 alpha^2 = alpha^4 - 2(x^2+y^2) alpha^2 + (x^2-y^2)^2
With w = 2 alpha z:
   E_{x,y}:  w^2 = alpha^4 - 2*(x^2+y^2)*alpha^2 + (x^2-y^2)^2

Genus 1 (smooth quartic) - an elliptic curve whose rational points
parametrize exactly the Euler z's for (x,y).

Sample:
  - smallest 50 (x,y) from DB (for speed)
  - 38 twin pairs (rank >= 1 expected)
  - 50 random medium size
Output: rank statistics.

Usage:  sage euler_xy_rank.sage
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
    """Returns E_{x,y} or None on error."""
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
    print("Loading sample...")

    # 1. Smallest 50 (x_prim, y_prim) (sorted by max(x,y))
    cur.execute("""
        SELECT x_prim, y_prim, GREATEST(x_prim, y_prim) AS mxy
        FROM pub.master_hits
        GROUP BY x_prim, y_prim
        ORDER BY mxy ASC
        LIMIT 50
    """)
    smallest = [(int(x), int(y)) for x, y, _ in cur.fetchall()]

    # 2. Twin pairs: all (x,y) with >=2 different z's
    cur.execute("""
        SELECT x_prim, y_prim FROM pub.master_hits
        GROUP BY x_prim, y_prim HAVING COUNT(DISTINCT z_prim) > 1
    """)
    twins_xy = [(int(x), int(y)) for x, y in cur.fetchall()]

    # 3. Medium size - 30 samples
    cur.execute("""
        SELECT x_prim, y_prim FROM pub.master_hits
        WHERE x_prim BETWEEN 100 AND 10000
        GROUP BY x_prim, y_prim
        ORDER BY RANDOM()
        LIMIT 30
    """)
    medium = [(int(x), int(y)) for x, y in cur.fetchall()]

    conn.close()

    print(f"Sample: {len(smallest)} smallest, {len(twins_xy)} twins, "
          f"{len(medium)} medium")

    # Canonicalize x <= y and dedup
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

    # Compute ranks
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

    print("\n=== Rank distribution ===")
    for (label, r), c in sorted(rank_dist.items()):
        print(f"  {label} rank={r}: {c}")

    print("\n=== Twin detail ===")
    for x, y, rl, ru in twin_results:
        print(f"  ({x},{y}): rk = {rl}-{ru}")


if __name__ == "__main__":
    main()
