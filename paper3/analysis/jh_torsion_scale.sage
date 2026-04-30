"""
Scaling of the torsion trick to all (m,n) up to M_MAX.

For each fiber:
  1. Compute rk(E_uV) and rk(E_3).
  2. If at least one is 0: rigorous, close Conjecture B.
  3. Collect statistics.

Usage:  sage jh_torsion_scale.sage [M_MAX]
"""
from sage.all import *
import sys
import signal
from math import gcd as pygcd
from time import time
from collections import Counter

pari.allocatemem(int(2e9))

M_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 50
TIMEOUT = int(45)


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(int(0))


def quartic_to_E(P_quartic):
    R = PolynomialRing(QQ, ['s', 'w'])
    s, w = R.gens()
    T = P_quartic.parent().gen()
    eqn = w**2 - P_quartic.subs({T: s})
    coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
    return EllipticCurve(coeffs)


def safe_rank_lo(E):
    """Return lower bound of rank (or None on error)."""
    try:
        result = with_timeout(TIMEOUT, lambda: E.pari_curve().ellrank())
        return int(result[0]), int(result[1])
    except Exception:
        return None, None


def torsion_size(E):
    try:
        return with_timeout(TIMEOUT, lambda: int(E.torsion_subgroup().order()))
    except Exception:
        return None


def main():
    print(f"Scaling torsion trick to (m,n) up to M_MAX={M_MAX}")
    print(f"We are looking for fibers with rk(E_uV) = 0 OR rk(E_3) = 0.\n")

    print(f"{'(m,n)':>10} {'rk_uV':>6} {'rk_3':>5} {'|tor_uV|':>9} {'|tor_3|':>8} "
          f"{'method':>8} {'bound':>6} {'OK?':>4}")

    n_total = 0
    n_proven = 0
    n_uncertain = 0
    proven_fibers = []
    by_torsion_size = Counter()

    for m in range(2, M_MAX + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1: continue
            if (m - n) % 2 != 1: continue
            n_total += 1

            U2 = m*m - n*n
            V2 = 2*m*n
            W2 = m*m + n*n

            # E_uV: V^2 = (V2^2 u^2 + 4(U2^2-V2^2))(W2^2 u^2 - 4 V2^2)
            Ru = PolynomialRing(QQ, 'U')
            U = Ru.gen()
            quartic_uV = (V2**2 * U**2 + 4*(U2**2 - V2**2)) * (W2**2 * U**2 - 4*V2**2)
            # E_3: V^2 = (V2^2 w^2 + 4 U2^2)(W2^2 w^2 + 4 U2^2)
            Rw = PolynomialRing(QQ, 'W')
            W = Rw.gen()
            quartic_3 = (V2**2 * W**2 + 4*U2**2) * (W2**2 * W**2 + 4*U2**2)

            try:
                E_uV = quartic_to_E(quartic_uV)
                E_3 = quartic_to_E(quartic_3)
            except Exception:
                continue

            rk_uV_lo, rk_uV_hi = safe_rank_lo(E_uV)
            rk_3_lo, rk_3_hi = safe_rank_lo(E_3)

            method = "-"
            tor_uV = None
            tor_3 = None
            bound = None
            ok = "-"

            # We pick the one whose rk is exactly 0 (lo == hi == 0)
            if rk_3_lo == 0 and rk_3_hi == 0:
                tor_3 = torsion_size(E_3)
                if tor_3 is not None:
                    method = "E_3"
                    bound = 2 * tor_3
                    # Empirically we know |H(Q)| >= 8 (6 affine + 2 inf, all trivial)
                    if bound == 8:
                        ok = "yes"
                        n_proven += 1
                        proven_fibers.append((m, n, "E_3", tor_3))
                        by_torsion_size[tor_3] += 1
                    else:
                        ok = f"<={bound}"
                        n_uncertain += 1
            elif rk_uV_lo == 0 and rk_uV_hi == 0:
                tor_uV = torsion_size(E_uV)
                if tor_uV is not None:
                    method = "E_uV"
                    bound = 2 * tor_uV
                    if bound == 8:
                        ok = "yes"
                        n_proven += 1
                        proven_fibers.append((m, n, "E_uV", tor_uV))
                        by_torsion_size[tor_uV] += 1
                    else:
                        ok = f"<={bound}"
                        n_uncertain += 1
            else:
                method = "-"

            tor_uV_str = str(tor_uV) if tor_uV else "?"
            tor_3_str = str(tor_3) if tor_3 else "?"
            print(f"{f'({m},{n})':>10} {rk_uV_lo if rk_uV_lo is not None else '?':>6} "
                  f"{rk_3_lo if rk_3_lo is not None else '?':>5} {tor_uV_str:>9} "
                  f"{tor_3_str:>8} {method:>8} {bound if bound else '?':>6} {ok:>4}")
            sys.stdout.flush()

    print(f"\n========== SUMMARY ==========")
    print(f"Total fibers checked: {n_total}")
    print(f"Conjecture B RIGOROUSLY PROVED: {n_proven}")
    print(f"Uncertain (bound > 8): {n_uncertain}")
    print(f"Distribution by torsion size:")
    for sz, c in sorted(by_torsion_size.items()):
        print(f"  |tors| = {sz}: {c} fibers")

    if proven_fibers:
        print(f"\nProved fibers (m, n):")
        for m, n, method, tor in proven_fibers:
            print(f"  ({m},{n}) via {method} (|tors|={tor})")


if __name__ == "__main__":
    main()
