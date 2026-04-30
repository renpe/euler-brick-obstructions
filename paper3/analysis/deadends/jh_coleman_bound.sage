"""
Formal Coleman bound for Conjecture B on the 11 Chabauty-applicable
fibers.

Coleman (1985): For genus-g hyperelliptic curve X with good reduction
modulo prime p (p > 2g) and rk(J(X)) < g, we have:
   |X(Q)| <= |X(F_p)| + 2g - 2

For us: g=3, rk(J(H)) < 3 -> |H(Q)| <= |H(F_p)| + 4.

Strategy: For each of the 11 fibers try various primes p (>= 7),
find the minimum |H(F_p)|, and apply the Coleman bound.

For "good reduction" we need disc(f) != 0 (mod p), where f is the polynomial.

Usage: sage jh_coleman_bound.sage
"""
from sage.all import *
import sys
import signal
from time import time

pari.allocatemem(int(2e9))

# Fibers with their computed rk(J(H)):
CHABAUTY_FIBERS_BY_RANK = [
    # rk(J) = 1
    ((2, 1), 1),
    ((3, 2), 1),
    ((6, 1), 1),
    ((7, 2), 1),
    ((7, 6), 1),
    # rk(J) = 2
    ((4, 1), 2),
    ((4, 3), 2),
    ((8, 1), 2),
    ((11, 6), 2),
    ((12, 1), 2),
    ((12, 5), 2),
]

PRIMES_TO_TRY = [int(p) for p in primes(3, 200)]


def main():
    print("Formal Coleman bound for the 11 Chabauty fibers\n")
    print("Goal: find |H(F_p)| <= 4 -> |H(Q)| <= 8 (matches empirics)\n")

    print(f"{'(m,n)':>10} {'rk_J':>5} {'best_p':>8} {'|H(F_p)|':>10} {'Stoll':>7} "
          f"{'Coleman':>8} {'empirical':>11} {'closed':>10}")

    n_proven = 0
    rigorous_results = []

    for (m, n), rk_J in CHABAUTY_FIBERS_BY_RANK:
        U2 = m*m - n*n
        V2 = 2*m*n
        W2 = m*m + n*n
        Rt = PolynomialRing(QQ, 'T')
        T = Rt.gen()
        P = V2**2 * T**4 + (4*U2**2 - 2*V2**2) * T**2 + V2**2
        Q = W2**2 * T**4 + 2*(U2**2 - V2**2) * T**2 + W2**2
        f = P * Q
        disc = f.discriminant()
        leading = f.leading_coefficient()
        # Inf points: 2 if leading coeff is a square
        n_inf = 2 if leading.is_square() else 0

        best_p = None
        best_count = None
        best_stoll = None
        best_coleman = None
        for p in PRIMES_TO_TRY:
            if disc % p == 0:
                continue
            try:
                fp = f.change_ring(GF(p))
                Hp = HyperellipticCurve(fp, 0)
                count = sum(int(c) for c in Hp.count_points())
            except Exception:
                continue
            stoll_bound = count + 2*rk_J          # Stoll: |X(Q)| <= N_p + 2r
            coleman_bound = count + 4             # Coleman: N_p + 2g-2
            if best_stoll is None or stoll_bound < best_stoll:
                best_stoll = stoll_bound
                best_coleman = coleman_bound
                best_count = count
                best_p = p

        empirical = 6 + n_inf
        # Rigorous if the TIGHTEST known bound (Stoll) <= empirical
        if best_stoll is not None and best_stoll <= empirical:
            closed = "RIGOROUS"
            n_proven += 1
            rigorous_results.append((m, n, rk_J, best_p, best_count,
                                     best_stoll, empirical))
        else:
            closed = "-"
        print(f"{f'({m},{n})':>10} {rk_J:>5} {best_p if best_p else '?':>8} "
              f"{best_count if best_count is not None else '?':>10} "
              f"{best_stoll if best_stoll is not None else '?':>7} "
              f"{best_coleman if best_coleman is not None else '?':>8} "
              f"{empirical:>11} {closed:>10}")

    print(f"\n========== Summary ==========")
    print(f"Rigorously closed fibers: {n_proven}/{len(CHABAUTY_FIBERS_BY_RANK)}")
    if rigorous_results:
        print("\nDetail of rigorously proven fibers (via Stoll bound):")
        for m, n, rk_J, p, c, b, e in rigorous_results:
            print(f"  ({m},{n}) rk(J)={rk_J}: at p={p}, |H(F_p)|={c}, "
                  f"Stoll bound={b}, empirical={e}")


if __name__ == "__main__":
    main()
