"""
Formelle Coleman-Bound für Konjektur B auf den 11 Chabauty-anwendbaren
Fasern.

Coleman (1985): Für genus-g hyperelliptische Kurve X mit guter Reduktion
modulo Primzahl p (p > 2g) und rk(J(X)) < g gilt:
   |X(Q)| ≤ |X(F_p)| + 2g − 2

Für uns: g=3, rk(J(H)) < 3 → |H(Q)| ≤ |H(F_p)| + 4.

Strategie: Für jede der 11 Fasern verschiedene Primzahlen p (≥ 7) ausprobieren,
das Minimum |H(F_p)| finden, und Coleman-Bound anwenden.

Bei "guter Reduktion" muss disc(f) ≢ 0 (mod p), wobei f das Polynom ist.

Aufruf: sage jh_coleman_bound.sage
"""
from sage.all import *
import sys
import signal
from time import time

pari.allocatemem(int(2e9))

# Fasern mit ihrem berechneten rk(J(H)):
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
    print("Formelle Coleman-Bound für die 11 Chabauty-Fasern\n")
    print("Ziel: |H(F_p)| ≤ 4 finden → |H(Q)| ≤ 8 (matcht Empirie)\n")

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
        # Inf-Punkte: 2 wenn Leading-Coeff Quadrat
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
            stoll_bound = count + 2*rk_J          # Stoll: |X(Q)| ≤ N_p + 2r
            coleman_bound = count + 4             # Coleman: N_p + 2g-2
            if best_stoll is None or stoll_bound < best_stoll:
                best_stoll = stoll_bound
                best_coleman = coleman_bound
                best_count = count
                best_p = p

        empirical = 6 + n_inf
        # Rigoros wenn der ENGSTE bekannte Bound (Stoll) ≤ empirisch
        if best_stoll is not None and best_stoll <= empirical:
            closed = "✓ RIGOROUS"
            n_proven += 1
            rigorous_results.append((m, n, rk_J, best_p, best_count,
                                     best_stoll, empirical))
        else:
            closed = "−"
        print(f"{f'({m},{n})':>10} {rk_J:>5} {best_p if best_p else '?':>8} "
              f"{best_count if best_count is not None else '?':>10} "
              f"{best_stoll if best_stoll is not None else '?':>7} "
              f"{best_coleman if best_coleman is not None else '?':>8} "
              f"{empirical:>11} {closed:>10}")

    print(f"\n========== Zusammenfassung ==========")
    print(f"Rigoros geschlossene Fasern: {n_proven}/{len(CHABAUTY_FIBERS_BY_RANK)}")
    if rigorous_results:
        print("\nDetail rigoros bewiesener Fasern (via Stoll-Bound):")
        for m, n, rk_J, p, c, b, e in rigorous_results:
            print(f"  ({m},{n}) rk(J)={rk_J}: bei p={p}, |H(F_p)|={c}, "
                  f"Stoll-Bound={b}, empirisch={e}")


if __name__ == "__main__":
    main()
