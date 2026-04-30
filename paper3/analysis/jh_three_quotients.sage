"""
Vollständige Analyse: für ALLE 518 Fasern bei M_MAX=50 die Ränge von
E_PQ, E_uV, E_3 bestimmen und die Abdeckung der Torsions-Trick-Methode
quantifizieren.

Eine Faser ist mit Torsions-Trick rigoros lösbar wenn:
  rk(E_PQ) = 0 ODER rk(E_uV) = 0 ODER rk(E_3) = 0
und |E_q(Q)_tors| = 4 (für den rk=0-Quotienten).

Wenn ALLE drei Ränge ≥ 1 sind: Chabauty greift NUR wenn Total-Rang < 3.

Aufruf: sage jh_three_quotients.sage [M_MAX]
"""
from sage.all import *
import sys
import signal
from math import gcd as pygcd
from collections import Counter
from time import time

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


def safe_rank(E):
    try:
        result = with_timeout(TIMEOUT, lambda: E.pari_curve().ellrank())
        return int(result[0]), int(result[1])
    except Exception:
        return None, None


def safe_torsion(E):
    try:
        return with_timeout(TIMEOUT, lambda: int(E.torsion_subgroup().order()))
    except Exception:
        return None


def main():
    print(f"Drei-Quotienten-Analyse für (m,n) bis M_MAX={M_MAX}\n")

    n_total = 0
    n_torsion_trick = 0   # mindestens ein rk=0 quotient + |tors|=4
    n_chabauty_only = 0    # total rk < 3 aber kein rk=0
    n_uncovered = 0        # total rk ≥ 3
    by_pattern = Counter()
    proven_torsion = []
    chabauty_only = []
    uncovered = []

    for m in range(2, M_MAX + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1: continue
            if (m - n) % 2 != 1: continue
            n_total += 1

            U2 = m*m - n*n
            V2 = 2*m*n
            W2 = m*m + n*n

            # Drei Quartiken
            Rs = PolynomialRing(QQ, 'S')
            S = Rs.gen()
            P_s = V2**2 * S**2 + (4*U2**2 - 2*V2**2) * S + V2**2
            Q_s = W2**2 * S**2 + 2*(U2**2 - V2**2) * S + W2**2
            quartic_PQ = P_s * Q_s

            Ru = PolynomialRing(QQ, 'U')
            U = Ru.gen()
            quartic_uV = (V2**2 * U**2 + 4*(U2**2 - V2**2)) * (W2**2 * U**2 - 4*V2**2)

            Rw = PolynomialRing(QQ, 'W')
            W = Rw.gen()
            quartic_3 = (V2**2 * W**2 + 4*U2**2) * (W2**2 * W**2 + 4*U2**2)

            try:
                E_PQ = quartic_to_E(quartic_PQ)
                E_uV = quartic_to_E(quartic_uV)
                E_3 = quartic_to_E(quartic_3)
            except Exception:
                continue

            rk_PQ_lo, rk_PQ_hi = safe_rank(E_PQ)
            rk_uV_lo, rk_uV_hi = safe_rank(E_uV)
            rk_3_lo, rk_3_hi = safe_rank(E_3)

            if any(r is None for r in [rk_PQ_lo, rk_uV_lo, rk_3_lo]):
                continue

            # Pattern-Klassifikation
            pattern = (rk_PQ_lo, rk_uV_lo, rk_3_lo)
            by_pattern[pattern] += 1

            # Welcher Quotient hat rk=0 mit eindeutigem rank (lo == hi)?
            method = None
            quotient = None
            if rk_3_lo == 0 and rk_3_hi == 0:
                method = "E_3"
                quotient = E_3
            elif rk_uV_lo == 0 and rk_uV_hi == 0:
                method = "E_uV"
                quotient = E_uV
            elif rk_PQ_lo == 0 and rk_PQ_hi == 0:
                method = "E_PQ"
                quotient = E_PQ

            if method is not None:
                tors = safe_torsion(quotient)
                if tors == 4:
                    n_torsion_trick += 1
                    proven_torsion.append((m, n, method, tors))
                else:
                    # rk=0 aber |tors|≠4: feinere Analyse nötig
                    pass
            else:
                # Alle Ränge ≥ 1
                total_rk_lo = rk_PQ_lo + rk_uV_lo + rk_3_lo
                total_rk_hi = rk_PQ_hi + rk_uV_hi + rk_3_hi
                if total_rk_hi < 3:
                    n_chabauty_only += 1
                    chabauty_only.append((m, n, pattern))
                else:
                    n_uncovered += 1
                    uncovered.append((m, n, pattern, total_rk_lo, total_rk_hi))

    print(f"\n========== ZUSAMMENFASSUNG ==========")
    print(f"Gesamt geprüfte Fasern: {n_total}")
    print(f"  Torsions-Trick (|tors|=4):           {n_torsion_trick}")
    print(f"  Nur via Chabauty (Total-Rang < 3):   {n_chabauty_only}")
    print(f"  Unabgedeckt (Total-Rang ≥ 3):        {n_uncovered}")
    print(f"  Differenz (rk=0 aber |tors|≠4):     {n_total - n_torsion_trick - n_chabauty_only - n_uncovered}")

    print(f"\nVerteilung nach (rk_PQ, rk_uV, rk_3):")
    for pat, c in sorted(by_pattern.items()):
        marker = ""
        if 0 in pat:
            marker += " [TRICK]"
        elif sum(pat) < 3:
            marker += " [CHAB]"
        else:
            marker += " [HARD]"
        print(f"  {pat}: {c}{marker}")

    if uncovered:
        print(f"\nUnabgedeckte Fasern (zuerst 20):")
        for m, n, pat, lo, hi in uncovered[:20]:
            print(f"  ({m},{n}): rk = {pat}, total = {lo}-{hi}")


if __name__ == "__main__":
    main()
