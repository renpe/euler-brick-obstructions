"""
Versuch, den Rang von J(H_{m,n}) für Stichprobe (m,n) zu schätzen.

Strategie:
  - H ist Genus 3, daher J(H) hat dim 3.
  - H → E_PQ via s = t² (Genus-1-Quotient): Y² = P(s)·Q(s).
  - Damit ist E_PQ ein Faktor: rk(E_PQ) ≤ rk(J(H)).
  - Wir kennen bereits rk(E_{m,n}), rk(E'_{m,n}), rk(E_PQ).
  - Falls die alle bekannt und summe < 3 → Chabauty ggf. anwendbar.

Hinweis: J(H) zerlegt sich nicht unbedingt vollständig in elliptische
Faktoren. Es gibt einen Prym-Anteil von Dim 2.

Aufruf: sage cuboid_jacobian_rank.sage [M_MAX]
"""
from sage.all import *
import sys
import signal
from math import gcd as pygcd
from time import time

pari.allocatemem(int(2e9))

M_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
TIMEOUT = int(60)


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


def main():
    print(f"M_MAX={M_MAX}\n")
    print(f"{'(m,n)':>10} {'rk_E':>5} {'rk_Ep':>6} {'rk_EPQ':>7} {'lower_bnd_J(H)':>15} {'time':>6}")

    rank_dist = []
    for m in range(2, M_MAX + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1: continue
            if (m - n) % 2 != 1: continue
            t0 = time()
            U2 = m*m - n*n
            V2 = 2*m*n
            W2 = m*m + n*n

            Rt = PolynomialRing(QQ, 'T')
            T = Rt.gen()
            P_t = V2**2 * T**4 + (4*U2**2 - 2*V2**2) * T**2 + V2**2
            Q_t = W2**2 * T**4 + 2*(U2**2 - V2**2) * T**2 + W2**2

            # E_PQ via s=t²: P(s)*Q(s) als Quartik in s
            Rs = PolynomialRing(QQ, 'S')
            S = Rs.gen()
            P_s = V2**2 * S**2 + (4*U2**2 - 2*V2**2) * S + V2**2
            Q_s = W2**2 * S**2 + 2*(U2**2 - V2**2) * S + W2**2
            PQ = P_s * Q_s   # Grad 4 in s

            try:
                E   = quartic_to_E(P_t)
                Ep  = quartic_to_E(Q_t)
                Epq = quartic_to_E(PQ)  # Genus-1-Quotient
            except Exception as ex:
                print(f"{f'({m},{n})':>10}  fail: {ex}")
                continue
            rE = safe_rank(E)
            rEp = safe_rank(Ep)
            rPQ = safe_rank(Epq)
            elapsed = time() - t0

            def fmt(rk):
                if rk[0] is None: return "?"
                if rk[0] == rk[1]: return f"{rk[0]}"
                return f"{rk[0]}-{rk[1]}"

            # Untere Schranke für rk(J(H)): mindestens rk(E_PQ).
            # Aber WICHTIG: rk(E) und rk(E') kommen aus der Volle C, nicht aus H.
            # Wir vermerken sie nur informativ.
            lb = fmt(rPQ)
            print(f"{f'({m},{n})':>10} {fmt(rE):>5} {fmt(rEp):>6} {fmt(rPQ):>7} {lb:>15} {elapsed:>6.1f}")
            sys.stdout.flush()
            rank_dist.append((m, n, rE, rEp, rPQ))

    print(f"\n========== Verteilung von rk(E_PQ) ==========")
    from collections import Counter
    counts = Counter()
    for _, _, _, _, rPQ in rank_dist:
        if rPQ[0] is not None:
            counts[rPQ[0]] += 1
        else:
            counts['?'] += 1
    for r in sorted(counts.keys(), key=lambda x: (isinstance(x, str), x)):
        print(f"  rk(E_PQ) = {r}: {counts[r]} Fasern")


if __name__ == "__main__":
    main()
