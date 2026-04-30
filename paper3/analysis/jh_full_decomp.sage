"""
Vollständige Rang-Bestimmung von J(H_{m,n}) durch Zerlegung in 3 elliptische
Faktoren.

H_{m,n}: v² = P(t²)·Q(t²) hat Automorphismen:
  σ₁: t → −t  (gerade Quartik in t)
  σ₂: t → 1/t (P, Q sind Palindrome in t)

Quotienten:
  E_PQ via s=t²:        Y² = P(s)·Q(s)                     (Grad 4 in s)
  E_uV via u=t+1/t:     V² = (V₂²u² + 4(U₂²−V₂²)) · (W₂²u² − 4V₂²)
  E_3   via Komposition σ₁σ₂: t → −1/t

J(H) ~ E_PQ × E_uV × E_3

Wir berechnen rk aller drei via PARI ellrank, summieren.

Aufruf: sage jh_full_decomp.sage [M_MAX]
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


def fmt(rk):
    if rk[0] is None: return "?"
    if rk[0] == rk[1]: return str(rk[0])
    return f"{rk[0]}-{rk[1]}"


def main():
    print(f"M_MAX={M_MAX}\n")
    print(f"{'(m,n)':>10} {'rk_PQ':>6} {'rk_uV':>6} {'rk_3':>6} {'rk_J(H)':>8} {'g':>2} "
          f"{'Chab?':>5}")

    chab_count = 0
    n_total = 0

    for m in range(2, M_MAX + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1: continue
            if (m - n) % 2 != 1: continue
            n_total += 1
            U2 = m*m - n*n
            V2 = 2*m*n
            W2 = m*m + n*n

            # Faktor 1: E_PQ : Y² = P(s)Q(s) (Quartik in s)
            Rs = PolynomialRing(QQ, 'S')
            S = Rs.gen()
            P_s = V2**2 * S**2 + (4*U2**2 - 2*V2**2) * S + V2**2
            Q_s = W2**2 * S**2 + 2*(U2**2 - V2**2) * S + W2**2
            PQ = P_s * Q_s   # Grad 4 in s

            # Faktor 2: E_uV : V² = (V₂²u² + 4(U₂²−V₂²)) · (W₂²u² − 4V₂²)
            Ru = PolynomialRing(QQ, 'U')
            U = Ru.gen()
            uV_quartic = (V2**2 * U**2 + 4*(U2**2 - V2**2)) * (W2**2 * U**2 - 4*V2**2)

            # Faktor 3: aus σ₁σ₂ (t → −1/t). Setzen wir w = t − 1/t,
            # dann t² + 1/t² = w² + 2, t⁴ + 1 = t²·(t² + 1/t²) = t²·(w²+2).
            # P(t²)/t² = V₂²(t² + 1/t²) + (4U₂² − 2V₂²) = V₂²(w² + 2) + 4U₂² − 2V₂²
            #          = V₂²w² + 4U₂²
            # Q(t²)/t² = W₂²(w² + 2) + 2(U₂² − V₂²) = W₂²w² + 2(U₂² − V₂² + W₂²)
            #          = W₂²w² + 2(2U₂²) = W₂²w² + 4U₂²
            # Also: P(t²)Q(t²)/t⁴ = (V₂²w² + 4U₂²)(W₂²w² + 4U₂²)
            # E_3 : Y² = (V₂²w² + 4U₂²)(W₂²w² + 4U₂²)
            Rw = PolynomialRing(QQ, 'W')
            W = Rw.gen()
            third_quartic = (V2**2 * W**2 + 4*U2**2) * (W2**2 * W**2 + 4*U2**2)

            try:
                E_PQ  = quartic_to_E(PQ)
                E_uV  = quartic_to_E(uV_quartic)
                E_3   = quartic_to_E(third_quartic)
            except Exception as ex:
                print(f"{f'({m},{n})':>10} fail: {ex}")
                continue

            rk_PQ = safe_rank(E_PQ)
            rk_uV = safe_rank(E_uV)
            rk_3  = safe_rank(E_3)

            # Summe-Bound
            if all(r[0] is not None for r in [rk_PQ, rk_uV, rk_3]):
                lo = rk_PQ[0] + rk_uV[0] + rk_3[0]
                hi = rk_PQ[1] + rk_uV[1] + rk_3[1]
                jh_str = f"{lo}-{hi}" if lo != hi else str(lo)
                # Chabauty: rk(J) < genus(H) = 3
                chab = "✓" if hi < 3 else "✗"
                if hi < 3:
                    chab_count += 1
            else:
                jh_str = "?"
                chab = "?"

            print(f"{f'({m},{n})':>10} {fmt(rk_PQ):>6} {fmt(rk_uV):>6} {fmt(rk_3):>6} "
                  f"{jh_str:>8} {'3':>2} {chab:>5}")
            sys.stdout.flush()

    print(f"\n========== Zusammenfassung ==========")
    print(f"Gesamt geprüfte Fasern: {n_total}")
    print(f"Chabauty anwendbar (rk(J)<3): {chab_count} → Konjektur B explizit beweisbar")


if __name__ == "__main__":
    main()
