"""
Full rank determination of J(H_{m,n}) by decomposition into 3 elliptic
factors.

H_{m,n}: v^2 = P(t^2) * Q(t^2) has automorphisms:
  sigma_1: t -> -t  (even quartic in t)
  sigma_2: t -> 1/t (P, Q are palindromes in t)

Quotients:
  E_PQ via s=t^2:        Y^2 = P(s) * Q(s)                     (degree 4 in s)
  E_uV via u=t+1/t:      V^2 = (V2^2 u^2 + 4(U2^2 - V2^2)) * (W2^2 u^2 - 4 V2^2)
  E_3   via composition sigma_1 sigma_2: t -> -1/t

J(H) ~ E_PQ x E_uV x E_3

We compute rk of all three via PARI ellrank, sum.

Usage: sage jh_full_decomp.sage [M_MAX]
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

            # Factor 1: E_PQ : Y^2 = P(s)Q(s) (quartic in s)
            Rs = PolynomialRing(QQ, 'S')
            S = Rs.gen()
            P_s = V2**2 * S**2 + (4*U2**2 - 2*V2**2) * S + V2**2
            Q_s = W2**2 * S**2 + 2*(U2**2 - V2**2) * S + W2**2
            PQ = P_s * Q_s   # degree 4 in s

            # Factor 2: E_uV : V^2 = (V2^2 u^2 + 4(U2^2 - V2^2)) * (W2^2 u^2 - 4 V2^2)
            Ru = PolynomialRing(QQ, 'U')
            U = Ru.gen()
            uV_quartic = (V2**2 * U**2 + 4*(U2**2 - V2**2)) * (W2**2 * U**2 - 4*V2**2)

            # Factor 3: from sigma_1 sigma_2 (t -> -1/t). Setting w = t - 1/t,
            # then t^2 + 1/t^2 = w^2 + 2, t^4 + 1 = t^2 * (t^2 + 1/t^2) = t^2 * (w^2+2).
            # P(t^2)/t^2 = V2^2 (t^2 + 1/t^2) + (4 U2^2 - 2 V2^2) = V2^2 (w^2 + 2) + 4 U2^2 - 2 V2^2
            #            = V2^2 w^2 + 4 U2^2
            # Q(t^2)/t^2 = W2^2 (w^2 + 2) + 2 (U2^2 - V2^2) = W2^2 w^2 + 2 (U2^2 - V2^2 + W2^2)
            #            = W2^2 w^2 + 2 (2 U2^2) = W2^2 w^2 + 4 U2^2
            # So: P(t^2) Q(t^2) / t^4 = (V2^2 w^2 + 4 U2^2) (W2^2 w^2 + 4 U2^2)
            # E_3 : Y^2 = (V2^2 w^2 + 4 U2^2) (W2^2 w^2 + 4 U2^2)
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

            # Sum bound
            if all(r[0] is not None for r in [rk_PQ, rk_uV, rk_3]):
                lo = rk_PQ[0] + rk_uV[0] + rk_3[0]
                hi = rk_PQ[1] + rk_uV[1] + rk_3[1]
                jh_str = f"{lo}-{hi}" if lo != hi else str(lo)
                # Chabauty: rk(J) < genus(H) = 3
                chab = "yes" if hi < 3 else "no"
                if hi < 3:
                    chab_count += 1
            else:
                jh_str = "?"
                chab = "?"

            print(f"{f'({m},{n})':>10} {fmt(rk_PQ):>6} {fmt(rk_uV):>6} {fmt(rk_3):>6} "
                  f"{jh_str:>8} {'3':>2} {chab:>5}")
            sys.stdout.flush()

    print(f"\n========== Summary ==========")
    print(f"Total fibers checked: {n_total}")
    print(f"Chabauty applicable (rk(J)<3): {chab_count} -> Conjecture B explicitly provable")


if __name__ == "__main__":
    main()
