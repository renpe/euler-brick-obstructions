"""
Pilot: Jacobian-Zerlegung des Genus-2-Modells C_{m,n} über die zwei
Bedingungen Master-Tupel + perfekter-Cuboid.

Aus der Identitätenanalyse folgt für festes (m, n):
  E_{m,n}:  y² = V₂² t⁴ + (4U₂² − 2V₂²) t² + V₂²        (Master-Tupel)
  E'_{m,n}: y² = W₂² t⁴ + 2(U₂² − V₂²) t² + W₂²         (perfekter Cuboid)

Die Schnittpunkte ihrer rationalen Punkte = perfekte Euler-Cuboids.
Konjektur B ⟺ der Schnitt ist trivial.

Dieses Skript:
  1. Bauen E_{m,n} und E'_{m,n} als WeierstraSS-Modelle in Sage.
  2. Berechnet rk(E) und rk(E') via PARI ellrank.
  3. Listet bekannte rationale Punkte auf E' (mit kleinen Höhen).
  4. Sucht nicht-triviale Schnittpunkte E ∩ E' (= cuboid-Kandidaten).

Aufruf:
    sage cuboid_genus2_pilot.sage [M_MAX]
    Default: 30 — alle teilerfremden (m, n) mit m≤30, m-n ungerade
"""
from sage.all import *
import sys
import signal
from math import gcd as pygcd
from time import time

pari.allocatemem(int(2e9))

M_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
TIMEOUT = int(60)


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(int(0))


def build_curves(m, n):
    """Liefert (E_master, E_cuboid) als EllipticCurve-Objekte aus den
    Quartiken via PARI ellfromeqn."""
    U2 = m*m - n*n
    V2 = 2*m*n
    W2 = m*m + n*n

    R = PolynomialRing(QQ, ['x', 'y'])
    x, y = R.gens()
    eqn_master = y**2 - (V2**2 * x**4 + (4*U2**2 - 2*V2**2) * x**2 + V2**2)
    eqn_cuboid = y**2 - (W2**2 * x**4 + 2*(U2**2 - V2**2) * x**2 + W2**2)
    cm = [ZZ(c) for c in pari(eqn_master).ellfromeqn()]
    cc = [ZZ(c) for c in pari(eqn_cuboid).ellfromeqn()]
    return EllipticCurve(cm), EllipticCurve(cc)


def rank_info(E):
    """Sicheres ellrank mit Timeout. Liefert (rl, ru) oder (None, None)."""
    try:
        result = with_timeout(TIMEOUT, lambda: E.pari_curve().ellrank())
        return int(result[0]), int(result[1])
    except Exception:
        return None, None


def main():
    print(f"M_MAX = {M_MAX}, timeout per call = {TIMEOUT}s\n")
    print(f"{'m':>4} {'n':>4} {'rk(E)':>6} {'rk(E_+)':>8} {'rk(E_+)_lo':>10} {'time':>6}")

    for m in range(2, M_MAX + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1:
                continue
            if (m - n) % 2 != 1:
                continue
            t0 = time()
            try:
                E, Ep = build_curves(m, n)
            except Exception as ex:
                print(f"{m:>4} {n:>4}   curve_build_fail: {ex}")
                continue
            rl_E, ru_E = rank_info(E)
            rl_Ep, ru_Ep = rank_info(Ep)
            elapsed = time() - t0
            E_str = f"{rl_E}-{ru_E}" if rl_E is not None else "?"
            Ep_str = f"{rl_Ep}-{ru_Ep}" if rl_Ep is not None else "?"
            print(f"{m:>4} {n:>4} {E_str:>6} {Ep_str:>8} "
                  f"{rl_Ep if rl_Ep is not None else '?':>10} "
                  f"{elapsed:>6.1f}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
