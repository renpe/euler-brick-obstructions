"""
Versuch, rk(J(H_{m,n})) für eine Auswahl Fasern direkt zu bestimmen.
Sage-Methoden: analytic_rank, rank_bound, etc.
"""
from sage.all import *
import sys
import signal

pari.allocatemem(int(2e9))


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(int(0))


def main():
    cases = [(2, 1), (3, 2), (5, 2), (7, 4), (8, 1), (11, 4)]
    for m, n in cases:
        U2 = m*m - n*n
        V2 = 2*m*n
        W2 = m*m + n*n
        Rt = PolynomialRing(QQ, 'T')
        T = Rt.gen()
        P = V2**2 * T**4 + (4*U2**2 - 2*V2**2) * T**2 + V2**2
        Q = W2**2 * T**4 + 2*(U2**2 - V2**2) * T**2 + W2**2
        f = P * Q
        H = HyperellipticCurve(f, 0)
        J = H.jacobian()

        print(f"=== (m,n)=({m},{n}) ===")
        print(f"  H: y² = {f}")
        print(f"  Genus: {H.genus()}")
        # Check available methods on J
        methods = [m for m in dir(J) if not m.startswith('_') and 'rank' in m.lower()]
        print(f"  Methoden mit 'rank': {methods}")

        # Try various rank computations
        for method_name in ['rank', 'analytic_rank', 'rank_bound', 'mwrank']:
            method = getattr(J, method_name, None)
            if method is None:
                continue
            try:
                result = with_timeout(60, lambda: method())
                print(f"  J.{method_name}() = {result}")
            except Exception as ex:
                print(f"  J.{method_name}() failed: {ex}")

        # Try absolute decomposition via PARI
        try:
            print(f"  L-Polynom (modulo p, kleiner Test):")
            for p in [3, 5, 7, 11]:
                try:
                    Hp = HyperellipticCurve(f.change_ring(GF(p)), 0)
                    n_pts = Hp.count_points()
                    print(f"    #H(F_{p}) = {n_pts}")
                except Exception:
                    pass
        except Exception:
            pass
        print()


if __name__ == "__main__":
    main()
