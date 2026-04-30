"""
Verification: what is the ACTUAL genus of the cuboid curve?

We had claimed: J(C_{m,n}) = E_{m,n} x E'_{m,n}, hence genus 2.
On a more careful check (Riemann-Hurwitz), everything points to
genus 5 (cover of a genus-1 via t^2 = s, with 8 ramification
points instead of the 2 I had assumed).

We set up the hyperelliptic form directly:
  v^2 = P(t^2) * Q(t^2)  - where v = y * w
P, Q are quadratics in s = t^2, so P*Q is degree 8 in t.
This is a hyperelliptic model GR of the curve.
"""
from sage.all import *
import sys

pari.allocatemem(int(2e9))


def main():
    m_n_pairs = [(5, 2), (8, 3), (9, 8), (11, 4)]
    for m, n in m_n_pairs:
        U2 = m*m - n*n
        V2 = 2*m*n
        W2 = m*m + n*n
        Rt = PolynomialRing(QQ, 'T')
        T = Rt.gen()
        P = V2**2 * T**4 + (4*U2**2 - 2*V2**2) * T**2 + V2**2
        Q = W2**2 * T**4 + 2*(U2**2 - V2**2) * T**2 + W2**2
        product = P * Q
        # v^2 = product, hyperelliptic curve
        # Calc genus via Hyperelliptic
        H = HyperellipticCurve(product, 0)
        g = H.genus()
        # Also: compute jacobian dimension if available
        print(f"(m,n)=({m},{n}):")
        print(f"  P(t) = {P}")
        print(f"  Q(t) = {Q}")
        print(f"  P*Q (degree {product.degree()}, squarefree? {product.is_squarefree()})")
        print(f"  Hyperelliptic curve v^2 = P(t)*Q(t) has Genus = {g}")
        print(f"  -> Jacobian dim = {g}")
        # Decomposition test: does the Jacobian decompose?
        try:
            J = H.jacobian()
            print(f"  Jacobian: {J}")
        except Exception:
            pass
        print()


if __name__ == "__main__":
    main()
