"""
Verifizierung: was ist der WIRKLICHE Genus der Cuboid-Kurve?

Wir hatten behauptet: J(C_{m,n}) = E_{m,n} × E'_{m,n}, also Genus 2.
Bei genauerem Nachrechnen (Riemann-Hurwitz) deutet aber alles auf
Genus 5 hin (Cover von einem Genus-1 via t² = s, mit 8 Verzweigungs-
punkten statt der von mir angenommenen 2).

Wir setzen die Hyperelliptische Form direkt auf:
  v² = P(t²) · Q(t²)  — wobei v = y·w
P, Q sind Quadratiken in s = t², also P·Q ist Grad 8 in t.
Das ist ein hyperelliptisches Modell GR der Kurve.
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
        # v² = product, hyperelliptic curve
        # Calc genus via Hyperelliptic
        H = HyperellipticCurve(product, 0)
        g = H.genus()
        # Also: compute jacobian dimension if available
        print(f"(m,n)=({m},{n}):")
        print(f"  P(t) = {P}")
        print(f"  Q(t) = {Q}")
        print(f"  P*Q (degree {product.degree()}, squarefree? {product.is_squarefree()})")
        print(f"  Hyperelliptic curve v² = P(t)·Q(t) has Genus = {g}")
        print(f"  → Jacobian dim = {g}")
        # Decomposition test: does the Jacobian decompose?
        try:
            J = H.jacobian()
            print(f"  Jacobian: {J}")
        except Exception:
            pass
        print()


if __name__ == "__main__":
    main()
