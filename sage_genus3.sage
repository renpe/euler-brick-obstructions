"""
sage_genus3.sage -- Genus-3 curve analysis for the Euler brick problem

Supports: Section 4 (genus-3 curve, Jacobian decomposition, quotients)

The complete condition for a perfect Euler brick is:

  w^2 = nu^8 + A*nu^4 + 1     (genus-3 curve!)

with A = 2 - 4*c^2, c = (s^4-6s^2+1)/(1+s^2)^2, s = a/b.

This is the ORIGINAL quartic w^2 = mu^4 + A*mu^2 + 1 with mu = nu^2.
The genus-3 curve arises by "unfolding" the square condition.

By Faltings' theorem every curve of genus >= 2 over Q has only finitely
many rational points.  For genus 3, Chabauty-Coleman can be applied when
the rank of the Jacobian is < 3.

Obvious points: (nu, w) = (0, +-1) and nu = +-1 if A = 0 (degenerate).

For each coprime pair (a, b) the script:
  1. builds the hyperelliptic curve w^2 = nu^8 + A*nu^4 + 1,
  2. searches for non-trivial rational points (integers and small fractions),
  3. back-substitutes any hit to check whether f1 and f2 are squares,
  4. prints a verification section explaining the link between the
     genus-3 curve and the original Euler brick problem.
"""

from sage.all import *

test_pairs = [
    (1, 2), (1, 3), (1, 4), (1, 5), (1, 7),
    (2, 3), (2, 5), (2, 7),
    (3, 4), (3, 5), (3, 7),
    (4, 5), (4, 7),
    (5, 7), (5, 8), (5, 9),
    (7, 9), (7, 11),
    (8, 9),
]

print("=" * 80)
print("GENUS-3 CURVE: w^2 = nu^8 + A*nu^4 + 1")
print("=" * 80)
print()

for a_val, b_val in test_pairs:
    if gcd(a_val, b_val) != 1:
        continue
    s_val = QQ(a_val) / QQ(b_val)
    c_val = (s_val^4 - 6*s_val^2 + 1) / (1 + s_val^2)^2
    A_val = 2 - 4*c_val^2

    # Construct the genus-3 curve w^2 = f(nu)
    R.<nu> = QQ[]
    f = nu^8 + A_val*nu^4 + 1

    # Scale to obtain integer coefficients
    # A_val is rational, so multiply f by an appropriate power of the denominator
    den = denominator(A_val)
    # w^2 = nu^8 + (p/q)*nu^4 + 1
    # Multiply by q^4: (q^2*w)^2 = q^4*nu^8 + p*q^3*nu^4 + q^4
    # Set W = q^2*w, N = q*nu: this gets complicated.
    # Simpler: work over QQ

    try:
        H = HyperellipticCurve(f)
        g = H.genus()

        print("s=%s (a=%d,b=%d): A=%s" % (s_val, a_val, b_val, A_val))
        print("  Genus = %d" % g)

        # Obvious points
        # nu=0: w^2=1, so (0,1) and (0,-1)
        # Search for additional small rational points
        pts = []
        for nu_test in range(-20, 21):
            nu_r = QQ(nu_test)
            w2 = nu_r^8 + A_val*nu_r^4 + 1
            if w2 >= 0 and is_square(numerator(w2)) and is_square(denominator(w2)):
                w_val = sqrt(numerator(w2)) / sqrt(denominator(w2))
                if nu_test != 0:
                    pts.append((nu_r, w_val))

        # Also test fractions
        for nu_num in range(1, 10):
            for nu_den in range(1, 10):
                if gcd(nu_num, nu_den) != 1:
                    continue
                for sign in [1, -1]:
                    nu_r = QQ(sign*nu_num) / QQ(nu_den)
                    w2 = nu_r^8 + A_val*nu_r^4 + 1
                    if w2 >= 0 and is_square(numerator(w2)) and is_square(denominator(w2)):
                        w_val = sqrt(numerator(w2)) / sqrt(denominator(w2))
                        pts.append((nu_r, w_val))

        if pts:
            print("  Non-trivial rational points found:")
            for nu_r, w_val in pts:
                # Back-substitution: mu = nu^2, lambda = nu
                # NOTE: lambda^2 = mu = nu^2, so lambda = nu
                # The Euler brick parameter would then be lambda = nu, s = a/b
                lam = nu_r  # lambda = nu (since mu = nu^2 = lambda^2)
                m_over_n = lam  # lambda = m/n
                print("    (nu, w) = (%s, %s) -> lambda = %s = m/n" % (nu_r, w_val, lam))

                # Verify: are f1 and f2 squares?
                t = s_val  # = a/b
                u = lam    # = m/n (= lambda)
                f1_norm = u^4 + 2*c_val*u^2 + 1
                f2_norm = u^4 - 2*c_val*u^2 + 1
                print("    f1_norm = %s (square: %s)" % (f1_norm, is_square(numerator(f1_norm)) and is_square(denominator(f1_norm)) if f1_norm > 0 else False))
                print("    f2_norm = %s (square: %s)" % (f2_norm, is_square(numerator(f2_norm)) and is_square(denominator(f2_norm)) if f2_norm > 0 else False))
        else:
            print("  No non-trivial rational points for nu in [-20..20] and fractions up to 9/9")

        # Jacobian rank (if possible)
        # For genus 3 this is difficult, but let us try
        # Sage can estimate the Jacobian rank for some hyperelliptic curves

        print()

    except Exception as e:
        print("s=%s (a=%d,b=%d): ERROR: %s" % (s_val, a_val, b_val, str(e)[:60]))
        print()

# Verification: confirm that the genus-3 curve is the correct one
print("=" * 80)
print("VERIFICATION: Link between genus-3 curve and original problem")
print("=" * 80)
print()
print("The quartic w^2 = mu^4 + A*mu^2 + 1 has genus 1.")
print("The condition mu = lambda^2 (square) sets mu = nu^2.")
print("Substituting: w^2 = nu^8 + A*nu^4 + 1 has genus 3.")
print()
print("A rational point (nu, w) with nu != 0 gives:")
print("  lambda = nu  (since lambda^2 = mu = nu^2)")
print("  m/n = lambda = nu")
print()
print("Together with s = a/b this yields the complete parameter set.")
print()

# Test the OBVIOUS points (0, +-1)
print("Obvious points:")
print("  (nu, w) = (0, 1): lambda = 0 -> m = 0 (degenerate)")
print("  (nu, w) = (0,-1): lambda = 0 -> m = 0 (degenerate)")
print()

# Are there points at infinity?
print("Points at nu = infinity:")
print("  w^2 ~ nu^8 -> w ~ nu^4, so (1 : nu^4 : 0) in projective space")
print("  But nu = infinity -> lambda = infinity -> n = 0 (degenerate)")
