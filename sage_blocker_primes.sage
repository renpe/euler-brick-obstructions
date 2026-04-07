"""
Blocker-prime analysis for the product f1 * f2.

For every coprime pair (a, b) with a < b and every coprime pair (m, n) with
m < n, we form the four face-diagonal factors

    F1- = S*T - 4*a*b*m*n,     F1+ = S*T + 4*a*b*m*n,
    F2- = a^2*(m-n)^2 + b^2*(m+n)^2,
    F2+ = a^2*(m+n)^2 + b^2*(m-n)^2,

where S = a^2 + b^2 and T = m^2 + n^2, and set

    prod = F1- * F1+ * F2- * F2+.

This product is NEVER a perfect square.  Therefore there always exists a
prime p whose p-adic valuation in the product is odd.  We call the smallest
such prime the "blocker prime" for that tuple (a, b, m, n).

The script collects blocker primes over a search range and tabulates them
by residue class mod 4.

Supports the claim in Section 7, item 3 of the paper:

    "The smallest blocker satisfies p = 1 (mod 4) in 88.4 % of cases
     and p = 2 in 11.6 %; no prime p = 3 (mod 4) ever appears."

The factorisation in Z[i] (Gaussian integers) explains why only primes
that split in Z[i] (i.e. p = 2 or p = 1 mod 4) can occur as blockers.
"""

from math import gcd as mgcd, isqrt
from collections import Counter

K.<i> = QuadraticField(-1)

def blocker_prime(N):
    """Return the smallest prime p with v_p(N) odd, or None if N is a perfect square."""
    if N == 0:
        return None
    N = abs(int(N))
    if N == 1:
        return None
    f = factor(N)
    for p, e in f:
        if e % 2 == 1:
            return int(p)
    return None  # N is a perfect square

bound = 40
blocker_counts = Counter()
blocker_by_class = Counter()
total = 0
squares = 0

print("Searching blocker primes for f1*f2, bound=%d" % bound)
print()

details = []

for a in range(1, bound):
    for b in range(a+1, bound):
        if mgcd(a, b) != 1:
            continue
        S = a^2 + b^2
        for m in range(1, bound):
            for n in range(m+1, bound):
                if mgcd(m, n) != 1:
                    continue

                T = m^2 + n^2

                F1m = S*T - 4*a*b*m*n
                F1p = S*T + 4*a*b*m*n
                F2m = a^2*(m-n)^2 + b^2*(m+n)^2
                F2p = a^2*(m+n)^2 + b^2*(m-n)^2

                prod = F1m * F1p * F2m * F2p
                total += 1

                bp = blocker_prime(prod)
                if bp is None:
                    squares += 1
                    print("*** SQUARE at (%d,%d,%d,%d)! ***" % (a,b,m,n))
                else:
                    blocker_counts[bp] += 1
                    blocker_by_class[bp % 4] += 1
                    if len(details) < 30:
                        details.append((a, b, m, n, bp, int(prod)))

    if a % 10 == 0:
        print("  a=%d/%d, %d tested, %d squares" % (a, bound-1, total, squares))

print()
print("=" * 60)
print("RESULT (bound=%d)" % bound)
print("=" * 60)
print("Total tested: %d" % total)
print("Of which square: %d" % squares)
print()

print("Most frequent blocker primes:")
for p, cnt in blocker_counts.most_common(20):
    pct = 100.0 * cnt / total
    print("  p=%d (%d mod 4): %d (%.1f%%)" % (p, p % 4, cnt, pct))

print()
print("Blockers by residue class mod 4:")
for r in sorted(blocker_by_class.keys()):
    cnt = blocker_by_class[r]
    pct = 100.0 * cnt / total
    print("  p = %d mod 4: %d (%.1f%%)" % (r, cnt, pct))

print()
print("First 30 examples:")
for a, b, m, n, bp, prod in details:
    # Factorise the blocker contribution across the four factors
    F1m = (a^2+b^2)*(m^2+n^2) - 4*a*b*m*n
    F1p = (a^2+b^2)*(m^2+n^2) + 4*a*b*m*n
    F2m = a^2*(m-n)^2 + b^2*(m+n)^2
    F2p = a^2*(m+n)^2 + b^2*(m-n)^2

    v1m = valuation(F1m, bp)
    v1p = valuation(F1p, bp)
    v2m = valuation(F2m, bp)
    v2p = valuation(F2p, bp)

    print("  (%d,%d,%d,%d): blocker p=%d, v_p(F1-)=%d, v_p(F1+)=%d, v_p(F2-)=%d, v_p(F2+)=%d, sum=%d" %
          (a, b, m, n, bp, v1m, v1p, v2m, v2p, v1m+v1p+v2m+v2p))

# Check: Is the blocker always a divisor of c_numerator?
print()
print("=" * 60)
print("IS THE BLOCKER A DIVISOR OF c_numerator?")
print("=" * 60)
print()

c_divides = 0
c_not_divides = 0
for a, b, m, n, bp, prod in details:
    c_z = a^4 - 6*a^2*b^2 + b^4
    if c_z % bp == 0:
        c_divides += 1
    else:
        c_not_divides += 1

print("  Blocker divides c_numerator: %d / %d" % (c_divides, len(details)))
print("  Blocker does NOT divide c_numerator: %d / %d" % (c_not_divides, len(details)))

# Exhaustive check
c_divides_total = 0
for a in range(1, min(bound, 20)):
    for b in range(a+1, min(bound, 20)):
        if mgcd(a, b) != 1:
            continue
        for m in range(1, min(bound, 20)):
            for n in range(m+1, min(bound, 20)):
                if mgcd(m, n) != 1:
                    continue
                S = a^2 + b^2
                T = m^2 + n^2
                prod = (S*T - 4*a*b*m*n) * (S*T + 4*a*b*m*n) * \
                       (a^2*(m-n)^2 + b^2*(m+n)^2) * (a^2*(m+n)^2 + b^2*(m-n)^2)
                bp = blocker_prime(prod)
                if bp is not None:
                    c_z = a^4 - 6*a^2*b^2 + b^4
                    d_z = m^4 - 6*m^2*n^2 + n^4  # analogous for (m,n)
                    if c_z % bp == 0:
                        c_divides_total += 1

print()
print("(Exhaustive check up to bound=20)")
