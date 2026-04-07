"""
Product analysis: Is f1 * f2 ever a perfect square?

Paper claim supported
--------------------
"For all 1 <= a < b <= 1000 and 1 <= m < n <= 1000 with
 gcd(a,b) = gcd(m,n) = 1, a-b and m-n odd:
 f1 * f2 is never a perfect square."

What this script does
---------------------
Given the Euler brick parametrisation

    A = 2*(a^2 - b^2)*m*n
    B = 4*a*b*m*n
    C = (a^2 + b^2)*(m^2 - n^2)

it defines

    f1 = A^2 + C^2
    f2 = B^2 + C^2

and exhaustively checks whether the product f1*f2 is ever a perfect
square across all coprime pairs (a,b) and (m,n) within a given bound.

By the Brahmagupta-Fibonacci identity the product decomposes as

    f1*f2 = (A*B - C^2)^2 + C^2*(A + B)^2

Why it matters
--------------
If f1*f2 is never a perfect square, then f1 and f2 cannot both be
perfect squares simultaneously.  Since a perfect Euler brick requires
both f1 and f2 to be perfect squares (they correspond to the face
diagonals), this would rule out its existence for the tested parameter
range -- supporting the structural proof that H_total is never a
perfect square.

Usage
-----
    python produkt_analyse.py [bound] [workers]

    bound   -- upper limit for parameters a, b, m, n (default 80)
    workers -- number of parallel processes        (default 8)
"""

from math import gcd, isqrt
from concurrent.futures import ProcessPoolExecutor
import sys


def ist_quadrat(n):
    """Return True if n is a perfect square, False otherwise."""
    if n <= 0:
        return n == 0
    r = isqrt(n)
    return r * r == n


def suche_a(args):
    """Search all valid (b, m, n) for a fixed value of a.

    For each coprime pair (a, b) with b > a and each coprime pair
    (m, n) with n > m, compute A, B, C from the parametrisation,
    then check whether f1*f2 is a perfect square.

    Returns a list of hits: (a, b, m, n, |A|, B, |C|,
                             f1_is_square, f2_is_square, product).
    """
    a, bound = args
    treffer = []  # "treffer" = hits / matches
    for b in range(a + 1, bound + 1):
        if gcd(a, b) != 1:
            continue

        # Precompute quantities that depend only on (a, b)
        S = a * a + b * b          # a^2 + b^2
        diff = a * a - b * b       # a^2 - b^2  (negative since b > a)

        for m in range(1, bound):
            for n in range(m + 1, bound + 1):
                if gcd(m, n) != 1:
                    continue

                # Euler brick parametrisation
                A = 2 * diff * m * n           # = 2*(a^2 - b^2)*m*n
                B = 4 * a * b * m * n          # = 4*a*b*m*n
                C = S * (m * m - n * n)        # = (a^2+b^2)*(m^2-n^2)

                # The two factors whose product we test
                f1 = A * A + C * C             # A^2 + C^2
                f2 = B * B + C * C             # B^2 + C^2
                produkt = f1 * f2              # the product to test

                # Check each factor and the product for perfect-squareness
                sq_prod = ist_quadrat(produkt)
                sq_f1 = ist_quadrat(f1)
                sq_f2 = ist_quadrat(f2)

                if sq_prod:
                    treffer.append((a, b, m, n, abs(A), B, abs(C),
                                    sq_f1, sq_f2, produkt))

    return treffer


def main():
    # Parse optional command-line arguments: bound and worker count
    bound = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8

    print(f"Searching f1*f2 = perfect square: b>a>0, n>m>0, gcd=1, bound={bound}")
    print()

    # Build one task per value of a; each task searches all (b, m, n)
    args_list = [(a, bound) for a in range(1, bound)]

    all_hits = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for i, result in enumerate(executor.map(suche_a, args_list)):
            all_hits.extend(result)
            # Progress report every 20 values of a
            if (i + 1) % 20 == 0:
                print(f"  a={i+1}/{bound-1}, {len(all_hits)} product hits so far")

    all_hits.sort()

    # --- Tally up the results ---
    n_prod = len(all_hits)
    n_f1 = sum(1 for t in all_hits if t[7])
    n_f2 = sum(1 for t in all_hits if t[8])
    n_both = sum(1 for t in all_hits if t[7] and t[8])
    n_neither = sum(1 for t in all_hits if not t[7] and not t[8])

    # --- Print summary ---
    print(f"\n{'='*70}")
    print(f"RESULT (bound={bound})")
    print(f"{'='*70}")
    print(f"f1*f2 = perfect square:       {n_prod}")
    print(f"  of which f1 = square:       {n_f1}")
    print(f"  of which f2 = square:       {n_f2}")
    print(f"  of which BOTH = square:     {n_both}")
    print(f"  of which NEITHER = square:  {n_neither}")
    print()

    if n_prod == 0:
        print("*** f1*f2 is NEVER a perfect square! ***")
        print("This is consistent with non-existence of perfect Euler bricks.")
    else:
        print(f"f1*f2 is a perfect square in {n_prod} cases.")
        print()

        # Table header for detailed output
        print(f"{'a':>3} {'b':>3} {'m':>3} {'n':>3} | {'A':>10} {'B':>10} {'C':>10} | "
              f"{'f1=□':>5} {'f2=□':>5} | sqrt(f1*f2)")
        print("-" * 80)
        for a, b, m, n, A, B, C, sq1, sq2, prod in all_hits[:50]:
            r = isqrt(prod)
            f1_str = "yes" if sq1 else "no"
            f2_str = "yes" if sq2 else "no"
            print(f"{a:>3} {b:>3} {m:>3} {n:>3} | {A:>10} {B:>10} {C:>10} | "
                  f"{f1_str:>5} {f2_str:>5} | {r}")

        # Extra section: cases where f1*f2 is square but neither f1 nor f2 is
        if n_neither > 0:
            print(f"\n--- Hits where f1*f2 = square but NEITHER factor alone: ---")
            for a, b, m, n, A, B, C, sq1, sq2, prod in all_hits:
                if not sq1 and not sq2:
                    r = isqrt(prod)
                    print(f"  ({a},{b},{m},{n}): A={A}, B={B}, C={C}, "
                          f"sqrt(f1*f2)={r}")


if __name__ == "__main__":
    main()
