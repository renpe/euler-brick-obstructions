"""
Saunderson Master-Hit generator.

Implements Theorem 5.1 of Peschmann, "Computational evidence and
structural obstructions for the perfect cuboid problem" (2026):
for coprime g > h > 0 with g - h odd,

    a = 4gh
    b = g^2 + h^2
    m = h(3g^2 - h^2)
    n = g(g^2 - 3h^2)     (with m, |n| swapped if n < 0)

gives a Master-Hit (a, b, m, n).

Usage:
    python saunderson_generator.py [max_generator]
"""

from __future__ import annotations

import sys
from math import gcd, isqrt


def saunderson_master_hit(g: int, h: int) -> tuple[int, int, int, int] | None:
    """Return the Master-Hit (a, b, m, n) produced by Saunderson's
    formula for generators (g, h), or None if the parameters are
    invalid (not coprime, wrong parity, or produce degenerate output).
    """
    if g <= h <= 0:
        return None
    if gcd(g, h) != 1:
        return None
    if (g - h) % 2 == 0:
        return None

    a = 4 * g * h
    b = g * g + h * h
    m = h * (3 * g * g - h * h)
    n = g * (g * g - 3 * h * h)

    # When n is negative (for small g/h with h relatively large),
    # swap m and |n| to obtain the standard Euclid form m > n > 0.
    if n < 0:
        m, n = n * -1, m
        m, n = max(m, n), min(m, n)
        if m <= n:
            return None
    if m <= n or n <= 0:
        return None
    if gcd(m, n) != 1 or (m - n) % 2 == 0:
        return None

    return (a, b, m, n)


def verify_master_condition(a: int, b: int, m: int, n: int) -> bool:
    """Return True iff the Master condition (V1*U2)^2 + (U1*V2)^2 = square
    holds for (a, b, m, n), where U1 = a^2 - b^2, V1 = 2ab, etc.
    """
    u1 = a * a - b * b
    v1 = 2 * a * b
    u2 = m * m - n * n
    v2 = 2 * m * n
    master = (v1 * u2) ** 2 + (u1 * v2) ** 2
    root = isqrt(master)
    return root * root == master


def enumerate_saunderson_hits(max_generator: int):
    """Yield all Master-Hits produced by Saunderson generators (g, h)
    with g <= max_generator.
    """
    for g in range(2, max_generator + 1):
        for h in range(1, g):
            hit = saunderson_master_hit(g, h)
            if hit is None:
                continue
            yield (g, h), hit


def main() -> int:
    max_generator = int(sys.argv[1]) if len(sys.argv) > 1 else 20

    print(f"Enumerating Saunderson Master-Hits for g <= {max_generator}")
    print()
    print(f"  {'g':>3s} {'h':>3s} | {'a':>6s} {'b':>6s} {'m':>7s} {'n':>7s}  verified")
    print(f"  {'-' * 3} {'-' * 3} | {'-' * 6} {'-' * 6} {'-' * 7} {'-' * 7}  --------")

    count = 0
    for (g, h), (a, b, m, n) in enumerate_saunderson_hits(max_generator):
        ok = verify_master_condition(a, b, m, n)
        flag = "yes" if ok else "FAIL"
        print(f"  {g:3d} {h:3d} | {a:6d} {b:6d} {m:7d} {n:7d}  {flag}")
        count += 1

    print()
    print(f"Total: {count} Master-Hits generated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
