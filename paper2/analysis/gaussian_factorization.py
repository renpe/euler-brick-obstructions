"""
Gaussian-integer factorization analysis of f1 = N(W1*U2 + i*U1*V2).

Implements the three-case classification of Section 4 of Peschmann (2026):

  Case A: p divides gcd(W1*U2, U1*V2)
          -> both conjugate Gaussian primes pi, pi_bar divide z
          -> v_p(f1) = v_pi(z) + v_pi_bar(z) (even in all observed data)

  Case B: p divides exactly one of W1*U2, U1*V2
          -> one Gaussian prime divides z, the other does not
          -> v_p(f1) = k, where k may be even or odd (a blocker if odd)

  Case C: p divides neither W1*U2 nor U1*V2 (split prime, "from the sum")
          -> still one Gaussian conjugate divides z
          -> v_p(f1) = k; a blocker when k is odd

Usage:
    python gaussian_factorization.py a b m n
"""

from __future__ import annotations

import sys
from math import gcd, isqrt

try:
    from sympy import factorint
except ImportError:
    factorint = None


def gaussian_prime_decomposition(p: int) -> tuple[int, int] | None:
    """Return (a, b) with a^2 + b^2 = p and a > b >= 0, for p = 2 or p ≡ 1 (mod 4).
    Returns None for inert primes (p ≡ 3 (mod 4), p > 2).
    """
    if p == 2:
        return (1, 1)
    if p % 4 == 3:
        return None
    for a in range(1, isqrt(p) + 1):
        bsq = p - a * a
        if bsq <= 0:
            break
        b = isqrt(bsq)
        if b * b == bsq:
            if a >= b:
                return (a, b)
            return (b, a)
    return None


def gaussian_valuations(z_re: int, z_im: int,
                        pi_re: int, pi_im: int) -> int:
    """Return v_pi(z_re + i*z_im), where pi = pi_re + i*pi_im (a Gaussian prime).
    """
    norm = pi_re * pi_re + pi_im * pi_im
    count = 0
    re, im = z_re, z_im
    while True:
        # (re + i*im) / (pi_re + i*pi_im) = ((re + i*im)(pi_re - i*pi_im)) / norm
        new_re = re * pi_re + im * pi_im
        new_im = im * pi_re - re * pi_im
        if new_re % norm != 0 or new_im % norm != 0:
            break
        re = new_re // norm
        im = new_im // norm
        count += 1
    return count


def classify_prime(a: int, b: int, m: int, n: int, p: int) -> dict:
    """Classify the role of prime p in f1 for the Master-Hit (a, b, m, n).

    Returns a dict with keys:
        case:       'A', 'B', or 'C'
        exponent:   total v_p(f1)
        v_pi:       valuation at pi (0 if inert)
        v_pi_bar:   valuation at conjugate (0 if inert)
        is_blocker: exponent is odd
    """
    u1 = a * a - b * b
    v1 = 2 * a * b
    w1 = a * a + b * b
    u2 = m * m - n * n
    v2 = 2 * m * n

    z_re = w1 * u2
    z_im = u1 * v2

    divides_re = (z_re % p == 0)
    divides_im = (z_im % p == 0)

    if p % 4 == 3:
        # Inert prime: v_p(f1) is twice v_p on either summand when both
        # are divisible simultaneously, otherwise f1 ≡ sum-of-squares mod p.
        case = "inert"
    elif divides_re and divides_im:
        case = "A"
    elif divides_re != divides_im:
        case = "B"
    else:
        case = "C"

    if p == 2 or p % 4 != 3:
        gp = gaussian_prime_decomposition(p)
        if gp is None:
            v_pi = v_pi_bar = 0
            exponent = 0
        else:
            pi_re, pi_im = gp
            v_pi = gaussian_valuations(z_re, z_im, pi_re, pi_im)
            v_pi_bar = gaussian_valuations(z_re, z_im, pi_re, -pi_im)
            exponent = v_pi + v_pi_bar
    else:
        # Inert: treat p as (p, 0) with norm p^2.
        v_pi = v_pi_bar = 0
        exponent = 0  # computed via sympy fallback below

    return {
        "prime": p,
        "case": case,
        "v_pi": v_pi,
        "v_pi_bar": v_pi_bar,
        "exponent": exponent,
        "is_blocker": exponent % 2 == 1,
        "divides_z_re": divides_re,
        "divides_z_im": divides_im,
    }


def analyse_f1(a: int, b: int, m: int, n: int) -> list[dict]:
    """Factor f1 and classify each prime."""
    if factorint is None:
        raise RuntimeError("sympy is required for f1 factorization")
    u1 = a * a - b * b
    v1 = 2 * a * b
    w1 = a * a + b * b
    u2 = m * m - n * n
    v2 = 2 * m * n
    f1 = (w1 * u2) ** 2 + (u1 * v2) ** 2
    factors = factorint(f1)
    rows = []
    for p, e in sorted(factors.items()):
        row = classify_prime(a, b, m, n, int(p))
        row["exponent_total"] = int(e)
        rows.append(row)
    return rows


def main() -> int:
    if len(sys.argv) != 5:
        print("Usage: python gaussian_factorization.py a b m n")
        return 1

    a, b, m, n = (int(x) for x in sys.argv[1:5])

    if gcd(a, b) != 1 or (a - b) % 2 == 0:
        print(f"Invalid Euclid pair (a, b) = ({a}, {b})")
        return 1
    if gcd(m, n) != 1 or (m - n) % 2 == 0:
        print(f"Invalid Euclid pair (m, n) = ({m}, {n})")
        return 1

    rows = analyse_f1(a, b, m, n)

    print(f"Prime factorization of f1 for (a,b,m,n) = ({a},{b},{m},{n}):")
    print()
    print(f"  {'prime':>12s}  {'exp':>4s}  {'case':>5s}  "
          f"{'v_pi':>5s}  {'v_pi_bar':>9s}  blocker")
    print(f"  {'-' * 12}  {'-' * 4}  {'-' * 5}  "
          f"{'-' * 5}  {'-' * 9}  -------")
    for row in rows:
        flag = "YES" if row["is_blocker"] else "no"
        case = row["case"]
        vp = row["v_pi"] if case != "inert" else "-"
        vpb = row["v_pi_bar"] if case != "inert" else "-"
        exp = row["exponent_total"]
        print(f"  {row['prime']:>12d}  {exp:>4d}  {case:>5s}  "
              f"{str(vp):>5s}  {str(vpb):>9s}  {flag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
