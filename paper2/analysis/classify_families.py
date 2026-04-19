"""
Parametric family classification for Master-Hits.

Checks whether a given Master-Hit (a, b, m, n) is produced by the
Saunderson closed-form generator (Theorem 5.1) or by any of the
three Himane two-triple constructions (Section 5).

Runs a bulk coverage check against the SQLite database produced by
`export_sqlite.py` (data/euler_cuboid.sqlite).

Usage:
    python classify_families.py                     # summary over DB
    python classify_families.py <a> <b> <m> <n>     # single tuple
"""

from __future__ import annotations

import sqlite3
import sys
from math import gcd, isqrt
from pathlib import Path


# ---------------------------------------------------------------------------
# Euclid-pair utilities
# ---------------------------------------------------------------------------

def is_valid_euclid_pair(p: int, q: int) -> bool:
    return p > q > 0 and gcd(p, q) == 1 and (p - q) % 2 == 1


def primitive_brick_edges(a: int, b: int, m: int, n: int) -> tuple[int, int, int]:
    """Return the three edges (U1*U2, V1*U2, U1*V2) divided by their gcd."""
    u1 = a * a - b * b
    v1 = 2 * a * b
    u2 = m * m - n * n
    v2 = 2 * m * n
    e1, e2, e3 = u1 * u2, v1 * u2, u1 * v2
    g = gcd(gcd(abs(e1), abs(e2)), abs(e3))
    if g == 0:
        raise ValueError("degenerate brick")
    return tuple(sorted((abs(e1) // g, abs(e2) // g, abs(e3) // g)))


# ---------------------------------------------------------------------------
# Saunderson
# ---------------------------------------------------------------------------

def matches_saunderson(a: int, b: int, m: int, n: int) -> tuple[int, int] | None:
    """If (a, b, m, n) is produced by Saunderson's formula, return (g, h).
    Otherwise None.
    """
    for g in range(2, 2 * max(m, a) + 1):
        for h in range(1, g):
            if gcd(g, h) != 1 or (g - h) % 2 == 0:
                continue
            if (4 * g * h, g * g + h * h) != (a, b):
                continue
            m_s = h * (3 * g * g - h * h)
            n_s = g * (g * g - 3 * h * h)
            if n_s < 0:
                m_s, n_s = n_s * -1, m_s
                m_s, n_s = max(m_s, n_s), min(m_s, n_s)
            if (m, n) == (m_s, n_s):
                return (g, h)
    return None


# ---------------------------------------------------------------------------
# Himane (Theorems 1, 2, 3)
# ---------------------------------------------------------------------------

def iterate_primitive_triples(max_generator: int):
    for g in range(2, max_generator + 1):
        for h in range(1, g):
            if gcd(g, h) != 1 or (g - h) % 2 == 0:
                continue
            u = g * g - h * h
            v = 2 * g * h
            w = g * g + h * h
            yield (u, v, w)


def generate_himane_bricks(max_generator: int = 80, max_scale: int = 12):
    """Generate primitive edge-sorted brick tuples from Himane's three
    theorems. Returns a dict {brick: (theorem, metadata)}.
    """
    triples = list(iterate_primitive_triples(max_generator))
    bricks: dict[tuple[int, int, int], tuple[str, dict]] = {}

    def register(brick, theorem, meta):
        if brick not in bricks:
            bricks[brick] = (theorem, meta)

    for u1, v1, w1 in triples:
        for u2, v2, w2 in triples:
            for k in range(1, max_scale + 1):
                # Theorem 1: v0*u1 = u0*v2  ->  u0 = u1*k, v0 = v2*k
                u0 = u1 * k
                v0 = v2 * k
                a = u0 * u2
                b = v0 * u1
                c = v0 * v1
                if a > 0 and b > 0 and c > 0:
                    brick = _primitive_triple(a, b, c)
                    if brick is not None:
                        register(brick, "himane_thm1",
                                 {"triples": ((u1, v1, w1), (u2, v2, w2)),
                                  "k": k})

                # Theorem 2: u0*u2 = v0*u1  ->  u0 = u1*k, v0 = u2*k
                u0 = u1 * k
                v0 = u2 * k
                a = u0 * u2
                b = v0 * v1
                c = u0 * v2
                if a > 0 and b > 0 and c > 0:
                    brick = _primitive_triple(a, b, c)
                    if brick is not None:
                        register(brick, "himane_thm2",
                                 {"triples": ((u1, v1, w1), (u2, v2, w2)),
                                  "k": k})

                # Theorem 3: v0*v1 = u0*v2  ->  u0 = v1*k, v0 = v2*k
                u0 = v1 * k
                v0 = v2 * k
                a = u0 * u2
                b = v0 * u1
                c = u0 * v2
                if a > 0 and b > 0 and c > 0:
                    brick = _primitive_triple(a, b, c)
                    if brick is not None:
                        register(brick, "himane_thm3",
                                 {"triples": ((u1, v1, w1), (u2, v2, w2)),
                                  "k": k})
    return bricks


def _primitive_triple(a: int, b: int, c: int) -> tuple[int, int, int] | None:
    g = gcd(gcd(abs(a), abs(b)), abs(c))
    if g == 0:
        return None
    return tuple(sorted((abs(a) // g, abs(b) // g, abs(c) // g)))


def matches_himane(a: int, b: int, m: int, n: int, bricks) -> tuple[str, dict] | None:
    """If (a, b, m, n) produces a brick found in the Himane catalogue,
    return (theorem_name, metadata). Otherwise None.
    """
    brick = primitive_brick_edges(a, b, m, n)
    return bricks.get(brick)


# ---------------------------------------------------------------------------
# Database coverage
# ---------------------------------------------------------------------------

def classify_database(db_path: Path, max_generator: int = 80,
                      max_scale: int = 12) -> dict[str, int]:
    bricks = generate_himane_bricks(max_generator, max_scale)
    counts = {"saunderson": 0, "himane_thm1": 0, "himane_thm2": 0,
              "himane_thm3": 0, "sporadic": 0, "total": 0}

    conn = sqlite3.connect(db_path)
    for a, b, m, n in conn.execute("SELECT a, b, m, n FROM master_hits"):
        counts["total"] += 1
        if matches_saunderson(a, b, m, n) is not None:
            counts["saunderson"] += 1
            continue
        hit = matches_himane(a, b, m, n, bricks)
        if hit is not None:
            counts[hit[0]] += 1
            continue
        counts["sporadic"] += 1
    conn.close()
    return counts


def main() -> int:
    if len(sys.argv) == 5:
        a, b, m, n = (int(x) for x in sys.argv[1:5])
        if not (is_valid_euclid_pair(a, b) and is_valid_euclid_pair(m, n)):
            print(f"Invalid Euclid pairs: ({a},{b}), ({m},{n})")
            return 1
        saund = matches_saunderson(a, b, m, n)
        if saund is not None:
            print(f"Saunderson with generators (g,h) = {saund}")
            return 0
        bricks = generate_himane_bricks()
        hit = matches_himane(a, b, m, n, bricks)
        if hit is not None:
            theorem, meta = hit
            print(f"Himane match: {theorem} with {meta}")
            return 0
        print("No known parametric family matches this tuple.")
        return 0

    db_path = Path(__file__).parent.parent / "data" / "euler_cuboid.sqlite"
    if not db_path.exists():
        print(f"Database not found: {db_path}")
        print("Run export_sqlite.py first.")
        return 1

    print(f"Classifying Master-Hits in {db_path} against parametric families...")
    counts = classify_database(db_path)

    print()
    total = counts["total"]
    for family in ("saunderson", "himane_thm1", "himane_thm2",
                   "himane_thm3", "sporadic"):
        pct = 100 * counts[family] / total if total else 0
        print(f"  {family:>15s}: {counts[family]:>6d}  ({pct:5.2f}%)")
    print(f"  {'total':>15s}: {total:>6d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
