"""
Shared DB helpers for the pub schema.

Imported by all scripts in scripts/paper4/. Encapsulates:
  - Connection setup with configurable DB
  - Insert of a Master-Hit including hit_groups (provenance)
  - Lookup of group_id for a provenance

Never INSERT INTO pub.master_hits directly - always go through
insert_master_hit() so that provenance is maintained consistently.
"""
from __future__ import annotations

import os
from math import gcd, isqrt
from typing import Optional

import psycopg

# Load .env (repo root scripts/ or paper4/) if python-dotenv is installed.
try:
    from dotenv import load_dotenv  # type: ignore
    _here = os.path.dirname(os.path.abspath(__file__))
    for _candidate in (
        os.path.join(_here, "..", "..", ".env"),  # scripts/ (repo root)
        os.path.join(_here, "..", ".env"),         # paper4/
    ):
        if os.path.isfile(_candidate):
            load_dotenv(_candidate)
            break
except ImportError:
    pass

DEFAULT_DB = os.environ.get(
    "EULER_DB",
    "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler",
)

# Binary paths used by other modules. Override via .env or environment.
YAFU_BIN = os.environ.get("YAFU_BIN", "yafu")
SAGE_BIN = os.environ.get("SAGE_BIN", "sage")
CADO_NFS = os.environ.get("CADO_NFS", "cado-nfs.py")
CADO_VENV_PY = os.environ.get("CADO_VENV_PY", "python3")


def connect(dsn: Optional[str] = None, autocommit: bool = False):
    return psycopg.connect(dsn or DEFAULT_DB, autocommit=autocommit)


# ---------------------------------------------------------------
# Group lookup / creation
# ---------------------------------------------------------------

def get_group_id(cur, short_name: str) -> int:
    """Returns id of generator_groups.short_name. Raises if not found."""
    cur.execute(
        "SELECT id FROM pub.generator_groups WHERE short_name = %s",
        (short_name,),
    )
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"generator_groups.short_name='{short_name}' not found.")
    return row[0]


def ensure_mw_group(cur, m: int, n: int) -> int:
    """Look up MW-{m}-{n}, creating it (and pub.fibers row) if missing.
    Returns the group_id. cur must support commit/rollback at caller.
    """
    short = f"MW-{m}-{n}"
    cur.execute("SELECT id FROM pub.generator_groups WHERE short_name=%s", (short,))
    row = cur.fetchone()
    if row is not None:
        return row[0]

    # Need parent_id for "Mordell-Weil"
    parent_id = get_group_id(cur, "Mordell-Weil")

    # Ensure fibers row exists
    cur.execute(
        "INSERT INTO pub.fibers (m, n) VALUES (%s, %s) ON CONFLICT (m, n) DO NOTHING RETURNING id",
        (m, n),
    )
    row = cur.fetchone()
    if row is None:
        cur.execute("SELECT id FROM pub.fibers WHERE m=%s AND n=%s", (m, n))
        fiber_id = cur.fetchone()[0]
    else:
        fiber_id = row[0]

    cur.execute(
        """INSERT INTO pub.generator_groups
           (short_name, parent_id, category, description, fiber_id)
           VALUES (%s, %s, 'provenance', %s, %s)
           RETURNING id""",
        (short, parent_id, f"Mordell-Weil generator on E_{{{m},{n}}}", fiber_id),
    )
    return cur.fetchone()[0]


# ---------------------------------------------------------------
# Master-Hit insertion with provenance
# ---------------------------------------------------------------

INSERT_HIT_SQL = """
INSERT INTO pub.master_hits
    (a, b, m, n, x, y, z, g_scale, x_prim, y_prim, z_prim, f1, num_blockers)
VALUES
    (%(a)s, %(b)s, %(m)s, %(n)s, %(x)s, %(y)s, %(z)s,
     %(g_scale)s, %(x_prim)s, %(y_prim)s, %(z_prim)s,
     %(f1)s, %(num_blockers)s)
ON CONFLICT (a, b, m, n) DO NOTHING
RETURNING id
"""


def compute_hit_fields(a: int, b: int, m: int, n: int) -> Optional[dict]:
    """Computes all derived fields for a Master-Hit. Returns None if M is not a square."""
    U1 = a*a - b*b
    V1 = 2*a*b
    U2 = m*m - n*n
    V2 = 2*m*n
    M_val = (V1*U2)**2 + (U1*V2)**2
    q = isqrt(M_val)
    if q*q != M_val:
        return None
    U1U2 = U1*U2
    f1 = M_val + U1U2*U1U2
    x = U1*U2; y = V1*U2; z = U1*V2
    g = gcd(gcd(int(x), int(y)), int(z))
    return dict(
        a=int(a), b=int(b), m=int(m), n=int(n),
        x=int(x), y=int(y), z=int(z),
        g_scale=int(g),
        x_prim=int(x)//g, y_prim=int(y)//g, z_prim=int(z)//g,
        f1=int(f1),
        num_blockers=None,
    )


def insert_master_hit(
    cur,
    a: int, b: int, m: int, n: int,
    provenance_group_id: int,
    details: Optional[dict] = None,
) -> Optional[int]:
    """Insert a Master-Hit into pub.master_hits and link it to the given
    provenance group via pub.hit_groups (is_primary=true).

    Returns the inserted hit_id, or None if the hit already existed.
    Caller is responsible for commit/rollback.
    """
    fields = compute_hit_fields(a, b, m, n)
    if fields is None:
        return None
    cur.execute(INSERT_HIT_SQL, fields)
    row = cur.fetchone()
    if row is None:
        return None
    hit_id = row[0]
    cur.execute(
        """INSERT INTO pub.hit_groups (hit_id, group_id, is_primary, details)
           VALUES (%s, %s, true, %s)""",
        (hit_id, provenance_group_id, psycopg.types.json.Jsonb(details) if details else None),
    )
    return hit_id


# ---------------------------------------------------------------
# Factorisation pipeline helper
# ---------------------------------------------------------------

def store_prim_brick_factors(
    cur,
    hit_id: int,
    x_prim: int, y_prim: int, z_prim: int,
    factors: dict[int, int],
) -> int:
    """Write factorisation of f1_prim = x_prim²+y_prim²+z_prim² into
    pub.prim_brick_factors and update num_blockers in pub.master_hits.

    factors: {prime: exponent} where Π p^e == x_prim²+y_prim²+z_prim².

    Returns number of blockers (odd-exponent primes).

    Note: primitive bricks are unique per (x_prim, y_prim, z_prim) — we
    use ON CONFLICT DO NOTHING because if two hits ever produced the same
    primitive brick, both would yield the same factorisation. (Empirically
    1.28M unique primitives in the DB → no conflicts in practice.)
    """
    n_blockers = 0
    for p, e in factors.items():
        is_blocker = (e % 2 == 1)
        if is_blocker:
            n_blockers += 1
        cur.execute(
            """INSERT INTO pub.prim_brick_factors
                   (x_prim, y_prim, z_prim, prime, exponent, is_blocker, prime_mod4)
               VALUES (%s, %s, %s, %s, %s, %s, %s)
               ON CONFLICT (x_prim, y_prim, z_prim, prime) DO UPDATE
                 SET exponent = EXCLUDED.exponent,
                     is_blocker = EXCLUDED.is_blocker,
                     prime_mod4 = EXCLUDED.prime_mod4""",
            (int(x_prim), int(y_prim), int(z_prim),
             int(p), int(e), is_blocker, int(p) % 4),
        )
    cur.execute(
        "UPDATE pub.master_hits SET num_blockers = %s WHERE id = %s",
        (n_blockers, hit_id),
    )
    return n_blockers
