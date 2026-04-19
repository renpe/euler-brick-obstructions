"""
Export the relevant tables for Peschmann (2026) Paper 2 from the
internal PostgreSQL database into a portable SQLite file.

This is a one-time developer script; the generated
``data/euler_cuboid.sqlite`` is bundled with the scripts so that the
analysis code runs without access to the original PostgreSQL instance.

Usage:
    python export_sqlite.py
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

try:
    import psycopg
except ImportError:
    psycopg = None


PG_CONNINFO = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"


def export(target: Path) -> None:
    if psycopg is None:
        raise RuntimeError("psycopg is required to read the source PostgreSQL database")

    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()

    sqlite = sqlite3.connect(target)
    sqlite.executescript("""
        CREATE TABLE master_hits (
            a INTEGER NOT NULL,
            b INTEGER NOT NULL,
            m INTEGER NOT NULL,
            n INTEGER NOT NULL,
            master_sqrt TEXT NOT NULL,
            f1 TEXT NOT NULL,
            PRIMARY KEY (a, b, m, n)
        );

        CREATE TABLE f1_factors (
            a INTEGER NOT NULL,
            b INTEGER NOT NULL,
            m INTEGER NOT NULL,
            n INTEGER NOT NULL,
            prime TEXT NOT NULL,
            exponent INTEGER NOT NULL,
            is_blocker INTEGER NOT NULL,
            PRIMARY KEY (a, b, m, n, prime)
        );

        CREATE TABLE chabauty_scan (
            m INTEGER NOT NULL,
            n INTEGER NOT NULL,
            has_finite_point INTEGER NOT NULL,
            points_raw TEXT,
            PRIMARY KEY (m, n)
        );
    """)
    sqlite.commit()

    with psycopg.connect(PG_CONNINFO) as pg:
        _export_master_hits(pg, sqlite)
        _export_f1_factors(pg, sqlite)
        _export_chabauty_scan(pg, sqlite)

    sqlite.commit()
    sqlite.close()


def _export_master_hits(pg, sqlite: sqlite3.Connection) -> None:
    cur = pg.cursor()
    cur.execute("SELECT a, b, m, n, q, f1 FROM master_hits")
    rows = cur.fetchall()
    sqlite.executemany(
        "INSERT INTO master_hits (a, b, m, n, master_sqrt, f1) VALUES (?, ?, ?, ?, ?, ?)",
        [(int(a), int(b), int(m), int(n), str(q), str(f1))
         for a, b, m, n, q, f1 in rows],
    )
    print(f"  master_hits: {len(rows)} rows")


def _export_f1_factors(pg, sqlite: sqlite3.Connection) -> None:
    cur = pg.cursor()
    cur.execute("""
        SELECT h.a, h.b, h.m, h.n, f.prime, f.exponent, f.is_blocker
        FROM master_hits h
        JOIN f1_factors f ON f.hit_id = h.id
    """)
    rows = cur.fetchall()
    sqlite.executemany(
        "INSERT INTO f1_factors (a, b, m, n, prime, exponent, is_blocker) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        [(int(a), int(b), int(m), int(n), str(p), int(e), 1 if blk else 0)
         for a, b, m, n, p, e, blk in rows],
    )
    print(f"  f1_factors: {len(rows)} rows")


def _export_chabauty_scan(pg, sqlite: sqlite3.Connection) -> None:
    cur = pg.cursor()
    cur.execute("""
        SELECT m, n, has_finite_point, points_raw
        FROM chabauty_vollscan
        WHERE status = 'done'
    """)
    rows = cur.fetchall()
    sqlite.executemany(
        "INSERT INTO chabauty_scan (m, n, has_finite_point, points_raw) "
        "VALUES (?, ?, ?, ?)",
        [(int(m), int(n), 1 if has else 0, pts)
         for m, n, has, pts in rows],
    )
    print(f"  chabauty_scan: {len(rows)} rows")


def main() -> int:
    target = Path(__file__).parent / "data" / "euler_cuboid.sqlite"
    print(f"Exporting to {target}")
    export(target)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
