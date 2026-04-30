"""
Import the rigorous rank results from pari_rank_scan.sage into pub.fibers.

Reads scripts/paper4/data/pari_ranks.csv and updates the corresponding
rows in pub.fibers (creating them if missing). The CSV contains:

    m, n, conductor_digits, rank_lower, rank_upper, sha2_dim, n_gens, error

We populate:
    pub.fibers.rank_known   = rank_lower (the proven lower bound)
    pub.fibers.rank_proven  = (rank_lower == rank_upper)

Idempotent: re-running over the same CSV updates the same rows.

Usage:
    python3 import_pari_ranks.py
"""
import os
import sys
import csv
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db

CSV_PATH = "scripts/paper4/data/pari_ranks.csv"


def main():
    if not os.path.exists(CSV_PATH):
        print(f"CSV not found: {CSV_PATH}")
        return 1

    conn = pub_db.connect(autocommit=False)
    cur = conn.cursor()

    n_updated = n_inserted = n_err = n_skipped = 0
    t0 = time()

    with open(CSV_PATH) as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) < 8: continue
            m, n, cond_d, rl, ru, sha2, ngens, err = row
            if err and err.strip():
                n_err += 1
                continue
            if not rl or not ru:
                n_skipped += 1
                continue
            m, n = int(m), int(n)
            rl, ru = int(rl), int(ru)
            proven = (rl == ru)
            cur.execute("""
                INSERT INTO pub.fibers (m, n, rank_known, rank_proven)
                VALUES (%s, %s, %s, %s)
                ON CONFLICT (m, n) DO UPDATE
                  SET rank_known = EXCLUDED.rank_known,
                      rank_proven = EXCLUDED.rank_proven
                RETURNING (xmax = 0) AS inserted
            """, (m, n, rl, proven))
            row_result = cur.fetchone()
            if row_result and row_result[0]:
                n_inserted += 1
            else:
                n_updated += 1

    conn.commit()
    print(f"Done in {time()-t0:.1f}s")
    print(f"  inserted: {n_inserted}")
    print(f"  updated:  {n_updated}")
    print(f"  errors:   {n_err}")
    print(f"  skipped:  {n_skipped}")
    conn.close()


if __name__ == "__main__":
    sys.exit(main() or 0)
