"""
Verify the odd-exponent blocker theorem (Theorem 4.1 of Paper 4):

    Every fully factored Master-Hit must have num_blockers >= 1.

A Master-Hit with num_blockers = 0 would be a counterexample.

Usage:
    python3 theorem_check.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def main():
    conn = pub_db.connect()
    cur = conn.cursor()

    cur.execute("SELECT count(*) FROM pub.master_hits WHERE num_blockers = 0")
    n_violations = cur.fetchone()[0]

    cur.execute("SELECT count(*) FROM pub.master_hits WHERE num_blockers IS NOT NULL")
    n_factored = cur.fetchone()[0]

    cur.execute("SELECT count(*) FROM pub.master_hits WHERE num_blockers IS NULL")
    n_unfactored = cur.fetchone()[0]

    cur.execute("""
        SELECT min(num_blockers), max(num_blockers), avg(num_blockers)
        FROM pub.master_hits WHERE num_blockers IS NOT NULL
    """)
    nb_min, nb_max, nb_avg = cur.fetchone()

    print(f"Factored Master-Hits:   {n_factored}")
    print(f"Unfactored Master-Hits: {n_unfactored}")
    print(f"")
    print(f"Theorem violations (num_blockers = 0): {n_violations}")
    print(f"")
    print(f"Blocker count statistics (factored only):")
    print(f"  min: {nb_min}, max: {nb_max}, avg: {nb_avg:.2f}")

    if n_violations > 0:
        print("\nVIOLATIONS DETECTED — perfect cuboid candidates:")
        cur.execute("""
            SELECT id, a, b, m, n, x, y, z, f1
            FROM pub.master_hits WHERE num_blockers = 0 LIMIT 20
        """)
        for row in cur:
            print(f"  {row}")
    else:
        print("\nTheorem holds on all factored Master-Hits.")

    conn.close()
    return 1 if n_violations > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
