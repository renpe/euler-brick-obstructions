"""
Diagnostic of the (g_+, g_-) fibration over pub.master_hits.

Background: the structure theorem says g_+ * g_- = g_scale. Hence all
primitive bricks (g_scale=1) lie in the same fiber (1,1) - the fibration
only separates non-primitive tuples. This diagnostic checks whether the
fibration reveals algebraically interesting classes, in particular:

  - What is the distribution of (g_+, g_-) pairs?
  - Does it correlate with num_blockers?
  - Do the known 1-blocker hits (Saunderson candidates) lie in
    a characteristic fiber?

Usage:
    python3 gpm_fibration.py
"""
from __future__ import annotations

import os
import sys
from collections import Counter, defaultdict
from math import gcd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def gpm(a: int, b: int, m: int, n: int) -> tuple[int, int]:
    """Returns (g_plus, g_minus) for a master tuple."""
    A = a*m + b*n
    B = a*n + b*m
    C = a*m - b*n
    D = a*n - b*m
    g_plus = gcd(A, B)
    g_minus = gcd(abs(C), abs(D))
    return g_plus, g_minus


def main():
    conn = pub_db.connect()
    cur = conn.cursor()
    cur.execute("""
        SELECT id, a, b, m, n, g_scale, num_blockers
        FROM pub.master_hits
        ORDER BY id
    """)
    rows = cur.fetchall()
    print(f"Loaded {len(rows)} master hits.\n")

    # --- 1. (g_+, g_-) distribution + sanity check g_+ * g_- == g_scale ----
    pair_count = Counter()
    pair_to_blockers = defaultdict(list)   # (gp, gm) -> [num_blockers, ...]
    pair_to_hit_ids = defaultdict(list)
    sanity_violations = 0
    for hid, a, b, m, n, g_scale, nb in rows:
        gp, gm = gpm(int(a), int(b), int(m), int(n))
        if gp * gm != int(g_scale):
            sanity_violations += 1
        pair_count[(gp, gm)] += 1
        pair_to_blockers[(gp, gm)].append(nb)
        pair_to_hit_ids[(gp, gm)].append(int(hid))
    print(f"Sanity check g_+ * g_- == g_scale: violations = {sanity_violations}\n")

    # --- 2. Top fibers by size ------------------------------------------
    print("=== Top 20 (g_+, g_-) fibers by hit count ===")
    print(f"{'g_+':>6} {'g_-':>6} {'#hits':>10} {'min_b':>6} {'max_b':>6} "
          f"{'mean_b':>8} {'unfact':>7}")
    for (gp, gm), c in pair_count.most_common(20):
        bs = pair_to_blockers[(gp, gm)]
        nfact = [b for b in bs if b is not None]
        nun = len(bs) - len(nfact)
        if nfact:
            mn = min(nfact); mx = max(nfact); mean = sum(nfact)/len(nfact)
        else:
            mn = mx = mean = float('nan')
        print(f"{gp:>6} {gm:>6} {c:>10} {mn:>6} {mx:>6} {mean:>8.2f} {nun:>7}")

    # --- 3. Saunderson localization: 1-blocker hits ---------------------
    print("\n=== Saunderson localization (num_blockers = 1) ===")
    saunderson_pairs = Counter()
    for (gp, gm), bs in pair_to_blockers.items():
        n1 = sum(1 for b in bs if b == 1)
        if n1 > 0:
            saunderson_pairs[(gp, gm)] = n1
    print(f"In total {sum(saunderson_pairs.values())} 1-blocker hits "
          f"in {len(saunderson_pairs)} different fibers.\n")
    print(f"{'g_+':>6} {'g_-':>6} {'#1-blocker':>11}")
    for (gp, gm), n in saunderson_pairs.most_common(20):
        print(f"{gp:>6} {gm:>6} {n:>11}")

    # --- 4. Primitive vs. non-primitive --------------------------------
    n_prim = pair_count[(1, 1)]
    n_nonprim = sum(c for k, c in pair_count.items() if k != (1, 1))
    print(f"\nPrimitive (g_+, g_-) = (1,1):   {n_prim:>10}")
    print(f"Non-primitive:               {n_nonprim:>10}")
    print(f"Number of distinct fibers:   {len(pair_count):>10}")

    # --- 5. Blocker distribution per "class" ---------------------------
    # Class 1: (1,1) primitive
    # Class 2: (k,1) - Saunderson candidate (g_+ scaled, g_- = 1)
    # Class 3: (1,k) - dual Saunderson candidate
    # Class 4: both >1
    print("\n=== Blocker distribution per fiber class ===")
    classes = {
        "(1,1) primitive":     lambda gp, gm: gp == 1 and gm == 1,
        "(k,1) k>1 'Saund+'":  lambda gp, gm: gp > 1 and gm == 1,
        "(1,k) k>1 'Saund-'":  lambda gp, gm: gp == 1 and gm > 1,
        "both >1":             lambda gp, gm: gp > 1 and gm > 1,
    }
    print(f"{'class':<24} {'#hits':>8} {'min_b':>6} {'max_b':>6} "
          f"{'#1-block':>9} {'unfact':>7}")
    for label, pred in classes.items():
        bs = []
        n1 = nun = nhit = 0
        for (gp, gm), bb in pair_to_blockers.items():
            if not pred(gp, gm):
                continue
            nhit += pair_count[(gp, gm)]
            for b in bb:
                if b is None:
                    nun += 1
                else:
                    bs.append(b)
                    if b == 1:
                        n1 += 1
        if bs:
            print(f"{label:<24} {nhit:>8} {min(bs):>6} {max(bs):>6} "
                  f"{n1:>9} {nun:>7}")
        else:
            print(f"{label:<24} {nhit:>8} {'-':>6} {'-':>6} {n1:>9} {nun:>7}")

    conn.close()


if __name__ == "__main__":
    main()
