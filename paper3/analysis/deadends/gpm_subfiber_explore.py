"""
Erkundung: wie verteilen sich Hits aus einer (g₊, g₋)-Faser über (m, n)?

Wir wählen z.B. (g₊=9, g₋=1) — die größte Saund⁺-Faser mit 14086 Hits — und
schauen, ob die Punkte in wenigen (m,n)-Subfasern konzentriert sind oder
breit gestreut. Falls konzentriert: wir können pro (g₊, g₋, m, n) eine
elliptische (Sub-)Kurve aufstellen und Mordell-Weil rechnen.

Aufruf:
    python3 gpm_subfiber_explore.py [g_plus] [g_minus]
    Default: 9 1
"""
from __future__ import annotations

import os
import sys
from collections import Counter
from math import gcd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def main():
    target_gp = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    target_gm = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    conn = pub_db.connect()
    cur = conn.cursor()
    cur.execute("""
        SELECT id, a, b, m, n, num_blockers
        FROM pub.master_hits
    """)
    rows = cur.fetchall()

    in_fiber = []
    for hid, a, b, m, n, nb in rows:
        a, b, m, n = int(a), int(b), int(m), int(n)
        A = a*m + b*n; B = a*n + b*m
        C = a*m - b*n; D = a*n - b*m
        gp = gcd(A, B); gm = gcd(abs(C), abs(D))
        if gp == target_gp and gm == target_gm:
            in_fiber.append((int(hid), a, b, m, n, nb))

    print(f"(g₊={target_gp}, g₋={target_gm}): {len(in_fiber)} Hits.\n")

    # --- Verteilung über (m, n) -------------------------------------
    mn_counter = Counter((m, n) for _, _, _, m, n, _ in in_fiber)
    print(f"Anzahl distinkter (m, n)-Subfasern: {len(mn_counter)}")
    print(f"Maximale Hits pro (m, n):           {max(mn_counter.values())}")
    print(f"Median Hits pro (m, n):             "
          f"{sorted(mn_counter.values())[len(mn_counter)//2]}")
    print(f"Subfasern mit nur 1 Hit:            "
          f"{sum(1 for v in mn_counter.values() if v == 1)}")
    print(f"Subfasern mit ≥10 Hits:             "
          f"{sum(1 for v in mn_counter.values() if v >= 10)}")
    print()

    print("=== Top 15 (m, n)-Subfasern in dieser Faser ===")
    print(f"{'m':>5} {'n':>5} {'#hits':>6} {'min_b':>6} {'max_b':>6}")
    for (m, n), c in mn_counter.most_common(15):
        bs = [int(nb) for hid, a, b, mm, nn, nb in in_fiber
              if mm == m and nn == n and nb is not None]
        if bs:
            print(f"{m:>5} {n:>5} {c:>6} {min(bs):>6} {max(bs):>6}")
        else:
            print(f"{m:>5} {n:>5} {c:>6} {'-':>6} {'-':>6}")

    # --- Verteilung über (a, b) -------------------------------------
    print()
    ab_counter = Counter((a, b) for _, a, b, _, _, _ in in_fiber)
    print(f"Anzahl distinkter (a, b)-Werte:     {len(ab_counter)}")
    print(f"Top (a, b) mit den meisten Hits (= verschiedene m,n):")
    for (a, b), c in ab_counter.most_common(10):
        print(f"  (a,b)=({a},{b}): {c} verschiedene (m,n)")

    # --- Single-Blocker innerhalb dieser Faser ----------------------
    sb = [(hid, a, b, m, n) for hid, a, b, m, n, nb in in_fiber if nb == 1]
    print(f"\n=== Single-Blocker in (g₊={target_gp}, g₋={target_gm}) ===")
    print(f"Anzahl: {len(sb)}")
    for hid, a, b, m, n in sb[:20]:
        print(f"  hit_id={hid}: (a,b,m,n)=({a},{b},{m},{n})")

    conn.close()


if __name__ == "__main__":
    main()
