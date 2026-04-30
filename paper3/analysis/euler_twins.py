"""
Detailanalyse der Euler-Brick-Zwillinge.

Empirisch: in 1.28M primitiven Bricks gibt es genau 38 Paare (x_prim, y_prim),
die mit ≥2 verschiedenen z_prim einen Euler-Brick bilden. Diese „Zwillinge"
sind die einzigen Punkte, an denen die (x,y) → z-Faserung mehrdeutig wird.

Wir suchen nach algebraischen Mustern:
  - Verhältnis z₁/z₂, z₁·z₂, z₁+z₂, z₂−z₁
  - ggT-Strukturen zwischen z₁ und z₂
  - Verhältnis f1(brick1) / f1(brick2)
  - Beziehung zu (a,b,m,n)-Parametern: liegen die zwei Bricks in derselben
    Master-Faser oder verschiedenen?
  - Existieren Tripletts (3+ z's für dasselbe (x,y))?
"""
from __future__ import annotations

import os
import sys
from collections import defaultdict
from math import gcd, isqrt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def is_square(n):
    if n < 0:
        return False
    r = isqrt(n)
    return r*r == n


def main():
    conn = pub_db.connect()
    cur = conn.cursor()
    print("Lade alle Master-Hits mit (a,b,m,n) und primitive Bricks...")
    cur.execute("""
        SELECT id, a, b, m, n, x_prim, y_prim, z_prim, g_scale
        FROM pub.master_hits
    """)
    rows = cur.fetchall()
    print(f"Geladen: {len(rows)} Hits.\n")

    # Sortiere (x, y) und sammle alle z's
    xy_to_entries = defaultdict(list)
    for hid, a, b, m, n, x, y, z, g in rows:
        x, y, z = int(x), int(y), int(z)
        a, b = sorted([x, y])  # vereinheitliche
        xy_to_entries[(a, b)].append({
            'hit_id': int(hid),
            'a': int(a), 'b': int(b), 'm': int(m), 'n': int(n),
            'x': x, 'y': y, 'z': z, 'g_scale': int(g),
        })

    twins = {xy: e for xy, e in xy_to_entries.items() if len({en['z'] for en in e}) > 1}
    triplets = {xy: e for xy, e in xy_to_entries.items() if len({en['z'] for en in e}) > 2}
    print(f"Zwillinge (≥2 verschiedene z): {len(twins)}")
    print(f"Tripletts (≥3 verschiedene z): {len(triplets)}")
    if triplets:
        print("★★★ Triplett(s) gefunden! ★★★")
        for xy, ee in triplets.items():
            print(f"  ({xy}): zs = {sorted({e['z'] for e in ee})}")
    print()

    # Für die 38 Zwillinge: detaillierte Analyse
    print(f"=== Detail-Analyse der {len(twins)} Zwillingspaare ===\n")

    twin_records = []
    for (x_prim, y_prim), entries in sorted(twins.items()):
        # Eindeutige (x,y,z)-Bricks (ein z kann mehrere (a,b,m,n)-Tupel haben)
        z_to_entry = {}
        for e in entries:
            if e['z'] not in z_to_entry:
                z_to_entry[e['z']] = e
        zs = sorted(z_to_entry.keys())
        if len(zs) != 2:
            continue
        z1, z2 = zs
        e1, e2 = z_to_entry[z1], z_to_entry[z2]
        f1_brick1 = x_prim**2 + y_prim**2 + z1**2
        f1_brick2 = x_prim**2 + y_prim**2 + z2**2
        twin_records.append({
            'x': x_prim, 'y': y_prim, 'z1': z1, 'z2': z2,
            'gcd_z1z2': gcd(z1, z2),
            'z2_over_z1': z2 / z1,
            'z2_minus_z1': z2 - z1,
            'z2_plus_z1': z2 + z1,
            'z2_times_z1': z2 * z1,
            'f1_1': f1_brick1, 'f1_2': f1_brick2,
            'f1_ratio': f1_brick2 / f1_brick1,
            'mn_1': (e1['m'], e1['n']), 'mn_2': (e2['m'], e2['n']),
            'g_scale_1': e1['g_scale'], 'g_scale_2': e2['g_scale'],
            'ab_1': (e1['a'], e1['b']), 'ab_2': (e2['a'], e2['b']),
        })

    # Tabellen-Ausgabe
    print(f"{'idx':>3} {'x':>10} {'y':>10} {'z1':>10} {'z2':>10} "
          f"{'gcd':>6} {'z2/z1':>8} {'mn1':>10} {'mn2':>10}")
    for i, r in enumerate(twin_records):
        print(f"{i:>3} {r['x']:>10} {r['y']:>10} {r['z1']:>10} {r['z2']:>10} "
              f"{r['gcd_z1z2']:>6} {r['z2_over_z1']:>8.3f} "
              f"{str(r['mn_1']):>10} {str(r['mn_2']):>10}")

    # Suche nach algebraischer Beziehung z1, z2 ↔ x, y
    print("\n=== Algebraische Tests ===")
    # Test 1: z1 + z2 = c * (x + y)?  oder z1*z2 = c * x*y?
    print("\nTest: ist z₁·z₂ irgendein 'einfacher' Wert in x,y?")
    print(f"{'idx':>3} {'z1·z2':>20} {'x²·y²':>20} {'x²+y²':>15} {'(x²+y²)²':>20} "
          f"{'z1·z2 / x²y²':>15}")
    for i, r in enumerate(twin_records[:15]):
        x, y = r['x'], r['y']
        z1z2 = r['z1'] * r['z2']
        x2y2 = x*x*y*y
        ratio = z1z2 / x2y2 if x2y2 else 0
        print(f"{i:>3} {z1z2:>20} {x2y2:>20} {x*x+y*y:>15} "
              f"{(x*x+y*y)**2:>20} {ratio:>15.6f}")

    # Test: z1 + z2 vs. einfache Polynome in x,y
    print("\nTest: ist z₁+z₂ einfach?")
    print(f"{'idx':>3} {'z1+z2':>15} {'2·x':>10} {'2·y':>10} {'x+y':>10} "
          f"{'x²+y²':>15} {'sqrt(z1·z2)':>15}")
    for i, r in enumerate(twin_records[:15]):
        x, y = r['x'], r['y']
        sum_z = r['z1'] + r['z2']
        # Suchen "einfacher" Ausdruck: ist sum_z teilbar durch x oder y?
        flag = ""
        if sum_z % x == 0:
            flag += " /x"
        if sum_z % y == 0:
            flag += " /y"
        if (x*x + y*y) > 0 and sum_z % (x*x + y*y) == 0:
            flag += " /(x²+y²)"
        sqrt_prod = isqrt(r['z1'] * r['z2'])
        is_sq = (sqrt_prod * sqrt_prod == r['z1'] * r['z2'])
        sq_str = "□" if is_sq else " "
        print(f"{i:>3} {sum_z:>15} {2*x:>10} {2*y:>10} {x+y:>10} "
              f"{x*x+y*y:>15} {sqrt_prod:>15}{sq_str}{flag}")

    # f1-Quotient: ist f1_2 / f1_1 = einfache Rationalzahl?
    print("\nTest: ist f1_2 / f1_1 ein 'schönes' Verhältnis?")
    print(f"{'idx':>3} {'f1_1':>15} {'f1_2':>15} {'ratio':>10} {'p/q vereinfacht':>20}")
    from fractions import Fraction
    for i, r in enumerate(twin_records[:20]):
        f = Fraction(r['f1_2'], r['f1_1'])
        print(f"{i:>3} {r['f1_1']:>15} {r['f1_2']:>15} "
              f"{r['f1_ratio']:>10.4f} {f.numerator}/{f.denominator}")

    # Brahmagupta-Test: gibt es eine quadratische Gleichung in z, deren
    # zwei Lösungen z1 und z2 sind, mit Koeffizienten in x,y?
    print("\nTest: Quadratische Gleichung az² + bz + c = 0 mit Roots z1, z2")
    print("       z1+z2 = -b/a, z1·z2 = c/a")
    print(f"{'idx':>3} {'z1+z2':>15} {'z1·z2':>20} {'(z1+z2)²−4·z1·z2':>20}")
    for i, r in enumerate(twin_records[:15]):
        sumz = r['z1'] + r['z2']
        prodz = r['z1'] * r['z2']
        disc = sumz*sumz - 4*prodz
        # ist disc = (z2-z1)² → ja immer per Vieta
        print(f"{i:>3} {sumz:>15} {prodz:>20} {disc:>20}  "
              f"= (z₂−z₁)² = {(r['z2']-r['z1'])**2}")

    conn.close()


if __name__ == "__main__":
    main()
