"""
Migrationsskript: public.master_hits → pub.master_hits.

Aufgaben:
  1. Kopiert alle Master-Hits in pub.master_hits (schlanke Spalten).
  2. Mappt search_bound auf eine primäre Provenance:
       200/500/2000/2300 → Exhaustive-Bound-N
       -1               → Saunderson-Generator
       -2               → Rathbun-Search
       -3               → MW-{m}-{n}  (legt MW-Untergruppen on the fly an)
  3. Klassifiziert nach Familien (Saunderson, Lenhart, Himane T1/T2/T3).
  4. Kopiert die f1_factors-Tabelle.
  5. Bei MW-Hits: legt fibers-Eintrag mit Weierstrass/Konduktor/Rang/Torsion an,
     soweit vorhanden (sonst nur (m,n) leer).

Vorbedingung: pub_schema.sql muss bereits ausgeführt sein.

Aufruf:
  .venv-linux/bin/python pub_migrate.py
"""
import sys
from math import isqrt, gcd
from time import time

import psycopg

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"

# ---------- Familien-Generation (wie in extract_family_mn.py) ----------

def is_square(n):
    if n <= 0: return n == 0
    r = isqrt(n)
    return r * r == n

def primitive_brick(e1, e2, e3):
    g = gcd(gcd(abs(e1), abs(e2)), abs(e3))
    if g == 0: return None
    return tuple(sorted([abs(e1)//g, abs(e2)//g, abs(e3)//g]))

def gen_saunderson(max_g):
    bricks = set()
    for g in range(2, max_g+1):
        for h in range(1, g):
            if gcd(g, h) != 1 or (g-h) % 2 == 0: continue
            u, v, w = g*g-h*h, 2*g*h, g*g+h*h
            uu, vv = min(u, v), max(u, v)
            if 3*uu*uu <= vv*vv: continue
            x = uu*(3*vv*vv - uu*uu); y = vv*(3*uu*uu - vv*vv); z = 4*uu*vv*w
            if x <= 0 or y <= 0: continue
            br = primitive_brick(x, y, z)
            if br: bricks.add(br)
    return bricks

def gen_lenhart(max_w):
    bricks = set()
    for w in range(1, max_w+1):
        target = 5*w*w
        for u in range(w+1, isqrt(target)+1):
            v_sq = target - u*u
            if v_sq <= w*w: continue
            if not is_square(v_sq): continue
            v = isqrt(v_sq)
            if v <= w: continue
            a = (u*u-w*w)*(v*v-w*w); b = 4*u*v*w*w
            for c in (2*u*w*(v*v-w*w), 2*v*w*(u*u-w*w)):
                if a > 0 and b > 0 and c > 0:
                    br = primitive_brick(a, b, c)
                    if br and is_square(a*a+b*b) and is_square(a*a+c*c) and is_square(b*b+c*c):
                        bricks.add(br)
    return bricks

def gen_pyth_triples(max_gen):
    return [(g*g-h*h, 2*g*h, g*g+h*h)
            for g in range(2, max_gen+1)
            for h in range(1, g)
            if gcd(g, h) == 1 and (g-h) % 2 == 1]

def gen_himane_per_theorem(triples, max_scale):
    """Returns dict {'T1': set, 'T2': set, 'T3': set}."""
    out = {'T1': set(), 'T2': set(), 'T3': set()}
    for u1, v1, w1 in triples:
        for u2, v2, w2 in triples:
            for k in range(1, max_scale+1):
                # T1
                u0, v0 = u1*k, v2*k
                a, b, c = u0*u2, v0*u1, v0*v1
                if a > 0 and b > 0 and c > 0 and is_square(a*a+c*c):
                    br = primitive_brick(a, b, c)
                    if br: out['T1'].add(br)
                # T2
                u0, v0 = u1*k, u2*k
                a, b, c = u0*u2, v0*v1, u0*v2
                if a > 0 and b > 0 and c > 0 and is_square(b*b+c*c):
                    br = primitive_brick(a, b, c)
                    if br: out['T2'].add(br)
                # T3
                u0, v0 = v1*k, v2*k
                a, b, c = u0*u2, v0*u1, v0*v1
                if a > 0 and b > 0 and c > 0 and is_square(a*a+b*b):
                    br = primitive_brick(a, b, c)
                    if br: out['T3'].add(br)
    return out


def primitive_brick_from_xyz(x, y, z):
    g = gcd(gcd(int(abs(x)), int(abs(y))), int(abs(z)))
    if g == 0: return None
    return tuple(sorted([int(abs(x))//g, int(abs(y))//g, int(abs(z))//g]))


def main():
    t0 = time()
    conn = psycopg.connect(DB, autocommit=False)
    cur = conn.cursor()

    print("Generiere Familien-Bricks...", flush=True)
    saunderson = gen_saunderson(500)
    lenhart    = gen_lenhart(300)
    himane     = gen_himane_per_theorem(gen_pyth_triples(80), 15)
    print(f"  Saunderson: {len(saunderson)}, Lenhart: {len(lenhart)}, "
          f"Himane T1/T2/T3: {len(himane['T1'])}/{len(himane['T2'])}/{len(himane['T3'])}",
          flush=True)

    # Group-IDs nachschlagen
    cur.execute("SELECT short_name, id FROM pub.generator_groups")
    group_id = dict(cur.fetchall())
    mw_parent_id = group_id['Mordell-Weil']

    # ------------------------------------------------------------------
    # Schritt 1: master_hits kopieren
    # ------------------------------------------------------------------
    print("\nKopiere master_hits → pub.master_hits...", flush=True)
    cur.execute("TRUNCATE pub.master_hits CASCADE")
    cur.execute("""
        INSERT INTO pub.master_hits
            (id, a, b, m, n, x, y, z, g_scale, x_prim, y_prim, z_prim, f1, num_blockers)
        SELECT id, a, b, m, n, x, y, z, g_scale, x_prim, y_prim, z_prim, f1,
               CASE WHEN num_blockers >= 0 THEN num_blockers ELSE NULL END
        FROM public.master_hits
    """)
    n_hits = cur.rowcount
    print(f"  {n_hits} Hits kopiert")
    conn.commit()

    # Sequence resetten, damit Inserts mit autoinkrementierender id nicht kollidieren
    cur.execute("SELECT setval(pg_get_serial_sequence('pub.master_hits','id'), MAX(id)) FROM pub.master_hits")
    conn.commit()

    # ------------------------------------------------------------------
    # Schritt 2: f1_factors kopieren
    # ------------------------------------------------------------------
    print("\nKopiere f1_factors → pub.f1_factors...", flush=True)
    cur.execute("TRUNCATE pub.f1_factors")
    cur.execute("""
        INSERT INTO pub.f1_factors (hit_id, prime, exponent, is_blocker, prime_mod4)
        SELECT hit_id, prime, exponent, is_blocker, prime_mod4
        FROM public.f1_factors
    """)
    n_fact = cur.rowcount
    print(f"  {n_fact} Faktor-Zeilen kopiert")
    conn.commit()

    # ------------------------------------------------------------------
    # Schritt 3: MW-Untergruppen anlegen
    # ------------------------------------------------------------------
    print("\nLege MW-Untergruppen an...", flush=True)
    cur.execute("""
        SELECT DISTINCT m, n FROM public.master_hits WHERE search_bound = -3
    """)
    mw_pairs = [(int(m), int(n)) for m, n in cur.fetchall()]
    for m, n in mw_pairs:
        short = f"MW-{m}-{n}"
        cur.execute("""
            INSERT INTO pub.fibers (m, n) VALUES (%s, %s)
            ON CONFLICT (m, n) DO NOTHING
            RETURNING id
        """, (m, n))
        row = cur.fetchone()
        if row is None:
            cur.execute("SELECT id FROM pub.fibers WHERE m=%s AND n=%s", (m, n))
            fiber_id = cur.fetchone()[0]
        else:
            fiber_id = row[0]
        cur.execute("""
            INSERT INTO pub.generator_groups (short_name, parent_id, category, description, fiber_id)
            VALUES (%s, %s, 'provenance', %s, %s)
            ON CONFLICT (short_name) DO NOTHING
        """, (short, mw_parent_id, f"Mordell-Weil generator on E_{{{m},{n}}}", fiber_id))
    print(f"  {len(mw_pairs)} MW-Untergruppen angelegt", flush=True)
    conn.commit()

    # Group-IDs nochmal nachschlagen (mit den neuen MW-X-Y)
    cur.execute("SELECT short_name, id FROM pub.generator_groups")
    group_id = dict(cur.fetchall())

    # ------------------------------------------------------------------
    # Schritt 4: hit_groups – primary provenance
    # ------------------------------------------------------------------
    print("\nWeise primäre Provenance zu...", flush=True)
    cur.execute("TRUNCATE pub.hit_groups")

    # 200/500/2000/2300 → Exhaustive-Bound-N
    for bound in (200, 500, 2000, 2300):
        gid = group_id.get(f'Exhaustive-Bound-{bound}')
        if gid is None: continue
        cur.execute("""
            INSERT INTO pub.hit_groups (hit_id, group_id, is_primary)
            SELECT id, %s, true FROM public.master_hits WHERE search_bound = %s
        """, (gid, bound))
        print(f"  Exhaustive-Bound-{bound}: {cur.rowcount} Hits")

    # -1 → Saunderson-Generator
    gid = group_id['Saunderson-Generator']
    cur.execute("""
        INSERT INTO pub.hit_groups (hit_id, group_id, is_primary)
        SELECT id, %s, true FROM public.master_hits WHERE search_bound = -1
    """, (gid,))
    print(f"  Saunderson-Generator: {cur.rowcount} Hits")

    # -2 → Rathbun-Search
    gid = group_id['Rathbun-Search']
    cur.execute("""
        INSERT INTO pub.hit_groups (hit_id, group_id, is_primary)
        SELECT id, %s, true FROM public.master_hits WHERE search_bound = -2
    """, (gid,))
    print(f"  Rathbun-Search: {cur.rowcount} Hits")

    # -3 → MW-{m}-{n}
    n_mw_inserted = 0
    for m, n in mw_pairs:
        gid = group_id[f'MW-{m}-{n}']
        cur.execute("""
            INSERT INTO pub.hit_groups (hit_id, group_id, is_primary)
            SELECT id, %s, true FROM public.master_hits
            WHERE search_bound = -3 AND m = %s AND n = %s
        """, (gid, m, n))
        n_mw_inserted += cur.rowcount
    print(f"  MW-*-*: {n_mw_inserted} Hits über {len(mw_pairs)} Fasern")
    conn.commit()

    # ------------------------------------------------------------------
    # Schritt 5: Familien-Klassifikation
    # ------------------------------------------------------------------
    print("\nFamilien-Klassifikation läuft...", flush=True)
    cur.execute("SELECT id, x, y, z FROM pub.master_hits")
    n_in_family = {'Saunderson': 0, 'Lenhart': 0, 'Himane-T1': 0, 'Himane-T2': 0, 'Himane-T3': 0}
    sporadic_count = 0
    sporadic_id = group_id['Sporadic']

    family_map = {
        'Saunderson': (saunderson, group_id['Saunderson']),
        'Lenhart':    (lenhart,    group_id['Lenhart']),
        'Himane-T1':  (himane['T1'], group_id['Himane-T1']),
        'Himane-T2':  (himane['T2'], group_id['Himane-T2']),
        'Himane-T3':  (himane['T3'], group_id['Himane-T3']),
    }

    inserts = []  # (hit_id, group_id)
    n_processed = 0
    for hit_id, x, y, z in cur.fetchall():
        br = primitive_brick_from_xyz(x, y, z)
        matched = False
        for fname, (bricks, gid) in family_map.items():
            if br in bricks:
                inserts.append((hit_id, gid))
                n_in_family[fname] += 1
                matched = True
        if not matched:
            inserts.append((hit_id, sporadic_id))
            sporadic_count += 1
        n_processed += 1
        if n_processed % 50000 == 0:
            print(f"  {n_processed} Hits klassifiziert (so far)...", flush=True)

    cur2 = conn.cursor()
    cur2.executemany(
        "INSERT INTO pub.hit_groups (hit_id, group_id, is_primary) VALUES (%s, %s, false) ON CONFLICT DO NOTHING",
        inserts
    )
    conn.commit()
    print("  Familien-Tags zugewiesen:")
    for fname, n in n_in_family.items():
        print(f"    {fname}: {n}")
    print(f"    Sporadic: {sporadic_count}")

    # ------------------------------------------------------------------
    # Schritt 6: Statistik in db_metadata
    # ------------------------------------------------------------------
    cur.execute("SELECT count(*) FROM pub.master_hits")
    total = cur.fetchone()[0]
    cur.execute("SELECT count(*) FROM pub.master_hits WHERE num_blockers IS NOT NULL")
    factored = cur.fetchone()[0]
    cur.execute("INSERT INTO pub.db_metadata (key, value) VALUES ('hit_count', %s) ON CONFLICT (key) DO UPDATE SET value=EXCLUDED.value, updated_at=now()", (str(total),))
    cur.execute("INSERT INTO pub.db_metadata (key, value) VALUES ('factored_count', %s) ON CONFLICT (key) DO UPDATE SET value=EXCLUDED.value, updated_at=now()", (str(factored),))
    cur.execute("INSERT INTO pub.db_metadata (key, value) VALUES ('mw_fibers', %s) ON CONFLICT (key) DO UPDATE SET value=EXCLUDED.value, updated_at=now()", (str(len(mw_pairs)),))
    cur.execute("INSERT INTO pub.db_metadata (key, value) VALUES ('snapshot_date', now()::text) ON CONFLICT (key) DO UPDATE SET value=EXCLUDED.value, updated_at=now()")
    conn.commit()

    print(f"\nFertig in {(time()-t0)/60:.1f} min")
    print(f"  Total Hits in pub: {total}")
    print(f"  Davon faktorisiert: {factored}")
    conn.close()


if __name__ == "__main__":
    main()
