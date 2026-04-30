"""
Classify pub.master_hits against the classical families (Saunderson,
Lenhart, Himane T1/T2/T3). Idempotent: only hits without any
family tag in pub.hit_groups are processed.

Run after every mw_dispatcher.py / mw_rerun.py to assign family tags
to newly-generated bricks.

Usage:
    python3 classify_families.py
"""
import os
import sys
from math import gcd, isqrt
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


# ---------------- Family enumeration -----------------

def is_square(n):
    if n <= 0: return n == 0
    r = isqrt(n)
    return r * r == n


def primitive_brick(e1, e2, e3):
    g = gcd(gcd(abs(e1), abs(e2)), abs(e3))
    if g == 0: return None
    return tuple(sorted([abs(e1)//g, abs(e2)//g, abs(e3)//g]))


def gen_saunderson(max_g=500):
    bricks = set()
    for g in range(2, max_g+1):
        for h in range(1, g):
            if gcd(g, h) != 1 or (g-h) % 2 == 0: continue
            u, v, w = g*g-h*h, 2*g*h, g*g+h*h
            uu, vv = min(u, v), max(u, v)
            if 3*uu*uu <= vv*vv: continue
            x = uu*(3*vv*vv - uu*uu)
            y = vv*(3*uu*uu - vv*vv)
            z = 4*uu*vv*w
            if x <= 0 or y <= 0: continue
            br = primitive_brick(x, y, z)
            if br: bricks.add(br)
    return bricks


def gen_lenhart(max_w=300):
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


def gen_pyth_triples(max_gen=80):
    return [(g*g-h*h, 2*g*h, g*g+h*h)
            for g in range(2, max_gen+1)
            for h in range(1, g)
            if gcd(g, h) == 1 and (g-h) % 2 == 1]


def gen_himane(triples, max_scale=15):
    out = {'T1': set(), 'T2': set(), 'T3': set()}
    for u1, v1, w1 in triples:
        for u2, v2, w2 in triples:
            for k in range(1, max_scale+1):
                u0, v0 = u1*k, v2*k
                a, b, c = u0*u2, v0*u1, v0*v1
                if a > 0 and b > 0 and c > 0 and is_square(a*a+c*c):
                    br = primitive_brick(a, b, c)
                    if br: out['T1'].add(br)
                u0, v0 = u1*k, u2*k
                a, b, c = u0*u2, v0*v1, u0*v2
                if a > 0 and b > 0 and c > 0 and is_square(b*b+c*c):
                    br = primitive_brick(a, b, c)
                    if br: out['T2'].add(br)
                u0, v0 = v1*k, v2*k
                a, b, c = u0*u2, v0*u1, v0*v1
                if a > 0 and b > 0 and c > 0 and is_square(a*a+b*b):
                    br = primitive_brick(a, b, c)
                    if br: out['T3'].add(br)
    return out


def main():
    print("Generating family bricks...", flush=True)
    saunderson = gen_saunderson()
    lenhart = gen_lenhart()
    himane = gen_himane(gen_pyth_triples())
    print(f"  Saunderson: {len(saunderson)}, Lenhart: {len(lenhart)}, "
          f"Himane T1/T2/T3: {len(himane['T1'])}/{len(himane['T2'])}/{len(himane['T3'])}",
          flush=True)

    conn = pub_db.connect(autocommit=False)
    cur = conn.cursor()

    family_map = {
        'Saunderson': (saunderson, pub_db.get_group_id(cur, 'Saunderson')),
        'Lenhart':    (lenhart,    pub_db.get_group_id(cur, 'Lenhart')),
        'Himane-T1':  (himane['T1'], pub_db.get_group_id(cur, 'Himane-T1')),
        'Himane-T2':  (himane['T2'], pub_db.get_group_id(cur, 'Himane-T2')),
        'Himane-T3':  (himane['T3'], pub_db.get_group_id(cur, 'Himane-T3')),
    }
    sporadic_id = pub_db.get_group_id(cur, 'Sporadic')

    # Find hits with no family tag
    cur.execute("""
        SELECT mh.id, mh.x, mh.y, mh.z FROM pub.master_hits mh
        WHERE NOT EXISTS (
            SELECT 1 FROM pub.hit_groups hg
            JOIN pub.generator_groups g ON g.id = hg.group_id
            WHERE hg.hit_id = mh.id AND g.category = 'family'
        )
    """)
    rows = cur.fetchall()
    print(f"  {len(rows)} hits without family tag", flush=True)
    if not rows:
        return

    inserts = []
    counts = {k: 0 for k in family_map}
    counts['Sporadic'] = 0
    t0 = time()

    for hit_id, x, y, z in rows:
        br = primitive_brick(int(x), int(y), int(z))
        matched = False
        for fname, (bricks, gid) in family_map.items():
            if br in bricks:
                inserts.append((hit_id, gid))
                counts[fname] += 1
                matched = True
        if not matched:
            inserts.append((hit_id, sporadic_id))
            counts['Sporadic'] += 1

    cur2 = conn.cursor()
    cur2.executemany(
        """INSERT INTO pub.hit_groups (hit_id, group_id, is_primary)
           VALUES (%s, %s, false) ON CONFLICT DO NOTHING""",
        inserts,
    )
    conn.commit()
    print(f"\nClassified in {time()-t0:.1f}s:")
    for k, v in counts.items():
        print(f"  {k}: {v}")
    conn.close()


if __name__ == "__main__":
    main()
