"""
Systematische Suche nach neuen hochrangigen (m,n)-Fasern.

Iteriert (m, n) mit gcd=1, (m-n) ungerade über m ∈ [2, M_MAX].
Berechnet analytic_rank für E_{m,n}. Bei Rang ≥ 2 in CSV speichern.

Skip:
  - Familien-Fasern (Saunderson/Himane/Lenhart) aus family_mn_skip.json
  - (88, 7)
  - alle (m, n), die schon in master_hits vorkommen

Output: rank_scan.csv (m, n, analytic_rank)
"""
from sage.all import *
from sage.all import alarm, cancel_alarm
import psycopg
import json
import csv
from math import gcd as pygcd
from time import time
import os
import sys

pari.allocatemem(2 * 10**9)

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"
# Argumente: M_MAX [WORKER_ID] [NUM_WORKERS] [M_MIN]
M_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
WORKER_ID = int(sys.argv[2]) if len(sys.argv) > 2 else 0
NUM_WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 1
M_MIN = int(sys.argv[4]) if len(sys.argv) > 4 else 2
OUT_FILE = f"rank_scan_w{WORKER_ID}of{NUM_WORKERS}_max{M_MAX}.csv"
RANK_TIMEOUT = 30  # Sekunden harter Timeout für analytic_rank


def compute_rank(m_p, n_p):
    """Returns analytic_rank or None on failure/timeout."""
    U2 = m_p*m_p - n_p*n_p
    V2 = 2*m_p*n_p
    A = V2*V2
    B = 4*U2*U2 - 2*V2*V2
    C = V2*V2
    R.<x, y> = QQ[]
    eqn = y**2 - (A*x**4 + B*x**2 + C)
    try:
        alarm(RANK_TIMEOUT)
        coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
        E = EllipticCurve(coeffs)
        ar = int(E.analytic_rank())
        cancel_alarm()
        return ar
    except AlarmInterrupt:
        cancel_alarm()
        return None
    except Exception:
        try: cancel_alarm()
        except Exception: pass
        return None


def main():
    print(f"M_MAX = {M_MAX}", flush=True)

    # Skip-Liste direkt aus pub.hit_groups lesen
    print("Lade Familien-Fasern aus pub.hit_groups...", flush=True)
    conn = psycopg.connect(DB)
    cur = conn.cursor()
    cur.execute("""
        SELECT DISTINCT mh.m, mh.n
        FROM pub.master_hits mh
        JOIN pub.hit_groups hg ON hg.hit_id = mh.id
        JOIN pub.generator_groups g ON g.id = hg.group_id
        WHERE g.short_name IN ('Saunderson', 'Lenhart', 'Himane-T1', 'Himane-T2', 'Himane-T3')
    """)
    skip_mn = set((int(m), int(n)) for m, n in cur.fetchall())
    skip_mn.add((88, 7))  # explicit, paper §3 example fibre
    print(f"  {len(skip_mn)} family/skip (m,n)")

    # DB-bekannte (m, n)
    print("Lade DB-(m,n) aus pub.master_hits...", flush=True)
    cur.execute("SELECT DISTINCT m, n FROM pub.master_hits")
    db_mn = set((int(m), int(n)) for m, n in cur.fetchall())
    conn.close()
    print(f"  {len(db_mn)} DB-(m,n)")

    # Faser-Kandidaten: jeder Worker nimmt m % NUM_WORKERS == WORKER_ID
    print(f"Worker {WORKER_ID}/{NUM_WORKERS}, m ∈ [{M_MIN}, {M_MAX}], "
          f"nimmt m mit m %% {NUM_WORKERS} == {WORKER_ID}", flush=True)
    todo = []
    for m_p in range(M_MIN, M_MAX + 1):
        if m_p % NUM_WORKERS != WORKER_ID: continue
        for n_p in range(1, m_p):
            if pygcd(m_p, n_p) != 1: continue
            if (m_p - n_p) % 2 != 1: continue
            if (m_p, n_p) in skip_mn: continue
            if (m_p, n_p) in db_mn: continue
            todo.append((m_p, n_p))
    print(f"  {len(todo)} (m,n)-Kandidaten zu testen")

    # Resume: bereits im CSV vorhandene Paare überspringen
    seen = set()
    if os.path.exists(OUT_FILE):
        with open(OUT_FILE) as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2 and row[0] != 'm':
                    seen.add((int(row[0]), int(row[1])))
        print(f"  Resume: {len(seen)} bereits gescannt")

    todo = [(m, n) for m, n in todo if (m, n) not in seen]
    print(f"  {len(todo)} verbleibend\n", flush=True)

    if not todo:
        print("Nichts zu tun.")
        return

    # CSV öffnen, falls neu Header schreiben
    is_new = not os.path.exists(OUT_FILE) or os.path.getsize(OUT_FILE) == 0
    fout = open(OUT_FILE, 'a', buffering=1)
    if is_new:
        fout.write("m,n,analytic_rank\n")

    t0 = time()
    n_high = 0   # rank >= 2
    n_done = 0
    last_report = t0

    for m_p, n_p in todo:
        ar = compute_rank(m_p, n_p)
        n_done += 1
        if ar is not None:
            fout.write(f"{m_p},{n_p},{ar}\n")
            if ar >= 2:
                n_high += 1
                print(f"  ★ ({m_p},{n_p}): analytic_rank = {ar}", flush=True)
        now = time()
        if now - last_report > 10:
            elapsed = now - t0
            rate = n_done / elapsed if elapsed > 0 else 0
            eta = (len(todo) - n_done) / rate if rate > 0 else 0
            print(f"  {n_done}/{len(todo)} | rank≥2: {n_high} | {rate:.2f}/s | "
                  f"ETA {eta/60:.1f} min", flush=True)
            last_report = now

    fout.close()
    print(f"\nFERTIG in {(time()-t0)/60:.1f} min")
    print(f"  Hochrangige (≥2): {n_high}")


if __name__ == "__main__":
    main()
