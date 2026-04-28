"""
Rigorous rank computation for the elliptic fibres E_{m,n}.

Phase A: every (m, n) for which pub.master_hits already has at least
         one record (i.e. all 411 MW-active fibres plus the classical
         search hits — about 55,000 total).
Phase B: every (m, n) with gcd(m, n) = 1, m - n odd, m <= M_MAX,
         that is NOT in Phase A. Default M_MAX = 500.

For each fibre we call PARI's `ellrank` (2-descent + 4-descent), which
returns rigorous lower and upper bounds on rk E_{m,n}(Q). Output goes
to scripts/paper4/data/pari_ranks.csv with columns:

    m, n, conductor_digits, rank_lower, rank_upper, sha2_dim, n_gens

Usage:
    sage pari_rank_scan.sage [PHASE] [M_MAX]
    PHASE in {A, B, both}; default both
    M_MAX default 500

Each fibre is given a 60s hard timeout via signal.SIGALRM.
"""
from sage.all import *
import psycopg
import csv
import os
import sys
import signal
from math import gcd as pygcd
from time import time

pari.allocatemem(2 * 10**9)

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"
PHASE = sys.argv[1] if len(sys.argv) > 1 else "both"
M_MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 500
OUT = "scripts/paper4/data/pari_ranks.csv"
TIMEOUT = 60   # seconds per fibre


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(seconds)
    try:
        return fn()
    finally:
        signal.alarm(0)


def compute_rank(m, n):
    """Returns dict with rank info or None on failure/timeout."""
    U2 = m*m - n*n
    V2 = 2*m*n
    A = V2*V2; B = 4*U2*U2 - 2*V2*V2; C = V2*V2
    R = PolynomialRing(QQ, ['x', 'y'])
    x, y = R.gens()
    eqn = y**2 - (A*x**4 + B*x**2 + C)
    try:
        coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
        E = EllipticCurve(coeffs)
        cond = E.conductor()
    except Exception as ex:
        return {'error': f'curve: {ex}'}
    try:
        result = with_timeout(TIMEOUT, lambda: E.pari_curve().ellrank())
    except TimeoutError:
        return {'error': 'timeout', 'conductor_digits': len(str(cond))}
    except Exception as ex:
        return {'error': f'ellrank: {ex}', 'conductor_digits': len(str(cond))}
    rl = int(result[0])
    ru = int(result[1])
    sha2 = int(result[2]) if len(result) > 2 else None
    gens = result[3] if len(result) > 3 else []
    return {
        'conductor_digits': len(str(cond)),
        'rank_lower': rl,
        'rank_upper': ru,
        'sha2_dim': sha2,
        'n_gens': len(gens) if isinstance(gens, list) else None,
    }


def load_phase_a():
    conn = psycopg.connect(DB)
    cur = conn.cursor()
    cur.execute("SELECT DISTINCT m, n FROM pub.master_hits")
    pairs = [(int(m), int(n)) for m, n in cur.fetchall()]
    pairs.sort(key=lambda mn: mn[0] + mn[1])
    conn.close()
    return pairs


def load_phase_b(m_max, exclude):
    todo = []
    for m in range(2, m_max + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1: continue
            if (m - n) % 2 != 1: continue
            if (m, n) in exclude: continue
            todo.append((m, n))
    return todo


def load_already_done(out_path):
    done = set()
    if os.path.exists(out_path):
        with open(out_path) as f:
            for row in csv.reader(f):
                if len(row) >= 2 and row[0] != 'm':
                    done.add((int(row[0]), int(row[1])))
    return done


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    is_new = not os.path.exists(OUT) or os.path.getsize(OUT) == 0
    fout = open(OUT, 'a', buffering=1)
    if is_new:
        fout.write("m,n,conductor_digits,rank_lower,rank_upper,sha2_dim,n_gens,error\n")

    done = load_already_done(OUT)
    print(f"Resume: {len(done)} fibres already done.", flush=True)

    pairs = []
    if PHASE in ("A", "both"):
        a = load_phase_a()
        pairs += [(m, n, "A") for (m, n) in a if (m, n) not in done]
        print(f"Phase A: {len(a)} DB-fibres ({len(pairs)} not yet done)", flush=True)
    if PHASE in ("B", "both"):
        a_set = set((m, n) for (m, n, _) in pairs) | done
        b = load_phase_b(M_MAX, a_set)
        pairs += [(m, n, "B") for (m, n) in b]
        print(f"Phase B: {len(b)} new fibres up to m={M_MAX}", flush=True)

    print(f"Total to scan: {len(pairs)}", flush=True)

    n_done = 0; n_err = 0; n_high = 0
    t0 = time()
    last = t0
    for (m, n, phase) in pairs:
        info = compute_rank(m, n) or {}
        cond = info.get('conductor_digits', '')
        rl = info.get('rank_lower', '')
        ru = info.get('rank_upper', '')
        sha2 = info.get('sha2_dim', '')
        ngens = info.get('n_gens', '')
        err = info.get('error', '')
        fout.write(f"{m},{n},{cond},{rl},{ru},{sha2},{ngens},{err}\n")
        n_done += 1
        if err:
            n_err += 1
        elif isinstance(rl, int) and rl >= 4:
            n_high += 1
            print(f"  ★ ({m},{n}) phase {phase}: rank={rl}-{ru}", flush=True)
        now = time()
        if now - last > 30:
            elapsed = now - t0
            rate = n_done / elapsed if elapsed > 0 else 0
            eta = (len(pairs) - n_done) / rate if rate > 0 else 0
            print(f"  {n_done}/{len(pairs)} | rank>=4: {n_high} | err: {n_err} "
                  f"| {rate:.2f}/s | ETA {eta/60:.1f} min", flush=True)
            last = now

    fout.close()
    print(f"\nFinished in {(time()-t0)/60:.1f} min")
    print(f"  scanned: {n_done}, errors: {n_err}, rank>=4: {n_high}")


if __name__ == "__main__":
    main()
