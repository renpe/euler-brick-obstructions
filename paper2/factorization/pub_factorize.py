"""
Factorise f1_prim = x_prim² + y_prim² + z_prim² (primitive Brick-Raumdiagonale²)
für jeden Master-Hit. Verwendet die algebraische Identität
    f1_prim = L_prim · R_prim
mit L_prim = (C/g_-)² + (D/g_-)², R_prim = (A/g_+)² + (B/g_+)²
(primitive Summen zweier Quadrate, halbe Stellenzahl).

Pipeline pro Faktor:
  Stage 1: trial division up to 10^6 (Python)
  Stage 2: sympy.factorint with 3-second timeout
  Stage 3: YAFU subprocess with 5-minute hard timeout

Targets every Master-Hit with num_blockers IS NULL.
Schreibt in pub.prim_brick_factors und updated num_blockers.

Usage:
    python3 pub_factorize.py [workers=6]
"""
import os
import re
import signal
import subprocess
import sys
from multiprocessing import Pool
from time import time

import psycopg
from sympy import factorint

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db

YAFU_BIN = pub_db.YAFU_BIN
# Strategy: many workers, few threads per YAFU. Most f1 are factored
# entirely by trial-division + sympy, so a worker rarely actually
# launches YAFU. With 30 workers × 1 thread we get 30 cores busy in
# the trial/sympy phase. When YAFU is actually needed (large residual),
# 30 single-threaded YAFU calls run in parallel — slower per call but
# higher overall throughput.
YAFU_THREADS = 1
SYMPY_TIMEOUT = 3
YAFU_TIMEOUT = 300       # 5 Minuten — schneller Durchlauf, partials später
TRIAL_LIMIT = 10**6


class FactorTimeout(Exception):
    pass


def _alarm(sig, frame):
    raise FactorTimeout()


def trial_div(n, limit):
    factors = {}
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2
    p = 3
    while p <= limit and p*p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 2
    return factors, n


def sympy_factor(n, seconds):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(seconds)
    try:
        return dict(factorint(n))
    finally:
        signal.alarm(0)


_FACTOR_RE = re.compile(r'^P\d+\s*=\s*(\d+)$', re.MULTILINE)


def yafu_factor(n, threads=YAFU_THREADS, timeout=YAFU_TIMEOUT):
    inp = f"factor({n})\n"
    tmpdir = f"/tmp/yafu_{os.getpid()}_{n % 100000}"
    os.makedirs(tmpdir, exist_ok=True)
    try:
        result = subprocess.run(
            [YAFU_BIN, "-threads", str(threads)],
            input=inp, capture_output=True, text=True,
            cwd=tmpdir, timeout=timeout,
        )
    finally:
        try:
            for f in os.listdir(tmpdir):
                os.remove(os.path.join(tmpdir, f))
            os.rmdir(tmpdir)
        except Exception:
            pass

    out = result.stdout
    primes = [int(m) for m in _FACTOR_RE.findall(out)]
    factors = {}
    rem = n
    for p in primes:
        e = 0
        while rem % p == 0:
            rem //= p
            e += 1
        if e > 0:
            factors[p] = e
    if rem != 1:
        return None
    return factors


def full_factor(n):
    """Return (factors_dict, status) with status in {'full', 'partial'}."""
    factors, rest = trial_div(n, TRIAL_LIMIT)
    if rest == 1:
        return factors, 'full'
    try:
        sub = sympy_factor(rest, SYMPY_TIMEOUT)
        for p, e in sub.items():
            factors[p] = factors.get(p, 0) + e
        return factors, 'full'
    except FactorTimeout:
        pass
    except Exception:
        pass
    try:
        sub = yafu_factor(rest)
        if sub is not None:
            for p, e in sub.items():
                factors[p] = factors.get(p, 0) + e
            return factors, 'full'
    except subprocess.TimeoutExpired:
        pass
    except Exception:
        pass
    return factors, 'partial'


def process_hit(row):
    """Faktorisiert f1_prim = x_prim² + y_prim² + z_prim² via algebraischer
    Identität f1_prim = L_prim · R_prim, mit
        L_prim = ((am-bn)/g_-)² + ((an-bm)/g_-)²    (g_- = gcd(C,D))
        R_prim = ((am+bn)/g_+)² + ((an+bm)/g_+)²    (g_+ = gcd(A,B))
    Jeder Faktor wird einzeln durch die trial-div → sympy → YAFU Pipeline
    gejagt; die Faktor-Dicts werden multiplikativ vereinigt. Status 'full'
    nur wenn beide Teilfaktorisierungen 'full' sind.
    """
    hit_id, a, b, m, n, x_prim, y_prim, z_prim = row
    hit_id = int(hit_id)
    a = int(a); b = int(b); m = int(m); n = int(n)
    x_prim = int(x_prim); y_prim = int(y_prim); z_prim = int(z_prim)
    try:
        f1_prim = x_prim*x_prim + y_prim*y_prim + z_prim*z_prim
        # algebraische Aufteilung mit g_+, g_-
        from math import gcd
        A = a*m + b*n; B = a*n + b*m
        C = a*m - b*n; D = a*n - b*m
        gp = gcd(A, B)
        gm = gcd(abs(C), abs(D))
        # L_prim und R_prim als primitive Summen zweier Quadrate
        Cp = C // gm; Dp = D // gm
        Ap = A // gp; Bp = B // gp
        L_prim = Cp*Cp + Dp*Dp
        R_prim = Ap*Ap + Bp*Bp
        # Sanity-Check: L_prim · R_prim == f1_prim (algebraische Identität)
        assert L_prim * R_prim == f1_prim, f"identity violated for hit {hit_id}"
        fL, sL = full_factor(L_prim)
        fR, sR = full_factor(R_prim)
        factors = dict(fL)
        for p, e in fR.items():
            factors[p] = factors.get(p, 0) + e
        status = 'full' if (sL == 'full' and sR == 'full') else 'partial'
    except Exception as e:
        return (hit_id, None, f"ERROR: {e}")
    return (hit_id, (x_prim, y_prim, z_prim, factors, status), None)


def flush(conn, results):
    cur = conn.cursor()
    full = partial = err = 0
    for hit_id, payload, _ in results:
        if payload is None:
            err += 1
            continue
        x_prim, y_prim, z_prim, factors, status = payload
        if status == 'full':
            pub_db.store_prim_brick_factors(cur, hit_id, x_prim, y_prim, z_prim, factors)
            full += 1
        else:
            # Partial: store the trial-division factors but leave num_blockers NULL
            for p, e in factors.items():
                cur.execute(
                    """INSERT INTO pub.prim_brick_factors
                           (x_prim, y_prim, z_prim, prime, exponent, is_blocker, prime_mod4)
                       VALUES (%s, %s, %s, %s, %s, %s, %s)
                       ON CONFLICT (x_prim, y_prim, z_prim, prime) DO NOTHING""",
                    (x_prim, y_prim, z_prim, int(p), int(e), e % 2 == 1, int(p) % 4),
                )
            partial += 1
    conn.commit()
    return full, partial, err


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    conn = pub_db.connect(autocommit=False)
    cur = conn.cursor()
    # KORREKTUR (2026-04-27): Faktorisiere ALLE Tupel, nicht nur g_scale=1.
    # Grund: die Standard-Parametrisierung (a,b,m,n) erzeugt nicht alle
    # primitiven Euler-Bricks direkt — z.B. Saunderson-Bricks erscheinen nur
    # als skalierte Versionen (g_scale > 1). Für Vollständigkeit der
    # Konjektur B müssen wir auch nicht-primitive Hits faktorisieren, denn
    # deren f₁ = g_scale² · f₁_prim enthält die Information über den
    # primitiven Wurzel-Brick (oft außerhalb der Standard-Familie).
    # Faktorisiere f1_prim = x_prim²+y_prim²+z_prim² statt f1.
    # Das ist das Wesentliche für Konjektur B (g_scale² ist immer Quadrat
    # und trägt nichts zur Blocker-Struktur bei). Faktorisiert wird via
    # L_prim · R_prim = f1_prim mit primitiven Summen zweier Quadrate.
    cur.execute("""
        SELECT id, a, b, m, n, x_prim, y_prim, z_prim FROM pub.master_hits
        WHERE num_blockers IS NULL
        ORDER BY (x_prim*x_prim + y_prim*y_prim + z_prim*z_prim)
    """)
    tasks = cur.fetchall()
    total = len(tasks)
    print(f"Unfactored hits: {total} | YAFU workers={workers} ({YAFU_THREADS} threads each)",
          flush=True)
    if total == 0:
        return

    t0 = time()
    n_full = n_partial = n_err = 0
    batch = []
    done = 0
    last_report = time()

    with Pool(processes=workers) as pool:
        for result in pool.imap_unordered(process_hit, tasks, chunksize=1):
            batch.append(result)
            done += 1
            if len(batch) >= 50:
                f, p, e = flush(conn, batch)
                n_full += f; n_partial += p; n_err += e
                batch = []
                now = time()
                if now - last_report > 10:
                    elapsed = now - t0
                    rate = done / elapsed if elapsed > 0 else 0
                    eta = (total - done) / rate if rate > 0 else 0
                    print(f"  {done}/{total} | full={n_full} partial={n_partial} err={n_err} "
                          f"| {rate:.1f}/s | ETA {eta/60:.1f} min", flush=True)
                    last_report = now
    if batch:
        f, p, e = flush(conn, batch)
        n_full += f; n_partial += p; n_err += e

    print(f"\nFertig in {(time()-t0)/60:.1f} min")
    print(f"  Full={n_full}, Partial={n_partial}, Err={n_err}")
    conn.close()


if __name__ == "__main__":
    main()
