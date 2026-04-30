"""
mwrank pass: for fibers where ellrank+kolyvagin returned an ambiguous
[0, k>0] for at least one quotient, try Sage's mwrank backend with
two_descent(second_limit=15). This sometimes resolves the rank exactly,
either confirming rank=0 (proves the fiber) or rank>=1 (rules out the
quotient).

Reads (m,n) pairs from CSV, writes JSONL with mwrank-resolved ranks
(only for the quotients that were originally ambiguous).

Usage:
    sage -python jh_mwrank_pass_worker.py <input.csv> <output.jsonl>
"""
from sage.all import *
import sys
import os
import json
import signal
import warnings
from math import gcd as pygcd
from time import time

warnings.filterwarnings('ignore')
pari.allocatemem(int(2e9))

TIMEOUT_MWRANK = int(os.environ.get('TIMEOUT_MWRANK', '300'))
SECOND_LIMIT = int(os.environ.get('MWRANK_SECOND_LIMIT', '12'))


def with_timeout(seconds, fn):
    def h(sig, frame):
        raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(0)


def quartic_to_E(P_quartic):
    R = PolynomialRing(QQ, ['s', 'w'])
    s, w = R.gens()
    T = P_quartic.parent().gen()
    eqn = w**2 - P_quartic.subs({T: s})
    coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
    return EllipticCurve(coeffs)


def mwrank_resolve(E):
    """Run mwrank's two_descent and return ([rank, rank_bound], certain)."""
    mw = E.mwrank_curve()
    mw.two_descent(second_limit=SECOND_LIMIT, verbose=False)
    return [int(mw.rank()), int(mw.rank_bound())], bool(mw.certain())


def build_quartics(m, n):
    U2 = m*m - n*n
    V2 = 2*m*n
    W2 = m*m + n*n

    Rs = PolynomialRing(QQ, 'S')
    Sv = Rs.gen()
    quartic_PQ = (V2**2 * Sv**2 + (4*U2**2 - 2*V2**2) * Sv + V2**2) * \
                 (W2**2 * Sv**2 + 2*(U2**2 - V2**2) * Sv + W2**2)

    Ru = PolynomialRing(QQ, 'U')
    Uv = Ru.gen()
    quartic_uV = (V2**2 * Uv**2 + 4*(U2**2 - V2**2)) * (W2**2 * Uv**2 - 4*V2**2)

    Rw = PolynomialRing(QQ, 'W')
    Wv = Rw.gen()
    quartic_3 = (V2**2 * Wv**2 + 4*U2**2) * (W2**2 * Wv**2 + 4*U2**2)

    return {'E_PQ': quartic_PQ, 'E_uV': quartic_uV, 'E_3': quartic_3}


def main():
    if len(sys.argv) != 3:
        print("Usage: sage -python jh_mwrank_pass_worker.py <input.csv> <output.jsonl>",
              file=sys.stderr)
        sys.exit(2)
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Input format: each line "m,n,quotient" — only retest the specified quotient
    pairs = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            m, n, q = int(parts[0]), int(parts[1]), parts[2]
            pairs.append((m, n, q))

    def to_py(o):
        if isinstance(o, dict):
            return {to_py(k): to_py(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [to_py(x) for x in o]
        if o is None or isinstance(o, (bool, str)):
            return o
        try: return int(o)
        except (TypeError, ValueError): pass
        try: return float(o)
        except (TypeError, ValueError): pass
        return str(o)

    t0 = time()
    with open(output_file, 'w') as out:
        for i, (m, n, q) in enumerate(pairs):
            result = {'m': m, 'n': n, 'quotient': q}
            try:
                quartics = build_quartics(m, n)
                E = quartic_to_E(quartics[q])
                t_call = time()
                try:
                    (rk, certain) = with_timeout(
                        TIMEOUT_MWRANK,
                        lambda: mwrank_resolve(E))
                    result['mwrank_rank'] = rk[0]
                    result['mwrank_rank_bound'] = rk[1]
                    result['mwrank_certain'] = certain
                    result['mwrank_time'] = float(time() - t_call)
                except TimeoutError:
                    result['error'] = f'timeout after {TIMEOUT_MWRANK}s'
                except Exception as e:
                    result['error'] = f'{type(e).__name__}: {str(e)[:120]}'
            except Exception as e:
                result['error'] = f'setup: {type(e).__name__}: {str(e)[:120]}'
            result['t_elapsed'] = float(time() - t0)
            out.write(json.dumps(to_py(result)) + '\n')
            out.flush()
            if (i + 1) % 5 == 0:
                print(f"  worker progress: {i+1}/{len(pairs)} ({100.0*(i+1)/len(pairs):.1f}%), "
                      f"elapsed {time()-t0:.0f}s",
                      file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
