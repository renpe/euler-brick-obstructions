"""
Worker for the full torsion-intersection scan over the three quotients
E_PQ, E_uV, E_3 of H_{m,n}.

For each (m, n), compute rank bounds and torsion size of all three
elliptic factors. A fiber is rigorously proven (Conjecture B excluded
on the entire fiber) when one of the three quotients satisfies:

    rank bounds [r_low, r_high] = [0, 0]
    AND
    2 * |tors| - r_q == 8

where r_q is the number of rational ramification points of H -> E_q
(determined geometrically; r_3 = 0, r_uV = 4, r_PQ = 4).

This gives:
    E_3:  proven when |tors| = 4
    E_uV: proven when |tors| = 6
    E_PQ: proven when |tors| = 6

The eight trivial points on H always exist (lower bound 8), so 2|tors|-r_q = 8
is tight: it forces |H(Q)| = 8.

Usage:
    sage jh_torsion_full_worker.sage <input.csv> <output.jsonl>
"""
from sage.all import *
import sys
import json
import signal
from math import gcd as pygcd
from time import time

pari.allocatemem(int(2e9))

import os
TIMEOUT_RANK = int(os.environ.get('TIMEOUT_RANK', '180'))
TIMEOUT_TORS = 30
TIMEOUT_ARUB = int(os.environ.get('TIMEOUT_ARUB', '300'))
ELLRANK_EFFORT = int(os.environ.get('ELLRANK_EFFORT', '2'))


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


def safe_rank(E):
    """Returns (rank_bounds, L_ratio_nonzero, is_semistable, resolution_method).

    rank_bounds is [r_lo, r_hi] from PARI ellrank, possibly tightened to [0, 0]
    when L(E, 1) != 0 forces rank=0 via the modularity theorem and Kolyvagin's
    theorem.

    L_ratio_nonzero is 1 if L(E, 1) / Omega_E != 0, 0 if zero, None if not computed.
    is_semistable is True iff the conductor is squarefree (Manin constant <= 2
    by Mazur, hence L_ratio() is provably correct).

    resolution_method is one of:
      'ellrank'              -- ellrank gave a sharp bound [0,0] on its own.
      'modular_symbol'       -- ambiguous ellrank, but L_ratio != 0 and the curve
                                 is semistable, so L(E,1) != 0 unconditionally;
                                 rank=0 follows by Kolyvagin's theorem.
      'modular_symbol_cond'  -- ambiguous ellrank, L_ratio != 0 but the curve is
                                 not semistable (squarefree conductor fails);
                                 the L_ratio computation depends on the Manin
                                 constant being <= 2.
      'ambiguous'            -- ambiguous ellrank, L_ratio = 0 or unavailable.
      'failed'               -- ellrank itself failed.
    """
    try:
        result = with_timeout(
            TIMEOUT_RANK,
            lambda: E.pari_curve().ellrank(ELLRANK_EFFORT),
        )
        rk_lo, rk_hi = int(result[0]), int(result[1])
    except Exception:
        return [None, None], None, None, 'failed'

    L_ratio_nonzero = None
    is_semistable = None
    method = 'ellrank' if rk_lo == rk_hi else 'ambiguous'

    if rk_lo == 0 and rk_hi > 0:
        try:
            ratio = with_timeout(
                TIMEOUT_ARUB,
                lambda: E.lseries().L_ratio(),
            )
            L_ratio_nonzero = 0 if ratio == 0 else 1
            cond = E.conductor()
            is_semistable = bool(cond == cond.squarefree_part())
            if L_ratio_nonzero == 1:
                rk_hi = 0  # rank=0 by Kolyvagin
                method = 'modular_symbol' if is_semistable \
                         else 'modular_symbol_cond'
        except Exception:
            pass

    return [rk_lo, rk_hi], L_ratio_nonzero, is_semistable, method


def safe_torsion_size(E):
    try:
        return int(with_timeout(TIMEOUT_TORS,
                                lambda: E.torsion_subgroup().order()))
    except Exception:
        return None


# Number of rational ramification points of H -> E_q on E_q (i.e. the count
# of branch points that are themselves rational, contributing 1 instead of 2
# preimages each to the cardinality bound).
RAMIFICATION_POINTS = {'E_3': 0, 'E_uV': 4, 'E_PQ': 4}

# For "proven via this quotient" with the naive bound 2|tors| - r_q = 8,
# we need |tors| = (8 + r_q) / 2.
TORSION_TARGET = {'E_3': 4, 'E_uV': 6, 'E_PQ': 6}


def count_H_preimages_via_E_uV(m, n):
    """Refined count of |H(Q)| using the quotient E_uV via PARI hyperellratpoints.

    For each rational point on the E_uV quartic
        V^2 = (V2^2 u^2 + 4(U2^2 - V2^2)) (W2^2 u^2 - 4 V2^2)
    we count H-preimages:
      * u = +-2 (ramification of sigma_2): 1 H-preimage each
      * affine non-ram (u, V) with V != 0: 2 H-preimages iff u^2-4 is a
        non-negative rational square (lift via t^2 - u t + 1 = 0)
      * each of the two points at infinity (rational since leading is V2^2 W2^2):
        2 H-preimages (sigma_2 swaps t=0 with t=infty in the H-orbit)

    Returns (h_count, n_affine_pts, leading_is_square) where:
      h_count       : total H-preimage count
      n_affine_pts  : number of distinct affine rational points found on the quartic
      leading_is_square : whether the two points at infinity are rational

    Caller must verify that n_affine_pts + (2 if leading_is_square else 0)
    equals |E_uV(Q)| = |tors_uV| before trusting h_count as a sharp upper
    bound. If enumeration is incomplete, h_count is only a lower bound.
    """
    U2 = m*m - n*n
    V2 = 2*m*n
    W2 = m*m + n*n

    Ru = PolynomialRing(QQ, 'U')
    U = Ru.gen()
    quartic = (V2**2 * U**2 + 4*(U2**2 - V2**2)) * (W2**2 * U**2 - 4*V2**2)

    # PARI's hyperellratpoints with a high height bound.
    pari_q = pari(quartic)
    raw = pari.hyperellratpoints(pari_q, 100000)
    affine_pts = set()
    for p in raw:
        u0 = QQ(p[0])
        v0 = QQ(p[1])
        affine_pts.add((u0, v0))

    h_count = 0
    for (u0, V0) in affine_pts:
        if u0 == 2 or u0 == -2:
            h_count += 1
        else:
            disc = u0**2 - 4
            if disc >= 0 and disc.is_square():
                h_count += 2

    leading = quartic.list()[-1]
    leading_is_sq = leading > 0 and leading.is_square()
    if leading_is_sq:
        h_count += 4  # two points at infinity, each with 2 H-preimages

    return int(h_count), int(len(affine_pts)), bool(leading_is_sq)


def analyze_fiber(m, n):
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

    result = {'m': m, 'n': n}
    quotients = [
        ('E_PQ', quartic_PQ),
        ('E_uV', quartic_uV),
        ('E_3',  quartic_3),
    ]

    proven_via = None
    proven_method = None
    for label, quartic in quotients:
        try:
            E = quartic_to_E(quartic)
        except Exception:
            result[f'rk_{label}'] = [None, None]
            result[f'L_ratio_nonzero_{label}'] = None
            result[f'semistable_{label}'] = None
            result[f'rk_method_{label}'] = 'failed'
            result[f'tors_{label}'] = None
            continue
        rk, L_ratio_nonzero, is_semistable, method = safe_rank(E)
        result[f'rk_{label}'] = rk
        result[f'L_ratio_nonzero_{label}'] = L_ratio_nonzero
        result[f'semistable_{label}'] = is_semistable
        result[f'rk_method_{label}'] = method
        if rk[0] == 0 and rk[1] == 0:
            tors = safe_torsion_size(E)
            result[f'tors_{label}'] = tors
            if tors == TORSION_TARGET[label] and proven_via is None:
                proven_via = label
                proven_method = 'naive'
        else:
            result[f'tors_{label}'] = None

    # Refinement pass: if not yet proven and E_uV has rk=[0,0] with |tors|=8,
    # try the explicit lift count.
    if proven_via is None and result.get('rk_E_uV') == [0, 0] \
            and result.get('tors_E_uV') == 8:
        try:
            h_count, n_affine, leading_sq = with_timeout(
                120, lambda: count_H_preimages_via_E_uV(m, n))
            result['refined_E_uV_h_count'] = h_count
            result['refined_E_uV_n_affine'] = n_affine
            result['refined_E_uV_leading_sq'] = leading_sq
            n_inf = 2 if leading_sq else 0
            total_pts = n_affine + n_inf
            result['refined_E_uV_total_pts'] = total_pts
            # Only accept refinement when enumeration is provably complete
            # (matches the torsion order) AND lift count equals trivial 8.
            if total_pts == 8 and h_count == 8:
                proven_via = 'E_uV'
                proven_method = 'refined'
        except Exception:
            result['refined_E_uV_h_count'] = None
            result['refined_E_uV_n_affine'] = None
            result['refined_E_uV_leading_sq'] = None
            result['refined_E_uV_total_pts'] = None

    result['proven_via'] = proven_via
    result['proven_method'] = proven_method
    result['proven'] = proven_via is not None
    return result


def main():
    if len(sys.argv) != 3:
        print("Usage: sage jh_torsion_full_worker.sage <input.csv> <output.jsonl>",
              file=sys.stderr)
        sys.exit(2)
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    pairs = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            m, n = line.split(',')
            pairs.append((int(m), int(n)))

    def to_py(obj):
        if isinstance(obj, dict):
            return {to_py(k): to_py(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [to_py(x) for x in obj]
        if obj is None or isinstance(obj, (bool, str)):
            return obj
        try:
            return int(obj)
        except (TypeError, ValueError):
            pass
        try:
            return float(obj)
        except (TypeError, ValueError):
            pass
        return str(obj)

    t0 = time()
    with open(output_file, 'w') as out:
        for i, (m, n) in enumerate(pairs):
            try:
                result = analyze_fiber(m, n)
            except Exception as e:
                result = {'m': m, 'n': n, 'error': str(e)}
            result['t_elapsed'] = time() - t0
            result = to_py(result)
            out.write(json.dumps(result) + '\n')
            out.flush()
            if (i + 1) % 10 == 0:
                pct = float((i + 1)) / float(len(pairs)) * 100.0
                el = float(time() - t0)
                print(f"  worker progress: {i+1}/{len(pairs)} "
                      f"({pct:.1f}%), elapsed {el:.0f}s",
                      file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
