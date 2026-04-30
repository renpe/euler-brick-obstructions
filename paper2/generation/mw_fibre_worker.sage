"""
Single-fiber worker. Arguments: m_p n_p MAX_SCALAR [DESCENT_LIMIT]
Output (stdout, JSON):
  {"status": "ok"|"fail"|"nogens", "gens_count": k, "candidates": [[a,b], ...]}

Called as a subprocess from the dispatcher, which sets a hard timeout.
DB inserts are done by the dispatcher, not the worker.
"""
from sage.all import *
import json
import sys
from math import gcd as pygcd
import psycopg

pari.allocatemem(2 * 10**9)

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"


def main():
    m_p = int(sys.argv[1])
    n_p = int(sys.argv[2])
    max_scalar = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    descent_limit = int(sys.argv[4]) if len(sys.argv) > 4 else 15

    U2 = m_p*m_p - n_p*n_p
    V2 = 2*m_p*n_p
    A = V2*V2
    B = 4*U2*U2 - 2*V2*V2
    C = V2*V2
    gamma = V2

    R.<x, y> = QQ[]
    eqn = y**2 - (A*x**4 + B*x**2 + C)
    try:
        coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
        E = EllipticCurve(coeffs)
    except Exception as ex:
        print(json.dumps({"status": "fail", "reason": f"curve: {ex}"}))
        return

    # Try E.gens() without timeout (subprocess timeout applies)
    gens = None
    try:
        gens = list(E.gens(proof=False, descent_second_limit=descent_limit))
    except RuntimeError as e:
        import re
        msg = str(e)
        mm = re.search(r'generators found=\((.*?)\)\.', msg)
        if mm:
            try:
                pts_str = '[' + mm.group(1) + ']'
                pts = eval(pts_str)
                gens_list = []
                for p in pts:
                    try:
                        P = E.point(p)
                        gens_list.append(P)
                    except Exception:
                        pass
                gens = gens_list if gens_list else None
            except Exception:
                gens = None
    except Exception:
        gens = None

    # Fallback: lift DB points
    if gens is None or len(gens) == 0:
        try:
            conn = psycopg.connect(DB)
            cur = conn.cursor()
            cur.execute("SELECT a, b FROM pub.master_hits WHERE m=%s AND n=%s", (m_p, n_p))
            db_ab = [(int(a), int(b)) for a, b in cur.fetchall()]
            conn.close()
        except Exception:
            db_ab = []
        R_QQ.<T> = QQ[]
        fT = A*T**4 + B*T**2 + C
        pts = []
        for a, b in db_ab:
            t = QQ(a) / QQ(b)
            s_sq = fT.subs(T=t)
            if not s_sq.is_square(): continue
            s_pos = sqrt(s_sq)
            for s_sig in (s_pos, -s_pos):
                X_val = 2*gamma*(s_sig + gamma) / t**2
                try:
                    lifts = E.lift_x(X_val, all=True)
                    if lifts:
                        pts.append(lifts[0])
                        break
                except Exception:
                    pass
        if len(pts) > 4: pts = pts[:4]
        gens = pts
        if not gens:
            print(json.dumps({"status": "nogens", "gens_count": 0, "candidates": []}))
            return

    # Scalar multiplication
    tors = list(E.torsion_subgroup())
    seen_X = set()
    candidates = []
    from itertools import product as iproduct
    ranges = [range(-max_scalar, max_scalar+1)] * len(gens)
    for cvec in iproduct(*ranges):
        if all(c == 0 for c in cvec): continue
        P = E(0)
        for g, c in zip(gens, cvec):
            P = P + c*g
        for T in tors:
            Q = P + T
            if Q == E(0): continue
            X_val = Q.xy()[0]
            if X_val in seen_X: continue
            seen_X.add(X_val)
            denom = X_val*X_val - 4*A*gamma*gamma
            if denom == 0: continue
            t_sq = 4*gamma*gamma*(X_val + B) / denom
            if t_sq < 0 or not t_sq.is_square(): continue
            t = sqrt(t_sq)
            if t == 0: continue
            a = abs(t.numerator())
            b = abs(t.denominator())
            if a <= b: continue
            if (a - b) % 2 == 0: continue
            if pygcd(int(a), int(b)) != 1: continue
            U1 = a*a - b*b
            V1_ = 2*a*b
            mv = (V1_*U2)**2 + (U1*V2)**2
            if not mv.is_square(): continue
            candidates.append([int(a), int(b)])

    print(json.dumps({
        "status": "ok",
        "gens_count": len(gens),
        "candidates": candidates,
    }))


if __name__ == "__main__":
    main()
