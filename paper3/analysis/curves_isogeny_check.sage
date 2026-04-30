"""
Structure test: for a sample of master tuples, build the five elliptic curves
   E_{m,n}, E'_{m,n}, E_{X,Y}, E_{X,Z}, E_{Y,Z}
and check whether they are isogenous (j-invariant, conductor, isogeny class).

Hypothesis: if e.g. E_{X,Y} and E_{m,n} are always isogenous, then the
"third curve" is just a travesty of the master fibration - no new lever.

Usage:
    sage curves_isogeny_check.sage [N_SAMPLES]
    Default: 10
"""
from sage.all import *
import psycopg
import sys
import signal

pari.allocatemem(int(2e9))

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
TIMEOUT = int(60)


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(int(0))


def quartic_to_E(P_quartic, var_pair=('s', 'w')):
    """Convert y^2 = P(T) (quartic) to EllipticCurve."""
    R = PolynomialRing(QQ, list(var_pair))
    s, w = R.gens()
    T = P_quartic.parent().gen()
    eqn = w**2 - P_quartic.subs({T: s})
    coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
    return EllipticCurve(coeffs)


def build_master_curve(m, n):
    U2 = m*m - n*n
    V2 = 2*m*n
    Rt = PolynomialRing(QQ, 'T')
    T = Rt.gen()
    P = V2**2 * T**4 + (4*U2**2 - 2*V2**2) * T**2 + V2**2
    return quartic_to_E(P)


def build_cuboid_curve(m, n):
    U2 = m*m - n*n
    V2 = 2*m*n
    W2 = m*m + n*n
    Rt = PolynomialRing(QQ, 'T')
    T = Rt.gen()
    P = W2**2 * T**4 + 2*(U2**2 - V2**2) * T**2 + W2**2
    return quartic_to_E(P)


def build_xy_curve(x, y):
    Ax = x*x + y*y
    Bx = (x*x - y*y)**2
    Rt = PolynomialRing(QQ, 'T')
    T = Rt.gen()
    P = T**4 - 2*Ax*T**2 + Bx
    return quartic_to_E(P)


def safe_j(E):
    try:
        return E.j_invariant()
    except Exception:
        return None


def safe_conductor(E):
    try:
        return with_timeout(TIMEOUT, lambda: E.conductor())
    except Exception:
        return None


def isogeny_class(E):
    try:
        return with_timeout(TIMEOUT, lambda: tuple(sorted(EE.j_invariant()
                                                          for EE in E.isogeny_class())))
    except Exception:
        return None


def main():
    conn = psycopg.connect(DB)
    cur = conn.cursor()
    print(f"Loading {N} smallest master tuples...")
    cur.execute("""
        SELECT a, b, m, n, x, y, z
        FROM pub.master_hits
        ORDER BY GREATEST(x, y, z) ASC
        LIMIT %s
    """, (N,))
    rows = cur.fetchall()
    conn.close()
    print(f"Loaded: {len(rows)}\n")

    n_jXY_eq_jE = 0
    n_jXZ_eq_jE = 0
    n_jYZ_eq_jE = 0
    n_jXY_eq_jEp = 0
    n_jXZ_eq_jEp = 0
    n_jYZ_eq_jEp = 0
    n_some_isog = 0
    n_total = 0

    for a, b, m, n, x, y, z in rows:
        a, b, m, n, x, y, z = int(a), int(b), int(m), int(n), int(x), int(y), int(z)
        try:
            E   = build_master_curve(m, n)
            Ep  = build_cuboid_curve(m, n)
            ExY = build_xy_curve(x, y)
            ExZ = build_xy_curve(x, z)
            EyZ = build_xy_curve(y, z)
        except Exception as ex:
            print(f"  ({a},{b},{m},{n}): build_fail: {ex}")
            continue

        n_total += 1
        jE   = safe_j(E)
        jEp  = safe_j(Ep)
        jXY  = safe_j(ExY)
        jXZ  = safe_j(ExZ)
        jYZ  = safe_j(EyZ)

        flags = []
        if jXY == jE:  n_jXY_eq_jE  += 1; flags.append("XY=E")
        if jXZ == jE:  n_jXZ_eq_jE  += 1; flags.append("XZ=E")
        if jYZ == jE:  n_jYZ_eq_jE  += 1; flags.append("YZ=E")
        if jXY == jEp: n_jXY_eq_jEp += 1; flags.append("XY=E'")
        if jXZ == jEp: n_jXZ_eq_jEp += 1; flags.append("XZ=E'")
        if jYZ == jEp: n_jYZ_eq_jEp += 1; flags.append("YZ=E'")

        # Isogeny test (j-comparison only on the fixed curve, isogeny_class
        # gives all class members)
        try:
            iso_E = with_timeout(TIMEOUT, lambda: set(EE.j_invariant() for EE in E.isogeny_class()))
        except Exception:
            iso_E = None

        if iso_E is not None:
            for label, jj in [("XY", jXY), ("XZ", jXZ), ("YZ", jYZ), ("E'", jEp)]:
                if jj is not None and jj in iso_E:
                    flags.append(f"{label}~E")
                    if label != "E'":
                        n_some_isog += 1

        print(f"  ({a},{b},{m},{n}) -> ({x},{y},{z})  flags={flags}")
        print(f"     j(E)={jE}")
        print(f"     j(E')={jEp}")
        print(f"     j(E_XY)={jXY}")
        print(f"     j(E_XZ)={jXZ}")
        print(f"     j(E_YZ)={jYZ}")
        sys.stdout.flush()

    print(f"\n========== Summary ==========")
    print(f"Total checked: {n_total}")
    print(f"j(E_XY) == j(E):   {n_jXY_eq_jE}")
    print(f"j(E_XZ) == j(E):   {n_jXZ_eq_jE}")
    print(f"j(E_YZ) == j(E):   {n_jYZ_eq_jE}")
    print(f"j(E_XY) == j(E'):  {n_jXY_eq_jEp}")
    print(f"j(E_XZ) == j(E'):  {n_jXZ_eq_jEp}")
    print(f"j(E_YZ) == j(E'):  {n_jYZ_eq_jEp}")
    print(f"E_xy isogenous to E (somehow): {n_some_isog}")


if __name__ == "__main__":
    main()
