"""
Structure test: is the (g_+ = 9, g_- = 1) layer inside the elliptic curve
E_{2368,1207} a subgroup / a coset?

Procedure:
  1. Load all 733 (a,b) hits from pub.master_hits with (m,n)=(2368,1207),
     g_+=9, g_-=1.
  2. Lift them to points on E.
  3. For a sample (P_i, P_j) compute Q = P_i + P_j on E.
  4. Lift Q back to (a, b) and compute (g_+, g_-) for (a, b, m, n).
  5. Collect the distribution - is it concentrated on one fiber?

Results:
  - "closed" -> the (9,1) layer is a subgroup.
  - "always (g_+', g_-') != (9,1) but constant" -> coset.
  - "scatters across many fibers" -> no group structure.

Usage:
    sage gpm_subfiber_mw_test.sage [m] [n] [g_plus] [g_minus] [n_pairs]
    Default: 2368 1207 9 1 200
"""
from sage.all import *
import psycopg
import sys
import os
from collections import Counter
from math import gcd as pygcd

pari.allocatemem(2 * 10**9)

DB = "host=192.168.178.63 port=5432 dbname=euler user=euler password=euler"

m_p = int(sys.argv[1]) if len(sys.argv) > 1 else 2368
n_p = int(sys.argv[2]) if len(sys.argv) > 2 else 1207
target_gp = int(sys.argv[3]) if len(sys.argv) > 3 else 9
target_gm = int(sys.argv[4]) if len(sys.argv) > 4 else 1
n_pairs = int(sys.argv[5]) if len(sys.argv) > 5 else 200

print(f"Test: (m,n)=({m_p},{n_p}), target fiber (g_+,g_-)=({target_gp},{target_gm}),"
      f" n_pairs={n_pairs}")

# --- Build the elliptic curve ----------------------------------------
U2 = m_p*m_p - n_p*n_p
V2 = 2*m_p*n_p
A = V2*V2
B = 4*U2*U2 - 2*V2*V2
C = V2*V2
gamma = V2

R = PolynomialRing(QQ, ['x', 'y'])
x, y = R.gens()
eqn = y**2 - (A*x**4 + B*x**2 + C)
coeffs = [ZZ(c) for c in pari(eqn).ellfromeqn()]
E = EllipticCurve(coeffs)
print(f"E: {E}")

# --- Forward map: (a,b) -> Point on E --------------------------------
RT = PolynomialRing(QQ, 'T')
T = RT.gen()
fT = A*T**4 + B*T**2 + C


def ab_to_point(a, b):
    """Lift (a,b) -> point on E. Returns None on failure."""
    t = QQ(a) / QQ(b)
    s_sq = fT.subs(T=t)
    if not s_sq.is_square():
        return None
    s_pos = sqrt(s_sq)
    # try both signs
    for s_sig in (s_pos, -s_pos):
        X_val = 2*gamma*(s_sig + gamma) / t**2
        try:
            lifts = E.lift_x(X_val, all=True)
            if lifts:
                return lifts[0]
        except Exception:
            pass
    return None


# --- Backward map: Point on E -> (a,b) -------------------------------
def point_to_ab_one(P):
    """Given point on E, recover candidate (a,b). Returns None if not in (a,b)-form."""
    if P == E(0):
        return None
    X_val = P.xy()[0]
    denom = X_val*X_val - 4*A*gamma*gamma
    if denom == 0:
        return None
    t_sq = 4*gamma*gamma*(X_val + B) / denom
    if t_sq < 0 or not t_sq.is_square():
        return None
    t_val = sqrt(t_sq)
    if t_val == 0:
        return None
    a = abs(t_val.numerator())
    b = abs(t_val.denominator())
    return (int(a), int(b))


# Precompute: all torsion points of the curve (small, usually <= 4 points)
TORSION = list(E.torsion_subgroup())
print(f"|E_torsion(Q)| = {len(TORSION)}")


def point_to_ab_all_sigs(P, m, n):
    """Collects all (g_+, g_-) signatures reachable through torsion translates
    of P. This is the full information modulo the trivial ambiguity of the
    (a,b) representation."""
    sigs = set()
    for T in TORSION:
        Q = P + T
        ab = point_to_ab_one(Q)
        if ab is None:
            continue
        a, b = ab
        if a == 0 or b == 0:
            continue
        sig = gpm_of(a, b, m, n)
        if sig is not None:
            sigs.add(sig)
    return sigs


def point_to_ab_canonical(P, m, n):
    """From all torsion translates of P, pick the (a,b) representative
    with minimal naive height max(|a|, |b|). Returns (a,b,sig) or None."""
    best = None
    best_height = None
    for T in TORSION:
        Q = P + T
        ab = point_to_ab_one(Q)
        if ab is None:
            continue
        a, b = ab
        if a == 0 or b == 0:
            continue
        h = max(int(a), int(b))
        if best_height is None or h < best_height:
            best_height = h
            best = (int(a), int(b), gpm_of(int(a), int(b), m, n))
    return best


def gpm_of(a, b, m, n):
    """Compute (g+, g-) for a master tuple. None on degenerate."""
    if a == 0 or b == 0:
        return None
    A = a*m + b*n
    B = a*n + b*m
    C = a*m - b*n
    D = a*n - b*m
    gp = pygcd(int(A), int(B))
    gm = pygcd(int(abs(C)), int(abs(D)))
    return (gp, gm)


# --- Load the 733 hits -----------------------------------------------
conn = psycopg.connect(DB)
cur = conn.cursor()
cur.execute(
    "SELECT a, b FROM pub.master_hits WHERE m=%s AND n=%s",
    (m_p, n_p),
)
all_ab = [(int(a), int(b)) for a, b in cur.fetchall()]
conn.close()

target_ab = [(a, b) for (a, b) in all_ab if gpm_of(a, b, m_p, n_p) == (target_gp, target_gm)]
print(f"Hits in subfiber: {len(target_ab)} of {len(all_ab)} total in (m,n)")

# --- Lift to points --------------------------------------------------
points = []
fail = 0
for a, b in target_ab:
    P = ab_to_point(a, b)
    if P is None:
        fail += 1
        continue
    points.append((P, (a, b)))
print(f"Lifted {len(points)} points (failed {fail}).")

if len(points) < 2:
    print("Not enough points to test addition.")
    sys.exit(0)

# --- Pairwise sum test ----------------------------------------------
import random
random.seed(int(42))
sample_idx = list(range(len(points)))
random.shuffle(sample_idx)

# --- Sanity check: do torsion translates of the INPUT points lead them
# out of their fiber? If so, (g_+, g_-) is not a property of the point
# itself but only of a representative choice.
print("\n=== Sanity: canonical signature (smallest naive height) of the INPUT points ===")
input_canonical_sigs = Counter()
input_match_target = 0
for P, ab in points:
    res = point_to_ab_canonical(P, m_p, n_p)
    if res is None:
        continue
    a, b, sig = res
    input_canonical_sigs[sig] += 1
    if sig == (target_gp, target_gm):
        input_match_target += 1
print(f"Canonical signatures among the {len(points)} inputs (all with original sig "
      f"({target_gp},{target_gm})):")
for sig, c in input_canonical_sigs.most_common(10):
    marker = "  <-TARGET" if sig == (target_gp, target_gm) else ""
    print(f"  {sig}: {c}{marker}")
print(f"-> Inputs whose canonical choice is again ({target_gp},{target_gm}): "
      f"{input_match_target}/{len(points)}")

# --- Sum test: canonical signature per sum --------------------------
sum_canon = Counter()
n_tested = 0
n_no_recovery = 0
for k in range(min(n_pairs, len(points)*(len(points)-1)//2)):
    i, j = random.sample(sample_idx, 2)
    P1, _ = points[i]
    P2, _ = points[j]
    Q = P1 + P2
    n_tested += 1
    res = point_to_ab_canonical(Q, m_p, n_p)
    if res is None:
        n_no_recovery += 1
        continue
    a, b, sig = res
    sum_canon[sig] += 1

print(f"\n=== Sum test ({n_tested} pairs, canonical signature) ===")
print(f"Recovery errors: {n_no_recovery}")
for sig, c in sum_canon.most_common(20):
    marker = "  <-TARGET" if sig == (target_gp, target_gm) else ""
    print(f"  {sig}: {c}{marker}")

# --- Difference test ------------------------------------------------
diff_canon = Counter()
for k in range(min(n_pairs, len(points)*(len(points)-1)//2)):
    i, j = random.sample(sample_idx, 2)
    P1, _ = points[i]
    P2, _ = points[j]
    Q = P1 - P2
    res = point_to_ab_canonical(Q, m_p, n_p)
    if res is None:
        continue
    a, b, sig = res
    diff_canon[sig] += 1

print(f"\n=== Difference test (canonical signature) ===")
for sig, c in diff_canon.most_common(20):
    print(f"  {sig}: {c}")

# --- 2P (doubling) - does (9,1) close on itself? ---------------------
print(f"\n=== Doubling test: 2*P_i for the first {min(100, len(points))} inputs ===")
double_canon = Counter()
for P, _ in points[:100]:
    Q = P + P
    res = point_to_ab_canonical(Q, m_p, n_p)
    if res is None:
        continue
    a, b, sig = res
    double_canon[sig] += 1
for sig, c in double_canon.most_common(10):
    print(f"  {sig}: {c}")
