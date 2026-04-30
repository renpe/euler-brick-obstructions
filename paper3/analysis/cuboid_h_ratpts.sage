"""
Search for rational points on the hyperelliptic genus-3 curve
   H_{m,n}: v^2 = P(t^2) * Q(t^2)
for a sample (m,n).

Each point (t_0, v_0) gives v_0^2 = P(t_0^2) * Q(t_0^2). For both P AND Q
to individually be squares, we need additional conditions - for a perfect
cuboid this would automatically be the case.

For each found point we check:
  - Is P(t_0^2) itself a square? (master tuple)
  - Is Q(t_0^2) itself a square? (cuboid condition)
  - If both: * perfect cuboid candidate (mod trivial t in {0, +-1, +-infty}).

PARI's hyperellratpoints finds all (t, v) with height <= B efficiently.

Usage:
    sage cuboid_h_ratpts.sage [M_MAX] [B]
    Default: 20 1000
"""
from sage.all import *
import sys
import signal
from math import gcd as pygcd
from time import time

pari.allocatemem(int(2e9))

M_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
B = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
TIMEOUT = int(60)


def with_timeout(seconds, fn):
    def h(sig, frame): raise TimeoutError()
    signal.signal(signal.SIGALRM, h)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(int(0))


def main():
    print(f"M_MAX={M_MAX}, height bound B={B}\n")
    print(f"{'(m,n)':>10} {'#pts':>5} {'*cuboid':>8} {'time':>6}")

    cuboid_findings = []
    n_total = 0
    n_with_pts = 0

    for m in range(2, M_MAX + 1):
        for n in range(1, m):
            if pygcd(m, n) != 1: continue
            if (m - n) % 2 != 1: continue
            n_total += 1
            U2 = m*m - n*n
            V2 = 2*m*n
            W2 = m*m + n*n

            Rt = PolynomialRing(QQ, 'T')
            T = Rt.gen()
            P = V2**2 * T**4 + (4*U2**2 - 2*V2**2) * T**2 + V2**2
            Q = W2**2 * T**4 + 2*(U2**2 - V2**2) * T**2 + W2**2
            f = P * Q  # polynomial in T, degree 8

            t0 = time()
            try:
                # PARI hyperellratpoints
                pts = with_timeout(
                    TIMEOUT,
                    lambda: pari(f).hyperellratpoints(B),
                )
                pts = [(QQ(p[0]), QQ(p[1])) for p in pts]
            except Exception as ex:
                print(f"{f'({m},{n})':>10}   FAIL: {ex}")
                continue
            elapsed = time() - t0

            n_cub = 0
            cubs_here = []
            for (tval, vval) in pts:
                # Trivial: t=0 (gives P(0)=V2^2, Q(0)=W2^2, both squares, but
                # corresponds to (a,b)=(0,1) - degenerate brick z=0)
                if tval == 0:
                    continue
                P_val = P.subs({T: tval})
                Q_val = Q.subs({T: tval})
                p_sq = P_val.is_square()
                q_sq = Q_val.is_square()
                if p_sq and q_sq:
                    # Trivial when t = +-1: then a=b, U1=0, X=0 -> degenerate
                    if abs(tval) == 1:
                        continue
                    n_cub += 1
                    cubs_here.append((tval, vval, P_val, Q_val))

            if pts:
                n_with_pts += 1
            print(f"{f'({m},{n})':>10} {len(pts):>5} {n_cub:>8} {elapsed:>6.1f}")
            sys.stdout.flush()

            if cubs_here:
                cuboid_findings.append((m, n, cubs_here))
                for tval, vval, P_val, Q_val in cubs_here:
                    print(f"    *** (m,n)=({m},{n}): t={tval}, "
                          f"P={P_val}=({sqrt(P_val)})^2, "
                          f"Q={Q_val}=({sqrt(Q_val)})^2")

    print(f"\n========== Summary ==========")
    print(f"Total fibers checked: {n_total}")
    print(f"With rational H-points: {n_with_pts}")
    print(f"Perfect cuboid candidates (non-trivial): "
          f"{sum(len(c) for _,_,c in cuboid_findings)}")
    if cuboid_findings:
        print("\n* ALARM *")
        for m, n, cs in cuboid_findings:
            print(f"  (m,n)=({m},{n}):")
            for tval, vval, P_val, Q_val in cs:
                print(f"    t={tval}, v={vval}")


if __name__ == "__main__":
    main()
