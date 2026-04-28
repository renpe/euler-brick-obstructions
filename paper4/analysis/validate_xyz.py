"""
Internal-consistency check on every record in pub.master_hits.

For every Master-Hit (a, b, m, n) we verify:

  (1)  x = U_1 U_2,   y = V_1 U_2,   z = U_1 V_2
       where U_1 = a^2 - b^2, V_1 = 2 a b,
             U_2 = m^2 - n^2, V_2 = 2 m n;

  (2)  g_scale = gcd(x, y, z);

  (3)  x_prim * g_scale = x,
       y_prim * g_scale = y,
       z_prim * g_scale = z.

Each check is a finite, exact integer identity. A successful run
proves that every record in pub.master_hits is internally consistent
(in particular, that g_scale = 1 is equivalent to the brick being
primitive in the integer-gcd sense).

Usage:
    python3 validate_xyz.py
"""
import os
import sys
from math import gcd
from time import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "_common"))
import pub_db


def main():
    conn = pub_db.connect()
    cur = conn.cursor("validate_xyz")  # server-side cursor for streaming
    cur.itersize = 10000
    cur.execute(
        "SELECT id, a, b, m, n, x, y, z, g_scale, x_prim, y_prim, z_prim "
        "FROM pub.master_hits"
    )
    n = 0
    err_xyz = 0
    err_g = 0
    err_prim = 0
    examples = []
    t0 = time()
    for row in cur:
        hid, a, b, m, nn, x, y, z, gs, xp, yp, zp = (int(v) for v in row)
        U1 = a*a - b*b
        V1 = 2*a*b
        U2 = m*m - nn*nn
        V2 = 2*m*nn
        ex = U1*U2
        ey = V1*U2
        ez = U1*V2
        if (x, y, z) != (ex, ey, ez):
            err_xyz += 1
            if len(examples) < 5:
                examples.append(("xyz mismatch", hid, (a, b, m, nn), (x, y, z), (ex, ey, ez)))
        g_expected = gcd(gcd(abs(x), abs(y)), abs(z))
        if g_expected != gs:
            err_g += 1
            if len(examples) < 5:
                examples.append(("g_scale mismatch", hid, (a, b, m, nn), gs, g_expected))
        if xp * gs != x or yp * gs != y or zp * gs != z:
            err_prim += 1
            if len(examples) < 5:
                examples.append(("primitive mismatch", hid, (a, b, m, nn), (xp, yp, zp), gs))
        n += 1
    elapsed = time() - t0
    print(f"Validated {n} Master-Hits in {elapsed:.1f} s.")
    print(f"  edges (x,y,z) match U_i V_i U_2 V_2:        {err_xyz} errors")
    print(f"  g_scale = gcd(x, y, z):                     {err_g} errors")
    print(f"  (x_prim, y_prim, z_prim) * g_scale = (x,y,z): {err_prim} errors")
    if examples:
        print("\nFirst few violations:")
        for e in examples:
            print(f"  {e}")
        conn.close()
        return 1
    else:
        print("\nAll records internally consistent.")
        conn.close()
        return 0


if __name__ == "__main__":
    sys.exit(main())
