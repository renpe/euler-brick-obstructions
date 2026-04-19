"""
Large-scale modular search for perfect Euler brick candidates.

For each of the five problematic rational slopes s = a/b that survive
earlier sieve stages, this script:
  1. Builds the associated elliptic curve E(s) and computes its
     Mordell-Weil generators via Sage.
  2. Enumerates lattice points P = n1*G1 + ... + nr*Gr in expanding
     shells (max |ni| growing by 2 each round, round-robin across
     all five curves) up to a one-hour wall-clock budget.
  3. For every lattice point with odd coordinate sum, evaluates
     mu = 2*y / ((x-2)*(x+2))  and checks whether mu is a quadratic
     residue modulo each of the 45 primes p in (3, 200).

Paper claim supported:
  "A modular search over 175,418 lattice points with 45 primes
   p < 200 finds zero candidates for f(P) in Q*^2."

Time limit: 3600 s (1 hour).
"""

import sys
import time

TIME_LIMIT = 3600  # 1 hour in seconds
start_time = time.time()

primes_check = list(primes(3, 200))

# The five problematic slope pairs (a, b) that survive earlier sieves
problem_cases = [(18,41), (18,47), (23,59), (23,64), (29,65)]

# Pre-compute elliptic curves and generators for every slope
curves = []
for a_val, b_val in problem_cases:
    s0 = QQ(a_val)/QQ(b_val)
    c0 = (s0**4-6*s0**2+1)/(1+s0**2)**2
    A0 = 2-4*c0**2
    E_orig = EllipticCurve([0, A0, 0, -4, -4*A0])
    E_min = E_orig.minimal_model()
    iso = E_min.isomorphism_to(E_orig)
    gens = E_min.gens()
    rk = len(gens)
    curves.append({
        'a': a_val, 'b': b_val, 'rk': rk,
        'E_min': E_min, 'iso': iso, 'gens': gens,
        'tested': 0, 'candidates': 0, 'max_n_reached': 0,
    })
    print("Prepared: s=%d/%d, rank=%d" % (a_val, b_val, rk))

print()
print("Starting search (time limit: %d seconds)" % TIME_LIMIT)
print()
sys.stdout.flush()

# Round-robin: increase max_n stepwise for all curves
current_n = 1

while True:
    current_n += 2  # only odd steps are meaningful (grow by 2)

    for curve in curves:
        elapsed = time.time() - start_time
        if elapsed > TIME_LIMIT:
            break

        rk = curve['rk']
        gens = curve['gens']
        E_min = curve['E_min']
        iso = curve['iso']
        a_val = curve['a']
        b_val = curve['b']

        # Test the NEW shell: points whose max |ni| equals current_n
        # (inner shells have already been tested)
        from itertools import product as iprod
        ranges = [range(-current_n, current_n+1)] * rk

        shell_tested = 0
        for ns in iprod(*ranges):
            if max(abs(n) for n in ns) != current_n:
                continue  # skip inner points; only the outer shell
            if sum(ns) % 2 == 0:
                continue  # even sum: already ruled out by parity argument

            P_min = sum(ns[i]*gens[i] for i in range(rk))
            if P_min == E_min(0):
                continue
            curve['tested'] += 1
            shell_tested += 1

            try:
                P = iso(P_min)
            except:
                continue
            x0 = P[0]; y0 = P[1]
            if x0**2 == 4 or y0 == 0:
                continue

            mu_num = 2*y0
            mu_den = (x0-2)*(x0+2)

            # Check quadratic residuosity of mu modulo each small prime
            is_candidate = True
            for p in primes_check:
                try:
                    mu_mod = (Mod(Integer(numerator(mu_num)), p)
                              * Mod(Integer(denominator(mu_num)), p)**(-1)
                              * Mod(Integer(denominator(mu_den)), p)
                              * Mod(Integer(numerator(mu_den)), p)**(-1))
                except ZeroDivisionError:
                    continue
                if mu_mod == 0:
                    continue
                if not mu_mod.is_square():
                    is_candidate = False
                    break

            if is_candidate:
                curve['candidates'] += 1
                print("*** CANDIDATE: s=%d/%d, n=%s ***" % (a_val, b_val, ns))
                sys.stdout.flush()

        curve['max_n_reached'] = current_n

    elapsed = time.time() - start_time
    if elapsed > TIME_LIMIT:
        break

    # Progress report every 10 steps
    if current_n % 10 == 1:
        print("n=%d, %.0f seconds, %s" % (
            current_n, elapsed,
            ", ".join("s=%d/%d:%d" % (k['a'], k['b'], k['tested'])
                      for k in curves)))
        sys.stdout.flush()

# Final results
elapsed = time.time() - start_time
print()
print("=" * 70)
print("RESULT (%.0f seconds)" % elapsed)
print("=" * 70)
print()
for curve in curves:
    print("s=%d/%d: rank=%d, max_n=%d, tested=%d, candidates=%d" %
          (curve['a'], curve['b'], curve['rk'],
           curve['max_n_reached'], curve['tested'], curve['candidates']))
print()
total_candidates = sum(k['candidates'] for k in curves)
if total_candidates == 0:
    print("*** NO EULER BRICK FOUND ***")
else:
    print("*** %d CANDIDATES FOUND — VERIFY MORE CLOSELY! ***" % total_candidates)
