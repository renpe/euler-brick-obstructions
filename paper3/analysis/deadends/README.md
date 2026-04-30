# Deadends -- documented dead-end approaches

These scripts belong to approaches we pursued and identified as not
viable. They are kept here because their outcomes are interesting for
the paper as "why this seemingly natural approach does not work".

## (g_+, g_-) fibration as an MW lever

- `gpm_subfiber_explore.py` -- distribution of hits across (m, n)
  subfibers within a (g_+, g_-) class.
- `gpm_subfiber_mw_test.sage` -- test whether the (g_+, g_-) layer in
  E_{m, n}(Q) is a subgroup/coset.

**Outcome:** no. The (g_+, g_-) labelling is a representative
convention, not a point invariant. Behind it sits only a 2-element
class mod squares, which is smeared by non-trivial torsion
(|tor| = 8 in the top fiber).

## E_{x, y} construction as a "third" elliptic curve

- `euler_xy_rank.sage` -- rank determination of the x, y curve.
- `euler_xy_cuboid_intersect.sage` -- intersection test with the
  cuboid condition.
- `euler_xy_scale.py` -- 100K-sample brute-force search.

**Outcome:** the curve E_{x, y} is algebraically identical to E_{m, n}
after a permutation of parameters (see `curves_isogeny_check.sage` in
the parent directory: j-invariants match). E_{x, y} therefore yields
no independent Mordell-Weil lever.

## Coleman-Stoll bound too loose for our genus-3 curves

- `jh_coleman_bound.sage` -- attempt to use |H(F_p)| + 2g - 2 or
  Stoll's |H(F_p)| + 2r as a rigorous bound.

**Outcome:** even with primes up to 200, |H(F_p)| >= 8 always (the
8 trivial Q-points reduce generically to 8 distinct F_p-points).
The bound is therefore >= 10 (rk = 1) or >= 12 (rk = 2) -- too loose
to force |H(Q)| = 8.

**What worked instead:** the torsion trick
(`jh_torsion_scale.sage` in the parent directory) -- directly via
|E_q(Q)| = 4 rather than via Stoll.

## Refined torsion trick with bug

- `jh_torsion_refined.sage` -- attempt to cover the |tors| != 4 cases
  by brute-force enumeration on the quartic.

**Outcome:** bug -- double-counting via the (w, +/-V) symmetry plus
points at infinity not captured. The correct approach is to use
Sage's `torsion_subgroup()` directly.

## Sage API test for genus-3 Jacobian rank

- `jacobian_rank_attempt.sage` -- attempt to use
  `J = HyperellipticCurve.jacobian()` for rank determination.

**Outcome:** Sage 10.7 has no `rank()` or `analytic_rank()`
implemented for genus-3 Jacobians. We get further via the explicit
decomposition J(H) ~ E_PQ x E_uV x E_3, not via the Jacobian class
directly.
