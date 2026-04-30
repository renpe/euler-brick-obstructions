# Computational scripts for Paper 3

*A torsion-intersection refinement: rigorous proof of the perfect-cuboid conjecture on an explicit family of master-tuple fibers*,
R. Peschmann (in preparation).

## Relation to Paper 1

Paper 3 is a **companion paper** to Peschmann (2026), arXiv:2604.09328
(*"Quartic reductions and elliptic obstructions for perfect Euler
bricks"*, Paper 1). Paper 1 establishes the genus-3 reduction
`C_A: w^2 = lambda^8 + A lambda^4 + 1`, the Klein-4 involutions with
three elliptic quotients `E_A x E_A' x E_A''`, and obstructions via
the Kummer character and 2-descent on `E_A`. Paper 3 builds on this
and contributes two specific new results:

1. **Reduction theorem (Theorem 4):** every primitive Euler brick
   descends from a master tuple. This strengthens Paper 1's
   Lemma 3.1 (forward implication only).

2. **Torsion-intersection argument on E_A''** (= E_3 in our notation):
   if rk(E_A'') = 0 with |tors| = 4, then H(Q) = {8 trivial points},
   so no perfect cuboid exists on the (m, n) fiber. Paper 1 only
   uses specialisations of E_A and E_A'.

Concrete output: **117 (m, n) fibers** at M_MAX = 50 on which
Conjecture B is rigorously proven (disjoint from Paper 1's 42 + 54
specialisations).

## Main results

### Theorem 1 (g_+ . g_- structure theorem)
For every master tuple: g_+ . g_- = gcd(U_1, U_2) = g_scale and
gcd(g_+, g_-) = 1.

### Theorem 2 (polynomial factorisation)
f_1 = L . R with L, R primitive sums of two squares.

### Theorem 4 (reduction to master tuples)
Every primitive Euler brick (X, Y, Z) arises from a master tuple via
the standard construction, scaled by g = gcd(U_1, U_2).

### Main theorem (Theorem 4.4)
For (m, n) with rk(E_3(Q)) = 0 and |E_3(Q)_tors| = 4:
|H_{m, n}(Q)| = 8 (all trivial). No perfect cuboid exists on this
fiber.

## Status

- Theorems 1, 2, 4: rigorously proven (Paper 3).
- Genus-3 + 3-quotient decomposition: cited from Paper 1 (Sect. 3, 4).
- Main theorem on 117 explicit fibers at M_MAX = 50.
- ~61 fibers with |tors| != 4: refinement possible.
- ~340 fibers with total rank >= 3: quadratic Chabauty required.
- Infinitely many (m, n) with max(m, n) > 50: open research problem.

See `CONJECTURE_B_ROADMAP.md`.

## Requirements

- Same as Paper 2: PostgreSQL, Python 3.10+, Sage 10.7+, PARI/GP 2.15+
- `pub_db.py` is shared via a symlink from `paper2/_common/`.

## Structure (public, in `scripts/paper3/`)

```
scripts/paper3/
├── README.md                            (this file)
├── _common/
│   └── pub_db.py → ../../paper2/_common/pub_db.py  (symlink)
├── data/
│   └── proven_fibers.csv                (117 fibers as CSV)
└── analysis/
    ├── 15 main-result scripts
    └── deadends/
        ├── README.md                    (why these approaches fail)
        └── 8 dead-end scripts
```

The manuscript material (LaTeX source, proof markdowns, compiled PDF)
lives entirely outside this repository and is published independently
of the code. This README only describes the publicly deployable
reproduction scripts.

## Reproducing the main results

### 117 fibers, rigorous (main theorem)

```bash
sage scripts/paper3/analysis/jh_torsion_scale.sage 50
```

Output: 518 fibers checked, 117 rigorously proven, 61 inconclusive
(|tors| != 4).

### Genus-3 form verification (Theorem 3 / recap of Paper 1)

```bash
sage scripts/paper3/analysis/cuboid_genus_check.sage
sage scripts/paper3/analysis/jh_three_quotients.sage 50
```

### Empirical height search (B = 10^6)

```bash
sage scripts/paper3/analysis/jh_chabauty_finalize.sage 1000000
```

11 fibers x 30 seconds ~= 5 min total; all find only the 6 trivial
affine points.

## Cross-references

- **Paper 1** (arXiv:2604.09328): genus-3 setup, Jacobian
  decomposition, Kummer character, 2-descent.
- **Paper 2** (in preparation): Master-Hit database, MW construction
  of new bricks, blocker analysis.
- Paper 3 uses `pub_db.py` and the `pub` schema from Paper 2.

## Data releases (off-repo)

The compact script outputs (`data/proven_fibers.csv`,
`data/scan_full_M100.jsonl`) are committed to this repository. Larger
data -- in particular the PostgreSQL dump of the `pub` schema, which
this paper consumes via Paper 2 -- is too large to ship here and is
hosted off-repo at

**<https://renepeschmann.de/research>**

with the same `zenodo/`-style layout used for Paper 2.
