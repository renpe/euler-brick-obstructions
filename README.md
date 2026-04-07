# Computational scripts for "Reductions and obstructions for the perfect Euler brick problem"

These scripts support the computational claims in the paper by R. Peschmann (2026).

## Requirements

- **PARI/GP** >= 2.15 (for `ellrank`, Simon 2-descent)
- **SageMath** >= 9.0 (for elliptic curve construction, point searches)
- **Python** >= 3.8 with `sympy` (for `product_square_check.py`)

## Scripts and paper claims

### Core verifications

| Script | Paper reference | What it verifies |
|--------|---------------|------------------|
| `product_square_check.py` | Section 7, item (1) | f1*f2 is never a perfect square for a,b,m,n up to 1000 |
| `pari_check_rank0.gp` | Remark 4.3 | 42 (resp. 54) certified rank-0 specialisations of E_A (resp. E'_A) |
| `sage_Eprime_root_diffs.sage` | Proposition 6.4 | Root differences of E'_A; c-primes are harmless for 2-Selmer |
| `sage_modular_search.sage` | Section 7, item (2) | Modular search: 175,418 lattice points, 0 candidates |
| `sage_blocker_primes.sage` | Section 7, item (3) | Blocker prime statistics (88.4% split, 11.6% p=2, 0% inert) |

### Supporting analyses

| Script | Paper reference | What it computes |
|--------|---------------|------------------|
| `sage_kummer_character.sage` | Theorem 5.4 | Kummer character chi_f on torsion of E_A |
| `sage_genus3.sage` | Section 4 | Genus-3 curve C_A and its elliptic quotients |
| `sage_EA_rank_evidence.sage` | Remark 4.3 | Rank-0 specialisations of E_A (Silverman evidence) |
| `sage_Eprime_rank_evidence.sage` | Remark 4.3 | Rank-0 specialisations of E'_A (Silverman evidence) |

## Running

```bash
# PARI/GP rank certification
gp -q pari_check_rank0.gp

# Product square check (Python, no Sage needed)
python3 product_square_check.py

# Sage scripts
sage sage_Eprime_root_diffs.sage
sage sage_modular_search.sage
sage sage_blocker_primes.sage
sage sage_kummer_character.sage
sage sage_genus3.sage
sage sage_EA_rank_evidence.sage
sage sage_Eprime_rank_evidence.sage
```

## Author

René Peschmann
