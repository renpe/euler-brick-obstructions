# Computational scripts for Paper 4

*Odd-exponent blockers and a Mordell-Weil construction of Euler bricks*,
R. Peschmann (in preparation).

These scripts produce and verify the database underlying Paper 4. The
public DB schema (`pub`) and the Mordell-Weil generator are the two
subjects of the paper, and every claim in §3, §5, and §6 is reproducible
from the scripts in this directory.

## Requirements

- **PostgreSQL** ≥ 14
- **Python** ≥ 3.10 with the dependencies in `requirements.txt`
- **SageMath** ≥ 10.5 with PARI/GP ≥ 2.15 (for elliptic-curve generation)
- **YAFU** ≥ 2.13 (for $f_1$ factorisation)
- **CADO-NFS** (optional, for partial residuals ≥ 84 digits)

Configuration is via environment variables; copy the repo-root
`.env.example` to `.env` and adjust:

```
EULER_DB=host=... port=... dbname=... user=... password=...
YAFU_BIN=/path/to/yafu
SAGE_BIN=/path/to/sage
CADO_NFS=/path/to/cado-nfs.py
CADO_VENV_PY=/path/to/cado-nfs.venv/bin/python
```

`pub_db.py` loads `.env` automatically via `python-dotenv`.
Binary variables default to the unqualified names (`yafu`, `sage`,
`cado-nfs.py`, `python3`); set them only if the binaries are not on
`PATH`.

## Structure

```
paper4/
├── README.md                         (this file)
├── _common/
│   └── pub_db.py                     # Helpers: insert_master_hit, ensure_mw_group
├── migration/
│   ├── pub_schema.sql                # DDL for the publication schema 'pub'
│   └── pub_migrate.py                # Migrate public.master_hits → pub.master_hits
├── generation/
│   ├── mw_dispatcher.py              # Subprocess-isolated MW generator
│   ├── mw_fibre_worker.sage          # Single-fibre worker (JSON output)
│   ├── mw_rerun.py                   # Rerun under-covered high-rank fibres
│   ├── search_new_fibres.sage        # Systematic (m,n) rank-scan
│   └── auto_scan_escalate.sh         # Phase-escalation wrapper for the scan
├── factorization/
│   └── pub_factorize.py              # YAFU pipeline targeting pub schema
└── analysis/
    ├── classify_families.py          # Re-classify against Saunderson/Lenhart/Himane
    ├── perfect_check.py              # Verify x²+y²+z² ≠ □ on all hits
    ├── theorem_check.py              # Verify num_blockers > 0 on all factored hits
    └── validate_xyz.py               # Internal consistency: x,y,z,g_scale,x_prim,…
```

## Reproduction workflow

The scripts are designed to run in this order. Each step is idempotent:
restarting after a crash will skip already-completed work.

### Phase 1 — Schema setup

```bash
PGPASSWORD=euler psql -h 192.168.178.63 -U euler -d euler \
    -f migration/pub_schema.sql
```

This creates the `pub` schema with six tables (`master_hits`, `fibers`,
`generator_groups`, `hit_groups`, `f1_factors`, `db_metadata`) and the
view `pub.hit_taxonomy`. Top-level family and provenance groups are
inserted statically.

### Phase 2 — Initial migration from the legacy schema

```bash
.venv-linux/bin/python migration/pub_migrate.py
```

Copies all Master-Hits from `public.master_hits` into `pub.master_hits`,
maps `search_bound` to provenance groups, classifies each hit against
the Saunderson, Lenhart, and Himane families, and creates the 411
`MW-{m}-{n}` provenance subgroups. Takes about 1 minute on the full DB.

### Phase 3 — Mordell-Weil generation (optional, extends the DB)

```bash
.venv-linux/bin/python generation/mw_dispatcher.py 4 4 7
```

Runs the subprocess-isolated MW generator over fibres in `pub.master_hits`
having between 4 and 7 hits. Each fibre is processed by a separate Sage
subprocess with a 60-second hard timeout, so SIGSEGVs in `eclib/mwrank`
do not affect the dispatcher. Output goes to `pub.master_hits` plus
`pub.hit_groups` with provenance `MW-{m}-{n}`.

For deeper exploration:

```bash
nohup ./generation/auto_scan_escalate.sh > /tmp/scan.log 2>&1 &
```

This runs three phases of `search_new_fibres.sage` (m≤200, m≤500,
m≤1000) with four parallel workers, then automatically replays
`mw_rerun.py` on every (m,n) found with `analytic_rank ≥ 2`.

### Phase 4 — Factorisation of $f_1$

```bash
nohup .venv-linux/bin/python factorization/pub_factorize.py 6 \
    > /tmp/pub_factorize.log 2>&1 &
```

Six parallel YAFU workers factorise $f_1$ for every hit with
`num_blockers IS NULL`. Trial division → sympy.factorint (3 s timeout)
→ YAFU (10 min timeout). Updates `pub.master_hits.num_blockers` and
populates `pub.f1_factors`.

### Phase 5 — Verification

```bash
.venv-linux/bin/python analysis/theorem_check.py
.venv-linux/bin/python analysis/perfect_check.py
```

Both must return zero counterexamples, in keeping with §4 and §5.4 of
the paper.

## Convention: provenance and family tags

Every hit in `pub.master_hits` carries (via `pub.hit_groups`)

- exactly one **primary provenance** tag from
  `Exhaustive-Bound-{200,500,2000,2300}`, `Rathbun-Search`,
  `Saunderson-Generator`, or `MW-{m}-{n}`;
- zero or more **family** tags from `Saunderson`, `Lenhart`,
  `Himane-T1`, `Himane-T2`, `Himane-T3`, with the fallback `Sporadic`
  if no family matches.

Use `pub.hit_taxonomy` (a view) for convenient retrieval:

```sql
SELECT id, a, b, m, n, families, provenance
FROM pub.hit_taxonomy
WHERE 'Saunderson' = ANY(families) AND num_blockers >= 0;
```

## License

The scripts are distributed under the licence stated in the top-level
[`LICENSE`](../../LICENSE).

## Author

René Peschmann
