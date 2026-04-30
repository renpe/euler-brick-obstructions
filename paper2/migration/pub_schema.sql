-- ============================================================
-- Publication schema for the Master-Hit database
-- ============================================================
-- Creates a dedicated schema 'pub' in the DB 'euler' that is lean,
-- documentable, and publication-ready. The old 'public.master_hits'
-- remains untouched as the working DB.
--
-- Apply:
--   PGPASSWORD=euler psql -h 192.168.178.63 -U euler -d euler -f pub_schema.sql
-- ============================================================

CREATE SCHEMA IF NOT EXISTS pub;
SET search_path TO pub, public;

-- ============================================================
-- master_hits: lean main table
-- ============================================================
DROP TABLE IF EXISTS pub.master_hits CASCADE;
CREATE TABLE pub.master_hits (
    id           bigserial PRIMARY KEY,
    -- Euclid parameters (a > b > 0, m > n > 0, gcd=1, a-b and m-n odd)
    a            numeric NOT NULL,
    b            numeric NOT NULL,
    m            numeric NOT NULL,
    n            numeric NOT NULL,
    -- Body cuboid edges
    x            numeric NOT NULL,    -- U1*U2
    y            numeric NOT NULL,    -- V1*U2
    z            numeric NOT NULL,    -- U1*V2
    -- Primitive edges
    g_scale      numeric NOT NULL,    -- gcd(x, y, z)
    x_prim       numeric NOT NULL,
    y_prim       numeric NOT NULL,
    z_prim       numeric NOT NULL,
    -- Master-Hit invariants
    q            numeric NOT NULL,    -- third face diagonal
    f1           numeric NOT NULL,    -- (W1*U2)^2 + (U1*V2)^2
    -- Blocker count: NULL = unfactored, >=0 = number of odd-exponent primes
    num_blockers integer,
    UNIQUE (a, b, m, n)
);

CREATE INDEX idx_pub_mh_mn         ON pub.master_hits (m, n);
CREATE INDEX idx_pub_mh_blockers   ON pub.master_hits (num_blockers);
CREATE INDEX idx_pub_mh_xyz        ON pub.master_hits (x_prim, y_prim, z_prim);

COMMENT ON TABLE  pub.master_hits IS 'One row per body cuboid Master-Hit.';
COMMENT ON COLUMN pub.master_hits.f1 IS '(W1*U2)^2 + (U1*V2)^2; the space-diagonal obstruction.';
COMMENT ON COLUMN pub.master_hits.num_blockers IS 'Count of primes with odd exponent in f1; NULL until factorisation completes.';

-- ============================================================
-- fibers: one row per (m,n) for MW-active fibers
-- ============================================================
DROP TABLE IF EXISTS pub.fibers CASCADE;
CREATE TABLE pub.fibers (
    id              serial PRIMARY KEY,
    m               numeric NOT NULL,
    n               numeric NOT NULL,
    -- Weierstrass coefficients [a1, a2, a3, a4, a6]
    weierstrass     numeric[],
    conductor       numeric,
    rank_known      integer,           -- known rank (possibly only a lower bound)
    rank_proven     boolean DEFAULT false,
    torsion         text,              -- e.g. "Z/2Z x Z/4Z"
    generators      jsonb,             -- list of [X, Y] pairs
    UNIQUE (m, n)
);

CREATE INDEX idx_pub_fibers_rank ON pub.fibers (rank_known);

COMMENT ON TABLE  pub.fibers IS 'Elliptic fibre E_{m,n} of the Master-Hit variety, one row per (m,n) covered by Mordell-Weil generation.';

-- ============================================================
-- generator_groups: hierarchy of all families and provenance tags
-- ============================================================
DROP TABLE IF EXISTS pub.generator_groups CASCADE;
CREATE TABLE pub.generator_groups (
    id              serial PRIMARY KEY,
    short_name      text UNIQUE NOT NULL,    -- e.g. 'Saunderson', 'MW-88-7'
    parent_id       integer REFERENCES pub.generator_groups(id),
    category        text NOT NULL CHECK (category IN ('family', 'provenance')),
    description     text,
    citation        text,
    -- Optional: for MW subgroups, reference to fibers
    fiber_id        integer REFERENCES pub.fibers(id)
);

CREATE INDEX idx_pub_gg_parent ON pub.generator_groups (parent_id);
CREATE INDEX idx_pub_gg_category ON pub.generator_groups (category);

COMMENT ON TABLE pub.generator_groups IS
'Hierarchical taxonomy. category=family: algebraic family (what the brick IS). category=provenance: discovery origin (HOW the brick was found).';

-- ============================================================
-- hit_groups: n:m link Hit <-> Groups
-- ============================================================
DROP TABLE IF EXISTS pub.hit_groups CASCADE;
CREATE TABLE pub.hit_groups (
    hit_id          bigint NOT NULL REFERENCES pub.master_hits(id) ON DELETE CASCADE,
    group_id        integer NOT NULL REFERENCES pub.generator_groups(id),
    is_primary      boolean DEFAULT false,   -- True for the original source
    details         jsonb,                   -- e.g. {"scalar_vector": [1,0,-1], "torsion_index": 0}
    PRIMARY KEY (hit_id, group_id)
);

CREATE INDEX idx_pub_hg_group     ON pub.hit_groups (group_id);
CREATE INDEX idx_pub_hg_primary   ON pub.hit_groups (hit_id) WHERE is_primary;

-- ============================================================
-- f1_factors: prime factorization of f1
-- ============================================================
DROP TABLE IF EXISTS pub.f1_factors CASCADE;
CREATE TABLE pub.f1_factors (
    hit_id          bigint NOT NULL REFERENCES pub.master_hits(id) ON DELETE CASCADE,
    prime           numeric NOT NULL,
    exponent        integer NOT NULL,
    is_blocker      boolean NOT NULL,        -- exponent odd
    prime_mod4      integer,
    PRIMARY KEY (hit_id, prime)
);

CREATE INDEX idx_pub_f1f_blocker ON pub.f1_factors (hit_id) WHERE is_blocker;
CREATE INDEX idx_pub_f1f_prime   ON pub.f1_factors (prime);

-- ============================================================
-- db_metadata: version, snapshot date, statistics
-- ============================================================
DROP TABLE IF EXISTS pub.db_metadata CASCADE;
CREATE TABLE pub.db_metadata (
    key             text PRIMARY KEY,
    value           text NOT NULL,
    updated_at      timestamptz DEFAULT now()
);

INSERT INTO pub.db_metadata (key, value) VALUES
    ('schema_version', '1.0'),
    ('paper_reference', 'Peschmann, R.: An odd-exponent blocker theorem and a Mordell-Weil construction of body cuboids, 2026'),
    ('source_database', 'euler.public.master_hits');

-- ============================================================
-- Top-level groups (static definitions)
-- ============================================================

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation) VALUES
    -- Family Top-Level
    ('Classical',         NULL, 'family',     'Classical infinite parametric families', NULL),
    ('Modern',            NULL, 'family',     'Modern parametric families (post-2000)', NULL),
    ('Sporadic',          NULL, 'family',     'Master-Hits not covered by any known family', NULL);

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Saunderson', id, 'family', 'Saunderson 1740: a=4gh, b=g^2+h^2, m=h(3g^2-h^2), n=g(g^2-3h^2)',
       'Guy, Unsolved Problems in Number Theory, 3rd ed.'
FROM pub.generator_groups WHERE short_name='Classical';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Euler', id, 'family', 'Euler family from x^2+y^2=2z^2', NULL
FROM pub.generator_groups WHERE short_name='Classical';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Lenhart', id, 'family', 'Lenhart family from u^2+v^2=5w^2', NULL
FROM pub.generator_groups WHERE short_name='Classical';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Himane', id, 'family', 'Himane 2024: three theorems from pairs of Pythagorean triples', 'arXiv:2405.13061'
FROM pub.generator_groups WHERE short_name='Modern';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Himane-T1', id, 'family', 'Himane Theorem 1', 'arXiv:2405.13061'
FROM pub.generator_groups WHERE short_name='Himane';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Himane-T2', id, 'family', 'Himane Theorem 2', 'arXiv:2405.13061'
FROM pub.generator_groups WHERE short_name='Himane';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Himane-T3', id, 'family', 'Himane Theorem 3', 'arXiv:2405.13061'
FROM pub.generator_groups WHERE short_name='Himane';

-- Provenance Top-Level
INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation) VALUES
    ('Exhaustive-Search', NULL, 'provenance', 'Classical brute-force search up to a parameter bound', NULL),
    ('Rathbun-Search',    NULL, 'provenance', 'Catalogue from Rathbun (2017)',                       'arXiv:1705.03515'),
    ('Saunderson-Generator', NULL, 'provenance', 'Direct closed-form Saunderson parametrisation', NULL),
    ('Mordell-Weil',      NULL, 'provenance', 'Mordell-Weil generator on elliptic fibre E_{m,n}',    NULL);

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Exhaustive-Bound-200', id, 'provenance', 'Exhaustive search up to bound 200', NULL
FROM pub.generator_groups WHERE short_name='Exhaustive-Search';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Exhaustive-Bound-500', id, 'provenance', 'Exhaustive search up to bound 500', NULL
FROM pub.generator_groups WHERE short_name='Exhaustive-Search';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Exhaustive-Bound-2000', id, 'provenance', 'Exhaustive search up to bound 2000', NULL
FROM pub.generator_groups WHERE short_name='Exhaustive-Search';

INSERT INTO pub.generator_groups (short_name, parent_id, category, description, citation)
SELECT 'Exhaustive-Bound-2300', id, 'provenance', 'Exhaustive search up to bound 2300', NULL
FROM pub.generator_groups WHERE short_name='Exhaustive-Search';

-- ============================================================
-- View: hit_taxonomy for easy analysis
-- ============================================================
CREATE OR REPLACE VIEW pub.hit_taxonomy AS
SELECT
    h.id,
    h.a, h.b, h.m, h.n,
    h.x_prim, h.y_prim, h.z_prim,
    h.num_blockers,
    -- Aggregated family tags
    (SELECT array_agg(g.short_name ORDER BY g.short_name)
     FROM pub.hit_groups hg
     JOIN pub.generator_groups g ON g.id = hg.group_id
     WHERE hg.hit_id = h.id AND g.category = 'family') AS families,
    -- Primary provenance tag
    (SELECT g.short_name
     FROM pub.hit_groups hg
     JOIN pub.generator_groups g ON g.id = hg.group_id
     WHERE hg.hit_id = h.id AND g.category = 'provenance' AND hg.is_primary
     LIMIT 1) AS provenance
FROM pub.master_hits h;

COMMENT ON VIEW pub.hit_taxonomy IS 'Convenience view: each Master-Hit with families[] and primary provenance.';
