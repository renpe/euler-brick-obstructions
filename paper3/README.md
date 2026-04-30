# Computational scripts for Paper 3

*A torsion-intersection refinement: rigorous proof of the perfect-cuboid conjecture on an explicit family of master-tuple fibers*,
R. Peschmann (in preparation).

## Verhältnis zu Paper 1

Paper 3 ist ein **Companion Paper** zu Peschmann (2026), arXiv:2604.09328 (*„Quartic reductions and elliptic obstructions for perfect Euler bricks"*, Paper 1). Paper 1 etabliert die Genus-3-Reduktion `C_A: w² = λ⁸ + Aλ⁴ + 1`, die Klein-4-Involutionen mit drei elliptischen Quotienten `E_A × E_A' × E_A''`, sowie Obstruktionen via Kummer-Charakter und 2-Descent auf `E_A`. Paper 3 baut darauf auf und liefert zwei spezifische neue Beiträge:

1. **Reduktionssatz (Theorem 4):** jeder primitive Euler-Brick stammt von einem Master-Tupel ab. Verstärkt Paper 1's Lemma 3.1 (nur Vorwärts-Implikation).

2. **Torsions-Schnitt-Argument auf E_A''** (= E_3 in unserer Notation): wenn rk(E_A'')=0 mit |tors|=4, dann ist H(ℚ) = {8 triviale Punkte}, ergo kein perfekter Cuboid auf der (m,n)-Faser. Paper 1 nutzt nur Spezialisierungen von E_A und E_A'.

Konkrete Liste: **117 (m,n)-Fasern** bei M_MAX=50, auf denen Konjektur B rigoros bewiesen ist (disjunkt von Paper 1's 42+54 Spezialisierungen).

## Hauptergebnisse

### Theorem 1 (g₊·g₋-Struktursatz)
Für jedes Master-Tupel: g₊·g₋ = gcd(U₁, U₂) = g_scale, gcd(g₊, g₋) = 1.

### Theorem 2 (Polynomielle Faktorisierung)
f₁ = L·R mit L, R primitive Summen zweier Quadrate.

### Theorem 4 (Reduktion auf Master-Tupel)
Jeder primitive Euler-Brick (X, Y, Z) entsteht aus einem Master-Tupel via Standardkonstruktion + Skalierung mit g = gcd(U₁, U₂).

### Hauptsatz (Theorem 4.4)
Für (m,n) mit rk(E_3(ℚ))=0 und |E_3(ℚ)_tors|=4: |H_{m,n}(ℚ)|=8 (alle trivial). Kein perfekter Cuboid auf dieser Faser.

## Status

- ✅ Theoreme 1, 2, 4: rigoros bewiesen (Paper 3).
- ✅ Genus-3 + 3-Quotienten-Zerlegung: zitiert aus Paper 1 (§3, §4).
- ✅ Hauptsatz auf 117 expliziten Fasern bei M_MAX=50.
- ⚠️ ~61 Fasern mit |tors|≠4: Verfeinerung möglich.
- ❌ ~340 Fasern mit Total-Rang ≥ 3: Quadratic Chabauty erforderlich.
- ❌ Unendlich viele (m,n) mit max(m,n) > 50: Forschungsproblem.

Siehe `CONJECTURE_B_ROADMAP.md`.

## Anforderungen

- Wie Paper 2: PostgreSQL, Python 3.10+, Sage 10.7+, PARI/GP 2.15+
- `pub_db.py` wird via Symlink aus paper2/_common/ geteilt.

## Struktur (öffentlich, in `scripts/paper3/`)

```
scripts/paper3/
├── README.md                            (dieses File)
├── _common/
│   └── pub_db.py → ../../paper2/_common/pub_db.py  (symlink)
├── data/
│   └── proven_fibers.csv                (117 Fasern als CSV)
└── analysis/
    ├── 15 Hauptergebnis-Skripte
    └── deadends/
        ├── README.md                    (warum diese Ansätze nicht funktionieren)
        └── 8 Sackgassen-Skripte
```

Das Manuskript-Material (LaTeX-Quelle, Beweis-Markdowns, kompiliertes PDF) liegt komplett außerhalb dieses Repositories und wird unabhängig vom Code veröffentlicht. Dieser README beschreibt nur die öffentlich-deploybaren Reproduktions-Skripte.

## Reproduktion der Hauptergebnisse

### 117 Fasern rigoros (Hauptsatz)

```bash
sage scripts/paper3/analysis/jh_torsion_scale.sage 50
```

Output: 518 Fasern geprüft, 117 rigoros bewiesen, 61 unsicher (|tors|≠4).

### Genus-3-Form-Verifikation (Theorem 3 / Recap aus Paper 1)

```bash
sage scripts/paper3/analysis/cuboid_genus_check.sage
sage scripts/paper3/analysis/jh_three_quotients.sage 50
```

### Empirische Höhensuche (B = 10⁶)

```bash
sage scripts/paper3/analysis/jh_chabauty_finalize.sage 1000000
```

11 Fasern × 30 Sekunden = ~5 min, alle finden nur die 6 trivialen affinen Punkte.

## Cross-References

- **Paper 1** (arXiv:2604.09328): Genus-3-Setup, Jacobian-Zerlegung, Kummer-Charakter, 2-Descent.
- **Paper 2** (in Vorbereitung): Master-Hit-Datenbank, MW-Konstruktion neuer Bricks, Blocker-Analyse.
- Paper 3 nutzt `pub_db.py` und das `pub`-Schema aus Paper 2.

## Memory

Relevant Memory-Einträge (für Sessions):
- `gpm_struktursatz` — Theorem 1
- `f1_polynomfaktor` — Theorem 2
- `cuboid_genus2_zerlegung` — Theorem 3 (zitiert aus Paper 1)
- `drei_kurven_kollaps` — Hilfslemma
- `parametrisierung_vollstaendigkeit` — Vorform von Theorem 4
- `paper3_setup` — Paper-5-Positionierung relativ zu Paper 1 und 4
