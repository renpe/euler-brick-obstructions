# Deadends — dokumentierte Sackgassen

Diese Skripte gehören zu Ansätzen, die wir verfolgt und als nicht zielführend identifiziert haben. Sie sind hier dokumentiert, weil ihre Ergebnisse für das Paper als „warum dieser naheliegende Ansatz nicht funktioniert" interessant sind.

## (g₊, g₋)-Faserung als MW-Hebel

- `gpm_subfiber_explore.py` — Verteilung der Hits über (m,n)-Subfasern in einer (g₊, g₋)-Klasse.
- `gpm_subfiber_mw_test.sage` — Test, ob (g₊, g₋)-Schicht in E_{m,n}(ℚ) eine Untergruppe/Coset ist.

**Ergebnis:** Nein. Die (g₊, g₋)-Etikettierung ist eine Repräsentanten-Konvention, kein Punkt-Invariant. Hinter ihr steckt nur eine 2-elementige Klasse mod Quadrate, die durch nicht-triviale Torsion (|tor|=8 in Top-Faser) verschmiert wird.

Memory: `gpm_fibration` (Negativ-Erkenntnis Abschnitt).

## E_{x,y}-Konstruktion als „dritte" elliptische Kurve

- `euler_xy_rank.sage` — Rang-Bestimmung der x,y-Kurve.
- `euler_xy_cuboid_intersect.sage` — Schnitt-Test mit Cuboid-Bedingung.
- `euler_xy_scale.py` — 100K-Stichprobe Brute-Force-Suche.

**Ergebnis:** Die Kurve E_{x,y} ist algebraisch identisch zu E_{m,n} nach einer Permutation der Parameter (siehe `curves_isogeny_check.sage` im Hauptverzeichnis: j-Invarianten stimmen überein). Damit liefert E_{x,y} keinen unabhängigen Mordell-Weil-Hebel.

Memory: `euler_xy_kurve` und `drei_kurven_kollaps`.

## Coleman-Stoll-Bound zu lose für unsere Genus-3-Kurven

- `jh_coleman_bound.sage` — Versuch, |H(F_p)| + 2g−2 oder Stoll's |H(F_p)| + 2r als rigorose Schranke zu nutzen.

**Ergebnis:** Auch mit Primzahlen bis 200 ist |H(F_p)| ≥ 8 immer (denn die 8 trivialen Q-Punkte reduzieren generisch zu 8 verschiedenen F_p-Punkten). Damit ist die Bound ≥ 10 (rk=1) bzw. ≥ 12 (rk=2) — zu lose um |H(Q)| = 8 zu erzwingen.

**Was stattdessen funktionierte:** Torsions-Trick (`jh_torsion_scale.sage` im Hauptverzeichnis) — direkt über |E_q(Q)| = 4 statt über Stoll.

## Verfeinerter Torsions-Trick mit Bug

- `jh_torsion_refined.sage` — Versuch, |tors|≠4-Fälle abzudecken durch Brute-Force-Enumeration auf der Quartik.

**Ergebnis:** Bug — Doppelzählung über (w, ±V)-Symmetrie + nicht erfasste Punkte im Unendlichen. Korrekt wäre: Sage's `torsion_subgroup()` direkt nutzen.

## Sage-API-Test für Genus-3-Jacobian-Rang

- `jacobian_rank_attempt.sage` — Versuch, `J = HyperellipticCurve.jacobian()` für Rang-Bestimmung zu nutzen.

**Ergebnis:** Sage 10.7 hat keine `rank()` oder `analytic_rank()` für Genus-3-Jacobians implementiert. Wir kommen über die explizite J(H) ~ E_PQ × E_uV × E_3 Zerlegung weiter, nicht über die Jacobian-Klasse direkt.
