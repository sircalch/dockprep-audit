# Pilot selection -- FROZEN 2026-08-21, EXPANDED 2026-08-22, RE-FROZEN 2026-08-23 (audit stage only)

**This cohort is frozen.** No case may be added, removed, or reassigned to a
different stratum without a new, dated entry in the Substitutions section
below explaining why. The authoritative frozen record is
`benchmark/pilot_manifest_frozen.csv`, checksummed at
`benchmark/pilot_manifest_frozen.sha256.txt` (SHA-256:
`ab37330651f5787bc25b4a7a301f7b7c4a98bfdd580139aa027898d6e5d4abfe` --
updated 2026-08-23, second re-freeze same day: `1SN5` got a receptor-chain
override, A;B;C;D instead of the geometric A;C, matching its
`REMARK 350`-declared tetrameric biological unit -- same pattern as the
`3FNU` override, see "Substitutions" below). The first 2026-08-23 re-freeze
(checksum `9313e6a7...6dbba`) substituted 9 cases after the technical smoke
test found genuine preparation-pipeline limitations, see "Substitutions"
below; the target-independence backfill noted below as a gap was also
completed in that re-freeze, now that `data.rcsb.org` is reachable again.
Cohort expanded from 12 to **60 cases (15 per stratum)** on 2026-08-22, Fase
E of `PROJECT-ROADMAP.md` section 9. See "Fase E expansion" below for the
original 48 new cases and how they were found/verified.

**Target-independence backfill completed 2026-08-23:** the `target_uniprot`
column is now populated for all 60 rows (was empty for the 48 new cases due
to a sustained RCSB outage during the 2026-08-22 freeze). No UniProt target
appears more than twice in the cohort except one incidental false positive
(`Q15788`/NCOA1, a shared coactivator peptide chain, not a drug target --
see "Substitutions" below for detail).

Regenerate both files with `python scripts/freeze_pilot_manifest.py` only
after a deliberate, documented cohort change -- never to silently refresh
provenance.

This is still a pilot for **audit-rule validation**, not a completed docking
benchmark, and must not be described as such in a manuscript. All files,
checksums and current automatic classifications are in
`pilot-inventory/cohort_candidates.csv`.

| PDB ID | Planned stratum | Why included | Status |
|---|---|---|---|
| 1M17 | Alternate location | EGFR complex; alternate-location rule triggered | Frozen |
| 1T46 | Alternate location | c-KIT-imatinib complex; site-local alternate-location rule triggered | Frozen |
| 4RJ3 | Alternate location | CDK2-inhibitor complex; site-local alternate-location rule triggered | Frozen |
| 5A2S | Metal/cofactor | HDAC inhibitor complex; zinc and sodium present | Frozen |
| 1CBX | Metal/cofactor | Carboxypeptidase A-L-benzylsuccinate complex; site-local zinc present | Frozen |
| 4EXS | Metal/cofactor | NDM-1-captopril complex; site-local zinc present | Frozen |
| 1OHR | Water policy | HIV protease-nelfinavir complex; water finding, no high-severity rule | Frozen |
| 3FNU | Water policy | Histo-aspartic protease-inhibitor complex; water finding, no high-severity rule | Frozen |
| 4GID | Water policy | Beta-secretase-inhibitor complex; water finding, no high-severity rule | Frozen |
| 1A28 | Low-risk control | Progesterone receptor complex; no high-severity rule | Frozen |
| 1QKT | Low-risk control | Estrogen receptor complex; no high-severity rule | Frozen |
| 1RBP | Low-risk control | Serum retinol-binding protein complex; no high-severity rule | Frozen |

## Checks completed before freezing (2026-08-21)

1. ✅ Biological assembly and intended receptor chain confirmed for all 12
   cases (`PROJECT-ROADMAP.md` section 15) -- 11/12 straightforward, `3FNU`
   resolved by explicit user decision (A+B, matching the real biological
   dimer; see `benchmark/receptor_chain_overrides.csv`).
2. ✅ Reference ligand identity and pose confirmed for all 12 cases via
   wwPDB validation metrics (RSCC/RSR, Mogul geometry outliers, occupancy;
   `PROJECT-ROADMAP.md` section 15). RSCC/RSR unavailable for 3 pre-1998
   depositions (`1CBX`, `1OHR`, `1RBP`) -- documented limitation, not fixable.
3. ✅ Water/metal binding-site relevance verified spatially (not just
   whole-structure presence) by `scripts/verify_pilot_eligibility.py`.
4. ✅ Full identifier registry of earlier cohorts searched (132 exclusion
   records, section 3 of `PROJECT-ROADMAP.md`), not only IDs named in the
   submitted manuscripts.
5. ✅ One substitution recorded before any docking run (`1CPS` -> `1CBX`,
   see below); no other exclusions pending.

Preparation-policy decisions (altLoc handling, water/metal retention,
protonation, box, seeds) are explicitly NOT part of this freeze -- they are
Fase B and are decided per-policy, not baked into the cohort.

## Biological-target independence (pendiente #3, decided 2026-08-21)

The previous-study exclusion registry (`benchmark/previous-study-exclusions.csv`,
132 records) is keyed by PDB ID only, so it cannot catch a *different* PDB
entry of the *same* biological target already used in a prior manuscript.
Checked by comparing UniProt accessions (unambiguous regardless of which
crystal form/PDB ID was deposited): **6 of the 12 frozen cases share a
UniProt accession with an excluded record**, all six from the same prior
source, `dude_receptor_prep_audit` (the 102-target DUD-E receptor-prep audit
submitted to *J. Mol. Model.*) -- not from the alternate-conformations
manuscript.

| This case | Target (UniProt) | Shares target with (excluded) |
|---|---|---|
| `1A28` | P06401 (progesterone receptor) | `3KBA` |
| `1M17` | P00533 (EGFR) | `2RGP` |
| `1QKT` | P03372 (estrogen receptor) | `1SJ0` |
| `1T46` | P10721 (KIT) | `3G0E` |
| `4GID` | P56817 (BACE1) | `3L5D` |
| `4RJ3` | P24941 (CDK2) | `1H00` |

**Decision: keep all 12 cases as frozen; document, do not replace.**
Reasoning: (1) these are different PDB entries -- different crystal, often
different resolution or ligand -- not reused data; (2) the prior study's
question (audit-rule prevalence across 102 DUD-E targets) is categorically
different from this study's question (does preparation policy change
redocking reproducibility); (3) requiring fully novel biological targets
for a risk-stratified docking pilot is close to infeasible in practice --
visible altLoc/water/metal richness concentrates in exactly the
well-characterized, frequently-crystallized targets the field already
studies most. Replacing 6 of 12 already-frozen, already-verified cases at
this stage was judged a worse trade-off than transparent disclosure. This
must be stated as an explicit limitation in the manuscript, with the table
above (or an equivalent) reproduced there. Full detail (including all
UniProt accessions and the specific excluded PDB/source per case) is
recorded in `benchmark/pilot_manifest_frozen.csv` (`target_uniprot`,
`shares_target_with_prior_pdb`, `shares_target_prior_study` columns),
generated reproducibly by `scripts/freeze_pilot_manifest.py`.

**Same check extended to the 48 Fase E cases, 2026-08-24:** the target-independence
backfill completed 2026-08-23 (once `data.rcsb.org` was reachable again) made
this the same one-line query as above, run against the full 60-case
manifest. The raw `target_uniprot` overlap flagged **17 of the 48** new
cases, but 3 of those are the `Q15788` (NCOA1 coactivator peptide) false
positive already noted under "Substitutions" above -- verified directly by
fetching the excluded record's own UniProt set (`2GTK` = `{Q15788,
P37231}`, i.e. PPAR-gamma + its bound coactivator peptide): `3BEJ` (FXR,
Q96RI1), `3KMR` (RARalpha, P10276) and `5HJS` (LXRbeta, Q13133) share only
the incidental coactivator chain with `2GTK`, not its real drug target, so
they are not genuine overlaps. That leaves **14 of 48 new cases with a
genuine target overlap**, all from `dude_receptor_prep_audit`, none from the
alternate-conformations manuscript. Combined with the original 6/12,
**20 of 60 frozen cases (33%)** share a real biological target with that
prior study.

| This case | Target (UniProt) | Shares target with (excluded) |
|---|---|---|
| `1G5Y` | P19793 (RXR-alpha) | `1MV9` |
| `1GS4` | P10275 (androgen receptor) | `2AM9` |
| `1HRN` | P00797 (renin) | `3G6Z` |
| `1O86` | P12821 (ACE) | `3BKL` |
| `1UZE` | P12821 (ACE) | `3BKL` |
| `1Z9Y` | P00918 (carbonic anhydrase II) | `1BCD` |
| `2OZ7` | P10275 (androgen receptor) | `2AM9` |
| `3GZ9` | Q03181 (PPAR-delta) | `2ZNP` |
| `3IBI` | P00918 (carbonic anhydrase II) | `1BCD` |
| `3KM4` | P00797 (renin) | `3G6Z` |
| `3NOA` | P37231 (PPAR-gamma) | `2GTK` |
| `4DJW` | P56817 (BACE1) | `3L5D` |
| `4RSY` | P09960 (LTA4 hydrolase) | `3CHP` |
| `6D8X` | P37231 (PPAR-gamma) | `2GTK` |

(`3BEJ`, `3KMR`, `5HJS` deliberately excluded from this table -- coactivator-only
false positive, see above.)

Same reasoning as the original 6/12 applies without modification: different
PDB entries, a categorically different research question, and the same
practical constraint that structurally-rich targets concentrate in the
field's most-studied proteins. **Decision: keep all 48 cases as frozen;
document, do not replace** -- consistent with not having replaced the
original 6, and with the substitution decisions elsewhere in this file
(which were made for genuine tool/chemistry failures, never to chase target
novelty).

## Fase E expansion (2026-08-22): 12 -> 60 cases

48 new cases were added, 12 per stratum, using the same audit engine and the
same spatial-relevance rule as the original 12 (altLoc/metal/bridging-water
within 6 A of the declared ligand; bridging water additionally requires
<=3.0 A to both ligand and receptor). Found via RCSB full-text/keyword
search across many protein families per stratum, then filtered
programmatically -- not hand-picked to fit a narrative. Full search log and
reasoning is in `PROJECT-ROADMAP.md` section 9 (Fase E). All 48 pass
`scripts/verify_pilot_eligibility.py` with zero manual review needed, and
none appear in the 132-record previous-study exclusion registry.

**Cap of 2 cases per biological target (UniProt), decided during the
search:** an initial "kinase inhibitor" search returned 6 of 13 alternate-
location candidates as the same protein (human PKA) -- accepting that would
have made the stratum a de facto single-target test dressed up as 12
independent cases. The cap was applied retroactively across all four
strata; some cases share a target with another case at exactly 2 (e.g.
`1O86`/`1UZE` are both human ACE, `4G9L`/`4JA1` are both MMP-3), documented
per-case in `benchmark/expansion_cases.csv` and inferrable from
`pilot_manifest_frozen.csv`'s `title` column.

**One quality-control catch worth recording:** two candidates found via an
"HIV-2 protease inhibitor" search (`1WBK`, `1WBM`) turned out, per their own
PDB `TITLE`/`COMPND` records, to be **HIV-1** protease (RCSB's full-text
search had matched a passing mention, not the actual target) -- caught by
reading the file header before accepting either, not by trusting the search
query's framing. Only `1WBK` was kept, pairing with the original `1OHR`
(also HIV-1 protease) at the 2-per-target cap.

New cases by stratum: `alternate_location` -- `6FTF`, `3KXG`, `3MWU`,
`7ORS`, `3P0M`, `2I4H`, `5E0J`, `5TG1`, `1E6U`, `1E7S`, `3D14`, `5F4N`.
`metal_or_cofactor` -- `3LXE`, `6ZR9`, `1Z9Y`, `3IBI`, `1O86`, `1UZE`,
`1DTH`, `4G9L`, `4JA1`, `6TMN`, `1THL`, `4RSY`. `water_policy` -- `1EPQ`,
`1ENT`, `3KM4`, `4DJW`, `6ASH`, `1PSO`, `1QRP`, `1WKR`, `1XDH`, `1CVZ`,
`2F25`, `1WBK`. `low_risk_control` -- `1GS4`, `2OZ7`, `3NOA`, `6D8X`,
`3GZ9`, `1SN5`, `1G5Y`, `1IE9`, `3FLI`, `3L1B`, `5HJS`, `3KMR`.

**This list is as originally found on 2026-08-22; 9 of these entries
(`6FTF`, `7ORS`, `1THL`, `1PSO`, `1WKR`, `1XDH`, `1EPQ`, `1ENT`, `3FLI`) were
substituted on 2026-08-23 after failing the technical smoke test -- see
"Substitutions" below for the current 9 replacement PDB IDs and why.**

**Technical smoke test completed 2026-08-23 for all 48 new-cohort slots**
(`scripts/run_expansion_smoke.py` plus manual reruns for the 9 substituted
cases): default policy (no water, exhaustiveness 8, seed 42) runs end-to-end
without a crash for all 48 current slots. Outcomes are diagnostic only (this
is a pipeline smoke test, not the Fase C scientific comparison), and are not
yet run at the full multi-seed/multi-policy Fase C setting.

**Ligand identity/pose (wwPDB validation metrics) and biological-assembly
confirmation completed 2026-08-23** for all 48 -- see
`benchmark/expansion-validation/wwpdb_validation.csv` and the `1SN5`
override entry above. Still pending: the full Fase B policy pipeline
(altLoc/His-tautomer/protonation/box policy decisions, as opposed to the
technical-default smoke test) and Fase C multi-seed runs.

## Substitutions

- **9 Fase E cases substituted (2026-08-23):** the technical smoke test run
  on all 48 new cases found 9 with a genuine preparation-pipeline limitation
  (not a project bug), each verified individually before substitution, same
  standard as `1CPS` below:
  - `1PSO`, `1WKR`, `1XDH` (`water_policy`): all three are pepstatin-bound
    aspartic proteases whose inhibitor contains the non-standard statine
    residue (`STA`), which breaks Meeko's polymer bond-padding logic
    (`Expected 2 paddings ... but got 1`) -- a known gap in Meeko's residue
    template coverage for peptidomimetics, not fixable by this project.
  - `1EPQ`, `1ENT` (`water_policy`), `6FTF`, `7ORS` (`alternate_location`),
    `1THL` (`metal_or_cofactor`), `3FLI` (`low_risk_control`): all six have
    a ligand for which no protonation variant proposed by dimorphite-dl could
    be built as a valid RDKit molecule (nitrogen valence errors) -- a real
    chemistry-complexity limit, same category as `1CPS`'s companion `7CI`
    ligand failure.
  - A first replacement attempt for `1EPQ`, `1MEM` (cathepsin K + a vinyl-
    sulfone inhibitor), also failed technical smoke testing: its ligand is
    covalently bonded to the protein backbone, hitting the identical Meeko
    padding-template limitation as the `STA` cases. Discarded in favor of
    `1EED`.
  - Every replacement below was accepted only after passing BOTH formal
    eligibility (`scripts/verify_pilot_eligibility.py`) AND the full
    technical pipeline (`scripts/smoke_redock_case.py`: extraction,
    protonation, Vina redocking) without error -- a `scoring_fail` or
    `sampling_fail` docking outcome is fine, only a crash disqualifies.

| Replaced (reason) | Replacement | Ligand (chain/resseq) | Why chosen |
|---|---|---|---|
| `6FTF` (ligand chemistry) | `3PNA` | CMP / A / 250 | PKA RIα docking domain + cAMP; altLoc 5.46 Å from ligand; new target (P00514) |
| `7ORS` (ligand chemistry) | `5K8S` | CMP / A / 501 | *P. falciparum* PKA-R domain + cAMP; altLoc 4.24 Å from ligand; new target (Q7KQK0) |
| `1THL` (ligand chemistry) | `1KJO` | THR / A / 1317 | Thermolysin + Z-L-Thr product; Zn 2.05 Å from ligand; same family as `6TMN` (P00800, 2/2 at cap) |
| `1PSO` (STA/pepstatin) | `1EPP` | 1Z1 / E / 333 | Endothiapepsin + small-molecule inhibitor; 9 waters within 4 Å (P11838) |
| `1WKR` (STA/pepstatin) | `1PPM` | 0P1 / E / 324 | Fungal aspartic protease transition-state mimic, small molecule; 10 waters within 4 Å (P00798) |
| `1XDH` (STA/pepstatin) | `1HRN` | 03D / A / 391 | Human renin, high-res, small-molecule inhibitor; 7 waters within 4 Å; same target as `3KM4` (P00797, 2/2 at cap) |
| `1EPQ` (ligand chemistry; `1MEM` also tried, covalent ligand, failed) | `1EED` | 0EO / P / 327 | Endothiapepsin + cyclohexyl renin-inhibitor-derived compound; 10 waters within 4 Å; same target as `1EPP` (P11838, 2/2 at cap) |
| `1ENT` (ligand chemistry) | `1WBM` | BLL / B / 1100 | HIV-1 protease symmetric inhibitor BEA450; 10 waters within 4 Å; same target as `1WBK`/`1OHR` (Q8Q3H0/P03366 family, 2/2 at cap) |
| `3FLI` (ligand chemistry) | `3BEJ` | MUF / A / 473 | FXR LBD + agonist MFA-1 + coactivator peptide; no altLoc/metal within 6 Å of ligand; FXR itself is a new target (Q96RI1) |

**Note on `3BEJ`'s coactivator chain:** `3BEJ`'s `target_uniprot` includes
`Q15788` (NCOA1/SRC-1), a coactivator peptide co-crystallized alongside the
nuclear-receptor ligand-binding domain -- the same peptide chain also appears
in `3KMR` (RARα) and `5HJS` (LXRβ), pushing its raw occurrence count to 3.
This is not a real target-duplication: the coactivator peptide is not the
drug target in any of the three structures, and the three actual receptor
targets (FXR/Q96RI1, RARα/P10276, LXRβ/Q13133) are all distinct. The 2-per-
target cap is applied to the drug target, not to incidental shared
secondary chains.

Working files for this substitution round: `benchmark/replacement-round2/`
(candidate search, eligibility, and smoke-test logs).

- **`1SN5` receptor-chain override, A;B;C;D (2026-08-23):** found during the
  wwPDB-validation/chain-assembly pass over all 48 new cases (results in
  `benchmark/expansion-validation/wwpdb_validation.csv`). `1SN5`'s
  `REMARK 350` declares `AUTHOR DETERMINED BIOLOGICAL UNIT: TETRAMERIC`, but
  the 6 Å geometric rule alone found only chains A;C near the bound `T3`
  (triiodothyronine) ligand -- chains B and D sit at 6.61/6.86 Å, just past
  the cutoff, not a distant crystallographic copy (contrast: two other
  cases flagged by the same pass, `1E6U`/`1E7S`, had their "extra copy"
  measured at >12 Å and correctly left as single-chain). Transthyretin's
  T3/T4 binding site is documented in the literature to sit in the central
  channel formed between the tetramer's two dimers, so B/D plausibly
  contribute to the real binding pocket -- same reasoning as the `3FNU`
  override below. Applied in `benchmark/receptor_chain_overrides.csv`,
  re-verified (60/60 pass, raw geometric finding A;C preserved in
  `receptor_chains_geometric` for transparency) and re-run through the
  technical smoke test (`sampling_fail`, no preparation error).

- **1CPS -> 1CBX (2026-08-21):** `1CPS` was dropped after a technical smoke
  test (see `PROJECT-ROADMAP.md` section 15) found that its residue
  `TYR A:204` has a carbonyl oxygen simultaneously ~1.23 A from C, ~1.28 A
  from CA and ~1.69 A from N of the same residue -- geometrically
  incompatible with a real carbonyl, and a hard blocker for automated
  receptor preparation. This is a genuine deposited-structure issue, not a
  parsing bug, and was not patched with a repair heuristic. `1CBX`
  (carboxypeptidase A + L-benzylsuccinate, zinc at 2.305 A from the ligand)
  was already a vetted round-2 replacement candidate that had passed spatial
  eligibility; it was re-verified against the current 12-case cohort and
  additionally passed the full technical pipeline (extraction, PDBQT
  preparation, Vina redocking, RMSD 0.682 A).

## Current limitation

The classifications in this pilot identify structural features in the source
PDB file. They do not infer protonation, biological relevance, binding-site
occupancy or a preferred repair. Those are deliberately separate decisions.
