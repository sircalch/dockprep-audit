# Manuscript skeleton — Article 1 (target: Journal of Molecular Modeling)

Status: outline only, no prose written yet. Every bracketed [DATA: ...] note
points to the exact frozen file/figure to pull the number or claim from --
nothing here should be written from memory once drafting starts.

**Template source:** `manuscript/reference/JMM-Manuscript-revised.docx` (the
author's own prior accepted-track submission to this same journal, the
DUD-E 102-receptor prep-audit study `dude_receptor_prep_audit` referenced
throughout PROJECT-ROADMAP.md). Converted to
`manuscript/reference/jmm_manuscript_reference.md` for reuse. Match its
structure, section labels, citation style, and epistemic tone exactly
unless there's a specific reason to deviate -- it is real evidence of what
this journal/editor already accepted from this author, not a guess.

## Working title

"DockPrep Audit: auditing receptor structural risk factors before docking,
and when they predict preparation-policy sensitivity" (placeholder --
revisit once Results are drafted; should name the actual finding, not just
the tool). Style precedent from the reference paper's title: states the
method class + the two things it does ("provenance-aware workflow for
strict receptor-preparation audits AND reference-pose recovery") -- mirror
that "does X and predicts Y" shape rather than a generic tool-name title.

## Authors / affiliations

**RESOLVED from the reference manuscript** -- same author team, likely same
order unless the user says otherwise:

1. Andrés Monreal Hernández (corresponding), Universidad Estatal de Sonora,
   Hermosillo, Sonora, Mexico. ORCID: 0009-0009-1207-8597.
   andres.monreal@ues.mx
2. Sara Lizbeth Franco Amaya, Doctorado en Nanotecnología, Universidad de
   Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0005-0272-0241
3. Carlos Ivanhoe Martínez Osorio, Doctorado en Ciencia de Materiales,
   Universidad de Sonora, Hermosillo, Sonora, Mexico. ORCID:
   0009-0003-7872-4965

[TODO: confirm this is still the intended team/order for THIS paper, and
whether contribution roles shift given this project's actual division of
work.]

## Abstract (~250 words, write last)

Must state: (1) the problem -- receptor preparation choices are often made
silently; (2) the tool -- DockPrep Audit detects altLoc/water/metal risk
factors pre-preparation; (3) the test -- 60-case stratified redocking
benchmark, 2 preparation policies x 3 seeds; (4) the finding -- the
water_policy stratum shows a distributed water-retention benefit (5/15
cases help clearly vs 1-2/15 in other strata) [DATA: PROJECT-ROADMAP.md
section 9, "Fase D sobre los 60 casos"]; (5) the implication -- audit-flagged
risk categories are weakly but measurably predictive of when preparation
policy matters.

## Keywords

molecular docking; receptor preparation; reproducibility; AutoDock Vina;
structural audit; alternate conformations; crystallographic water; benchmark

---

## 1. Introduction

- Motivation: docking reproducibility depends on receptor preparation
  choices (altLoc handling, water retention, metal treatment, protonation)
  that are rarely documented explicitly in published docking studies.
- Prior work this differentiates from [DATA: PROJECT-ROADMAP.md section 3]:
  (a) the author's own DUD-E 102-receptor prep-audit study
  (`dude_receptor_prep_audit`, submitted to this same journal);
  (b) the alternate-conformations/pose-recovery manuscript on the prior
  cohort. State explicitly how this study's question differs (see section 3
  of the roadmap for the exact differentiation language already agreed).
- The research question, verbatim (already fixed since project start):
  "Do structural features detectable before receptor preparation predict
  when different preparation policies change docking reproducibility?"
- Contribution statement: (1) an open-source rule-based audit engine
  (DockPrep Audit v0.1.0) that flags altLoc/water/metal risk factors from a
  raw PDB file before any preparation step; (2) a 60-case, 4-stratum
  (15/stratum) redocking benchmark with frozen provenance; (3) evidence that
  the audit-flagged risk category is weakly predictive of preparation-policy
  sensitivity, strongest in the water-policy stratum.

## 2. Materials and Methods

### 2.1 The DockPrep Audit engine

- Rule-based detection: ALTLOC_PRESENT, WATERS_PRESENT, METAL_PRESENT
  [DATA: benchmark/pilot-inventory/reports/*.json schema; describe severity
  levels high/medium/low].
- Explicitly NOT a preparation tool -- audits, does not fix or choose.

### 2.2 Cohort construction and freezing

- Original 12-case pilot (3 per stratum), expanded to 60 (15 per stratum) in
  Fase E [DATA: benchmark/PILOT_SELECTION.md, full search methodology and
  2-per-UniProt-target cap rationale].
- Eligibility rule: spatial relevance within 6 A of the declared ligand
  instance (altLoc/metal), bridging-water criterion <=3.0 A to both ligand
  and receptor [DATA: scripts/verify_pilot_eligibility.py].
- Freezing/provenance: SHA-256 checksums of source PDB files, frozen
  manifest regenerated only on documented, dated changes
  [DATA: benchmark/pilot_manifest_frozen.csv, current checksum
  `ab37330651f5787bc25b4a7a301f7b7c4a98bfdd580139aa027898d6e5d4abfe`].
- 10 total substitutions across the project's life (1CPS->1CBX plus 9 Fase E
  cases) and 2 documented receptor-chain overrides (3FNU, 1SN5) where the
  6 A geometric rule disagreed with the REMARK 350 biological assembly
  [DATA: benchmark/PILOT_SELECTION.md "Substitutions" section -- state the
  reason for each substitution category: STA/pepstatin Meeko-template gap,
  and ligand-chemistry protonation failures].
- Biological-target independence disclosure: 20/60 cases (33%) share a
  UniProt target with the author's own prior DUD-E study; state explicitly
  as a limitation, not hidden [DATA: benchmark/PILOT_SELECTION.md
  "Biological-target independence" section, full table].

### 2.3 Preparation policies (frozen before results were observed)

- altLoc: highest occupancy retained, tie-break to conformer 'A'
  [DATA: PROJECT-ROADMAP.md section 15, altLoc policy discussion].
- Water: two compared policies -- `conservative_water` (bridging water
  retained, <=3.0 A to both ligand and receptor) vs `simplified_no_water`
  (all water removed). Note the originally-planned 3rd policy was merged
  with policy 2 after finding them computationally identical in this
  pipeline (no separate open-flow difference existed).
- Metal ions: always retained, native Meeko templates (Zn/Mg/Ca/Mn/Fe); His
  tautomer assigned geometrically when metal-coordinating (8/12 originally
  needed HID, not Meeko's default HIE).
- Ligand protonation: dimorphite-dl at physiological pH 7.4 (replaces
  `scrubber`, incompatible with Python 3.12).
- Docking box: cubic, 2.9 x ligand radius of gyration, centered on ligand
  centroid (Feinstein & Brylinski 2015) -- cite this reference properly in
  Methods.
- Vina settings: exhaustiveness 32, 3 predeclared seeds (42, 123, 2024) per
  case x policy, `--cpu 0` (multi-threaded, non-deterministic by Vina's own
  documented behavior -- state this explicitly as a reason for 3 seeds).

### 2.4 Redocking validation and metrics

- Heavy-atom RMSD to the crystallographic ligand pose.
- Exact-coordinate atom matching between reference ligand and Meeko's
  PDBQT output (tolerance 0.01 A) -- state WHY: Meeko's torsion-tree writer
  reorders atoms unpredictably, so the `REMARK SMILES IDX` hint cannot be
  trusted for correspondence.
- 3-way outcome classification: success (top-1 pose <=2.0 A), scoring_fail
  (some pose <=2.0 A but not ranked first), sampling_fail (no pose
  <=2.0 A) -- standard in redocking literature, cite an Astex-Diverse-style
  precedent.
- Both top-1-pose and best-of-9-poses metrics reported; explain why: the
  12-case pilot found the water-policy effect on `1M17` was invisible under
  the binary top-1 threshold and only appeared when using best-of-N
  [DATA: PROJECT-ROADMAP.md section 9, Fase D 12-case investigation].

## 3. Results

### 3.1 Cohort characteristics

- Finding frequency and co-occurrence across the 60 cases
  [FIGURE: fig_finding_frequency.png].
- Note explicitly that "WATERS_PRESENT" here is the generic any-water rule
  (59/60 cases), distinct from the site-local bridging-water criterion used
  to stratify `water_policy` -- do not conflate the two in the text.

### 3.2 Pipeline robustness

- 360/360 Fase C runs completed without a preparation error across the
  frozen 60-case cohort [DATA: benchmark/phase-c-60/phase_c_summary.csv].
- Report the smoke-test failure/substitution process as a methods-honesty
  note, not hidden: 12/48 new cases failed a first pass (3 retriable, 9
  genuine tool/chemistry limits requiring substitution) -- this is evidence
  the audit + substitution protocol is auditable, not evidence of a flawed
  cohort.

### 3.3 Per-case redocking outcomes

- [FIGURE: fig_heatmap_case_policy.png] -- full 60 x 2 heatmap.
- [FIGURE: fig_rmsd_distribution.png] -- violin distributions per
  stratum x policy.

### 3.4 Stratum-level water-policy effect (the central finding)

- n=15/stratum aggregate table [DATA: benchmark/phase-c-60/phase_d_by_stratum.csv]:
  water_policy succ(best) 0.644 (conservative) vs 0.556 (simplified), the
  only stratum with a meaningful gap; other 3 strata within +-0.044.
- Per-case breakdown: 5/15 water_policy cases show a clear benefit (>0.3 A
  improvement: 1CVZ, 1EPP, 1WBK, 4GID, 6ASH), 2/15 show the opposite
  (1EED, 2F25), 8/15 show no effect -- state plainly this is NOT
  deterministic.
- Contrast with the other 3 strata (1-2/15 "helps" cases each) to support
  the claim that the effect concentrates in the audit-flagged category,
  not uniformly across the cohort.
- Explicitly reference the already-refuted simple hypothesis from the
  12-case pilot (H-bond distance <=3.2 A to a ligand heteroatom does NOT
  predict which cases respond) so the paper does not overclaim a clean
  predictive rule -- report the refutation, not just the positive finding.

### 3.5 Example geometries

- [FIGURE: fig_binding_site_examples.png] -- one worked example per finding
  type (1M17 altLoc, 1OHR bridging water -- the classic HIV-protease flap
  water, 1CBX Zn coordination), real coordinates, not schematic.

## 4. Discussion

- Interpret the water-policy finding: audit-detectable bridging water is a
  weak but real predictor of preparation-policy sensitivity, concentrated
  in ~1/3 of the flagged stratum -- useful as a triage signal, not a
  deterministic rule.
- Explicitly discuss why the effect is not universal within the flagged
  stratum -- open question, not resolved in this study (candidate for
  follow-up: cavity geometry after water removal, docking-score-function
  desolvation modeling limitations already flagged as Limitation #1 in
  PROJECT-ROADMAP.md section 17).
- Position relative to the tool landscape: DockPrep Audit is not a docking
  engine and does not compete with Schrodinger/AMDOCK-style tools; it
  targets the preparation-decision transparency gap upstream of docking
  (see the open-access positioning discussion already had with the user --
  MolProbity/PDB2PQR are the closer conceptual neighbors, not docking GUIs).

## 5. Limitations

Pull directly, near-verbatim, from PROJECT-ROADMAP.md section 17 (already
written and reviewed):
1. Water model has no explicit desolvation term (rigid TIP3P charge only).
2. Metals lack explicit tetrahedral coordination geometry (no AD4Zn).
3. 20/60 cases (33%) share a biological target with the author's prior
   DUD-E study -- documented, not silently avoided.
4. n=15/stratum is larger than the original pilot but still not a formal
   statistical test; report as descriptive/exploratory evidence.
5. Structural-validation coverage gap correlates with deposition era, not a
   clean pre-1998 cutoff (11/48 new cases lack RSCC/RSR, verified directly
   against missing structure-factor files at RCSB).
6. Vina is not fully deterministic even with fixed seeds under multi-thread
   (`--cpu 0`) -- mitigated with 3 seeds/case, not eliminated.

## 6. Conclusion

One paragraph: restate the research question, the answer (weak, real,
stratum-concentrated predictive signal for the water-policy category), and
the practical recommendation (audit before preparing; treat a water-policy
flag as a signal to compare both policies, not to assume one is right).

## Data and code availability

- GitHub repo: [TODO: URL once the separate public repo exists]
- Zenodo DOI: [TODO: after archiving the frozen version]
- Frozen manifest checksum: `ab37330651f5787bc25b4a7a301f7b7c4a98bfdd580139aa027898d6e5d4abfe`
- No PDB coordinates redistributed; only identifiers/URLs/checksums per
  RCSB redistribution practice.

## Declarations

- Funding: [TODO]
- Conflicts of interest: [TODO]
- Author contributions: [TODO]

---

## Open items before drafting prose

1. **RESOLVED (format/style):** matched to `manuscript/reference/jmm_manuscript_reference.md`:
   - Abstract has exactly two labeled parts, `Context:` and `Methods:` --
     no separate `Results:`/`Conclusions:` label. Results/Discussion content
     gets folded into `Context:`'s closing sentences instead (see reference
     abstract: it states the headline numbers inline under `Context:`).
     Adapt: our abstract should state the water-policy finding inline the
     same way, not as a separate labeled clause.
   - 5 keywords, semicolon-separated, lowercase except proper nouns.
   - Numbered sections `1 Introduction`, `2 Materials and methods` (with
     `2.1`, `2.2`, ... subsections), `3 Results`, `4 Discussion`,
     `5 Conclusions` -- no separate numbered Limitations section in the
     reference (limitations are folded into Discussion as explicit
     sentences, e.g. "insufficient to benchmark docking broadly", "ADP
     template automatically constructed... retained as an explicit
     limitation"). **Decision needed:** keep our Limitations as its own
     subsection (more scannable, matches PROJECT-ROADMAP.md section 17
     already being separately maintained) or fold into Discussion to match
     the reference exactly -- lean toward keeping it separate since we have
     6 distinct limitations, more than the reference's inline handful.
   - `# Figure captions` as its own section after Conclusions, each caption
     as "Fig. N <description, 1-2 sentences, no trailing period>" followed
     by the image.
   - `# References`, numbered in order of first citation (Vancouver-style),
     full author lists, journal abbreviations, DOI on every entry.
   - `# Statements and declarations` as the final section, exact
     subsection headers: `## Funding`, `## Competing interests`,
     `## Author contributions`, `## Data availability`. Use these four,
     same order, same header wording.
   - Tone: precise, bounded claims; the reference repeatedly draws an
     explicit "evidence boundary" (its own 3.3 subsection) stating what the
     data does NOT show. Our paper should do the same for the water-policy
     finding: state plainly it is a stratum-level descriptive pattern, not
     a validated predictive rule, and is not evidence about any specific
     drug target's biology.
2. **RESOLVED (main-text/supplementary split decided 2026-08-27):** no
   explicit word/figure-count limit found in the reference file itself, but
   it calibrates length (~2500 words, 2 figures, 2 tables, 15 references).
   Our 7 figures split as:
   - **Main text (4):** `fig_finding_frequency.png` (cohort characterization,
     Results 3.1), `fig_heatmap_case_policy.png` and
     `fig_rmsd_distribution.png` (the core water-policy finding, Results
     3.3-3.4), `fig_binding_site_examples.png` (mechanistic illustration,
     Results 3.5 -- strong for reviewers, real coordinates not schematic).
   - **Online Resource 1 (3):** `fig_graphical_summary.png`,
     `fig_architecture.png`, `fig_decision_tree.png` -- these explain how
     the tool works, not what was found; useful for reproducibility but not
     load-bearing for the scientific narrative, matching the reference
     paper's own pattern of putting protocol/verification detail in
     supplementary material rather than the main figures.
   - Two tables planned for the main text, mirroring the reference's
     Table 1/Table 2 pattern: (a) stratum-level water-policy summary
     (from `benchmark/phase-c-60/phase_d_by_stratum.csv`), (b) per-case
     breakdown of the 5 helps/2 hurts/8 no-effect water_policy cases.
3. Confirm co-author list/order for THIS paper -- reference gives the
   default (Monreal Hernández corresponding + Franco Amaya + Martínez
   Osorio), but contribution roles will differ since the actual work split
   on DockPrep Audit hasn't been the same as the DUD-E audit paper.
4. **RESOLVED:** numbered Vancouver-style citations, in order of first
   appearance, as shown in the reference's References section.
5. **RESOLVED, and PROJECT-ROADMAP.md needs a small consistency fix:** Fase
   C already ran both policies x 3 seeds x exhaustiveness 32 on all 60
   frozen cases (360/360 runs, section 9 "Fase D sobre los 60 casos"). The
   "Trabajo restante" sentence near the end of that same subsection still
   says "Fase C multi-semilla sobre los 60 casos" as pending -- that's
   stale, written before the run completed later the same session. Fix
   that line in PROJECT-ROADMAP.md before drafting Methods, so the
   manuscript's Methods section doesn't accidentally understate what was
   actually done.
6. **DECIDED (2026-08-27): no AI-assistance disclosure paragraph.** The
   reference manuscript includes one (naming OpenAI Codex), but the user
   confirmed every result in this project was manually validated and does
   not want an equivalent statement here. Do not add a Methods 2.4-style
   disclosure when drafting.
