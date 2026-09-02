# DockPrep Audit: study design for a software-and-methods article

## Working question

Do objectively detectable receptor-structure features (alternate locations,
retained waters, metal ions, non-standard residues, incomplete annotations and
duplicate atom records) predict differences in the reproducibility of a
molecular-docking workflow?

The study does **not** assume that one preparation choice is universally
correct. It measures when an explicit choice is necessary and whether results
change across documented alternatives.

## Separation from previous manuscripts

The evaluation cohort must be independent of the 30-PDB alternate-conformer
study and the 102-target DUD-E workflow audit already submitted. No aggregate
result from those cohorts will be reused as a new result. They may be cited as
motivation only.

## Cohort plan

Create a preregistered, public PDB cohort of 60 protein–ligand complexes:

| Stratum | Structures | Purpose |
|---|---:|---|
| Alternate locations | 15 | Test explicit conformer-selection decisions. |
| Waters in/near binding site | 15 | Test retain-versus-remove policy effects. |
| Metals/cofactors | 15 | Test documented treatment of coordination environments. |
| Low-risk controls | 15 | Establish a baseline with no high-severity audit finding. |

Eligibility: experimental protein–ligand complex, a resolvable reference
ligand, binding-site coordinates usable for redocking, an open PDB record, and
no overlap with the submitted-study cohorts. Selection and exclusions will be
stored in `benchmark/cohort.csv` before docking begins.

## Experimental conditions

For every eligible complex, run fixed-grid redocking under three fully
documented preparation policies:

1. **Reference-preserving:** retain source annotations and explicitly select
   the stated alternate conformer when one exists.
2. **Conventional simplified:** apply a clearly specified solvent/heteroatom
   removal policy.
3. **Tool-default:** use a named open preparation pipeline with its version and
   defaults recorded.

The docking engine, search box, exhaustiveness, seed, ligand protonation
assumptions and scoring settings are held fixed within each comparison. Every
run produces a provenance manifest.

## Primary outcomes

- Redocking success rate, using ligand heavy-atom RMSD at or below 2 Å.
- Difference in RMSD between preparation policies.
- Rank/pose agreement across policies.
- Prevalence of each audit finding and its association with divergent output.

## Secondary outcomes

- Human-review agreement for each audit rule on a blinded subset.
- Time required to identify and document a preparation decision.
- False-positive and false-negative rule counts against the curated annotations.

## Statistical plan

Report effect sizes with confidence intervals rather than relying solely on
threshold tests. Pairwise preparation comparisons use paired analyses within
complex. The low-risk control stratum is analysed separately to determine
whether differences are concentrated in flagged structures.

## Proposed figures

1. Graphical abstract: PDB record → audit → explicit decision → docking result.
2. Study-flow diagram: eligibility, stratification, exclusions and final cohort.
3. Architecture and provenance diagram of DockPrep Audit.
4. Frequency and co-occurrence of audit findings.
5. Heatmap of complex × audit finding × preparation policy.
6. Redocking success and RMSD distribution across policies.
7. Three molecular case studies (alternate location, water, metal).
8. Agreement between automated findings and blinded expert review.
9. Decision tree translating results into reproducible preparation guidance.

## Reproducibility package

Release the cohort manifest, raw audit reports, preparation manifests,
command logs, docking outputs, figure scripts and an environment lockfile. A
versioned Zenodo archive will receive a DOI at each public release.

## Milestones

1. Freeze cohort eligibility and create the public manifest.
2. Implement the remaining audit rules and visual report summaries.
3. Curate and annotate the 60 structures.
4. Execute the preregistered redocking benchmark.
5. Write the article after the results are fixed.
