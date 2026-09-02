# DockPrep Audit

**DockPrep Audit** is an open-source, local-first quality-control tool for receptor structures used in molecular docking. It does not alter a structure or claim that it is biologically correct. Instead, it makes preparation choices visible before a docking run.

## What v0.2 audits

Whole-structure findings (any occurrence anywhere in the receptor):

- Alternate atom locations (`ALTLOC_PRESENT`): a receptor requires an explicit conformer policy.
- Water molecules, non-standard residues, and metal atoms (`WATERS_PRESENT`, `METAL_PRESENT`): potential protocol-dependent decisions.
- Missing occupancy and missing element annotations.
- Duplicate atom records within a residue.

Site-local findings (scoped to the declared ligand's binding site, added in v0.2 — require a `ligand` argument):

- `SITE_BRIDGING_WATER_PRESENT`: a water within hydrogen-bond range of both the ligand and the receptor.
- `SITE_ALTLOC_PRESENT`: a receptor alternate-location conformer near the binding site (excludes disordered waters and the ligand's own atoms).
- `SITE_METAL_PRESENT`: a metal ion near the binding site.

The output is a machine-readable JSON report and a portable HTML report. The structure stays local; no data are uploaded.

## Companion study

`manuscript/draft.md` and `benchmark/` contain a 60-case redocking benchmark testing whether these findings actually predict preparation-policy sensitivity (water retention, altLoc resolution, metal retention), including a 76-case independent validation cohort and a formal predictive evaluation. `PROJECT-ROADMAP.md` section 19 documents the full methodology and results (items F1-F9).

## Quick start

```powershell
py -m pip install -e .
dockprep-audit audit path\to\receptor.pdb --output audit-results
```

Open `audit-results/report.html` in a browser. Exit status is `0` when no high-severity finding is present and `2` otherwise.

## Scientific scope

This is a preparation-audit tool, not a docking engine and not an automatic repair program. Its intended contribution is traceability: before any repair, deletion, protonation, or conformer selection, the user receives an explicit report of the decisions that can affect a docking workflow.

## Status

1. ~~Validate rules against curated DUD-E/PDB receptor examples.~~ Done (102-target audit, 60-case redocking benchmark).
2. ~~Quantify rule prevalence and inter-protocol disagreement.~~ Done (`manuscript/draft.md` Sections 3.5-3.10).
3. Add a graphical interface and exportable provenance manifest. Not started.
4. Archive a release in Zenodo and submit the manuscript. Manuscript drafted, pending final author/DOI details (`manuscript/draft.md`, `[DATA:]` markers); Zenodo archive and GitHub repo not yet created.

## License

MIT. See [LICENSE](LICENSE).
