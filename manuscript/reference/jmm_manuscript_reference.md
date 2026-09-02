**A provenance-aware workflow for strict receptor-preparation audits and
reference-pose recovery in molecular docking**

Andrés Monreal Hernández1\*, Sara Lizbeth Franco Amaya2, and Carlos
Ivanhoe Martínez Osorio3

1 Universidad Estatal de Sonora, Hermosillo, Sonora, Mexico. ORCID:
0009-0009-1207-8597

2 Doctorado en Nanotecnología, Universidad de Sonora, Hermosillo,
Sonora, Mexico. ORCID: 0009-0005-0272-0241

3 Doctorado en Ciencia de Materiales, Universidad de Sonora, Hermosillo,
Sonora, Mexico. ORCID: 0009-0003-7872-4965

\*Corresponding author: andres.monreal@ues.mx

# Abstract

Context: Molecular-docking workflows can become difficult to evaluate
when coordinate provenance, receptor-preparation choices, and
pose-selection rules are not retained. This work presents a
provenance-aware workflow for recording strict receptor-preparation
outcomes and bounded reference-pose recovery. The workflow evaluated 102
DUD-E receptor entries. Direct preparation of the original DUD-E PDB
representation failed for all 102 entries because blank element fields
were recorded. In separate official RCSB conditions, strict preparation
produced 26/102 successes from mmCIF and 27/102 successes from legacy
PDB. The two RCSB conditions agreed for 101 target-level outcomes; the
only difference was DPP4/2I78, a timeout in the frozen mmCIF audit that
completed in subsequent executions. Reference-pose campaigns retained
all nine Vina poses for BRAF 3D4Q/SM5 and KIF11 3CJO/K30. The top-score
and lowest-RMSD values were reported separately.

Methods: Python scripts recorded source URLs, checksums, command
outcomes, and normalized error classes. Strict Meeko preparation was run
without repair, alternate-location selection, template addition, or
post-outcome exclusion. For the two reference cases, receptor/cofactor
policies, ligand-component qualification, AutoDock Vina 1.2.7
configuration, and identity-mapped aligned heavy-atom RMSD were
specified before execution. Deterministic SVG figures and CSV
verification scripts support inspection of the versioned evidence.

Keywords: molecular docking; receptor preparation; reproducibility;
provenance; Meeko; AutoDock Vina.

# 1 Introduction

Molecular docking is frequently used to generate structural hypotheses
in computer-aided molecular research. Its practical interpretation
depends on coordinate provenance, treatment of the receptor and ligand,
and the rules used to retain and compare poses \[1-4\]. Benchmark
resources such as DUD-E and structural archives such as the RCSB Protein
Data Bank make it possible to define source inputs explicitly \[1,2\].
However, a workflow can still be difficult to audit if it does not
retain the preprocessing outcome and the reason for every failure.

This study addresses a narrower question than prospective virtual
screening: whether a declared direct receptor-preparation route can
process a fixed set of source structures, and how its outcomes can be
preserved for review. The work does not estimate affinity, biological
activity, selectivity, safety, or therapeutic effect. Instead, it
records a strict direct-preparation audit for 102 DUD-E targets, two
official RCSB source-and-reader conditions, and two reference-pose
recovery cases.

The workflow was designed so that negative results remain available
rather than being removed by undocumented repair. Retaining all output
poses is also important because the lowest-energy pose and the pose with
the lowest geometry-based RMSD need not be the same \[12,13\]. The
objective was therefore to provide a reusable, evidence-bounded software
workflow for source tracking, preparation compatibility, and
reference-pose reporting.

# 2 Materials and methods

## 2.1 Target registry and coordinate provenance

The frozen registry contains 102 DUD-E target/PDB identifier pairs. The
first condition retrieved the original receptor.pdb file distributed by
DUD-E for each target. A second condition used the corresponding
official RCSB PDBx/mmCIF file. A predeclared comparator retrieved the
official RCSB legacy-PDB file for the same 102 identifiers. For every
available source, the workflow retained the source URL, byte count, and
SHA-256 checksum in non-coordinate result tables. Raw coordinate files
were retained locally for computation but were not redistributed in the
repository.

## 2.2 Strict preparation audit

Each source was processed once by the declared Meeko pathway \[7\].
DUD-E and RCSB legacy-PDB files used the direct PDB reader; RCSB mmCIF
files used the declared ProDy reader \[6\]. Gemmi supported structural
processing in the workflow \[5\]. A direct success required return code
zero and a nonempty PDBQT output. All other outcomes were retained as
direct failures with a normalized error class. The protocol prohibited
coordinate repair, element correction, alternate-location selection,
residue deletion, manual template addition, and post-outcome exclusion.
The RCSB legacy-PDB comparison was selected before execution and was
evaluated as a distinct source-and-reader condition rather than as a
universal reference standard.

## 2.3 Reference-pose recovery cases

Two cases from the RCSB mmCIF direct-success subset were evaluated: BRAF
3D4Q with ligand SM5 and KIF11 3CJO with ligand K30 \[10,11\]. Ligand
component identity and a deposited ligand instance were qualified before
ligand preparation. The BRAF receptor consisted of author chain A
polymer. The KIF11 receptor retained author chain A polymer, deposited
ADP, and Mg2+; K30 and waters were excluded. The automatic Meeko
construction of an ADP template was retained as an explicit limitation.

AutoDock Vina 1.2.7 was run with CPU 1, seed 20260809, exhaustiveness
16, num_modes 9, minimum RMSD 1.0 Å, and energy range 3.0 kcal/mol
\[3,4\]. All nine poses from each case were retained. Pose recovery used
an identity-mapped, aligned heavy-atom RMSD calculation. Top-score and
lowest-RMSD poses were identified independently.

## 2.4 Reproducibility controls and computational assistance

The project uses versioned Python scripts, CSV tables, configuration
files, and deterministic SVG renderers. The verification script checks
denominators, outcome classes, pose-table consistency, and SHA-256
manifests. An independent direct DUD-E rerun was byte-identical to the
frozen result table. Independent reference-pose campaigns were
byte-identical to their frozen pose tables. The full mmCIF rerun agreed
for 101/102 targets; DPP4/2I78 differed because a frozen 300-s timeout
later completed successfully. Three predeclared follow-up attempts for
2I78 completed in 147.002-160.443 s and produced byte-identical PDBQT
files. The non-coordinate derived tables, protocols, and verification
scripts are supplied as Online Resource 1.

OpenAI Codex (GPT-5) was used as computational assistance during
software and document organization. The human authors reviewed the
computational evidence, figures, references, and manuscript material and
remain responsible for the methods, interpretation, and final submitted
text.

# 3 Results

## 3.1 Strict preparation outcomes

The original DUD-E PDB condition produced no nonempty PDBQT outputs: all
102 attempts failed, and the retained audit recorded blank PDB element
fields in all observed coordinate records. The RCSB mmCIF condition
produced 26 direct successes and 76 failures. The RCSB legacy-PDB
condition produced 27 direct successes and 75 failures (Table 1; Fig.
1).

**Table 1 Strict direct receptor-preparation outcomes by
source-and-reader condition**

  --------------------------------------------------------------------------
  **Condition**     **Direct          **Direct          **Scope**
                    successes**       failures**        
  ----------------- ----------------- ----------------- --------------------
  Original DUD-E    0/102             102/102           Frozen original
  PDB, direct                                           receptor
  reader                                                representation

  Official RCSB     26/102            76/102            Frozen
  mmCIF, ProDy                                          source-alternative
  reader                                                audit

  Official RCSB     27/102            75/102            Predeclared
  legacy PDB,                                           source-and-reader
  direct reader                                         comparator
  --------------------------------------------------------------------------

The RCSB mmCIF and legacy-PDB conditions agreed for 26 joint successes
and 75 joint failures. For the jointly failed records, all normalized
error classes agreed: 28 alternate-location cases, 29 residue-template
mismatches, 15 other preparation failures, and three element errors. The
single discordant outcome was DPP4/2I78, which was a timeout in the
frozen mmCIF table and a success in the later legacy-PDB condition. This
record is reported as execution-sensitive under the declared conditions,
not as an error correction to the frozen release.

## 3.2 Reference-pose recovery

All nine generated poses were retained for each reference case (Table 2;
Fig. 2). For BRAF 3D4Q/SM5, the top-score pose had a mapped heavy-atom
RMSD of 1.134 Å and the lowest observed RMSD was 1.083 Å. For KIF11
3CJO/K30, the top-score pose had RMSD 1.477 Å and the lowest observed
RMSD was 1.316 Å.

Table 2 Retained reference-pose recovery outcomes under declared
configurations

  -------------------------------------------------------------------------------
  **Target    **Retained   **Top-score    **RMSD of   **Lowest    **Declared
  and         poses**      value          top-score   RMSD (Å)**  limitation**
  ligand**                 (kcal/mol)**   pose (Å)**              
  ----------- ------------ -------------- ----------- ----------- ---------------
  BRAF 3D4Q / 9            -10.929        1.134       1.083       No post hoc
  SM5                                                             RMSD threshold
                                                                  was introduced

  KIF11 3CJO  9            -10.994        1.477       1.316       ADP template
  / K30                                                           was
                                                                  automatically
                                                                  constructed by
                                                                  Meeko
  -------------------------------------------------------------------------------

## 3.3 Evidence boundary

The audit and reference-pose outcomes are computational observations for
the specified sources, versions, and settings. They provide no estimate
of binding affinity or of biological, pharmacological, toxicological, or
clinical effects. The BRAF and KIF11 cases are deposited-structure
examples, not evidence of cancer-related efficacy.

# 4 Discussion

The results show why source identity and preparation behavior should be
reported before interpreting downstream docking calculations. The
original DUD-E representation could not be processed by the declared
direct route because its element fields were blank. In contrast, the two
official RCSB conditions produced 26 and 27 direct successes. Their
101/102 target-level agreement, including all normalized error classes
among shared failures, indicates that the observed compatibility pattern
persisted across these two specified source-and-reader conditions. It
does not demonstrate that the file formats are generally
interchangeable, because the source representation and reader both
differed.

The DPP4/2I78 result illustrates a separate reproducibility boundary. A
300-s timeout belongs to the frozen release record, while the subsequent
full rerun and three focused repetitions completed. Reporting both
results avoids silently replacing a versioned negative observation with
a later positive one. This practice is consistent with retaining
complete computational provenance and treating divergence as a result to
investigate rather than material to remove \[8,9\].

The BRAF and KIF11 campaigns demonstrate how a deposited reference
ligand, explicit receptor policy, fixed docking configuration, and
all-pose retention can be documented in a compact workflow. The
distinction between top-score and lowest-RMSD values is deliberate:
score ranking alone should not be presented as a geometric validation
result \[12,13\]. The two cases are insufficient to benchmark docking
broadly, but they provide inspectable reference-pose records that can be
extended through separately predeclared campaigns.

# 5 Conclusions

This software workflow preserves strict receptor-preparation outcomes,
provenance records, failure classes, and bounded reference-pose recovery
results. The versioned evidence identifies a direct-compatibility
limitation in the original DUD-E receptor files, documents two RCSB
source-and-reader conditions, and retains all poses for two reference
cases. Its appropriate use is reproducibility-oriented workflow
auditing; it must not be used to infer biological or therapeutic
effects.

# Figure captions

Fig. 1 Strict receptor-preparation outcomes across original DUD-E PDB,
official RCSB mmCIF, and official RCSB legacy-PDB conditions. Bars show
direct successes out of 102 entries. The RCSB conditions agreed for 101
target-level outcomes; DPP4/2I78 differed by a frozen timeout versus a
later success

![](media/image1.png){width="6.5in" height="2.364374453193351in"}

Fig. 2 Reference-pose recovery panel for BRAF 3D4Q/SM5 and KIF11
3CJO/K30. All nine Vina poses per case were retained. Top-score and
lowest-RMSD poses are reported under separate rules

![](media/image2.png){width="6.5in" height="2.4239588801399825in"}

# References

1\. Mysinger MM, Carchia M, Irwin JJ, Shoichet BK (2012) Directory of
Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better
Benchmarking. J Med Chem 55:6582-6594. https://doi.org/10.1021/jm300687e

2\. Burley SK, Bhikadiya C, Bi C, Bittrich S, Chao H, Chen L, Craig PA,
Crichlow GV, Dalenberg K, Duarte JM, Dutta S, Fayazi M, Feng Z, Flatt
JW, Ganesan S, Ghosh S, Goodsell DS, Green RK, Guranovic V, Henry J,
Hudson BP, Khokhriakov I, Lawson CL, Liang Y, Lowe R, Peisach E,
Persikova I, Piehl DW, Rose Y, Sali A, Segura J, Sekharan M, Shao C,
Vallat B, Voigt M, Webb B, Westbrook JD, Whetstone S, Young JY, Zalevsky
A, Zardecki C (2023) RCSB Protein Data Bank (RCSB.org): delivery of
experimentally-determined PDB structures alongside one million computed
structure models of proteins from artificial intelligence/machine
learning. Nucleic Acids Res 51:D488-D508.
https://doi.org/10.1093/nar/gkac1077

3\. Trott O, Olson AJ (2010) AutoDock Vina: improving the speed and
accuracy of docking with a new scoring function, efficient optimization,
and multithreading. J Comput Chem 31:455-461.
https://doi.org/10.1002/jcc.21334

4\. Eberhardt J, Santos-Martins D, Tillack AF, Forli S (2021) AutoDock
Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python
Bindings. J Chem Inf Model 61:3891-3898.
https://doi.org/10.1021/acs.jcim.1c00203

5\. Wojdyr M (2022) GEMMI: A library for structural biology. J Open
Source Softw 7:4200. https://doi.org/10.21105/joss.04200

6\. Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K,
Doruker P, Li H, Bahar I (2021) ProDy 2.0: increased scale and scope
after 10 years of protein dynamics modelling with Python. Bioinformatics
37:3657-3659. https://doi.org/10.1093/bioinformatics/btab187

7\. Forli Lab (2025) Meeko 0.7.1. Python Package Index.
https://pypi.org/project/meeko/0.7.1/

8\. Wilkinson MD, Dumontier M, Aalbersberg IJ, Appleton G, Axton M, Baak
A, Blomberg N, Boiten J-W, da Silva Santos LB, Bourne PE, Bouwman J,
Brookes AJ, Clark T, Crosas M, Dillo I, Dumon O, Edmunds S, Evelo CT,
Finkers R, Gonzalez-Beltran A, Gray AJG, Groth P, Goble C, Grethe JS,
Heringa J, \'t Hoen PAC, Hooft R, Kuhn T, Kok R, Kok J, Lusher SJ,
Martone ME, Mons A, Packer AL, Persson B, Rocca-Serra P, Roos M, van
Schaik R, Sansone S-A, Schultes E, Sengstag T, Slater T, Strawn G,
Swertz MA, Thompson M, van der Lei J, van Mulligen E, Velterop J,
Waagmeester A, Wittenburg P, Wolstencroft K, Zhao J, Mons B (2016) The
FAIR Guiding Principles for scientific data management and stewardship.
Sci Data 3:160018. https://doi.org/10.1038/sdata.2016.18

9\. Sandve GK, Nekrutenko A, Taylor J, Hovig E (2013) Ten Simple Rules
for Reproducible Computational Research. PLoS Comput Biol 9:e1003285.
https://doi.org/10.1371/journal.pcbi.1003285

10\. Hansen JD, Grina J, Newhouse B, Welch M, Topalov G, Littman N,
Callejo M, Gloor S, Martinson M, Laird E, Brandhuber BJ, Vigers G,
Morales T, Woessner R, Randolph N, Lyssikatos J, Olivero A (2008) Potent
and selective pyrazole-based inhibitors of B-Raf kinase. Bioorg Med Chem
Lett 18:4692-4695. https://doi.org/10.1016/j.bmcl.2008.07.002

11\. Cox CD, Coleman PJ, Breslin MJ, Whitman DB, Garbaccio RM, Fraley
ME, Buser CA, Walsh ES, Hamilton K, Schaber MD, Lobell RB, Tao W, Davide
JP, Diehl RE, Abrams MT, South VJ, Huber HE, Torrent M, Prueksaritanont
T, Li C, Slaughter DE, Mahan E, Fernandez-Metzler C, Yan Y, Kuo LC, Kohl
NE, Hartman GD (2008) Kinesin spindle protein (KSP) inhibitors. 9.
Discovery of MK-0731 for the treatment of taxane-refractory cancer. J
Med Chem 51:4239-4252. https://doi.org/10.1021/jm800386y

12\. Mukherjee S, Balius TE, Rizzo RC (2010) Docking validation
resources: protein family and ligand flexibility experiments. J Chem Inf
Model 50:1986-2000. https://doi.org/10.1021/ci1001982

13\. Ramírez D, Caballero J (2018) Is It Reliable to Take the Molecular
Docking Top Scoring Position as the Best Solution without Considering
Available Structural Data? Molecules 23:1038.
https://doi.org/10.3390/molecules23051038

14\. Allen WJ, Rizzo RC (2014) Implementation of the Hungarian algorithm
to account for ligand symmetry and similarity in structure-based design.
J Chem Inf Model 54:518-529. https://doi.org/10.1021/ci400534h

15\. Stanzione F, Giangreco I, Cole JC (2021) Use of molecular docking
computational tools in drug discovery. Prog Med Chem 60:273-343.
https://doi.org/10.1016/bs.pmch.2021.01.004

# Statements and declarations

## Funding

The authors declare that no funds, grants, or other support were
received during the preparation of this manuscript. The work used
available local computing resources.

## Competing interests

The authors have no relevant financial or non-financial interests to
disclose.

## Author contributions

Sara Lizbeth Franco Amaya: Conceptualization; Writing - review and
editing. Carlos Ivanhoe Martínez Osorio: Methodology. Andrés Monreal
Hernández: Data curation; Formal analysis; Investigation; Software;
Visualization; Writing - original draft. All authors reviewed the
evidence, contributed to manuscript review, and must approve the final
submitted version.

## Data availability

Versioned code, derived non-coordinate result tables, protocols,
validation artefacts, and deterministic figures are available from the
public project repository. Supporting non-coordinate derived tables,
protocols, and verification scripts are supplied as Online Resource 1.
The frozen v0.2.0 evidence archive is available at
https://doi.org/10.5281/zenodo.21866318. Original coordinates and local
PDBQT/log files are not redistributed; source identifiers, URLs, and
checksums are provided to support retrieval and verification.
