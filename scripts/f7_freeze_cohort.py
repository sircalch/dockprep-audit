"""F7 step 2: from the screened candidates (f7_screen_candidates.py), select
and freeze the external validation cohort.

Pre-registered selection rule, fixed before this script is run and not
adjusted afterward: take every QUALIFIES candidate, sorted by PDB ID for a
deterministic order (no outcome-based cherry-picking, since no docking has
been run yet at this stage), applying a cap of at most 2 cases per
ligand-component-id (a cheap proxy for "same target/series", avoiding the
UniProt lookup cost the main cohort used) to avoid a handful of chemical
series dominating the cohort the way the main cohort avoided one dominant
biological target (Section 2.2 of the manuscript). Stops once the target N
is reached or candidates are exhausted.

Writes benchmark/f7-external-validation/manifest.csv (frozen: pdb_id,
ligand_component_id/chain/resseq, receptor_chains, source_sha256) and
copies each selected structure's raw PDB alongside a checksum file, mirroring
the provenance discipline already used for the main 60-case cohort
(Section 2.2).
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "benchmark" / "f7-external-validation"
RAW_DIR = OUT_DIR / "raw-pdb"
TARGET_N = 80
MAX_PER_LIGAND = 2


def main() -> None:
    rows = list(csv.DictReader(open(OUT_DIR / "screening_results.csv", encoding="utf-8")))
    qualifying = [r for r in rows if r["status"] == "QUALIFIES"]
    qualifying.sort(key=lambda r: r["pdb_id"])
    print(f"{len(qualifying)} qualifying candidates available")

    selected = []
    per_ligand_count: dict[str, int] = {}
    for r in qualifying:
        lig = r["ligand_component_id"]
        if per_ligand_count.get(lig, 0) >= MAX_PER_LIGAND:
            continue
        selected.append(r)
        per_ligand_count[lig] = per_ligand_count.get(lig, 0) + 1
        if len(selected) >= TARGET_N:
            break

    print(f"{len(selected)} selected after per-ligand cap (max {MAX_PER_LIGAND} per component id)")

    manifest_rows = []
    for r in selected:
        pdb_id = r["pdb_id"]
        raw_path = RAW_DIR / f"{pdb_id}.pdb"
        sha256 = hashlib.sha256(raw_path.read_bytes()).hexdigest()
        manifest_rows.append({
            "pdb_id": pdb_id,
            "ligand_component_id": r["ligand_component_id"],
            "ligand_chain": r["ligand_chain"],
            "ligand_resseq": r["ligand_resseq"],
            "receptor_chains": r["ligand_chain"],  # single-protein-chain candidates by search construction
            "source_url": f"https://files.rcsb.org/download/{pdb_id}.pdb",
            "source_sha256": sha256,
            "findings": r["findings"],
        })

    with (OUT_DIR / "manifest.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
        w.writeheader()
        w.writerows(manifest_rows)
    print(f"Wrote {OUT_DIR / 'manifest.csv'} ({len(manifest_rows)} cases)")

    # A single combined-manifest checksum, mirroring the main cohort's practice
    combined = "\n".join(f"{r['pdb_id']},{r['source_sha256']}" for r in manifest_rows)
    manifest_sha256 = hashlib.sha256(combined.encode("utf-8")).hexdigest()
    (OUT_DIR / "manifest.sha256.txt").write_text(manifest_sha256 + "\n", encoding="utf-8")
    print(f"Manifest checksum: {manifest_sha256}")


if __name__ == "__main__":
    main()
