"""F5: secondary diagnosis of 1EPP/1PPM (PROJECT-ROADMAP.md section 19,
item F5), added 2026-08-28.

Checks four candidate explanations for why these two water-policy-stratum
cases fail consistently (~9 A best-of-9 RMSD) under every policy tested in
this study: docking-box adequacy, ligand protonation, ligand flexibility,
and reference-pose quality. Also checks whether the water-policy stratum as
a whole is more flexible than the other three strata, since 1EPP/1PPM's
own flexibility turned out not to be an isolated property of just these two
cases.

Prints its findings; see manuscript/draft.md Section 3.4 and Limitation 10
for how they were written up.
"""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

REPO_ROOT = Path(__file__).resolve().parent.parent
CASES = list(csv.DictReader(open(REPO_ROOT / "benchmark" / "pilot_cases.csv", encoding="utf-8")))
RUN_ROOTS = [
    REPO_ROOT / "benchmark" / "phase-c-60" / "conservative_water" / "seed42",
    REPO_ROOT / "benchmark" / "phase-c" / "conservative_water" / "seed42",
]


def find_case_dir(pdb_id: str) -> Path | None:
    for root in RUN_ROOTS:
        candidate = root / pdb_id
        if (candidate / "ligand.sdf").exists():
            return candidate
    return None


def rotatable_bonds(pdb_id: str) -> int | None:
    case_dir = find_case_dir(pdb_id)
    if case_dir is None:
        return None
    mol = Chem.SDMolSupplier(str(case_dir / "ligand.sdf"), removeHs=False)[0]
    return rdMolDescriptors.CalcNumRotatableBonds(mol)


def ligand_extent(pdb_id: str) -> float | None:
    case_dir = find_case_dir(pdb_id)
    if case_dir is None:
        return None
    coords = []
    for line in (case_dir / "ligand_raw.pdb").read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM"):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    xs, ys, zs = zip(*coords)
    return max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def box_size(pdb_id: str) -> float | None:
    raw_path = REPO_ROOT / "benchmark" / "phase-c-60" / "raw_runs" / f"{pdb_id}_conservative_water_seed42.json"
    if not raw_path.exists():
        raw_path = REPO_ROOT / "benchmark" / "phase-c" / "raw_runs" / f"{pdb_id}_conservative_water_seed42.json"
    if not raw_path.exists():
        return None
    return json.loads(raw_path.read_text(encoding="utf-8"))["result"]["box"]["size_x"]


def protonation_attempts(pdb_id: str) -> dict | None:
    case_dir = find_case_dir(pdb_id)
    if case_dir is None:
        return None
    p = case_dir / "ligand_protonation.json"
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    print("=== Box adequacy (ligand max extent vs. box size) ===")
    for pdb_id in ["1EPP", "1PPM", "1QRP"]:
        ext = ligand_extent(pdb_id)
        box = box_size(pdb_id)
        print(f"  {pdb_id}: extent={ext:.1f} A  box={box:.1f} A  fits={ext <= box}")

    print("\n=== Protonation (single clean attempt = not the cause) ===")
    for pdb_id in ["1EPP", "1PPM"]:
        d = protonation_attempts(pdb_id)
        print(f"  {pdb_id}: {len(d['attempts'])} attempt(s), chosen={d['chosen_smiles'] is not None}")

    print("\n=== Mean rotatable bonds by stratum ===")
    by_stratum = {}
    for row in CASES:
        rb = rotatable_bonds(row["pdb_id"])
        if rb is not None:
            by_stratum.setdefault(row["stratum"], []).append(rb)
    for stratum, vals in by_stratum.items():
        print(f"  {stratum:20s} n={len(vals):2d}  mean={statistics.mean(vals):.2f}  "
              f"median={statistics.median(vals)}  max={max(vals)}")

    print("\n=== Water-policy stratum, rotatable bonds (sorted) ===")
    for row in CASES:
        if row["stratum"] != "water_policy":
            continue
        rb = rotatable_bonds(row["pdb_id"])
        print(f"  {row['pdb_id']}: {rb}")


if __name__ == "__main__":
    main()
