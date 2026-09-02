"""F4: symmetry-corrected RMSD (PROJECT-ROADMAP.md section 19, item F4),
added 2026-08-28.

Recomputes RMSD for every already-completed redocking run (Fase C's 360,
F2's 90, F3's 90 -- 540 total, no new docking needed) using molecular
automorphism search instead of the fixed exact-coordinate index map alone,
so a ligand with locally symmetric groups (e.g. a carboxylate's two
oxygens, a phenyl ring) is not penalized for an arbitrary atom-labeling
choice that has no physical meaning.

Method, deliberately NOT a rigid-body realignment (which would defeat the
purpose of a docking-pose RMSD -- a badly placed but correctly shaped pose
must still score badly):
  1. ref_coords <-> ligand.pdbqt index map is built exactly as in
     smoke_redock_case.py (exact-coordinate match at the pre-docking pose,
     0.01 A tolerance) -- this establishes physical atom identity, not
     symmetry.
  2. ligand.sdf (built during preparation from the same deposited
     coordinates, via AssignBondOrdersFromTemplate + AddHs, which appends
     new H atoms after the original heavy atoms without reordering them)
     gives the heavy-atom molecular graph in the SAME order as ref_coords.
  3. mol.GetSubstructMatches(mol, uniquify=False, useChirality=True) on the
     heavy-atom-only graph enumerates every valid automorphism (atom
     permutation that preserves the graph, respecting stereochemistry).
  4. For each docked pose, the symmetry-corrected RMSD is the minimum over
     all automorphisms of the plain (unaligned) RMSD between ref_coords and
     the pose's coordinates permuted by that automorphism.

Writes benchmark/symmetry-rmsd/{phase_c60,f2_altloc,f3_metal}_symmetry_rmsd.csv
with both the original and symmetry-corrected RMSD per (case, policy, seed,
rank), plus a summary of how many cases/outcomes changed.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from rdkit import Chem  # noqa: E402

from smoke_redock_case import (  # noqa: E402
    build_reference_to_pdbqt_index_map,
    count_models,
    read_all_atom_coords_pdbqt,
    read_heavy_atom_coords_pdb,
    rmsd,
)

SUCCESS_THRESHOLD_A = 2.0


def get_automorphisms(sdf_path: Path, n_heavy: int) -> list[tuple[int, ...]]:
    """Heavy-atom automorphisms of the ligand, in the SDF's own atom order
    (which matches ref_coords' order for the first n_heavy atoms -- see
    module docstring point 2)."""
    mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"RDKit could not read {sdf_path}")
    mol_noH = Chem.RemoveHs(mol, sanitize=False)
    if mol_noH.GetNumAtoms() != n_heavy:
        raise ValueError(
            f"{sdf_path}: heavy-atom count mismatch after RemoveHs "
            f"({mol_noH.GetNumAtoms()} vs expected {n_heavy})"
        )
    matches = mol_noH.GetSubstructMatches(mol_noH, uniquify=False, useChirality=True, maxMatches=100000)
    if not matches:
        matches = (tuple(range(n_heavy)),)
    return matches


def symmetry_rmsd(ref_coords, pose_coords, automorphisms) -> tuple[float, tuple[int, ...]]:
    """Minimum RMSD over all automorphisms; automorphism maps output-position
    -> input-position (RDKit substruct-match convention), so permuted[i] =
    pose_coords[automorphism[i]]."""
    best = (float("inf"), None)
    for perm in automorphisms:
        permuted = [pose_coords[p] for p in perm]
        d = rmsd(ref_coords, permuted)
        if d < best[0]:
            best = (d, perm)
    return best


def find_case_files(run_dir: Path):
    ligand_raw = run_dir / "ligand_raw.pdb"
    ligand_sdf = run_dir / "ligand.sdf"
    ligand_pdbqt = run_dir / "ligand.pdbqt"
    redock_out = run_dir / "redock_out.pdbqt"
    if all(p.exists() for p in (ligand_raw, ligand_sdf, ligand_pdbqt, redock_out)):
        return ligand_raw, ligand_sdf, ligand_pdbqt, redock_out
    return None


def resolve_run_dir(roots: list[Path], policy_dir: str, seed: int, pdb_id: str) -> Path | None:
    """The original 12-case pilot cohort's intermediate files live under
    benchmark/phase-c/ (run before the 48-case expansion, kept under
    phase-c-60/); every other case/comparison lives under its own single
    root. Search all given roots in order and return the first that has a
    complete file set."""
    for root in roots:
        candidate = root / policy_dir / f"seed{seed}" / pdb_id
        if find_case_files(candidate) is not None:
            return candidate
    return None


def classify(rmsds_by_rank):
    if rmsds_by_rank[0] <= SUCCESS_THRESHOLD_A:
        return "success"
    if any(r <= SUCCESS_THRESHOLD_A for r in rmsds_by_rank):
        return "scoring_fail"
    return "sampling_fail"


def process_run_dir(run_dir: Path, automorphism_cache: dict) -> dict | None:
    files = find_case_files(run_dir)
    if files is None:
        return None
    ligand_raw, ligand_sdf, ligand_pdbqt, redock_out = files

    ref_coords = read_heavy_atom_coords_pdb(ligand_raw)
    ligand_pdbqt_coords = read_all_atom_coords_pdbqt(ligand_pdbqt, model=None)
    index_map = build_reference_to_pdbqt_index_map(ref_coords, ligand_pdbqt_coords)

    cache_key = str(ligand_sdf)
    if cache_key not in automorphism_cache:
        automorphism_cache[cache_key] = get_automorphisms(ligand_sdf, len(ref_coords))
    automorphisms = automorphism_cache[cache_key]

    n_poses = count_models(redock_out)
    orig_rmsds, sym_rmsds = [], []
    for rank in range(1, n_poses + 1):
        pose_all_atoms = read_all_atom_coords_pdbqt(redock_out, model=rank)
        pose_coords = [pose_all_atoms[j] for j in index_map]
        orig_rmsds.append(round(rmsd(ref_coords, pose_coords), 3))
        sym_d, _ = symmetry_rmsd(ref_coords, pose_coords, automorphisms)
        sym_rmsds.append(round(sym_d, 3))

    return {
        "n_heavy_atoms": len(ref_coords),
        "n_automorphisms": len(automorphisms),
        "orig_rmsd_by_rank": orig_rmsds,
        "sym_rmsd_by_rank": sym_rmsds,
        "orig_top1": orig_rmsds[0],
        "orig_best": min(orig_rmsds),
        "sym_top1": sym_rmsds[0],
        "sym_best": min(sym_rmsds),
        "orig_outcome": classify(orig_rmsds),
        "sym_outcome": classify(sym_rmsds),
    }


CASES = [r["pdb_id"] for r in csv.DictReader(open(REPO_ROOT / "benchmark" / "pilot_cases.csv", encoding="utf-8"))]
SEEDS = [42, 123, 2024]
OUT_DIR = REPO_ROOT / "benchmark" / "symmetry-rmsd"


def run_comparison(name: str, roots: list[Path], policies: list[str]) -> list[dict]:
    automorphism_cache: dict[str, list] = {}
    rows = []
    missing = []
    for pdb_id in CASES:
        for policy in policies:
            for seed in SEEDS:
                run_dir = resolve_run_dir(roots, policy, seed, pdb_id)
                if run_dir is None:
                    missing.append((pdb_id, policy, seed))
                    continue
                try:
                    result = process_run_dir(run_dir, automorphism_cache)
                except Exception as exc:
                    print(f"  ERROR {pdb_id} {policy} seed={seed}: {exc}")
                    missing.append((pdb_id, policy, seed))
                    continue
                rows.append({"pdb_id": pdb_id, "policy": policy, "seed": seed, **result,
                             "orig_rmsd_by_rank": ";".join(str(x) for x in result["orig_rmsd_by_rank"]),
                             "sym_rmsd_by_rank": ";".join(str(x) for x in result["sym_rmsd_by_rank"])})
    print(f"[{name}] {len(rows)} runs processed, {len(missing)} missing/failed")
    if missing:
        print(f"  missing: {missing[:10]}{' ...' if len(missing) > 10 else ''}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_DIR / f"{name}_symmetry_rmsd.csv"
    if rows:
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"  wrote {out_csv}")
    return rows


def summarize(name: str, rows: list[dict]) -> None:
    n_sym_gt = sum(1 for r in rows if r["n_automorphisms"] > 1)
    delta_top1 = [r["orig_top1"] - r["sym_top1"] for r in rows]
    delta_best = [r["orig_best"] - r["sym_best"] for r in rows]
    outcome_changed = sum(1 for r in rows if r["orig_outcome"] != r["sym_outcome"])
    big_top1 = [(r["pdb_id"], r["policy"], r["seed"], round(d, 3)) for r, d in zip(rows, delta_top1) if d > 0.3]
    print(f"\n=== {name} summary ({len(rows)} runs) ===")
    print(f"  runs with >1 automorphism (symmetry could matter): {n_sym_gt}/{len(rows)}")
    print(f"  mean delta top1 (orig-sym): {sum(delta_top1)/len(delta_top1):.4f} A "
          f"(max {max(delta_top1):.3f}, n>0.3A: {sum(1 for d in delta_top1 if d > 0.3)})")
    print(f"  mean delta best (orig-sym): {sum(delta_best)/len(delta_best):.4f} A "
          f"(max {max(delta_best):.3f}, n>0.3A: {sum(1 for d in delta_best if d > 0.3)})")
    print(f"  outcome category changed (success/scoring_fail/sampling_fail): {outcome_changed}/{len(rows)}")
    if big_top1:
        print(f"  runs with >0.3A top1 correction: {big_top1[:15]}")


def main() -> None:
    water_rows = run_comparison(
        "phase_c60",
        [REPO_ROOT / "benchmark" / "phase-c-60", REPO_ROOT / "benchmark" / "phase-c"],
        ["conservative_water", "simplified_no_water"],
    )
    altloc_rows = run_comparison(
        "f2_altloc",
        [REPO_ROOT / "benchmark" / "phase-f2-altloc"],
        ["highest_occupancy", "lowest_occupancy"],
    )
    metal_rows = run_comparison(
        "f3_metal",
        [REPO_ROOT / "benchmark" / "phase-f3-metal"],
        ["retain", "remove"],
    )
    if water_rows:
        summarize("Water policy (Fase C, 60 cases)", water_rows)
    if altloc_rows:
        summarize("AltLoc policy (F2, 15 cases)", altloc_rows)
    if metal_rows:
        summarize("Metal policy (F3, 15 cases)", metal_rows)


if __name__ == "__main__":
    main()
