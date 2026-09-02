"""F7 step 5: symmetry-corrected RMSD + aggregation for the external
validation cohort, and comparison against the main 60-case cohort's
water-policy stratum finding (Section 3.5) -- the discovery/validation
check the whole of F7 exists to run. No definition (success threshold,
Wilcoxon test, symmetry-correction method, per-case effect threshold) is
adjusted from what Sections 2.4/3.5 already fixed before this cohort was
even screened.
"""

from __future__ import annotations

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import wilcoxon

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from compute_symmetry_rmsd import get_automorphisms, symmetry_rmsd  # noqa: E402
from smoke_redock_case import (  # noqa: E402
    build_reference_to_pdbqt_index_map,
    count_models,
    read_all_atom_coords_pdbqt,
    read_heavy_atom_coords_pdb,
    rmsd,
)

OUT_DIR = REPO_ROOT / "benchmark" / "f7-external-validation"
RUNS_DIR = OUT_DIR / "runs"
SEEDS = [42, 123, 2024]
SUCCESS_THRESHOLD_A = 2.0
EFFECT_THRESHOLD_A = 0.3


def find_case_files(run_dir: Path):
    files = [run_dir / n for n in ("ligand_raw.pdb", "ligand.sdf", "ligand.pdbqt", "redock_out.pdbqt")]
    return tuple(files) if all(p.exists() for p in files) else None


def classify(rmsds_by_rank):
    if rmsds_by_rank[0] <= SUCCESS_THRESHOLD_A:
        return "success"
    if any(r <= SUCCESS_THRESHOLD_A for r in rmsds_by_rank):
        return "scoring_fail"
    return "sampling_fail"


def process_run(run_dir: Path, automorphism_cache: dict) -> dict | None:
    files = find_case_files(run_dir)
    if files is None:
        return None
    ligand_raw, ligand_sdf, ligand_pdbqt, redock_out = files

    ref_coords = read_heavy_atom_coords_pdb(ligand_raw)
    ligand_pdbqt_coords = read_all_atom_coords_pdbqt(ligand_pdbqt, model=None)
    try:
        index_map = build_reference_to_pdbqt_index_map(ref_coords, ligand_pdbqt_coords)
    except ValueError:
        return None

    cache_key = str(ligand_sdf)
    if cache_key not in automorphism_cache:
        try:
            automorphism_cache[cache_key] = get_automorphisms(ligand_sdf, len(ref_coords))
        except Exception:
            automorphism_cache[cache_key] = None
    automorphisms = automorphism_cache[cache_key]
    if automorphisms is None:
        return None

    n_poses = count_models(redock_out)
    sym_rmsds = []
    for rank in range(1, n_poses + 1):
        pose_all_atoms = read_all_atom_coords_pdbqt(redock_out, model=rank)
        pose_coords = [pose_all_atoms[j] for j in index_map]
        d, _ = symmetry_rmsd(ref_coords, pose_coords, automorphisms)
        sym_rmsds.append(round(d, 3))

    return {"sym_top1": sym_rmsds[0], "sym_best": min(sym_rmsds), "outcome": classify(sym_rmsds)}


def main() -> None:
    manifest_rows = list(csv.DictReader(open(OUT_DIR / "manifest.csv", encoding="utf-8")))
    cases = [r["pdb_id"] for r in manifest_rows]

    automorphism_cache: dict[str, list] = {}
    case_rows = []
    skipped = []
    for pdb_id in cases:
        for policy in ["conservative_water", "simplified_no_water"]:
            top1s, bests, outcomes = [], [], []
            for seed in SEEDS:
                run_dir = RUNS_DIR / policy / f"seed{seed}" / pdb_id
                result = process_run(run_dir, automorphism_cache)
                if result is None:
                    continue
                top1s.append(result["sym_top1"])
                bests.append(result["sym_best"])
                outcomes.append(result["outcome"])
            if len(top1s) < 3:
                skipped.append((pdb_id, policy, len(top1s)))
                continue
            case_rows.append({
                "pdb_id": pdb_id, "policy": policy,
                "median_top1": statistics.median(top1s), "median_best": statistics.median(bests),
                "unstable": len(set(outcomes)) > 1,
            })

    if skipped:
        print(f"{len(skipped)} (case, policy) pairs skipped (incomplete/unusable data): {skipped}")

    by_case = defaultdict(dict)
    for r in case_rows:
        by_case[r["pdb_id"]][r["policy"]] = r
    complete_cases = [pdb for pdb, d in by_case.items()
                      if "conservative_water" in d and "simplified_no_water" in d]
    print(f"\n{len(complete_cases)} cases with both policies complete (of {len(cases)} in manifest)")

    cons_best = [by_case[c]["conservative_water"]["median_best"] for c in complete_cases]
    simp_best = [by_case[c]["simplified_no_water"]["median_best"] for c in complete_cases]
    cons_top1 = [by_case[c]["conservative_water"]["median_top1"] for c in complete_cases]
    simp_top1 = [by_case[c]["simplified_no_water"]["median_top1"] for c in complete_cases]

    succ_cons_best = sum(1 for x in cons_best if x <= SUCCESS_THRESHOLD_A) / len(cons_best)
    succ_simp_best = sum(1 for x in simp_best if x <= SUCCESS_THRESHOLD_A) / len(simp_best)
    succ_cons_top1 = sum(1 for x in cons_top1 if x <= SUCCESS_THRESHOLD_A) / len(cons_top1)
    succ_simp_top1 = sum(1 for x in simp_top1 if x <= SUCCESS_THRESHOLD_A) / len(simp_top1)

    print(f"\n=== F7 external validation cohort, symmetry-corrected RMSD (n={len(complete_cases)}) ===")
    print(f"  Conservative: succ(top1)={succ_cons_top1:.3f}  RMSD(top1)={statistics.mean(cons_top1):.3f}  "
          f"succ(best9)={succ_cons_best:.3f}  RMSD(best9)={statistics.mean(cons_best):.3f}")
    print(f"  Simplified:   succ(top1)={succ_simp_top1:.3f}  RMSD(top1)={statistics.mean(simp_top1):.3f}  "
          f"succ(best9)={succ_simp_best:.3f}  RMSD(best9)={statistics.mean(simp_best):.3f}")

    stat, p = wilcoxon(cons_best, simp_best)
    mean_diff = statistics.mean(c - s for c, s in zip(cons_best, simp_best))
    print(f"\nPaired Wilcoxon (best-of-9 RMSD, conservative vs simplified, n={len(complete_cases)}): "
          f"W={stat:.1f}, p={p:.4f}, mean paired diff (cons-simp)={mean_diff:+.3f} A")

    print(f"\nMain-cohort comparison (Table 3, water_policy, n=11): succ(best9) 0.697 vs 0.576, "
          f"RMSD(best9) 2.719 vs 2.952 A, Wilcoxon p=0.21")

    deltas = []
    for c in complete_cases:
        cons = by_case[c]["conservative_water"]["median_best"]
        simp = by_case[c]["simplified_no_water"]["median_best"]
        deltas.append((c, cons, simp, simp - cons))
    deltas.sort(key=lambda x: -abs(x[3]))
    n_help = sum(1 for _, _, _, d in deltas if d > EFFECT_THRESHOLD_A)
    n_hurt = sum(1 for _, _, _, d in deltas if d < -EFFECT_THRESHOLD_A)
    n_none = len(deltas) - n_help - n_hurt
    print(f"\n{n_help} help / {n_hurt} hurt / {n_none} no-effect (of {len(complete_cases)}, "
          f">{EFFECT_THRESHOLD_A} A threshold)")
    print("\nTop 10 |delta| cases:")
    for c, cons, simp, d in deltas[:10]:
        tag = "water helps" if d > EFFECT_THRESHOLD_A else ("water hurts" if d < -EFFECT_THRESHOLD_A else "no effect")
        print(f"  {c}  cons={cons:.3f}  simp={simp:.3f}  delta(simp-cons)={d:+.3f}  {tag}")

    with (OUT_DIR / "f7_by_case_symmetry.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "policy", "median_top1", "median_best", "unstable"])
        w.writeheader()
        w.writerows(case_rows)
    print(f"\nWrote {OUT_DIR / 'f7_by_case_symmetry.csv'}")


if __name__ == "__main__":
    main()
