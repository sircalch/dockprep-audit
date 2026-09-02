"""F9 aggregation: symmetry-corrected RMSD for the 4 oriented-water runs,
compared against those same 4 cases' existing conservative_water (bare
oxygen) and simplified_no_water (no water) results already computed in
benchmark/symmetry-rmsd/phase_c60_symmetry_rmsd.csv (F4). No new definition
introduced -- same symmetry-correction method, same success threshold.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from compute_symmetry_rmsd import process_run_dir, resolve_run_dir  # noqa: E402

CASES = ["6ASH", "1WBK", "1CVZ", "4GID"]
SEEDS = [42, 123, 2024]
OUT_DIR = REPO_ROOT / "benchmark" / "symmetry-rmsd"
F9_RUNS_ROOT = REPO_ROOT / "benchmark" / "phase-f9-oriented-water" / "conservative_oriented"
MAIN_CSV = OUT_DIR / "phase_c60_symmetry_rmsd.csv"


def load_main_cohort_bests() -> dict:
    rows = list(csv.DictReader(open(MAIN_CSV, encoding="utf-8")))
    by = {}
    for r in rows:
        if r["pdb_id"] in CASES:
            by.setdefault((r["pdb_id"], r["policy"]), []).append(float(r["sym_best"]))
    return by


def main() -> None:
    automorphism_cache: dict = {}
    new_rows = []
    for pdb_id in CASES:
        for seed in SEEDS:
            run_dir = F9_RUNS_ROOT / f"seed{seed}" / pdb_id
            result = process_run_dir(run_dir, automorphism_cache)
            if result is None:
                print(f"MISSING: {pdb_id} seed={seed}")
                continue
            new_rows.append({"pdb_id": pdb_id, "policy": "conservative_oriented", "seed": seed, **result,
                              "orig_rmsd_by_rank": ";".join(str(x) for x in result["orig_rmsd_by_rank"]),
                              "sym_rmsd_by_rank": ";".join(str(x) for x in result["sym_rmsd_by_rank"])})

    out_csv = OUT_DIR / "f9_oriented_water_symmetry_rmsd.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(new_rows[0].keys()))
        w.writeheader()
        w.writerows(new_rows)
    print(f"Wrote {out_csv} ({len(new_rows)} runs)\n")

    main_bests = load_main_cohort_bests()

    print(f"{'case':6s} {'bare-O (cons.)':>15s} {'no-water (simp.)':>17s} {'oriented (F9)':>14s}")
    for pdb_id in CASES:
        bare_o = statistics.median(main_bests.get((pdb_id, "conservative_water"), []))
        no_water = statistics.median(main_bests.get((pdb_id, "simplified_no_water"), []))
        oriented_vals = [r["sym_best"] for r in new_rows if r["pdb_id"] == pdb_id]
        oriented = statistics.median(oriented_vals) if oriented_vals else float("nan")
        print(f"{pdb_id:6s} {bare_o:15.3f} {no_water:17.3f} {oriented:14.3f}")

    print("\nDeltas relative to no-water baseline (positive = water helps):")
    for pdb_id in CASES:
        no_water = statistics.median(main_bests.get((pdb_id, "simplified_no_water"), []))
        bare_o = statistics.median(main_bests.get((pdb_id, "conservative_water"), []))
        oriented_vals = [r["sym_best"] for r in new_rows if r["pdb_id"] == pdb_id]
        oriented = statistics.median(oriented_vals) if oriented_vals else float("nan")
        print(f"  {pdb_id}: bare-O delta={no_water - bare_o:+.3f} A   "
              f"oriented delta={no_water - oriented:+.3f} A")


if __name__ == "__main__":
    main()
