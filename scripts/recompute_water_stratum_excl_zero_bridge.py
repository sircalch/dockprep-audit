"""F1 follow-up (2026-08-28): SITE_BRIDGING_WATER_PRESENT, added to the audit
engine to close the software/paper gap, revealed that 3 of the 15
water_policy-stratum cases (1HRN, 1PPM, 3FNU) had bridging_waters_kept=0
under the conservative policy -- the cohort-selection screen used a looser
4.0 A water-to-ligand-only criterion (scripts/verify_pilot_eligibility.py)
than the strict, symmetric Eq. 1 criterion (<=3.0 A to both ligand and
receptor) that actually decides which water is retained. For these 3 cases
the conservative and simplified policies are physically identical inputs,
so any RMSD difference is Vina's own multi-threaded non-determinism
(Limitation 6), not a water-policy effect.

This script recomputes the water_policy row of Table 3 and the per-case
Table 4 breakdown with those 3 cases excluded (n=12), using the same
methodology as scripts/run_phase_d.py (equal-weight mean of per-case
values; paired Wilcoxon signed-rank test on best-of-9 RMSD).
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

from scipy.stats import wilcoxon

REPO_ROOT = Path(__file__).resolve().parent.parent
BY_CASE_CSV = REPO_ROOT / "benchmark" / "phase-c-60" / "phase_d_by_case.csv"

ZERO_BRIDGE_CASES = {"1HRN", "1PPM", "3FNU"}


def load():
    rows = list(csv.DictReader(open(BY_CASE_CSV, encoding="utf-8")))
    by_case = {}
    for r in rows:
        if r["stratum"] != "water_policy":
            continue
        by_case.setdefault(r["pdb_id"], {})[r["policy"]] = r
    return by_case


def stratum_row(cases, policy, label):
    rows = [c[policy] for c in cases.values()]
    succ_top1 = [float(r["success_fraction_top1"]) for r in rows]
    succ_best = [float(r["success_fraction_best"]) for r in rows]
    top1 = [float(r["median_top1_rmsd_A"]) for r in rows]
    best = [float(r["median_best_pose_rmsd_A"]) for r in rows]
    print(f"{label:10s} {policy:22s} n={len(rows):2d}  "
          f"succ(top1)={statistics.mean(succ_top1):.3f}  RMSD(top1)={statistics.mean(top1):.3f}  "
          f"succ(best9)={statistics.mean(succ_best):.3f}  RMSD(best9)={statistics.mean(best):.3f}")


def main():
    all_cases = load()
    excl_cases = {k: v for k, v in all_cases.items() if k not in ZERO_BRIDGE_CASES}

    print(f"All water_policy cases (n={len(all_cases)}):")
    stratum_row(all_cases, "conservative_water", "n=15")
    stratum_row(all_cases, "simplified_no_water", "n=15")

    print()
    print(f"Excluding zero-bridging-water cases {sorted(ZERO_BRIDGE_CASES)} (n={len(excl_cases)}):")
    stratum_row(excl_cases, "conservative_water", "n=12")
    stratum_row(excl_cases, "simplified_no_water", "n=12")

    cons_best = [float(v["conservative_water"]["median_best_pose_rmsd_A"]) for v in excl_cases.values()]
    simp_best = [float(v["simplified_no_water"]["median_best_pose_rmsd_A"]) for v in excl_cases.values()]
    diffs = [c - s for c, s in zip(cons_best, simp_best)]
    stat, p = wilcoxon(cons_best, simp_best)
    print()
    print(f"Paired Wilcoxon (best-of-9 RMSD, conservative vs simplified, n={len(cons_best)}): "
          f"W={stat:.1f}, p={p:.4f}, mean paired diff={statistics.mean(diffs):+.3f} A")

    print()
    print("Per-case best-of-9 RMSD, n=12 (excludes 1HRN, 1PPM, 3FNU):")
    for pdb, v in sorted(excl_cases.items(), key=lambda kv: float(kv[1]["simplified_no_water"]["median_best_pose_rmsd_A"]) - float(kv[1]["conservative_water"]["median_best_pose_rmsd_A"])):
        c = float(v["conservative_water"]["median_best_pose_rmsd_A"])
        s = float(v["simplified_no_water"]["median_best_pose_rmsd_A"])
        print(f"  {pdb}  cons={c:.3f}  simp={s:.3f}  delta(simp-cons)={s - c:+.3f}")


if __name__ == "__main__":
    main()
