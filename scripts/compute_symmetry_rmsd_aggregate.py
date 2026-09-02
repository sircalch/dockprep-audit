"""F4 aggregation: recompute the manuscript's key stratum/comparison
statistics using symmetry-corrected RMSD (scripts/compute_symmetry_rmsd.py
output) instead of the exact-coordinate-index RMSD, to check whether
symmetry correction changes any reported conclusion.

1OHR is excluded from the water-policy comparison: its ligand.sdf has only
40 heavy atoms against 44 in the deposited structure, because the chosen
protonation-variant SMILES (from fetch_ligand_smiles('1UN')) itself
represents a 40-atom molecule -- an RCSB chemical-component data
discrepancy, not a pipeline bug. The original, symmetry-agnostic RMSD
(reported throughout the manuscript) is unaffected, since it never used
ligand.sdf's bond graph, only the deposited PDB coordinates matched
directly against ligand.pdbqt.
"""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

from scipy.stats import wilcoxon

REPO_ROOT = Path(__file__).resolve().parent.parent
SYM_DIR = REPO_ROOT / "benchmark" / "symmetry-rmsd"

STRATUM_BY_CASE = {}
for row in csv.DictReader(open(REPO_ROOT / "benchmark" / "pilot_cases.csv", encoding="utf-8")):
    STRATUM_BY_CASE[row["pdb_id"]] = row["stratum"]

ZERO_BRIDGE_CASES = {"1HRN", "1PPM", "3FNU"}
EXCLUDED_BUG_ALTLOC = {"1GS4", "1SN5"}
SDF_MISMATCH_CASES = {"1OHR"}


def load(name):
    return list(csv.DictReader(open(SYM_DIR / f"{name}_symmetry_rmsd.csv", encoding="utf-8")))


def median_per_case(rows, policy_key, value_key, exclude=frozenset()):
    by_case_policy = defaultdict(list)
    for r in rows:
        if r["pdb_id"] in exclude:
            continue
        by_case_policy[(r["pdb_id"], r[policy_key])].append(float(r[value_key]))
    return {k: statistics.median(v) for k, v in by_case_policy.items()}


def compare_metric(label, orig_medians, sym_medians, cases, policy_a, policy_b, threshold=2.0):
    a_orig = [orig_medians[(c, policy_a)] for c in cases]
    b_orig = [orig_medians[(c, policy_b)] for c in cases]
    a_sym = [sym_medians[(c, policy_a)] for c in cases]
    b_sym = [sym_medians[(c, policy_b)] for c in cases]

    def stats(a, b):
        succ_a = sum(1 for x in a if x <= threshold) / len(a)
        succ_b = sum(1 for x in b if x <= threshold) / len(b)
        mean_a, mean_b = statistics.mean(a), statistics.mean(b)
        stat, p = wilcoxon(a, b)
        return succ_a, succ_b, mean_a, mean_b, stat, p

    so_a, so_b, mo_a, mo_b, wo, po = stats(a_orig, b_orig)
    ss_a, ss_b, ms_a, ms_b, ws, ps = stats(a_sym, b_sym)

    print(f"\n--- {label} (n={len(cases)}) ---")
    print(f"  {'':14s} {'succ_A':>8s} {'succ_B':>8s} {'RMSD_A':>8s} {'RMSD_B':>8s} {'W':>7s} {'p':>7s}")
    print(f"  {'original':14s} {so_a:8.3f} {so_b:8.3f} {mo_a:8.3f} {mo_b:8.3f} {wo:7.2f} {po:7.4f}")
    print(f"  {'symmetry':14s} {ss_a:8.3f} {ss_b:8.3f} {ms_a:8.3f} {ms_b:8.3f} {ws:7.2f} {ps:7.4f}")


def main():
    water_rows = load("phase_c60")
    altloc_rows = load("f2_altloc")
    metal_rows = load("f3_metal")

    # --- Water: best-of-9, per stratum ---
    orig_best = median_per_case(water_rows, "policy", "orig_best")
    sym_best = median_per_case(water_rows, "policy", "sym_best")
    for stratum in ["water_policy", "alternate_location", "metal_or_cofactor", "low_risk_control"]:
        cases = [c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == stratum]
        exclude = ZERO_BRIDGE_CASES | SDF_MISMATCH_CASES if stratum == "water_policy" else SDF_MISMATCH_CASES
        cases = [c for c in cases if c not in exclude]
        compare_metric(f"Water policy, best-of-9, {stratum}", orig_best, sym_best, cases,
                        "conservative_water", "simplified_no_water")

    # --- AltLoc: best-of-9, n=15 ---
    a_orig_best = median_per_case(altloc_rows, "policy", "orig_best", exclude=EXCLUDED_BUG_ALTLOC)
    a_sym_best = median_per_case(altloc_rows, "policy", "sym_best", exclude=EXCLUDED_BUG_ALTLOC)
    altloc_cases = [c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == "alternate_location"]
    compare_metric("AltLoc policy, best-of-9", a_orig_best, a_sym_best, altloc_cases,
                    "highest_occupancy", "lowest_occupancy")

    # --- Metal: best-of-9, n=15 ---
    m_orig_best = median_per_case(metal_rows, "policy", "orig_best")
    m_sym_best = median_per_case(metal_rows, "policy", "sym_best")
    metal_cases = [c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == "metal_or_cofactor"]
    compare_metric("Metal policy, best-of-9", m_orig_best, m_sym_best, metal_cases, "retain", "remove")


if __name__ == "__main__":
    main()
