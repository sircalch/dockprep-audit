"""F4 promotion: build every symmetry-corrected number the manuscript needs
(Tables 3-8 equivalents, Wilcoxon tests, per-case delta breakdowns), mirroring
the exact methodology already used for the original (exact-coordinate-index)
numbers in run_phase_d.py / run_phase_f2_aggregate.py / run_phase_f3_aggregate.py,
but reading sym_top1/sym_best from benchmark/symmetry-rmsd/*.csv instead of
recomputing RMSD.

1OHR is excluded throughout (see compute_symmetry_rmsd_aggregate.py docstring
for why: ligand.sdf built from a 40-heavy-atom SMILES against a 44-atom
deposited ligand, an RCSB data discrepancy, not a pipeline bug -- it has no
usable automorphism data). It keeps its original (non-symmetry-corrected)
value in the low_risk... no, water_policy stratum reporting, flagged
explicitly, exactly like the 1HRN/1PPM/3FNU zero-bridge exclusion.
"""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

from scipy.stats import wilcoxon

REPO_ROOT = Path(__file__).resolve().parent.parent
SYM_DIR = REPO_ROOT / "benchmark" / "symmetry-rmsd"
OUT_DIR = SYM_DIR

STRATUM_BY_CASE = {}
for row in csv.DictReader(open(REPO_ROOT / "benchmark" / "pilot_cases.csv", encoding="utf-8")):
    STRATUM_BY_CASE[row["pdb_id"]] = row["stratum"]

ZERO_BRIDGE_CASES = {"1HRN", "1PPM", "3FNU"}
EXCLUDED_BUG_ALTLOC = {"1GS4", "1SN5"}
SDF_MISMATCH_CASES = {"1OHR"}
SUCCESS_THRESHOLD = 2.0
EFFECT_THRESHOLD = 0.3


def load(name):
    return list(csv.DictReader(open(SYM_DIR / f"{name}_symmetry_rmsd.csv", encoding="utf-8")))


def per_case_medians(rows, exclude=frozenset()):
    """Returns {(pdb_id, policy): {"top1": median, "best": median, "outcome_by_seed": [...]}}"""
    by = defaultdict(lambda: {"top1": [], "best": [], "outcomes": []})
    for r in rows:
        if r["pdb_id"] in exclude:
            continue
        key = (r["pdb_id"], r["policy"])
        by[key]["top1"].append(float(r["sym_top1"]))
        by[key]["best"].append(float(r["sym_best"]))
        by[key]["outcomes"].append(r["sym_outcome"])
    out = {}
    for key, d in by.items():
        out[key] = {
            "median_top1": statistics.median(d["top1"]),
            "median_best": statistics.median(d["best"]),
            "success_top1": sum(1 for x in d["top1"] if x <= SUCCESS_THRESHOLD) / len(d["top1"]),
            "success_best": sum(1 for x in d["best"] if x <= SUCCESS_THRESHOLD) / len(d["best"]),
            "unstable": len(set(d["outcomes"])) > 1,
        }
    return out


def stratum_row(case_medians, cases, policy):
    vals = [case_medians[(c, policy)] for c in cases]
    return {
        "n": len(vals),
        "succ_top1": statistics.mean(v["success_top1"] for v in vals),
        "rmsd_top1": statistics.mean(v["median_top1"] for v in vals),
        "succ_best": statistics.mean(v["success_best"] for v in vals),
        "rmsd_best": statistics.mean(v["median_best"] for v in vals),
        "unstable": sum(1 for v in vals if v["unstable"]),
    }


def wilcoxon_test(case_medians, cases, policy_a, policy_b):
    a = [case_medians[(c, policy_a)]["median_best"] for c in cases]
    b = [case_medians[(c, policy_b)]["median_best"] for c in cases]
    stat, p = wilcoxon(a, b)
    mean_diff = statistics.mean(x - y for x, y in zip(a, b))
    return stat, p, mean_diff


def delta_table(case_medians, cases, policy_a, policy_b, label_a, label_b, sign_a_minus_b=True):
    rows = []
    for c in cases:
        va = case_medians[(c, policy_a)]["median_best"]
        vb = case_medians[(c, policy_b)]["median_best"]
        delta = vb - va if sign_a_minus_b else va - vb
        rows.append((c, va, vb, delta))
    rows.sort(key=lambda r: -abs(r[3]))
    return rows


def print_stratum_table(title, case_medians, strata_policy_pairs):
    print(f"\n=== {title} ===")
    for stratum, policy, cases in strata_policy_pairs:
        r = stratum_row(case_medians, cases, policy)
        print(f"  {stratum:22s} {policy:22s} n={r['n']:2d}  succ(top1)={r['succ_top1']:.3f}  "
              f"RMSD(top1)={r['rmsd_top1']:.3f}  succ(best9)={r['succ_best']:.3f}  "
              f"RMSD(best9)={r['rmsd_best']:.3f}  unstable={r['unstable']}/{r['n']}")


def main():
    water_rows = load("phase_c60")
    altloc_rows = load("f2_altloc")
    metal_rows = load("f3_metal")

    # ---------------- Water: Table 3 equivalent (4 strata) ----------------
    water_exclude_general = SDF_MISMATCH_CASES
    water_medians = per_case_medians(water_rows, exclude=water_exclude_general)
    strata_pairs = []
    for stratum in ["alternate_location", "metal_or_cofactor", "water_policy", "low_risk_control"]:
        cases = sorted(c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == stratum and c not in water_exclude_general)
        if stratum == "water_policy":
            cases_excl_zero = [c for c in cases if c not in ZERO_BRIDGE_CASES]
        else:
            cases_excl_zero = cases
        for policy in ["conservative_water", "simplified_no_water"]:
            strata_pairs.append((stratum, policy, cases_excl_zero if stratum == "water_policy" else cases))
    print_stratum_table("Table 3 (symmetry-corrected)", water_medians, strata_pairs)

    water_policy_cases = sorted(c for c in STRATUM_BY_CASE
                                 if STRATUM_BY_CASE[c] == "water_policy" and c not in ZERO_BRIDGE_CASES | water_exclude_general)
    print(f"\nwater_policy cases used (n={len(water_policy_cases)}): {water_policy_cases}")
    w_stat, w_p, w_diff = wilcoxon_test(water_medians, water_policy_cases, "conservative_water", "simplified_no_water")
    print(f"Wilcoxon water (conservative vs simplified), n={len(water_policy_cases)}: W={w_stat:.1f} p={w_p:.4f} mean_diff(cons-simp)={w_diff:+.3f}")

    print("\nTable 4 (symmetry-corrected) candidates, water_policy, |delta|>0.3:")
    for c, va, vb, delta in delta_table(water_medians, water_policy_cases, "conservative_water", "simplified_no_water", "cons", "simp"):
        tag = "water helps" if delta > EFFECT_THRESHOLD else ("water hurts" if delta < -EFFECT_THRESHOLD else "no effect")
        print(f"  {c}  cons={va:.3f}  simp={vb:.3f}  delta(simp-cons)={delta:+.3f}  {tag}")

    # also report other 3 strata Wilcoxon for completeness
    for stratum in ["alternate_location", "metal_or_cofactor", "low_risk_control"]:
        cases = sorted(c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == stratum and c not in water_exclude_general)
        st, p, d = wilcoxon_test(water_medians, cases, "conservative_water", "simplified_no_water")
        print(f"Wilcoxon water manipulation, {stratum}: W={st:.1f} p={p:.4f}")

    # ---------------- AltLoc: Table 5/6 equivalent ----------------
    altloc_medians = per_case_medians(altloc_rows, exclude=EXCLUDED_BUG_ALTLOC)
    altloc_cases = sorted(c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == "alternate_location")
    print_stratum_table("Table 5 (symmetry-corrected)", altloc_medians,
                        [("alternate_location", "highest_occupancy", altloc_cases),
                         ("alternate_location", "lowest_occupancy", altloc_cases)])
    a_stat, a_p, a_diff = wilcoxon_test(altloc_medians, altloc_cases, "highest_occupancy", "lowest_occupancy")
    print(f"Wilcoxon altLoc (high vs low), n={len(altloc_cases)}: W={a_stat:.1f} p={a_p:.4f} mean_diff(high-low)={a_diff:+.3f}")
    print("\nTable 6 (symmetry-corrected) candidates, |delta|>0.3:")
    for c, va, vb, delta in delta_table(altloc_medians, altloc_cases, "highest_occupancy", "lowest_occupancy", "high", "low"):
        tag = "high-occ helps" if delta > EFFECT_THRESHOLD else ("high-occ hurts" if delta < -EFFECT_THRESHOLD else "no effect")
        print(f"  {c}  high={va:.3f}  low={vb:.3f}  delta(low-high)={delta:+.3f}  {tag}")

    # ---------------- Metal: Table 7/8 equivalent ----------------
    metal_medians = per_case_medians(metal_rows)
    metal_cases = sorted(c for c in STRATUM_BY_CASE if STRATUM_BY_CASE[c] == "metal_or_cofactor")
    print_stratum_table("Table 7 (symmetry-corrected)", metal_medians,
                        [("metal_or_cofactor", "retain", metal_cases),
                         ("metal_or_cofactor", "remove", metal_cases)])
    m_stat, m_p, m_diff = wilcoxon_test(metal_medians, metal_cases, "retain", "remove")
    print(f"Wilcoxon metal (retain vs remove), n={len(metal_cases)}: W={m_stat:.1f} p={m_p:.4f} mean_diff(retain-remove)={m_diff:+.3f}")
    print("\nTable 8 (symmetry-corrected) candidates, |delta|>0.3:")
    for c, va, vb, delta in delta_table(metal_medians, metal_cases, "retain", "remove", "retain", "remove"):
        tag = "retain helps" if delta < -EFFECT_THRESHOLD else ("retain hurts" if delta > EFFECT_THRESHOLD else "no effect")
        print(f"  {c}  retain={va:.3f}  remove={vb:.3f}  delta(remove-retain)={delta:+.3f}  {tag}")

    # ---------------- 1EPP/1PPM specific values (Section 3.4 prose) ----------------
    print("\n1EPP/1PPM symmetry-corrected best-of-9 (water_policy prose in 3.4):")
    for c in ["1EPP", "1PPM"]:
        cw = water_medians[(c, "conservative_water")]
        sw = water_medians[(c, "simplified_no_water")]
        print(f"  {c}: cons_top1={cw['median_top1']:.3f} cons_best={cw['median_best']:.3f} "
              f"simp_top1={sw['median_top1']:.3f} simp_best={sw['median_best']:.3f}")

    # 6ASH (3.4 example) and 1A28 (3.4 example)
    print("\n6ASH / 1A28 symmetry-corrected (Section 3.4 prose examples):")
    for c in ["6ASH", "1A28"]:
        cw = water_medians[(c, "conservative_water")]
        sw = water_medians[(c, "simplified_no_water")]
        print(f"  {c}: cons_best={cw['median_best']:.3f} simp_best={sw['median_best']:.3f} "
              f"cons_top1={cw['median_top1']:.3f} simp_top1={sw['median_top1']:.3f}")


if __name__ == "__main__":
    main()
