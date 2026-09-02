"""F2 aggregation: real altLoc-policy comparison (PROJECT-ROADMAP.md section
19, item F2), added 2026-08-28.

Mirrors run_phase_d.py's methodology exactly (median-of-3-seeds per case,
mean-of-per-case-medians for the group aggregate, "unstable" = top-1 outcome
category differs across seeds) but for the two altloc policies
(highest_occupancy vs lowest_occupancy; see run_phase_f2_altloc.py) on the
15 SITE_ALTLOC_PRESENT cases (the full alternate_location stratum),
instead of the two water policies on all 60.

Also reports the paired Wilcoxon signed-rank test (best-of-9 RMSD, n=15)
and the per-case Delta breakdown, in the same format as Table 3/4 of the
manuscript and scripts/recompute_water_stratum_excl_zero_bridge.py, so the
result can be dropped into Results/Discussion directly once the raw runs
finish.

1GS4 and 1SN5 are reported separately, if cached runs for them exist, as a
side note rather than part of the primary group: an earlier dockprep_audit
bug (fixed 2026-08-28, v0.2.0) had counted a disordered water (1GS4) and
the ligand's own alternate conformer (1SN5) as SITE_ALTLOC_PRESENT hits,
which is why they were in an earlier version of the run list. Both are
correctly excluded post-fix (no protein conformer sits within 6 A of the
ligand for either case), but their two policy runs are not both expected
to be identical: 1GS4's receptor is byte-identical between policies (the
water altloc never enters the prepared receptor); 1SN5's is not, because a
second T3 ligand copy elsewhere in the file is retained as part of the
receptor and its altloc is resolved differently -- see
run_phase_f2_altloc.py's SITE_ALTLOC_CASES comment for the full account.
This section reports both cases' actual values rather than asserting an
expected outcome for either.

Safe to run before benchmark/phase-f2-altloc/raw_runs/ is fully populated:
reports on whatever cases have all 6 runs (2 policies x 3 seeds) completed
so far, and lists which cases are still incomplete.
"""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

from scipy.stats import wilcoxon

REPO_ROOT = Path(__file__).resolve().parent.parent
PHASE_F2_DIR = REPO_ROOT / "benchmark" / "phase-f2-altloc"
SITE_ALTLOC_CASES = [
    "1E6U", "1E7S", "1M17", "1T46", "2I4H", "3D14", "3KXG",
    "3MWU", "3P0M", "3PNA", "4RJ3", "5E0J", "5F4N", "5K8S", "5TG1",
]
EXCLUDED_BUG_CASES = ["1GS4", "1SN5"]
ALTLOC_POLICIES = ["highest_occupancy", "lowest_occupancy"]
SEEDS = [42, 123, 2024]
EFFECT_THRESHOLD_A = 0.3


def load_raw_runs(raw_runs_dir: Path) -> list[dict]:
    records = []
    for p in sorted(raw_runs_dir.glob("*.json")):
        try:
            records.append(json.loads(p.read_text(encoding="utf-8")))
        except Exception:
            continue
    return records


def main() -> None:
    raw_dir = PHASE_F2_DIR / "raw_runs"
    if not raw_dir.exists():
        print(f"{raw_dir} does not exist yet -- run scripts/run_phase_f2_altloc.py first.")
        return

    records = load_raw_runs(raw_dir)
    print(f"Loaded {len(records)} raw run records from {raw_dir} "
          f"(expect {len(SITE_ALTLOC_CASES) * len(ALTLOC_POLICIES) * len(SEEDS)} when complete)")

    all_known_cases = SITE_ALTLOC_CASES + EXCLUDED_BUG_CASES
    by_case_policy: dict[tuple[str, str], list[dict]] = defaultdict(list)
    n_prep_fail = 0
    for r in records:
        if r["pdb_id"] not in all_known_cases:
            continue
        by_case_policy[(r["pdb_id"], r["altloc_policy"])].append(r)
        if r["returncode"] != 0:
            n_prep_fail += 1

    case_rows = []
    for (pdb_id, policy), recs in sorted(by_case_policy.items()):
        results = [r["result"] for r in recs if r.get("result")]
        top1_rmsds = [res["top_pose_rmsd_A"] for res in results]
        best_rmsds = [res.get("best_pose_rmsd_A", res["top_pose_rmsd_A"]) for res in results]
        outcomes = [res["outcome"] for res in results]
        unstable = len(set(outcomes)) > 1 if len(outcomes) > 1 else False
        case_rows.append({
            "pdb_id": pdb_id, "altloc_policy": policy,
            "n_seeds": len(recs), "n_ok": len(results),
            "n_preparation_failures": len(recs) - len(results),
            "median_top1_rmsd_A": round(statistics.median(top1_rmsds), 3) if top1_rmsds else "",
            "median_best_pose_rmsd_A": round(statistics.median(best_rmsds), 3) if best_rmsds else "",
            "success_fraction_top1": round(sum(1 for o in outcomes if o == "success") / len(outcomes), 3) if outcomes else "",
            "success_fraction_best": round(sum(1 for r in best_rmsds if r <= 2.0) / len(best_rmsds), 3) if best_rmsds else "",
            "unstable_across_seeds": unstable,
            "outcome_by_seed": ";".join(outcomes),
        })

    PHASE_F2_DIR.mkdir(parents=True, exist_ok=True)
    with (PHASE_F2_DIR / "phase_f2_by_case.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(case_rows[0].keys()) if case_rows else [])
        if case_rows:
            w.writeheader()
            w.writerows(case_rows)
    print(f"Wrote {PHASE_F2_DIR / 'phase_f2_by_case.csv'} ({len(case_rows)} case x policy rows)")

    by_case = defaultdict(dict)
    for row in case_rows:
        by_case[row["pdb_id"]][row["altloc_policy"]] = row

    all_complete = {
        pdb: d for pdb, d in by_case.items()
        if all(p in d and d[p]["n_ok"] == 3 for p in ALTLOC_POLICIES)
    }
    complete_cases = {pdb: d for pdb, d in all_complete.items() if pdb in SITE_ALTLOC_CASES}
    bug_check_cases = {pdb: d for pdb, d in all_complete.items() if pdb in EXCLUDED_BUG_CASES}
    incomplete_cases = sorted(set(SITE_ALTLOC_CASES) - set(complete_cases))

    print(f"\n{len(complete_cases)}/{len(SITE_ALTLOC_CASES)} cases have all 6 runs (2 policies x 3 seeds) complete.")
    if incomplete_cases:
        print(f"Still incomplete: {incomplete_cases}")

    if bug_check_cases:
        print(f"\n--- Side note (excluded from primary n=15 group; see EXCLUDED_BUG_CASES docstring): "
              f"{sorted(bug_check_cases)} ---")
        for pdb, d in sorted(bug_check_cases.items()):
            h = d["highest_occupancy"]["median_best_pose_rmsd_A"]
            l = d["lowest_occupancy"]["median_best_pose_rmsd_A"]
            note = "identical (expected: no altloc enters the prepared receptor)" if h == l \
                else f"differs by {l - h:+.3f} A (expected for 1SN5: a second ligand copy's altloc IS resolved differently; not a primary-group finding)"
            print(f"  {pdb}  high_occ={h:.3f}  low_occ={l:.3f}  {note}")

    if not complete_cases:
        print("\nNo complete cases yet -- nothing to aggregate.")
        return

    print("\n--- Group summary (mean of per-case medians, n={}) ---".format(len(complete_cases)))
    group_rows = []
    for policy in ALTLOC_POLICIES:
        rows = [d[policy] for d in complete_cases.values()]
        top1 = [r["median_top1_rmsd_A"] for r in rows]
        best = [r["median_best_pose_rmsd_A"] for r in rows]
        succ_top1 = [r["success_fraction_top1"] for r in rows]
        succ_best = [r["success_fraction_best"] for r in rows]
        n_unstable = sum(1 for r in rows if r["unstable_across_seeds"])
        group_rows.append({
            "altloc_policy": policy, "n_cases": len(rows),
            "mean_success_fraction_top1": round(statistics.mean(succ_top1), 3),
            "mean_success_fraction_best": round(statistics.mean(succ_best), 3),
            "mean_median_top1_rmsd_A": round(statistics.mean(top1), 3),
            "mean_median_best_rmsd_A": round(statistics.mean(best), 3),
            "n_cases_unstable_across_seeds": n_unstable,
        })
        print(f"{policy:20s} n={len(rows):2d}  succ(top1)={group_rows[-1]['mean_success_fraction_top1']}  "
              f"succ(best9)={group_rows[-1]['mean_success_fraction_best']}  "
              f"RMSD(top1)={group_rows[-1]['mean_median_top1_rmsd_A']}  "
              f"RMSD(best9)={group_rows[-1]['mean_median_best_rmsd_A']}  unstable={n_unstable}/{len(rows)}")

    with (PHASE_F2_DIR / "phase_f2_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(group_rows[0].keys()))
        w.writeheader()
        w.writerows(group_rows)
    print(f"\nWrote {PHASE_F2_DIR / 'phase_f2_summary.csv'}")

    high = [d["highest_occupancy"]["median_best_pose_rmsd_A"] for d in complete_cases.values()]
    low = [d["lowest_occupancy"]["median_best_pose_rmsd_A"] for d in complete_cases.values()]
    if len(complete_cases) >= 2:
        stat, p = wilcoxon(high, low)
        mean_diff = statistics.mean(h - l for h, l in zip(high, low))
        print(f"\nPaired Wilcoxon (best-of-9 RMSD, highest_occupancy vs lowest_occupancy, n={len(complete_cases)}): "
              f"W={stat:.1f}, p={p:.4f}, mean paired diff (high-low)={mean_diff:+.3f} A")
    else:
        print("\nNeed at least 2 complete cases for the Wilcoxon test.")

    print(f"\nPer-case best-of-9 RMSD, sorted by |delta| descending (threshold {EFFECT_THRESHOLD_A} A):")
    deltas = []
    for pdb, d in complete_cases.items():
        h = d["highest_occupancy"]["median_best_pose_rmsd_A"]
        l = d["lowest_occupancy"]["median_best_pose_rmsd_A"]
        deltas.append((pdb, h, l, l - h))
    deltas.sort(key=lambda x: -abs(x[3]))
    n_help, n_hurt, n_none = 0, 0, 0
    for pdb, h, l, delta in deltas:
        if delta > EFFECT_THRESHOLD_A:
            tag = "high-occ helps"
            n_help += 1
        elif delta < -EFFECT_THRESHOLD_A:
            tag = "high-occ hurts"
            n_hurt += 1
        else:
            tag = "no effect"
            n_none += 1
        print(f"  {pdb}  high_occ={h:.3f}  low_occ={l:.3f}  delta(low-high)={delta:+.3f}  {tag}")
    print(f"\n{n_help} help / {n_hurt} hurt / {n_none} no-effect (of {len(complete_cases)} complete cases, "
          f">{EFFECT_THRESHOLD_A} A threshold)")

    print(f"\nTotal preparation failures across all runs: {n_prep_fail}")


if __name__ == "__main__":
    main()
