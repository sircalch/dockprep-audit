"""F3 aggregation: real metal-policy comparison (PROJECT-ROADMAP.md section
19, item F3), added 2026-08-28.

Mirrors run_phase_f2_aggregate.py's methodology exactly (median-of-3-seeds
per case, mean-of-per-case-medians for the group aggregate, "unstable" =
top-1 outcome category differs across seeds) but for the two metal policies
(retain vs remove; see run_phase_f3_metal.py) on the 15 SITE_METAL_PRESENT
cases (the full metal_or_cofactor stratum).

Also reports the paired Wilcoxon signed-rank test (best-of-9 RMSD, n=15)
and the per-case Delta breakdown, in the same format as Table 3/4/5/6 of
the manuscript, so the result can be dropped into Results/Discussion
directly once the raw runs finish.

Safe to run before benchmark/phase-f3-metal/raw_runs/ is fully populated:
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
PHASE_F3_DIR = REPO_ROOT / "benchmark" / "phase-f3-metal"
SITE_METAL_CASES = [
    "5A2S", "1CBX", "4EXS", "3LXE", "6ZR9", "1Z9Y", "3IBI", "1O86",
    "1UZE", "1DTH", "4G9L", "4JA1", "6TMN", "1KJO", "4RSY",
]
METAL_POLICIES = ["retain", "remove"]
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
    raw_dir = PHASE_F3_DIR / "raw_runs"
    if not raw_dir.exists():
        print(f"{raw_dir} does not exist yet -- run scripts/run_phase_f3_metal.py first.")
        return

    records = load_raw_runs(raw_dir)
    print(f"Loaded {len(records)} raw run records from {raw_dir} "
          f"(expect {len(SITE_METAL_CASES) * len(METAL_POLICIES) * len(SEEDS)} when complete)")

    by_case_policy: dict[tuple[str, str], list[dict]] = defaultdict(list)
    n_prep_fail = 0
    for r in records:
        if r["pdb_id"] not in SITE_METAL_CASES:
            continue
        by_case_policy[(r["pdb_id"], r["metal_policy"])].append(r)
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
            "pdb_id": pdb_id, "metal_policy": policy,
            "n_seeds": len(recs), "n_ok": len(results),
            "n_preparation_failures": len(recs) - len(results),
            "median_top1_rmsd_A": round(statistics.median(top1_rmsds), 3) if top1_rmsds else "",
            "median_best_pose_rmsd_A": round(statistics.median(best_rmsds), 3) if best_rmsds else "",
            "success_fraction_top1": round(sum(1 for o in outcomes if o == "success") / len(outcomes), 3) if outcomes else "",
            "success_fraction_best": round(sum(1 for r in best_rmsds if r <= 2.0) / len(best_rmsds), 3) if best_rmsds else "",
            "unstable_across_seeds": unstable,
            "outcome_by_seed": ";".join(outcomes),
        })

    PHASE_F3_DIR.mkdir(parents=True, exist_ok=True)
    with (PHASE_F3_DIR / "phase_f3_by_case.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(case_rows[0].keys()) if case_rows else [])
        if case_rows:
            w.writeheader()
            w.writerows(case_rows)
    print(f"Wrote {PHASE_F3_DIR / 'phase_f3_by_case.csv'} ({len(case_rows)} case x policy rows)")

    by_case = defaultdict(dict)
    for row in case_rows:
        by_case[row["pdb_id"]][row["metal_policy"]] = row

    complete_cases = {
        pdb: d for pdb, d in by_case.items()
        if all(p in d and d[p]["n_ok"] == 3 for p in METAL_POLICIES)
    }
    incomplete_cases = sorted(set(SITE_METAL_CASES) - set(complete_cases))

    print(f"\n{len(complete_cases)}/{len(SITE_METAL_CASES)} cases have all 6 runs (2 policies x 3 seeds) complete.")
    if incomplete_cases:
        print(f"Still incomplete: {incomplete_cases}")

    if not complete_cases:
        print("\nNo complete cases yet -- nothing to aggregate.")
        return

    print("\n--- Group summary (mean of per-case medians, n={}) ---".format(len(complete_cases)))
    group_rows = []
    for policy in METAL_POLICIES:
        rows = [d[policy] for d in complete_cases.values()]
        top1 = [r["median_top1_rmsd_A"] for r in rows]
        best = [r["median_best_pose_rmsd_A"] for r in rows]
        succ_top1 = [r["success_fraction_top1"] for r in rows]
        succ_best = [r["success_fraction_best"] for r in rows]
        n_unstable = sum(1 for r in rows if r["unstable_across_seeds"])
        group_rows.append({
            "metal_policy": policy, "n_cases": len(rows),
            "mean_success_fraction_top1": round(statistics.mean(succ_top1), 3),
            "mean_success_fraction_best": round(statistics.mean(succ_best), 3),
            "mean_median_top1_rmsd_A": round(statistics.mean(top1), 3),
            "mean_median_best_rmsd_A": round(statistics.mean(best), 3),
            "n_cases_unstable_across_seeds": n_unstable,
        })
        print(f"{policy:10s} n={len(rows):2d}  succ(top1)={group_rows[-1]['mean_success_fraction_top1']}  "
              f"succ(best9)={group_rows[-1]['mean_success_fraction_best']}  "
              f"RMSD(top1)={group_rows[-1]['mean_median_top1_rmsd_A']}  "
              f"RMSD(best9)={group_rows[-1]['mean_median_best_rmsd_A']}  unstable={n_unstable}/{len(rows)}")

    with (PHASE_F3_DIR / "phase_f3_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(group_rows[0].keys()))
        w.writeheader()
        w.writerows(group_rows)
    print(f"\nWrote {PHASE_F3_DIR / 'phase_f3_summary.csv'}")

    retain = [d["retain"]["median_best_pose_rmsd_A"] for d in complete_cases.values()]
    remove = [d["remove"]["median_best_pose_rmsd_A"] for d in complete_cases.values()]
    if len(complete_cases) >= 2:
        stat, p = wilcoxon(retain, remove)
        mean_diff = statistics.mean(r - x for r, x in zip(retain, remove))
        print(f"\nPaired Wilcoxon (best-of-9 RMSD, retain vs remove, n={len(complete_cases)}): "
              f"W={stat:.1f}, p={p:.4f}, mean paired diff (retain-remove)={mean_diff:+.3f} A")
    else:
        print("\nNeed at least 2 complete cases for the Wilcoxon test.")

    print(f"\nPer-case best-of-9 RMSD, sorted by |delta| descending (threshold {EFFECT_THRESHOLD_A} A):")
    deltas = []
    for pdb, d in complete_cases.items():
        r = d["retain"]["median_best_pose_rmsd_A"]
        x = d["remove"]["median_best_pose_rmsd_A"]
        deltas.append((pdb, r, x, x - r))
    deltas.sort(key=lambda t: -abs(t[3]))
    n_help, n_hurt, n_none = 0, 0, 0
    for pdb, r, x, delta in deltas:
        if delta > EFFECT_THRESHOLD_A:
            tag = "retain helps"
            n_help += 1
        elif delta < -EFFECT_THRESHOLD_A:
            tag = "retain hurts"
            n_hurt += 1
        else:
            tag = "no effect"
            n_none += 1
        print(f"  {pdb}  retain={r:.3f}  remove={x:.3f}  delta(remove-retain)={delta:+.3f}  {tag}")
    print(f"\n{n_help} help / {n_hurt} hurt / {n_none} no-effect (of {len(complete_cases)} complete cases, "
          f">{EFFECT_THRESHOLD_A} A threshold)")

    print(f"\nTotal preparation failures across all runs: {n_prep_fail}")


if __name__ == "__main__":
    main()
