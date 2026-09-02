"""Fase D driver: stratum-level aggregation of Fase C results.

Generalizes the ad-hoc analysis done by hand for the 12-case pilot (see
PROJECT-ROADMAP.md section 9, Fase D) so it can be re-run identically on the
60-case cohort. Two aggregation levels are computed per case x policy,
because the 12-case pilot found that collapsing to a binary success/fail
threshold can hide a real, reproducible water-policy effect (1M17: best pose
2.17-2.38 A with bridging water vs. 4.29-4.78 A without, across all 3 seeds
-- invisible in the top-1-pose binary outcome). Report both, not just one:

  - top1: median top-1-pose RMSD and success fraction across seeds
           (top-1 = the pose Vina itself ranks first by score)
  - best:  median best-of-all-poses RMSD across seeds
           (best = the closest of the n poses Vina returned, regardless of
           Vina's own ranking -- this is what caught the 1M17 effect)

A case is flagged "unstable" if its top1 outcome category (success /
scoring_fail / sampling_fail) differs across the 3 seeds under the SAME
policy -- this is expected under Vina's documented multi-threaded
non-determinism, not a pipeline bug, but it means a single-seed result for
that case should not be trusted.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def load_raw_runs(raw_runs_dir: Path) -> list[dict]:
    records = []
    for p in sorted(raw_runs_dir.glob("*.json")):
        try:
            records.append(json.loads(p.read_text(encoding="utf-8")))
        except Exception:
            continue
    return records


def load_strata(cases_csv: Path) -> dict[str, str]:
    with cases_csv.open(newline="", encoding="utf-8") as f:
        return {row["pdb_id"]: row["stratum"] for row in csv.DictReader(f)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate Fase C results by stratum (Fase D).")
    parser.add_argument("--phase-c-dir", type=Path, default=REPO_ROOT / "benchmark" / "phase-c-60")
    parser.add_argument("--cases", type=Path, default=REPO_ROOT / "benchmark" / "pilot_cases.csv")
    parser.add_argument("--out-dir", type=Path, default=None,
                         help="Default: same as --phase-c-dir")
    args = parser.parse_args()
    out_dir = args.out_dir or args.phase_c_dir

    strata = load_strata(args.cases)
    records = load_raw_runs(args.phase_c_dir / "raw_runs")
    print(f"Loaded {len(records)} raw run records from {args.phase_c_dir / 'raw_runs'}")

    by_case_policy: dict[tuple[str, str], list[dict]] = defaultdict(list)
    n_prep_fail = 0
    for r in records:
        if r["pdb_id"] not in strata:
            continue
        by_case_policy[(r["pdb_id"], r["policy"])].append(r)
        if r["returncode"] != 0:
            n_prep_fail += 1

    case_rows = []
    for (pdb_id, policy), recs in sorted(by_case_policy.items()):
        results = [r["result"] for r in recs if r.get("result")]
        top1_rmsds = [res["top_pose_rmsd_A"] for res in results]
        best_rmsds = [res.get("best_pose_rmsd_A", res["top_pose_rmsd_A"]) for res in results]
        outcomes = [res["outcome"] for res in results]
        n_seeds_expected = len(recs)
        unstable = len(set(outcomes)) > 1 if len(outcomes) > 1 else False
        case_rows.append({
            "pdb_id": pdb_id, "stratum": strata[pdb_id], "policy": policy,
            "n_seeds": n_seeds_expected, "n_ok": len(results),
            "n_preparation_failures": n_seeds_expected - len(results),
            "median_top1_rmsd_A": round(statistics.median(top1_rmsds), 3) if top1_rmsds else "",
            "median_best_pose_rmsd_A": round(statistics.median(best_rmsds), 3) if best_rmsds else "",
            "success_fraction_top1": round(sum(1 for o in outcomes if o == "success") / len(outcomes), 3) if outcomes else "",
            "success_fraction_best": round(sum(1 for r in best_rmsds if r <= 2.0) / len(best_rmsds), 3) if best_rmsds else "",
            "unstable_across_seeds": unstable,
            "outcome_by_seed": ";".join(outcomes),
        })

    with (out_dir / "phase_d_by_case.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(case_rows[0].keys()))
        w.writeheader()
        w.writerows(case_rows)
    print(f"Wrote {out_dir / 'phase_d_by_case.csv'} ({len(case_rows)} case x policy rows)")

    by_stratum_policy: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in case_rows:
        by_stratum_policy[(row["stratum"], row["policy"])].append(row)

    stratum_rows = []
    for (stratum, policy), rows in sorted(by_stratum_policy.items()):
        top1 = [r["median_top1_rmsd_A"] for r in rows if r["median_top1_rmsd_A"] != ""]
        best = [r["median_best_pose_rmsd_A"] for r in rows if r["median_best_pose_rmsd_A"] != ""]
        succ_top1 = [r["success_fraction_top1"] for r in rows if r["success_fraction_top1"] != ""]
        succ_best = [r["success_fraction_best"] for r in rows if r["success_fraction_best"] != ""]
        n_unstable = sum(1 for r in rows if r["unstable_across_seeds"])
        stratum_rows.append({
            "stratum": stratum, "policy": policy, "n_cases": len(rows),
            "mean_success_fraction_top1": round(statistics.mean(succ_top1), 3) if succ_top1 else "",
            "mean_success_fraction_best": round(statistics.mean(succ_best), 3) if succ_best else "",
            "mean_median_top1_rmsd_A": round(statistics.mean(top1), 3) if top1 else "",
            "mean_median_best_rmsd_A": round(statistics.mean(best), 3) if best else "",
            "n_cases_unstable_across_seeds": n_unstable,
        })

    with (out_dir / "phase_d_by_stratum.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(stratum_rows[0].keys()))
        w.writeheader()
        w.writerows(stratum_rows)
    print(f"Wrote {out_dir / 'phase_d_by_stratum.csv'} ({len(stratum_rows)} stratum x policy rows)")

    print(f"\nTotal preparation failures across all runs: {n_prep_fail}")
    print("\n--- Stratum summary (mean of per-case medians) ---")
    for row in stratum_rows:
        print(f"{row['stratum']:20s} {row['policy']:22s} n={row['n_cases']:3d}  "
              f"succ(top1)={row['mean_success_fraction_top1']}  succ(best)={row['mean_success_fraction_best']}  "
              f"RMSD(top1)={row['mean_median_top1_rmsd_A']}  RMSD(best)={row['mean_median_best_rmsd_A']}  "
              f"unstable={row['n_cases_unstable_across_seeds']}")


if __name__ == "__main__":
    main()
