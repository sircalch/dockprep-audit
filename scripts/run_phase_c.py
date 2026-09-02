"""Fase C driver: run the 2 real preparation policies (conservative water vs
simplified/open, decided 2026-08-21 -- policies 2 and 3 are computationally
identical in this pipeline, so only 2 are run, not 3, per that decision)
across the 12 frozen pilot cases, each with 3 predeclared seeds (42, 123,
2024) at exhaustiveness 32, per PROJECT-ROADMAP.md pendiente #11/#12.

Writes one JSON result per (case, policy, seed) plus a consolidated CSV
summary (median top-1 RMSD and outcome per case x policy).
"""

from __future__ import annotations

import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CASES = ["1M17", "1T46", "4RJ3", "5A2S", "1CBX", "4EXS", "1OHR", "3FNU", "4GID", "1A28", "1QKT", "1RBP"]
POLICIES = {"conservative_water": "conservative", "simplified_no_water": "none"}
SEEDS = [42, 123, 2024]
EXHAUSTIVENESS = 32
OUT_DIR = REPO_ROOT / "benchmark" / "phase-c"


def run_one(pdb_id: str, policy_name: str, water_policy: str, seed: int) -> dict:
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "smoke_redock_case.py"), pdb_id,
         "--out-dir", str(OUT_DIR / policy_name / f"seed{seed}"), "--water-policy", water_policy,
         "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(seed)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    record = {"pdb_id": pdb_id, "policy": policy_name, "seed": seed, "returncode": result.returncode}
    if result.returncode == 0:
        try:
            record["result"] = json.loads(result.stdout[result.stdout.index("{"):])
        except Exception:
            record["result"] = None
            record["stdout_tail"] = result.stdout[-2000:]
    else:
        record["stderr_tail"] = result.stderr[-2000:]
    (OUT_DIR / "raw_runs").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "raw_runs" / f"{pdb_id}_{policy_name}_seed{seed}.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8"
    )
    return record


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_records = []
    total = len(CASES) * len(POLICIES) * len(SEEDS)
    done = 0
    for pdb_id in CASES:
        for policy_name, water_policy in POLICIES.items():
            for seed in SEEDS:
                done += 1
                print(f"[{done}/{total}] {pdb_id} {policy_name} seed={seed}", flush=True)
                rec = run_one(pdb_id, policy_name, water_policy, seed)
                all_records.append(rec)
                status = "OK" if rec["returncode"] == 0 else "FAIL"
                rmsd = (rec.get("result") or {}).get("top_pose_rmsd_A") if rec.get("result") else None
                outcome = (rec.get("result") or {}).get("outcome") if rec.get("result") else None
                print(f"    -> {status} rmsd={rmsd} outcome={outcome}", flush=True)

    summary_rows = []
    for pdb_id in CASES:
        for policy_name in POLICIES:
            recs = [r for r in all_records if r["pdb_id"] == pdb_id and r["policy"] == policy_name]
            rmsds = [r["result"]["top_pose_rmsd_A"] for r in recs if r.get("result")]
            outcomes = [r["result"]["outcome"] for r in recs if r.get("result")]
            n_fail_prep = sum(1 for r in recs if r["returncode"] != 0)
            summary_rows.append({
                "pdb_id": pdb_id, "policy": policy_name,
                "n_runs": len(recs), "n_preparation_failures": n_fail_prep,
                "median_top1_rmsd_A": round(statistics.median(rmsds), 3) if rmsds else "",
                "rmsd_by_seed": ";".join(str(r) for r in rmsds),
                "outcome_by_seed": ";".join(outcomes),
            })

    with (OUT_DIR / "phase_c_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nWrote {OUT_DIR / 'phase_c_summary.csv'}")


if __name__ == "__main__":
    main()
