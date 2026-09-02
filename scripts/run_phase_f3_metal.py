"""F3 driver: real metal-policy comparison (PROJECT-ROADMAP.md section 19,
item F3), added 2026-08-28.

Runs the two metal policies (retain = project default, geometrically
assigned coordinating-HIS tautomers; remove = metal ions dropped entirely,
HIS tautomers fall back to Meeko's default; see extract_receptor_atoms() in
smoke_redock_case.py) x 3 predeclared seeds, restricted to the 15 of 60
cases confirmed by the audit engine's SITE_METAL_PRESENT finding to have a
metal ion within 6 A of the declared ligand -- the site-local criterion
(exactly the metal_or_cofactor stratum, verified 2026-08-28), following the
same F1/F2 lesson that only a site-local feature can plausibly affect the
redocking result.

water_policy is held fixed at "conservative" and altloc_policy at
"highest_occupancy" for every run in this comparison, isolating the
metal-policy manipulation the same way F2 isolated the altLoc manipulation.

Writes one JSON result per (case, metal_policy, seed) under
benchmark/phase-f3-metal/raw_runs/, plus a consolidated per-case CSV.
Caches already-completed runs, so this script can be safely re-run/resumed.
"""

from __future__ import annotations

import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

SITE_METAL_CASES = [
    "5A2S", "1CBX", "4EXS", "3LXE", "6ZR9", "1Z9Y", "3IBI", "1O86",
    "1UZE", "1DTH", "4G9L", "4JA1", "6TMN", "1KJO", "4RSY",
]
METAL_POLICIES = ["retain", "remove"]
SEEDS = [42, 123, 2024]
EXHAUSTIVENESS = 32
WATER_POLICY = "conservative"
ALTLOC_POLICY = "highest_occupancy"
OUT_DIR = REPO_ROOT / "benchmark" / "phase-f3-metal"
PYTHON = REPO_ROOT / ".venv" / "Scripts" / "python.exe"


def run_one(pdb_id: str, metal_policy: str, seed: int) -> dict:
    out_subdir = OUT_DIR / metal_policy / f"seed{seed}"
    result = subprocess.run(
        [str(PYTHON if PYTHON.exists() else sys.executable),
         str(REPO_ROOT / "scripts" / "smoke_redock_case.py"), pdb_id,
         "--out-dir", str(out_subdir), "--water-policy", WATER_POLICY,
         "--altloc-policy", ALTLOC_POLICY, "--metal-policy", metal_policy,
         "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(seed)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    record = {"pdb_id": pdb_id, "metal_policy": metal_policy, "seed": seed, "returncode": result.returncode}
    if result.returncode == 0:
        try:
            record["result"] = json.loads(result.stdout[result.stdout.index("{"):])
        except Exception:
            record["result"] = None
            record["stdout_tail"] = result.stdout[-2000:]
    else:
        record["stderr_tail"] = result.stderr[-2000:]
    (OUT_DIR / "raw_runs").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "raw_runs" / f"{pdb_id}_{metal_policy}_seed{seed}.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8"
    )
    return record


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    already_done = {}
    if (OUT_DIR / "raw_runs").exists():
        for p in (OUT_DIR / "raw_runs").glob("*.json"):
            try:
                rec = json.loads(p.read_text(encoding="utf-8"))
                already_done[(rec["pdb_id"], rec["metal_policy"], rec["seed"])] = rec
            except Exception:
                pass

    all_records = []
    total = len(SITE_METAL_CASES) * len(METAL_POLICIES) * len(SEEDS)
    done = 0
    for pdb_id in SITE_METAL_CASES:
        for metal_policy in METAL_POLICIES:
            for seed in SEEDS:
                done += 1
                key = (pdb_id, metal_policy, seed)
                if key in already_done:
                    print(f"[{done}/{total}] {pdb_id} {metal_policy} seed={seed} -> cached", flush=True)
                    all_records.append(already_done[key])
                    continue
                print(f"[{done}/{total}] {pdb_id} {metal_policy} seed={seed}", flush=True)
                rec = run_one(pdb_id, metal_policy, seed)
                all_records.append(rec)
                status = "OK" if rec["returncode"] == 0 else "FAIL"
                rmsd = (rec.get("result") or {}).get("top_pose_rmsd_A") if rec.get("result") else None
                outcome = (rec.get("result") or {}).get("outcome") if rec.get("result") else None
                print(f"    -> {status} rmsd={rmsd} outcome={outcome}", flush=True)

    summary_rows = []
    for pdb_id in SITE_METAL_CASES:
        for metal_policy in METAL_POLICIES:
            recs = [r for r in all_records if r["pdb_id"] == pdb_id and r["metal_policy"] == metal_policy]
            top1_rmsds = [r["result"]["top_pose_rmsd_A"] for r in recs if r.get("result")]
            best_rmsds = [r["result"]["best_pose_rmsd_A"] for r in recs if r.get("result")]
            outcomes = [r["result"]["outcome"] for r in recs if r.get("result")]
            n_fail_prep = sum(1 for r in recs if r["returncode"] != 0)
            summary_rows.append({
                "pdb_id": pdb_id, "metal_policy": metal_policy,
                "n_runs": len(recs), "n_preparation_failures": n_fail_prep,
                "median_top1_rmsd_A": round(statistics.median(top1_rmsds), 3) if top1_rmsds else "",
                "median_best_pose_rmsd_A": round(statistics.median(best_rmsds), 3) if best_rmsds else "",
                "success_fraction_top1": round(sum(1 for o in outcomes if o == "success") / len(outcomes), 3) if outcomes else "",
                "success_fraction_best": round(sum(1 for r in best_rmsds if r <= 2.0) / len(best_rmsds), 3) if best_rmsds else "",
                "outcome_by_seed": ";".join(outcomes),
            })

    with (OUT_DIR / "phase_f3_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nWrote {OUT_DIR / 'phase_f3_summary.csv'}")
    n_fail = sum(1 for r in all_records if r["returncode"] != 0)
    print(f"{len(all_records) - n_fail}/{len(all_records)} runs completed without a preparation error.")


if __name__ == "__main__":
    main()
