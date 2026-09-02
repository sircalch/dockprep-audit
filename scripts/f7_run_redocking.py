"""F7 step 3: redock the frozen external validation cohort (manifest.csv)
with the same two water policies, exhaustiveness, and seeds as the main
60-case cohort's Fase C -- no protocol adjustment after seeing screening
results, per the discovery/validation design (module docstring of
f7_freeze_cohort.py).

Mirrors run_phase_f2_altloc.py's driver structure exactly. Writes one JSON
result per (case, policy, seed) under benchmark/f7-external-validation/raw_runs/,
cached so this script is safely resumable.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "benchmark" / "f7-external-validation"
RAW_DIR = OUT_DIR / "raw-pdb"
PYTHON = REPO_ROOT / ".venv" / "Scripts" / "python.exe"

POLICIES = {"conservative_water": "conservative", "simplified_no_water": "none"}
SEEDS = [42, 123, 2024]
EXHAUSTIVENESS = 32


def run_one(pdb_id: str, policy_name: str, water_policy: str, seed: int, case_row: dict) -> dict:
    out_subdir = OUT_DIR / "runs" / policy_name / f"seed{seed}"
    result = subprocess.run(
        [str(PYTHON if PYTHON.exists() else sys.executable),
         str(REPO_ROOT / "scripts" / "smoke_redock_case.py"), pdb_id,
         "--cases", str(OUT_DIR / "manifest_as_cases.csv"),
         "--eligibility", str(OUT_DIR / "manifest_as_eligibility.csv"),
         "--raw-dir", str(RAW_DIR),
         "--out-dir", str(out_subdir), "--water-policy", water_policy,
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


def prepare_case_and_eligibility_csvs(manifest_rows: list[dict]) -> None:
    """smoke_redock_case.py expects separate --cases (ligand identity) and
    --eligibility (receptor_chains) CSVs; the F7 manifest already has both,
    just split into the two column sets it expects."""
    with (OUT_DIR / "manifest_as_cases.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "ligand_component_id", "ligand_chain", "ligand_resseq"])
        w.writeheader()
        for r in manifest_rows:
            w.writerow({"pdb_id": r["pdb_id"], "ligand_component_id": r["ligand_component_id"],
                        "ligand_chain": r["ligand_chain"], "ligand_resseq": r["ligand_resseq"]})
    with (OUT_DIR / "manifest_as_eligibility.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "receptor_chains"])
        w.writeheader()
        for r in manifest_rows:
            w.writerow({"pdb_id": r["pdb_id"], "receptor_chains": r["receptor_chains"]})


def main() -> None:
    manifest_rows = list(csv.DictReader(open(OUT_DIR / "manifest.csv", encoding="utf-8")))
    prepare_case_and_eligibility_csvs(manifest_rows)
    cases = [r["pdb_id"] for r in manifest_rows]

    already_done = {}
    raw_runs_dir = OUT_DIR / "raw_runs"
    if raw_runs_dir.exists():
        for p in raw_runs_dir.glob("*.json"):
            try:
                rec = json.loads(p.read_text(encoding="utf-8"))
                already_done[(rec["pdb_id"], rec["policy"], rec["seed"])] = rec
            except Exception:
                pass

    all_records = []
    total = len(cases) * len(POLICIES) * len(SEEDS)
    done = 0
    by_pdb = {r["pdb_id"]: r for r in manifest_rows}
    for pdb_id in cases:
        for policy_name, water_policy in POLICIES.items():
            for seed in SEEDS:
                done += 1
                key = (pdb_id, policy_name, seed)
                if key in already_done:
                    print(f"[{done}/{total}] {pdb_id} {policy_name} seed={seed} -> cached", flush=True)
                    all_records.append(already_done[key])
                    continue
                print(f"[{done}/{total}] {pdb_id} {policy_name} seed={seed}", flush=True)
                rec = run_one(pdb_id, policy_name, water_policy, seed, by_pdb[pdb_id])
                all_records.append(rec)
                status = "OK" if rec["returncode"] == 0 else "FAIL"
                rmsd = (rec.get("result") or {}).get("top_pose_rmsd_A") if rec.get("result") else None
                outcome = (rec.get("result") or {}).get("outcome") if rec.get("result") else None
                print(f"    -> {status} rmsd={rmsd} outcome={outcome}", flush=True)

    n_fail = sum(1 for r in all_records if r["returncode"] != 0)
    print(f"\n{len(all_records) - n_fail}/{len(all_records)} runs completed without a preparation error.")
    if n_fail:
        fails = [(r["pdb_id"], r["policy"], r["seed"]) for r in all_records if r["returncode"] != 0]
        print(f"Failures: {fails}")


if __name__ == "__main__":
    main()
