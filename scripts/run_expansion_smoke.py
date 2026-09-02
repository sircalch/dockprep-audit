"""Technical smoke test for the 48 new Fase E cases (pendiente in
PROJECT-ROADMAP.md section 9): one run per case, default policy
(water_policy=none, exhaustiveness=8, seed=42) -- the same purpose as the
original 12-case smoke test: catch preparation failures (e.g. the 1CPS
TYR:204 impossible-geometry case) before committing to the full multi-seed,
multi-policy Fase C compute run. This is diagnostic, not a scientific
result.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "benchmark" / "expansion-smoke"
CASES_CSV = REPO_ROOT / "benchmark" / "expansion_cases_for_eligibility.csv"
ELIGIBILITY_CSV = REPO_ROOT / "benchmark" / "expansion-eligibility" / "eligibility_results.csv"


def main() -> None:
    with (REPO_ROOT / "benchmark" / "expansion_cases.csv").open(newline="", encoding="utf-8") as f:
        cases = [row["pdb_id"] for row in csv.DictReader(f)]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_path = OUT_DIR / "smoke_log.jsonl"
    summary = []
    for i, pdb_id in enumerate(cases, 1):
        print(f"[{i}/{len(cases)}] {pdb_id}", flush=True)
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "smoke_redock_case.py"), pdb_id,
             "--cases", str(CASES_CSV), "--eligibility", str(ELIGIBILITY_CSV),
             "--out-dir", str(OUT_DIR)],
            cwd=REPO_ROOT, capture_output=True, text=True,
        )
        record = {"pdb_id": pdb_id, "returncode": result.returncode}
        if result.returncode == 0:
            try:
                record["result"] = json.loads(result.stdout[result.stdout.index("{"):])
            except Exception:
                record["result"] = None
                record["stdout_tail"] = result.stdout[-1500:]
        else:
            record["stderr_tail"] = result.stderr[-1500:]
        with log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(record) + "\n")
        summary.append(record)
        status = "OK" if record["returncode"] == 0 else "FAIL"
        rmsd = (record.get("result") or {}).get("top_pose_rmsd_A") if record.get("result") else None
        outcome = (record.get("result") or {}).get("outcome") if record.get("result") else None
        print(f"    -> {status} rmsd={rmsd} outcome={outcome}", flush=True)

    n_ok = sum(1 for r in summary if r["returncode"] == 0)
    print(f"\n{n_ok}/{len(summary)} cases prepared and docked without error.")
    fails = [r["pdb_id"] for r in summary if r["returncode"] != 0]
    if fails:
        print("Preparation/docking failures:", fails)


if __name__ == "__main__":
    main()
