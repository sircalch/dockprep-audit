"""F9 (PROJECT-ROADMAP.md section 19, optional item): robustness check of
the water-policy finding (Section 3.5) against a richer explicit-water
representation. The 'conservative' bridging-water policy used throughout
Fase C-F models each kept water as a bare, undirected oxygen point charge
(TIP3P O charge, no hydrogens). This script re-runs the 4 cases with the
strongest water-helps effect in the main 60-case cohort (6ASH, 1WBK, 1CVZ,
4GID; Table 4) under 'conservative_oriented' -- the same bridging-water
selection, but each kept water is a rigid O+H+H TIP3P body with hydrogens
oriented toward the same ligand/receptor contacts that qualified it as
bridging in the first place (smoke_redock_case.py,
_oriented_hydrogen_positions).

Same 3 seeds, same exhaustiveness (32), same box policy as the main cohort
run (run_phase_c_60.py) -- nothing about the redocking protocol changes,
only the water representation.
"""

from __future__ import annotations

import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CASES = ["6ASH", "1WBK", "1CVZ", "4GID"]
SEEDS = [42, 123, 2024]
EXHAUSTIVENESS = 32
OUT_DIR = REPO_ROOT / "benchmark" / "phase-f9-oriented-water"


def run_one(pdb_id: str, seed: int) -> dict:
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "smoke_redock_case.py"), pdb_id,
         "--out-dir", str(OUT_DIR / "conservative_oriented" / f"seed{seed}"),
         "--water-policy", "conservative_oriented",
         "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(seed)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    record = {"pdb_id": pdb_id, "policy": "conservative_oriented", "seed": seed,
              "returncode": result.returncode}
    if result.returncode == 0:
        try:
            record["result"] = json.loads(result.stdout[result.stdout.index("{"):])
        except Exception:
            record["result"] = None
            record["stdout_tail"] = result.stdout[-2000:]
    else:
        record["stderr_tail"] = result.stderr[-2000:]
    (OUT_DIR / "raw_runs").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "raw_runs" / f"{pdb_id}_conservative_oriented_seed{seed}.json").write_text(
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
                already_done[(rec["pdb_id"], rec["seed"])] = rec
            except Exception:
                pass

    all_records = []
    total = len(CASES) * len(SEEDS)
    done = 0
    for pdb_id in CASES:
        for seed in SEEDS:
            done += 1
            key = (pdb_id, seed)
            if key in already_done:
                print(f"[{done}/{total}] {pdb_id} seed={seed} -> cached", flush=True)
                all_records.append(already_done[key])
                continue
            print(f"[{done}/{total}] {pdb_id} seed={seed}", flush=True)
            rec = run_one(pdb_id, seed)
            all_records.append(rec)
            status = "OK" if rec["returncode"] == 0 else "FAIL"
            rmsd = (rec.get("result") or {}).get("top_pose_rmsd_A") if rec.get("result") else None
            n_waters = (rec.get("result") or {}).get("bridging_waters_kept") if rec.get("result") else None
            print(f"    -> {status} rmsd={rmsd} waters={n_waters}", flush=True)

    n_fail = sum(1 for r in all_records if r["returncode"] != 0)
    print(f"\n{len(all_records) - n_fail}/{len(all_records)} runs completed without a preparation error.")
    if n_fail:
        print("Failures:", [(r["pdb_id"], r["seed"]) for r in all_records if r["returncode"] != 0])


if __name__ == "__main__":
    main()
