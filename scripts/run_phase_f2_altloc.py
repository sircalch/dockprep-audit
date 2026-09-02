"""F2 driver: real altLoc-policy comparison (PROJECT-ROADMAP.md section 19,
item F2), added 2026-08-28.

Runs the two altloc policies (highest_occupancy = project default,
lowest_occupancy = the symmetric-opposite manipulation; see
choose_altloc_conformers() in smoke_redock_case.py) x 3 predeclared seeds,
restricted to the 17 of 60 cases confirmed by the audit engine's
SITE_ALTLOC_PRESENT finding to have an alternate-location conformer within
6 A of the declared ligand -- the site-local criterion, not the whole-
structure ALTLOC_PRESENT finding (22/60) or the alternate_location stratum
label alone (15/60), following the same F1 lesson that only a site-local
feature can plausibly affect the redocking result. The 17 cases are 15 of
the 15 alternate_location-stratum cases plus 2 low_risk_control cases
(1GS4, 1SN5) that happen to also have a site-local altLoc despite being
selected for that stratum on other grounds.

water_policy is held fixed at "conservative" for every run in this
comparison, isolating the altLoc-policy manipulation rather than crossing
it with the water-policy manipulation (which would need 4x the runs for a
comparison this study is not designed to make simultaneously).

Writes one JSON result per (case, altloc_policy, seed) under
benchmark/phase-f2-altloc/raw_runs/, plus a consolidated per-case CSV.
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
CASES_CSV = REPO_ROOT / "benchmark" / "pilot_cases.csv"
MANIFEST_CSV = REPO_ROOT / "benchmark" / "pilot_manifest_frozen.csv"

SITE_ALTLOC_CASES = [
    "1E6U", "1E7S", "1M17", "1T46", "2I4H", "3D14", "3KXG",
    "3MWU", "3P0M", "3PNA", "4RJ3", "5E0J", "5F4N", "5K8S", "5TG1",
]
# 1GS4 and 1SN5 were included in an earlier version of this list (before
# 2026-08-28) due to a dockprep_audit bug that counted a disordered WATER
# oxygen (1GS4: HOH A2026) and the ligand's own alternate conformer (1SN5:
# T3 C601) as SITE_ALTLOC_PRESENT hits. Neither is a receptor
# conformer-selection question, so both are correctly excluded from the
# site-local finding after the fix (dockprep_audit v0.2.0,
# src/dockprep_audit/audit.py); SITE_ALTLOC_PRESENT now recovers exactly
# the 15 alternate_location-stratum cases, matching SITE_METAL_PRESENT's
# clean 15/15 recovery of metal_or_cofactor.
#
# Their two altloc-policy runs are NOT necessarily identical, though, and
# should not be assumed to be: 1GS4's receptor is confirmed byte-identical
# between policies (its only altloc, the water, never enters
# extract_receptor_atoms()'s output). 1SN5 is NOT identical -- its raw
# structure contains a second T3 copy (chain D, resseq 602) that
# extract_receptor_atoms() retains as part of the receptor (as an
# HETATM-promoted-to-ATOM residue, independently of the declared C601
# ligand instance and independently of the 6 A site-local radius), and that
# copy's own altloc conformer IS resolved differently by the two policies,
# producing a real, non-trivial difference in the redocked RMSD (observed:
# ~2.4 A at seed 42). Why this second copy is retained as "receptor" at all
# is a separate, not-yet-investigated question -- flagged in
# PROJECT-ROADMAP.md rather than resolved here. Either way, 1SN5 is
# correctly excluded from the primary n=15 SITE_ALTLOC_PRESENT group (no
# protein conformer sits within 6 A of the docked ligand instance), so this
# does not affect Table 3/4-equivalent F2 statistics -- it only means the
# "should be a trivial no-op" framing that applies to 1GS4 does not also
# apply to 1SN5, and scripts/run_phase_f2_aggregate.py's confirmatory-check
# section reports both cases' actual values rather than assuming either
# outcome.
ALTLOC_POLICIES = ["highest_occupancy", "lowest_occupancy"]
SEEDS = [42, 123, 2024]
EXHAUSTIVENESS = 32
WATER_POLICY = "conservative"
OUT_DIR = REPO_ROOT / "benchmark" / "phase-f2-altloc"
PYTHON = REPO_ROOT / ".venv" / "Scripts" / "python.exe"


def run_one(pdb_id: str, altloc_policy: str, seed: int) -> dict:
    out_subdir = OUT_DIR / altloc_policy / f"seed{seed}"
    result = subprocess.run(
        [str(PYTHON if PYTHON.exists() else sys.executable),
         str(REPO_ROOT / "scripts" / "smoke_redock_case.py"), pdb_id,
         "--out-dir", str(out_subdir), "--water-policy", WATER_POLICY,
         "--altloc-policy", altloc_policy,
         "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(seed)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    record = {"pdb_id": pdb_id, "altloc_policy": altloc_policy, "seed": seed, "returncode": result.returncode}
    if result.returncode == 0:
        try:
            record["result"] = json.loads(result.stdout[result.stdout.index("{"):])
        except Exception:
            record["result"] = None
            record["stdout_tail"] = result.stdout[-2000:]
    else:
        record["stderr_tail"] = result.stderr[-2000:]
    (OUT_DIR / "raw_runs").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "raw_runs" / f"{pdb_id}_{altloc_policy}_seed{seed}.json").write_text(
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
                already_done[(rec["pdb_id"], rec["altloc_policy"], rec["seed"])] = rec
            except Exception:
                pass

    all_records = []
    total = len(SITE_ALTLOC_CASES) * len(ALTLOC_POLICIES) * len(SEEDS)
    done = 0
    for pdb_id in SITE_ALTLOC_CASES:
        for altloc_policy in ALTLOC_POLICIES:
            for seed in SEEDS:
                done += 1
                key = (pdb_id, altloc_policy, seed)
                if key in already_done:
                    print(f"[{done}/{total}] {pdb_id} {altloc_policy} seed={seed} -> cached", flush=True)
                    all_records.append(already_done[key])
                    continue
                print(f"[{done}/{total}] {pdb_id} {altloc_policy} seed={seed}", flush=True)
                rec = run_one(pdb_id, altloc_policy, seed)
                all_records.append(rec)
                status = "OK" if rec["returncode"] == 0 else "FAIL"
                rmsd = (rec.get("result") or {}).get("top_pose_rmsd_A") if rec.get("result") else None
                outcome = (rec.get("result") or {}).get("outcome") if rec.get("result") else None
                print(f"    -> {status} rmsd={rmsd} outcome={outcome}", flush=True)

    summary_rows = []
    for pdb_id in SITE_ALTLOC_CASES:
        for altloc_policy in ALTLOC_POLICIES:
            recs = [r for r in all_records if r["pdb_id"] == pdb_id and r["altloc_policy"] == altloc_policy]
            top1_rmsds = [r["result"]["top_pose_rmsd_A"] for r in recs if r.get("result")]
            best_rmsds = [r["result"]["best_pose_rmsd_A"] for r in recs if r.get("result")]
            outcomes = [r["result"]["outcome"] for r in recs if r.get("result")]
            n_fail_prep = sum(1 for r in recs if r["returncode"] != 0)
            summary_rows.append({
                "pdb_id": pdb_id, "altloc_policy": altloc_policy,
                "n_runs": len(recs), "n_preparation_failures": n_fail_prep,
                "median_top1_rmsd_A": round(statistics.median(top1_rmsds), 3) if top1_rmsds else "",
                "median_best_pose_rmsd_A": round(statistics.median(best_rmsds), 3) if best_rmsds else "",
                "success_fraction_top1": round(sum(1 for o in outcomes if o == "success") / len(outcomes), 3) if outcomes else "",
                "success_fraction_best": round(sum(1 for r in best_rmsds if r <= 2.0) / len(best_rmsds), 3) if best_rmsds else "",
                "outcome_by_seed": ";".join(outcomes),
            })

    with (OUT_DIR / "phase_f2_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nWrote {OUT_DIR / 'phase_f2_summary.csv'}")
    n_fail = sum(1 for r in all_records if r["returncode"] != 0)
    print(f"{len(all_records) - n_fail}/{len(all_records)} runs completed without a preparation error.")


if __name__ == "__main__":
    main()
