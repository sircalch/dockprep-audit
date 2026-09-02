"""F7 step 4: one substitution pass for cases that failed technical
validation on the first redocking pass (mirrors the main 60-case cohort's
Table S2 substitution protocol exactly -- same rule: a preparation-stage
crash disqualifies a candidate, a scoring/sampling failure does not, and
every replacement must independently pass the same screen as everyone
else, here: SITE_BRIDGING_WATER_PRESENT, from the same pre-frozen,
already-screened candidate pool in screening_results.csv, never a newly
searched candidate chosen after seeing this failure).

Run once, after run_phase's first pass is complete and all failures are
known -- not iteratively per failure, to avoid any appearance of tuning
substitutions to a desired result.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "benchmark" / "f7-external-validation"
RAW_DIR = OUT_DIR / "raw-pdb"
MAX_PER_LIGAND = 2


def find_failed_cases() -> list[str]:
    raw_runs = list((OUT_DIR / "raw_runs").glob("*.json"))
    by_case_ok = {}
    for p in raw_runs:
        rec = json.loads(p.read_text(encoding="utf-8"))
        pdb_id = rec["pdb_id"]
        by_case_ok.setdefault(pdb_id, []).append(rec["returncode"] == 0)
    failed = sorted(pdb for pdb, oks in by_case_ok.items() if not all(oks))
    return failed


def main() -> None:
    manifest_rows = list(csv.DictReader(open(OUT_DIR / "manifest.csv", encoding="utf-8")))
    manifest_by_pdb = {r["pdb_id"]: r for r in manifest_rows}
    current_ids = set(manifest_by_pdb)
    per_ligand_count = {}
    for r in manifest_rows:
        per_ligand_count[r["ligand_component_id"]] = per_ligand_count.get(r["ligand_component_id"], 0) + 1

    failed = find_failed_cases()
    print(f"{len(failed)} cases failed technical validation: {failed}")

    screening_rows = list(csv.DictReader(open(OUT_DIR / "screening_results.csv", encoding="utf-8")))
    qualifying = [r for r in screening_rows if r["status"] == "QUALIFIES"]
    qualifying.sort(key=lambda r: r["pdb_id"])
    pool = [r for r in qualifying if r["pdb_id"] not in current_ids]

    substitutes = []
    for r in pool:
        if len(substitutes) >= len(failed):
            break
        lig = r["ligand_component_id"]
        if per_ligand_count.get(lig, 0) >= MAX_PER_LIGAND:
            continue
        substitutes.append(r)
        per_ligand_count[lig] = per_ligand_count.get(lig, 0) + 1

    print(f"{len(substitutes)} substitutes selected: {[r['pdb_id'] for r in substitutes]}")

    # Remove failed cases from the manifest, add substitutes, keep frozen provenance for both.
    kept_rows = [manifest_by_pdb[p] for p in manifest_by_pdb if p not in failed]
    sub_manifest_rows = []
    for r in substitutes:
        pdb_id = r["pdb_id"]
        raw_path = RAW_DIR / f"{pdb_id}.pdb"
        if not raw_path.exists():
            import urllib.request
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", raw_path)
        sha256 = hashlib.sha256(raw_path.read_bytes()).hexdigest()
        sub_manifest_rows.append({
            "pdb_id": pdb_id,
            "ligand_component_id": r["ligand_component_id"],
            "ligand_chain": r["ligand_chain"],
            "ligand_resseq": r["ligand_resseq"],
            "receptor_chains": r["ligand_chain"],
            "source_url": f"https://files.rcsb.org/download/{pdb_id}.pdb",
            "source_sha256": sha256,
            "findings": r["findings"],
        })

    new_manifest_rows = kept_rows + sub_manifest_rows
    new_manifest_rows.sort(key=lambda r: r["pdb_id"])

    with (OUT_DIR / "manifest.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(new_manifest_rows[0].keys()))
        w.writeheader()
        w.writerows(new_manifest_rows)
    print(f"Wrote updated manifest.csv ({len(new_manifest_rows)} cases: "
          f"{len(kept_rows)} original + {len(sub_manifest_rows)} substitutes)")

    with (OUT_DIR / "substitutions.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["replaced_case", "cause", "replacement"])
        w.writeheader()
        for failed_id, sub in zip(failed, substitutes):
            w.writerow({"replaced_case": failed_id, "cause": "preparation-stage failure (Meeko/dimorphite-dl)",
                        "replacement": sub["pdb_id"]})
    print(f"Wrote substitutions.csv ({len(failed)} replacements documented)")

    combined = "\n".join(f"{r['pdb_id']},{r['source_sha256']}" for r in new_manifest_rows)
    manifest_sha256 = hashlib.sha256(combined.encode("utf-8")).hexdigest()
    (OUT_DIR / "manifest.sha256.txt").write_text(manifest_sha256 + "\n", encoding="utf-8")
    print(f"New manifest checksum: {manifest_sha256}")


if __name__ == "__main__":
    main()
