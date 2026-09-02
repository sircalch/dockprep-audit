"""Download a declared PDB pilot set and create an auditable inventory.

This script deliberately does not choose a "correct" preparation. It retains
the source file checksum and records only rule-based features detected by
DockPrep Audit. Candidate selection must be documented before docking begins.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from dockprep_audit import audit_pdb

EXCLUDED_PREVIOUS_STUDIES = {"1STP", "1B9V", "1HVR", "1IEP", "3PTB", "3CJO", "3D4Q", "2I78"}


def category(report: dict) -> str:
    codes = {finding["code"] for finding in report["findings"]}
    if "ALTLOC_PRESENT" in codes:
        return "alternate_location"
    if "METAL_PRESENT" in codes:
        return "metal_or_cofactor"
    if "WATERS_PRESENT" in codes:
        return "water_present"
    return "low_risk_candidate"


def download(pdb_id: str, destination: Path) -> tuple[str, str]:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    payload = urlopen(url, timeout=30).read()
    destination.write_bytes(payload)
    return url, hashlib.sha256(payload).hexdigest()


def entry_metadata(pdb_id: str) -> dict[str, str]:
    """Fetch concise RCSB metadata used for transparent eligibility screening."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    record = json.loads(urlopen(url, timeout=30).read().decode("utf-8"))
    resolution = record.get("rcsb_entry_info", {}).get("resolution_combined", [])
    methods = [item.get("method", "") for item in record.get("exptl", [])]
    return {
        "experimental_method": "; ".join(methods),
        "resolution_angstrom": "; ".join(str(value) for value in resolution),
        "title": record.get("struct", {}).get("title", "").replace("\n", " "),
    }


def hetero_residues(path: Path) -> str:
    residues = {
        line[17:20].strip()
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("HETATM") and line[17:20].strip() not in {"HOH", "WAT", "H2O", "DOD"}
    }
    return "; ".join(sorted(residues))


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a checksum-tracked PDB pilot inventory.")
    parser.add_argument("ids", nargs="+", help="PDB identifiers to retrieve")
    parser.add_argument("--output", type=Path, default=Path("benchmark/pilot-inventory"))
    parser.add_argument("--exclusions", type=Path, default=Path("benchmark/previous-study-exclusions.csv"))
    args = parser.parse_args()

    excluded = set(EXCLUDED_PREVIOUS_STUDIES)
    if args.exclusions.is_file():
        with args.exclusions.open(newline="", encoding="utf-8") as handle:
            excluded.update(row["pdb_id"].upper() for row in csv.DictReader(handle))

    output = args.output
    raw = output / "raw-pdb"
    raw.mkdir(parents=True, exist_ok=True)
    records = []
    for raw_id in args.ids:
        pdb_id = raw_id.upper()
        record = {
            "pdb_id": pdb_id, "selected": False, "reason": "candidate", "source_url": "", "sha256": "",
            "category": "", "experimental_method": "", "resolution_angstrom": "", "hetero_residues": "", "title": "", "status": "",
        }
        if pdb_id in excluded:
            record.update(status="excluded", reason="present in a previous-study registry")
            records.append(record)
            continue
        try:
            path = raw / f"{pdb_id}.pdb"
            url, digest = download(pdb_id, path)
            report = audit_pdb(path)
            metadata = entry_metadata(pdb_id)
            (output / "reports").mkdir(exist_ok=True)
            (output / "reports" / f"{pdb_id}.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
            record.update(source_url=url, sha256=digest, category=category(report), hetero_residues=hetero_residues(path), status="retrieved", **metadata)
        except (HTTPError, URLError, TimeoutError) as exc:
            record.update(status="unavailable", reason=str(exc))
        records.append(record)

    fields = ["pdb_id", "selected", "reason", "source_url", "sha256", "category", "experimental_method", "resolution_angstrom", "hetero_residues", "title", "status"]
    with (output / "cohort_candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    print(f"Wrote {len(records)} candidate records to {output / 'cohort_candidates.csv'}")


if __name__ == "__main__":
    main()
