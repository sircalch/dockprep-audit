"""F7 step 1: screen RCSB candidates for the external validation cohort
(PROJECT-ROADMAP.md section 19, item F7), added 2026-08-28.

Downloads each candidate's raw PDB (RCSB search API result, resolution
<=2.2 A, X-ray, single protein chain, has characterized binding-affinity
data), excludes anything already used in this study's frozen 60-case
cohort or the two prior studies' 132-identifier registry, identifies the
largest plausible small-molecule ligand per structure (excluding water,
common crystallization additives, and bare ions), and runs
dockprep_audit.audit_pdb() with that ligand to check SITE_BRIDGING_WATER_PRESENT
-- the exact same finding used to define the water_policy stratum in the
main 60-case cohort (Section 2.1/2.2). Enrichment for this finding is
deliberate (F7's stated purpose is validating the water-policy result
specifically), not an attempt to replicate the full four-stratum design.

Writes benchmark/f7-external-validation/screening_results.csv (one row per
candidate attempted) and prints a running qualifying-candidate count.
"""

from __future__ import annotations

import csv
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "src"))

import dockprep_audit as audit  # noqa: E402

OUT_DIR = REPO_ROOT / "benchmark" / "f7-external-validation"
RAW_DIR = OUT_DIR / "raw-pdb"
RAW_DIR.mkdir(parents=True, exist_ok=True)

WATER_RESIDUES = {"HOH", "WAT", "H2O", "DOD"}
BUFFER_BLOCKLIST = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "ACT", "TRS", "BME", "DMS",
    "CIT", "IMD", "MPD", "BOG", "CO3", "NO3", "ACY", "FMT", "EOH", "MES",
    "TAM", "BTB", "P6G", "1PE", "15P", "PGE", "MRD", "IPA", "NH4", "UNX",
    "UNL", "EPE", "CAC", "HED", "PGO", "OXL", "SIN", "MLA", "PLM", "MYR",
    "OCT", "DTT", "DTU", "TCE", "PEO", "PE4", "6JZ", "SPD", "SPM", "PUT",
}
ION_BLOCKLIST = {
    "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CD", "CU", "CO", "NI",
    "HG", "AL", "BA", "SR", "CS", "LI", "RB", "PB", "AG", "AU", "IOD", "BR",
}
MIN_LIGAND_HEAVY_ATOMS = 8


def download_pdb(pdb_id: str) -> Path | None:
    out_path = RAW_DIR / f"{pdb_id}.pdb"
    if out_path.exists():
        return out_path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
        out_path.write_bytes(data)
        return out_path
    except urllib.error.HTTPError:
        return None
    except Exception:
        return None


def find_best_ligand(pdb_path: Path) -> dict | None:
    """Largest non-water, non-buffer, non-ion HETATM component; returns
    {"component_id", "chain", "resseq", "n_heavy"} or None."""
    residues: dict[tuple[str, str, str], int] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname in WATER_RESIDUES or resname in BUFFER_BLOCKLIST or resname in ION_BLOCKLIST:
            continue
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        key = (resname, chain, resseq)
        residues[key] = residues.get(key, 0) + 1

    candidates = [(k, v) for k, v in residues.items() if v >= MIN_LIGAND_HEAVY_ATOMS]
    if not candidates:
        return None
    (resname, chain, resseq), n_heavy = max(candidates, key=lambda kv: kv[1])
    return {"component_id": resname, "chain": chain, "resseq": resseq, "n_heavy": n_heavy}


def main() -> None:
    candidates = (OUT_DIR / "candidates.txt").read_text(encoding="utf-8").splitlines()
    exclusions = set()
    for row in csv.DictReader(open(REPO_ROOT / "benchmark" / "previous-study-exclusions.csv", encoding="utf-8")):
        exclusions.add(row["pdb_id"])
    for row in csv.DictReader(open(REPO_ROOT / "benchmark" / "pilot_cases.csv", encoding="utf-8")):
        exclusions.add(row["pdb_id"])

    print(f"{len(candidates)} raw candidates, {len(exclusions)} excluded identifiers")

    results = []
    n_qualifying = 0
    for i, pdb_id in enumerate(candidates, 1):
        if pdb_id in exclusions:
            results.append({"pdb_id": pdb_id, "status": "excluded_prior_study"})
            continue

        pdb_path = download_pdb(pdb_id)
        if pdb_path is None:
            results.append({"pdb_id": pdb_id, "status": "download_failed"})
            continue

        ligand = find_best_ligand(pdb_path)
        if ligand is None:
            results.append({"pdb_id": pdb_id, "status": "no_plausible_ligand"})
            continue

        try:
            report = audit.audit_pdb(str(pdb_path), ligand={
                "component_id": ligand["component_id"], "chain": ligand["chain"], "resseq": ligand["resseq"],
            })
        except Exception as exc:
            results.append({"pdb_id": pdb_id, "status": f"audit_error: {exc}"})
            continue

        codes = {f["code"] for f in report["findings"]}
        has_bridging_water = "SITE_BRIDGING_WATER_PRESENT" in codes
        status = "QUALIFIES" if has_bridging_water else "no_bridging_water"
        if has_bridging_water:
            n_qualifying += 1
        results.append({
            "pdb_id": pdb_id, "status": status,
            "ligand_component_id": ligand["component_id"], "ligand_chain": ligand["chain"],
            "ligand_resseq": ligand["resseq"], "ligand_n_heavy": ligand["n_heavy"],
            "findings": ";".join(sorted(codes)),
        })

        if i % 50 == 0:
            print(f"  [{i}/{len(candidates)}] {n_qualifying} qualifying so far", flush=True)

    fieldnames = ["pdb_id", "status", "ligand_component_id", "ligand_chain", "ligand_resseq",
                  "ligand_n_heavy", "findings"]
    with (OUT_DIR / "screening_results.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    print(f"\nDone. {n_qualifying} qualifying candidates (SITE_BRIDGING_WATER_PRESENT) "
          f"out of {len(candidates)} raw candidates.")


if __name__ == "__main__":
    main()
