"""Freeze the pilot cohort manifest: consolidate provenance (URL, checksum),
receptor chain, ligand identity, and inclusion basis for each case into one
versioned, checksummed file. Once frozen, no case in this file may be added,
removed, or reassigned to a different stratum without a new, dated entry in
benchmark/PILOT_SELECTION.md explaining why -- that is what "frozen" means
here: the cohort composition is locked before any docking result is
observed (Fase A of PROJECT-ROADMAP.md section 9).

This script only consolidates already-produced audit artifacts; it makes no
new structural or scientific judgment calls of its own.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.request
from pathlib import Path
from urllib.error import URLError

REPO_ROOT = Path(__file__).resolve().parent.parent
CASES = REPO_ROOT / "benchmark" / "pilot_cases.csv"
CANDIDATES = REPO_ROOT / "benchmark" / "pilot-inventory" / "cohort_candidates.csv"
ELIGIBILITY = REPO_ROOT / "benchmark" / "pilot-eligibility" / "eligibility_results.csv"
PREVIOUS_EXCLUSIONS = REPO_ROOT / "benchmark" / "previous-study-exclusions.csv"
OUT_CSV = REPO_ROOT / "benchmark" / "pilot_manifest_frozen.csv"
OUT_SHA256 = REPO_ROOT / "benchmark" / "pilot_manifest_frozen.sha256.txt"

FIELDS = [
    "pdb_id", "stratum", "status",
    "source_url", "source_sha256", "experimental_method", "resolution_angstrom", "title",
    "ligand_component_id", "ligand_chain", "ligand_resseq",
    "receptor_chains", "receptor_chains_basis", "receptor_chains_override_reason",
    "target_uniprot", "shares_target_with_prior_pdb", "shares_target_prior_study",
]

RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
TARGET_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { uniprot_ids }
    }
  }
}
"""


def fetch_uniprot_ids(pdb_ids: list[str]) -> dict[str, set[str]]:
    """Independence check (pendiente #3, seccion 8), decided 2026-08-21: the
    previous-study exclusion registry (benchmark/previous-study-exclusions.csv)
    is keyed by PDB ID only, so it cannot catch a *different* PDB entry of
    the *same* biological target (protein) already used in a prior
    manuscript. UniProt accession is the unambiguous way to check that,
    independent of which specific crystal form/PDB ID was deposited.
    """
    result: dict[str, set[str]] = {}
    for i in range(0, len(pdb_ids), 50):
        batch = pdb_ids[i:i + 50]
        payload = json.dumps({"query": TARGET_QUERY, "variables": {"ids": batch}}).encode()
        req = urllib.request.Request(RCSB_GRAPHQL_URL, data=payload, headers={"Content-Type": "application/json"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    data = json.load(resp)
                break
            except (URLError, OSError) as exc:
                if attempt == 4:
                    raise
                wait = 5 * (attempt + 1)
                print(f"  retry {attempt + 1}/5 after {exc} (waiting {wait}s)")
                time.sleep(wait)
        for entry in data.get("data", {}).get("entries") or []:
            if entry is None:
                continue
            uniprots = set()
            for pe in entry.get("polymer_entities") or []:
                ids = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []
                uniprots.update(ids)
            result[entry["rcsb_id"]] = uniprots
    return result


def main() -> None:
    with CASES.open(newline="", encoding="utf-8") as f:
        cases = {row["pdb_id"].upper(): row for row in csv.DictReader(f)}
    with CANDIDATES.open(newline="", encoding="utf-8") as f:
        candidates = {row["pdb_id"].upper(): row for row in csv.DictReader(f)}
    with ELIGIBILITY.open(newline="", encoding="utf-8") as f:
        eligibility = {row["pdb_id"].upper(): row for row in csv.DictReader(f)}
    with PREVIOUS_EXCLUSIONS.open(newline="", encoding="utf-8") as f:
        excluded_rows = list(csv.DictReader(f))

    excluded_ids = sorted({row["pdb_id"].upper() for row in excluded_rows})
    excluded_study_by_id = {row["pdb_id"].upper(): row["source_study"] for row in excluded_rows}
    try:
        uniprots = fetch_uniprot_ids(list(cases.keys()) + excluded_ids)
    except (URLError, OSError) as exc:
        print(f"WARNING: RCSB target-independence lookup unavailable ({exc}). "
              "Freezing WITHOUT target_uniprot/shares_target columns populated -- "
              "rerun this script later to backfill that transparency data; it does not "
              "block freezing the core provenance (URL/checksum/chain/ligand).")
        uniprots = {}

    rows = []
    for pdb_id, case in cases.items():
        cand = candidates[pdb_id]
        elig = eligibility[pdb_id]
        override_reason = elig.get("receptor_chains_override_reason", "").strip()

        my_uniprots = uniprots.get(pdb_id, set())
        overlap_pdb, overlap_study = "", ""
        for ex_id in excluded_ids:
            if my_uniprots & uniprots.get(ex_id, set()):
                overlap_pdb, overlap_study = ex_id, excluded_study_by_id[ex_id]
                break  # one documented example is enough for transparency; not exhaustive

        rows.append({
            "pdb_id": pdb_id,
            "stratum": case["stratum"],
            "status": "frozen",
            "source_url": cand["source_url"],
            "source_sha256": cand["sha256"],
            "experimental_method": cand["experimental_method"],
            "resolution_angstrom": cand["resolution_angstrom"],
            "title": cand["title"],
            "ligand_component_id": case["ligand_component_id"],
            "ligand_chain": case["ligand_chain"],
            "ligand_resseq": case["ligand_resseq"],
            "receptor_chains": elig["receptor_chains"],
            "receptor_chains_basis": "user_override" if override_reason else "geometric_6A_rule",
            "receptor_chains_override_reason": override_reason,
            "target_uniprot": ";".join(sorted(my_uniprots)),
            "shares_target_with_prior_pdb": overlap_pdb,
            "shares_target_prior_study": overlap_study,
        })

    rows.sort(key=lambda r: r["pdb_id"])

    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    digest = hashlib.sha256(OUT_CSV.read_bytes()).hexdigest()
    OUT_SHA256.write_text(f"{digest}  {OUT_CSV.name}\n", encoding="utf-8")

    print(f"Froze {len(rows)} cases into {OUT_CSV}")
    print(f"SHA-256: {digest}")


if __name__ == "__main__":
    main()
