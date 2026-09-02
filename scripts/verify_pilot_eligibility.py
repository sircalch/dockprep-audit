"""Verify declared ligands and binding-site relevance for the pilot cohort."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

WATERS = {"HOH", "WAT", "H2O", "DOD"}
METALS = {"AL", "CA", "CD", "CO", "CU", "FE", "HG", "K", "MG", "MN", "NA", "NI", "ZN"}


def parse_atom(line: str) -> dict:
    return {
        "record": line[0:6].strip(), "name": line[12:16].strip(), "altloc": line[16:17].strip(),
        "resname": line[17:20].strip(), "chain": line[21:22].strip() or "_", "resseq": line[22:26].strip(),
        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
        "element": line[76:78].strip().upper(),
    }


def distance(a: dict, b: dict) -> float:
    return math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))


def min_distance(group_a: list[dict], group_b: list[dict]) -> float | None:
    if not group_a or not group_b:
        return None
    return min(distance(a, b) for a in group_a for b in group_b)


def main() -> None:
    parser = argparse.ArgumentParser(description="Check pilot ligands and site-local structural features.")
    parser.add_argument("--cases", type=Path, default=Path("benchmark/pilot_cases.csv"))
    parser.add_argument("--raw", type=Path, default=Path("benchmark/pilot-inventory/raw-pdb"))
    parser.add_argument("--output", type=Path, default=Path("benchmark/pilot-eligibility"))
    parser.add_argument("--overrides", type=Path, default=Path("benchmark/receptor_chain_overrides.csv"),
                         help="Explicit, documented overrides of the 6 A geometric receptor_chains rule "
                              "(e.g. real biological dimers whose interface matters even when this "
                              "specific bound ligand instance sits closer to one monomer). Never applied "
                              "silently: the raw geometric finding is kept in receptor_chains_geometric.")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    overrides = {}
    if args.overrides.exists():
        with args.overrides.open(newline="", encoding="utf-8") as handle:
            overrides = {row["pdb_id"].upper(): row for row in csv.DictReader(handle)}

    results = []
    details = {}
    with args.cases.open(newline="", encoding="utf-8") as handle:
        cases = list(csv.DictReader(handle))
    for case in cases:
        pdb_id = case["pdb_id"].upper()
        path = args.raw / f"{pdb_id}.pdb"
        atoms = [parse_atom(line) for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith(("ATOM  ", "HETATM"))]
        matching_ligands = [a for a in atoms if a["record"] == "HETATM" and a["resname"] == case["ligand_component_id"]]
        declared_chain = case.get("ligand_chain", "").strip()
        declared_resseq = case.get("ligand_resseq", "").strip()
        ligand = [
            a for a in matching_ligands
            if (not declared_chain or a["chain"] == declared_chain) and (not declared_resseq or a["resseq"] == declared_resseq)
        ]
        protein = [a for a in atoms if a["record"] == "ATOM"]
        altloc = [a for a in protein if a["altloc"]]
        metals = [a for a in atoms if a["element"] in METALS]
        waters = [a for a in atoms if a["resname"] in WATERS]
        ligand_instances = sorted({f"{a['chain']}:{a['resseq']}" for a in ligand})
        receptor_chains_geometric = sorted({a["chain"] for a in protein if any(distance(a, lig) <= 6.0 for lig in ligand)}) if ligand else []
        override = overrides.get(pdb_id)
        if override:
            receptor_chains = override["receptor_chains"].split(";")
            receptor_chains_override_reason = override["reason"]
        else:
            receptor_chains = receptor_chains_geometric
            receptor_chains_override_reason = ""
        close_waters = sorted({f"{a['chain']}:{a['resseq']}" for a in waters if any(distance(a, lig) <= 4.0 for lig in ligand)}) if ligand else []
        altloc_distance = min_distance(altloc, ligand)
        metal_distance = min_distance(metals, ligand)
        issues = []
        if not ligand:
            issues.append("declared ligand not found")
        if len(ligand_instances) != 1:
            issues.append(f"expected one ligand instance; found {len(ligand_instances)}")
        if not receptor_chains:
            issues.append("no receptor chain within 6 A of ligand")
        site_feature = {
            "alternate_location": altloc_distance is not None and altloc_distance <= 6.0,
            "metal_or_cofactor": metal_distance is not None and metal_distance <= 6.0,
            "water_policy": bool(close_waters),
            "low_risk_control": not (altloc_distance is not None and altloc_distance <= 6.0) and not (metal_distance is not None and metal_distance <= 6.0),
        }[case["stratum"]]
        if not site_feature:
            issues.append("declared stratum feature is not binding-site relevant")
        result = {
            **case, "ligand_instances": len(ligand_instances), "ligand_atoms": len(ligand),
            "receptor_chains": ";".join(receptor_chains),
            "receptor_chains_geometric": ";".join(receptor_chains_geometric),
            "receptor_chains_override_reason": receptor_chains_override_reason,
            "binding_site_waters_4A": len(close_waters),
            "nearest_altloc_A": "" if altloc_distance is None else f"{altloc_distance:.3f}",
            "nearest_metal_A": "" if metal_distance is None else f"{metal_distance:.3f}",
            "eligibility_status": "pass" if not issues else "review", "issues": "; ".join(issues),
        }
        results.append(result)
        details[pdb_id] = {"ligand_instances": ligand_instances, "receptor_chains": receptor_chains,
                            "receptor_chains_geometric": receptor_chains_geometric, "binding_site_waters": close_waters, "issues": issues}

    fields = list(results[0])
    with (args.output / "eligibility_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    (args.output / "site_details.json").write_text(json.dumps(details, indent=2), encoding="utf-8")
    print(f"Verified {len(results)} cases: {sum(row['eligibility_status'] == 'pass' for row in results)} pass, {sum(row['eligibility_status'] == 'review' for row in results)} review")


if __name__ == "__main__":
    main()
