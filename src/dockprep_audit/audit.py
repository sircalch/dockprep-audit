"""Rule-based PDB audit that preserves source structure content."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

WATER_RESIDUES = {"HOH", "WAT", "H2O", "DOD"}
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",
}
METALS = {"AL", "CA", "CD", "CO", "CU", "FE", "HG", "K", "MG", "MN", "NA", "NI", "ZN"}

# Site-local finding thresholds. WATER_BRIDGE_DISTANCE_A matches Eq. 1 of the
# companion manuscript exactly (symmetric ligand/receptor bridging-water
# criterion already used to stratify the 60-case benchmark); SITE_RADIUS_A
# matches the cohort's own spatial-relevance inclusion rule (Section 2.2).
WATER_BRIDGE_DISTANCE_A = 3.0
SITE_RADIUS_A = 6.0


def _atom(line: str) -> dict[str, Any]:
    """Parse the fixed-width PDB fields needed by the audit."""
    try:
        coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        coord = None
    return {
        "record": line[0:6].strip(),
        "atom_name": line[12:16].strip(),
        "altloc": line[16:17].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21:22].strip() or "_",
        "resseq": line[22:26].strip(),
        "icode": line[26:27].strip(),
        "occupancy": line[54:60].strip(),
        "element": line[76:78].strip().upper(),
        "coord": coord,
    }


def _min_dist(coords_a: list[tuple[float, float, float]], coords_b: list[tuple[float, float, float]]) -> float:
    return min((math.dist(a, b) for a in coords_a for b in coords_b), default=float("inf"))


def _finding(code: str, severity: str, message: str, count: int, examples: list[str]) -> dict[str, Any]:
    return {"code": code, "severity": severity, "message": message, "count": count, "examples": examples[:10]}


def audit_pdb(
    path: str | Path,
    ligand: dict[str, str] | None = None,
    receptor_chains: set[str] | None = None,
) -> dict[str, Any]:
    """Audit a PDB receptor and return a serialisable provenance report.

    ``ligand``, if given, identifies the declared ligand instance as
    ``{"component_id": ..., "chain": ..., "resseq": ...}`` (component_id is
    the 3-letter HETATM residue name). When present, three additional
    site-local findings are evaluated relative to that ligand: whether a
    water molecule bridges the ligand and the receptor (Eq. 1 of the
    companion manuscript), and whether an altLoc or metal finding lies
    within the same site-local radius used to build the benchmark cohort.
    These are distinct from the whole-structure ALTLOC_PRESENT/
    WATERS_PRESENT/METAL_PRESENT findings, which fire on any occurrence
    anywhere in the file regardless of ligand proximity.

    ``receptor_chains``, if given, restricts which polymer chains count as
    "the receptor" for the water-bridging and site-radius checks -- matching
    the declared biological assembly for the case (Section 2.2) rather than
    every polymer chain physically present in the deposited file, which may
    include crystallographic copies outside the modeled complex. Ignored if
    ``ligand`` is not also given.
    """
    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = [_atom(line) for line in lines if line.startswith(("ATOM  ", "HETATM"))]

    # Alternate locations on a water oxygen (disordered solvent) are a
    # water-handling question (WATERS_PRESENT), not the receptor
    # conformer-selection question this finding is about (Section 2.3);
    # excluding them here avoids double-counting the same disorder under
    # two unrelated finding codes.
    altlocs = [a for a in atoms if a["altloc"] and a["resname"] not in WATER_RESIDUES]
    waters = [a for a in atoms if a["resname"] in WATER_RESIDUES]
    nonstandard = [a for a in atoms if a["record"] == "ATOM" and a["resname"] not in STANDARD_RESIDUES]
    metals = [a for a in atoms if a["element"] in METALS]
    missing_occupancy = [a for a in atoms if not a["occupancy"]]
    missing_element = [a for a in atoms if not a["element"]]

    atom_keys = [(a["chain"], a["resseq"], a["icode"], a["resname"], a["atom_name"], a["altloc"]) for a in atoms]
    duplicate_keys = [key for key, count in Counter(atom_keys).items() if count > 1]

    findings = []
    if altlocs:
        residues = sorted({f"{a['resname']} {a['chain']}{a['resseq']}{a['icode']}" for a in altlocs})
        findings.append(_finding("ALTLOC_PRESENT", "high", "Alternate atom locations require a documented conformer-selection policy before docking.", len(altlocs), residues))
    if waters:
        residues = sorted({f"{a['resname']} {a['chain']}{a['resseq']}" for a in waters})
        findings.append(_finding("WATERS_PRESENT", "medium", "Retaining or removing waters is protocol-dependent; record the decision.", len(residues), residues))
    if nonstandard:
        residues = sorted({f"{a['resname']} {a['chain']}{a['resseq']}" for a in nonstandard})
        findings.append(_finding("NONSTANDARD_RESIDUE", "medium", "Non-standard polymer residues may need explicit treatment in preparation tools.", len(residues), residues))
    if metals:
        residues = sorted({f"{a['element']} {a['resname']} {a['chain']}{a['resseq']}" for a in metals})
        findings.append(_finding("METAL_PRESENT", "medium", "Metal coordination and charge treatment should be documented.", len(residues), residues))
    if missing_occupancy:
        findings.append(_finding("MISSING_OCCUPANCY", "low", "Some atoms have no occupancy value in the PDB record.", len(missing_occupancy), []))
    if missing_element:
        findings.append(_finding("MISSING_ELEMENT", "low", "Some atoms have no element annotation; downstream format conversion may infer it ambiguously.", len(missing_element), []))
    if duplicate_keys:
        examples = [f"{key[3]} {key[0]}{key[1]} atom {key[4]} altLoc {key[5] or '-'}" for key in duplicate_keys]
        findings.append(_finding("DUPLICATE_ATOM_RECORD", "high", "Duplicate atom records were found within one residue and conformer.", len(duplicate_keys), examples))

    if ligand is not None:
        lig_chain = ligand.get("chain", "")
        lig_resseq = str(ligand.get("resseq", "")).strip()
        lig_component = ligand.get("component_id", "")
        ligand_atoms = [
            a for a in atoms
            if a["resname"] == lig_component and a["chain"] == (lig_chain or a["chain"]) and a["resseq"] == lig_resseq
            and a["coord"] is not None
        ]
        if not ligand_atoms:
            findings.append(_finding(
                "SITE_LIGAND_NOT_FOUND", "low",
                f"Declared ligand {lig_component} {lig_chain}{lig_resseq} was not found in this file; "
                "site-local findings (SITE_BRIDGING_WATER_PRESENT, SITE_ALTLOC_PRESENT, SITE_METAL_PRESENT) "
                "were not evaluated.",
                0, [],
            ))
        else:
            ligand_coords = [a["coord"] for a in ligand_atoms]

            def _in_receptor_chains(a: dict[str, Any]) -> bool:
                return receptor_chains is None or a["chain"] in receptor_chains

            receptor_coords = [
                a["coord"] for a in atoms
                if a["record"] == "ATOM" and a["coord"] is not None and _in_receptor_chains(a)
            ]
            waters_in_scope = [a for a in waters if _in_receptor_chains(a)]

            bridging_waters = []
            for a in waters_in_scope:
                if a["coord"] is None:
                    continue
                if _min_dist([a["coord"]], ligand_coords) <= WATER_BRIDGE_DISTANCE_A \
                        and _min_dist([a["coord"]], receptor_coords) <= WATER_BRIDGE_DISTANCE_A:
                    bridging_waters.append(f"{a['resname']} {a['chain']}{a['resseq']}")
            if bridging_waters:
                residues = sorted(set(bridging_waters))
                findings.append(_finding(
                    "SITE_BRIDGING_WATER_PRESENT", "medium",
                    f"A water molecule lies within {WATER_BRIDGE_DISTANCE_A:.1f} Å of both the declared "
                    "ligand and the receptor (Eq. 1); retaining or removing it is protocol-dependent and "
                    "has been shown to change the redocking result in a minority of such cases.",
                    len(residues), residues,
                ))

            def _is_ligand_atom(a: dict[str, Any]) -> bool:
                return a["resname"] == lig_component and a["chain"] == (lig_chain or a["chain"]) and a["resseq"] == lig_resseq

            site_altlocs = [
                a for a in altlocs
                if a["coord"] is not None and _in_receptor_chains(a) and not _is_ligand_atom(a)
                and _min_dist([a["coord"]], ligand_coords) <= SITE_RADIUS_A
            ]
            if site_altlocs:
                residues = sorted({f"{a['resname']} {a['chain']}{a['resseq']}{a['icode']}" for a in site_altlocs})
                findings.append(_finding(
                    "SITE_ALTLOC_PRESENT", "high",
                    f"An alternate-location conformer lies within {SITE_RADIUS_A:.0f} Å of the declared "
                    "ligand; the conformer-selection policy directly affects the local binding-site geometry.",
                    len(residues), residues,
                ))

            site_metals = [
                a for a in metals
                if a["coord"] is not None and _in_receptor_chains(a)
                and _min_dist([a["coord"]], ligand_coords) <= SITE_RADIUS_A
            ]
            if site_metals:
                residues = sorted({f"{a['element']} {a['resname']} {a['chain']}{a['resseq']}" for a in site_metals})
                findings.append(_finding(
                    "SITE_METAL_PRESENT", "medium",
                    f"A metal ion lies within {SITE_RADIUS_A:.0f} Å of the declared ligand; coordination "
                    "and charge treatment directly affect the local binding-site geometry.",
                    len(residues), residues,
                ))

    severities = Counter(finding["severity"] for finding in findings)
    return {
        "schema_version": "0.1",
        "tool": {"name": "DockPrep Audit", "version": "0.2.0"},
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(source), "format": "PDB", "lines": len(lines), "atom_records": len(atoms)},
        "summary": {
            "status": "review_required" if severities["high"] else "review_recommended" if findings else "no_rules_triggered",
            "findings": len(findings),
            "high": severities["high"], "medium": severities["medium"], "low": severities["low"],
        },
        "findings": findings,
    }
