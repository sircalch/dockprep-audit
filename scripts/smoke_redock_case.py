"""Reproducible technical smoke test: extract receptor/ligand, prepare PDBQT with
Meeko (no ProDy), redock with AutoDock Vina, and report heavy-atom RMSD to the
crystallographic pose.

This validates the PIPELINE, not a scientific preparation policy. It uses a
single technical default (drop template-mismatched residues via
--allow_bad_res, keep protein ATOM records only, assign ligand bond orders
from the RCSB chemical-component SMILES). Preparation policy comparisons
belong to Fase B, not to this script.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from urllib.error import URLError

import dimorphite_dl
from rdkit import Chem
from rdkit.Chem import AllChem

PHYSIOLOGICAL_PH = 7.4

REPO_ROOT = Path(__file__).resolve().parent.parent
VINA_BIN = REPO_ROOT / "tools" / "vina" / "vina_1.2.7_win.exe"
MK_PREPARE_RECEPTOR = REPO_ROOT / ".venv" / "Scripts" / "mk_prepare_receptor.exe"
MK_PREPARE_LIGAND = REPO_ROOT / ".venv" / "Scripts" / "mk_prepare_ligand.exe"


def fetch_ligand_smiles(component_id: str) -> str:
    """RCSB's data.rcsb.org endpoint had a sustained outage during this
    session (confirmed with direct curl, not just from this script) --
    retry with backoff rather than failing every case in a batch on a
    transient/prolonged external issue unrelated to the structures
    themselves.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{component_id}"
    last_exc = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.load(resp)
            break
        except (URLError, OSError) as exc:
            last_exc = exc
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))
    else:
        raise last_exc
    for entry in data["pdbx_chem_comp_descriptor"]:
        if entry["type"] == "SMILES":
            return entry["descriptor"]
    raise ValueError(f"No SMILES descriptor found for component {component_id}")


BACKBONE_ATOM_NAMES = {"N", "CA", "C"}
WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}
# Metal ion policy (pendiente #7, seccion 8), decided 2026-08-21: keep metal
# ions structurally (they are rarely-if-ever optional context for a
# metalloenzyme active site, unlike water, whose relevance is genuinely
# case-by-case) using Meeko's own native ion templates -- which already
# carry the correct formal charge (Zn2+, Mg2+, Ca2+, Mn2+, Fe3+) -- rather
# than inventing ad hoc parameters. This is a shared, structural-fidelity
# rule applied identically across all three Fase B preparation policies,
# the same reasoning as the altLoc policy: no defensible "preparation
# philosophy" argues for silently deleting a coordinating active-site ion.
# Confirmed present in meeko/data/residue_chem_templates.json; NA/K/CU/CO/
# NI/CD have no native template there and are not used by any metal in
# this cohort (all three metal_or_cofactor cases use ZN).
MEEKO_METAL_ION_RESNAMES = {"ZN", "MG", "CA", "MN", "FE"}


def find_polymer_hetatm_residues(raw_pdb_text: str, chains: set[str]) -> set[tuple[str, str]]:
    """Identify HETATM residues that are actually part of the polymer chain
    backbone (modified/non-standard amino acids such as ALY, MSE, SEP, ...)
    rather than a bound ligand, water, or metal.

    PDB marks any non-standard residue as HETATM regardless of whether it
    sits in the middle of the polypeptide. A residue with N/CA/C backbone
    atom names present is, for any practical purpose, a polymer residue --
    a small-molecule ligand does not use protein backbone atom naming.
    Dropping such a residue (as a plain 'keep only ATOM records' filter
    would) deletes real backbone atoms and creates an artificial chain gap,
    breaking inter-residue bond perception for its neighbors too (found in
    4RJ3: ALY A:33, 3.49 A from the ligand, sits between LEU A:32 and
    LYS A:34 in SEQRES with a complete backbone).
    """
    backbone_atoms_by_residue: dict[tuple[str, str], set[str]] = {}
    for line in raw_pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname in WATER_RESNAMES:
            continue
        chain = line[21:22]
        if chains and chain not in chains:
            continue
        resseq = line[22:27]  # includes insertion code
        key = (chain, resseq)
        name = line[12:16].strip()
        if name in BACKBONE_ATOM_NAMES:
            backbone_atoms_by_residue.setdefault(key, set()).add(name)
    return {key for key, names in backbone_atoms_by_residue.items() if BACKBONE_ATOM_NAMES <= names}


def find_meeko_supported_metal_residues(raw_pdb_text: str, chains: set[str]) -> set[tuple[str, str]]:
    """Identify metal-ion HETATM residues that Meeko has a native residue
    template for (see MEEKO_METAL_ION_RESNAMES). These are single-atom
    residues with no backbone, so find_polymer_hetatm_residues() does not
    (and should not) catch them; without this, they are silently dropped by
    an 'ATOM records only' filter exactly like ALY was before that fix.
    """
    keys = set()
    for line in raw_pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname not in MEEKO_METAL_ION_RESNAMES:
            continue
        chain = line[21:22]
        if chains and chain not in chains:
            continue
        keys.add((chain, line[22:27]))
    return keys


METAL_COORDINATION_DISTANCE_A = 2.6


def assign_metal_coordinating_histidine_tautomers(raw_pdb_text: str, chains: set[str]) -> tuple[dict[tuple[str, str], str], list[dict]]:
    """Receptor protonation policy, part 1 (pendiente #8, seccion 8), decided
    2026-08-21: a HIS whose ND1 or NE2 nitrogen directly coordinates a metal
    ion (<= METAL_COORDINATION_DISTANCE_A) has its tautomer FIXED by that
    geometry, not left to Meeko's blind default. Meeko has no crystal
    structure to consult when a plain 'HIS' resname is ambiguous between
    HID/HIE/HIP (all three share the same heavy atoms), so it always tries
    HIE first (confirmed in residue_chem_templates.json's 'ambiguous' list:
    HIS -> ['HIE', 'HID', 'HIP', ...], first match wins). But whichever
    nitrogen is engaged in metal coordination must be the deprotonated
    lone-pair donor: if ND1 coordinates, the residue must be HIE (proton on
    NE2); if NE2 coordinates, it must be HID (proton on ND1). Checked
    empirically across the 3 metal_or_cofactor cases: Meeko's HIE-first
    default would have been chemically wrong for 8 of 12 zinc-coordinating
    histidines (5A2S, 1CBX, 4EXS combined, both chains).

    This is a shared, structural-fidelity rule applied identically across
    all three Fase B preparation policies, same reasoning as altLoc and
    metal retention -- not a policy comparison axis.

    Returns (tautomer_by_residue, decisions_log) where tautomer_by_residue
    maps (chain, resseq) -> 'HIE' or 'HID' only for residues whose tautomer
    was fixed by metal coordination; all other HIS residues are left alone
    (Meeko's HIE-first default applies, undisturbed).
    """
    metal_coords = []
    his_nitrogens: dict[tuple[str, str], dict[str, tuple[float, float, float]]] = {}
    for line in raw_pdb_text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        chain = line[21:22]
        if chains and chain not in chains:
            continue
        resname = line[17:20].strip()
        coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if line.startswith("HETATM") and resname in MEEKO_METAL_ION_RESNAMES:
            metal_coords.append(coord)
        elif resname == "HIS":
            name = line[12:16].strip()
            if name in ("ND1", "NE2"):
                his_nitrogens.setdefault((chain, line[22:27]), {})[name] = coord

    tautomer_by_residue = {}
    decisions_log = []
    for key, nitrogens in his_nitrogens.items():
        coordinating = {
            name: min((math.dist(coord, m) for m in metal_coords), default=math.inf)
            for name, coord in nitrogens.items()
        }
        nd1_dist, ne2_dist = coordinating.get("ND1", math.inf), coordinating.get("NE2", math.inf)
        if nd1_dist > METAL_COORDINATION_DISTANCE_A and ne2_dist > METAL_COORDINATION_DISTANCE_A:
            continue  # not metal-coordinating; leave to Meeko's default
        tautomer = "HIE" if nd1_dist <= ne2_dist else "HID"
        tautomer_by_residue[key] = tautomer
        decisions_log.append({
            "chain": key[0], "resseq": key[1].strip(),
            "nd1_metal_dist_A": round(nd1_dist, 3) if nd1_dist != math.inf else None,
            "ne2_metal_dist_A": round(ne2_dist, 3) if ne2_dist != math.inf else None,
            "assigned_tautomer": tautomer,
        })
    return tautomer_by_residue, decisions_log


ALTLOC_POLICIES = ("highest_occupancy", "lowest_occupancy")


def choose_altloc_conformers(raw_pdb_text: str, chains: set[str],
                              policy: str = "highest_occupancy") -> tuple[dict[tuple[str, str], str], list[dict]]:
    """Pick one altLoc conformer per residue.

    ``policy="highest_occupancy"`` (the project's original, default policy,
    decided 2026-08-21, pendiente #5 seccion 8) keeps the conformer with the
    highest reported occupancy, breaking an exact tie to 'A' as a declared
    default rather than an arbitrary one. This is a structural-fidelity
    decision (which conformation most likely existed in the crystal), not a
    preparation-philosophy choice, so both policies below apply identically
    regardless of water-handling policy.

    ``policy="lowest_occupancy"`` (F2, added 2026-08-28 to test whether the
    audit-flagged altLoc risk factor predicts *when* the conformer-selection
    choice itself changes the redocking result, mirroring the water-policy
    manipulation of Section 2.3/Eq. 1) keeps the conformer with the lowest
    reported occupancy instead -- the symmetric opposite manipulation, with
    an exact tie broken to the last declared letter rather than 'A'.

    Returns (winner_by_residue, decisions_log) where winner_by_residue maps
    (chain, resseq) -> chosen altloc letter, and decisions_log records, per
    residue, the occupancies seen and whether the choice came from occupancy
    or from the tie-break, for transparency (never a silent choice).
    """
    if policy not in ALTLOC_POLICIES:
        raise ValueError(f"Unknown altloc policy: {policy!r} (expected one of {ALTLOC_POLICIES})")

    occupancy_by_residue: dict[tuple[str, str], dict[str, float]] = {}
    resname_by_residue: dict[tuple[str, str], str] = {}
    for line in raw_pdb_text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        altloc = line[16:17]
        if altloc == " ":
            continue
        chain = line[21:22]
        if chains and chain not in chains:
            continue
        resseq = line[22:27]
        key = (chain, resseq)
        occ = float(line[54:60])
        residue_occupancies = occupancy_by_residue.setdefault(key, {})
        residue_occupancies[altloc] = max(residue_occupancies.get(altloc, 0.0), occ)
        resname_by_residue[key] = line[17:20].strip()

    winner_by_residue = {}
    decisions_log = []
    for key, occupancies in occupancy_by_residue.items():
        if policy == "highest_occupancy":
            target_occ = max(occupancies.values())
            tied = sorted(letter for letter, occ in occupancies.items() if occ == target_occ)
            chosen = tied[0] if tied[0] == "A" or "A" not in tied else "A"
            basis = "occupancy" if len(tied) == 1 else ("tie_break_A" if chosen == "A" else "tie_break_first_letter")
        else:
            target_occ = min(occupancies.values())
            tied = sorted(letter for letter, occ in occupancies.items() if occ == target_occ)
            chosen = tied[-1]
            basis = "occupancy" if len(tied) == 1 else "tie_break_last_letter"
        winner_by_residue[key] = chosen
        decisions_log.append({
            "chain": key[0], "resseq": key[1].strip(), "resname": resname_by_residue[key],
            "occupancies": occupancies, "chosen_altloc": chosen, "basis": basis, "policy": policy,
        })
    return winner_by_residue, decisions_log


METAL_POLICIES = ("retain", "remove")


def extract_receptor_atoms(raw_pdb: Path, chains: set[str], out_path: Path,
                            altloc_log_path: Path | None = None,
                            his_tautomer_log_path: Path | None = None,
                            altloc_policy: str = "highest_occupancy",
                            metal_policy: str = "retain") -> None:
    """Keep protein backbone/sidechain records for the given chains,
    including HETATM-recorded modified residues that are polymer-continuous
    (see find_polymer_hetatm_residues) -- excluding them would be a bug, not
    a policy choice, since it deletes real backbone atoms. altLoc conformers
    are resolved by choose_altloc_conformers() (default: occupancy-based,
    tie-break A; see that function for the F2 lowest_occupancy alternative),
    and metal-coordinating HIS tautomers by
    assign_metal_coordinating_histidine_tautomers() -- both project-wide
    policies shared by all three Fase B preparation policies, decided
    2026-08-21.
    """
    if metal_policy not in METAL_POLICIES:
        raise ValueError(f"Unknown metal policy: {metal_policy!r} (expected one of {METAL_POLICIES})")

    raw_text = raw_pdb.read_text(encoding="utf-8", errors="replace")
    polymer_hetatm_residues = find_polymer_hetatm_residues(raw_text, chains)
    metal_ion_residues = find_meeko_supported_metal_residues(raw_text, chains) if metal_policy == "retain" else set()
    kept_hetatm_residues = polymer_hetatm_residues | metal_ion_residues
    winner_by_residue, decisions_log = choose_altloc_conformers(raw_text, chains, policy=altloc_policy)
    if altloc_log_path is not None:
        altloc_log_path.write_text(json.dumps(decisions_log, indent=2), encoding="utf-8")
    # The geometric HIS-tautomer-fixing rule (assign_metal_coordinating_histidine_tautomers)
    # only makes physical sense when the coordinating metal is actually present in the
    # receptor; under the "remove" metal policy there is nothing to coordinate, so HIS
    # tautomers fall back to Meeko's own default rather than being forced by a metal
    # that is no longer part of the model.
    if metal_policy == "retain":
        his_tautomer_by_residue, his_decisions_log = assign_metal_coordinating_histidine_tautomers(raw_text, chains)
    else:
        his_tautomer_by_residue, his_decisions_log = {}, []
    if his_tautomer_log_path is not None:
        his_tautomer_log_path.write_text(json.dumps(his_decisions_log, indent=2), encoding="utf-8")

    lines = []
    for line in raw_text.splitlines():
        is_atom = line.startswith("ATOM  ")
        is_kept_hetatm = line.startswith("HETATM") and (line[21:22], line[22:27]) in kept_hetatm_residues
        if not (is_atom or is_kept_hetatm):
            continue
        if chains and line[21:22] not in chains:
            continue
        altloc = line[16:17]
        if altloc != " ":
            winner = winner_by_residue.get((line[21:22], line[22:27]))
            if altloc != winner:
                continue
            line = line[:16] + " " + line[17:]
        if is_kept_hetatm:
            line = "ATOM  " + line[6:]
        if line[17:20] == "HIS":
            tautomer = his_tautomer_by_residue.get((line[21:22], line[22:27]))
            if tautomer is not None:
                line = line[:17] + tautomer + line[20:]
        lines.append(line)
    out_path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


WATER_BRIDGE_DISTANCE_A = 3.0


def find_bridging_waters(raw_pdb_text: str, ligand_atoms: list[tuple[float, float, float]],
                          chains: set[str]) -> list[dict]:
    """Water policy, Policia 1 ("referencia conservadora"), decided 2026-08-21
    (pendiente #6, seccion 8): a water is kept only if it plausibly bridges
    the ligand and the receptor -- its oxygen within WATER_BRIDGE_DISTANCE_A
    of BOTH a ligand atom and a receptor atom. This is the "loose water"
    definition used in docking studies (e.g. GOLD/Verdonk-style water
    filtering), not the looser "any water within 4 A of the ligand" screen
    the audit engine already uses to flag candidates for human review --
    that screen is deliberately permissive (for surfacing), this one is
    deliberately strict (for deciding what to keep in a receptor model).

    Policies 2 (simplified conventional) and 3 (open-tool default) both
    drop all waters unconditionally -- already the pipeline's existing
    behavior, requiring no extra code.
    """
    protein_atoms = []
    water_atoms_by_residue: dict[tuple[str, str], list[tuple[float, float, float]]] = {}
    water_meta_by_residue: dict[tuple[str, str], dict] = {}
    for line in raw_pdb_text.splitlines():
        if line.startswith("ATOM  "):
            if chains and line[21:22] not in chains:
                continue
            protein_atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif line.startswith("HETATM") and line[17:20].strip() in WATER_RESNAMES:
            chain = line[21:22]
            if chains and chain not in chains:
                continue
            key = (chain, line[22:27])
            coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            water_atoms_by_residue.setdefault(key, []).append(coord)
            water_meta_by_residue[key] = {"chain": chain, "resseq": line[22:26].strip()}

    bridging = []
    for key, coords in water_atoms_by_residue.items():
        ligand_hits = [(math.dist(w, l), l) for w in coords for l in ligand_atoms
                       if math.dist(w, l) <= WATER_BRIDGE_DISTANCE_A]
        if not ligand_hits:
            continue
        protein_hits = [(math.dist(w, p), p) for w in coords for p in protein_atoms
                         if math.dist(w, p) <= WATER_BRIDGE_DISTANCE_A]
        if not protein_hits:
            continue
        meta = dict(water_meta_by_residue[key])
        oxygen_coord = coords[0]
        meta["oxygen_coord"] = oxygen_coord
        # F9 (2026-09-01): record the specific ligand/protein atoms that
        # qualified this water as bridging -- its own nearest contacts --
        # so a richer water representation can orient explicit hydrogens
        # toward them, instead of leaving the water an undirected point charge.
        # Anchored specifically to oxygen_coord (not just any conformer in
        # coords): a disordered water can have multiple HETATM records for
        # the same residue, and ligand_hits/protein_hits above deliberately
        # pool contacts across all of them to decide bridging status (that
        # qualification logic is unchanged and must stay that way -- it is
        # what Sections 3.5-3.10 and F1-F8 were computed against). But
        # picking the globally nearest contact across every conformer, while
        # placing only oxygen_coord's specific conformer, could orient a
        # hydrogen toward an atom that has nothing to do with this water's
        # actual position. Restricting to hits whose water atom is
        # oxygen_coord itself keeps the two consistent; falling back to the
        # global nearest only if oxygen_coord itself has no in-range contact
        # of that kind (shouldn't happen given the bridging check above, but
        # keeps this defined rather than raising on an unexpected input).
        own_ligand_hits = [(math.dist(oxygen_coord, l), l) for l in ligand_atoms
                            if math.dist(oxygen_coord, l) <= WATER_BRIDGE_DISTANCE_A]
        own_protein_hits = [(math.dist(oxygen_coord, p), p) for p in protein_atoms
                             if math.dist(oxygen_coord, p) <= WATER_BRIDGE_DISTANCE_A]
        meta["nearest_ligand_atom"] = min(own_ligand_hits or ligand_hits, key=lambda t: t[0])[1]
        meta["nearest_protein_atom"] = min(own_protein_hits or protein_hits, key=lambda t: t[0])[1]
        bridging.append(meta)
    return bridging


def _splice_static_atoms_into_pdbqt(receptor_pdbqt: Path, build_lines) -> None:
    """Shared by append_bridging_waters_to_pdbqt and
    append_oriented_bridging_waters_to_pdbqt: read the prepared receptor,
    hand `build_lines` the next free atom serial number to generate new
    HETATM lines from, and splice them in before a trailing 'TORSDOF 0' line
    (or at the end if there isn't one). Kept as one function shared by both
    water policies so a future fix to the splice point only has to be made
    once.
    """
    lines = receptor_pdbqt.read_text(encoding="utf-8").splitlines()
    next_serial = sum(1 for l in lines if l.startswith(("ATOM", "HETATM"))) + 1
    new_lines = build_lines(next_serial)
    if lines and lines[-1].strip() == "TORSDOF 0":
        lines = lines[:-1] + new_lines + [lines[-1]]
    else:
        lines = lines + new_lines
    receptor_pdbqt.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_bridging_waters_to_pdbqt(receptor_pdbqt: Path, bridging_waters: list[dict]) -> None:
    """Append kept bridging waters to the prepared receptor PDBQT as static,
    non-bonded oxygen atoms (AutoDock atom type OA, H-bond acceptor), using
    the TIP3P water-oxygen partial charge (-0.834 e) since no other
    receptor-wide charge model is being applied to them. Meeko has no
    residue template for HOH (waters aren't a polymer residue), so they
    cannot go through the normal Polymer/template pipeline -- they are
    appended directly, which is standard practice for treating a handful of
    ordered waters as part of a rigid AutoDock/Vina receptor.
    """
    def build(serial: int) -> list[str]:
        new_lines = []
        for w in bridging_waters:
            x, y, z = w["oxygen_coord"]
            new_lines.append(
                f"HETATM{serial:5d}  O   HOH {w['chain']}{w['resseq']:>4s}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    -0.834 OA"
            )
            serial += 1
        return new_lines

    _splice_static_atoms_into_pdbqt(receptor_pdbqt, build)


TIP3P_OH_LENGTH_A = 0.9572
TIP3P_HOH_ANGLE_DEG = 104.52


def _oriented_hydrogen_positions(oxygen: tuple[float, float, float],
                                  target_a: tuple[float, float, float],
                                  target_b: tuple[float, float, float]
                                  ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Place two O-H vectors at the idealized TIP3P bond length and H-O-H
    angle, oriented in the plane spanned by the water's two bridging
    contacts (nearest_ligand_atom, nearest_protein_atom) -- the same two
    atoms that already qualified this water as "bridging" in
    find_bridging_waters. This is a geometric best-guess, not a resolved
    crystallographic position (waters are essentially never resolved with
    hydrogens at the diffraction resolution typical of these structures);
    it replaces the prior undirected point-oxygen model with a directional
    dipole aimed at the same contacts the water-policy definition itself
    already privileges, rather than introducing a new, separately-tuned
    criterion for orientation. F9, PROJECT-ROADMAP.md section 19.
    """
    def sub(u, v):
        return (u[0] - v[0], u[1] - v[1], u[2] - v[2])

    def norm(v):
        n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
        return (v[0] / n, v[1] / n, v[2] / n) if n > 1e-9 else (0.0, 0.0, 1.0)

    def dot(u, v):
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    def add(u, v):
        return (u[0] + v[0], u[1] + v[1], u[2] + v[2])

    def scale(u, s):
        return (u[0] * s, u[1] * s, u[2] * s)

    d1 = norm(sub(target_a, oxygen))
    d2 = norm(sub(target_b, oxygen))
    bisector = norm(add(d1, d2))
    n_raw = sub(d1, d2)
    if math.dist(n_raw, (0.0, 0.0, 0.0)) < 1e-6:
        # target directions coincide (degenerate); pick an arbitrary
        # perpendicular so the two hydrogens don't collapse onto each other.
        arbitrary = (1.0, 0.0, 0.0) if abs(d1[0]) < 0.9 else (0.0, 1.0, 0.0)
        n_perp = norm(sub(arbitrary, scale(bisector, dot(arbitrary, bisector))))
    else:
        n_raw = norm(n_raw)
        n_perp = norm(sub(n_raw, scale(bisector, dot(n_raw, bisector))))

    half_angle = math.radians(TIP3P_HOH_ANGLE_DEG / 2.0)
    h1_dir = add(scale(bisector, math.cos(half_angle)), scale(n_perp, math.sin(half_angle)))
    h2_dir = add(scale(bisector, math.cos(half_angle)), scale(n_perp, -math.sin(half_angle)))
    h1 = add(oxygen, scale(norm(h1_dir), TIP3P_OH_LENGTH_A))
    h2 = add(oxygen, scale(norm(h2_dir), TIP3P_OH_LENGTH_A))
    return h1, h2


def append_oriented_bridging_waters_to_pdbqt(receptor_pdbqt: Path, bridging_waters: list[dict]) -> None:
    """F9 richer water representation: same rigid, static-atom treatment as
    append_bridging_waters_to_pdbqt, but each water is a 3-atom O+H+H rigid
    body (idealized TIP3P geometry, oriented per _oriented_hydrogen_positions)
    carrying the full TIP3P charge set (O=-0.834e, H=+0.417e each) instead of
    a bare oxygen at the oxygen's full -0.834e. Hydrogens use AutoDock atom
    type HD (polar/donor hydrogen).
    """
    def build(serial: int) -> list[str]:
        new_lines = []
        for w in bridging_waters:
            ox = w["oxygen_coord"]
            h1, h2 = _oriented_hydrogen_positions(ox, w["nearest_ligand_atom"], w["nearest_protein_atom"])
            chain, resseq = w["chain"], w["resseq"]
            new_lines.append(
                f"HETATM{serial:5d}  O   HOH {chain}{resseq:>4s}    "
                f"{ox[0]:8.3f}{ox[1]:8.3f}{ox[2]:8.3f}  1.00  0.00    -0.834 OA"
            )
            serial += 1
            for h_name, h_coord in (("H1", h1), ("H2", h2)):
                new_lines.append(
                    f"HETATM{serial:5d}  {h_name:<4s}HOH {chain}{resseq:>4s}    "
                    f"{h_coord[0]:8.3f}{h_coord[1]:8.3f}{h_coord[2]:8.3f}  1.00  0.00    +0.417 HD"
                )
                serial += 1
        return new_lines

    _splice_static_atoms_into_pdbqt(receptor_pdbqt, build)


def extract_ligand_atoms(raw_pdb: Path, component_id: str, chain: str, resseq: str, out_path: Path) -> None:
    """Extract the declared ligand instance, resolving altLoc with the same
    policy as the receptor (pendiente #5, seccion 8: highest occupancy,
    tie-break A) -- found missing 2026-08-22 when `3P0M` and `1SN5` (Fase E
    expansion cases) failed at the RMSD atom-matching step. Their ligand
    itself has altLoc conformers (e.g. `1SN5`'s T3 has three: A/B/C at
    0.40/0.40/0.20 occupancy) and this function was including all of them
    undifferentiated, giving RDKit/Meeko a chemically invalid double- or
    triple-exposure of the same atoms instead of one real conformer. This
    was never caught by the original 12-case pilot because none of those
    ligands happened to have altLoc.
    """
    matching = [
        line for line in raw_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("HETATM") and line[17:20].strip() == component_id
        and (not chain or line[21:22] == chain) and (not resseq or line[22:26].strip() == resseq)
    ]
    if not matching:
        raise ValueError(f"No HETATM lines matched component={component_id} chain={chain} resseq={resseq}")

    occupancy_by_altloc: dict[str, float] = {}
    for line in matching:
        altloc = line[16:17]
        if altloc == " ":
            continue
        occ = float(line[54:60])
        occupancy_by_altloc[altloc] = max(occupancy_by_altloc.get(altloc, 0.0), occ)

    lines = []
    if occupancy_by_altloc:
        best_occ = max(occupancy_by_altloc.values())
        tied = sorted(l for l, occ in occupancy_by_altloc.items() if occ == best_occ)
        winner = "A" if "A" in tied else tied[0]
        for line in matching:
            altloc = line[16:17]
            if altloc != " " and altloc != winner:
                continue
            if altloc == winner:
                line = line[:16] + " " + line[17:]
            lines.append(line)
    else:
        lines = matching

    out_path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def protonate_ligand_smiles(smiles: str) -> tuple[str, dict]:
    """Ligand protonation policy (pendiente #9, seccion 8), decided
    2026-08-21: assign the dominant ionization state at physiological pH
    (7.4) with dimorphite-dl, rather than docking whatever neutral state
    the RCSB chemical-component SMILES happens to encode.

    This matters concretely: BZS (1CBX's ligand, L-benzylsuccinic acid) is
    deposited as the neutral diacid (formal charge 0), but at pH 7.4
    (pKa of a carboxylic acid ~4-5) both carboxyls are expected to be
    ionized (net charge -2) -- docking the neutral form would be docking a
    species that is not the dominant one in solution, right where it
    matters most (a succinate-mimetic binding a catalytic zinc).

    The originally recommended tool for this step (scrubber, per the
    official Vina zinc-docking guide) would not install on Python 3.12
    (ModuleNotFoundError: 'imp', a module removed in 3.12; its setup.py
    predates that removal). dimorphite-dl is the underlying, actively
    maintained, lighter-weight protonation engine and installs cleanly.

    Returns (chosen_smiles, decision_log). If dimorphite-dl proposes more
    than one variant at this pH (a real possibility when two ionizable
    groups have close pKas), all are logged and the first is chosen as a
    declared, documented tie-break -- never a silent one.
    """
    variants = dimorphite_dl.protonate_smiles(smiles, ph_min=PHYSIOLOGICAL_PH, ph_max=PHYSIOLOGICAL_PH,
                                               precision=0.5)
    if not variants:
        raise ValueError(f"dimorphite-dl returned no protonation state for SMILES: {smiles}")
    log = {
        "input_smiles": smiles, "ph": PHYSIOLOGICAL_PH, "variants": variants,
        "basis": "only_variant" if len(variants) == 1 else "multiple_variants_tried_in_order",
    }
    return variants, log


def ligand_pdb_to_sdf(ligand_pdb: Path, component_id: str, sdf_path: Path,
                       protonation_log_path: Path | None = None) -> None:
    """Build the docking-ready ligand SDF from the physiologically protonated
    template. When dimorphite-dl returns more than one candidate ionization
    state (real for multi-site ligands whose pKas straddle 7.4, e.g. 1OHR's
    nelfinavir-like inhibitor with both a phenol and a tertiary amine), each
    is tried in order against the raw crystallographic heavy-atom
    connectivity; the first that successfully assigns bond orders and
    sanitizes is kept. This is not picking whichever answer is convenient --
    it is a documented fallback for candidates that fail to even build as a
    valid molecule (confirmed for 1OHR: dimorphite-dl's own first variant hit
    an internal RDKit valence error), not a re-ranking of chemically valid
    candidates by their docking outcome.
    """
    deposited_smiles = fetch_ligand_smiles(component_id)
    variants, protonation_log = protonate_ligand_smiles(deposited_smiles)
    raw = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False, sanitize=True)
    if raw is None:
        raise ValueError(f"RDKit could not parse ligand PDB {ligand_pdb}")

    mol = None
    attempts = []
    for variant_smiles in variants:
        try:
            template = Chem.MolFromSmiles(variant_smiles)
            if template is None:
                attempts.append({"smiles": variant_smiles, "ok": False, "error": "RDKit could not parse SMILES"})
                continue
            candidate = AllChem.AssignBondOrdersFromTemplate(template, raw)
            candidate = Chem.AddHs(candidate, addCoords=True)
            Chem.SanitizeMol(candidate)
            attempts.append({"smiles": variant_smiles, "ok": True, "error": None})
            mol = candidate
            break
        except Exception as exc:
            attempts.append({"smiles": variant_smiles, "ok": False, "error": str(exc)})
    protonation_log["attempts"] = attempts
    protonation_log["chosen_smiles"] = attempts[-1]["smiles"] if mol is not None else None
    if protonation_log_path is not None:
        protonation_log_path.write_text(json.dumps(protonation_log, indent=2), encoding="utf-8")
    if mol is None:
        raise ValueError(f"No protonation variant for {component_id} could be built: {attempts}")

    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()


RG_BOX_MULTIPLIER = 2.9


def compute_ligand_box(ligand_pdb: Path, padding: float | None = None) -> dict:
    """Docking box policy (pendiente #10, seccion 8), decided 2026-08-21: a
    cubic box of size 2.9 x the ligand's radius of gyration (Rg), centered
    on its geometric center -- the box size empirically found to maximize
    redocking/pocket-recovery accuracy (Feinstein & Brylinski 2015, PMC4468813:
    RMSD 4.9->4.0 A, binding-residue recovery 0.78->0.92 when moving from an
    ad hoc box to this rule). Replaces the earlier 'ligand bounding box +
    padding' heuristic used during pipeline validation (padding kept as an
    optional override, e.g. for reproducing old runs, but is not the policy).
    """
    coords = []
    for line in ligand_pdb.read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM"):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    if padding is not None:
        xs, ys, zs = zip(*coords)
        return {
            "center_x": (min(xs) + max(xs)) / 2, "center_y": (min(ys) + max(ys)) / 2,
            "center_z": (min(zs) + max(zs)) / 2,
            "size_x": max(xs) - min(xs) + padding, "size_y": max(ys) - min(ys) + padding,
            "size_z": max(zs) - min(zs) + padding,
        }
    rg = math.sqrt(sum((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2 for c in coords) / n)
    size = RG_BOX_MULTIPLIER * rg
    return {
        "center_x": cx, "center_y": cy, "center_z": cz,
        "size_x": size, "size_y": size, "size_z": size,
        "radius_of_gyration_A": round(rg, 3),
    }


VINA_BOX_KEYS = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")


def write_box_file(box: dict, path: Path) -> None:
    path.write_text("\n".join(f"{k} = {box[k]:.3f}" for k in VINA_BOX_KEYS) + "\n", encoding="utf-8")


def _wait_for_file(path: Path, timeout_s: float = 5.0) -> None:
    """Guard against a transient file-visibility race observed on this
    OneDrive-synced checkout: a just-written file is occasionally not yet
    visible to a subprocess spawned immediately afterward (mk_prepare_receptor
    raising FileNotFoundError on a file confirmed present moments later).
    Purely a wait -- no science-affecting behavior.
    """
    deadline = time.time() + timeout_s
    while not path.exists() and time.time() < deadline:
        time.sleep(0.1)


def run(cmd: list[str], cwd: Path) -> str:
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")
    return result.stdout + result.stderr


def parse_ignored_residues(mk_prepare_receptor_output: str) -> list[str]:
    """Extract the residue keys Meeko silently dropped via --allow_bad_res so
    they can be recorded instead of disappearing into console output. This is
    a transparency requirement, not a scientific decision about which
    residues SHOULD be dropped -- see roadmap section 8, pendiente #1/#7.
    """
    match = re.search(r"Template matching failed for: (\[[^\]]*\])", mk_prepare_receptor_output)
    if not match:
        return []
    return ast.literal_eval(match.group(1))


def read_heavy_atom_coords_pdb(path: Path) -> list[tuple[float, float, float]]:
    coords = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM"):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


def read_all_atom_coords_pdbqt(path: Path, model: int | None = 1) -> list[tuple[float, float, float]]:
    """Read every atom's coordinates (heavy and hydrogen) in file order.

    ligand.pdbqt has no MODEL/ENDMDL wrapper (pass model=None); a docked
    output file has one MODEL block per pose.
    """
    coords = []
    in_model = model is None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            in_model = int(line.split()[1]) == model
        elif line.startswith("ENDMDL"):
            if in_model and model is not None:
                break
        elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


def build_reference_to_pdbqt_index_map(ref_coords: list[tuple[float, float, float]],
                                        ligand_pdbqt_coords: list[tuple[float, float, float]],
                                        tolerance: float = 0.01) -> list[int]:
    """Map each reference (crystal) heavy atom to its position in ligand.pdbqt
    by exact coordinate match.

    Meeko's PDBQT torsion-tree writer reorders atoms relative to the input
    molecule (confirmed by inspection: this reordering is a no-op for a
    rigid ligand like progesterone but scrambles order for a flexible one
    like retinol), so file position cannot be trusted, and the 'REMARK
    SMILES IDX' hint turned out to reference the SMILES string's own atom
    order rather than the input order, giving wrong correspondences too.
    Coordinates are the one thing Meeko does not touch: mk_prepare_ligand
    keeps the exact input 3-D coordinates for heavy atoms, so matching by
    position (not identity/order) is exact and verified to 0.000 A here.
    """
    used: set[int] = set()
    index_map = []
    for ref_atom in ref_coords:
        best_j, best_d = None, tolerance
        for j, cand in enumerate(ligand_pdbqt_coords):
            if j in used:
                continue
            d = math.dist(ref_atom, cand)
            if d <= best_d:
                best_d, best_j = d, j
        if best_j is None:
            raise ValueError(f"No ligand.pdbqt atom within {tolerance} A of reference atom {ref_atom}")
        used.add(best_j)
        index_map.append(best_j)
    return index_map


def rmsd(a: list[tuple[float, float, float]], b: list[tuple[float, float, float]]) -> float:
    if len(a) != len(b):
        raise ValueError(f"Atom count mismatch: reference={len(a)} pose={len(b)}")
    sq = sum((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 for p, q in zip(a, b))
    return math.sqrt(sq / len(a))


def count_models(pdbqt_path: Path) -> int:
    model_ids = [int(line.split()[1]) for line in pdbqt_path.read_text(encoding="utf-8").splitlines()
                 if line.startswith("MODEL")]
    return max(model_ids, default=1)


def classify_redocking(rmsds_by_rank: list[float], success_threshold: float = 2.0) -> str:
    """Three-way outcome used in redocking benchmarks (e.g. the Astex Diverse
    set evaluations): a scoring function can fail to RANK a near-native pose
    first even when the search DID sample one. Distinguishing this from a
    true sampling failure matters for this project because the two point to
    different causes -- one is a scoring-function limitation independent of
    receptor prep, the other may reflect a genuinely broken binding site
    (e.g. a removed structural water/metal the ligand needed).

    - success:      rank-1 pose is within the threshold.
    - scoring_fail:  some pose is within the threshold, but not rank 1.
    - sampling_fail: no pose sampled is within the threshold.
    """
    if rmsds_by_rank[0] <= success_threshold:
        return "success"
    if any(r <= success_threshold for r in rmsds_by_rank):
        return "scoring_fail"
    return "sampling_fail"


def run_case(pdb_id: str, component_id: str, chain: str, resseq: str, receptor_chains: str,
             raw_dir: Path, out_dir: Path, padding: float | None, exhaustiveness: int, seed: int,
             water_policy: str = "none", altloc_policy: str = "highest_occupancy",
             metal_policy: str = "retain") -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_pdb = raw_dir / f"{pdb_id}.pdb"

    receptor_raw = out_dir / "receptor_raw.pdb"
    chains = {c for c in receptor_chains.split(";") if c} if receptor_chains else set()
    extract_receptor_atoms(raw_pdb, chains, receptor_raw, altloc_log_path=out_dir / "altloc_decisions.json",
                           his_tautomer_log_path=out_dir / "his_tautomer_decisions.json",
                           altloc_policy=altloc_policy, metal_policy=metal_policy)
    _wait_for_file(receptor_raw)

    ligand_raw = out_dir / "ligand_raw.pdb"
    extract_ligand_atoms(raw_pdb, component_id, chain, resseq, ligand_raw)
    _wait_for_file(ligand_raw)

    receptor_pdbqt = out_dir / "receptor.pdbqt"
    receptor_prep_log = run([
        str(MK_PREPARE_RECEPTOR), "--read_pdb", str(receptor_raw), "-o", str(out_dir / "receptor"),
        "-p", "--allow_bad_res",
    ], cwd=out_dir)
    (out_dir / "receptor_prep.log").write_text(receptor_prep_log, encoding="utf-8")
    ignored_residues = parse_ignored_residues(receptor_prep_log)

    bridging_waters = []
    if water_policy in ("conservative", "conservative_oriented"):
        raw_text = raw_pdb.read_text(encoding="utf-8", errors="replace")
        ligand_coords_for_water = read_heavy_atom_coords_pdb(ligand_raw)
        bridging_waters = find_bridging_waters(raw_text, ligand_coords_for_water, chains)
        if bridging_waters:
            if water_policy == "conservative_oriented":
                append_oriented_bridging_waters_to_pdbqt(receptor_pdbqt, bridging_waters)
            else:
                append_bridging_waters_to_pdbqt(receptor_pdbqt, bridging_waters)
    elif water_policy != "none":
        raise ValueError(
            f"Unknown water_policy: {water_policy!r} "
            "(expected 'none', 'conservative', or 'conservative_oriented')")
    (out_dir / "bridging_waters.json").write_text(
        json.dumps([{**w, "oxygen_coord": list(w["oxygen_coord"])} for w in bridging_waters], indent=2),
        encoding="utf-8",
    )

    ligand_sdf = out_dir / "ligand.sdf"
    ligand_pdb_to_sdf(ligand_raw, component_id, ligand_sdf,
                      protonation_log_path=out_dir / "ligand_protonation.json")

    ligand_pdbqt = out_dir / "ligand.pdbqt"
    run([str(MK_PREPARE_LIGAND), "-i", str(ligand_sdf), "-o", str(ligand_pdbqt)], cwd=out_dir)

    box = compute_ligand_box(ligand_raw, padding)
    box_file = out_dir / "ligand_box.txt"
    write_box_file(box, box_file)

    redock_out = out_dir / "redock_out.pdbqt"
    vina_stdout = out_dir / "vina.log"
    result = subprocess.run(
        [str(VINA_BIN), "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
         "--config", str(box_file), "--exhaustiveness", str(exhaustiveness), "--seed", str(seed),
         "--out", str(redock_out)],
        cwd=out_dir, capture_output=True, text=True,
    )
    vina_stdout.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0 or not redock_out.exists():
        raise RuntimeError(f"Vina failed for {pdb_id}: {result.stderr[-2000:]}")

    ref_coords = read_heavy_atom_coords_pdb(ligand_raw)
    ligand_pdbqt_coords = read_all_atom_coords_pdbqt(ligand_pdbqt, model=None)
    index_map = build_reference_to_pdbqt_index_map(ref_coords, ligand_pdbqt_coords)

    n_poses = count_models(redock_out)
    rmsds_by_rank = []
    for rank in range(1, n_poses + 1):
        pose_all_atoms = read_all_atom_coords_pdbqt(redock_out, model=rank)
        pose_coords = [pose_all_atoms[j] for j in index_map]
        rmsds_by_rank.append(round(rmsd(ref_coords, pose_coords), 3))

    best_rank = min(range(len(rmsds_by_rank)), key=lambda i: rmsds_by_rank[i]) + 1

    return {
        "pdb_id": pdb_id,
        "ligand_heavy_atoms": len(ref_coords),
        "n_poses": n_poses,
        "rmsd_by_rank_A": rmsds_by_rank,
        "top_pose_rmsd_A": rmsds_by_rank[0],
        "best_pose_rmsd_A": rmsds_by_rank[best_rank - 1],
        "best_pose_rank": best_rank,
        "outcome": classify_redocking(rmsds_by_rank),
        "receptor_residues_ignored": ignored_residues,
        "water_policy": water_policy,
        "bridging_waters_kept": len(bridging_waters),
        "altloc_policy": altloc_policy,
        "metal_policy": metal_policy,
        "box": box,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb_id", help="PDB id from benchmark/pilot_cases.csv, e.g. 1A28")
    parser.add_argument("--cases", type=Path, default=REPO_ROOT / "benchmark" / "pilot_cases.csv")
    parser.add_argument("--eligibility", type=Path,
                         default=REPO_ROOT / "benchmark" / "pilot-eligibility" / "eligibility_results.csv")
    parser.add_argument("--raw-dir", type=Path, default=REPO_ROOT / "benchmark" / "pilot-inventory" / "raw-pdb")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "benchmark" / "pilot-extraction")
    parser.add_argument("--padding", type=float, default=None,
                         help="Angstrom padding around ligand bbox (legacy heuristic, e.g. --padding 20). "
                              "Default (unset) uses the decided box policy: cubic, 2.9x radius of gyration.")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--water-policy", choices=["none", "conservative", "conservative_oriented"],
                         default="none",
                         help="'none' = drop all waters (Fase B policies 2/3, current default). "
                              "'conservative' = keep bridging waters only (Fase B policy 1, "
                              "decided 2026-08-21): water <=3.0 A of both ligand and receptor, "
                              "modeled as a bare oxygen point charge. 'conservative_oriented' = same "
                              "bridging-water selection, but each water is a rigid O+H+H TIP3P body "
                              "with hydrogens oriented toward its own bridging contacts (F9, "
                              "added 2026-09-01).")
    parser.add_argument("--altloc-policy", choices=list(ALTLOC_POLICIES), default="highest_occupancy",
                         help="'highest_occupancy' = keep the most-occupied conformer, tie-break 'A' "
                              "(project default, decided 2026-08-21). 'lowest_occupancy' = keep the "
                              "least-occupied conformer, tie-break last letter (F2, added 2026-08-28, "
                              "the symmetric-opposite manipulation used to test altLoc-policy sensitivity).")
    parser.add_argument("--metal-policy", choices=list(METAL_POLICIES), default="retain",
                         help="'retain' = keep metal ions with geometrically assigned coordinating-HIS "
                              "tautomers (project default). 'remove' = drop metal ions entirely and let "
                              "HIS tautomers fall back to Meeko's default (F3, added 2026-08-28, the "
                              "manipulation used to test metal-policy sensitivity).")
    args = parser.parse_args()

    with args.cases.open(newline="", encoding="utf-8") as handle:
        cases = {row["pdb_id"].upper(): row for row in csv.DictReader(handle)}
    with args.eligibility.open(newline="", encoding="utf-8") as handle:
        eligibility = {row["pdb_id"].upper(): row for row in csv.DictReader(handle)}

    pdb_id = args.pdb_id.upper()
    case = cases[pdb_id]
    elig = eligibility[pdb_id]

    result = run_case(
        pdb_id=pdb_id,
        component_id=case["ligand_component_id"],
        chain=case.get("ligand_chain", "").strip(),
        resseq=case.get("ligand_resseq", "").strip(),
        receptor_chains=elig["receptor_chains"],
        raw_dir=args.raw_dir,
        out_dir=args.out_dir / pdb_id,
        padding=args.padding,
        exhaustiveness=args.exhaustiveness,
        seed=args.seed,
        water_policy=args.water_policy,
        altloc_policy=args.altloc_policy,
        metal_policy=args.metal_policy,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
