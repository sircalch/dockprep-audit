"""F6: verify the docking-box formula (PROJECT-ROADMAP.md section 19, item
F6), added 2026-08-28.

Confirms L = 2.9 x Rg (Section 2.3, Eq. 2) faithfully reproduces Feinstein &
Brylinski 2015 (J Cheminform 7:18, PMC4468813) -- their paper states
verbatim "the dimensions of the search space are 2.9 times larger than the
radius of gyration", validated on ligands of 6-100 heavy atoms with a cubic
box, with no explicit caveat for elongated or very small ligands (checked
via live fetch 2026-08-28) -- and empirically checks, across all 60 frozen
cases, whether that box actually contains the deposited ligand's
axis-aligned extent, and whether elongation (a shape the Rg-based formula
does not explicitly account for) predicts an inadequate box.
"""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CASES = list(csv.DictReader(open(REPO_ROOT / "benchmark" / "pilot_cases.csv", encoding="utf-8")))
RUN_ROOTS = [
    REPO_ROOT / "benchmark" / "phase-c-60" / "conservative_water" / "seed42",
    REPO_ROOT / "benchmark" / "phase-c" / "conservative_water" / "seed42",
]
RAW_RUN_DIRS = [
    REPO_ROOT / "benchmark" / "phase-c-60" / "raw_runs",
    REPO_ROOT / "benchmark" / "phase-c" / "raw_runs",
]


def find_case_dir(pdb_id: str) -> Path | None:
    for root in RUN_ROOTS:
        candidate = root / pdb_id
        if (candidate / "ligand_raw.pdb").exists():
            return candidate
    return None


def ligand_extents(pdb_id: str) -> tuple[float, float, float] | None:
    case_dir = find_case_dir(pdb_id)
    if case_dir is None:
        return None
    coords = []
    for line in (case_dir / "ligand_raw.pdb").read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM"):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    xs, ys, zs = zip(*coords)
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def box_size_and_rg(pdb_id: str) -> tuple[float, float, int] | None:
    for raw_dir in RAW_RUN_DIRS:
        raw_path = raw_dir / f"{pdb_id}_conservative_water_seed42.json"
        if raw_path.exists():
            result = json.loads(raw_path.read_text(encoding="utf-8"))["result"]
            return result["box"]["size_x"], result["box"]["radius_of_gyration_A"], result["ligand_heavy_atoms"]
    return None


def main() -> None:
    rows = []
    for row in CASES:
        pdb_id = row["pdb_id"]
        ext = ligand_extents(pdb_id)
        box_info = box_size_and_rg(pdb_id)
        if ext is None or box_info is None:
            print(f"  SKIP {pdb_id}: missing data")
            continue
        box, rg, n_heavy = box_info
        max_ext = max(ext)
        min_ext = min(ext)
        elongation = max_ext / min_ext if min_ext > 0 else float("inf")
        margin = box - max_ext
        rows.append({
            "pdb_id": pdb_id, "stratum": row["stratum"], "n_heavy": n_heavy, "rg": rg, "box": box,
            "ext_x": ext[0], "ext_y": ext[1], "ext_z": ext[2], "max_ext": max_ext,
            "elongation": elongation, "margin": margin, "fits": margin >= 0,
        })

    rows.sort(key=lambda r: r["margin"])

    print(f"=== Box adequacy, all {len(rows)} cases (sorted by tightest margin first) ===")
    print(f"{'case':8s} {'n_heavy':>7s} {'Rg':>6s} {'box':>6s} {'max_ext':>7s} {'margin':>7s} {'elong':>6s} {'fits':>5s}")
    for r in rows:
        print(f"{r['pdb_id']:8s} {r['n_heavy']:7d} {r['rg']:6.2f} {r['box']:6.2f} {r['max_ext']:7.2f} "
              f"{r['margin']:7.2f} {r['elongation']:6.2f} {str(r['fits']):>5s}")

    n_fail = sum(1 for r in rows if not r["fits"])
    print(f"\n{n_fail}/{len(rows)} cases have ligand extent exceeding the box "
          f"({n_fail/len(rows)*100:.1f}%)")

    print("\n=== Correlation: elongation vs. margin (tight boxes concentrated among elongated ligands?) ===")
    elong_sorted = sorted(rows, key=lambda r: -r["elongation"])
    print("Top 10 most elongated ligands:")
    for r in elong_sorted[:10]:
        print(f"  {r['pdb_id']}: elongation={r['elongation']:.2f}  margin={r['margin']:.2f}  fits={r['fits']}")

    print("\n=== Small-ligand check: n_heavy < 15 ===")
    small = [r for r in rows if r["n_heavy"] < 15]
    for r in sorted(small, key=lambda r: r["n_heavy"]):
        print(f"  {r['pdb_id']}: n_heavy={r['n_heavy']}  Rg={r['rg']:.2f}  box={r['box']:.2f}  "
              f"margin={r['margin']:.2f}  elongation={r['elongation']:.2f}")

    print(f"\nMean margin: {statistics.mean(r['margin'] for r in rows):.2f} A, "
          f"median: {statistics.median(r['margin'] for r in rows):.2f} A, "
          f"min: {min(r['margin'] for r in rows):.2f} A")


if __name__ == "__main__":
    main()
