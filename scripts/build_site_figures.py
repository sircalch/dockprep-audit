"""Roadmap fig 7: 3D binding-site geometry for one example case per finding
type (altLoc, bridging water, metal), built directly from the real
deposited coordinates in benchmark/pilot-inventory/raw-pdb/ -- no molecular
rendering engine (PyMOL requires native DLLs not available in this plain
pip environment; see PROJECT-ROADMAP.md for that tradeoff), just matplotlib
3D with bonds inferred by a covalent-distance cutoff. Same palette/type as
the rest of the figure set.

Cases (all already verified/discussed elsewhere in the project):
  - 1M17, ASP A831: the two altLoc conformers (0.50/0.50 occupancy) diverge
    from CB onward -- the geometric basis for the occupancy-tiebreak policy.
  - 1OHR, bridging water HOH A303: the classic HIV-protease "flap water",
    2.82 A from a ligand heteroatom and 2.88 A from Ile50 -- literally
    bridging ligand and receptor.
  - 1CBX, ZN A309: tetrahedral-ish coordination by His69/Glu72/His196 and
    the ligand's carboxylate (BZS), the reason metals are never dropped.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = REPO_ROOT / "benchmark" / "pilot-inventory" / "raw-pdb"
FIG_DIR = REPO_ROOT / "benchmark" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")

BLUE = "#2a78d6"
ORANGE = "#eb6834"
CRITICAL = "#d03b3b"
VIOLET = "#4a3aa7"
YELLOW = "#eda100"
INK = "#0b0b0b"
INK_MUTED = "#898781"
SURFACE = "#fcfcfb"

ELEMENT_COLOR = {"C": "#5a5a56", "N": BLUE, "O": CRITICAL, "S": YELLOW, "ZN": VIOLET}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9.5, "text.color": INK,
    "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def parse_atoms(pdb_id: str, matches) -> list[dict]:
    """matches: list of (resname, chain, resi_or_None, altloc_filter_or_None)."""
    out = []
    for line in (RAW_DIR / f"{pdb_id}.pdb").read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        rn, ch, ri, alt = line[17:20].strip(), line[21], line[22:26].strip(), line[16]
        for m_rn, m_ch, m_ri, m_alt in matches:
            if rn == m_rn and ch == m_ch and (m_ri is None or ri == m_ri) and (m_alt is None or alt in (" ", m_alt)):
                name = line[12:16].strip()
                element = "".join(c for c in name if c.isalpha())[:2]
                if element.upper() not in ("ZN",):
                    element = name[0] if name[0].isalpha() else name[1]
                pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                out.append({"name": name, "pos": pos, "resname": rn, "chain": ch, "resi": ri, "alt": alt,
                            "element": "ZN" if rn == "ZN" else element})
                break
    return out


def infer_bonds(atoms: list[dict], cutoff: float = 1.75) -> list[tuple[int, int]]:
    bonds = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if atoms[i]["element"] == "ZN" or atoms[j]["element"] == "ZN":
                continue
            d = np.linalg.norm(atoms[i]["pos"] - atoms[j]["pos"])
            if d <= cutoff:
                bonds.append((i, j))
    return bonds


def draw_sticks(ax, atoms, bonds, color_override=None, lw=2.6, s=55, alpha=1.0):
    positions = np.array([a["pos"] for a in atoms])
    for i, j in bonds:
        p, q = positions[i], positions[j]
        mid = (p + q) / 2
        c1 = color_override or ELEMENT_COLOR.get(atoms[i]["element"], INK_MUTED)
        c2 = color_override or ELEMENT_COLOR.get(atoms[j]["element"], INK_MUTED)
        ax.plot(*zip(p, mid), color=c1, linewidth=lw, solid_capstyle="round", alpha=alpha, zorder=2)
        ax.plot(*zip(mid, q), color=c2, linewidth=lw, solid_capstyle="round", alpha=alpha, zorder=2)
    for a in atoms:
        color = color_override or ELEMENT_COLOR.get(a["element"], INK_MUTED)
        ax.scatter(*a["pos"], color=color, s=s, edgecolors="white", linewidths=0.4, alpha=alpha, zorder=3)


LABEL_BOX = dict(facecolor="white", alpha=0.82, edgecolor="none", pad=1.2)


def dashed_distance(ax, p1, p2, label, color=INK_MUTED, frac=0.5, offset=(0.0, 0.0, 0.0), label_pos=None):
    ax.plot(*zip(p1, p2), color=color, linewidth=1.3, linestyle=(0, (3, 2)), zorder=1)
    pt = np.array(label_pos) if label_pos is not None else np.array(p1) + (np.array(p2) - np.array(p1)) * frac
    pt = pt + np.array(offset)
    ax.text(pt[0], pt[1], pt[2], label, fontsize=8, color=color, ha="center", zorder=5, bbox=LABEL_BOX)


def label_atom(ax, pos, text, color, offset=(0.0, 0.0, 0.4), fontsize=8.2, fontweight="normal"):
    p = np.array(pos) + np.array(offset)
    ax.text(p[0], p[1], p[2], text, fontsize=fontsize, color=color, ha="center", zorder=5,
            fontweight=fontweight, bbox=LABEL_BOX)


def set_equal_aspect(ax, all_points, elev=20, azim=-60, pad=1.12):
    pts = np.array(all_points)
    center = pts.mean(axis=0)
    span = max((pts.max(axis=0) - pts.min(axis=0)).max() / 2, 1.2) * pad
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()


def panel_altloc(ax):
    a_conf = parse_atoms("1M17", [("ASP", "A", "831", "A")])
    b_conf = parse_atoms("1M17", [("ASP", "A", "831", "B")])
    bonds_a = infer_bonds(a_conf)
    bonds_b = infer_bonds(b_conf)
    draw_sticks(ax, a_conf, bonds_a, color_override=BLUE, s=45)
    draw_sticks(ax, b_conf, bonds_b, color_override=ORANGE, s=45)

    od2_a = next(a["pos"] for a in a_conf if a["name"] == "OD2")
    od2_b = next(a["pos"] for a in b_conf if a["name"] == "OD2")
    d = np.linalg.norm(od2_a - od2_b)
    dashed_distance(ax, od2_a, od2_b, f"{d:.1f} Å", color=INK, offset=(0, 0, 0.6))

    set_equal_aspect(ax, [a["pos"] for a in a_conf + b_conf], elev=15, azim=25, pad=1.25)
    ax.set_title("ASP A:831 (1M17)\nconformeros A (azul, 0.50) / B (naranja, 0.50)", fontsize=9.3, fontweight="bold", pad=2)


def panel_water(ax):
    lig = parse_atoms("1OHR", [("1UN", "A", "201", None)])
    ile = parse_atoms("1OHR", [("ILE", "A", "50", None)])
    water = parse_atoms("1OHR", [("HOH", "A", "303", None)])

    lig_hetero = [a for a in lig if a["element"] in ("N", "O")]
    w_pos = water[0]["pos"]
    lig_near = min(lig_hetero, key=lambda a: np.linalg.norm(a["pos"] - w_pos))
    ile_near = min((a for a in ile if a["element"] in ("N", "O")), key=lambda a: np.linalg.norm(a["pos"] - w_pos))

    draw_sticks(ax, lig, infer_bonds(lig), color_override=None, s=40)
    draw_sticks(ax, ile, infer_bonds(ile), color_override=INK_MUTED, s=40)
    ax.scatter(*w_pos, color=CRITICAL, s=140, edgecolors="white", linewidths=0.6, zorder=4)
    label_atom(ax, w_pos, "H2O A:303", CRITICAL, offset=(0, 0, -0.8), fontweight="bold")

    d1 = np.linalg.norm(w_pos - lig_near["pos"])
    d2 = np.linalg.norm(w_pos - ile_near["pos"])
    dashed_distance(ax, w_pos, lig_near["pos"], f"{d1:.2f} Å", color=INK, frac=0.35, offset=(0, 0, 0.5))
    dashed_distance(ax, w_pos, ile_near["pos"], f"{d2:.2f} Å", color=INK, frac=0.65, offset=(0, 0, -0.3))

    set_equal_aspect(ax, [a["pos"] for a in lig + ile] + [w_pos], elev=12, azim=-40, pad=1.2)
    ax.set_title("Agua puente HOH A:303 (1OHR)\nligando 1UN ↔ Ile50 (flap de la proteasa)", fontsize=9.3, fontweight="bold", pad=2)


def panel_metal(ax):
    lig = parse_atoms("1CBX", [("BZS", "A", "500", None)])
    his69 = parse_atoms("1CBX", [("HIS", "A", "69", None)])
    glu72 = parse_atoms("1CBX", [("GLU", "A", "72", None)])
    his196 = parse_atoms("1CBX", [("HIS", "A", "196", None)])
    zn = parse_atoms("1CBX", [("ZN", "A", "309", None)])
    zn_pos = zn[0]["pos"]

    for res in (his69, glu72, his196):
        draw_sticks(ax, res, infer_bonds(res), color_override=INK_MUTED, s=38)
    draw_sticks(ax, lig, infer_bonds(lig), color_override=None, s=42)
    ax.scatter(*zn_pos, color=VIOLET, s=180, edgecolors="white", linewidths=0.6, zorder=4)

    coordinators = [
        ("1", "His69", next(a["pos"] for a in his69 if a["name"] == "ND1")),
        ("2", "Glu72", next(a["pos"] for a in glu72 if a["name"] == "OE2")),
        ("3", "His196", next(a["pos"] for a in his196 if a["name"] == "ND1")),
        ("4", "ligando O1", next(a["pos"] for a in lig if a["name"] == "O1")),
    ]
    dists = []
    for num, name, pos in coordinators:
        d = np.linalg.norm(zn_pos - pos)
        dists.append((num, name, d))
        direction = (np.array(pos) - zn_pos)
        direction = direction / np.linalg.norm(direction)
        label_pos = np.array(pos) + direction * 1.1
        ax.plot(*zip(zn_pos, pos), color=VIOLET, linewidth=1.3, linestyle=(0, (3, 2)), zorder=1)
        ax.text(*label_pos, num, fontsize=9, color="white", ha="center", va="center", zorder=5,
                fontweight="bold", bbox=dict(boxstyle="circle,pad=0.25", facecolor=VIOLET, edgecolor="none"))

    set_equal_aspect(ax, [a["pos"] for a in his69 + glu72 + his196 + lig] + [zn_pos], elev=35, azim=140, pad=1.4)
    ax.set_title("Zn(II) A:309 (1CBX)\nHis69/Glu72/His196 + carboxilato del ligando", fontsize=9.3, fontweight="bold", pad=2)
    legend = "   ".join(f"{n}={name} ({d:.2f} Å)" for n, name, d in dists)
    ax.text2D(0.5, -0.03, legend, transform=ax.transAxes, fontsize=7.6, color=VIOLET, ha="center")


def main() -> None:
    fig = plt.figure(figsize=(13.5, 5.2))
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax3 = fig.add_subplot(1, 3, 3, projection="3d")

    panel_altloc(ax1)
    panel_water(ax2)
    panel_metal(ax3)

    fig.suptitle("Figura — Geometría real del sitio de unión para cada tipo de hallazgo",
                 fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.02)
    fig.text(0.02, 0.965,
             "Coordenadas depositadas reales (no esquemáticas); enlaces inferidos por distancia de covalencia (≤1.75 Å).",
             fontsize=9, color=INK_MUTED, style="italic")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_binding_site_examples.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_binding_site_examples.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_binding_site_examples.png/.pdf")


if __name__ == "__main__":
    main()
