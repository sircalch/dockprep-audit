"""Roadmap fig 7 (v3): real ball-and-stick 3D rendering (PyVista/VTK) of the
binding-site geometry, with labels placed as a 2D callout overlay (matplotlib)
instead of embedded 3D text -- so every label sits in clear whitespace with a
thin leader line back to its atom, never covering part of the structure.
See build_site_figures.py (matplotlib-only v1) and PROJECT-ROADMAP.md for the
earlier iterations and why PyMOL itself was not usable in this environment.

Same 3 cases, same coordinates/distances as before:
  - 1M17, ASP A831: two altLoc conformers (0.50/0.50), OD2-OD2 diverge 4.3 A.
  - 1OHR, HOH A303: classic HIV-protease "flap water", 2.82 A from the
    ligand, 2.88 A from Ile50.
  - 1CBX, ZN A309: coordination by His69/Glu72/His196 + ligand carboxylate.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = REPO_ROOT / "benchmark" / "pilot-inventory" / "raw-pdb"
FIG_DIR = REPO_ROOT / "benchmark" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

pv.OFF_SCREEN = True
fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")
plt.rcParams.update({"font.family": "Arial"})

BLUE = "#2a78d6"
ORANGE = "#eb6834"
CRITICAL = "#d03b3b"
VIOLET = "#4a3aa7"
YELLOW = "#eda100"
CARBON = "#8a8a86"
INK = "#0b0b0b"

ELEMENT_COLOR = {"C": CARBON, "N": BLUE, "O": CRITICAL, "S": YELLOW, "ZN": VIOLET}
ELEMENT_RADIUS = {"C": 0.30, "N": 0.28, "O": 0.28, "S": 0.34, "ZN": 0.48}

WIN = 1700


def parse_atoms(pdb_id: str, matches) -> list[dict]:
    out = []
    for line in (RAW_DIR / f"{pdb_id}.pdb").read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        rn, ch, ri, alt = line[17:20].strip(), line[21], line[22:26].strip(), line[16]
        for m_rn, m_ch, m_ri, m_alt in matches:
            if rn == m_rn and ch == m_ch and (m_ri is None or ri == m_ri) and (m_alt is None or alt in (" ", m_alt)):
                name = line[12:16].strip()
                element = "".join(c for c in name if c.isalpha())[:2]
                if element.upper() != "ZN":
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
            if np.linalg.norm(atoms[i]["pos"] - atoms[j]["pos"]) <= cutoff:
                bonds.append((i, j))
    return bonds


def add_atoms(plotter, atoms, color_override=None, radius_scale=1.0):
    for a in atoms:
        color = color_override or ELEMENT_COLOR.get(a["element"], CARBON)
        r = ELEMENT_RADIUS.get(a["element"], 0.3) * radius_scale
        sphere = pv.Sphere(radius=r, center=a["pos"], theta_resolution=40, phi_resolution=40)
        plotter.add_mesh(sphere, color=color, smooth_shading=True, specular=0.4, specular_power=20)


def add_bonds(plotter, atoms, bonds, color_override=None, radius=0.09):
    for i, j in bonds:
        p, q = atoms[i]["pos"], atoms[j]["pos"]
        mid = (p + q) / 2
        c1 = color_override or ELEMENT_COLOR.get(atoms[i]["element"], CARBON)
        c2 = color_override or ELEMENT_COLOR.get(atoms[j]["element"], CARBON)
        for start, end, color in ((p, mid, c1), (mid, q, c2)):
            height = np.linalg.norm(end - start)
            if height < 1e-6:
                continue
            direction = (end - start) / height
            cyl = pv.Cylinder(center=(start + end) / 2, direction=direction, radius=radius, height=height,
                               resolution=28)
            plotter.add_mesh(cyl, color=color, smooth_shading=True)


def add_dashed_3d(plotter, p1, p2, color, n_dashes=10, line_width=3):
    p1, p2 = np.array(p1), np.array(p2)
    for k in range(n_dashes):
        if k % 2 == 1:
            continue
        t0, t1 = k / n_dashes, (k + 0.6) / n_dashes
        seg = pv.Line(p1 + (p2 - p1) * t0, p1 + (p2 - p1) * t1)
        plotter.add_mesh(seg, color=color, line_width=line_width)


def project(plotter, world_pos) -> np.ndarray:
    coord = vtkCoordinate()
    coord.SetCoordinateSystemToWorld()
    coord.SetValue(*world_pos)
    x, y = coord.GetComputedDisplayValue(plotter.renderer)
    return np.array([x, WIN - y])  # flip to image (row-0-at-top) pixel space


def render_panel(build_fn) -> tuple[np.ndarray, dict]:
    """build_fn(plotter) populates the scene and returns a dict of named 3D
    anchor points; render_panel screenshots the clean scene (no embedded
    text) and projects those anchors to 2D image pixel coordinates."""
    plotter = pv.Plotter(off_screen=True, window_size=(WIN, WIN), border=False)
    plotter.set_background("white")
    anchors, camera_fn = build_fn(plotter)
    camera_fn(plotter)
    plotter.enable_ssao(radius=1.5)
    plotter.enable_anti_aliasing("ssaa")
    plotter.render()
    img = plotter.screenshot(return_img=True)
    anchors_2d = {name: project(plotter, pos) for name, pos in anchors.items()}
    plotter.close()
    return img, anchors_2d


def build_altloc(plotter):
    a_conf = parse_atoms("1M17", [("ASP", "A", "831", "A")])
    b_conf = parse_atoms("1M17", [("ASP", "A", "831", "B")])
    add_atoms(plotter, a_conf, color_override=BLUE, radius_scale=0.85)
    add_bonds(plotter, a_conf, infer_bonds(a_conf), color_override=BLUE, radius=0.075)
    add_atoms(plotter, b_conf, color_override=ORANGE, radius_scale=0.85)
    add_bonds(plotter, b_conf, infer_bonds(b_conf), color_override=ORANGE, radius=0.075)
    od2_a = next(a["pos"] for a in a_conf if a["name"] == "OD2")
    od2_b = next(a["pos"] for a in b_conf if a["name"] == "OD2")
    add_dashed_3d(plotter, od2_a, od2_b, "black")

    def cam(p):
        p.camera_position = "iso"
        p.camera.azimuth = 25
        p.camera.elevation = 10
        p.reset_camera()
        p.camera.zoom(1.3)

    return {"dist_od2": (od2_a + od2_b) / 2}, cam


def build_water(plotter):
    lig = parse_atoms("1OHR", [("1UN", "A", "201", None)])
    ile = parse_atoms("1OHR", [("ILE", "A", "50", None)])
    water = parse_atoms("1OHR", [("HOH", "A", "303", None)])
    w_pos = water[0]["pos"]
    lig_hetero = [a for a in lig if a["element"] in ("N", "O")]
    lig_near = min(lig_hetero, key=lambda a: np.linalg.norm(a["pos"] - w_pos))
    ile_near = min((a for a in ile if a["element"] in ("N", "O")), key=lambda a: np.linalg.norm(a["pos"] - w_pos))

    add_atoms(plotter, lig, radius_scale=0.85)
    add_bonds(plotter, lig, infer_bonds(lig))
    add_atoms(plotter, ile, color_override=CARBON, radius_scale=0.85)
    add_bonds(plotter, ile, infer_bonds(ile), color_override=CARBON)
    add_atoms(plotter, water, color_override=CRITICAL, radius_scale=1.3)
    add_dashed_3d(plotter, w_pos, lig_near["pos"], "black")
    add_dashed_3d(plotter, w_pos, ile_near["pos"], "black")

    def cam(p):
        p.camera_position = "iso"
        p.camera.azimuth = -30
        p.reset_camera()
        p.camera.zoom(0.92)

    return {
        "water": w_pos,
        "dist_lig": (w_pos + lig_near["pos"]) / 2,
        "dist_ile": (w_pos + ile_near["pos"]) / 2,
    }, cam


def build_metal(plotter):
    lig = parse_atoms("1CBX", [("BZS", "A", "500", None)])
    his69 = parse_atoms("1CBX", [("HIS", "A", "69", None)])
    glu72 = parse_atoms("1CBX", [("GLU", "A", "72", None)])
    his196 = parse_atoms("1CBX", [("HIS", "A", "196", None)])
    zn = parse_atoms("1CBX", [("ZN", "A", "309", None)])
    zn_pos = zn[0]["pos"]

    for res in (his69, glu72, his196):
        add_atoms(plotter, res, color_override=CARBON, radius_scale=0.8)
        add_bonds(plotter, res, infer_bonds(res), color_override=CARBON, radius=0.07)
    add_atoms(plotter, lig, radius_scale=0.85)
    add_bonds(plotter, lig, infer_bonds(lig))
    add_atoms(plotter, zn, color_override=VIOLET, radius_scale=1.2)

    coordinators = [
        ("his69", next(a["pos"] for a in his69 if a["name"] == "ND1")),
        ("glu72", next(a["pos"] for a in glu72 if a["name"] == "OE2")),
        ("his196", next(a["pos"] for a in his196 if a["name"] == "ND1")),
        ("ligand_o1", next(a["pos"] for a in lig if a["name"] == "O1")),
    ]
    anchors = {"zn": zn_pos}
    for key, pos in coordinators:
        add_dashed_3d(plotter, zn_pos, pos, VIOLET)
        anchors[f"dist_{key}"] = zn_pos + (pos - zn_pos) * 0.55

    def cam(p):
        p.camera_position = "iso"
        p.camera.azimuth = 60
        p.reset_camera()
        p.camera.zoom(0.82)

    return anchors, cam


def box_whiteness(img, cx, cy, box_w, box_h) -> float:
    """Fraction of near-white pixels in the image region centered at
    (cx, cy) -- used to test whether a candidate label position is clear
    of any rendered structure (or of a previously placed label)."""
    h, w = img.shape[0], img.shape[1]
    x0, x1 = int(cx - box_w / 2), int(cx + box_w / 2)
    y0, y1 = int(cy - box_h / 2), int(cy + box_h / 2)
    if x0 < 0 or y0 < 0 or x1 > w or y1 > h:
        return 0.0
    patch = img[y0:y1, x0:x1]
    return float(np.mean(np.all(patch > 235, axis=-1)))


def find_clear_spot(img, anchor_2d, text_len, occupied, min_r_frac=0.10, max_r_frac=0.42):
    """Search rings of candidate directions at increasing radius from the
    anchor; return the first candidate whose label-sized box is >=98% clear
    of both the rendered structure and any already-placed label box."""
    h, w = img.shape[0], img.shape[1]
    box_w, box_h = 15 * text_len + 20, 34
    n_angles = 16
    for r_frac in np.linspace(min_r_frac, max_r_frac, 9):
        r = r_frac * min(h, w)
        for k in range(n_angles):
            angle = 2 * np.pi * k / n_angles
            cx = anchor_2d[0] + r * np.cos(angle)
            cy = anchor_2d[1] + r * np.sin(angle)
            if box_whiteness(img, cx, cy, box_w, box_h) < 0.98:
                continue
            box = (cx - box_w / 2, cy - box_h / 2, cx + box_w / 2, cy + box_h / 2)
            if any(not (box[2] < o[0] or box[0] > o[2] or box[3] < o[1] or box[1] > o[3]) for o in occupied):
                continue
            occupied.append(box)
            return np.array([cx, cy])
    # fallback: farthest ring, ignore label-overlap constraint
    r = max_r_frac * min(h, w)
    return anchor_2d + np.array([r, -r])


def draw_callouts(ax, img, anchors_2d, labels):
    """labels: list of (anchor_key, text, color) or (anchor_key, text, color,
    (px, py)) with an explicit label pixel position. Explicit positions were
    hand-picked from the real projected anchor coordinates (printed via
    render_panel) after the automatic clear-spot search proved too fragile
    once >=4 anchors cluster tightly (metal panel) -- it kept finding
    technically-white but mutually overlapping slots. Anchor-only entries
    still go through find_clear_spot."""
    ax.imshow(img)
    ax.axis("off")
    occupied = []
    for entry in labels:
        key, text, color = entry[0], entry[1], entry[2]
        ax_pt = anchors_2d[key]
        if len(entry) == 4:
            label_pt = np.array(entry[3], dtype=float)
        else:
            label_pt = find_clear_spot(img, ax_pt, len(text), occupied)
        ax.annotate(
            text, xy=ax_pt, xytext=label_pt,
            fontsize=13, fontweight="bold", color=color, ha="center", va="center",
            arrowprops=dict(arrowstyle="-", color=color, linewidth=1.1, shrinkA=0, shrinkB=4),
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.92),
        )


def main() -> None:
    img_a, anc_a = render_panel(build_altloc)
    img_w, anc_w = render_panel(build_water)
    img_m, anc_m = render_panel(build_metal)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.0))

    draw_callouts(axes[0], img_a, anc_a, [
        ("dist_od2", "4.3 A", INK),
    ])
    axes[0].set_title("ASP A:831 (1M17)\nconformers A (blue) / B (orange), 0.50/0.50",
                       fontsize=11, fontweight="bold", color=INK, loc="left")

    draw_callouts(axes[1], img_w, anc_w, [
        ("water", "H2O A:303 (red sphere)", CRITICAL, (150, 1150)),
        ("dist_lig", "2.82 A", INK, (1350, 700)),
        ("dist_ile", "2.88 A", INK, (1000, 1350)),
    ])
    axes[1].set_title("Bridging water HOH A:303 (1OHR)\nligand 1UN <-> Ile50 (protease flap)",
                       fontsize=11, fontweight="bold", color=INK, loc="left")

    draw_callouts(axes[2], img_m, anc_m, [
        ("zn", "Zn(II) (violet sphere)", VIOLET, (1350, 900)),
        ("dist_his69", "His69 2.04 A", VIOLET, (150, 500)),
        ("dist_glu72", "Glu72 2.03 A", VIOLET, (1350, 300)),
        ("dist_his196", "His196 2.00 A", VIOLET, (1350, 1250)),
        ("dist_ligand_o1", "ligand O1 2.30 A", VIOLET, (150, 1250)),
    ])
    axes[2].set_title("Zn(II) A:309 (1CBX)\nHis69/Glu72/His196 + ligand carboxylate",
                       fontsize=11, fontweight="bold", color=INK, loc="left")

    fig.tight_layout()
    out_png = FIG_DIR / "fig_binding_site_examples_pv.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
