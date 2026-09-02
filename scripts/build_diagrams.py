"""Generate the 3 conceptual (non-data) figures from PROJECT-ROADMAP.md
section 10: graphical summary (fig 1), software/provenance architecture
(fig 3), and the reproducible-preparation decision tree (fig 9). Same
palette/typography as scripts/build_figures.py for a consistent figure set.

These are schematic, not data-driven -- built with matplotlib box/arrow
patches rather than an external diagramming tool, so they stay in the same
reproducible-by-code pipeline as the data figures (section 18).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = REPO_ROOT / "benchmark" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
VIOLET = "#4a3aa7"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
SURFACE = "#fcfcfb"
BOX_FILL = "#eef2f6"
GOOD = "#0ca30c"
CRITICAL = "#d03b3b"

BORDER = "#c7c7c2"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "text.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def box(ax, xy, w, h, text, accent=None, fontsize=8.3, fontweight="normal", textcolor=INK, title=None, title_size=8.6):
    """White box, thin neutral border, minimal corner rounding, optional thin
    left accent bar for category color -- restrained journal-figure style
    (accent tag + near-black text) rather than a filled pastel box."""
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
        linewidth=0.9, edgecolor=BORDER, facecolor="white", zorder=2,
    )
    ax.add_patch(patch)
    if accent:
        bar_w = min(0.045 * w, 0.05)
        bar = FancyBboxPatch((x, y), bar_w, h, boxstyle="round,pad=0,rounding_size=0.006",
                              linewidth=0, facecolor=accent, zorder=3)
        ax.add_patch(bar)
    text_x = x + w / 2 + (bar_w / 2 if accent else 0)
    if title:
        ax.text(text_x, y + h - h * 0.22, title, ha="center", va="center",
                fontsize=title_size, fontweight="bold", color=INK, zorder=4)
        ax.text(text_x, y + h * 0.38, text, ha="center", va="center",
                fontsize=fontsize, fontweight=fontweight, color=textcolor, zorder=4, linespacing=1.3)
    else:
        ax.text(text_x, y + h / 2, text, ha="center", va="center",
                fontsize=fontsize, fontweight=fontweight, color=textcolor, zorder=4, linespacing=1.3)
    return (x, y, w, h)


def arrow(ax, start, end, color=INK_MUTED, style="-|>", lw=0.9, connectionstyle="arc3,rad=0.0", label=None, label_offset=(0, 0.1)):
    a = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=8,
                         linewidth=lw, color=color, zorder=1, connectionstyle=connectionstyle,
                         shrinkA=0, shrinkB=0)
    ax.add_patch(a)
    if label:
        mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
        ax.text(mx + label_offset[0], my + label_offset[1], label, ha="center", va="bottom",
                fontsize=7.3, color=INK_MUTED, style="italic")


def fig_graphical_summary() -> None:
    fig, ax = plt.subplots(figsize=(12.5, 3.7))
    ax.set_xlim(0, 12.5); ax.set_ylim(0, 3.7); ax.axis("off")

    stages = [
        ("Estructura PDB\ndepositada", "fuente RCSB,\nchecksum registrado", INK_MUTED),
        ("DockPrep Audit", "altLoc · aguas ·\nmetal no estándar", BLUE),
        ("Decisión de\npolítica explícita", "registrada, no\nautomática", ORANGE),
        ("Preparación +\nredocking", "Meeko, luego\nAutoDock Vina", AQUA),
        ("Resultado\nclasificado", "success / scoring_fail\n/ sampling_fail", INK_MUTED),
    ]
    w, h, gap = 2.1, 1.35, 0.4
    x = 0.15
    y = 1.55
    positions = []
    for title, sub, accent in stages:
        pos = box(ax, (x, y - h / 2), w, h, sub, accent=accent, title=title, fontsize=7.8, title_size=8.8)
        positions.append(pos)
        x += w + gap

    for i in range(len(positions) - 1):
        x0, y0, w0, h0 = positions[i]
        x1, y1, w1, h1 = positions[i + 1]
        arrow(ax, (x0 + w0, y0 + h0 / 2), (x1, y1 + h1 / 2))

    ax.text(0.02, 3.35,
            "¿Las características estructurales detectables antes de preparar un receptor predicen cuándo\n"
            "distintas políticas de preparación cambian la reproducibilidad del docking?",
            fontsize=9.6, fontweight="bold", color=INK, ha="left", va="top")
    ax.text(0.02, 0.35,
            "Cada flecha es una decisión documentada, no un paso automático silencioso.",
            fontsize=8, color=INK_MUTED, ha="left", va="top", style="italic")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_graphical_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_graphical_summary.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_graphical_summary.png/.pdf")


def fig_architecture() -> None:
    fig, ax = plt.subplots(figsize=(9.2, 10.2))
    ax.set_xlim(0, 9.2); ax.set_ylim(0, 10.2); ax.axis("off")

    col_l, col_r = 0.4, 5.0
    w_wide, w_narrow = 3.7, 3.8

    chain = [
        ("Fuente pública", "files.rcsb.org/download/{id}.pdb", INK_MUTED, 9.0, 0.85),
        ("build_pilot_inventory.py", "checksum SHA-256 · category()", BLUE, 7.7, 0.85),
        ("dockprep_audit (v0.1.0)", "altLoc · aguas · metal → JSON", BLUE, 6.4, 0.85),
        ("verify_pilot_eligibility.py", "relevancia espacial ≤6 Å · overrides", ORANGE, 5.1, 0.85),
        ("freeze_pilot_manifest.py", "manifiesto congelado + checksum", ORANGE, 3.8, 0.85),
        ("smoke_redock_case.py", "Meeko (PDBQT) → Vina → RMSD\n(matching exacto de coordenadas)", AQUA, 2.3, 1.05),
        ("run_phase_c*.py / run_phase_d.py", "n semillas × políticas → agregado", AQUA, 0.9, 0.85),
    ]
    left_chain = []
    for title, sub, accent, y, h in chain:
        pos = box(ax, (col_l, y), w_wide, h, sub, accent=accent, title=title, fontsize=7.4, title_size=7.9)
        left_chain.append(pos)

    for i in range(len(left_chain) - 1):
        x0, y0, w0, h0 = left_chain[i]
        x1, y1, w1, h1 = left_chain[i + 1]
        arrow(ax, (x0 + w0 / 2, y0), (x1 + w1 / 2, y1 + h1))

    box(ax, (col_r, 5.7), w_narrow, 2.0,
        "altLoc → mayor ocupancia\naguas → puente ≤3.0 Å\nmetal → siempre conservar\nprotonación → pH 7.4",
        accent=YELLOW, title="Decisiones de política\n(sección 7/8, congeladas\nantes de ver resultados)",
        fontsize=7.4, title_size=7.7)
    arrow(ax, (col_l + w_wide, 6.4 + 0.42), (col_r, 6.7), color=INK_MUTED, connectionstyle="arc3,rad=-0.15")

    box(ax, (col_r, 3.75), w_narrow, 1.3,
        "receptor_chain_overrides.csv\n(3FNU, 1SN5 — REMARK 350 vs.\nregla geométrica de 6 Å)",
        accent=CRITICAL, title="Overrides documentados", fontsize=7.4, title_size=7.7)
    arrow(ax, (col_l + w_wide, 5.1 + 0.42), (col_r, 4.4), color=INK_MUTED, connectionstyle="arc3,rad=-0.1")

    box(ax, (col_r, 0.9), w_narrow, 2.0,
        "phase_c_summary.csv\nphase_d_by_case.csv\nphase_d_by_stratum.csv\nfiguras (build_figures.py)",
        accent=INK_MUTED, title="Salidas versionadas", fontsize=7.4, title_size=7.7)
    arrow(ax, (col_l + w_wide, 1.32), (col_r, 1.6), color=INK_MUTED, connectionstyle="arc3,rad=-0.15")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_architecture.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_architecture.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_architecture.png/.pdf")


def fig_decision_tree() -> None:
    fig, ax = plt.subplots(figsize=(12.5, 7.8))
    ax.set_xlim(0, 12.5); ax.set_ylim(0, 7.8); ax.axis("off")

    root_w = 2.9
    box(ax, (4.8, 6.6), root_w, 0.8, "detectado por la auditoría", accent=BLUE,
        title="Hallazgo estructural", fontsize=7.8, title_size=8.2, fontweight="bold")

    col_w = 2.5
    yn_w = 1.5
    branches = [
        (0.15, "Conformación alternativa\n(altLoc)", "¿Ocupancia desigual?",
         "conservar la de\nmayor ocupancia", "conservar 'A'\ndeclarado (empate)",
         "Registro por residuo en\naltloc_decisions.json"),
        (4.35, "Agua cristalográfica", "¿Puente ligando–\nreceptor ≤3.0 Å?",
         "conservar\n(política 'conservador')", "descartar\n(política 'simplificado')",
         "Ambas políticas se corren\ny comparan — no se elige\nuna sola de antemano"),
        (8.55, "Metal / cofactor o\nresiduo no estándar", "¿Meeko tiene\nplantilla nativa?",
         "conservar, carga nativa\n(Zn/Mg/Ca/Mn/Fe, ALY)", "falla documentada\n(STA, ligando covalente)",
         "Ej.: pepstatina (STA) sin\nplantilla → sustitución\ncon justificación"),
    ]

    for x0, title, question, yes_text, no_text, note in branches:
        col_c = x0 + col_w / 2
        arrow(ax, (6.25, 6.6), (col_c, 5.85), color=INK_MUTED, connectionstyle="arc3,rad=0.0")
        box(ax, (x0, 5.0), col_w, 0.75, title, accent=INK_MUTED, fontweight="bold", fontsize=7.9)
        arrow(ax, (col_c, 5.0), (col_c, 4.65))
        box(ax, (x0, 3.75), col_w, 0.75, question, accent=YELLOW, fontsize=7.7)

        yes_x = col_c - 0.9 - yn_w / 2
        no_x = col_c + 0.9 - yn_w / 2
        arrow(ax, (col_c - 0.35, 3.75), (yes_x + yn_w / 2, 2.75), color=GOOD, connectionstyle="arc3,rad=-0.25",
              label="sí", label_offset=(-0.3, 0.0))
        box(ax, (yes_x, 1.85), yn_w, 0.9, yes_text, accent=GOOD, fontsize=6.9)

        arrow(ax, (col_c + 0.35, 3.75), (no_x + yn_w / 2, 2.75), color=CRITICAL, connectionstyle="arc3,rad=0.25",
              label="no", label_offset=(0.3, 0.0))
        box(ax, (no_x, 1.85), yn_w, 0.9, no_text, accent=CRITICAL, fontsize=6.9)

        ax.text(col_c, 1.5, note, ha="center", va="top", fontsize=6.7, color=INK_MUTED, style="italic", linespacing=1.3)

    fig.text(0.01, 0.02, "Cada rama termina en una decisión registrada en un archivo versionado, nunca en una corrección silenciosa.",
              fontsize=7.6, color=INK_MUTED, style="italic")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_decision_tree.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_decision_tree.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_decision_tree.png/.pdf")


def main() -> None:
    fig_graphical_summary()
    fig_architecture()
    fig_decision_tree()
    print(f"\nAll diagrams written to {FIG_DIR}")


if __name__ == "__main__":
    main()
