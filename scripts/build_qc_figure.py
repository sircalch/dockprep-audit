"""Roadmap fig, Results 3.2 support: deposition year vs. RSCC availability
for the 48 Fase E expansion cases, showing the "not a clean 1998 cutoff"
finding directly rather than only in prose. Same palette/typography as the
rest of the figure set.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = REPO_ROOT / "benchmark" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")

BLUE = "#2a78d6"
CRITICAL = "#d03b3b"
INK = "#0b0b0b"
INK_MUTED = "#898781"
GRID = "#e1e0d9"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 10,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_MUTED, "ytick.color": INK_MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
})


def main() -> None:
    rows = list(csv.DictReader(open(REPO_ROOT / "benchmark" / "expansion-validation" / "wwpdb_validation.csv", encoding="utf-8")))
    have = [(int(r["deposition_year"]), float(r["rscc"])) for r in rows if r["rscc"].strip()]
    missing = [int(r["deposition_year"]) for r in rows if not r["rscc"].strip()]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    years_h, rscc_h = zip(*have)
    ax.scatter(years_h, rscc_h, s=34, color=BLUE, alpha=0.85, edgecolors="white", linewidths=0.5,
               label=f"RSCC available (n={len(have)})", zorder=3)
    ax.scatter(missing, [0.79] * len(missing), s=34, color=CRITICAL, marker="x", linewidths=1.6,
               label=f"No structure-factor file (n={len(missing)})", zorder=3)
    ax.axhline(0.8, color=INK_MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(1991.5, 0.805, "field reliability threshold (RSCC = 0.8)", fontsize=7.5, color=INK_MUTED, va="bottom")
    ax.axvline(1998, color=INK_MUTED, linewidth=0.8, linestyle=(0, (1, 2)))
    ax.text(1998.3, 0.98, "1998", fontsize=7.5, color=INK_MUTED, ha="left")

    ax.set_xlabel("Deposition year")
    ax.set_ylabel("RSCC")
    ax.set_ylim(0.75, 1.0)
    ax.legend(loc="lower right", fontsize=8.5, frameon=False)
    ax.set_title("wwPDB validation coverage, 48 expansion-cohort cases", fontsize=11, fontweight="bold", loc="left")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_qc_rscc_by_year.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_qc_rscc_by_year.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_qc_rscc_by_year.png/.pdf")


if __name__ == "__main__":
    main()
