"""Two more figures from the second review-batch, built from already-frozen
Fase D data (median best-of-9 RMSD per case x policy), no new experiments:

  1. fig_cumulative_success.png      - cumulative success-fraction curves
                                        (fraction of cases with RMSD <= x) vs.
                                        RMSD threshold, conservative vs.
                                        simplified, water_policy stratum and
                                        all 60 cases combined.
  2. fig_threshold_sensitivity.png   - the water-policy success-fraction gap
                                        (conservative - simplified) as a
                                        function of the success threshold,
                                        showing the 2.0 A threshold used in
                                        the paper is not cherry-picked.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = REPO_ROOT / "benchmark" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
BY_CASE_CSV = REPO_ROOT / "benchmark" / "phase-c-60" / "phase_d_by_case.csv"

fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")

BLUE = "#2a78d6"
AMBER = "#e2a63f"
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

THRESHOLD_USED = 2.0


def load():
    rows = list(csv.DictReader(open(BY_CASE_CSV, encoding="utf-8")))
    by_stratum_policy = defaultdict(lambda: defaultdict(list))
    for r in rows:
        rmsd = float(r["median_best_pose_rmsd_A"])
        by_stratum_policy[r["stratum"]][r["policy"]].append(rmsd)
        by_stratum_policy["__all__"][r["policy"]].append(rmsd)
    return by_stratum_policy


def save(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {name}.png/.pdf")


def cumulative_curve(rmsds, thresholds):
    rmsds = np.asarray(rmsds)
    return np.array([(rmsds <= t).mean() for t in thresholds])


def fig1_cumulative_success(data):
    thresholds = np.linspace(0, 10, 201)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)

    panels = [("water_policy", "Water-policy stratum (n=15)"), ("__all__", "All 60 cases combined")]
    for ax, (key, title) in zip(axes, panels):
        cons = cumulative_curve(data[key]["conservative_water"], thresholds)
        simp = cumulative_curve(data[key]["simplified_no_water"], thresholds)
        ax.plot(thresholds, cons, color=BLUE, linewidth=2, label="Conservative (water retained)")
        ax.plot(thresholds, simp, color=AMBER, linewidth=2, label="Simplified (water removed)")
        ax.fill_between(thresholds, cons, simp, color=INK_MUTED, alpha=0.12, linewidth=0)
        ax.axvline(THRESHOLD_USED, color=INK, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.set_xlabel("RMSD threshold (\u00c5)")
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1.0)

    axes[0].set_ylabel("Cumulative success fraction\n(cases with best-of-9 RMSD \u2264 threshold)")
    axes[0].text(THRESHOLD_USED + 0.15, 0.05, "2.0 \u00c5\n(used in\nthis study)", fontsize=7.5, color=INK_MUTED)
    axes[1].legend(loc="lower right", fontsize=8.5, frameon=False)
    fig.suptitle("Cumulative docking-success curves are consistently higher under the conservative policy",
                 fontsize=11, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "fig_cumulative_success")


def fig2_threshold_sensitivity(data):
    thresholds = np.arange(0.5, 5.01, 0.1)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    stratum_colors = {
        "water_policy": "#3fa66e",
        "alternate_location": BLUE,
        "metal_or_cofactor": AMBER,
        "low_risk_control": "#8a8f98",
    }
    stratum_labels = {
        "water_policy": "Water policy",
        "alternate_location": "Alternate location",
        "metal_or_cofactor": "Metal / cofactor",
        "low_risk_control": "Low-risk control",
    }
    for stratum in ["water_policy", "alternate_location", "metal_or_cofactor", "low_risk_control"]:
        cons = cumulative_curve(data[stratum]["conservative_water"], thresholds)
        simp = cumulative_curve(data[stratum]["simplified_no_water"], thresholds)
        gap = cons - simp
        lw = 2.4 if stratum == "water_policy" else 1.3
        alpha = 1.0 if stratum == "water_policy" else 0.65
        ax.plot(thresholds, gap, color=stratum_colors[stratum], linewidth=lw, alpha=alpha,
                label=stratum_labels[stratum], zorder=3 if stratum == "water_policy" else 2)

    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvline(THRESHOLD_USED, color=INK_MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(THRESHOLD_USED + 0.05, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 0.2,
            "2.0 \u00c5 (used in this study)", fontsize=7.5, color=INK_MUTED, va="top")

    ax.set_xlabel("Success threshold (\u00c5)")
    ax.set_ylabel("Success-fraction gap\n(conservative \u2212 simplified)")
    ax.set_title("Success-fraction gap by stratum, as a function of the success threshold",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    ax.set_xlim(thresholds.min(), thresholds.max())

    fig.tight_layout()
    save(fig, "fig_threshold_sensitivity")


def main():
    data = load()
    fig1_cumulative_success(data)
    fig2_threshold_sensitivity(data)


if __name__ == "__main__":
    main()
