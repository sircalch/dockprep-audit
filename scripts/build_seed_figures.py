"""Two more figures from the second review-batch, both requiring per-seed
raw-run data (already frozen in benchmark/phase-c-60/raw_runs, no new
experiments needed):

  1. fig_seed_stability.png   - per-case spread of best-of-9 RMSD across the
                                 3 predeclared seeds (42, 123, 2024), one line
                                 per case x policy, sorted by median RMSD,
                                 with the 15 seed-unstable-outcome cases
                                 flagged directly on the plot.
  2. fig_pose_rank_profiles.png - full 9-pose RMSD-by-rank profile (seed 42)
                                 for the 6 most water-sensitive cases (largest
                                 |Delta median best-of-9 RMSD| in the
                                 water_policy stratum), conservative vs.
                                 simplified, small-multiples layout.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = REPO_ROOT / "benchmark" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
RAW_RUNS = REPO_ROOT / "benchmark" / "phase-c-60" / "raw_runs"
BY_CASE_CSV = REPO_ROOT / "benchmark" / "phase-c-60" / "phase_d_by_case.csv"

fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")

BLUE = "#2a78d6"
AMBER = "#e2a63f"
GREEN = "#3fa66e"
GRAY = "#8a8f98"
CRITICAL = "#d03b3b"
INK = "#0b0b0b"
INK_MUTED = "#898781"
GRID = "#e1e0d9"

STRATUM_COLOR = {
    "alternate_location": BLUE,
    "metal_or_cofactor": AMBER,
    "water_policy": GREEN,
    "low_risk_control": GRAY,
}
STRATUM_LABEL = {
    "alternate_location": "Alternate location",
    "metal_or_cofactor": "Metal / cofactor",
    "water_policy": "Water policy",
    "low_risk_control": "Low-risk control",
}
STRATUM_ORDER = ["water_policy", "alternate_location", "metal_or_cofactor", "low_risk_control"]
SEEDS = [42, 123, 2024]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 10,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_MUTED, "ytick.color": INK_MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
})


def save(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {name}.png/.pdf")


def load_by_case():
    rows = list(csv.DictReader(open(BY_CASE_CSV, encoding="utf-8")))
    cases = defaultdict(dict)
    for r in rows:
        cases[r["pdb_id"]]["stratum"] = r["stratum"]
        cases[r["pdb_id"]][r["policy"]] = r
    return cases


def load_raw(pdb_id, policy, seed):
    path = RAW_RUNS / f"{pdb_id}_{policy}_seed{seed}.json"
    return json.load(open(path, encoding="utf-8"))["result"]


# ---------------------------------------------------------------------------
# 1. Seed stability: best-of-9 RMSD spread across 3 seeds, all 60 cases
# ---------------------------------------------------------------------------

def fig1_seed_stability(cases):
    rows = []
    for pdb, d in cases.items():
        stratum = d["stratum"]
        for policy in ["conservative_water", "simplified_no_water"]:
            vals = [load_raw(pdb, policy, s)["best_pose_rmsd_A"] for s in SEEDS]
            unstable = d[policy]["unstable_across_seeds"] == "True"
            rows.append({
                "pdb": pdb, "stratum": stratum, "policy": policy,
                "vals": vals, "median": float(np.median(vals)),
                "spread": max(vals) - min(vals), "unstable": unstable,
            })

    rows.sort(key=lambda r: (STRATUM_ORDER.index(r["stratum"]), r["pdb"], r["policy"]))

    fig, ax = plt.subplots(figsize=(7.8, 13))
    y = np.arange(len(rows))
    for yi, r in zip(y, rows):
        color = STRATUM_COLOR[r["stratum"]]
        marker_face = CRITICAL if r["unstable"] else color
        ax.plot([min(r["vals"]), max(r["vals"])], [yi, yi], color=color, linewidth=1.3, alpha=0.55, zorder=2)
        ax.scatter(r["vals"], [yi] * 3, color=marker_face, s=16, zorder=3,
                   edgecolors="white", linewidths=0.3)

    labels = [f"{r['pdb']}  ({'C' if r['policy']=='conservative_water' else 'S'})" for r in rows]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.3, fontfamily="monospace")
    ax.set_xlabel("Best-of-9 RMSD across 3 seeds (42, 123, 2024), \u00c5")
    ax.set_title("Per-case, per-policy seed stability, all 60 cases \u00d7 2 policies",
                 fontsize=11, fontweight="bold", loc="left")
    ax.set_ylim(-1, len(rows))
    ax.invert_yaxis()

    handles = [plt.Line2D([0], [0], marker="o", color=STRATUM_COLOR[s], linestyle="",
                          markersize=6, label=STRATUM_LABEL[s]) for s in STRATUM_ORDER]
    handles.append(plt.Line2D([0], [0], marker="o", color=CRITICAL, linestyle="",
                              markersize=6, label="Unstable outcome across seeds (n=15/120)"))
    ax.legend(handles=handles, loc="lower right", fontsize=7.5, frameon=False)

    prev = None
    for yi, r in zip(y, rows):
        if r["stratum"] != prev and prev is not None:
            ax.axhline(yi - 0.5, color=GRID, linewidth=0.8)
        prev = r["stratum"]

    fig.tight_layout()
    save(fig, "fig_seed_stability")


# ---------------------------------------------------------------------------
# 2. 9-pose rank profiles for the 6 most water-sensitive cases
# ---------------------------------------------------------------------------

def fig2_pose_rank_profiles(cases):
    water_cases = []
    for pdb, d in cases.items():
        if d["stratum"] != "water_policy":
            continue
        c = float(d["conservative_water"]["median_best_pose_rmsd_A"])
        s = float(d["simplified_no_water"]["median_best_pose_rmsd_A"])
        water_cases.append((pdb, abs(c - s), c, s))

    water_cases.sort(key=lambda x: -x[1])
    top6 = water_cases[:6]

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.4), sharey=False)

    for ax, (pdb, delta, c_med, s_med) in zip(axes.flat, top6):
        c_raw = load_raw(pdb, "conservative_water", 42)
        s_raw = load_raw(pdb, "simplified_no_water", 42)
        c_vals = c_raw["rmsd_by_rank_A"]
        s_vals = s_raw["rmsd_by_rank_A"]
        ax.plot(np.arange(1, len(c_vals) + 1), c_vals, color=BLUE, marker="o", markersize=4,
                linewidth=1.6, label="Conservative")
        ax.plot(np.arange(1, len(s_vals) + 1), s_vals, color=AMBER, marker="o", markersize=4,
                linewidth=1.6, label="Simplified")
        ax.axhline(2.0, color=INK_MUTED, linewidth=0.7, linestyle=(0, (3, 2)))
        ax.set_title(f"{pdb}  (\u0394median={delta:.2f} \u00c5)", fontsize=9.5, loc="left")
        ax.set_xticks([1, 3, 5, 7, 9])
        ax.set_xlabel("Vina pose rank", fontsize=8)
        ax.set_ylabel("RMSD (\u00c5)", fontsize=8)
        ax.tick_params(labelsize=7.5)

    axes.flat[0].legend(loc="upper right", fontsize=7.5, frameon=False)
    fig.suptitle("Full 9-pose RMSD-by-rank profiles, 6 most water-sensitive cases (seed 42)",
                 fontsize=12, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "fig_pose_rank_profiles")


def main():
    cases = load_by_case()
    fig1_seed_stability(cases)
    fig2_pose_rank_profiles(cases)


if __name__ == "__main__":
    main()
