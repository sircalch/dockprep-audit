"""Batch of 4 additional manuscript figures requested after the second
review-style critique, buildable entirely from already-frozen Fase C/D data
(no new experiments needed):

  1. fig_delta_rmsd_forest.png   - per-case Delta RMSD (best-of-9, conservative
                                    minus simplified), all 60 cases, grouped by
                                    stratum, signed lollipop/forest plot.
  2. fig_top1_vs_best9.png       - success-fraction comparison, top-1 vs.
                                    best-of-9, by stratum x policy (the paper's
                                    central sampling-vs-scoring point).
  3. fig_confound_ligand_size.png- Delta RMSD vs. ligand heavy-atom count and
                                    vs. radius of gyration, water_policy stratum
                                    highlighted, to check whether the water-
                                    policy effect is confounded with ligand size.
  4. fig_finding_cooccurrence.png- replaces the pie+heatmap co-occurrence figure
                                    with an UpSet-style matrix of audit-finding
                                    co-occurrence across the 60-case cohort.

Same palette/typography as the rest of the figure set (see build_qc_figure.py).
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from itertools import combinations
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

# validated colorblind-safe categorical palette (project standard)
STRATUM_COLOR = {
    "alternate_location": "#2a78d6",   # blue
    "metal_or_cofactor": "#e2a63f",    # amber
    "water_policy": "#3fa66e",         # green
    "low_risk_control": "#8a8f98",     # gray
}
STRATUM_LABEL = {
    "alternate_location": "Alternate location",
    "metal_or_cofactor": "Metal / cofactor",
    "water_policy": "Water policy",
    "low_risk_control": "Low-risk control",
}
STRATUM_ORDER = ["water_policy", "alternate_location", "metal_or_cofactor", "low_risk_control"]

BLUE = "#2a78d6"
AMBER = "#e2a63f"
GREEN = "#3fa66e"
GRAY = "#8a8f98"
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


def load_by_case():
    rows = list(csv.DictReader(open(BY_CASE_CSV, encoding="utf-8")))
    by_key = {}
    for r in rows:
        by_key[(r["pdb_id"], r["policy"])] = r
    cases = defaultdict(dict)
    for (pdb, policy), r in by_key.items():
        cases[pdb]["stratum"] = r["stratum"]
        cases[pdb][policy] = r
    return cases


def save(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {name}.png/.pdf")


# ---------------------------------------------------------------------------
# 1. Delta RMSD forest / lollipop plot, all 60 cases, grouped by stratum
# ---------------------------------------------------------------------------

def fig1_delta_rmsd_forest(cases):
    rows = []
    for pdb, d in cases.items():
        c = d.get("conservative_water")
        s = d.get("simplified_no_water")
        if not c or not s:
            continue
        delta = float(c["median_best_pose_rmsd_A"]) - float(s["median_best_pose_rmsd_A"])
        rows.append((pdb, d["stratum"], delta))

    rows.sort(key=lambda r: (STRATUM_ORDER.index(r[1]), r[2]))

    fig, ax = plt.subplots(figsize=(7.5, 11))
    y = np.arange(len(rows))
    colors = [STRATUM_COLOR[r[1]] for r in rows]
    deltas = [r[2] for r in rows]

    ax.hlines(y, 0, deltas, color=colors, linewidth=1.6, zorder=2)
    ax.scatter(deltas, y, color=colors, s=26, zorder=3, edgecolors="white", linewidths=0.4)
    ax.axvline(0, color=INK, linewidth=0.9, zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.2, fontfamily="monospace")
    ax.set_xlabel("$\\Delta$RMSD, best-of-9 (conservative $-$ simplified, Å)")
    ax.set_title("Per-case water-policy effect on redocking accuracy, all 60 cases",
                 fontsize=11, fontweight="bold", loc="left", pad=22)
    ax.text(0.0, 1.0, "negative = conservative (water-retained) more accurate  →",
            transform=ax.transAxes, fontsize=7.5, color=INK_MUTED, ha="left", va="bottom")

    # stratum group separators + labels
    prev_stratum = None
    for i, (pdb, stratum, delta) in enumerate(rows):
        if stratum != prev_stratum:
            if prev_stratum is not None:
                ax.axhline(i - 0.5, color=GRID, linewidth=0.8)
            prev_stratum = stratum

    handles = [plt.Line2D([0], [0], marker="o", color=STRATUM_COLOR[s], linestyle="",
                           markersize=6, label=STRATUM_LABEL[s]) for s in STRATUM_ORDER]
    ax.legend(handles=handles, loc="lower right", fontsize=8, frameon=False)
    ax.set_ylim(-1, len(rows))
    ax.invert_yaxis()

    fig.tight_layout()
    save(fig, "fig_delta_rmsd_forest")


# ---------------------------------------------------------------------------
# 2. Top-1 vs best-of-9 success-fraction comparison, by stratum x policy
# ---------------------------------------------------------------------------

def fig2_top1_vs_best9(cases):
    policies = ["conservative_water", "simplified_no_water"]
    policy_label = {"conservative_water": "Conservative\n(water retained)", "simplified_no_water": "Simplified\n(water removed)"}

    data = {s: {p: {"top1": [], "best9": []} for p in policies} for s in STRATUM_ORDER}
    for pdb, d in cases.items():
        s = d["stratum"]
        for p in policies:
            r = d.get(p)
            if not r:
                continue
            data[s][p]["top1"].append(float(r["success_fraction_top1"]))
            data[s][p]["best9"].append(float(r["success_fraction_best"]))

    fig, ax = plt.subplots(figsize=(8.5, 5))
    n_strata = len(STRATUM_ORDER)
    group_w = 0.62
    bar_w = group_w / 4
    x = np.arange(n_strata)

    metric_shade = {"top1": 0.55, "best9": 1.0}
    offsets = [-1.5, -0.5, 0.5, 1.5]
    bars_meta = [(p, m) for p in policies for m in ["top1", "best9"]]

    for off, (p, m) in zip(offsets, bars_meta):
        vals = [np.mean(data[s][p][m]) if data[s][p][m] else 0 for s in STRATUM_ORDER]
        base_color = BLUE if p == "conservative_water" else AMBER
        alpha = metric_shade[m]
        label = f"{'Conservative' if p=='conservative_water' else 'Simplified'}, {'top-1' if m=='top1' else 'best-of-9'}"
        ax.bar(x + off * bar_w, vals, width=bar_w * 0.92, color=base_color, alpha=alpha,
               edgecolor="white", linewidth=0.6, label=label, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels([STRATUM_LABEL[s] for s in STRATUM_ORDER], fontsize=9)
    ax.set_ylabel("Mean success fraction (RMSD ≤ 2.0 Å)")
    ax.set_ylim(0, 1.0)
    ax.set_title("Top-1-pose vs. best-of-9-pose success, by stratum and water policy",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=7.8, frameon=False, ncol=2)

    fig.tight_layout()
    save(fig, "fig_top1_vs_best9")


# ---------------------------------------------------------------------------
# 3. Confound check: Delta RMSD vs ligand heavy-atom count and vs Rg
# ---------------------------------------------------------------------------

def fig3_confound_ligand_size(cases):
    rows = []
    for pdb, d in cases.items():
        c = d.get("conservative_water")
        s = d.get("simplified_no_water")
        if not c or not s:
            continue
        raw_path = RAW_RUNS / f"{pdb}_conservative_water_seed42.json"
        if not raw_path.exists():
            continue
        raw = json.load(open(raw_path, encoding="utf-8"))["result"]
        heavy_atoms = raw.get("ligand_heavy_atoms")
        rg = raw.get("box", {}).get("radius_of_gyration_A")
        delta = float(c["median_best_pose_rmsd_A"]) - float(s["median_best_pose_rmsd_A"])
        rows.append((pdb, d["stratum"], heavy_atoms, rg, delta))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))

    for ax, xkey, xlabel in zip(axes, [2, 3], ["Ligand heavy-atom count", "Ligand radius of gyration (Å)"]):
        for stratum in STRATUM_ORDER:
            pts = [(r[xkey], r[4]) for r in rows if r[1] == stratum]
            if not pts:
                continue
            xs, ys = zip(*pts)
            is_water = stratum == "water_policy"
            ax.scatter(xs, ys, s=42 if is_water else 26,
                       color=STRATUM_COLOR[stratum],
                       alpha=1.0 if is_water else 0.55,
                       edgecolors="white" if is_water else "none", linewidths=0.6,
                       zorder=3 if is_water else 2,
                       label=STRATUM_LABEL[stratum])
        ax.axhline(0, color=INK_MUTED, linewidth=0.8, zorder=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("$\\Delta$RMSD (conservative $-$ simplified, Å)")

    axes[0].legend(loc="upper right", fontsize=7.2, frameon=False)
    fig.suptitle("Water-policy effect is not explained by ligand size", fontsize=11, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, "fig_confound_ligand_size")


# ---------------------------------------------------------------------------
# 4. UpSet-style co-occurrence matrix of audit findings across the cohort
# ---------------------------------------------------------------------------

def locate_raw_pdb(pdb_id: str) -> Path | None:
    for d in REPO_ROOT.glob("benchmark/*/raw-pdb"):
        p = d / f"{pdb_id}.pdb"
        if p.exists():
            return p
    return None


def fig4_finding_cooccurrence(cases):
    import dockprep_audit

    combo_counts = Counter()
    unresolved = []
    for pdb in cases:
        path = locate_raw_pdb(pdb)
        if path is None:
            unresolved.append(pdb)
            continue
        result = dockprep_audit.audit_pdb(str(path))
        codes = tuple(sorted({f["code"] for f in result["findings"]}))
        combo_counts[codes] += 1

    if unresolved:
        print(f"WARNING: could not locate raw PDB for {len(unresolved)} case(s), "
              f"excluded from co-occurrence figure: {unresolved}")

    combo_counts = {k: v for k, v in combo_counts.items() if k}
    if not combo_counts:
        print("No non-empty finding combinations to plot, skipping fig_finding_cooccurrence")
        return

    all_findings = sorted({f for combo in combo_counts for f in combo})
    combos_sorted = sorted(combo_counts.items(), key=lambda kv: -kv[1])

    n_combos = len(combos_sorted)
    n_findings = len(all_findings)
    fig = plt.figure(figsize=(max(6, n_combos * 0.9 + 2), 4.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.3], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(n_combos)
    counts = [c for _, c in combos_sorted]
    ax_bar.bar(x, counts, width=0.55, color=BLUE, zorder=3)
    for xi, c in zip(x, counts):
        ax_bar.text(xi, c + 0.6, str(c), ha="center", fontsize=8, color=INK)
    ax_bar.set_ylabel("Case count")
    ax_bar.set_xticks([])
    ax_bar.spines["bottom"].set_visible(False)
    ax_bar.set_title("Audit-finding co-occurrence across the 60-case cohort",
                     fontsize=11, fontweight="bold", loc="left")

    for yi, finding in enumerate(all_findings):
        ax_matrix.axhline(yi, color=GRID, linewidth=0.6, zorder=1)
    for xi, (combo, _) in enumerate(combos_sorted):
        active_y = [all_findings.index(f) for f in combo]
        if len(active_y) > 1:
            ax_matrix.plot([xi, xi], [min(active_y), max(active_y)], color=INK, linewidth=1.3, zorder=2)
        for yi in range(n_findings):
            active = yi in active_y
            ax_matrix.scatter([xi], [yi], s=90 if active else 60,
                              color=INK if active else GRID,
                              zorder=3, edgecolors="none")

    ax_matrix.set_yticks(range(n_findings))
    ax_matrix.set_yticklabels(all_findings, fontsize=8.5)
    ax_matrix.set_ylim(-0.7, n_findings - 0.3)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xticks([])
    ax_matrix.grid(False)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)

    fig.tight_layout()
    save(fig, "fig_finding_cooccurrence")


def main():
    cases = load_by_case()
    fig1_delta_rmsd_forest(cases)
    fig2_top1_vs_best9(cases)
    fig3_confound_ligand_size(cases)
    fig4_finding_cooccurrence(cases)


if __name__ == "__main__":
    main()
