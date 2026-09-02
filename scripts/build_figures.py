"""Generate publication-quality figures from frozen Fase C/D data.

Figures are built by code from already-frozen data (checksum
`ab37330651f5787bc25b4a7a301f7b7c4a98bfdd580139aa027898d6e5d4abfe`), per
PROJECT-ROADMAP.md section 18 ("las figuras se generen por codigo a partir
de datos congelados"). Colors follow the validated categorical/sequential
palette in the project's dataviz reference (blue/orange for the 2
preparation policies, matching the CVD-safe fixed order), not an arbitrary
choice.

Three figures are produced (roadmap section 10, items 4/5/6). Item 2
(cohort selection/exclusion funnel) is deliberately NOT included here: the
raw candidate-search bookkeeping across ~20 ad-hoc search batches was never
deduplicated into one authoritative count, and fabricating a precise-looking
funnel number from messy intermediate files would misrepresent uncertain
data as exact -- better to leave it out than publish an unverified figure.
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

fm.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
fm.fontManager.addfont("C:/Windows/Fonts/arialbd.ttf")

# Validated categorical palette (dataviz skill reference, light mode) --
# fixed order, not cycled: slot 1 blue, slot 2 orange, slot 3 aqua, slot 4 yellow.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
VIOLET = "#4a3aa7"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"

POLICY_COLOR = {"conservative_water": BLUE, "simplified_no_water": ORANGE}
POLICY_LABEL = {"conservative_water": "Con agua puente", "simplified_no_water": "Sin agua"}
STRATUM_COLOR = {
    "alternate_location": VIOLET, "metal_or_cofactor": AQUA,
    "water_policy": BLUE, "low_risk_control": INK_MUTED,
}
STRATUM_LABEL = {
    "alternate_location": "Conf. alternativa", "metal_or_cofactor": "Metal/cofactor",
    "water_policy": "Política de aguas", "low_risk_control": "Control bajo riesgo",
}
STRATUM_ORDER = ["alternate_location", "metal_or_cofactor", "water_policy", "low_risk_control"]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 10.5,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_SECONDARY, "ytick.color": INK_SECONDARY,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False,
})


def load_strata() -> dict[str, str]:
    with (REPO_ROOT / "benchmark" / "pilot_cases.csv").open(newline="", encoding="utf-8") as f:
        return {row["pdb_id"]: row["stratum"] for row in csv.DictReader(f)}


def fig_finding_frequency(strata: dict[str, str]) -> None:
    """Roadmap fig 4: frequency and co-occurrence of structural findings."""
    reports_dir = REPO_ROOT / "benchmark" / "pilot-inventory" / "reports"
    codes_by_case = {}
    for pdb_id in strata:
        data = json.loads((reports_dir / f"{pdb_id}.json").read_text(encoding="utf-8"))
        codes_by_case[pdb_id] = {f["code"] for f in data["findings"]}

    all_codes = ["ALTLOC_PRESENT", "WATERS_PRESENT", "METAL_PRESENT"]
    code_label = {"ALTLOC_PRESENT": "Conf. alternativa", "WATERS_PRESENT": "Aguas", "METAL_PRESENT": "Metal"}
    counts = {c: sum(1 for codes in codes_by_case.values() if c in codes) for c in all_codes}

    co = np.zeros((3, 3), dtype=int)
    for codes in codes_by_case.values():
        present = [c in codes for c in all_codes]
        for i in range(3):
            for j in range(3):
                if present[i] and present[j]:
                    co[i, j] += 1

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.2), gridspec_kw={"width_ratios": [1, 1.15]})

    bars = ax1.barh(
        [code_label[c] for c in all_codes], [counts[c] for c in all_codes],
        color=[BLUE, ORANGE, AQUA], height=0.55, zorder=3,
    )
    ax1.set_xlim(0, 60)
    ax1.set_xlabel("N.º de casos (de 60)")
    ax1.set_title("Frecuencia de hallazgos estructurales", fontsize=11, fontweight="bold", loc="left")
    ax1.grid(axis="y", visible=False)
    for bar, code in zip(bars, all_codes):
        ax1.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2, str(counts[code]),
                  va="center", ha="left", fontsize=9.5, color=INK)

    im = ax2.imshow(co, cmap="Blues", vmin=0, vmax=60)
    ax2.set_xticks(range(3)); ax2.set_xticklabels([code_label[c] for c in all_codes], rotation=30, ha="right")
    ax2.set_yticks(range(3)); ax2.set_yticklabels([code_label[c] for c in all_codes])
    ax2.set_title("Coocurrencia (n.º de casos con ambos)", fontsize=11, fontweight="bold", loc="left")
    for i in range(3):
        for j in range(3):
            color = "white" if co[i, j] > 30 else INK
            ax2.text(j, i, str(co[i, j]), ha="center", va="center", color=color, fontsize=10)
    ax2.grid(visible=False)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=8, color=INK_MUTED)
    cbar.outline.set_visible(False)

    fig.suptitle("Figura — Hallazgos estructurales en la cohorte de 60 casos", fontsize=12.5, fontweight="bold", x=0.02, ha="left", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_finding_frequency.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_finding_frequency.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_finding_frequency.png/.pdf")


def fig_heatmap(strata: dict[str, str]) -> None:
    """Roadmap fig 5: heatmap case x policy, colored by best-of-9-pose RMSD."""
    rows = list(csv.DictReader((REPO_ROOT / "benchmark" / "phase-c-60" / "phase_d_by_case.csv").open(encoding="utf-8")))
    by_case = defaultdict(dict)
    for r in rows:
        by_case[r["pdb_id"]][r["policy"]] = float(r["median_best_pose_rmsd_A"])

    ordered_cases = sorted(strata, key=lambda p: (STRATUM_ORDER.index(strata[p]), p))
    policies = ["conservative_water", "simplified_no_water"]
    matrix = np.array([[by_case[c][p] for p in policies] for c in ordered_cases])

    fig, ax = plt.subplots(figsize=(4.6, 13.5))
    im = ax.imshow(matrix, cmap="Blues_r", vmin=0, vmax=6, aspect="auto")
    ax.set_xticks(range(2)); ax.set_xticklabels([POLICY_LABEL[p] for p in policies], rotation=20, ha="right")
    ax.set_yticks(range(len(ordered_cases))); ax.set_yticklabels(ordered_cases, fontsize=7.2, fontfamily="monospace")
    ax.grid(visible=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    prev_stratum = None
    for i, c in enumerate(ordered_cases):
        s = strata[c]
        if s != prev_stratum:
            if prev_stratum is not None:
                ax.axhline(i - 0.5, color=INK, linewidth=1.1)
            prev_stratum = s
        for j in range(2):
            v = matrix[i, j]
            color = "white" if v < 2.0 else INK
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.3, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("RMSD mejor pose (Å) — más bajo = mejor", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5)
    cbar.outline.set_visible(False)

    ax.set_title("Figura — RMSD (mejor de 9 poses) por caso y política\n60 casos, agrupados por estrato de riesgo",
                  fontsize=11, fontweight="bold", loc="left", pad=12)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_heatmap_case_policy.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_heatmap_case_policy.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_heatmap_case_policy.png/.pdf")


def fig_rmsd_distribution(strata: dict[str, str]) -> None:
    """Roadmap fig 6: RMSD distribution by policy, faceted by stratum."""
    raw_dir = REPO_ROOT / "benchmark" / "phase-c-60" / "raw_runs"
    values = defaultdict(lambda: defaultdict(list))
    for p in raw_dir.glob("*.json"):
        rec = json.loads(p.read_text(encoding="utf-8"))
        if not rec.get("result"):
            continue
        pdb_id = rec["pdb_id"]
        if pdb_id not in strata:
            continue
        values[strata[pdb_id]][rec["policy"]].append(rec["result"]["best_pose_rmsd_A"])

    policies = ["conservative_water", "simplified_no_water"]
    fig, axes = plt.subplots(1, 4, figsize=(13, 4.6), sharey=True)
    for ax, stratum in zip(axes, STRATUM_ORDER):
        data = [values[stratum][p] for p in policies]
        parts = ax.violinplot(data, positions=[0, 1], widths=0.7, showmedians=True, showextrema=False)
        for body, policy in zip(parts["bodies"], policies):
            body.set_facecolor(POLICY_COLOR[policy])
            body.set_alpha(0.55)
            body.set_edgecolor(POLICY_COLOR[policy])
            body.set_linewidth(1.0)
        parts["cmedians"].set_color(INK)
        parts["cmedians"].set_linewidth(1.4)
        for i, p in enumerate(policies):
            xs = np.random.default_rng(42).normal(i, 0.05, size=len(data[i]))
            ax.scatter(xs, data[i], s=9, color=POLICY_COLOR[p], alpha=0.5, zorder=3, edgecolors="none")
        ax.axhline(2.0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Con\nagua", "Sin\nagua"], fontsize=9)
        ax.set_title(STRATUM_LABEL[stratum], fontsize=10.5, fontweight="bold")
        ax.set_ylim(-0.3, 11.5)
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("RMSD, mejor de 9 poses (Å)")
    axes[0].text(1.55, 2.15, "umbral de éxito (2.0 Å)", fontsize=7.5, color=INK_MUTED, ha="center")

    fig.suptitle("Figura — Distribución de RMSD por política de preparación, n=15/estrato",
                 fontsize=12.5, fontweight="bold", x=0.01, ha="left", y=1.03)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_rmsd_distribution.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_rmsd_distribution.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote fig_rmsd_distribution.png/.pdf")


def main() -> None:
    strata = load_strata()
    fig_finding_frequency(strata)
    fig_heatmap(strata)
    fig_rmsd_distribution(strata)
    print(f"\nAll figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
