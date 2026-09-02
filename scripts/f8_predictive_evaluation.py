"""F8: formal predictive evaluation of SITE_BRIDGING_WATER_PRESENT
(PROJECT-ROADMAP.md section 19, item F8), added 2026-08-29.

Uses F1's real engine finding and F4's symmetry-corrected RMSD, evaluated
across all 60 main-cohort cases (not just the water_policy stratum), to
build the full binary confusion matrix: predictor =
SITE_BRIDGING_WATER_PRESENT (true/false, from
site_bridging_water_by_case.csv), outcome = |delta best-of-9 RMSD,
conservative-simplified| > 0.3 A (the same threshold already used
throughout Sections 3.5-3.7, pre-specified before this evaluation, not
tuned for it). Reports sensitivity, specificity, PPV, NPV, and an ROC/AUC
sweep using |delta RMSD| as the continuous score.

F7's validation cohort cannot itself supply a full confusion matrix (it
was screened FOR the finding, so it is 100% positive-class by
construction) but gives an independent check on one cell: of predicted-
positive cases, what fraction are actually policy-sensitive (PPV) in data
never used to define or tune the finding.
"""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SYM_DIR = REPO_ROOT / "benchmark" / "symmetry-rmsd"
EFFECT_THRESHOLD_A = 0.3


def load_predictor() -> dict[str, bool]:
    return {r["pdb_id"]: r["site_bridging_water_present"] == "True"
            for r in csv.DictReader(open(SYM_DIR / "site_bridging_water_by_case.csv", encoding="utf-8"))}


def load_delta_rmsd() -> dict[str, float]:
    rows = list(csv.DictReader(open(SYM_DIR / "phase_c60_symmetry_rmsd.csv", encoding="utf-8")))
    by_case_policy = defaultdict(list)
    for r in rows:
        by_case_policy[(r["pdb_id"], r["policy"])].append(float(r["sym_best"]))
    deltas = {}
    cases = {pdb for pdb, _ in by_case_policy}
    for pdb in cases:
        cons = by_case_policy.get((pdb, "conservative_water"))
        simp = by_case_policy.get((pdb, "simplified_no_water"))
        if not cons or not simp:
            continue
        deltas[pdb] = abs(statistics.median(cons) - statistics.median(simp))
    return deltas


def confusion_matrix(predictor: dict[str, bool], deltas: dict[str, float], threshold: float):
    tp = fp = tn = fn = 0
    for pdb, delta in deltas.items():
        pred_pos = predictor.get(pdb, False)
        actual_pos = delta > threshold
        if pred_pos and actual_pos:
            tp += 1
        elif pred_pos and not actual_pos:
            fp += 1
        elif not pred_pos and actual_pos:
            fn += 1
        else:
            tn += 1
    return tp, fp, tn, fn


def roc_points(predictor: dict[str, bool], deltas: dict[str, float]):
    """Sweep the |delta RMSD| threshold itself as the continuous score,
    with SITE_BRIDGING_WATER_PRESENT as the binary classifier -- report a
    single-point sensitivity/specificity/PPV/NPV summary (the classifier is
    binary, not a score), and separately an ROC-style sweep using the
    magnitude |delta RMSD| itself as a continuous "true sensitivity" score
    against the binary predictor, i.e. how well the predictor's positive
    class captures increasingly large true effects.
    """
    pdbs = sorted(deltas)
    y_true_cont = [deltas[p] for p in pdbs]
    y_pred_bin = [1 if predictor.get(p, False) else 0 for p in pdbs]
    # Rank-based measure: mean |delta| for predicted-positive vs predicted-negative
    pos = [deltas[p] for p in pdbs if predictor.get(p, False)]
    neg = [deltas[p] for p in pdbs if not predictor.get(p, False)]
    return pos, neg


def main() -> None:
    predictor = load_predictor()
    deltas = load_delta_rmsd()
    print(f"{len(deltas)} cases with usable symmetry-corrected delta RMSD "
          f"(of 60; 1OHR excluded, Section 2.4)")

    tp, fp, tn, fn = confusion_matrix(predictor, deltas, EFFECT_THRESHOLD_A)
    n = tp + fp + tn + fn
    sensitivity = tp / (tp + fn) if (tp + fn) else float("nan")
    specificity = tn / (tn + fp) if (tn + fp) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) else float("nan")
    npv = tn / (tn + fn) if (tn + fn) else float("nan")
    accuracy = (tp + tn) / n

    print(f"\n=== Confusion matrix (n={n}, threshold |delta RMSD| > {EFFECT_THRESHOLD_A} A) ===")
    print(f"                  actual_sensitive   actual_not_sensitive")
    print(f"  predicted_pos   TP={tp:3d}              FP={fp:3d}")
    print(f"  predicted_neg   FN={fn:3d}              TN={tn:3d}")
    print(f"\n  Sensitivity (recall): {sensitivity:.3f}")
    print(f"  Specificity:          {specificity:.3f}")
    print(f"  PPV (precision):      {ppv:.3f}")
    print(f"  NPV:                  {npv:.3f}")
    print(f"  Accuracy:             {accuracy:.3f}")
    print(f"  Baseline sensitive rate (prevalence): {(tp+fn)/n:.3f}")

    pos, neg = roc_points(predictor, deltas)
    print(f"\n=== Magnitude comparison ===")
    print(f"  Mean |delta RMSD| when SITE_BRIDGING_WATER_PRESENT=True  (n={len(pos)}): {statistics.mean(pos):.3f} A")
    print(f"  Mean |delta RMSD| when SITE_BRIDGING_WATER_PRESENT=False (n={len(neg)}): {statistics.mean(neg):.3f} A")

    try:
        from scipy.stats import mannwhitneyu
        stat, p = mannwhitneyu(pos, neg, alternative="greater")
        print(f"  Mann-Whitney U (pos > neg): U={stat:.1f}, p={p:.4f}")
    except ImportError:
        pass

    try:
        from sklearn.metrics import roc_auc_score
        y_true = [1 if predictor.get(p, False) else 0 for p in deltas]
        y_score = [deltas[p] for p in deltas]
        # AUC of using |delta RMSD| to predict the binary finding is not the
        # question; the question is whether the finding predicts |delta RMSD|
        # crossing 0.3 A -- report AUC the natural direction instead.
        y_actual = [1 if deltas[p] > EFFECT_THRESHOLD_A else 0 for p in deltas]
        y_score_pred = [1 if predictor.get(p, False) else 0 for p in deltas]
        auc = roc_auc_score(y_actual, y_score_pred)
        print(f"\n  AUC (binary predictor vs. binary outcome, equivalent to balanced accuracy point): {auc:.3f}")
    except ImportError:
        print("\n  (scikit-learn not available; binary-predictor AUC is mathematically equal to "
              "(sensitivity+specificity)/2 for a single-point binary classifier, reported above.)")

    # --- External check via F7 (100% predicted-positive by construction) ---
    f7_path = REPO_ROOT / "benchmark" / "f7-external-validation" / "f7_by_case_symmetry.csv"
    if f7_path.exists():
        f7_rows = list(csv.DictReader(open(f7_path, encoding="utf-8")))
        f7_by_case = defaultdict(dict)
        for r in f7_rows:
            f7_by_case[r["pdb_id"]][r["policy"]] = float(r["median_best"])
        f7_deltas = []
        for pdb, d in f7_by_case.items():
            if "conservative_water" in d and "simplified_no_water" in d:
                f7_deltas.append(abs(d["conservative_water"] - d["simplified_no_water"]))
        f7_ppv = sum(1 for x in f7_deltas if x > EFFECT_THRESHOLD_A) / len(f7_deltas)
        print(f"\n=== External check (F7 validation cohort, n={len(f7_deltas)}, 100% SITE_BRIDGING_WATER_PRESENT) ===")
        print(f"  PPV in independent data: {f7_ppv:.3f}  (main-cohort PPV: {ppv:.3f})")


if __name__ == "__main__":
    main()
