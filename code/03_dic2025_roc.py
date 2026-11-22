# -*- coding: utf-8 -*-
"""
02_dic2025_roc.py -- ROC of Day-1 neutrophil markers for overt DIC (ISTH 2025)
=====================================================================
Outcome: `DIC_2025` (computed by 01_dic2025_scoring.py; total >= 5, complete cases only)
Predictors: the five Day-1 neutrophil markers (same set as Figure 2B)
    CD63 MFI / CD62L MFI / CD66b MFI / CD16 MFI / CD10- neutrophil percentage

Statistics follow `code/tables.py::table3` so the results are comparable with Table 3:
    - markers with AUROC < 0.5 are flipped before plotting;
    - 95% CI: patient-level bootstrap, B = 2000, seed = 20260507;
    - operating point: Youden optimum in the full sample, held fixed during resampling;
    - sensitivity / specificity intervals: Clopper-Pearson.

Usage:
    python code/02_dic2025_roc.py                # default outcome: overt DIC
    python code/02_dic2025_roc.py --outcome mortality   # alternative: 28-day mortality

Output:
    results/tables/DIC2025_roc_5biomarkers.csv
    figures/DIC2025_roc_5biomarkers_day1.png / .pdf
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from statsmodels.stats.proportion import proportion_confint

# ------------------------------------------------------------------ paths
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
for _ in range(6):
    if os.path.isdir(os.path.join(ROOT, "data")) and os.path.isdir(os.path.join(ROOT, "code")):
        break
    ROOT = os.path.dirname(ROOT)
DATA = os.path.join(ROOT, "data")
OUT_T = os.path.join(ROOT, "results", "tables")
OUT_F = os.path.join(ROOT, "figures")
os.makedirs(OUT_T, exist_ok=True)
os.makedirs(OUT_F, exist_ok=True)

# ------------------------------------------------------------------ markers
BM = [
    ("All_Neu(CD63_MFI)", "CD63 MFI", "CD63", "#B2182B"),
    ("All_Neu_Mean_(CD62L_MFI)", "CD62L MFI", "CD62L", "#990099"),
    ("All_Neu_Mean_(CD66b_MFI)", "CD66b MFI", "CD66b", "#339C62"),
    ("All_Neu_Mean_(CD16_MFI)", "CD16 MFI", "CD16", "#E69F00"),
    ("immature_Neu_percentage", "CD10$^{-}$ neutrophils (%)", "iNeu", "#1A4A8C"),
]

SEED = 20260507
N_BOOT = 2000


def _load_scorer():
    """Load 01_dic2025_scoring.py as a module."""
    path = os.path.join(HERE, "01_dic2025_scoring.py")
    spec = importlib.util.spec_from_file_location("dic2025", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _metrics(x, y, flipped, t_opt):
    """Operating-point metrics on the given arrays."""
    pred = ((-x if flipped else x) >= t_opt).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    npv = tn / (tn + fn) if (tn + fn) else np.nan
    f1 = 2 * ppv * sens / (ppv + sens) if (ppv + sens) else np.nan
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, Sens=sens, Spec=spec, PPV=ppv, NPV=npv, F1=f1)


def analyse(cohort, col, y):
    """AUROC / AUPRC / operating point / bootstrap CI for one marker."""
    d = cohort[[col]].copy()
    x = d[col].values.astype(float)
    fpr, tpr, thrs = roc_curve(y, x)
    flipped = auc(fpr, tpr) < 0.5
    fpr, tpr, thrs = roc_curve(y, -x if flipped else x)
    roc_auc = auc(fpr, tpr)
    t_opt = float(thrs[int(np.argmax(tpr - fpr))])
    cut = -t_opt if flipped else t_opt
    m = _metrics(x, y, flipped, t_opt)
    pr, rc, _ = precision_recall_curve(y, -x if flipped else x)
    auprc = auc(rc, pr)

    rng = np.random.default_rng(SEED)
    boot = {k: [] for k in ["AUC", "AUPRC", "Sens", "Spec", "PPV", "NPV", "F1"]}
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(x), len(x))
        xb, yb = x[idx], y[idx]
        if len(np.unique(yb)) < 2:
            continue
        fb, tb, _ = roc_curve(yb, -xb if flipped else xb)
        boot["AUC"].append(auc(fb, tb))
        pb, rb, _ = precision_recall_curve(yb, -xb if flipped else xb)
        boot["AUPRC"].append(auc(rb, pb))
        mb = _metrics(xb, yb, flipped, t_opt)          # fixed threshold
        for k in ["Sens", "Spec", "PPV", "NPV", "F1"]:
            boot[k].append(mb[k])
    ci = {k: tuple(float(v) for v in np.nanpercentile(vals, [2.5, 97.5]))
          for k, vals in boot.items()}
    for k in ["Sens", "Spec", "PPV", "NPV"]:
        num, den = {"Sens": (m["tp"], m["tp"] + m["fn"]), "Spec": (m["tn"], m["tn"] + m["fp"]),
                    "PPV": (m["tp"], m["tp"] + m["fp"]), "NPV": (m["tn"], m["tn"] + m["fn"])}[k]
        if den:
            ci[k] = tuple(float(v) for v in proportion_confint(int(num), int(den), method="beta"))
    m.update(AUC=roc_auc, AUPRC=auprc, ci=ci, cutoff=cut, flipped=flipped,
             direction="<=" if flipped else ">=", n=len(x), fpr=fpr, tpr=tpr)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outcome", choices=["dic", "mortality"], default="dic")
    args = ap.parse_args()

    d1 = pd.read_csv(os.path.join(DATA, "NE_baseline_D1.csv"), encoding="utf-8-sig")
    sc = _load_scorer()
    d1, _coag_src = sc.apply_backfill(d1)      # same Day-2 backfill as step 01
    res = sc.score_frame(d1)

    if args.outcome == "dic":
        keep = res["DIC_2025"].notna()
        cohort = d1.loc[keep.values].copy().reset_index(drop=True)
        y = res.loc[keep, "DIC_2025"].astype(int).values
        out_label = "overt DIC (ISTH 2025, score \u22655)"
        stem = "DIC2025_roc_5biomarkers_day1"
        series = [(c, lab, short, col) for c, lab, short, col in BM]
    else:
        mort = pd.to_numeric(d1["28-day mortality"], errors="coerce")
        keep = mort.notna()
        cohort = d1.loc[keep.values].copy().reset_index(drop=True)
        y = mort[keep].astype(int).values
        out_label = "28-day mortality"
        stem = "DIC2025_roc_mortality_day1"
        series = [(c, lab, short, col) for c, lab, short, col in BM]
        series.append(("DIC_total_2025", "DIC score (ISTH 2025)", "DICscore", "#4D4D4D"))
        cohort["DIC_total_2025"] = res.loc[keep, "DIC_total"].values

    for c, *_ in series:
        cohort[c] = pd.to_numeric(cohort[c], errors="coerce")

    # ---------- per-marker statistics ----------
    rows, results = [], {}
    for col, lab, short, color in series:
        sub = cohort[[col]].copy()
        ok = sub[col].notna().values
        r = analyse(cohort.loc[ok].reset_index(drop=True), col, y[ok])
        r["label"], r["color"], r["short"] = lab, color, short
        results[short] = r
        lo, hi = r["ci"]["AUC"]
        cut_disp = (int(round(float(r["cutoff"]))) if abs(r["cutoff"]) >= 100
                    else round(float(r["cutoff"]), 3))
        rows.append(dict(Biomarker=lab, N=r["n"], AUC=round(r["AUC"], 3),
                         AUC_95CI="%.3f\u2013%.3f" % (lo, hi),
                         AUPRC=round(r["AUPRC"], 3),
                         Cutoff=cut_disp,
                         Direction=r["direction"],
                         Sensitivity=round(r["Sens"], 3), Specificity=round(r["Spec"], 3),
                         PPV=round(r["PPV"], 3), NPV=round(r["NPV"], 3), F1=round(r["F1"], 3)))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT_T, stem + ".csv"), index=False, encoding="utf-8-sig")

    tab_print = tab.copy()
    tab_print["Cutoff"] = [("%d" % round(v)) if abs(v) >= 100 else ("%.3f" % v)
                           for v in tab["Cutoff"].astype(float)]

    print("=" * 70)
    print("Outcome: %s   (positive %d / negative %d)" % (out_label, int(y.sum()), int(len(y) - y.sum())))
    if args.outcome == "dic":
        print("PT control %.1f s; D-dimer ULN %.0f ug/L" % (sc.PT_CONTROL, sc.DD_ULN))
    print("=" * 70)
    print(tab_print.to_string(index=False))

    # ---------- plot (same layout as Figure 2B) ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.unicode_minus"] = False

    FS = {"tick": 10, "label": 14, "legend": 9.5, "annot": 11}
    fig, ax = plt.subplots(figsize=(7.6, 6.4))

    ax.plot([0, 1], [0, 1], color="#8A8A8A", lw=1.5, ls="--", label="No-discrimination line")
    for col, lab, short, color in series:
        r = results[short]
        lo, hi = r["ci"]["AUC"]
        ax.plot(r["fpr"], r["tpr"], color=color, lw=2.5,
                label="%s: %.3f (%.3f\u2013%.3f)" % (lab, r["AUC"], lo, hi))

    r = results["iNeu"]
    ax.text(0.97, 0.42,
            "CD10$^{-}$ neutrophil cutoff = %.1f%%\nSensitivity %.3f\nSpecificity %.3f"
            % (r["cutoff"] * 100, r["Sens"], r["Spec"]),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=FS["annot"], style="italic", color="#222222")
    ax.text(0.03, 0.97, "Outcome: %s\nn = %d (%d events)" % (out_label, len(y), int(y.sum())),
            transform=ax.transAxes, ha="left", va="top",
            fontsize=FS["annot"] - 0.5, color="#444444")

    ax.set_xlim([-0.02, 1.02]); ax.set_ylim([-0.02, 1.02])
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=FS["tick"])
    ax.tick_params(axis="y", labelsize=FS["tick"])
    ax.set_xlabel("1 - Specificity", fontsize=FS["label"], fontweight="bold")
    ax.set_ylabel("Sensitivity", fontsize=FS["label"], fontweight="bold")
    ax.legend(loc="lower right", fontsize=FS["legend"] * 0.9, frameon=True, fancybox=True,
              handlelength=1.5, labelspacing=0.32,
              facecolor="white", edgecolor="gray", borderpad=0.5)
    ax.grid(True, alpha=0.2, linestyle=":")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT_F, stem + "." + ext), dpi=300,
                    bbox_inches="tight", facecolor="white")
    print("\nwritten: %s.png / .pdf" % os.path.join(OUT_F, stem))
    return tab, results


if __name__ == "__main__":
    main()
