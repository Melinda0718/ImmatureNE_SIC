# -*- coding: utf-8 -*-
"""
01_dic2025_scoring.py -- recompute Day-1 overt DIC with the ISTH 2025 criteria
================================================================
Input : data/NE_baseline_D1.csv (99 patients)
Output: results/tables/DIC2025_D1_scores.csv   per-patient subscores, total and call
      results/tables/DIC2025_summary.csv     case count, score distribution, cross-tabs
      results/tables/DIC2025_sensitivity.csv sensitivity analysis (D-dimer ULN, PT control)

ISTH 2025 score (total >= 5 in a patient with a predisposing condition supports overt DIC)
    platelet count (x10^9/L)  >=100 -> 0; >=50 and <100 -> 1; <50 -> 2
    D-dimer / local upper limit of normal  <=3x -> 0; >3x and <=7x -> 2; >7x -> 3
    PT prolongation (seconds above control)  <3 -> 0; >=3 and <6 -> 1; >=6 -> 2
    fibrinogen (g/L)  >=1.0 -> 0; <1.0 -> 1
Two points where this implementation differs from common practice:
    - PT is scored by seconds of prolongation, not by INR thresholds;
    - D-dimer is scaled to the local assay upper limit of normal: exactly 3x scores 0,
      exactly 7x scores 2, and only above 7x does it score 3.

Assay parameters used for this cohort (see README_DIC2025.md for their derivation):
    PT control = 11.0 s     -- back-calculated from the INR values: 10.89 s (R^2 = 0.974),
                                  and 10.59 s from PTA; 11.0 s is used as the round value
    D-dimer ULN = 500 ug/L  -- the LIS report gives <=500 ug/L for the D-dimer assay
                                 (`D_dimer` in the data file is in ug/L)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ paths
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
for _ in range(6):
    if os.path.isdir(os.path.join(ROOT, "data")) and os.path.isdir(os.path.join(ROOT, "code")):
        break
    ROOT = os.path.dirname(ROOT)
DATA = os.path.join(ROOT, "data")
OUT_T = os.path.join(ROOT, "results", "tables")
os.makedirs(OUT_T, exist_ok=True)

# ------------------------------------------------------------------ parameters
PT_CONTROL = 11.0        # PT control in seconds; LIS reference 9.4-12.5 (midpoint 10.95)
DD_ULN = 500.0           # D-dimer upper limit of normal (ug/L)
CUT = 5                  # ISTH 2025: total >= 5 supports overt DIC

# ------------------------------------------------------------------ Day-2 backfill
# Two patients are missing the whole Day-1 coagulation panel. Day-2 fibrinolysis
# results are used for DIC scoring only; the data file is not modified and platelets stay Day 1.
# Source: ICU laboratory system screenshots of the routine fibrinolysis panel.
D2_BACKFILL = {
    "ZD-SIC-124": dict(date="2024-12-29", PT=10.6, INR=0.95, PTA=110.0, FIB=4.9,
                       APTT=35.1, TT=15.5, ATIII=55.0, FDP=4.30, D_dimer=1144.0),
    "ZD-SIC-135": dict(date="2025-03-05", PT=12.2, INR=1.10, PTA=87.0, FIB=4.0,
                       APTT=36.2, TT=19.7, ATIII=86.0, FDP=42.80, D_dimer=18127.0),
}
FILL_COLS = ["PT", "INR", "PTA", "FIB", "APTT", "TT", "ATIII", "FDP", "D_dimer"]


# ------------------------------------------------------------------ scoring
def pts_plt(v):
    """Platelet count (x10^9/L): >=100 -> 0, >=50 and <100 -> 1, <50 -> 2."""
    if pd.isna(v):
        return np.nan
    if v >= 100:
        return 0
    if v >= 50:
        return 1
    return 2


def pts_ddimer(v, uln=DD_ULN):
    """D-dimer / upper limit of normal: <=3x -> 0, >3x and <=7x -> 2, >7x -> 3."""
    if pd.isna(v):
        return np.nan
    r = v / uln
    if r <= 3:
        return 0
    if r <= 7:
        return 2
    return 3


def pts_pt(v, control=PT_CONTROL):
    """PT prolongation in seconds: <3 -> 0, >=3 and <6 -> 1, >=6 -> 2."""
    if pd.isna(v):
        return np.nan
    d = v - control
    if d < 3:
        return 0
    if d < 6:
        return 1
    return 2


def pts_fib(v):
    """Fibrinogen (g/L): >=1.0 -> 0, <1.0 -> 1."""
    if pd.isna(v):
        return np.nan
    return 0 if v >= 1.0 else 1


MAXPTS = {"PLT_pts": 2, "Ddimer_pts": 3, "PT_pts": 2, "FIB_pts": 1}


def score_frame(d1, pt_control=PT_CONTROL, dd_uln=DD_ULN):
    """Return the per-patient score table (including the missing-value columns)."""
    d = d1.copy()
    for c in ["PLT", "PT", "FIB", "D_dimer"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    out = pd.DataFrame({"Code": d["Code"]})
    for c in ["PLT", "PT", "FIB", "D_dimer"]:
        out[c] = d[c]
    out["PT_prolong_s"] = d["PT"] - pt_control
    out["Ddimer_xULN"] = d["D_dimer"] / dd_uln

    out["PLT_pts"] = d["PLT"].map(pts_plt)
    out["Ddimer_pts"] = d["D_dimer"].map(lambda v: pts_ddimer(v, dd_uln))
    out["PT_pts"] = d["PT"].map(lambda v: pts_pt(v, pt_control))
    out["FIB_pts"] = d["FIB"].map(pts_fib)

    sub = out[["PLT_pts", "Ddimer_pts", "PT_pts", "FIB_pts"]]
    out["DIC_total"] = sub.sum(axis=1, min_count=1)          # sum of available items
    out["n_missing"] = sub.isna().sum(axis=1)

    # upper bound assuming the lowest score for missing items
    ceil = sub.fillna(0).sum(axis=1)
    for k, mx in MAXPTS.items():
        ceil = ceil + np.where(out[k].isna(), mx, 0)
    out["max_possible"] = ceil

    # call: available items already >=5 -> yes; even the maximum with missing items <5 -> no; else indeterminate
    out["DIC_2025"] = np.select(
        [out["DIC_total"] >= CUT, out["max_possible"] < CUT],
        [1.0, 0.0], default=np.nan)
    out["evaluable"] = out["DIC_2025"].notna()

    out["SIC_D1"] = d["SIC_D1"]
    out["DIC_old"] = pd.to_numeric(d["DIC"], errors="coerce")
    return out


def apply_backfill(d1, backfill=None):
    """Fill Day-2 values into the empty Day-1 slots (never overwrites existing values).

    Returns (filled copy, coagulation source column). The data file is untouched.
    """
    backfill = D2_BACKFILL if backfill is None else backfill
    d = d1.copy()
    src = pd.Series("D1", index=d.index, name="coag_source")
    for code, vals in backfill.items():
        m = d["Code"] == code
        if not m.any():
            continue
        i = d.index[m][0]
        filled = []
        for c in FILL_COLS:
            if pd.isna(pd.to_numeric(d.at[i, c], errors="coerce")) and c in vals:
                d.at[i, c] = vals[c]
                filled.append(c)
        if filled:
            src.at[i] = f"Day-2 backfill ({vals['date']})"
    return d, src


# ------------------------------------------------------------------ main
def main():
    d1_raw = pd.read_csv(os.path.join(DATA, "NE_baseline_D1.csv"), encoding="utf-8-sig")
    d1, coag_src = apply_backfill(d1_raw)
    res = score_frame(d1)
    res["coag_source"] = coag_src.values

    cols = ["Code", "PLT", "PLT_pts", "PT", "PT_prolong_s", "PT_pts",
            "D_dimer", "Ddimer_xULN", "Ddimer_pts", "FIB", "FIB_pts",
            "DIC_total", "n_missing", "max_possible", "DIC_2025",
            "coag_source", "SIC_D1", "DIC_old"]
    res[cols].to_csv(os.path.join(OUT_T, "DIC2025_D1_scores.csv"),
                     index=False, encoding="utf-8-sig")

    n = len(res)
    n_pos = int((res["DIC_2025"] == 1).sum())
    n_neg = int((res["DIC_2025"] == 0).sum())
    n_ind = int(res["DIC_2025"].isna().sum())

    print("=" * 62)
    print("ISTH 2025 overt DIC - Day-1 recalculation")
    print("=" * 62)
    print(f"PT control = {PT_CONTROL:.1f} s; D-dimer ULN = {DD_ULN:.0f} ug/L")
    print(f"Cohort n = {n}; assessable {n_pos + n_neg}")
    print(f"Overt DIC (total >= {CUT}) = {n_pos} ({100 * n_pos / n:.1f}%)")
    print(f"No DIC = {n_neg}; indeterminate = {n_ind}")
    if n_ind:
        print("  indeterminate: " + ", ".join(res.loc[res["DIC_2025"].isna(), "Code"]))

    filled = res[res["coag_source"] != "D1"]
    if len(filled):
        print()
        print(f"[backfill] {len(filled)} patient(s) with a missing Day-1 panel:")
        for _, r in filled.iterrows():
            print("  %s  %s" % (r["Code"], r["coag_source"]))
            print("     PLT %.0f -> %.0f pts | PT %.1f s (%+.1f s) -> %.0f pts | FIB %.1f -> %.0f pts | D-dimer %.0f ug/L = %.2fxULN -> %.0f pts"
                  % (r["PLT"], r["PLT_pts"], r["PT"], r["PT_prolong_s"], r["PT_pts"],
                     r["FIB"], r["FIB_pts"], r["D_dimer"], r["Ddimer_xULN"], r["Ddimer_pts"]))
            print("     total %.0f -> %s" % (r["DIC_total"], "DIC" if r["DIC_2025"] == 1 else "no DIC"))

    print("\nTotal score distribution:")
    dist = res["DIC_total"].value_counts().sort_index()
    for k, v in dist.items():
        print(f"  {int(k) if float(k).is_integer() else k:>2} points: {v:>2} patients")

    # compare with the DIC column computed with ULN = 1000
    old_pos = res["DIC_old"] >= CUT
    both = res["DIC_2025"].notna()
    agree = int((res.loc[both, "DIC_2025"].astype(bool) == old_pos[both]).sum())
    print(f"\nAgreement with the DIC column computed at ULN = 1000: {agree}/{int(both.sum())}"
          f"previous positives {int(old_pos.sum())} -> now {n_pos}")
    chg = res.loc[both & (res["DIC_2025"].astype(bool) != old_pos)]
    if len(chg):
        print("  newly classified as DIC (previous total < 5):")
        for _, r in chg.iterrows():
            print("    %s: total %.0f (previous %s), D-dimer %.0f ug/L = %.2fxULN"
                  % (r["Code"], r["DIC_total"], int(r["DIC_old"]), r["D_dimer"], r["Ddimer_xULN"]))

    print("\nCross-tabulation with Day-1 SIC:")
    ct = pd.crosstab(res["SIC_D1"], res["DIC_2025"].fillna(-1),
                     rownames=["SIC_D1"], colnames=["DIC_2025"])
    print(ct.to_string().replace("-1.0", "indeterminate"))

    # summary table
    summ = pd.DataFrame([
        dict(item="n", value=n),
        dict(item="PT_control_s", value=PT_CONTROL),
        dict(item="D_dimer_ULN_ug_per_L", value=DD_ULN),
        dict(item="DIC_cutoff", value=CUT),
        dict(item="DIC_positive", value=n_pos),
        dict(item="DIC_negative", value=n_neg),
        dict(item="DIC_indeterminate", value=n_ind),
        dict(item="n_D2_backfilled", value=int((res["coag_source"] != "D1").sum())),
        dict(item="DIC_prevalence_pct", value=round(100 * n_pos / n, 1)),
        dict(item="median_score", value=float(res["DIC_total"].median())),
        dict(item="mean_score", value=round(float(res["DIC_total"].mean()), 2)),
    ])
    summ.to_csv(os.path.join(OUT_T, "DIC2025_summary.csv"), index=False, encoding="utf-8-sig")

    # parameter sensitivity
    rows = []
    for uln in [500.0, 800.0, 1000.0, 2000.0]:
        for ctrl in [10.9, 11.0, 12.0, 12.9]:
            r = score_frame(d1, pt_control=ctrl, dd_uln=uln)
            rows.append(dict(DD_ULN=uln, PT_control=ctrl,
                             n_DIC=int((r["DIC_2025"] == 1).sum()),
                             n_evaluable=int(r["DIC_2025"].notna().sum()),
                             prevalence_pct=round(100 * (r["DIC_2025"] == 1).sum() / n, 1)))
    sens = pd.DataFrame(rows)
    sens.to_csv(os.path.join(OUT_T, "DIC2025_sensitivity.csv"), index=False, encoding="utf-8-sig")
    print("\nParameter sensitivity (DIC cases):")
    print(sens.pivot(index="PT_control", columns="DD_ULN", values="n_DIC").to_string())
    return res


if __name__ == "__main__":
    main()
