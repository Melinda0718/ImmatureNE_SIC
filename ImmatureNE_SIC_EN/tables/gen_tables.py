# -*- coding: utf-8 -*-
"""
gen_tables.py - single data source for the basic statistical tables
====================================================

Input : `data/NE_baseline_D1.csv` / `_D3.csv` / `_D7.csv` / `NE_HC.csv`
Output:
  1. `outputs/tables/table_*.csv`            - intermediate tables (also printed in the notebook)
  2. `tables/tex/Table1.tex` … + `.pdf`  - typeset with the manuscript template, compiled by xelatex
  3. `tables/*.docx`                       - numbers replaced in place, formatting preserved (backup in `tables/_backup/`）

Command line:
    python tables/build_tables.py             # csv + tex + pdf + docx
    python tables/build_tables.py --no-docx
    python tables/build_tables.py --no-tex
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, chi2_contingency, fisher_exact

# ------------------------------------------------------------------ paths
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
for _ in range(5):
    if os.path.isdir(os.path.join(ROOT, "data")) and os.path.isdir(os.path.join(ROOT, "outputs")):
        break
    ROOT = os.path.dirname(ROOT)
DATA = os.path.join(ROOT, "data")
OUT_T = os.path.join(ROOT, "outputs", "tables")
TEX_DIR = os.path.join(HERE, "tex")
BACKUP_DIR = os.path.join(HERE, "_backup")

# ------------------------------------------------------------------ markers
BM = [
    ("All_Neu(CD63_MFI)", "CD63 MFI in total neutrophils", "CD63"),
    ("All_Neu_Mean_(CD62L_MFI)", "CD62L MFI in total neutrophils", "CD62L"),
    ("All_Neu_Mean_(CD66b_MFI)", "CD66b MFI in total neutrophils", "CD66b"),
    ("All_Neu_Mean_(CD16_MFI)", "CD16 MFI in total neutrophils", "CD16"),
    ("immature_Neu_percentage", "CD10- neutrophils, %", "iNeu"),
]
BM_COL = [b[0] for b in BM]
BM_NAME = [b[1] for b in BM]
BM_KEY = [b[2] for b in BM]
# short names used in headers (matching the hist tables)
BM_SHORT = ["CD63 MFI", "CD62L Mean MFI", "CD66b Mean MFI", "CD16 Mean MFI",
            "CD10- neutrophil percentage"]

SEC = {
    "demo": "Demographics",
    "sev": "Disease severity",
    "com": "Comorbidity diseases",
    "lab": "Laboratory tests",
    "inf": "Infection site",
    "path": "Pathogen",
    "sup": "Support treatment",
    "out": "Clinical outcomes",
}

T1_SPEC = {
    "demo": [("cont", "Age", "Age (years)", 0),
             ("bin", "Gender", "Gender (Male, %)", "Male"),
             ("cont", "BMI", "BMI (kg/m2)", 1)],
    "sev": [("cont", "APACHEII", "APACHE II score", 0),
            ("cont", "SOFA", "SOFA score", 0),
            ("bin", "Septic Shock", "Septic shock", 1)],
    "com": [("bin", "hypertension", "Hypertension", 1),
            ("bin", "diabetes", "Diabetes mellitus", 1),
            ("bin", "chronic_heart_disease", "Chronic heart disease", 1),
            ("bin", "chronic_pulmonary_disease", "Chronic pulmonary disease", 1),
            ("bin", "chronic_kidney_disease", "Chronic kidney disease", 1),
            ("bin", "chronic_liver_disease", "Chronic liver disease", 1),
            ("bin", "cerebral_vascular_accident", "Cerebrovascular disease", 1),
            ("bin", "tumor", "Malignancy", 1),
            ("bin", "Thyroid dysfunction", "Thyroid dysfunction", 1)],
    "lab": [("cont", "WBC", "WBC (10^9/L)", 0), ("cont", "NE count", "NE count (10^9/L)", 0),
            ("cont", "LY count", "LY count (10^9/L)", 0), ("cont", "MO count", "MO count (10^9/L)", 0),
            ("cont", "PLT", "Platelet count (10^9/L)", 0), ("cont", "PT", "PT (seconds)", 1),
            ("cont", "APTT", "APTT (seconds)", 1), ("cont", "FIB", "FIB (g/L)", 1),
            ("cont", "INR", "INR", 1), ("cont", "FDP", "FDP (mg/L)", 0),
            ("cont", "D_dimer", "D-dimer (ug/L)", 0), ("cont", "PTA", "PTA (%)", 0),
            ("cont", "ATIII", "ATIII (%)", 0), ("cont", "TT", "TT (seconds)", 1),
            ("cont", "pH", "pH", 2), ("cont", "HCO3_", "Bicarbonate (mmol/L)", 0),
            ("cont", "Ca2+", "Ca2+ (mmol/L)", 2), ("cont", "lac", "Lactate (mmol/L)", 1),
            ("cont", "PCT", "PCT (ng/mL)", 0), ("cont", "CRP", "CRP (mg/L)", 0),
            ("cont", "AST", "AST (U/L)", 0), ("cont", "ALT", "ALT (U/L)", 0),
            ("cont", "DBIL", "DBIL (\u03bcmol/L)", 0), ("cont", "TBIL", "TBIL (\u03bcmol/L)", 0),
            ("cont", "ALB", "Albumin (g/L)", 0), ("cont", "BUN", "BUN (mmol/L)", 0),
            ("cont", "Scr", "Creatinine (\u03bcmol/L)", 0), ("cont", "TM", "TM (ng/mL)", 0),
            ("cont", "TAT", "TAT (ng/mL)", 0), ("cont", "t_PAIC", "t-PAIC (ng/mL)", 0),
            ("cont", "PIC", "PIC (ug/L)", 0)],
    "inf": [("bin", "Infection site_lung", "Lung", 1),
            ("bin", "Infection site_urinary", "Urinary", 1),
            ("bin", "Infection site_digestive", "Digestive", 1),
            ("bin", "Infection site_soft tissue", "Soft tissue", 1),
            ("bin", "Infection site_other", "Other", 1)],
    "path": [("bin", "G+", "G+", 1), ("bin", "G-", "G-", 1), ("bin", "Virus", "Virus", 1),
             ("bin", "fungal", "Fungal", 1),
             ("bin", "No outcome for infection", "No pathogen identified", 1)],
    "sup": [("bin", "Vasoactive drug use", "Vasoactive drug use", 1),
            ("bin", "Heparin use", "Heparin use", 1),
            ("bin", "VT_D1", "Mechanical ventilation", 1),
            ("bin", "CRRT_D1", "CRRT", 1)],
    "out": [("cont", "LOS of Hospital", "LOS of hospital (days)", 0),
            ("cont", "LOS of ICU", "LOS of ICU (days)", 0),
            ("bin", "28-day mortality", "28-day mortality", 1),
            ("bin", "ICU mortality", "ICU mortality", 1),
            ("bin", "Hospital mortality", "Hospital mortality", 1)],
}

T1_ORDER = ["demo", "sev", "com", "lab", "inf", "path", "sup", "out"]


def load_data():
    d1 = pd.read_csv(os.path.join(DATA, "NE_baseline_D1.csv"), encoding="utf-8-sig")
    d3 = pd.read_csv(os.path.join(DATA, "NE_baseline_D3.csv"), encoding="utf-8-sig")
    d7 = pd.read_csv(os.path.join(DATA, "NE_baseline_D7.csv"), encoding="utf-8-sig")
    hc = pd.read_csv(os.path.join(DATA, "NE_HC.csv"), encoding="utf-8-sig")
    return d1, d3, d7, hc


# ------------------------------------------------------------------ formatting
def fnum(v, dec):
    return "%.*f" % (dec, v)


def med_iqr(s, dec):
    q = s.quantile([.5, .25, .75]).values
    return "%s (%s, %s)" % (fnum(q[0], dec), fnum(q[1], dec), fnum(q[2], dec))


def p_fmt(p):
    return "<0.001" if p < 0.001 else "%.3f" % p


def fmt_bm(v, dec=1):
    s = "%.*f" % (dec, v)
    if dec == 1 and s.endswith(".0"):
        s = s[:-2]
    return s


def iqr_pair(s, dec, sep="-"):
    q = s.quantile([.5, .25, .75]).values
    return "%s (%s%s%s)" % (fmt_bm(q[0], dec), fmt_bm(q[1], dec), sep, fmt_bm(q[2], dec))


@dataclass
class Table:
    key: str
    caption: str
    header: list
    rows: list
    note: str = ""
    abbrev: str = ""
    kind: str = "plain"          # plain | section
    landscape: bool = False
    longtable: bool = False
    colspec: str = ""
    file_stem: str = ""

    @property
    def df(self):
        return pd.DataFrame(self.rows, columns=self.header)


# ------------------------------------------------------------------ DeLong
def _compute_midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2


def _fast_delong(preds_sorted_transposed, n_pos):
    m = n_pos; n = preds_sorted_transposed.shape[1] - m
    pos = preds_sorted_transposed[:, :m]; neg = preds_sorted_transposed[:, m:]
    k = preds_sorted_transposed.shape[0]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r, :] = _compute_midrank(pos[r, :]); ty[r, :] = _compute_midrank(neg[r, :])
        tz[r, :] = _compute_midrank(preds_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n; v10 = 1.0 - tz[:, m:] / m
    sx = np.cov(v01); sy = np.cov(v10)
    return aucs, sx / m + sy / n


def delong_test(y, s1, s2):
    """Paired DeLong test: returns (auc1, auc2, p, z)."""
    from scipy import stats as _st
    y = np.asarray(y); s1 = np.asarray(s1, float); s2 = np.asarray(s2, float)
    order = np.argsort(-y, kind="mergesort")
    preds = np.vstack((s2[order], s1[order]))
    aucs, cov = _fast_delong(preds, int(np.sum(y == 1)))
    L = np.array([[1, -1]])
    z = float((L @ aucs).item() / np.sqrt((L @ cov @ L.T).item()))
    return float(aucs[0]), float(aucs[1]), 2 * _st.norm.sf(abs(z)), z


# ================================================================== Table 1
def table1():
    d1, _, _, _ = load_data()
    cohort = d1[d1["SIC_D1"].isin(["SIC", "NSIC"])].copy()
    nsic = cohort[cohort["SIC_D1"] == "NSIC"]
    sic = cohort[cohort["SIC_D1"] == "SIC"]

    def cont(var, dec):
        a = pd.to_numeric(nsic[var], errors="coerce").dropna()
        b = pd.to_numeric(sic[var], errors="coerce").dropna()
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        return [med_iqr(a, dec), med_iqr(b, dec), p_fmt(p)]

    def binary(var, pos=1):
        num = pd.to_numeric(cohort[var], errors="coerce")
        if num.notna().any():
            vals, posval = num.fillna(0), pos
        else:
            vals, posval = cohort[var].astype(str).str.strip(), str(pos)
        ps = int((vals[sic.index] == posval).sum())
        pn = int((vals[nsic.index] == posval).sum())
        tab = np.array([[len(nsic) - pn, pn], [len(sic) - ps, ps]])
        p = fisher_exact(tab).pvalue if tab.min() < 5 else chi2_contingency(tab, correction=True).pvalue
        return ["%d (%.1f%%)" % (pn, 100 * pn / len(nsic)),
                "%d (%.1f%%)" % (ps, 100 * ps / len(sic)), p_fmt(p)]

    rows = []
    for sec in T1_ORDER:
        rows.append([SEC[sec], "", "", ""])
        for item in T1_SPEC[sec]:
            if item[0] == "cont":
                rows.append([item[2]] + cont(item[1], item[3]))
            else:
                rows.append([item[2]] + binary(item[1], item[3]))

    return Table(
        key="T1",
        caption=("Table 1. Clinical characteristics and outcomes of patients with and without "
                 "sepsis-induced coagulopathy"),
        header=["Variable", "NSIC (n=%d)" % len(nsic), "SIC (n=%d)" % len(sic), "p-value"],
        rows=rows, kind="section", longtable=True, file_stem="Table1",
        colspec=(r">{\hsize=1.6\hsize}X"
                 r">{\hsize=0.8\hsize\centering\arraybackslash}X"
                 r">{\hsize=0.8\hsize\centering\arraybackslash}X"
                 r">{\hsize=0.8\hsize\centering\arraybackslash}X"),
        note=("Data are presented as median (interquartile range) for continuous variables and "
              "n (%) for categorical variables."),
        abbrev=("SIC, sepsis-induced coagulopathy; NSIC, non-sepsis-induced coagulopathy; BMI, body mass index; "
                "APACHE II, Acute Physiology and Chronic Health Evaluation II; SOFA, Sequential Organ Failure "
                "Assessment; WBC, white blood cell count; NE, neutrophil count; LY, lymphocyte count; MO, monocyte "
                "count; PT, prothrombin time; APTT, activated partial thromboplastin time; FIB, fibrinogen; INR, "
                "international normalized ratio; FDP, fibrin degradation products; PTA, prothrombin time activity; "
                "ATIII, antithrombin III; TT, thrombin time; PCT, procalcitonin; CRP, C-reactive protein; AST, "
                "aspartate aminotransferase; ALT, alanine aminotransferase; DBIL, direct bilirubin; TBIL, total "
                "bilirubin; BUN, blood urea nitrogen; TM, thrombomodulin; TAT, thrombin-antithrombin complex; "
                "t-PAIC, tissue plasminogen activator-inhibitor complex; PIC, plasmin-alpha2-plasmin inhibitor "
                "complex; CRRT, continuous renal replacement therapy; ICU, intensive care unit; LOS, length of stay."),
    )


# ================================================================== Table 2
def table2():
    d1, _, _, hc = load_data()
    cohort = d1[d1["SIC_D1"].isin(["SIC", "NSIC"])]
    sic = cohort[cohort["SIC_D1"] == "SIC"]
    nsic = cohort[cohort["SIC_D1"] == "NSIC"]

    rows = []
    for col, name in zip(BM_COL, BM_NAME):
        is_frac = "immature" in col
        if is_frac:
            f = lambda df: iqr_pair(pd.to_numeric(df[col], errors="coerce").dropna() * 100, 2)
        else:
            f = lambda df: iqr_pair(pd.to_numeric(df[col], errors="coerce").dropna(), 1)
        a = pd.to_numeric(sic[col], errors="coerce").dropna()
        b = pd.to_numeric(nsic[col], errors="coerce").dropna()
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        rows.append([name, f(hc), f(sic), f(nsic), p_fmt(p)])

    return Table(
        key="T2",
        caption="Table 2. Neutrophil phenotypes on Day 1 according to SIC status and in healthy controls",
        header=["Neutrophil Biomarker", "HC (n=%d)" % len(hc), "SIC (n=%d)" % len(sic),
                "NSIC (n=%d)" % len(nsic), "P value"],
        rows=rows, landscape=True, colspec="l c c c c", file_stem="Table2",
        note=("Data is presented as median (interquartile range). Patients with sepsis were classified as SIC "
              "or NSIC according to their SIC status on Day 1. P values represent comparisons between the SIC "
              "and NSIC groups. CD10- neutrophils are expressed as a percentage of total neutrophils."),
        abbrev=("HC, healthy controls; SIC, sepsis-induced coagulopathy; NSIC, non-sepsis-induced coagulopathy; "
                "MFI, mean fluorescence intensity."),
    )


# ================================================================== Table 3
def table3():
    """Diagnostic-performance table (DeLong is not a row of the formatted table but is returned for CSV/notebook use)."""
    from sklearn.metrics import roc_curve, auc, precision_recall_curve

    d1, _, _, _ = load_data()
    cohort = d1[d1["SIC_D1"].isin(["SIC", "NSIC"])].copy()
    rng = np.random.default_rng(20260507)
    N_BOOT = 2000
    out = {}
    for (col, name, _), short in zip(BM, BM_KEY):
        d = cohort[[col, "SIC_D1"]].dropna()
        x = d[col].values.astype(float)
        y = (d["SIC_D1"].values == "SIC").astype(int)
        fpr, tpr, thrs = roc_curve(y, x)
        flipped = auc(fpr, tpr) < 0.5
        fpr, tpr, thrs = roc_curve(y, -x if flipped else x)
        roc_auc = auc(fpr, tpr)
        t_opt = float(thrs[int(np.argmax(tpr - fpr))])
        cut = -t_opt if flipped else t_opt
        direction = "<=" if flipped else ">="

        def rule(xv):
            return ((-xv if flipped else xv) >= t_opt).astype(int)

        pred = rule(x)
        tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
        tn = int(((pred == 0) & (y == 0)).sum()); fn = int(((pred == 0) & (y == 1)).sum())
        sens, spec = tp / (tp + fn), tn / (tn + fp)
        ppv = tp / (tp + fp) if (tp + fp) else np.nan
        npv = tn / (tn + fn) if (tn + fn) else np.nan
        f1 = 2 * ppv * sens / (ppv + sens) if (ppv + sens) else np.nan
        pr_, rc_, _ = precision_recall_curve(y, -x if flipped else x)
        auprc = auc(rc_, pr_)

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
            pb_ = rule(xb)
            tpb = int(((pb_ == 1) & (yb == 1)).sum()); fpb = int(((pb_ == 1) & (yb == 0)).sum())
            tnb = int(((pb_ == 0) & (yb == 0)).sum()); fnb = int(((pb_ == 0) & (yb == 1)).sum())
            boot["Sens"].append(tpb / (tpb + fnb) if (tpb + fnb) else np.nan)
            boot["Spec"].append(tnb / (tnb + fpb) if (tnb + fpb) else np.nan)
            boot["PPV"].append(tpb / (tpb + fpb) if (tpb + fpb) else np.nan)
            boot["NPV"].append(tnb / (tnb + fnb) if (tnb + fnb) else np.nan)
            boot["F1"].append(2 * tpb / (2 * tpb + fpb + fnb) if (2 * tpb + fpb + fnb) else np.nan)
        ci = {k: (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5)) for k, v in boot.items()}
        out[short] = dict(AUC=roc_auc, AUPRC=auprc, Sens=sens, Spec=spec, PPV=ppv, NPV=npv, F1=f1,
                          ci=ci, cutoff=cut, direction=direction, n=len(x), flipped=flipped,
                          p=None, z=None, n_delong=None)

    # paired DeLong against CD10- (CSV/notebook only, not in the formatted table)
    d_i = cohort[["immature_Neu_percentage", "SIC_D1"]].dropna()
    for (col, name, _), key in zip(BM, BM_KEY):
        if key == "iNeu":
            continue
        d = cohort[[col, "immature_Neu_percentage", "SIC_D1"]].dropna()
        yv = (d["SIC_D1"].values == "SIC").astype(int)
        s1 = d["immature_Neu_percentage"].values.astype(float)
        s2 = d[col].values.astype(float)
        if out[key]["flipped"]:
            s2 = -s2
        a1, a2, pv, zz = delong_test(yv, s1, s2)
        out[key].update(p=pv, z=zz, n_delong=len(d),
                        delong=dict(auc_ref=a1, auc_marker=a2, p=pv, z=zz, n=len(d)))

    def ci3(short, key):
        v, (lo, hi) = out[short][key], out[short]["ci"][key]
        return "%.3f (%.3f-%.3f)" % (v, lo, hi)

    def thr(v, key):
        if key == "iNeu":
            return "%.1f%%" % (v * 100 if v <= 1.5 else v)
        return ("%.1f" % v)[:-2] if abs(v - round(v)) < 1e-9 else "%.1f" % v

    rows = []
    for metric in ["AUC", "AUPRC", "Sensitivity", "Specificity", "PPV", "NPV", "F1 Score", "Optimal threshold"]:
        row = [metric]
        for s, hdr in zip(BM_KEY, BM_SHORT):
            r = out[s]
            if metric in ("AUC", "AUPRC"):
                row.append(ci3(s, metric))
            elif metric == "Optimal threshold":
                row.append({">=": "\u2265", "<=": "\u2264"}[r["direction"]] + thr(r["cutoff"], s))
            else:
                key = {"Sensitivity": "Sens", "Specificity": "Spec"}.get(metric, metric.replace(" Score", ""))
                row.append(ci3(s, key))
        rows.append(row)

    t = Table(
        key="T3",
        caption="Table 3. Diagnostic performance of Day 1 neutrophil biomarkers for SIC",
        header=["Metric"] + ["%s (n=%d)" % (h, out[k]["n"]) for k, h in zip(BM_KEY, BM_SHORT)],
        rows=rows, landscape=True, colspec="l c c c c c", file_stem="Table3",
        note=("Values are presented as point estimates with 95% confidence intervals unless otherwise indicated. "
              "Diagnostic performance was evaluated for concurrent SIC on Day 1. The 95% confidence intervals were "
              "derived from 2,000 bootstrap resamples."),
        abbrev=("SIC, sepsis-induced coagulopathy; MFI, mean fluorescence intensity; AUC, area under the receiver "
                "operating characteristic curve; AUPRC, area under the precision-recall curve; PPV, positive "
                "predictive value; NPV, negative predictive value."),
    )
    return t, out


# ================================================================== sTable S1
def stable_s1():
    d1, d3, d7, hc = load_data()
    cohort = d1[d1["SIC_D1"].isin(["SIC", "NSIC"])]
    d3c = d3[d3["SIC_D3"].isin(["SIC", "NSIC"])]
    d7c = d7[d7["SIC_D7"].isin(["SIC", "NSIC"])]

    def vals(df, gcol, col, grp):
        s = pd.to_numeric(df.loc[df[gcol] == grp, col], errors="coerce").dropna()
        if "immature" in col:
            s = s * 100
        return s

    def cell(df, gcol, col, grp):
        s = vals(df, gcol, col, grp)
        dec = 2 if "immature" in col else 1
        return iqr_pair(s, dec), len(s)

    def pval(df, gcol, col):
        a = vals(df, gcol, col, "SIC"); b = vals(df, gcol, col, "NSIC")
        return mannwhitneyu(a, b, alternative="two-sided").pvalue

    rows = []
    n3n = n3s = n7n = n7s = 0
    for col, name in zip(BM_COL, BM_NAME):
        c1n, _ = cell(cohort, "SIC_D1", col, "NSIC")
        c1s, _ = cell(cohort, "SIC_D1", col, "SIC")
        c3n, n3n = cell(d3c, "SIC_D3", col, "NSIC")
        c3s, n3s = cell(d3c, "SIC_D3", col, "SIC")
        c7n, n7n = cell(d7c, "SIC_D7", col, "NSIC")
        c7s, n7s = cell(d7c, "SIC_D7", col, "SIC")
        hc_txt, _ = cell(hc.assign(SIC_D1="SIC"), "SIC_D1", col, "SIC")
        rows.append([name, hc_txt, c1n, c1s, p_fmt(pval(cohort, "SIC_D1", col)),
                     c3n, c3s, p_fmt(pval(d3c, "SIC_D3", col)),
                     c7n, c7s, p_fmt(pval(d7c, "SIC_D7", col))])

    header = ["Neutrophil Biomarker", "HC (n=%d)" % len(hc),
              "Day 1 NSIC (n=%d)" % len(cohort[cohort.SIC_D1 == "NSIC"]),
              "Day 1 SIC (n=%d)" % len(cohort[cohort.SIC_D1 == "SIC"]), "P value",
              "Day 3 NSIC (n=%d)" % n3n, "Day 3 SIC (n=%d)" % n3s, "P value",
              "Day 7 NSIC (n=%d)" % n7n, "Day 7 SIC (n=%d)" % n7s, "P value"]
    return Table(
        key="S1",
        caption=("sTable S1. Cross-sectional distribution of neutrophil phenotypes according to SIC "
                 "status on the same day"),
        header=header, rows=rows, landscape=True, colspec="l c c c c c c c c c c", file_stem="sTableS1",
        note=("Data is presented as median (interquartile range). Patients were grouped according to their "
              "SIC status assessed on the same day (Day 1, Day 3 or Day 7). This analysis is "
              "cross-sectional: the SIC and NSIC groups at Days 3 and 7 do not contain the same patients "
              "as at Day 1, and within-patient change is not described here; the longitudinal analysis "
              "based on Day 1 status is given in sTable S2 and Supplementary Figure S3."),
        abbrev=("HC, healthy controls; SIC, sepsis-induced coagulopathy; NSIC, non-sepsis-induced coagulopathy; "
                "MFI, median fluorescence intensity."),
    )


# ================================================================== sTable S2 (LMM)
def stable_s2(return_models=False):
    import warnings
    import statsmodels.formula.api as smf
    from scipy import stats
    warnings.filterwarnings("ignore")

    d1, d3, d7, _ = load_data()
    frames = {}
    for dd, t in [(d1, 1), (d3, 3), (d7, 7)]:
        dd = dd.copy(); dd["time"] = t; dd["PatientID"] = dd["Code"].astype(str)
        frames[t] = dd
    baseline = frames[1].set_index("PatientID")["SIC_D1"].map({"SIC": 1, "NSIC": 0})
    dfl = pd.concat([frames[t] for t in (1, 3, 7)], ignore_index=True)
    dfl["SIC_Group"] = dfl["PatientID"].map(baseline)
    dfl = dfl.dropna(subset=["SIC_Group"]); dfl["SIC_Group"] = dfl["SIC_Group"].astype(int)
    dfl["time_f"] = pd.Categorical(dfl["time"].map({1: "Day 1", 3: "Day 3", 7: "Day 7"}),
                                   categories=["Day 1", "Day 3", "Day 7"], ordered=True)

    def wald(res, which):
        names = list(res.params.index)
        if which == "time":
            terms = [n for n in names if n.startswith("time_f[")]
        elif which == "group":
            terms = [n for n in names if n == "SIC_Group"]
        else:
            terms = [n for n in names if ":" in n]
        idx = [i for i, n in enumerate(res.params.index) if n in terms]
        b = np.asarray(res.params, float)[idx]
        V = np.asarray(res.cov_params(), float)[np.ix_(idx, idx)]
        W = float(b @ np.linalg.pinv(V) @ b)
        return W, float(stats.chi2.sf(W, len(idx))), len(idx)

    rows, lmm_res = [], {}
    for col, name in zip(BM_COL, BM_NAME):
        sub = dfl[["PatientID", "time_f", "SIC_Group", col]].dropna()
        model = smf.mixedlm("Q('%s') ~ time_f * SIC_Group" % col, sub, groups=sub["PatientID"])
        try:
            res = model.fit(reml=True, method="lbfgs")
        except Exception:
            res = model.fit(reml=True, method="bfgs")
        lmm_res[col] = res
        rows.append([name, sub["PatientID"].nunique(), len(sub),
                     p_fmt(wald(res, "time")[1]), p_fmt(wald(res, "group")[1]),
                     p_fmt(wald(res, "interaction")[1])])

    tab = Table(
        key="S2",
        caption="sTable S2. Mixed-effects models for longitudinal changes in neutrophil phenotypes",
        header=["Biomarker", "Patients, n", "Observations, n", "P for time", "P for Day 1 SIC group",
                "P for time \u00d7 group interaction"],
        rows=rows, colspec="l c c c c c", file_stem="sTableS2",
        note=("Linear mixed-effects models with a random intercept per patient; fixed effects were time "
              "(Day 1/3/7), Day 1 SIC group, and their interaction. P values are from joint Wald tests."),
        abbrev="SIC, sepsis-induced coagulopathy; MFI, mean fluorescence intensity.",
    )
    if return_models:
        return tab, dfl, lmm_res, wald, BM_COL
    return tab


# ================================================================== CSV
def write_csv(t: Table):
    os.makedirs(OUT_T, exist_ok=True)
    path = os.path.join(OUT_T, "table_%s.csv" % t.key.lower())
    t.df.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def write_delong(out):
    """Store DeLong tests in a separate CSV (not in docx / tex / pdf)."""
    rows = []
    for k in BM_KEY:
        r = out[k]
        if k == "iNeu":
            continue
        d = r.get("delong")
        if not d:
            continue
        rows.append({"Biomarker": k, "n (paired)": d["n"],
                     "AUC (CD10-)": round(d["auc_ref"], 4),
                     "AUC (marker)": round(d["auc_marker"], 4),
                     "DeLong Z": round(d["z"], 4), "P (vs CD10-)": round(d["p"], 4)})
    df = pd.DataFrame(rows)
    path = os.path.join(OUT_T, "table_t3_delong.csv")
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path, df


# ================================================================== LaTeX
PREAMBLE = r"""\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{times}
\usepackage{booktabs}
\usepackage{array}
\usepackage{multirow}
\usepackage{amsmath}
\usepackage{pdflscape}
\usepackage{graphicx}
\usepackage{xltabular}
\usepackage{geometry}
\usepackage{makecell}
\usepackage{caption}
\usepackage{textcomp}

\geometry{%(geometry)s}

\begin{document}

%(body)s
\end{document}
"""


def tex_label(s):
    rep = [("10^9", r"10\textsuperscript{9}"), ("m2", r"m\textsuperscript{2}"),
           ("Ca2+", r"Ca\textsuperscript{2+}"), ("CD10-", r"CD10\textsuperscript{-}"),
           ("\u03bcmol", r"\textmu mol"), ("ug/L", r"\textmu g/L"),
           ("\u2265", r"$\geq$"), ("\u2264", r"$\leq$"), ("\u00d7", r"$\times$"),
           ("%", r"\%")]
    for a, b in rep:
        s = s.replace(a, b)
    return s


def tex_cell(s):
    s = str(s)
    if "(" in s and ")" in s and "-" in s:
        head, rest = s.split("(", 1)
        body, tail = rest.rsplit(")", 1)
        body = body.replace("-", "--")
        s = head + "(" + body + ")" + tail
    return tex_label(s)


def tex_header(t):
    out = []
    for h in t.header:
        hh = tex_label(h)
        hh = hh.replace("P value", r"\textit{P} value").replace("p-value", r"\textit{p}-value")
        out.append(r"\textbf{%s}" % hh)
    return " & ".join(out) + r" \\"


def render_tex(t: Table) -> str:
    geo = ("a4paper, landscape, left=1.2cm, right=1.2cm, top=1.5cm, bottom=1.5cm, footskip=10pt"
           if t.landscape else
           "a4paper, left=1.2cm, right=1.2cm, top=0.8cm, bottom=1.2cm, footskip=10pt")
    ncol = len(t.header)
    head = tex_header(t)

    body = []
    for row in t.rows:
        if t.kind == "section" and all(str(c) == "" for c in row[1:]):
            body.append(r"\addlinespace[0.3em]")
            body.append(r"\multicolumn{%d}{l}{\textbf{\textit{%s}}} \\" % (ncol, tex_label(row[0])))
            continue
        cells = [tex_label(row[0])] + [tex_cell(c) for c in row[1:]]
        prefix = r"\quad " if t.kind == "section" else ""
        body.append(prefix + " & ".join(cells) + r" \\")
    body_txt = "\n".join(body)

    if t.longtable:
        tab = (r"\begin{xltabular}{\textwidth}{%s}" "\n" % t.colspec +
               r"\caption*{\textbf{%s}} \\" "\n" % tex_label(t.caption) +
               r"\toprule[1pt]" "\n" + head + "\n" + r"\midrule[1pt]" "\n" +
               r"\endfirsthead" "\n\n" +
               r"\caption*{\textbf{%s} (continued)} \\" "\n" % tex_label(t.caption) +
               r"\toprule[1pt]" "\n" + head + "\n" + r"\midrule[1pt]" "\n" +
               r"\endhead" "\n\n" + r"\midrule" "\n" + r"\endfoot" "\n\n" +
               r"\endlastfoot" "\n\n" + body_txt + "\n" + r"\bottomrule[1pt]" "\n" +
               r"\end{xltabular}" "\n")
    else:
        tab = (r"\begin{tabular}{%s}" "\n" % t.colspec + r"\toprule[1pt]" "\n" + head + "\n" +
               r"\midrule[1pt]" "\n" + body_txt + "\n" + r"\bottomrule[1pt]" "\n" +
               r"\end{tabular}" "\n")
        if t.landscape and ncol >= 10:
            tab = (r"\resizebox{\dimexpr\textwidth+1.2cm\relax}{!}{%" "\n" + tab + "}\n")

    note_txt = ("\\vspace{2mm}\n\\begin{minipage}{\\textwidth}\n\\footnotesize\n"
                "\\textbf{Note:} " + tex_label(t.note) + "\n\n"
                "\\textbf{Abbreviations:} " + tex_label(t.abbrev) + "\n"
                "\\end{minipage}\n")

    if t.longtable:
        block = "\\begin{center}\n" + tab + "\n" + note_txt + "\\end{center}\n"
    else:
        block = ("\\begin{center}\n"
                 + ("\\vspace*{\\fill}\n" if t.landscape else "")
                 + "\\captionof*{table}{\\textbf{%s}}" % tex_label(t.caption) + "\n\n"
                 "\\renewcommand{\\arraystretch}{1.45}\n"
                 "\\setlength{\\tabcolsep}{4pt}\n"
                 "\\small\n\n" + tab + "\n" + note_txt
                 + ("\\vspace*{\\fill}\n" if t.landscape else "")
                 + "\\end{center}\n")

    if t.landscape:
        block = "\\clearpage\n" + block
    return PREAMBLE % {"geometry": geo, "body": block}


def compile_tex(tex_path):
    cwd = os.path.dirname(tex_path)
    name = os.path.basename(tex_path)
    for _ in range(2):
        subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error", name],
                       cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pdf = os.path.splitext(tex_path)[0] + ".pdf"
    return pdf if os.path.exists(pdf) else None


def write_tex(t: Table):
    os.makedirs(TEX_DIR, exist_ok=True)
    tex = os.path.join(TEX_DIR, t.file_stem + ".tex")
    with open(tex, "w", encoding="utf-8") as fh:
        fh.write(render_tex(t))
    pdf = compile_tex(tex)
    return tex, pdf


def render_supp_tex(tables):
    """Combined sTable rendering: several tables on one landscape page."""
    geo = "a4paper, landscape, left=1.0cm, right=1.0cm, top=1.0cm, bottom=1.0cm, footskip=10pt"
    blocks = []
    for i, t in enumerate(tables):
        ncol = len(t.header)
        head = tex_header(t)
        body = []
        for row in t.rows:
            cells = [tex_label(row[0])] + [tex_cell(c) for c in row[1:]]
            body.append(" & ".join(cells) + r" \\")
        tab_inner = (r"\begin{tabular}{%s}" "\n" % t.colspec + r"\toprule[1pt]" "\n" + head + "\n"
                     + r"\midrule[1pt]" "\n" + "\n".join(body) + "\n" + r"\bottomrule[1pt]" "\n"
                     + r"\end{tabular}" "\n")
        if ncol >= 10:
            tab_inner = (r"\resizebox{\dimexpr\textwidth\relax}{!}{%" "\n" + tab_inner + "}\n")
        blk = ("\n\\vspace{%s}\n" % ("0pt" if i == 0 else "10mm")
               + "\\begin{center}\n"
               + "\\captionof*{table}{\\textbf{%s}}" % tex_label(t.caption) + "\n\n"
               + "\\renewcommand{\\arraystretch}{1.35}\n"
               + "\\setlength{\\tabcolsep}{3pt}\n"
               + "\\footnotesize\n\n" + tab_inner + "\n"
               + "\\vspace{1.5mm}\n\\begin{minipage}{0.98\\textwidth}\n\\scriptsize\n"
               + "\\textbf{Note:} " + tex_label(t.note) + "\n\n"
               + "\\textbf{Abbreviations:} " + tex_label(t.abbrev) + "\n"
               + "\\end{minipage}\n\\end{center}\n")
        blocks.append(blk)
    return PREAMBLE % {"geometry": geo, "body": "\n".join(blocks)}


def write_supp_tex(tables, stem="sTables_S1_S2"):
    os.makedirs(TEX_DIR, exist_ok=True)
    tex = os.path.join(TEX_DIR, stem + ".tex")
    with open(tex, "w", encoding="utf-8") as fh:
        fh.write(render_supp_tex(tables))
    pdf = compile_tex(tex)
    return tex, pdf


# ================================================================== docx
def _norm(s):
    """Normalise row labels: superscript minus/square and redundant spaces."""
    return (str(s).replace("\u207b", "-").replace("\u00b2", "2").replace("\u00b9", "1")
            .replace("\n", " ").strip())


def _set_cell_text(cell, text):
    paras = cell.paragraphs
    if not paras:
        cell.text = text
        return
    p = paras[0]
    if p.runs:
        p.runs[0].text = text
        for r in p.runs[1:]:
            r.text = ""
    else:
        p.add_run(text)
    for extra in paras[1:]:
        for r in extra.runs:
            r.text = ""


def update_docx(t: Table, docx_path: str, table_index: int = 0):
    import docx
    doc = docx.Document(docx_path)
    tbl = doc.tables[table_index]
    rows = {_norm(r.cells[0].text): r for r in tbl.rows}
    n = 0
    for row in t.rows:
        if all(str(c) == "" for c in row[1:]):
            continue
        r = rows.get(_norm(row[0]))
        if r is None:
            continue
        for j, val in enumerate(row[1:], start=1):
            if j < len(r.cells):
                _set_cell_text(r.cells[j], str(val))
                n += 1
    doc.save(docx_path)
    return n


def update_word(tables):
    os.makedirs(BACKUP_DIR, exist_ok=True)
    jobs = [("Table 1.docx", "T1", 0), ("Table 2.docx", "T2", 0), ("Table 3.docx", "T3", 0),
            ("supplementary_tables.docx", "S1", 0), ("supplementary_tables.docx", "S2", 1)]
    by_key = {t.key: t for t in tables}
    backed = set()
    for fname, key, idx in jobs:
        path = os.path.join(HERE, fname)
        if not os.path.exists(path) or key not in by_key:
            continue
        if fname not in backed:
            shutil.copy2(path, os.path.join(BACKUP_DIR, fname))
            backed.add(fname)
        try:
            n = update_docx(by_key[key], path, table_index=idx)
            print("  docx %-22s %-3s updated %d cells" % (fname, key, n))
        except PermissionError:
            alt = os.path.splitext(path)[0] + "_updated.docx"
            if not os.path.exists(alt):
                shutil.copy2(path, alt)
            n = update_docx(by_key[key], alt, table_index=idx)
            print("  docx %-22s %-3s file locked -> saved as %s (%d cells)"
                  % (fname, key, os.path.basename(alt), n))


# ================================================================== entry point
def build_all(do_tex=True, do_docx=True):
    t1 = table1()
    t2 = table2()
    t3, t3_out = table3()
    s1 = stable_s1()
    s2 = stable_s2()
    tables = [t1, t2, t3, s1, s2]
    write_delong(t3_out)
    for t in tables:
        write_csv(t)
    if do_tex:
        for t in (t1, t2, t3):
            tex, pdf = write_tex(t)
            print("  tex  %-3s %s%s" % (t.key, os.path.relpath(tex, ROOT), "" if pdf else " [PDF failed]"))
        tex, pdf = write_supp_tex([s1, s2])
        print("  tex  S*  %s%s" % (os.path.relpath(tex, ROOT), "" if pdf else " [PDF failed]"))
    print("  csv  -> outputs/tables/table_*.csv")
    if do_docx:
        update_word(tables)
    return tables


if __name__ == "__main__":
    build_all(do_tex="--no-tex" not in sys.argv, do_docx="--no-docx" not in sys.argv)
