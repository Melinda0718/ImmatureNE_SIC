# -*- coding: utf-8 -*-
"""
gen_pred_tex.py - LaTeX / PDF for the prediction-model supplementary tables (sTable S3-S7)
=================================================================
Reads `outputs/tables/nested_cv_results.pkl` and `nested_cv_stage2.pkl`,
typeset with the manuscript template (xltabular + booktabs + caption* + Note/Abbreviations)
into one landscape document, compiled twice with `xelatex` -> `tables/tex/sTables_S3_S7.tex/.pdf`.

    python code/prediction_model/layout/gen_pred_tex.py
"""
import os
import pickle
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
for _ in range(6):
    if os.path.isdir(os.path.join(ROOT, "data")) and os.path.isdir(os.path.join(ROOT, "outputs")):
        break
    ROOT = os.path.dirname(ROOT)

sys.path.insert(0, os.path.join(ROOT, "tables"))
import gen_tables as gt  # noqa: E402  reuse preamble / escaping / compilation

OUT_T = os.path.join(ROOT, "outputs", "tables")
TEX = os.path.join(ROOT, "tables", "tex")

R = pickle.load(open(os.path.join(OUT_T, "nested_cv_results.pkl"), "rb"))
S2 = pickle.load(open(os.path.join(OUT_T, "nested_cv_stage2.pkl"), "rb"))
summ = R["summary"]

PARS = S2["pars_features"]
N_PARS = len(PARS)
PARS_HEAD = "Simplified model (%d predictors)" % N_PARS
PARS_SEC = "Simplified model (four-predictor specification, secondary analysis)"


def f3(v):
    return "%.3f (%.3f-%.3f)" % (v[0], v[1], v[2])


def fp(p):
    return "<0.001" if p < 0.001 else "%.3f" % p


def fp4(p):
    return "<0.0001" if p < 0.0001 else "%.4f" % p


def tab(key, caption, header, rows, note, colspec=None):
    return gt.Table(key=key, caption=caption, header=header, rows=rows,
                    note=note, abbrev="", colspec=colspec or ("l" + " c" * (len(header) - 1)),
                    file_stem=key)


# ---------------------------------------------------------------- S3
algo = S2["algo_res"]
rows3 = [[n, f3((r["pooled_point"]["auc"], *r["boot_ci"]["auc"])),
          f3((r["pooled_point"]["auprc"], *r["boot_ci"]["auprc"]))] for n, r in algo.items()]
S3 = tab("sTableS6", "sTable S6. Cross-validated performance of four candidate algorithms",
         ["Algorithm", "AUROC (95% CI)", "AUPRC (95% CI)"], rows3,
         "Algorithms were compared using the same outer splits; only hyperparameters were tuned "
         "within each outer training fold. Point estimates are from pooled out-of-fold predictions; "
         "95% confidence intervals are from 2,000 patient-level bootstrap resamples.")

# ---------------------------------------------------------------- S4
desc = S2["desc"]
KEYS = ["Clinical_Only", "Single_Immature", "Clinical_Immature", "Parsimonious_Pipeline"]
HEAD4 = ["Clinical model (7 predictors)", "CD10- neutrophils alone",
         "Combined model (8 predictors)", PARS_SEC]
rows4 = []
for m, lab in [("auc", "AUROC"), ("auprc", "AUPRC"), ("sensitivity", "Sensitivity"),
               ("specificity", "Specificity"), ("ppv", "PPV"), ("npv", "NPV"), ("f1", "F1 score")]:
    rows4.append([lab] + [f3((summ[k]["pooled_point"][m], *summ[k]["boot_ci"][m]))
                          for k in KEYS])
_c = S2["cmp_rows"]
rows4.append(["Delta AUROC"] + ["\u2014", "\u2014"] + [f3(_c[i]["d_auc"]) for i in (0, 1)])
rows4.append(["Continuous NRI"] + ["\u2014", "\u2014"] + [f3(_c[i]["nri"]) for i in (0, 1)])
rows4.append(["IDI"] + ["\u2014", "\u2014"] + [f3(_c[i]["idi"]) for i in (0, 1)])
rows4.append(["P for Delta AUROC"] + ["\u2014", "\u2014"] + [fp4(_c[i]["p_d_auc"]) for i in (0, 1)])
S4 = tab("sTableS3", "sTable S3. Predictive performance of the clinical, biomarker, combined, and simplified models",
         ["Metric"] + HEAD4, rows4,
         "All metrics were calculated on pooled patient-level out-of-fold predictions, the same "
         "probabilities used in the ROC, calibration and decision-curve analyses. The four "
         "predictors were those most frequently retained across the 100 outer folds; because this "
         "set was defined in the same cohort, the estimate should be read as the performance of "
         "the full modelling strategy (0.938 when predictors are re-selected within every fold). "
         "Calibration slope and intercept for the simplified model were 2.15 and 0.31 "
         "(unpenalised logistic recalibration); a slope substantially above 1 indicates that "
         "predicted probabilities were under-confident. Delta "
         "values, continuous NRI and IDI are relative to the clinical model and were obtained from "
         "2,000 paired patient-level bootstrap resamples; NRI and IDI are exploratory and unstable "
         "at this sample size.",
         colspec="l" + " c" * 4)

# ---------------------------------------------------------------- S5
rows5 = []
from statsmodels.stats.multitest import multipletests as _mt
_q3 = _mt([r["p_d_auc"] for r in S2["cmp_rows"]], method="fdr_bh")[1]
for _ri, r in enumerate(S2["cmp_rows"]):
    def _p(tag):
        if tag == "d_auc":
            return "q %s (P %s)" % (fp4(_q3[_ri]), fp4(r["p_" + tag]))
        return "P %s" % fp4(r["p_" + tag])
    _lab = r["comparison"].replace("Simplified model vs",
                                   "Simplified model (secondary analysis) vs")
    rows5.append([_lab,
                  "%s; %s" % (f3(r["d_auc"]), _p("d_auc")),
                  "%s; %s" % (f3(r["d_auprc"]), _p("d_auprc")),
                  "%s; %s" % (f3(r["d_brier"]), _p("d_brier")),
                  "%s; %s" % (f3(r["nri"]), _p("nri")),
                  "%s; %s" % (f3(r["idi"]), _p("idi"))])
S5 = tab("sTableS4", "sTable S4. Pairwise comparison of the four prediction models",
         ["Comparison", "Delta AUROC (95% CI)", "Delta AUPRC (95% CI)",
          "Delta Brier score (95% CI)", "Continuous NRI (95% CI)", "IDI (95% CI)"], rows5,
         "Contrasts were computed on pooled patient-level predictions; 95% confidence intervals "
         "were derived from 2,000 paired patient-level bootstrap resamples in which both models "
         "were scored on the same resampled patients. Positive Delta AUROC, Delta AUPRC, NRI and "
         "IDI favour the first-listed model; negative Delta Brier indicates better calibration. "
         "Fifteen tests were performed in total (five contrasts for each of three pairwise "
         "comparisons). Only Delta AUROC was formally tested, across the three contrasts, and the "
         "q values are Benjamini-Hochberg adjusted across these three tests. Delta AUPRC, "
         "Delta Brier score, continuous NRI and IDI are reported descriptively, without "
         "multiplicity adjustment, and are exploratory and unstable at this sample size.")

# ---------------------------------------------------------------- S6
mk = S2["markers"]
rows6 = [[n, f3(d["auc"]), f3(d["auprc"])] for n, d in mk.items()]
for k, lab in [("Clinical_Only", "Clinical model (7 predictors)"),
               ("Clinical_Immature", "Combined model (8 predictors)"),
               ("Parsimonious_Pipeline", PARS_SEC)]:
    s = summ[k]
    rows6.append([lab, f3((s["pooled_point"]["auc"], *s["boot_ci"]["auc"])),
                  f3((s["pooled_point"]["auprc"], *s["boot_ci"]["auprc"]))])
S6 = tab("sTableS5", "sTable S5. Discrimination of individual predictors and multivariable models",
         ["Predictor or model", "AUROC (95% CI)", "AUPRC (95% CI)"], rows6,
         "Individual predictors were evaluated on their raw Day 1 values using the same outer test "
         "folds; multivariable models are the nested cross-validated pooled out-of-fold estimates.")

# ---------------------------------------------------------------- S7 (descriptive)
sg = S2["subgroup"]
_models = sg["models"]
def _iqr(d, k):
    return "%.3f (%.3f-%.3f)" % (d[k], *d[k.replace("median", "iqr")])
rows7 = [["Patients, n (events)"] + ["%d (%d)" % (sg["n"], sg["events"])] * len(_models),
         ["Predicted risk, median (IQR) - events"] + [_iqr(d, "median_event") for d in _models.values()],
         ["Predicted risk, median (IQR) - no events"] + [_iqr(d, "median_nonevent") for d in _models.values()]]
S7 = tab("sTableS7", "sTable S7. Descriptive analysis of the subgroup without SIC on Day 1",
         ["Metric"] + [n for n in _models.keys()], rows7,
         "Descriptive analysis only; no discrimination statistics were estimated because of the "
         "very small number of incident events. Predicted risks are nested cross-validated pooled "
         "out-of-fold probabilities, summarised as median and interquartile range.",
         colspec="l" + " c" * len(_models))

tables = [S4, S5, S6, S3, S7]   # narrative order matching the Word document

# ---------------------------------------------------------------- rendering (one table per page)
geo = "a4paper, landscape, left=1.0cm, right=1.0cm, top=1.0cm, bottom=1.0cm, footskip=10pt"
blocks = []
for i, t in enumerate(tables):
    head = gt.tex_header(t)
    body = []
    for row in t.rows:
        cells = [gt.tex_label(row[0])] + [gt.tex_cell(c) for c in row[1:]]
        body.append(" & ".join(cells) + r" \\")
    tabular = (r"\begin{tabular}{%s}" "\n" % t.colspec + r"\toprule[1pt]" "\n" + head + "\n"
               + r"\midrule[1pt]" "\n" + "\n".join(body) + "\n" + r"\bottomrule[1pt]" "\n"
               + r"\end{tabular}" "\n")
    if len(t.header) >= 5:
        tabular = (r"\resizebox{\dimexpr\textwidth\relax}{!}{%" "\n" + tabular + "}\n")
    head_blocks = ("\n\\clearpage\n" if i else "\n") + "\\begin{center}\n"
    blocks.append(head_blocks
                  + "\\captionof*{table}{\\textbf{%s}}" % gt.tex_label(t.caption) + "\n\n"
                  + "\\renewcommand{\\arraystretch}{1.35}\n\\setlength{\\tabcolsep}{4pt}\n\\small\n\n"
                  + tabular + "\n"
                  + "\\vspace{1.5mm}\n\\begin{minipage}{0.98\\textwidth}\n\\footnotesize\n"
                  + "\\textbf{Note:} " + gt.tex_label(t.note) + "\n"
                  + "\\end{minipage}\n\\end{center}\n")

tex = gt.PREAMBLE % {"geometry": geo, "body": "".join(blocks)}
os.makedirs(TEX, exist_ok=True)
path = os.path.join(TEX, "sTables_S3_S7.tex")
open(path, "w", encoding="utf-8").write(tex)
pdf = gt.compile_tex(path)
print("written", path, "| pdf:", bool(pdf))
