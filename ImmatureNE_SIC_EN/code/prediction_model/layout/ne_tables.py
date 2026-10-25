# -*- coding: utf-8 -*-
"""
ne_tables.py  (table builder)
================================================================================
Updates the prediction-model tables in `tables/supplementary_tables.docx` in place
(sTable S3-S7) and their notes: formatting preserved, values replaced by the nested-CV results.
Also updates `figures/Figure legend.docx` (Figures 3/4) and `figures/Supplementary Figure.docx`
(S6/S7) methodological wording.
================================================================================
"""
import os
import sys
import copy
import pickle
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import docx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
for _ in range(5):
    if os.path.isdir(os.path.join(ROOT, 'data')) and os.path.isdir(os.path.join(ROOT, 'outputs')):
        break
    ROOT = os.path.dirname(ROOT)
OUT_T = os.path.join(ROOT, 'outputs', 'tables')
TAB = os.path.join(ROOT, 'tables')
FIG = os.path.join(ROOT, 'figures')
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

with open(os.path.join(OUT_T, 'nested_cv_results.pkl'), 'rb') as f:
    R = pickle.load(f)
with open(os.path.join(OUT_T, 'nested_cv_stage2.pkl'), 'rb') as f:
    S2 = pickle.load(f)

summ = R['summary']
pars_features = S2['pars_features']
N_PARS = len(pars_features)
PARS_HEAD = f'Simplified model ({N_PARS} predictors)'
PARS_SEC = ('Simplified model (four-predictor specification, '
            'secondary analysis)')
pars_display = [ {'INR': 'INR', 'PLT': 'platelet count', 'lac': 'lactate',
                  'SOFA D1 Circulation': 'SOFA cardiovascular subscore',
                  'SOFA D1 Respiratory System': 'SOFA respiratory subscore',
                  'VT_D1': 'mechanical ventilation status',
                  'SOFA_Creatinine': 'SOFA renal subscore',
                  'immature_Neu_percentage': 'CD10\u207b neutrophil percentage'}[f]
                 for f in pars_features ]
# predictor order in the table notes: CD10 first
pars_display_note = ['CD10\u207b neutrophil percentage'] + \
    [d for d in pars_display if not d.startswith('CD10')]


def r3(x):
    return np.floor(float(x) * 1000 + 0.5 + 1e-9) / 1000


def fnum(x, plus=False):
    s = f'{r3(x):.3f}'
    return ('+' + s) if (plus and x >= 0) else s


def fci(m, lo, hi, plus=False):
    return f'{fnum(m, plus)} ({fnum(lo, plus)}, {fnum(hi, plus)})'


def fp(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return '-'
    return '<0.0001' if p < 0.0001 else f'{p:.4f}'.rstrip('0').rstrip('.') if p >= 0.1 else f'{p:.4f}'


def fp_cmp(p):
    if p < 0.0001:
        return 'P < 0.0001'
    return f'P = {p:.4f}'.replace('0.', '0.', 1)


def set_cell(cell, text):
    """Replace cell text while keeping the formatting of the first run."""
    para = cell.paragraphs[0]
    # drop extra paragraphs
    for p in cell.paragraphs[1:]:
        p._element.getparent().remove(p._element)
    runs = para.runs
    if not runs:
        para.add_run(text)
        return
    runs[0].text = text
    for r in runs[1:]:
        r.text = ''


def set_para(para, text):
    runs = para.runs
    if not runs:
        para.add_run(text)
        return
    runs[0].text = text
    for r in runs[1:]:
        r.text = ''


# ================================================================== supplementary_tables.docx
doc = docx.Document(os.path.join(ROOT, 'tables', 'supplementary_tables.docx'))

# ---------- locate tables by caption number (robust to reordering) ----------
from docx.oxml.ns import qn as _qn


def _table_map(document):
    body = document.element.body
    cur, out = None, {}
    for el in body.iterchildren():
        if el.tag == _qn('w:p'):
            txt = ''.join(t.text or '' for t in el.iter(_qn('w:t'))).strip()
            if txt.startswith('sTable S'):
                cur = txt.split('.')[0].replace('sTable S', '').strip()
        elif el.tag == _qn('w:tbl') and cur:
            out[cur] = _DocxTable(el, document)
    return out


from docx.table import Table as _DocxTable  # noqa: E402
TM = _table_map(doc)

# ---------- sTable S3 (tables[2]): four algorithms ----------
algo = S2['algo_res']
t = TM['6']
METS = ['auc', 'auprc', 'sensitivity', 'specificity', 'ppv', 'npv', 'f1']
for ri, name in enumerate(['Logistic Regression', 'Gradient Boosting',
                           'Linear SVM', 'Gaussian Naive Bayes'], start=1):
    res = algo[name]
    if name == 'Linear SVM':
        set_cell(t.rows[ri].cells[0], 'SVM')
    vals = []
    for m in METS:
        point = res['pooled_point'][m]
        lo, hi = res['boot_ci'][m]
        vals.append(fci(point, lo, hi))
    for ci_, v in enumerate(vals, start=1):
        set_cell(t.rows[ri].cells[ci_], v)

# ---------- sTable S4 (tables[3]): four models ----------
desc = S2['desc']
KEYS = ['Clinical_Only', 'Single_Immature', 'Clinical_Immature', 'Parsimonious_Pipeline']
t = TM['3']
set_cell(t.rows[0].cells[4], PARS_SEC)
ROW_MET = ['auc', 'auprc', 'sensitivity', 'specificity', 'ppv', 'npv', 'f1']
for ri, m in enumerate(ROW_MET, start=1):
    for cj, k in enumerate(KEYS, start=1):
        s = summ[k]
        point = s['pooled_point'][m]
        lo, hi = s['boot_ci'][m]
        set_cell(t.rows[ri].cells[cj], fci(point, lo, hi))
# dAUROC / NRI / IDI / p (rows 8-11; combined and simplified columns only)
# taken from cmp_rows (the same paired bootstrap draws), independent of desc
_cmp = S2['cmp_rows']
for cj, row in [(3, _cmp[0]), (4, _cmp[1])]:
    set_cell(t.rows[8].cells[cj], fci(*row['d_auc'], plus=True))
    set_cell(t.rows[9].cells[cj], fci(*row['nri'], plus=True))
    set_cell(t.rows[10].cells[cj], fci(*row['idi'], plus=True))
    set_cell(t.rows[11].cells[cj], fp_cmp(row['p_d_auc']))

# ---------- sTable S5 (tables[4]): paired comparisons ----------
t = TM['4']
cmp_rows = S2['cmp_rows']
ROWLAB = [
    f'Combined model (8 predictors) vs. Clinical model (7 predictors)',
    f'{PARS_SEC} vs. Clinical model (7 predictors)',
    f'{PARS_SEC} vs. Combined model (8 predictors)',
]
from statsmodels.stats.multitest import multipletests as _mt
_q3 = _mt([r['p_d_auc'] for r in cmp_rows], method='fdr_bh')[1]
def _pv(x):
    return 'P < 0.0001' if x < 0.0001 else 'P = %.4f' % x


for ri, (lab, row) in enumerate(zip(ROWLAB, cmp_rows), start=1):
    set_cell(t.rows[ri].cells[0], lab)
    set_cell(t.rows[ri].cells[1],
             f"{fci(*row['d_auc'], plus=True)}; q = {_q3[ri - 1]:.4f} "
             f"({_pv(row['p_d_auc'])})")
    set_cell(t.rows[ri].cells[2],
             f"{fci(*row['d_auprc'], plus=True)}; {_pv(row['p_d_auprc'])}")
    set_cell(t.rows[ri].cells[3],
             f"{fci(*row['d_brier'], plus=True)}; {_pv(row['p_d_brier'])}")
    set_cell(t.rows[ri].cells[4],
             f"{fci(*row['nri'], plus=True)}; {_pv(row['p_nri'])}")
    set_cell(t.rows[ri].cells[5],
             f"{fci(*row['idi'], plus=True)}; {_pv(row['p_idi'])}")

# ---------- sTable S6 (tables[5]): single markers + three models ----------
t = TM['5']
markers = S2['markers']


def _marker(name):
    """Accept either hyphen style for CD10 (superscript minus vs plain hyphen)."""
    for cand in (name, name.replace('\u207b', '-'), name.replace('-', '\u207b')):
        if cand in markers:
            return markers[cand]
    raise KeyError(name)


mrows = ['INR', 'Platelet count', 'CD10\u207b neutrophils, %', 'Day 1 SIC score']
for ri, mn in enumerate(mrows, start=1):
    mk = _marker(mn)
    set_cell(t.rows[ri].cells[1], fci(*mk['auc']))
    set_cell(t.rows[ri].cells[2], fci(*mk['auprc']))
for ri, k in [(5, 'Clinical_Only'), (6, 'Clinical_Immature'), (7, 'Parsimonious_Pipeline')]:
    s = summ[k]
    set_cell(t.rows[ri].cells[1], fci(s['pooled_point']['auc'], *s['boot_ci']['auc']))
    set_cell(t.rows[ri].cells[2], fci(s['pooled_point']['auprc'], *s['boot_ci']['auprc']))
set_cell(t.rows[7].cells[0], PARS_SEC)

# ---------- sTable S7 (tables[6]): NSIC subgroup (descriptive only) ----------
t = TM['7']
sub = S2['subgroup']
sub_keys = ['Clinical model (7 predictors)', 'Combined model (8 predictors)',
            f'Simplified model ({N_PARS} predictors)']
set_cell(t.rows[0].cells[3], PARS_HEAD + ' ')
set_cell(t.rows[1].cells[0], 'Patients, n (incident SIC events)')

# insert the non-events row before the events row
_row_ev = t.rows[2]
_new = t.add_row()
t._tbl.remove(_new._tr)
_row_ev._tr.addprevious(_new._tr)

set_cell(t.rows[2].cells[0], 'No Day-3 SIC, predicted risk, median (IQR)')
set_cell(t.rows[3].cells[0], 'Day-3 SIC, predicted risk, median (IQR)')
for cj, sk in enumerate(sub_keys, start=1):
    d = sub['models'][sk]
    set_cell(t.rows[1].cells[cj], f"{sub['n']} ({sub['events']})")
    set_cell(t.rows[2].cells[cj], '%.3f (%.3f-%.3f)' % (d['median_nonevent'], *d['iqr_nonevent']))
    set_cell(t.rows[3].cells[cj], '%.3f (%.3f-%.3f)' % (d['median_event'], *d['iqr_event']))

# ---------- table notes ----------
note_s3 = ('Note: Values are presented as cross-validated point estimates with 95% confidence '
           'intervals. All four algorithms were fitted to the same eight candidate predictors measured '
           'on Day 1 (INR, lactate, platelet count, SOFA cardiovascular, renal and respiratory '
           'subscores, mechanical ventilation status, and CD10- neutrophil percentage) to '
           'predict SIC status on Day 3. SIC on Day 3 was defined '
           'as the positive outcome. Model performance was evaluated using nested cross-validation '
           '(20 repeats of five-fold outer cross-validation; hyperparameter tuning was performed '
           'within each outer training fold). Point estimates were calculated from pooled '
           'out-of-fold predictions averaged across repeats for each patient, and 95% confidence '
           'intervals were derived from 2,000 patient-level bootstrap resamples of these pooled '
           'predictions. Multiple imputation and feature standardization were performed within '
           'each training fold only. Sensitivity, specificity, PPV, NPV, and F1 score were '
           'calculated at Youden-optimal thresholds derived from inner-cross-validation '
           'out-of-fold predictions within each training fold. AUROC was used as the primary metric for model comparison.')
note_s4 = ('Note: Values are presented as cross-validated point estimates with 95% confidence '
           'intervals. The positive outcome was SIC on Day 3. The clinical model included seven '
           'Day 1 predictors: INR, lactate, platelet count, SOFA cardiovascular, renal, and '
           'respiratory subscores, and mechanical ventilation status. The biomarker-only model '
           'included Day 1 CD10\u207b neutrophil percentage, and the combined model incorporated '
           'CD10\u207b neutrophil percentage in addition to the seven clinical predictors. The '
           f'simplified {N_PARS}-predictor model was defined as the predictor set most '
           'frequently retained across the 100 outer training folds, with features ranked by '
           'recursive feature elimination (RFE) inside each fold; it comprised '
           f'{", ".join(pars_display_note)}. Model performance was '
           'evaluated using nested cross-validation (20 repeats of five-fold outer '
           'cross-validation), in which the entire modelling pipeline\u2014multiple imputation, '
           'feature standardization, elastic-net hyperparameter tuning, feature-importance '
           'ranking, recursive feature elimination (RFE), and classification-threshold selection\u2014'
           'was re-executed independently within each outer training fold, so that feature '
           'selection was fully separated from performance evaluation; each outer test fold was '
           'predicted only once. Operating-point metrics (sensitivity, specificity, PPV, NPV, '
           'and F1 score) were calculated at Youden-optimal thresholds derived from inner-cross-validation '
           'out-of-fold predictions within each training fold. Point estimates were calculated from pooled out-of-fold '
           'predictions averaged across repeats for each patient, and 95% confidence intervals '
           'were derived from 2,000 patient-level bootstrap resamples of these pooled predictions. '
           '\u0394AUROC, continuous NRI, and IDI were calculated relative to the clinical model, '
           'using the pooled patient-level predictions so that every metric in this table shares '
           'the same estimator as the ROC, calibration, and decision-curve analyses in Figure 3. '
           'Confidence intervals for these contrasts were derived from 2,000 paired patient-level '
           'bootstrap resamples, in which both models were re-evaluated within the same resample. '
           'The four predictors were those most frequently retained across the 100 outer folds; '
           'because this set was defined using the same cohort, the resulting estimate is best '
           'interpreted as the performance of the full modelling strategy, and a strategy-level '
           'estimate that re-selects predictors within every fold yielded an AUROC of 0.938. '
           'Continuous NRI and IDI are reported for completeness and should be regarded as '
           'exploratory; both are unstable at this sample size.')

# ---- calibration slope / intercept (written by figure4() in notebook 01) ----
_calib_csv = os.path.join(OUT_T, 'calibration_slope_intercept.csv')
if os.path.exists(_calib_csv):
    import pandas as _pd
    _c = _pd.read_csv(_calib_csv)
    _row = _c.loc[_c['model'] == 'Parsimonious_Pipeline']
    if len(_row):
        note_s4 += (' Calibration slope and intercept for the simplified model were %.2f and %.2f, '
                    'respectively (unpenalised logistic recalibration of the pooled out-of-fold '
                    'predictions); a slope substantially above 1 indicates that predicted '
                    'probabilities were compressed towards the middle, i.e. under-confident.'
                    % (float(_row['slope'].iloc[0]), float(_row['intercept'].iloc[0])))

note_s5 = ('Note: Changes (\u0394) are expressed relative to the reference model listed second in '
           'each comparison. Positive values for \u0394AUROC, \u0394AUPRC, continuous NRI, and IDI '
           'favor the first model, whereas negative values for \u0394Brier score indicate improved '
           'performance. The clinical model included seven Day 1 clinical predictors; the combined '
           'model additionally included Day 1 CD10\u207b neutrophil percentage; and the '
           f'simplified model included {", ".join(pars_display_note)}. P values and 95% '
           'confidence intervals for \u0394AUROC, \u0394AUPRC, and \u0394Brier score were derived '
           'from 2,000 paired patient-level bootstrap resamples: in each resample both models were '
           'scored on the same resampled patients and the difference retained, which accounts for '
           'the paired design and for the correlation between cross-validation repeats sharing the '
           'same patients. Nadeau\u2013Bengio corrected resampled t-tests on repeat-level estimates '
           'are provided in the accompanying analysis file as a sensitivity analysis. '
           'Only \u0394AUROC was formally tested, across three contrasts; q values are '
           'Benjamini\u2013Hochberg adjusted across these three tests. \u0394AUPRC, \u0394Brier '
           'score, continuous NRI and IDI are reported descriptively, without multiplicity '
           'adjustment, and should be regarded as exploratory and unstable at this sample size.')
note_s6 = ('Note: Values are presented as cross-validated point estimates with 95% confidence '
           'intervals derived from 2,000 patient-level bootstrap resamples. All individual '
           'predictors and multivariable models were evaluated using the same nested '
           'cross-validation splits (20 repeats of five-fold). The positive outcome was SIC on '
           'Day 3. Day 1 SIC score was included as a baseline composite benchmark. For platelet '
           'count, the prediction direction was reversed because lower platelet counts indicate '
           'greater SIC risk; all other predictors were evaluated in their natural direction. The '
           'clinical model included seven Day 1 predictors: INR, lactate, platelet count, SOFA '
           'cardiovascular, renal, and respiratory subscores, and mechanical ventilation status. '
           'The combined model additionally included Day 1 CD10\u207b neutrophil percentage. The '
           f'simplified model included {", ".join(pars_display_note)}.')
note_s7 = ('Note: This exploratory analysis was restricted to the %d patients without SIC on Day 1, '
           'of whom %d developed SIC by Day 3. Predicted risks are nested cross-validated pooled '
           'out-of-fold probabilities and are summarised as median and interquartile range separately '
           'for patients who did and did not develop SIC on Day 3. Because only %d incident events '
           'occurred, no discrimination or calibration statistic, and no between-group test, was '
           'calculated; an AUC or a P value in this subgroup would be dominated by sampling variance '
           'and would not support inference. The values are therefore reported descriptively and are '
           'hypothesis-generating only.' % (sub['n'], sub['events'], sub['events']))

def _apply_note_by_key(document, mapping):
    """Locate the first Note of a block by sTable number, overwrite it and drop duplicates."""
    from docx.oxml.ns import qn as _q
    body = document.element.body
    cur, first = None, {}
    extra = []
    for el in body.iterchildren():
        if el.tag == _q('w:p'):
            txt = ''.join(t.text or '' for t in el.iter(_q('w:t'))).strip()
            if txt.startswith('sTable S'):
                cur = txt.split('.')[0].replace('sTable S', '').strip()
            elif txt.startswith('Note:') and cur:
                if cur not in first:
                    first[cur] = el
                else:
                    extra.append(el)
    for el in extra:
        body.remove(el)
    for key, txt in mapping.items():
        el = first.get(key)
        if el is None:
            continue
        from docx.text.paragraph import Paragraph
        set_para(Paragraph(el, document), txt)


_apply_note_by_key(doc, {'6': note_s3, '3': note_s4, '4': note_s5, '5': note_s6, '7': note_s7})

# if the simplified set is no longer four predictors, the headers were already replaced above.
# ---------- physical reordering of the tables ----------
def reorder_supp_tables(document, order=('1', '2', '3', '4', '5', '6', '7')):
    body = document.element.body
    kids = list(body.iterchildren())
    starts = []
    for idx, el in enumerate(kids):
        if el.tag == _qn('w:p'):
            txt = ''.join(t.text or '' for t in el.iter(_qn('w:t'))).strip()
            if txt.startswith('sTable S'):
                starts.append(idx)
    if not starts:
        return False
    blocks, keys = {}, []
    for i, st in enumerate(starts):
        en = starts[i + 1] if i + 1 < len(starts) else next(
            (j for j in range(st, len(kids)) if kids[j].tag.endswith('sectPr')), len(kids))
        txt = ''.join(t.text or '' for t in kids[st].iter(_qn('w:t'))).strip()
        key = txt.split('.')[0].replace('sTable S', '').strip()
        blocks[key] = kids[st:en]
        keys.append(key)
    want = [k for k in order if k in blocks] + [k for k in keys if k not in order]
    if want == keys:
        return False
    prefix = kids[:starts[0]]
    sect = [el for el in kids[starts[-1]:] if el.tag.endswith('sectPr')]
    for el in kids:
        body.remove(el)
    for el in prefix:
        body.append(el)
    for k in want:
        for el in blocks[k]:
            body.append(el)
    for el in sect:
        body.append(el)
    return True


if reorder_supp_tables(doc):
    print('Tables reordered: S1, S2, S3, S4, S5, S6, S7', flush=True)

doc.save(os.path.join(ROOT, 'tables', 'supplementary_tables.docx'))
print('Updated tables/supplementary_tables.docx', flush=True)


# ================================================================== figure legend update
def replace_in_doc(path_in, path_out, replacements):
    d = docx.Document(path_in)
    for para in d.paragraphs:
        full = para.text
        new = full
        for old, newt in replacements:
            if old in new:
                new = new.replace(old, newt)
        if new != full:
            set_para(para, new)
    d.save(path_out)
    print(f'Updated {os.path.basename(path_out)}', flush=True)


# Figure 3/4 legends (in figures/Figure legend.docx)
leg_repl = [
    ('Bars represent mean feature importance across repeated cross-validation runs, and error '
     'bars indicate ±1 standard deviation.',
     'Bars represent mean feature importance across repeated cross-validation runs on the full '
     'derivation dataset, and error bars indicate ±1 standard deviation.'),
    ('Points represent mean AUROC estimates from 20 repeats of five-fold cross-validation, with '
     'error bars indicating 95% confidence intervals.',
     'Points represent mean AUROC estimates from 20 repeats of five-fold cross-validation on the '
     'full derivation dataset, with error bars indicating 95% confidence intervals. This panel '
     'illustrates the derivation of the final model only; unbiased performance estimates were '
     'obtained by nested cross-validation (Figure 4).'),
    ('The four-predictor model was selected as the parsimonious model and included CD10\u207b '
     'neutrophil percentage, platelet count, SOFA cardiovascular subscore, and SOFA renal '
     'subscore.',
     f'The {N_PARS}-predictor model was selected as the parsimonious model and included '
     f'{", ".join(pars_display_note)}.'),
    ('and the four-predictor parsimonious model',
     f'and the {N_PARS}-predictor parsimonious model'),
    ('Model performance was evaluated using 20 repeats of five-fold cross-validation, with ROC '
     'curves constructed from out-of-fold predictions. AUROCs are presented with 95% confidence '
     'intervals derived from repeat-level estimates.',
     'Model performance was evaluated using nested cross-validation (20 repeats of five-fold '
     'outer cross-validation), with the entire feature-selection and tuning pipeline re-executed '
     'within each outer training fold; ROC curves were constructed from pooled out-of-fold '
     'predictions averaged across repeats for each patient. AUROCs are presented with 95% '
     'confidence intervals derived from 2,000 patient-level bootstrap resamples.'),
    ('Patients were grouped according to predicted risk, and the observed proportion of Day 3 '
     'SIC was plotted against the mean predicted probability within each risk group.',
     'Patients were grouped according to nested cross-validated predicted risk, and the observed '
     'proportion of Day 3 SIC was plotted against the mean predicted probability within each '
     'risk group.'),
    ('Net benefit was evaluated across a range of threshold probabilities using cross-validated '
     'out-of-fold predicted probabilities.',
     'Net benefit was evaluated across a range of threshold probabilities using nested '
     'cross-validated out-of-fold predicted probabilities.'),
]
replace_in_doc(os.path.join(ROOT, 'figures', 'Figure legend.docx'),
               os.path.join(ROOT, 'figures', 'Figure legend.docx'), leg_repl)

sleg_repl = [
    ('developed using the same seven Day 1 clinical predictors',
     'developed using the same eight candidate predictors (seven clinical variables plus '
     'CD10- neutrophil percentage)'),
    ('linear support vector machine', 'support vector machine'),
    ('Model performance was evaluated using 20 repeats of five-fold cross-validation, and ROC '
     'curves were based on out-of-fold predictions. AUROCs are presented with 95% confidence '
     'intervals derived from repeat-level cross-validation estimates.',
     'Model performance was evaluated using nested cross-validation (20 repeats of five-fold '
     'outer cross-validation with hyperparameter tuning confined to each outer training fold), '
     'and ROC curves were based on pooled out-of-fold predictions averaged across repeats for '
     'each patient. AUROCs are presented with 95% confidence intervals derived from 2,000 '
     'patient-level bootstrap resamples; shaded bands represent pointwise 95% confidence '
     'intervals across repeats.'),
    ('Model performance was evaluated using repeated cross-validation, and AUROCs are presented '
     'with 95% confidence intervals derived from repeat-level estimates.',
     'Model performance was evaluated using repeated cross-validation, and AUROCs are presented '
     'with 95% confidence intervals derived from 2,000 patient-level bootstrap resamples of '
     'pooled out-of-fold predictions.'),
    ('the four-predictor parsimonious model in the exploratory subgroup',
     f'the {N_PARS}-predictor parsimonious model in the exploratory subgroup'),
]
replace_in_doc(os.path.join(ROOT, 'figures', 'Supplementary Figure.docx'),
               os.path.join(ROOT, 'figures', 'Supplementary Figure.docx'), sleg_repl)

print('\nTables and figure legends updated.')
