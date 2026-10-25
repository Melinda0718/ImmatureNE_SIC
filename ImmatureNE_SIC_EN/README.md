# NE-99: Day-1 neutrophil phenotypes and Day-3 SIC prediction

Analysis package accompanying the manuscript on Day-1 neutrophil phenotypes (flow cytometry) and
the prediction of sepsis-induced coagulopathy (SIC) on Day 3 in sepsis.

## Layout

```
NE_99_261022_EN/
|-- code/
|   |-- classical_analysis/02_classical_statistics.ipynb   # Tables 1-3, sTables S1-S2, Figures 2, S2-S5
|   `-- prediction_model/
|       |-- 01_prediction_nestedCV.ipynb                   # nested CV, Figures 3-4, sTables S3-S7
|       |-- pipeline.md                                    # pipeline description (with diagram)
|       `-- layout/                                        # Word/LaTeX table builders
|-- data/           # de-identified analysis datasets (D1/D3/D7 + healthy controls)
|-- outputs/        # intermediate tables and single-panel figures
|-- tables/         # formatted Word/LaTeX tables
`-- figures/        # final figures (png / pdf / tif, 300 dpi)
```


## Requirements

Python 3.13 with the packages listed in `requirements.txt` (scikit-learn 1.6.1, statsmodels 0.14.6,
pandas 2.2.3, matplotlib, seaborn, python-docx). LaTeX output needs `xelatex`.

```bash
pip install -r requirements.txt
```

## Running

```bash
# 1) prediction model: nested cross-validation (about 30 min single core)
#    set NE_JOBS=1 on machines with few cores
jupyter nbconvert --to notebook --execute --inplace code/prediction_model/01_prediction_nestedCV.ipynb

# 2) classical statistics: tables and descriptive figures (about 1 min)
jupyter nbconvert --to notebook --execute --inplace code/classical_analysis/02_classical_statistics.ipynb

# 3) rebuild the formatted tables (csv + tex + pdf + docx)
python tables/build_tables.py
```

## Methodological notes

- The outer loop (5-fold x 20 repeats = 100 folds) is used for evaluation only; imputation,
  standardisation, feature selection, tuning and threshold derivation are repeated inside every
  outer training fold, and the test fold is predicted once.
- Performance is summarised from pooled out-of-fold predictions (one value per patient, averaged
  over the 20 repeats); 95% CIs come from 2,000 patient-level bootstrap resamples.
- Paired comparisons are bootstrapped within the same resamples; Nadeau-Bengio corrected t-tests
  are provided as a sensitivity analysis.
- The simplified model has two definitions: the strategy-level estimate (feature selection repeated
  in every fold, primary) and a fixed four-predictor specification (secondary, descriptive).

## Licence

Code: MIT. Data: available from the corresponding author on reasonable request, subject to the
ethics approval of the participating centre.

## Data availability

Patient-level data are **not** distributed with this repository. De-identified datasets are available
from the corresponding author upon reasonable request and after approval by the Institutional Review
Board of Zhongda Hospital, Southeast University (contact: [corresponding author], [email]).

Data are shared as csv files with study codes only (`P001`-`P099` for patients, `HC001`-`HC017` for
healthy controls); source identifiers and calendar dates are removed.

## Expected input files

Place the following files in `data/` before running the notebooks (column names follow the manuscript):

| File | Content |
|---|---|
| `data/NE_baseline_D1.csv` | 99 patients, Day-1 clinical and flow-cytometry data (outcome column `SIC`) |
| `data/NE_baseline_D3.csv` | Day-3 SIC status per patient |
| `data/NE_baseline_D7.csv` | Day-7 SIC status per patient |
| `data/NE_HC.csv` | Healthy controls (reference group) |

The patient identifier column must be named `Code`; it is used to link the three time points.
