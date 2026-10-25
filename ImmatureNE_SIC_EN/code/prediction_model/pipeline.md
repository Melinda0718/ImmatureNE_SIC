# Prediction model pipeline

Nested cross-validation for the Day-1 -> Day-3 SIC prediction model.

## Overview

```
data/NE_baseline_D1.csv + NE_baseline_D3.csv   (not distributed; see README)
|-- intersect on Code -> 99 patients (Day-3 SIC 40 / NSIC 59)
   |-- outer RepeatedStratifiedKFold (5 folds x 20 repeats = 100 folds)
      |-- training fold (~79 patients)
      |  |-- imputation (MICE) -> standardisation      (outer training data only)
      |  |-- fixed models (7 / 1 / 8 predictors): elastic-net logistic regression,
      |  |   hyperparameters tuned by inner 5-fold CV
      |  |-- simplified model: RFECV inside the training fold (8 -> 3 predictors),
      |  |   subset size chosen by inner-CV AUROC
      |  |-- Youden threshold from inner out-of-fold predictions
      |-- test fold (~20 patients): predicted once
   |-- pooled OOF: 20 repeats averaged per patient
      |-- 2,000 patient-level bootstrap -> 95% CI
      |-- paired bootstrap -> dAUROC / dAUPRC / dBrier / NRI / IDI
```

## Definitions

| Item | Value |
|---|---|
| Outcome | Day-3 ISTH SIC score >= 4 |
| Candidates | 7 clinical predictors + CD10- neutrophil percentage |
| Outer loop | 5-fold x 20 repeats, stratified, seed 42 |
| Inner loop | 5-fold (RFECV and hyperparameter search) |
| Selection | RFECV, subset sizes 8 to 3 |
| Final model | 4 most frequently selected predictors |
| CI | 2,000 patient-level bootstrap resamples |

## Reproducing

Run `01_prediction_nestedCV.ipynb` from top to bottom; Part 1 writes `outputs/tables/nested_cv_results.pkl`,
Part 2 writes `outputs/tables/nested_cv_stage2.pkl`, and Part 3 renders the figures into `figures/`.


> The input csv files are not part of this repository; request them from the corresponding author.
> Expected columns and file names are listed in the README.
