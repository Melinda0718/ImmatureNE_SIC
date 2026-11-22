# Circulating CD10-negative neutrophil percentage in sepsis-induced coagulopathy

Analysis code and de-identified data for the study:

> **Circulating CD10⁻ Neutrophil Percentage as a Novel Biomarker for Sepsis-Induced
> Coagulopathy Diagnosis: A Prospective Cohort Study**

Single-centre prospective observational cohort (Zhongda Hospital, Southeast University,
May 2024 – April 2025): 99 patients with sepsis (56 with SIC and 43 without SIC on Day 1)
plus 17 healthy volunteers. Whole-blood flow cytometry quantified the CD10⁻ neutrophil
percentage and CD63 / CD16 / CD62L / CD66b expression on Days 1, 3 and 7.

## Contents

| Path | Description |
|---|---|
| `data/NE_baseline_D1.csv` | Day 1 clinical and flow-cytometry variables, 99 patients × 110 variables |
| `data/NE_baseline_D3.csv` | Day 3 variables, 99 patients × 56 variables |
| `data/NE_baseline_D7.csv` | Day 7 variables, 99 patients × 56 variables |
| `data/NE_HC.csv` | Healthy controls, 17 × 26 variables |
| `data/DATA_DICTIONARY.md` | Variable list and definitions |
| `code/` | Analysis scripts (**see `code/README.md`**) |
| `results/tables/` | Intermediate CSV outputs of the table pipeline |

## De-identification

Patient records are relabelled `P001`–`P099`; the same labels are used across Day 1, 3
and 7 files. Healthy controls are `HC-001`–`HC-017`. The original laboratory codes, the
flow-cytometry file names and all other direct identifiers have been removed.

## Environment

```
python 3.13
pandas, numpy, scipy, statsmodels, scikit-learn, matplotlib
```

Install with:

```bash
pip install -r requirements.txt
```

## Licence

Code: MIT. Data: CC BY 4.0 (see `LICENSE`).

## Contact

Corresponding author: Xiaojing Wu (36410257@qq.com); Yingzi Huang (yz_huang@126.com).
