# Circulating CD10-negative neutrophil percentage in sepsis-induced coagulopathy

Analysis code for the study:

> **Circulating CD10⁻ Neutrophil Percentage as a Novel Biomarker for Sepsis-Induced
> Coagulopathy Diagnosis: A Prospective Cohort Study**

Single-centre prospective observational cohort (Zhongda Hospital, Southeast University,
May 2024 – April 2025): 99 patients with sepsis (56 with SIC and 43 without SIC on Day 1)
plus 17 healthy volunteers. Whole-blood flow cytometry quantified the CD10⁻ neutrophil
percentage and CD63 / CD16 / CD62L / CD66b expression on Days 1, 3 and 7.

## Contents

| Path | Description |
|---|---|
| `code/01_main_analysis.ipynb` | Main analysis: Table 1–3, sTable S1/S2/S8, Figure 2 panels, Supplementary Figures S2–S5 |
| `code/tables.py` | Module that computes and formats every table |
| `code/02_dic2025_scoring.py` | Day-1 overt DIC recomputed with the ISTH 2025 criteria |
| `code/03_dic2025_roc.py` | ROC of the five Day-1 markers for overt DIC |
| `code/README.md` | How to run the code |
| `data/DATA_DICTIONARY.md` | Variable names, units and coding |
| `results/tables/` | Intermediate CSV outputs of the table pipeline |

## Data availability

The patient-level dataset is **not distributed with this repository**. The de-identified data
are available from the corresponding author on reasonable request:

- Xiaojing Wu — 36410257@qq.com
- Yingzi Huang — yz_huang@126.com

`data/DATA_DICTIONARY.md` documents every variable so that the code can be applied to the data
once they are obtained.

## Environment

```
python 3.13
pandas, numpy, scipy, statsmodels, scikit-learn, matplotlib
```

Install with:

```bash
pip install -r requirements.txt
```

Then, from the repository root:

```bash
jupyter lab code/01_main_analysis.ipynb
python code/02_dic2025_scoring.py
python code/03_dic2025_roc.py
```

## Licence

Code: MIT (see `LICENSE`). No data are licensed here because no patient-level data are
distributed.

## Contact

Yanzi Guo (gyz2002@126.com) — please open an issue or write if you have any problem
running the code.
