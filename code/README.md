# Code

| File | What it does | Needs |
|---|---|---|
| `01_main_analysis.ipynb`* | Main analysis: Table 1-3, sTable S1/S2/S8, Figure 2 panels, Supplementary Figures S2-S5 | `data/`, `code/tables.py` |
| `02_dic2025_scoring.py` | Recomputes Day-1 overt DIC with the ISTH 2025 criteria | `data/` |
| `03_dic2025_roc.py` | ROC of the five Day-1 markers for overt DIC | `data/`, `code/02_dic2025_scoring.py` |

Run from the repository root:

```bash
python code/02_dic2025_scoring.py
python code/03_dic2025_roc.py
jupyter lab code/01_main_analysis.ipynb
```

`*` the patient-level data are not distributed; request them from the corresponding author
and place the CSVs in `data/`.

`01_main_analysis.ipynb` imports `tables.py`, the module that computes and formats every
table from `data/`. It is provided with the repository; the Word/LaTeX export helpers inside it
expect the manuscript templates and are skipped when those files are absent.
