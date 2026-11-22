# Data dictionary

Common identifiers and grouping variables
-----------------------------------------
| Variable | Meaning |
|---|---|
| `Code` | De-identified patient label, `P001`-`P099`; identical across Day 1/3/7 |
| `SIC` | Sepsis-induced coagulopathy at that time point (0/1) |
| `SIC_D1` / `SIC_D3` / `SIC_D7` | SIC status on that day |
| `SIC_Group` | Group label used in the manuscript (SIC / NSIC) |
| `SIC_score` | Total SIC score (PLT + INR + SOFA components) |
| `SIC_PLT_pts` / `SIC_INR_pts` / `SIC_SOFA_pts` | SIC subscores |
| `DIC` | Disseminated intravascular coagulation (ISTH overt DIC) |
| `ID` (healthy controls only) | `HC-001`-`HC-017` |

Clinical variables (`NE_baseline_D1.csv`)
-----------------------------------------
Demographics and severity: `Age`, `Gender`, `Weight`, `Height`, `BMI`, `APACHEII`, `SOFA`,
`Septic Shock`, comorbidities (`hypertension`, `diabetes`, `chronic_heart_disease`,
`chronic_pulmonary_disease`, `chronic_kidney_disease`, `chronic_liver_disease`,
`cerebral_vascular_accident`, `tumor`, `Thyroid dysfunction`), infection site indicators
(`Infection site_*`), pathogen categories (`G+`, `G-`, `Virus`, `fungal`),
treatments (`Vasoactive drug use`, `Heparin use`, `CRRT_D1`, `VT_D1`) and outcomes
(`LOS of Hospital`, `LOS of ICU`, `28-day mortality`, `ICU mortality`, `Hospital mortality`).

Laboratory variables: `WBC`, `NE count`, `LY count`, `MO count`, `NE/LY`, `PLT`, `PT`, `APTT`,
`FIB`, `INR`, `FDP`, `D_dimer`, `PTA`, `ATIII`, `TT`, `pH`, `HCO3_`, `Ca2+`, `lac`, `PCT`,
`CRP`, `AST`, `ALT`, `DBIL`, `TBIL`, `ALB`, `BUN`, `Scr`, `TM`, `TAT`, `t_PAIC`, `PIC`.

SOFA components: `SOFA D1 Respiratory System`, `SOFA D1 Platelets`, `SOFA D1 Bilirubin`,
`SOFA D1 Circulation`, `SOFA D1 GCS Score`, `SOFA D1 Creatinine`, `SOFA D1 Urine Output`,
`SOFA_Creatinine`.

Flow cytometry (Day 1/3/7)
--------------------------
For total neutrophils (`All_Neu_*`), the CD10-negative subset (`immature_Neu_*`) and the
CD10-positive subset (`mature_Neu_*`):

| Suffix | Meaning |
|---|---|
| `_percentage` | Subset as percentage of the parent gate |
| `_Mean_(CDxx_MFI)` | Mean fluorescence intensity of the marker in that subset |
| `_Median_(CDxx_MFI)` | Median fluorescence intensity of the marker in that subset |

Markers: CD63 (degranulation), CD62L (L-selectin), CD16 (FcγRIII), CD66b (activation).
`immature_Neu_percentage` is stored as a 0-1 fraction; multiply by 100 for the percentage.

Notes
-----
- Group assignment for Day 3 and Day 7 uses the **Day 1** SIC status (fixed cohort grouping).
- `SIC_D3` and `SIC_D7` are provided for completeness but were not used in the manuscript.
- No imputation was applied; missing values are left blank.
