# Physical Activity and High-Impact Chronic Pain — NHIS 2020

Analysis code for *Physical Activity and High-Impact Chronic Pain Among US Adults: A Cross-Sectional Dose–Response Analysis of the 2020 National Health Interview Survey* (Gannon Brochin).

## Data

`adult20csv.zip` is the 2020 NHIS Sample Adult public-use file, unmodified, from the National Center for Health Statistics:
https://www.cdc.gov/nchs/nhis/documentation/2020-nhis.html

The HICP activity-limitation item (`PAIWKLM3M_A`) was fielded only July–December 2020 (`INTV_QRT` 3–4), so the analytic sample is restricted to that period. Survey-weighted models treat it as a domain (`subset()`) of the full-year design.

## Reproduce

1. Open `PA and HICP.Rproj` (or set the working directory to this folder) and unzip the data:

   ```r
   unzip("adult20csv.zip")
   ```

2. Install packages (versions used are in `sessioninfo.txt`):

   ```r
   install.packages(c("tidyverse", "survey", "broom", "janitor", "emmeans", "officer",
                      "flextable", "tableone", "here"))
   ```

3. Run the main analysis, then the supplemental analyses:

   ```r
   source("nhis_pa_hicp_2020.R")
   source("nhis_pa_hicp_2020_supplement.R")
   ```

All results are written to `Output/`.

| Script | Produces |
|---|---|
| `nhis_pa_hicp_2020.R` | Table 1, Table 2, Figure 2, prevalence ratios, trend test, spline thresholds |
| `nhis_pa_hicp_2020_supplement.R` | Supplement 2 Tables S1–S8 as CSV files (`supplement_S*.csv`, `prevalence_ratios_model_A.csv`, `prevalence_ratios_model_B.csv`) |

Expected checkpoints: analytic n = 16,918; HICP cases = 1,392; weighted HICP prevalence = 7.45%.
