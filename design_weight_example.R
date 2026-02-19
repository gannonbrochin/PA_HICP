# =============================================================================
# Author: Gannon Brochin 
# Title: Applying NHIS Survey Design Weights in R
# Source: nhis_pa_hicp_2020.R (Brochin, 2025–2026, currently unpublished)
# Full source available at https://github.com/gannonbrochin/PA_HICP 
# For use with NHIS 2020 Adult Sample File 
# =============================================================================
# The NHIS uses a complex, multi-stage sampling design:
# Stage 1: The country is divided into geographic areas, and a subset of areas are selected 
#   (these become the primary sampling units / PSUs, stored in PPSU)
# Stage 2: Within those areas, households are sampled
# Stage 3: Within households, one adult is randomly selected to be interviewed

# To produce nationally representative estimates, analyses must account for:
#   - Sampling weights   (WTFA_A)  — adjust for unequal selection probability
#                                     and nonresponse
#   - Stratification     (PSTRAT)  — variance strata
#   - Clustering / PSUs  (PPSU)    — primary sampling units (geographic clusters)  
#
# Ignoring the design leads to biased point estimates (unweighted) and
# understated standard errors (no clustering/stratification adjustment).
# NOTE: Unweighted values are needed for items like sensitivity testing,
#       missingness descriptions, and sample descriptions so they should be calculated
#       and included in full analyses on top of these weighted values.
# =============================================================================

# NOTE: Uncomment the line below to install required packages if not already installed.
# install.packages(c("survey", "readr", "janitor"))

suppressPackageStartupMessages({
  library(survey)
  library(readr)
  library(janitor)
})

# ---------- 1. Load data ----------
# Replace with your NHIS CSV path
d <- read_csv("adult20.csv", show_col_types = FALSE) |> clean_names()

# ---------- 2. Handle lonely PSUs ----------
# Some strata may contain only one PSU after subsetting. This option applies
# a conservative degrees-of-freedom adjustment rather than throwing an error.
options(survey.lonely.psu = "adjust")

# ---------- 3. Declare the survey design ----------
# This tells the survey package HOW the sample was drawn so that all
# subsequent analyses correctly weight observations and estimate variances.
#
#   ids     = ~ppsu       → cluster identifier (primary sampling unit)
#   strata  = ~pstrat     → variance stratum
#   weights = ~wtfa_a     → sampling weight for the Sample Adult file
#   nest    = TRUE        → PSU IDs are only unique WITHIN strata
#
des <- svydesign(
  ids     = ~ppsu,
  strata  = ~pstrat,
  weights = ~wtfa_a,
  nest    = TRUE,
  data    = d
)

# ---------- 4. Example analyses using the design ----------

# (a) Weighted mean of a continuous variable (e.g., age)
svymean(~agep_a, design = des, na.rm = TRUE)

# (b) Weighted prevalence of a binary outcome
#     (e.g., proportion reporting chronic pain — paifrq3m_a codes 3 or 4)
d$chronic_pain <- ifelse(d$paifrq3m_a %in% c(3, 4), 1L,
                         ifelse(d$paifrq3m_a %in% c(1, 2), 0L, NA))
des <- update(des, chronic_pain = d$chronic_pain)
svymean(~chronic_pain, design = des, na.rm = TRUE)

# (c) Weighted prevalence by subgroup, with 95% CIs
svyby(~chronic_pain, ~sex_a, design = des, svymean,
      vartype = "ci", na.rm = TRUE)

# (d) Design-based regression (linear probability model)
fit <- svyglm(chronic_pain ~ agep_a + factor(sex_a),
              design = des, family = gaussian())
summary(fit)
