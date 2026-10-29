suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(broom)
  library(janitor)
  library(readr)
  library(emmeans)
  library(splines)
  library(ggplot2)
  library(officer)
  library(flextable)
  library(forcats)
  library(tableone)
  library(here)
})

dir.create(here("Output"), showWarnings = FALSE)
data_path <- here("adult20.csv")

d <- read_csv(data_path, show_col_types = FALSE) |> clean_names()

d <- d |>
  mutate(
    chronic_pain = case_when(paifrq3m_a %in% c(3,4) ~ 1L,
                             paifrq3m_a %in% c(1,2) ~ 0L,
                             TRUE ~ NA_integer_),
    hicp = case_when(
      chronic_pain == 1L & paiwklm3m_a %in% c(3,4) ~ 1L,
      chronic_pain == 0L ~ 0L,
      chronic_pain == 1L & paiwklm3m_a %in% c(1,2) ~ 0L,
      TRUE ~ NA_integer_
    ),
    hicp_window = intv_qrt %in% c(3, 4)  # PAIWKLM3M_A fielded Jul–Dec 2020 only
  )

d <- d |>
  mutate(
    modfreqw_raw = as.numeric(modfreqw_a),
    vigfreqw_raw = as.numeric(vigfreqw_a),
    strfreqw_raw = as.numeric(strfreqw_a)
  )

recode_freq <- function(x) {
  x <- as.numeric(x)
  case_when(
    x %in% c(94, 95)          ~ 0,  # never / unable
    x %in% c(96, 97, 98, 99)  ~ NA_real_,
    TRUE                       ~ x
  )
}

recode_min <- function(x) {
  x <- as.numeric(x)
  case_when(
    x >= 996    ~ NA_real_,
    !is.na(x)   ~ x,
    TRUE         ~ NA_real_
  )
}

d <- d |>
  mutate(
    modfreqw_a = recode_freq(modfreqw_a),
    vigfreqw_a = recode_freq(vigfreqw_a),
    strfreqw_a = recode_freq(strfreqw_a),
    modmin_a   = recode_min(modmin_a),
    vigmin_a   = recode_min(vigmin_a),
    modmin_a   = ifelse(!is.na(modfreqw_a) & modfreqw_a == 0 & is.na(modmin_a), 0, modmin_a),
    vigmin_a   = ifelse(!is.na(vigfreqw_a) & vigfreqw_a == 0 & is.na(vigmin_a), 0, vigmin_a),
    mod_min_wk = modfreqw_a * modmin_a,
    vig_min_wk = vigfreqw_a * vigmin_a,
    eq_min_wk  = mod_min_wk + 2 * vig_min_wk,
    meet_aerobic  = ifelse(!is.na(eq_min_wk) & eq_min_wk >= 150, 1L,
                           ifelse(!is.na(eq_min_wk), 0L, NA)),
    meet_strength = ifelse(!is.na(strfreqw_a) & strfreqw_a >= 2, 1L,
                           ifelse(!is.na(strfreqw_a), 0L, NA)),
    pa_cat = case_when(
      !is.na(eq_min_wk) & eq_min_wk == 0 ~ "No PA",
      !is.na(meet_aerobic) & meet_aerobic == 1 &
        !is.na(meet_strength) & meet_strength == 1 ~ "Meets aerobic + strength",
      !is.na(meet_aerobic) & meet_aerobic == 1 &
        !is.na(meet_strength) & meet_strength == 0 ~ "Meets aerobic only",
      !is.na(eq_min_wk) & eq_min_wk > 0 ~ "Some PA but below guidelines",
      TRUE ~ NA_character_
    ),
    pa_cat = factor(pa_cat, levels = c(
      "No PA",
      "Some PA but below guidelines",
      "Meets aerobic only",
      "Meets aerobic + strength"
    ))
  )

d <- d |>
  mutate(
    sex_a = factor(sex_a, levels = c(1,2), labels = c("Male","Female")),
    is_hisp = case_when(hisp_a == 1 ~ 1L, hisp_a == 2 ~ 0L, TRUE ~ NA_integer_),
    race_eth = case_when(
      is_hisp == 1 ~ "Hispanic",
      is_hisp == 0 & raceallp_a == 1 ~ "NH White",
      is_hisp == 0 & raceallp_a == 2 ~ "NH Black",
      is_hisp == 0 & raceallp_a == 4 ~ "NH Asian",
      is_hisp == 0 ~ "NH Other/Multiracial",
      TRUE ~ NA_character_
    ) |> factor(levels = c("NH White","NH Black","NH Asian","Hispanic","NH Other/Multiracial"))
  )

d <- d |>
  mutate(
    educ_cat = case_when(
      educ_a %in% c(0,1,2,3)   ~ "Less than HS",
      educ_a == 4               ~ "HS diploma/GED",
      educ_a %in% c(5,6,7,8)   ~ "Some college/Associate",
      educ_a %in% c(9,10)      ~ "Bachelor or higher",
      TRUE                      ~ NA_character_
    ) |> factor(levels = c("Less than HS","HS diploma/GED",
                            "Some college/Associate","Bachelor or higher")),
    bmi_cat = case_when(
      bmicat_a == 1 ~ "Underweight",
      bmicat_a == 2 ~ "Normal weight",
      bmicat_a == 3 ~ "Overweight",
      bmicat_a == 4 ~ "Obese",
      TRUE          ~ NA_character_
    ) |> factor(levels = c("Normal weight","Underweight","Overweight","Obese")),
    arthritis = case_when(
      arthev_a == 1 ~ 1L,
      arthev_a == 2 ~ 0L,
      TRUE          ~ NA_integer_
    ),
    dep_freq = case_when(
      depfreq_a %in% c(1,2)   ~ "Daily/Weekly",
      depfreq_a == 3           ~ "Monthly",
      depfreq_a == 4           ~ "A few times a year",
      depfreq_a == 5           ~ "Never",
      TRUE                     ~ NA_character_
    ) |> factor(levels = c("Never","A few times a year","Monthly","Daily/Weekly")),
    anx_freq = case_when(
      anxfreq_a %in% c(1,2)   ~ "Daily/Weekly",
      anxfreq_a == 3           ~ "Monthly",
      anxfreq_a == 4           ~ "A few times a year",
      anxfreq_a == 5           ~ "Never",
      TRUE                     ~ NA_character_
    ) |> factor(levels = c("Never","A few times a year","Monthly","Daily/Weekly")),
    income_cat = case_when(
      ratcat_a %in% c(1,2,3,4,5) ~ "<2.00",
      ratcat_a %in% c(6,7,8,9)   ~ "2.00-3.99",
      ratcat_a %in% c(10,11)     ~ "4.00-4.99",
      ratcat_a %in% c(12,13,14)  ~ "5.00+",
      TRUE                        ~ NA_character_
    ) |> factor(levels = c("<2.00","2.00-3.99","4.00-4.99","5.00+")),
    has_code95 = (modfreqw_raw == 95 | vigfreqw_raw == 95 | strfreqw_raw == 95),
    has_code95 = ifelse(is.na(has_code95), FALSE, has_code95)
  )

d <- d |>
  filter(!is.na(wtfa_a), !is.na(pstrat), !is.na(ppsu)) |>
  mutate(
    hicp1    = as.numeric(hicp == 1),
    pa_score = case_when(
      pa_cat == "No PA" ~ 0,
      pa_cat == "Some PA but below guidelines" ~ 1,
      pa_cat == "Meets aerobic only" ~ 2,
      pa_cat == "Meets aerobic + strength" ~ 3
    ),
    eq_capped = pmin(eq_min_wk, 1500),
    str_days  = strfreqw_a,
    age_grp   = case_when(
      agep_a < 40                ~ "18-39",
      agep_a >= 40 & agep_a < 65 ~ "40-64",
      agep_a >= 65               ~ "65+"
    ) |> factor(levels = c("18-39","40-64","65+")),
    in_base = hicp_window & !is.na(hicp) & !is.na(pa_cat) &
      !is.na(agep_a) & !is.na(sex_a) & !is.na(race_eth),
    in_exp  = in_base & !is.na(educ_cat) & !is.na(bmi_cat) & !is.na(arthritis) &
      !is.na(dep_freq) & !is.na(anx_freq) & !is.na(income_cat),
    in_no95 = in_base & !has_code95
  )

d_cc      <- d |> filter(in_base)
d_cc_exp  <- d |> filter(in_exp)
d_cc_no95 <- d |> filter(in_no95)

message("Jul–Dec sample adults = ", sum(d$hicp_window))
message("Base analytic n = ", nrow(d_cc))
message("Expanded analytic n = ", nrow(d_cc_exp))
message("Sensitivity (no code 95) n = ", nrow(d_cc_no95))
message("Code 95 respondents removed = ", sum(d_cc$has_code95))

options(survey.lonely.psu = "adjust")
des_full <- svydesign(ids = ~ppsu, strata = ~pstrat, weights = ~wtfa_a,
                      nest = TRUE, data = d)
des_base <- subset(des_full, in_base)  # domains of full design; annual WTFA_A, so totals are ~2x low
des_exp  <- subset(des_full, in_exp)
des_no95 <- subset(des_full, in_no95)

fmt_p <- function(p) {
  p <- as.numeric(p)
  ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
}

message("\n===== S1: Missingness Comparison =====")

d_full <- d |>
  filter(hicp_window) |>
  mutate(inclusion = ifelse(in_base, "Included", "Excluded"))

miss_vars <- c("agep_a")
miss_factor <- c("sex_a","race_eth")

d_full <- d_full |>
  mutate(hicp_avail = ifelse(!is.na(hicp), as.character(hicp), NA_character_))

tab_miss <- CreateTableOne(
  vars    = c("agep_a","sex_a","race_eth","hicp_avail"),
  strata  = "inclusion",
  data    = d_full,
  factorVars = c("sex_a","race_eth","hicp_avail"),
  includeNA  = FALSE,
  test       = TRUE,
  addOverall = FALSE
)

mat_miss <- print(tab_miss, showAllLevels = TRUE, quote = FALSE,
                  noSpaces = TRUE, printToggle = FALSE)
df_miss  <- data.frame(Characteristic = rownames(mat_miss), mat_miss,
                        row.names = NULL, check.names = FALSE)

df_miss$Characteristic <- gsub("agep_a \\(mean \\(SD\\)\\)", "Age, mean (SD)", df_miss$Characteristic)
df_miss$Characteristic <- gsub("sex_a \\(%\\)", "Sex (%)", df_miss$Characteristic)
df_miss$Characteristic <- gsub("race_eth \\(%\\)", "Race/Ethnicity (%)", df_miss$Characteristic)
df_miss$Characteristic <- gsub("hicp_avail \\(%\\)", "HICP status (%)", df_miss$Characteristic)

if ("p" %in% colnames(df_miss)) {
  df_miss$p <- sapply(df_miss$p, function(x) {
    if (is.na(x) || x == "") return("")
    val <- suppressWarnings(as.numeric(x))
    if (is.na(val)) {
      if (grepl("^<", x)) return(x)
      return(x)
    }
    fmt_p(val)
  })
}

write_csv(df_miss, here("Output", "supplement_S1_missingness_comparison.csv"))
message("S1 saved: supplement_S1_missingness_comparison.csv")

message("\n===== S2: Two-Model Confounding Assessment =====")

emm_options(df = survey::degf(des_exp))

fit_A <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth,
                design = des_exp, family = gaussian())
emm_A  <- summary(emmeans(fit_A, ~ pa_cat))

trend_A <- svyglm(hicp1 ~ pa_score + agep_a + sex_a + race_eth,
                  design = des_exp, family = gaussian())
trend_A_test <- survey::regTermTest(trend_A, ~ pa_score)
b_A  <- coef(summary(trend_A))["pa_score","Estimate"]
se_A <- coef(summary(trend_A))["pa_score","Std. Error"]

fit_B <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth +
                  educ_cat + bmi_cat + arthritis + dep_freq + anx_freq + income_cat,
                design = des_exp, family = gaussian())
emm_B  <- summary(emmeans(fit_B, ~ pa_cat,
                           nuisance = c("educ_cat","bmi_cat","dep_freq","anx_freq","income_cat")))

trend_B <- svyglm(hicp1 ~ pa_score + agep_a + sex_a + race_eth +
                    educ_cat + bmi_cat + arthritis + dep_freq + anx_freq + income_cat,
                  design = des_exp, family = gaussian())
trend_B_test <- survey::regTermTest(trend_B, ~ pa_score)
b_B  <- coef(summary(trend_B))["pa_score","Estimate"]
se_B <- coef(summary(trend_B))["pa_score","Std. Error"]

s2_tbl <- tibble(
  `PA category` = as.character(emm_A$pa_cat),
  `Model A: Adj. HICP %` = sprintf("%.1f", 100 * emm_A$emmean),
  `Model A: 95% CI`      = sprintf("%.1f-%.1f", 100 * emm_A$lower.CL, 100 * emm_A$upper.CL),
  `Model B: Adj. HICP %` = sprintf("%.1f", 100 * emm_B$emmean),
  `Model B: 95% CI`      = sprintf("%.1f-%.1f", 100 * emm_B$lower.CL, 100 * emm_B$upper.CL)
)

write_csv(s2_tbl, here("Output", "supplement_S2_two_model_comparison.csv"))

s2_trend <- tibble(
  Model = c("A (demographics)", "B (fully adjusted)"),
  `Per-category change (pp)` = sprintf("%.2f", c(100*b_A, 100*b_B)),
  `95% CI`                   = sprintf("%.2f to %.2f",
                                        c(100*(b_A - 1.96*se_A), 100*(b_B - 1.96*se_B)),
                                        c(100*(b_A + 1.96*se_A), 100*(b_B + 1.96*se_B))),
  `Trend p-value`            = c(
    fmt_p(trend_A_test$p),
    fmt_p(trend_B_test$p)
  )
)

write_csv(s2_trend, here("Output", "supplement_S2_trend_comparison.csv"))
message("S2 saved: supplement_S2_two_model_comparison.csv, supplement_S2_trend_comparison.csv")

attenuation_pct <- 100 * (1 - abs(b_B) / abs(b_A))
message(sprintf("Attenuation from Model A to B: %.1f%%", attenuation_pct))

message("\n===== S3: Sensitivity Excluding Code 95 =====")

emm_options(df = survey::degf(des_no95))

fit_no95 <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth,
                   design = des_no95, family = gaussian())
emm_no95 <- summary(emmeans(fit_no95, ~ pa_cat))

trend_no95 <- svyglm(hicp1 ~ pa_score + agep_a + sex_a + race_eth,
                     design = des_no95, family = gaussian())
trend_no95_test <- survey::regTermTest(trend_no95, ~ pa_score)
b_no95  <- coef(summary(trend_no95))["pa_score","Estimate"]
se_no95 <- coef(summary(trend_no95))["pa_score","Std. Error"]

n_code95_hicp <- sum(d_cc$has_code95 & d_cc$hicp1 == 1, na.rm = TRUE)
n_code95_total <- sum(d_cc$has_code95)
message(sprintf("Code 95 respondents: n = %d; of these, %d (%.1f%%) had HICP",
                n_code95_total, n_code95_hicp, 100*n_code95_hicp/n_code95_total))

s3_tbl <- tibble(
  `PA category`         = as.character(emm_no95$pa_cat),
  `Adj. HICP % (no 95)` = sprintf("%.1f", 100 * emm_no95$emmean),
  `95% CI`              = sprintf("%.1f-%.1f", 100 * emm_no95$lower.CL, 100 * emm_no95$upper.CL)
)

s3_trend <- tibble(
  Analysis = "Excluding code 95",
  `Per-category change (pp)` = sprintf("%.2f", 100*b_no95),
  `95% CI` = sprintf("%.2f to %.2f", 100*(b_no95 - 1.96*se_no95), 100*(b_no95 + 1.96*se_no95)),
  `Trend p-value` = fmt_p(trend_no95_test$p),
  `n removed` = n_code95_total,
  `Analytic n` = nrow(d_cc_no95)
)

write_csv(s3_tbl, here("Output", "supplement_S3_sensitivity_no95_prevalences.csv"))
write_csv(s3_trend, here("Output", "supplement_S3_sensitivity_no95_trend.csv"))
message("S3 saved.")

message("\n===== S4: Subgroup Interaction Analyses =====")

fit_age_int <- svyglm(hicp1 ~ pa_score * age_grp + sex_a + race_eth,
                      design = des_base, family = gaussian())
age_int_test <- survey::regTermTest(fit_age_int, ~ pa_score:age_grp)

age_slopes <- tibble()
for (ag in levels(d_cc$age_grp)) {
  des_sub <- subset(des_base, age_grp == ag)
  fit_sub <- svyglm(hicp1 ~ pa_score + sex_a + race_eth,
                    design = des_sub, family = gaussian())
  b  <- coef(summary(fit_sub))["pa_score","Estimate"]
  se <- coef(summary(fit_sub))["pa_score","Std. Error"]
  p  <- coef(summary(fit_sub))["pa_score","Pr(>|t|)"]
  age_slopes <- bind_rows(age_slopes, tibble(
    Subgroup = ag,
    `Per-category change (pp)` = sprintf("%.2f", 100*b),
    `95% CI` = sprintf("%.2f to %.2f", 100*(b - 1.96*se), 100*(b + 1.96*se)),
    `p-value` = fmt_p(p)
  ))
}

fit_sex_int <- svyglm(hicp1 ~ pa_score * sex_a + agep_a + race_eth,
                      design = des_base, family = gaussian())
sex_int_test <- survey::regTermTest(fit_sex_int, ~ pa_score:sex_a)

sex_slopes <- tibble()
for (sx in levels(d_cc$sex_a)) {
  des_sub <- subset(des_base, sex_a == sx)
  fit_sub <- svyglm(hicp1 ~ pa_score + agep_a + race_eth,
                    design = des_sub, family = gaussian())
  b  <- coef(summary(fit_sub))["pa_score","Estimate"]
  se <- coef(summary(fit_sub))["pa_score","Std. Error"]
  p  <- coef(summary(fit_sub))["pa_score","Pr(>|t|)"]
  sex_slopes <- bind_rows(sex_slopes, tibble(
    Subgroup = sx,
    `Per-category change (pp)` = sprintf("%.2f", 100*b),
    `95% CI` = sprintf("%.2f to %.2f", 100*(b - 1.96*se), 100*(b + 1.96*se)),
    `p-value` = fmt_p(p)
  ))
}

fit_race_int <- svyglm(hicp1 ~ pa_score * race_eth + agep_a + sex_a,
                       design = des_base, family = gaussian())
race_int_test <- survey::regTermTest(fit_race_int, ~ pa_score:race_eth)

race_slopes <- tibble()
for (rc in levels(d_cc$race_eth)) {
  if (sum(d_cc$race_eth == rc) < 100) next
  des_sub <- subset(des_base, race_eth == rc)
  fit_sub <- svyglm(hicp1 ~ pa_score + agep_a + sex_a,
                    design = des_sub, family = gaussian())
  b  <- coef(summary(fit_sub))["pa_score","Estimate"]
  se <- coef(summary(fit_sub))["pa_score","Std. Error"]
  p  <- coef(summary(fit_sub))["pa_score","Pr(>|t|)"]
  race_slopes <- bind_rows(race_slopes, tibble(
    Subgroup = rc,
    `Per-category change (pp)` = sprintf("%.2f", 100*b),
    `95% CI` = sprintf("%.2f to %.2f", 100*(b - 1.96*se), 100*(b + 1.96*se)),
    `p-value` = fmt_p(p)
  ))
}

s4_slopes <- bind_rows(
  age_slopes  |> mutate(Modifier = "Age group", .before = 1),
  sex_slopes  |> mutate(Modifier = "Sex", .before = 1),
  race_slopes |> mutate(Modifier = "Race/Ethnicity", .before = 1)
)

s4_interactions <- tibble(
  Modifier         = c("Age group","Sex","Race/Ethnicity"),
  `Interaction p`  = c(
    fmt_p(age_int_test$p),
    fmt_p(sex_int_test$p),
    fmt_p(race_int_test$p)
  )
)

write_csv(s4_slopes, here("Output", "supplement_S4_subgroup_slopes.csv"))
write_csv(s4_interactions, here("Output", "supplement_S4_interaction_tests.csv"))
message("S4 saved.")

message("\n===== S5: LPM Bounds Check =====")

fit_lpm_base <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth,
                       design = des_base, family = gaussian())
preds_lpm <- predict(fit_lpm_base, type = "response")
pred_vals <- as.numeric(preds_lpm)

s5_bounds <- tibble(
  Metric = c("Min predicted probability",
             "Max predicted probability",
             "N predictions < 0",
             "N predictions > 1",
             "Proportion in [0,1]"),
  Value  = c(
    sprintf("%.4f", min(pred_vals, na.rm = TRUE)),
    sprintf("%.4f", max(pred_vals, na.rm = TRUE)),
    as.character(sum(pred_vals < 0, na.rm = TRUE)),
    as.character(sum(pred_vals > 1, na.rm = TRUE)),
    sprintf("%.4f", mean(pred_vals >= 0 & pred_vals <= 1, na.rm = TRUE))
  )
)

write_csv(s5_bounds, here("Output", "supplement_S5_lpm_bounds.csv"))
message("S5 saved.")

message("\n===== S6: Spline Knot Locations =====")

ns_obj    <- ns(d_cc$eq_capped, df = 3)
int_knots <- attr(ns_obj, "knots")
bnd_knots <- attr(ns_obj, "Boundary.knots")

s6_knots <- tibble(
  `Knot type` = c("Boundary (lower)", "Internal (33rd percentile)",
                  "Internal (67th percentile)", "Boundary (upper)"),
  `Moderate-eq. min/wk` = c(bnd_knots[1], int_knots, bnd_knots[2])
)

write_csv(s6_knots, here("Output", "supplement_S6_spline_knots.csv"))
message("S6 saved.")

message("\n===== S7: Cap at 1500 min/wk Justification =====")

n_above_1500 <- sum(d_cc$eq_min_wk > 1500, na.rm = TRUE)
n_above_1000 <- sum(d_cc$eq_min_wk > 1000, na.rm = TRUE)
n_total_pa   <- sum(!is.na(d_cc$eq_min_wk))
max_val      <- max(d_cc$eq_min_wk, na.rm = TRUE)
p99          <- quantile(d_cc$eq_min_wk, 0.99, na.rm = TRUE)
p995         <- quantile(d_cc$eq_min_wk, 0.995, na.rm = TRUE)

s7_cap <- tibble(
  Metric = c("Total with valid moderate-eq. min/wk",
             "Maximum (min/wk)",
             "99th percentile (min/wk)",
             "99.5th percentile (min/wk)",
             "N > 1,000 min/wk",
             "% > 1,000 min/wk",
             "N > 1,500 min/wk",
             "% > 1,500 min/wk"),
  Value  = c(
    as.character(n_total_pa),
    as.character(max_val),
    sprintf("%.1f", p99),
    sprintf("%.1f", p995),
    as.character(n_above_1000),
    sprintf("%.2f%%", 100*n_above_1000/n_total_pa),
    as.character(n_above_1500),
    sprintf("%.2f%%", 100*n_above_1500/n_total_pa)
  )
)

write_csv(s7_cap, here("Output", "supplement_S7_cap_justification.csv"))
message("S7 saved.")

message("\n===== Prevalence Ratios (Modified Poisson) =====")

fit_pr <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth,
                 design = des_base, family = quasipoisson(link = "log"))

pr_coefs <- coef(summary(fit_pr))
pr_names <- rownames(pr_coefs)
pa_rows  <- grep("^pa_cat", pr_names)

pr_tbl <- tibble(
  `PA category` = c("No PA (ref)", sub("^pa_cat", "", pr_names[pa_rows])),
  PR            = c("1.00", sprintf("%.2f", exp(pr_coefs[pa_rows, "Estimate"]))),
  `95% CI`      = c("ref", sprintf("%.2f-%.2f",
                     exp(pr_coefs[pa_rows, "Estimate"] - 1.96 * pr_coefs[pa_rows, "Std. Error"]),
                     exp(pr_coefs[pa_rows, "Estimate"] + 1.96 * pr_coefs[pa_rows, "Std. Error"])))
)

write_csv(pr_tbl, here("Output", "prevalence_ratios_model_A.csv"))
message("Prevalence ratios (Model A) saved.")

fit_pr_B <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth +
                     educ_cat + bmi_cat + arthritis + dep_freq + anx_freq + income_cat,
                   design = des_exp, family = quasipoisson(link = "log"))

pr_coefs_B <- coef(summary(fit_pr_B))
pr_names_B <- rownames(pr_coefs_B)
pa_rows_B  <- grep("^pa_cat", pr_names_B)

pr_tbl_B <- tibble(
  `PA category` = c("No PA (ref)", sub("^pa_cat", "", pr_names_B[pa_rows_B])),
  PR            = c("1.00", sprintf("%.2f", exp(pr_coefs_B[pa_rows_B, "Estimate"]))),
  `95% CI`      = c("ref", sprintf("%.2f-%.2f",
                     exp(pr_coefs_B[pa_rows_B, "Estimate"] - 1.96 * pr_coefs_B[pa_rows_B, "Std. Error"]),
                     exp(pr_coefs_B[pa_rows_B, "Estimate"] + 1.96 * pr_coefs_B[pa_rows_B, "Std. Error"])))
)

write_csv(pr_tbl_B, here("Output", "prevalence_ratios_model_B.csv"))
message("Prevalence ratios (Model B) saved.")

message("\n========== ALL SUPPLEMENTAL ANALYSES COMPLETE ==========")
