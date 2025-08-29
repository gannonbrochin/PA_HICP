
# NHIS 2020: PA vs HICP — Trend test, dose–response figure, adjusted rates, and TableOne
# v14 (fix TableOne header rename; no education; outputs to working directory)
#
# Files written:
#   - trend_test_2020.txt
#   - fig_dose_response_hicp_pa_2020.png
#   - adjusted_rates_by_pa_2020.csv
#   - adjusted_rates_by_pa_2020.docx
#   - tableone_by_pa_2020.csv
#   - tableone_by_pa_2020.docx
#
# install.packages(c("tidyverse","survey","broom","janitor","readr",
#                    "emmeans","splines","ggplot2","officer","flextable","forcats","tableone"))

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
})

# ---------- paths ----------
data_path <- "adult20.csv"   # set to your 2020 Sample Adult CSV

# ---------- load & recode ----------
d <- read_csv(data_path, show_col_types = FALSE) |> clean_names()

need <- c("wtfa_a","pstrat","ppsu","srvy_yr","agep_a","sex_a","hisp_a","raceallp_a",
          "paifrq3m_a","paiwklm3m_a","modfreqw_a","modmin_a","vigfreqw_a","vigmin_a","strfreqw_a","educ_a")
has  <- need %in% names(d)
if (!all(has)) stop("Missing variables in 2020 file: ", paste(need[!has], collapse = ", "))

# Outcome: HICP
d <- d |>
  mutate(
    chronic_pain = case_when(paifrq3m_a %in% c(3,4) ~ 1L,
                             paifrq3m_a %in% c(1,2) ~ 0L,
                             TRUE ~ NA_integer_),
    hicp = case_when(chronic_pain==1L & paiwklm3m_a %in% c(3,4) ~ 1L,
                     !is.na(chronic_pain) & paiwklm3m_a %in% c(1,2) ~ 0L,
                     TRUE ~ NA_integer_)
  )

# Exposure: PA categories and continuous dose
clean_val <- function(x) ifelse(x %in% c(7,8,9,97,98,99), NA, x)
d <- d |>
  mutate(across(c(modfreqw_a, modmin_a, vigfreqw_a, vigmin_a, strfreqw_a), as.numeric),
         across(c(modfreqw_a, modmin_a, vigfreqw_a, vigmin_a, strfreqw_a), clean_val),
         modmin_a = ifelse(modmin_a >= 997, NA, modmin_a),
         vigmin_a = ifelse(vigmin_a >= 997, NA, vigmin_a),
         mod_min_wk = modfreqw_a * modmin_a,
         vig_min_wk = vigfreqw_a * vigmin_a,
         eq_min_wk  = mod_min_wk + 2*vig_min_wk,
         meet_aerobic  = ifelse(!is.na(eq_min_wk) & eq_min_wk >= 150, 1L, ifelse(!is.na(eq_min_wk), 0L, NA)),
         meet_strength = ifelse(!is.na(strfreqw_a) & strfreqw_a >= 2, 1L, ifelse(!is.na(strfreqw_a), 0L, NA)),
         pa_cat = case_when(
           !is.na(eq_min_wk) & eq_min_wk==0 ~ "No PA",
           !is.na(meet_aerobic) & meet_aerobic==1 & !is.na(meet_strength) & meet_strength==1 ~ "Meets aerobic + strength",
           !is.na(meet_aerobic) & meet_aerobic==1 & !is.na(meet_strength) & meet_strength==0 ~ "Meets aerobic only",
           !is.na(eq_min_wk) & eq_min_wk>0 ~ "Some PA but below guidelines",
           TRUE ~ NA_character_
         ),
         pa_cat = factor(pa_cat, levels = c("No PA","Some PA but below guidelines","Meets aerobic only","Meets aerobic + strength"))
  )

# Covariates (sex, race/eth only — education removed)
d <- d |>
  mutate(
    sex_a = factor(sex_a, levels=c(1,2), labels=c("Male","Female")),
    is_hisp = case_when(hisp_a==1 ~ 1L, hisp_a==2 ~ 0L, TRUE ~ NA_integer_),
    race_eth = case_when(
      is_hisp==1 ~ "Hispanic",
      is_hisp==0 & raceallp_a==1 ~ "NH White",
      is_hisp==0 & raceallp_a==2 ~ "NH Black",
      is_hisp==0 & raceallp_a==4 ~ "NH Asian",
      is_hisp==0 ~ "NH Other/Multiracial",
      TRUE ~ NA_character_
    ) |> factor(levels=c("NH White","NH Black","NH Asian","Hispanic","NH Other/Multiracial"))
  )

# Analytic complete-case data (education omitted from filter)
d_cc <- d |>
  filter(!is.na(wtfa_a), !is.na(pstrat), !is.na(ppsu),
         !is.na(hicp), !is.na(pa_cat),
         !is.na(agep_a), !is.na(sex_a), !is.na(race_eth)) |>
  mutate(
    hicp1 = as.numeric(hicp == 1),
    pa_cat = fct_drop(pa_cat),
    sex_a = fct_drop(sex_a),
    race_eth = fct_drop(race_eth)
  )

message("Analytic n = ", nrow(d_cc), "; HICP events = ", sum(d_cc$hicp==1, na.rm=TRUE))

# Survey design
options(survey.lonely.psu = "adjust")
des <- svydesign(ids=~ppsu, strata=~pstrat, weights=~wtfa_a, nest=TRUE, data=d_cc)

# ---------- (1) Ordered trend test across PA (NO education) ----------
d_cc <- d_cc |>
  mutate(
    pa_score = case_when(
      pa_cat == "No PA" ~ 0,
      pa_cat == "Some PA but below guidelines" ~ 1,
      pa_cat == "Meets aerobic only" ~ 2,
      pa_cat == "Meets aerobic + strength" ~ 3,
      TRUE ~ NA_real_
    )
  )
des <- update(des, pa_score = d_cc$pa_score)

fit_trend <- svyglm(hicp1 ~ pa_score + agep_a + sex_a + race_eth,
                    design = des, family = gaussian())
trend_test <- survey::regTermTest(fit_trend, ~ pa_score)

# Append per-category and 0→3 differences (percentage points) to the trend text
b  <- coef(summary(fit_trend))["pa_score","Estimate"]
se <- coef(summary(fit_trend))["pa_score","Std. Error"]
est_pp <- 100*b
ci_pp  <- 100*c(b - 1.96*se, b + 1.96*se)
est_pp_0to3 <- 100*(3*b)
ci_pp_0to3  <- 100*c(3*b - 1.96*(3*se), 3*b + 1.96*(3*se))

sink("trend_test_2020.txt")
cat("Design-based linear trend test of HICP across PA score (0..3) — no education in model\n")
print(trend_test)
cat("\nPer-category change (one-level increase in PA):\n")
cat(sprintf("Estimate = %.2f pp; 95%% CI = [%.2f, %.2f] pp\n", est_pp, ci_pp[1], ci_pp[2]))
cat("End-to-end change (No PA -> Meets aerobic + strength; 3 levels):\n")
cat(sprintf("Estimate = %.2f pp; 95%% CI = [%.2f, %.2f] pp\n", est_pp_0to3, ci_pp_0to3[1], ci_pp_0to3[2]))
sink()

# ---------- (2) Dose–response spline model and figure (NO education) ----------
d_cc <- d_cc |>
  mutate(eq_capped = pmin(eq_min_wk, 1500),
         str_days = strfreqw_a)

des <- update(des, eq_capped = d_cc$eq_capped, str_days = d_cc$str_days)

fit_spline <- svyglm(hicp1 ~ ns(eq_capped, df = 3) + str_days +
                       agep_a + sex_a + race_eth,
                     design = des, family = quasipoisson(link = "log"))

# Prediction grid (set covariates at typical/reference values)
grid <- expand.grid(eq_capped = seq(0, 1000, by = 25))
grid$str_days <- stats::median(d_cc$str_days, na.rm=TRUE)
grid$agep_a   <- stats::median(d_cc$agep_a, na.rm=TRUE)
grid$sex_a    <- factor(levels(d_cc$sex_a)[1], levels=levels(d_cc$sex_a))
grid$race_eth <- factor(levels(d_cc$race_eth)[1], levels=levels(d_cc$race_eth))

# Robust prediction with SEs
pred <- try(predict(fit_spline, newdata = grid, type = "link", se.fit = TRUE), silent = TRUE)
if (!inherits(pred, "try-error") && is.list(pred) && all(c("fit","se.fit") %in% names(pred))) {
  eta <- as.numeric(pred$fit)
  se  <- as.numeric(pred$se.fit)
} else {
  TT <- delete.response(terms(fit_spline))
  X  <- model.matrix(TT, grid)
  beta <- coef(fit_spline); V <- vcov(fit_spline)
  eta  <- as.vector(X %*% beta)
  se   <- sqrt(rowSums((X %*% V) * X))
}

grid$prev  <- exp(eta)
grid$lower <- exp(eta - 1.96*se)
grid$upper <- exp(eta + 1.96*se)

p <- ggplot(grid, aes(x = eq_capped, y = 100*prev)) +
  geom_ribbon(aes(ymin = 100*lower, ymax = 100*upper), alpha = 0.2) +
  geom_line(size = 1) +
  labs(x = "Moderate-equivalent minutes/week (capped at 1500)",
       y = "Adjusted HICP prevalence (%)",
       title = "Adjusted dose–response of HICP across leisure-time PA (NHIS 2020)") +
  theme_minimal(base_size = 12)
ggsave("fig_dose_response_hicp_pa_2020.png", p, width = 7, height = 5, dpi = 300)

# ---------- (3) Adjusted marginal rates by PA category (NO education) ----------
emm_options(df = survey::degf(des))
fit_lpm <- svyglm(hicp1 ~ pa_cat + agep_a + sex_a + race_eth,
                  design = des, family = gaussian())
emm_pa <- emmeans(fit_lpm, ~ pa_cat)
emm_tbl <- summary(emm_pa) |>
  as_tibble() |>
  transmute(
    `PA category` = pa_cat,
    `Adj. HICP %` = 100*emmean,
    `95% CI` = sprintf("%.1f–%.1f", 100*lower.CL, 100*upper.CL)
  )

write_csv(emm_tbl, "adjusted_rates_by_pa_2020.csv")

ft <- regulartable(emm_tbl) |> autofit()
doc <- read_docx()
doc <- body_add_par(doc, "Adjusted HICP prevalence by PA category (NHIS 2020) — no education in model", style = "table title")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "Notes: Design-based linear model (identity link) adjusted for age, sex, and race/ethnicity; education omitted.", style = "Normal")
print(doc, target = "adjusted_rates_by_pa_2020.docx")

# ---------- (4) TableOne (survey) WITHOUT education, analytic sample only; safe column renaming ----------
vars <- c("agep_a","sex_a","race_eth")
factorVars <- c("sex_a","race_eth")

tab1 <- svyCreateTableOne(vars = vars, strata = "pa_cat", data = des,
                          factorVars = factorVars, includeNA = FALSE, test = TRUE, addOverall = TRUE)

mat <- print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# Convert to data.frame
df_tab1 <- data.frame(Characteristic = rownames(mat), mat, row.names = NULL, check.names = FALSE)

# UNWEIGHTED Ns for headers
n_overall <- nrow(d_cc)
pa_levels <- levels(d_cc$pa_cat)
n_by_pa   <- as.integer(table(d_cc$pa_cat)[pa_levels])

# Safe rename: map existing column names to desired headers, preserving unmatched (e.g., 'p')
orig_cols <- colnames(mat)
new_cols  <- orig_cols

# Replace 'Overall' header if present
idx_overall <- which(orig_cols == "Overall")
if (length(idx_overall) == 1) {
  new_cols[idx_overall] <- paste0("Overall (N=", format(n_overall, big.mark=","), ")")
}

# Replace each PA level header if present
for (i in seq_along(pa_levels)) {
  lvl <- pa_levels[i]
  idx <- which(orig_cols == lvl)
  if (length(idx) == 1) {
    new_cols[idx] <- paste0(lvl, " (n=", format(n_by_pa[i], big.mark=","), ")")
  }
}

# Now set final names: Characteristic + mapped group names (plus any extra columns like 'p')
final_names <- c("Characteristic", new_cols)
stopifnot(length(final_names) == ncol(df_tab1))
stopifnot(!any(is.na(final_names)))
stopifnot(anyDuplicated(final_names) == 0)

colnames(df_tab1) <- final_names

# Save CSV and DOCX
write_csv(df_tab1, "tableone_by_pa_2020.csv")

ft_tab1 <- regulartable(df_tab1) |> autofit()
doc2 <- read_docx()
doc2 <- body_add_par(doc2, "TableOne (survey-weighted) by PA category — NHIS 2020 (no education)", style = "table title")
doc2 <- body_add_flextable(doc2, ft_tab1)

# Add note: Ns are unweighted for analytic sample; include total 2020 respondents (unweighted)
total_2020_unw <- nrow(d)
note_txt <- paste0(
  "Notes: Values are survey-weighted (WTFA_A) using the complete-case analytic sample (",
  format(n_overall, big.mark=","), " adults). Column headers show UNWEIGHTED counts. ",
  "For reference, total 2020 Sample Adult respondents (unweighted): ",
  format(total_2020_unw, big.mark=","), "."
)
doc2 <- body_add_par(doc2, note_txt, style = "Normal")

print(doc2, target = "tableone_by_pa_2020.docx")

message("✓ Wrote: ",
        "\n  - trend_test_2020.txt",
        "\n  - fig_dose_response_hicp_pa_2020.png",
        "\n  - adjusted_rates_by_pa_2020.csv",
        "\n  - adjusted_rates_by_pa_2020.docx",
        "\n  - tableone_by_pa_2020.csv",
        "\n  - tableone_by_pa_2020.docx")
