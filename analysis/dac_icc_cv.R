#!/usr/bin/env Rscript
# dac_icc_cv.R -- ICC and within-patient CV analysis of DAC CRF longitudinal biomarker data
#
# Reads dac_long.csv and dac_baseline.csv produced by dac_pull.py.
# Outputs tables and figures to results/dac_icc_cv/.
#
# Usage:
#   Rscript analysis/dac_icc_cv.R [--data-dir PATH] [--output-dir PATH]
#
# Methods:
#   ICC  -- LMM-based (lmer null model): ICC = patient variance / total variance.
#           Handles unbalanced visits. 95% CI via profile likelihood on variance components.
#           Sensitivity: classic ICC (psych::ICC) on V1-V2 pairs only.
#   CV   -- per-patient SD/mean x 100 (%) across all available visits.
#           Predictor analysis: correlation with baseline covariates.
#
# Requires: tidyverse, lme4, psych, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(patchwork)
  library(lme4)
  library(psych)
  library(gt)
})

# -- Paths ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) normalizePath(file.path(getwd(), "analysis"), mustWork = FALSE)
)
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

data_dir <- file.path(project_dir, "data")
out_dir  <- file.path(project_dir, "results", "dac_icc_cv")

for (i in seq_along(args)) {
  if (args[i] == "--data-dir"    && i < length(args)) data_dir <- args[i + 1]
  if (args[i] == "--output-dir"  && i < length(args)) out_dir  <- args[i + 1]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "theme_research.R"))
source(file.path(script_dir, "palettes.R"))

message("Reading data from: ", data_dir)
message("Writing results to: ", out_dir)

# -- Load data -----------------------------------------------------------------

long <- read_csv(file.path(data_dir, "dac_long.csv"), show_col_types = FALSE) |>
  mutate(
    visit = factor(visit, levels = c("V1", "V2", "V3", "V4")),
    visit_date = as.Date(visit_date),
    blood_date = as.Date(blood_date)
  )

baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("Female", "Male")))

# Longitudinal subcohort: >= 2 visits
long_ids <- long |>
  group_by(record_id) |>
  filter(n_distinct(visit) >= 2) |>
  pull(record_id) |>
  unique()

long_sub <- long |> filter(record_id %in% long_ids)

message(sprintf(
  "Longitudinal subcohort: %d patients, %d rows",
  n_distinct(long_sub$record_id), nrow(long_sub)
))

# -- Biomarker metadata --------------------------------------------------------

biomarkers <- list(
  ptau217  = list(label = "pTau217",          unit = "pg/mL",   log_scale = TRUE),
  nfl      = list(label = "NfL",              unit = "pg/mL",   log_scale = TRUE),
  gfap     = list(label = "GFAP",             unit = "pg/mL",   log_scale = TRUE),
  ab42_40  = list(label = "Ab42/40",          unit = "ratio",   log_scale = FALSE),
  lucent_ad = list(label = "Lucent AD",       unit = "score",   log_scale = FALSE)
)

# -- Helper: LMM-based ICC with 95% CI -----------------------------------------
# ICC(1) from one-way random effects model:
#   value ~ 1 + (1 | record_id)
#   ICC = var_patient / (var_patient + var_residual)
# 95% CI derived from profile likelihood confidence intervals on variance components.

lmm_icc <- function(df, biomarker_col) {
  dat <- df |>
    select(record_id, value = all_of(biomarker_col)) |>
    filter(!is.na(value))

  n_patients <- n_distinct(dat$record_id)
  n_obs      <- nrow(dat)

  # Need >=2 patients and >=2 total obs
  if (n_patients < 2 || n_obs < 4) {
    return(tibble(
      biomarker = biomarker_col, n_patients = n_patients, n_obs = n_obs,
      icc = NA_real_, icc_lo = NA_real_, icc_hi = NA_real_,
      var_patient = NA_real_, var_residual = NA_real_
    ))
  }

  fit <- lmer(value ~ 1 + (1 | record_id), data = dat, REML = TRUE)

  vc          <- as.data.frame(VarCorr(fit))
  var_patient <- vc$vcov[vc$grp == "record_id"]
  var_resid   <- vc$vcov[vc$grp == "Residual"]
  icc_est     <- var_patient / (var_patient + var_resid)

  # 95% CI via parametric bootstrap on ICC
  # Use boot_icc() below; fall back to NA on error
  ci <- tryCatch({
    boot_out <- bootMer(
      fit,
      FUN = function(m) {
        vc_b <- as.data.frame(VarCorr(m))
        vp   <- vc_b$vcov[vc_b$grp == "record_id"]
        vr   <- vc_b$vcov[vc_b$grp == "Residual"]
        vp / (vp + vr)
      },
      nsim = 500,
      use.u = FALSE,
      seed = 42
    )
    quantile(boot_out$t, probs = c(0.025, 0.975), na.rm = TRUE)
  }, error = function(e) {
    message("  Bootstrap failed for ", biomarker_col, ": ", conditionMessage(e))
    c(NA_real_, NA_real_)
  })

  tibble(
    biomarker   = biomarker_col,
    n_patients  = n_patients,
    n_obs       = n_obs,
    icc         = icc_est,
    icc_lo      = ci[[1]],
    icc_hi      = ci[[2]],
    var_patient = var_patient,
    var_residual = var_resid
  )
}

# -- Helper: classic ICC on V1-V2 pairs (sensitivity) -------------------------

classic_icc_v1v2 <- function(df, biomarker_col) {
  wide <- df |>
    filter(visit %in% c("V1", "V2")) |>
    select(record_id, visit, value = all_of(biomarker_col)) |>
    filter(!is.na(value)) |>
    pivot_wider(names_from = visit, values_from = value) |>
    filter(!is.na(V1), !is.na(V2))

  if (nrow(wide) < 3) {
    return(list(n_pairs = nrow(wide), icc3_single = NA, lo = NA, hi = NA))
  }

  mat   <- wide |> select(V1, V2) |> as.matrix()
  icc_r <- ICC(mat, lmer = FALSE)  # psych::ICC

  # ICC3 = two-way mixed, consistency, single measures
  icc3_row <- icc_r$results |> filter(type == "ICC3")

  list(
    n_pairs    = nrow(wide),
    icc3_single = round(icc3_row$ICC,  3),
    lo          = round(icc3_row$`lower bound`, 3),
    hi          = round(icc3_row$`upper bound`, 3)
  )
}

# -- Run LMM ICC for all biomarkers --------------------------------------------

message("\nComputing LMM-based ICC (may take ~30s for bootstrap)...")

icc_results <- map_dfr(names(biomarkers), function(bm) {
  message("  ", biomarkers[[bm]]$label)
  lmm_icc(long_sub, bm)
})

# Interpretation bands
icc_results <- icc_results |>
  mutate(
    label = map_chr(biomarker, ~ biomarkers[[.x]]$label),
    interpretation = case_when(
      icc >= 0.90 ~ "Excellent",
      icc >= 0.75 ~ "Good",
      icc >= 0.50 ~ "Moderate",
      icc >= 0.00 ~ "Poor",
      TRUE        ~ NA_character_
    )
  )

message("\nICC results:")
print(icc_results |> select(label, n_patients, n_obs, icc, icc_lo, icc_hi, interpretation))

# -- Sensitivity: classic ICC on V1-V2 pairs -----------------------------------

message("\nClassic ICC (V1-V2 pairs only)...")

classic_results <- map_dfr(names(biomarkers), function(bm) {
  res <- classic_icc_v1v2(long_sub, bm)
  tibble(
    biomarker   = bm,
    label       = biomarkers[[bm]]$label,
    n_pairs     = res$n_pairs,
    icc3        = res$icc3_single,
    icc3_lo     = res$lo,
    icc3_hi     = res$hi
  )
})

message("\nClassic ICC (V1-V2, ICC3 consistency):")
print(classic_results)

# -- Per-patient CV ------------------------------------------------------------

message("\nComputing per-patient CV...")

cv_long <- long_sub |>
  pivot_longer(
    cols = all_of(names(biomarkers)),
    names_to = "biomarker",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  group_by(record_id, biomarker) |>
  filter(n() >= 2) |>
  summarise(
    n_obs    = n(),
    mean_val = mean(value),
    sd_val   = sd(value),
    cv_pct   = sd_val / mean_val * 100,
    .groups  = "drop"
  ) |>
  mutate(label = map_chr(biomarker, ~ biomarkers[[.x]]$label))

# CV summary per biomarker
cv_summary <- cv_long |>
  group_by(biomarker, label) |>
  summarise(
    n_patients  = n(),
    median_cv   = median(cv_pct),
    q25_cv      = quantile(cv_pct, 0.25),
    q75_cv      = quantile(cv_pct, 0.75),
    pct_above_20 = mean(cv_pct > 20) * 100,
    .groups     = "drop"
  )

message("\nCV summary:")
print(cv_summary)

# -- Predictor analysis for CV -------------------------------------------------
# Join per-patient CV with baseline covariates, compute correlations + regression.

message("\nPreparing CV predictor analysis...")

# Focus on ptau217 and nfl (largest n); run separate models
predictor_vars <- c(
  "age_at_enrol", "sex", "htn", "dm2", "ckd", "egfr",
  "apoe_genotype", "pet_result", "ptau217_risk"
)

# Pull baseline fields + ptau217_risk from V1 visit (per-patient)
v1_ptau_risk <- long_sub |>
  filter(visit == "V1") |>
  select(record_id, ptau217_risk) |>
  distinct()

baseline_pred <- baseline |>
  select(record_id, all_of(intersect(predictor_vars, names(baseline)))) |>
  left_join(v1_ptau_risk, by = "record_id") |>
  mutate(
    sex_f   = as.integer(sex == "Female"),
    htn_f   = as.integer(htn == "Yes"),
    dm2_f   = as.integer(dm2 == "Yes"),
    ckd_f   = as.integer(ckd == "Yes"),
    apoe4   = as.integer(grepl("e4", apoe_genotype)),
    amyloid_pos = case_when(
      pet_result == "Positive" ~ 1L,
      pet_result == "Negative" ~ 0L,
      TRUE ~ NA_integer_
    ),
    ptau_high = as.integer(ptau217_risk == "High")
  )

# Merge CV with baseline predictors
cv_pred <- cv_long |>
  left_join(baseline_pred, by = "record_id")

# Spearman correlation of CV with continuous predictors
cv_cor <- function(cv_df, bm_label) {
  sub <- cv_df |> filter(label == bm_label)
  cont_vars <- c("age_at_enrol", "egfr")
  binary_vars <- c("sex_f", "htn_f", "dm2_f", "ckd_f", "apoe4", "amyloid_pos", "ptau_high")

  bind_rows(
    map_dfr(cont_vars, function(v) {
      vals <- sub[[v]]
      cv   <- sub$cv_pct
      ok   <- !is.na(vals) & !is.na(cv)
      if (sum(ok) < 5) return(tibble(variable = v, n = sum(ok), rho = NA, p = NA))
      ct <- cor.test(cv[ok], vals[ok], method = "spearman", exact = FALSE)
      tibble(variable = v, n = sum(ok), rho = ct$estimate, p = ct$p.value)
    }),
    map_dfr(binary_vars, function(v) {
      vals <- sub[[v]]
      cv   <- sub$cv_pct
      ok   <- !is.na(vals) & !is.na(cv) & vals %in% c(0L, 1L)
      if (sum(ok) < 5) return(tibble(variable = v, n = sum(ok), rho = NA, p = NA))
      # Skip zero-variance predictors (all same value)
      if (var(as.numeric(vals[ok])) == 0) return(tibble(variable = v, n = sum(ok), rho = NA, p = NA))
      ct <- cor.test(cv[ok], as.numeric(vals[ok]), method = "spearman", exact = FALSE)
      tibble(variable = v, n = sum(ok), rho = ct$estimate, p = ct$p.value)
    })
  ) |>
    arrange(p) |>
    mutate(biomarker = bm_label)
}

cor_ptau <- cv_cor(cv_pred, "pTau217")
cor_nfl  <- cv_cor(cv_pred, "NfL")

message("\nCV predictors -- pTau217 (Spearman rho):")
print(cor_ptau)
message("\nCV predictors -- NfL (Spearman rho):")
print(cor_nfl)

# -- Table 1: ICC summary ------------------------------------------------------

message("\nBuilding ICC summary table...")

icc_tbl <- icc_results |>
  select(label, n_patients, n_obs, icc, icc_lo, icc_hi, interpretation) |>
  left_join(
    classic_results |> select(label, n_pairs, icc3, icc3_lo, icc3_hi),
    by = "label"
  ) |>
  mutate(
    icc_fmt      = sprintf("%.2f (%.2f\u2013%.2f)", icc, icc_lo, icc_hi),
    icc3_fmt     = case_when(
      !is.na(icc3) ~ sprintf("%.2f (%.2f\u2013%.2f)", icc3, icc3_lo, icc3_hi),
      TRUE         ~ sprintf("N=%d (insufficient)", n_pairs)
    )
  )

gt_icc <- icc_tbl |>
  select(
    Biomarker = label,
    `N patients` = n_patients,
    `N obs` = n_obs,
    `ICC (LMM, 95% CI)` = icc_fmt,
    Interpretation = interpretation,
    `N V1-V2 pairs` = n_pairs,
    `ICC3 V1-V2 (95% CI)` = icc3_fmt
  ) |>
  gt() |>
  tab_header(
    title    = "Intraclass Correlation Coefficients for Plasma AD Biomarkers",
    subtitle = "DAC CRF cohort, longitudinal subcohort (N=51)"
  ) |>
  tab_spanner(
    label   = "LMM-based ICC (all visits)",
    columns = c("N patients", "N obs", "ICC (LMM, 95% CI)", "Interpretation")
  ) |>
  tab_spanner(
    label   = "Classic ICC3 (V1-V2 pairs only)",
    columns = c("N V1-V2 pairs", "ICC3 V1-V2 (95% CI)")
  ) |>
  tab_footnote(
    footnote = "LMM ICC: random-intercept model (lme4 REML); ICC = patient variance / total variance. 95% CI via parametric bootstrap (500 replicates).",
    locations = cells_column_spanners("LMM-based ICC (all visits)")
  ) |>
  tab_footnote(
    footnote = "Classic ICC3: two-way mixed, consistency, single measures (psych::ICC). Restricted to patients with both V1 and V2 biomarker values.",
    locations = cells_column_spanners("Classic ICC3 (V1-V2 pairs only)")
  ) |>
  tab_style(
    style     = cell_fill(color = "#fff3cd"),
    locations = cells_body(
      columns = Interpretation,
      rows    = Interpretation %in% c("Poor", "Moderate")
    )
  ) |>
  tab_style(
    style     = cell_fill(color = "#d4edda"),
    locations = cells_body(
      columns = Interpretation,
      rows    = Interpretation %in% c("Good", "Excellent")
    )
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_icc, file.path(out_dir, "table_icc_summary.html"))
gtsave(gt_icc, file.path(out_dir, "table_icc_summary.png"), expand = 20)
message("Saved: table_icc_summary")

# -- Table 2: CV summary -------------------------------------------------------

gt_cv <- cv_summary |>
  mutate(
    cv_fmt       = sprintf("%.1f (%.1f\u2013%.1f)", median_cv, q25_cv, q75_cv),
    pct_above_fmt = sprintf("%.0f%%", pct_above_20)
  ) |>
  select(
    Biomarker      = label,
    `N patients`   = n_patients,
    `Median CV % (IQR)` = cv_fmt,
    `>20% CV`      = pct_above_fmt
  ) |>
  gt() |>
  tab_header(
    title    = "Within-Patient Coefficient of Variation",
    subtitle = "Serial plasma AD biomarkers, DAC CRF cohort"
  ) |>
  tab_footnote(
    footnote = "CV = SD/mean \u00d7 100% computed across all available visits (\u22652 visits per patient). >20% CV threshold from Mayfield et al. preanalytical total allowable error."
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_cv, file.path(out_dir, "table_cv_summary.html"))
gtsave(gt_cv, file.path(out_dir, "table_cv_summary.png"), expand = 20)
message("Saved: table_cv_summary")

# -- Figure 1: CV distributions (violin + box + jitter) -----------------------

# Order biomarkers by median CV
bm_order <- cv_summary |>
  arrange(median_cv) |>
  pull(label)

cv_plot_df <- cv_long |>
  mutate(label = factor(label, levels = bm_order))

fig_cv <- ggplot(cv_plot_df, aes(x = label, y = cv_pct)) +
  geom_violin(fill = "#4575b4", alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "#4575b4") +
  geom_jitter(width = 0.08, size = 1.8, alpha = 0.55, color = "#4575b4") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#E41A1C", linewidth = 0.7) +
  annotate("text", x = 0.55, y = 21.5, label = "20% TEa", color = "#E41A1C",
           hjust = 0, size = 3.5) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Within-Patient Coefficient of Variation by Biomarker",
    subtitle = "Longitudinal subcohort (N>=51, >=2 visits); dashed line = 20% preanalytical TEa",
    x        = NULL,
    y        = "Within-patient CV (%)"
  ) +
  theme_research() +
  theme(legend.position = "none")

fsz <- get_figure_size("violin_single")
ggsave(
  file.path(out_dir, "fig_cv_distributions.pdf"),
  fig_cv, width = fsz["width"], height = fsz["height"]
)
ggsave(
  file.path(out_dir, "fig_cv_distributions.png"),
  fig_cv, width = fsz["width"], height = fsz["height"], dpi = 300
)
message("Saved: fig_cv_distributions")

# -- Figure 2: ICC forest plot -------------------------------------------------

icc_plot_df <- icc_results |>
  mutate(
    label = factor(label, levels = rev(map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)))
  ) |>
  filter(!is.na(icc))

# Interpretation bands as background
interp_bands <- tibble(
  ymin  = c(0, 0.5, 0.75, 0.9),
  ymax  = c(0.5, 0.75, 0.9, 1.0),
  fill  = c("#f8d7da", "#fff3cd", "#d4edda", "#cce5ff"),
  label = c("Poor (<0.50)", "Moderate (0.50-0.75)", "Good (0.75-0.90)", "Excellent (>0.90)")
)

fig_icc <- ggplot() +
  # Background bands
  geom_rect(
    data = interp_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = label),
    alpha = 0.35, inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = setNames(interp_bands$fill, interp_bands$label),
    name   = "Reliability"
  ) +
  # Error bars + points
  geom_errorbar(
    data = icc_plot_df,
    aes(x = label, ymin = icc_lo, ymax = icc_hi),
    width = 0.15, linewidth = 0.7, color = "#1d3557"
  ) +
  geom_point(
    data = icc_plot_df,
    aes(x = label, y = icc),
    size = 3.5, color = "#1d3557"
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Intraclass Correlation Coefficient by Biomarker",
    subtitle = "LMM-based ICC (95% CI), longitudinal subcohort",
    x        = NULL,
    y        = "ICC"
  ) +
  theme_research(legend_position = "right") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6)))

fsz_forest <- get_figure_size("forest")
ggsave(
  file.path(out_dir, "fig_icc_forest.pdf"),
  fig_icc, width = fsz_forest["width"], height = fsz_forest["height"] * 0.5
)
ggsave(
  file.path(out_dir, "fig_icc_forest.png"),
  fig_icc, width = fsz_forest["width"], height = fsz_forest["height"] * 0.5, dpi = 300
)
message("Saved: fig_icc_forest")

# -- Figure 3: CV predictor dot plot (pTau217 and NfL) -------------------------

cor_all <- bind_rows(cor_ptau, cor_nfl) |>
  filter(!is.na(rho)) |>
  mutate(
    sig      = p < 0.05,
    var_label = recode(variable,
      age_at_enrol = "Age",
      egfr         = "eGFR",
      sex_f        = "Sex (female)",
      htn_f        = "Hypertension",
      dm2_f        = "Diabetes (T2)",
      ckd_f        = "CKD",
      apoe4        = "APOE4 carrier",
      amyloid_pos  = "Amyloid+ (PET)",
      ptau_high    = "pTau217 High risk"
    )
  )

if (nrow(cor_all) > 0) {
  fig_cor <- ggplot(cor_all, aes(x = rho, y = reorder(var_label, rho), color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(
      aes(x = 0, xend = rho, y = reorder(var_label, rho), yend = reorder(var_label, rho)),
      linewidth = 0.5
    ) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "#999999"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05"),
      name   = NULL
    ) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    facet_wrap(~ biomarker, ncol = 2) +
    labs(
      title    = "Predictors of Within-Patient CV",
      subtitle = "Spearman correlation with per-patient CV%; pTau217 and NfL",
      x        = "Spearman rho",
      y        = NULL
    ) +
    theme_research()

  fsz_forest2 <- get_figure_size("forest")
  ggsave(
    file.path(out_dir, "fig_cv_predictors.pdf"),
    fig_cor, width = fsz_forest2["width"], height = fsz_forest2["height"] * 0.65
  )
  ggsave(
    file.path(out_dir, "fig_cv_predictors.png"),
    fig_cor, width = fsz_forest2["width"], height = fsz_forest2["height"] * 0.65, dpi = 300
  )
  message("Saved: fig_cv_predictors")
} else {
  message("Skipping fig_cv_predictors: insufficient data for correlations")
}

# -- Write ICC + CV CSVs for downstream use ------------------------------------

write_csv(icc_results, file.path(out_dir, "icc_lmm_results.csv"))
write_csv(classic_results, file.path(out_dir, "icc_classic_v1v2.csv"))
write_csv(cv_long,    file.path(out_dir, "cv_per_patient.csv"))
write_csv(cv_summary, file.path(out_dir, "cv_summary.csv"))
write_csv(bind_rows(cor_ptau, cor_nfl), file.path(out_dir, "cv_predictor_correlations.csv"))

message("\nDone. Outputs written to: ", out_dir)
