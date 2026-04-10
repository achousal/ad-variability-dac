#!/usr/bin/env Rscript
# dac_mixed_model.R -- Random slopes mixed model for individual biomarker trajectories
#
# Fits linear mixed models with random intercepts + slopes to estimate each
# patient's rate of change, shrunk toward the population mean. Answers whether
# the distribution of individual slopes is centered at zero (oscillatory/noise)
# or has tails suggesting a subgroup with genuine drift.
#
# Usage:
#   Rscript analysis/dac_mixed_model.R [--data-dir PATH] [--output-dir PATH]
#
# Models fitted per biomarker:
#   M0: value ~ 1 + (1 | record_id)                         [intercept only]
#   M1: value ~ days_from_v1 + (1 | record_id)              [fixed slope, random intercept]
#   M2: value ~ days_from_v1 + (1 + days_from_v1 | record_id) [random slopes]
#
# M2 is the primary model. Falls back to M1 if M2 fails to converge or
# produces singular fit (often happens with N~40 and 2 visits/patient).
#
# Requires: tidyverse, lme4, lmerTest, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(patchwork)
  library(lme4)
  library(gt)
})

# lmerTest: adds p-values to lmer summaries; load after lme4
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# -- Paths ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) normalizePath(file.path(getwd(), "analysis"), mustWork = FALSE)
)
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

data_dir <- file.path(project_dir, "data")
out_dir  <- file.path(project_dir, "results", "dac_mixed_model")

for (i in seq_along(args)) {
  if (args[i] == "--data-dir"   && i < length(args)) data_dir <- args[i + 1]
  if (args[i] == "--output-dir" && i < length(args)) out_dir  <- args[i + 1]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "theme_research.R"))
source(file.path(script_dir, "palettes.R"))

message("Reading data from: ", data_dir)
message("Writing results to: ", out_dir)

# -- Load data -----------------------------------------------------------------

long <- read_csv(file.path(data_dir, "dac_long.csv"), show_col_types = FALSE) |>
  mutate(
    visit      = factor(visit, levels = c("V1", "V2", "V3", "V4")),
    blood_date = as.Date(blood_date)
  )

# V1 amyloid risk for stratification
v1_amyloid <- long |>
  filter(visit == "V1", !is_redraw) |>
  select(record_id, ptau217_risk, lucent_interpretation) |>
  distinct(record_id, .keep_all = TRUE)

# Longitudinal subcohort: >=2 visits, non-redraw, has days_from_v1
long_sub <- long |>
  filter(!is_redraw, !is.na(days_from_v1)) |>
  group_by(record_id) |>
  filter(n_distinct(visit) >= 2) |>
  ungroup()

message(sprintf("Longitudinal subcohort: %d patients, %d rows",
                n_distinct(long_sub$record_id), nrow(long_sub)))

# -- Biomarker metadata --------------------------------------------------------

biomarkers <- list(
  ptau217   = list(label = "pTau217",   unit = "pg/mL", log_transform = TRUE),
  nfl       = list(label = "NfL",       unit = "pg/mL", log_transform = TRUE),
  lucent_ad = list(label = "Lucent AD", unit = "score",  log_transform = FALSE)
)

# -- Helper: fit mixed models and extract slopes --------------------------------

fit_trajectory_models <- function(long_df, bm, v1_amyloid_df) {
  meta <- biomarkers[[bm]]

  dat <- long_df |>
    select(record_id, visit, days_from_v1, value = all_of(bm)) |>
    filter(!is.na(value)) |>
    mutate(
      # Scale time to years for interpretable slope units
      years_from_v1 = days_from_v1 / 365.25,
      # Log-transform skewed analytes for better linearity assumption
      value_fit = if (meta$log_transform) log(value) else value
    )

  n_pts <- n_distinct(dat$record_id)
  message(sprintf("  %s: %d patients, %d obs", meta$label, n_pts, nrow(dat)))

  if (n_pts < 5) {
    message("  Skipping: insufficient N")
    return(NULL)
  }

  # --- M0: intercept only (null) ---
  m0 <- lmer(value_fit ~ 1 + (1 | record_id), data = dat, REML = TRUE)

  # --- M1: fixed slope, random intercept ---
  m1 <- lmer(value_fit ~ years_from_v1 + (1 | record_id), data = dat, REML = TRUE)

  # --- M2: random slopes ---
  m2 <- tryCatch(
    lmer(value_fit ~ years_from_v1 + (1 + years_from_v1 | record_id),
         data = dat, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 1e5))),
    error = function(e) {
      message("  M2 failed: ", conditionMessage(e), " -- using M1")
      NULL
    }
  )

  singular_m2 <- !is.null(m2) && isSingular(m2, tol = 1e-4)
  if (singular_m2) {
    message("  M2 singular fit -- reporting but flagging")
  }

  primary_model <- if (!is.null(m2)) m2 else m1
  primary_label <- if (!is.null(m2)) "M2 (random slopes)" else "M1 (fixed slope)"
  message(sprintf("  Primary model: %s%s", primary_label,
                  if (singular_m2) " [SINGULAR]" else ""))

  # --- Variance components ---
  vc_m0 <- as.data.frame(VarCorr(m0))
  vc_primary <- as.data.frame(VarCorr(primary_model))

  # --- LRT: M1 vs M0 (does time matter at all?) ---
  lrt_m1_m0 <- tryCatch(
    anova(m0, m1, refit = TRUE),
    error = function(e) NULL
  )

  # --- LRT: M2 vs M1 (do random slopes improve fit?) ---
  lrt_m2_m1 <- if (!is.null(m2)) {
    tryCatch(anova(m1, m2, refit = TRUE), error = function(e) NULL)
  } else NULL

  # --- Extract fixed effects ---
  fe <- fixef(primary_model)
  fe_ci <- tryCatch(
    confint(primary_model, method = "Wald", parm = "years_from_v1"),
    error = function(e) matrix(c(NA, NA), nrow = 1)
  )

  slope_fixed       <- fe[["years_from_v1"]]
  slope_fixed_lo    <- if (!is.null(fe_ci)) fe_ci[1] else NA
  slope_fixed_hi    <- if (!is.null(fe_ci)) fe_ci[2] else NA

  # Back-transform slope for log-transformed analytes:
  # log-scale slope (per year) -> approximate % change per year
  slope_pct_yr <- if (meta$log_transform) exp(slope_fixed) - 1 else NA_real_

  # --- Extract patient-level slopes (BLUPs) ---
  re <- ranef(primary_model)$record_id

  if ("years_from_v1" %in% colnames(re)) {
    slopes_df <- re |>
      tibble::rownames_to_column("record_id") |>
      mutate(record_id = as.integer(record_id)) |>
      rename(slope_re = years_from_v1, intercept_re = `(Intercept)`) |>
      mutate(
        # Individual slope = population slope + random deviation
        slope_individual = slope_fixed + slope_re
      )
  } else {
    # M1: only random intercepts; individual slope = fixed slope for everyone
    slopes_df <- re |>
      tibble::rownames_to_column("record_id") |>
      mutate(
        record_id        = as.integer(record_id),
        intercept_re     = `(Intercept)`,
        slope_re         = 0,
        slope_individual = slope_fixed
      )
  }

  # Attach baseline value and amyloid risk
  v1_vals <- dat |>
    filter(visit == "V1") |>
    select(record_id, baseline_value = value, baseline_fit = value_fit)

  slopes_df <- slopes_df |>
    left_join(v1_vals, by = "record_id") |>
    left_join(v1_amyloid_df |> select(record_id, ptau217_risk), by = "record_id") |>
    mutate(
      biomarker        = bm,
      label            = meta$label,
      log_transformed  = meta$log_transform,
      singular_m2      = isTRUE(singular_m2)
    )

  # --- Fitted trajectories for plotting ---
  fitted_df <- dat |>
    left_join(slopes_df |> select(record_id, intercept_re, slope_re), by = "record_id") |>
    mutate(
      fitted_value_fit = (fixef(primary_model)[["(Intercept)"]] + intercept_re) +
        (slope_fixed + slope_re) * years_from_v1,
      fitted_value = if (meta$log_transform) exp(fitted_value_fit) else fitted_value_fit
    )

  list(
    biomarker       = bm,
    label           = meta$label,
    model_m0        = m0,
    model_m1        = m1,
    model_m2        = m2,
    primary_model   = primary_model,
    primary_label   = primary_label,
    singular        = singular_m2,
    lrt_m1_m0       = lrt_m1_m0,
    lrt_m2_m1       = lrt_m2_m1,
    slope_fixed     = slope_fixed,
    slope_fixed_lo  = slope_fixed_lo,
    slope_fixed_hi  = slope_fixed_hi,
    slope_pct_yr    = slope_pct_yr,
    slopes_df       = slopes_df,
    fitted_df       = fitted_df,
    n_patients      = n_pts,
    vc_primary      = vc_primary,
    log_transform   = meta$log_transform
  )
}

# -- Run models ----------------------------------------------------------------

message("\nFitting mixed models...")
model_results <- map(names(biomarkers), function(bm) {
  message(sprintf("\n--- %s ---", biomarkers[[bm]]$label))
  fit_trajectory_models(long_sub, bm, v1_amyloid)
}) |>
  set_names(names(biomarkers)) |>
  compact()

# -- Model summary table -------------------------------------------------------

message("\nBuilding model summary table...")

summary_df <- map_dfr(model_results, function(res) {
  # LRT p-values
  lrt_time_p <- tryCatch(
    anova(res$model_m0, res$model_m1, refit = TRUE)[2, "Pr(>Chisq)"],
    error = function(e) NA_real_
  )
  lrt_slopes_p <- if (!is.null(res$lrt_m2_m1)) {
    tryCatch(res$lrt_m2_m1[2, "Pr(>Chisq)"], error = function(e) NA_real_)
  } else NA_real_

  # Variance components from primary model
  vc <- res$vc_primary
  var_int  <- vc$vcov[vc$grp == "record_id" & vc$var1 == "(Intercept)" & is.na(vc$var2)]
  var_slp  <- vc$vcov[vc$grp == "record_id" & vc$var1 == "years_from_v1" & is.na(vc$var2)]
  var_res  <- vc$vcov[vc$grp == "Residual"]

  if (length(var_slp) == 0) var_slp <- NA_real_

  tibble(
    label          = res$label,
    n_patients     = res$n_patients,
    model          = res$primary_label,
    singular       = res$singular,
    slope_fixed    = res$slope_fixed,
    slope_lo       = res$slope_fixed_lo,
    slope_hi       = res$slope_fixed_hi,
    slope_pct_yr   = res$slope_pct_yr,
    lrt_time_p     = lrt_time_p,
    lrt_slopes_p   = lrt_slopes_p,
    var_intercept  = var_int,
    var_slope      = var_slp,
    var_residual   = var_res
  )
})

message("\nModel summary:")
print(summary_df |> select(label, model, slope_fixed, slope_pct_yr, lrt_time_p, singular))

# -- GT table ------------------------------------------------------------------

gt_models <- summary_df |>
  mutate(
    slope_fmt = case_when(
      !is.na(slope_lo) & !is.na(slope_hi) ~
        sprintf("%.4f (%.4f, %.4f)", slope_fixed, slope_lo, slope_hi),
      TRUE ~ sprintf("%.4f", slope_fixed)
    ),
    pct_yr_fmt  = if_else(!is.na(slope_pct_yr),
                          sprintf("%.1f%%/yr", slope_pct_yr * 100), "—"),
    time_p_fmt  = sprintf("%.3f", lrt_time_p),
    slope_p_fmt = if_else(!is.na(lrt_slopes_p), sprintf("%.3f", lrt_slopes_p), "—"),
    singular_fmt = if_else(singular, "Yes", "No")
  ) |>
  select(
    Biomarker      = label,
    `N patients`   = n_patients,
    Model          = model,
    Singular       = singular_fmt,
    `Fixed slope (95% CI, log-scale/yr)` = slope_fmt,
    `% change/yr`  = pct_yr_fmt,
    `LRT: time p`  = time_p_fmt,
    `LRT: random slopes p` = slope_p_fmt
  ) |>
  gt() |>
  tab_header(
    title    = "Linear Mixed Model: Individual Biomarker Trajectories",
    subtitle = "DAC CRF longitudinal subcohort"
  ) |>
  tab_footnote(
    footnote = "Fixed slope = population-average rate of change (log-scale per year for pTau217/NfL). % change/yr = exp(slope) - 1, back-transformed. LRT time p: M1 vs M0. LRT random slopes p: M2 vs M1. Singular: random slope variance near zero."
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_models, file.path(out_dir, "table_model_summary.html"))
gtsave(gt_models, file.path(out_dir, "table_model_summary.png"), expand = 20)
message("Saved: table_model_summary")

# -- Figure 1: Distribution of individual slopes -------------------------------

slopes_all <- map_dfr(model_results, ~ .x$slopes_df) |>
  mutate(
    label        = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)),
    ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High"))
  )

risk_colors <- c("High" = "#E41A1C", "Intermediate" = "#FF7F00",
                 "Low" = "#4DAF4A", "NA" = "#999999")

fig_slope_dist <- ggplot(slopes_all,
                         aes(x = slope_individual, fill = ptau217_risk)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_histogram(bins = 15, alpha = 0.75, color = "white", position = "stack") +
  scale_fill_manual(values = risk_colors, na.value = "#cccccc",
                    name = "pTau217 risk") +
  facet_wrap(~ label, scales = "free_x", ncol = 3) +
  labs(
    title    = "Distribution of Individual Slopes (BLUPs)",
    subtitle = "Random slopes model; dashed = zero (no change); color = amyloid risk",
    x        = "Individual slope (log-scale units/year)",
    y        = "N patients"
  ) +
  theme_research()

fsz <- get_figure_size("violin_grouped")
ggsave(file.path(out_dir, "fig_slope_distribution.pdf"),
       fig_slope_dist, width = fsz["width"] * 1.4, height = fsz["height"])
ggsave(file.path(out_dir, "fig_slope_distribution.png"),
       fig_slope_dist, width = fsz["width"] * 1.4, height = fsz["height"], dpi = 300)
message("Saved: fig_slope_distribution")

# -- Figure 2: Slope vs baseline value ----------------------------------------

fig_slope_base <- ggplot(
  slopes_all |> filter(!is.na(baseline_value)),
  aes(x = baseline_value, y = slope_individual, color = ptau217_risk)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8,
              linetype = "solid", alpha = 0.15) +
  scale_color_manual(values = risk_colors, na.value = "#cccccc",
                     name = "pTau217 risk") +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  labs(
    title    = "Individual Slope vs. Baseline Value",
    subtitle = "Regression to the mean expected (higher baseline -> lower slope)",
    x        = "Baseline value (V1)",
    y        = "Individual slope (/year)"
  ) +
  theme_research()

fsz_s <- get_figure_size("scatter_multi")
ggsave(file.path(out_dir, "fig_slope_vs_baseline.pdf"),
       fig_slope_base, width = fsz_s["width"] * 0.8, height = fsz_s["height"] * 0.55)
ggsave(file.path(out_dir, "fig_slope_vs_baseline.png"),
       fig_slope_base, width = fsz_s["width"] * 0.8, height = fsz_s["height"] * 0.55,
       dpi = 300)
message("Saved: fig_slope_vs_baseline")

# -- Figure 3: Fitted trajectories (spaghetti + population line) ---------------

make_fitted_spaghetti <- function(bm) {
  res  <- model_results[[bm]]
  if (is.null(res)) return(NULL)
  meta <- biomarkers[[bm]]

  fd <- res$fitted_df |>
    left_join(v1_amyloid |> select(record_id, ptau217_risk), by = "record_id") |>
    mutate(ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High")))

  # Population mean line from fixed effects
  t_range <- range(fd$years_from_v1)
  pop_line <- tibble(
    years_from_v1 = seq(t_range[1], t_range[2], length.out = 50),
    fitted_value  = {
      fe <- fixef(res$primary_model)
      yhat_log <- fe[["(Intercept)"]] + fe[["years_from_v1"]] * seq(t_range[1], t_range[2], length.out = 50)
      if (meta$log_transform) exp(yhat_log) else yhat_log
    }
  )

  ggplot() +
    geom_line(
      data = fd,
      aes(x = years_from_v1, y = fitted_value, group = record_id,
          color = ptau217_risk),
      alpha = 0.35, linewidth = 0.45
    ) +
    geom_line(
      data = pop_line,
      aes(x = years_from_v1, y = fitted_value),
      color = "black", linewidth = 1.2, linetype = "solid"
    ) +
    scale_color_manual(values = risk_colors, na.value = "#cccccc",
                       name = "pTau217 risk") +
    labs(
      title    = meta$label,
      subtitle = sprintf("N=%d; black = population mean trajectory", res$n_patients),
      x        = "Years from V1",
      y        = sprintf("%s (%s)", meta$label, meta$unit)
    ) +
    theme_research()
}

fitted_plots <- compact(map(names(biomarkers), make_fitted_spaghetti))

if (length(fitted_plots) > 0) {
  fig_fitted <- wrap_plots(fitted_plots, ncol = length(fitted_plots)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  fsz_f <- get_figure_size("scatter_multi")
  ggsave(file.path(out_dir, "fig_fitted_trajectories.pdf"),
         fig_fitted, width = fsz_f["width"], height = fsz_f["height"] * 0.55)
  ggsave(file.path(out_dir, "fig_fitted_trajectories.png"),
         fig_fitted, width = fsz_f["width"], height = fsz_f["height"] * 0.55, dpi = 300)
  message("Saved: fig_fitted_trajectories")
}

# -- Figure 4: Slope by pTau217 risk (box + jitter) ----------------------------

fig_slope_risk <- slopes_all |>
  filter(!is.na(ptau217_risk)) |>
  ggplot(aes(x = ptau217_risk, y = slope_individual,
             fill = ptau217_risk, color = ptau217_risk)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "white",
               linewidth = 0.6, color = "grey30") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.75) +
  scale_color_manual(values = risk_colors, guide = "none") +
  scale_fill_manual(values  = risk_colors, guide = "none") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    title    = "Individual Slopes by Amyloid Risk",
    subtitle = "BLUPs from random slopes model; dashed = zero",
    x        = "pTau217 risk",
    y        = "Individual slope (/year)"
  ) +
  theme_research()

ggsave(file.path(out_dir, "fig_slope_by_risk.pdf"),
       fig_slope_risk, width = fsz["width"] * 1.4, height = fsz["height"])
ggsave(file.path(out_dir, "fig_slope_by_risk.png"),
       fig_slope_risk, width = fsz["width"] * 1.4, height = fsz["height"], dpi = 300)
message("Saved: fig_slope_by_risk")

# -- Write CSVs ----------------------------------------------------------------

write_csv(slopes_all,  file.path(out_dir, "slopes_per_patient.csv"))
write_csv(summary_df,  file.path(out_dir, "model_summary.csv"))

message("\nDone. Outputs written to: ", out_dir)
