#!/usr/bin/env Rscript
# dac_eda.R -- Exploratory longitudinal descriptive analysis of DAC CRF data
#
# Reads dac_long.csv and dac_baseline.csv produced by dac_pull.py.
# Outputs tables and figures to _code/results/dac_eda/.
#
# Usage:
#   Rscript dac_eda.R [--data-dir PATH] [--output-dir PATH]
#
# Requires: tidyverse, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(gt)
})

# -- Paths --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    # fallback: assume running from _code/R/
    normalizePath(file.path(getwd(), "_code", "R"), mustWork = FALSE)
  }
)
code_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

data_dir <- file.path(code_dir, "data")
out_dir  <- file.path(code_dir, "results", "dac_eda")

# Parse simple --data-dir / --output-dir args
for (i in seq_along(args)) {
  if (args[i] == "--data-dir" && i < length(args))   data_dir <- args[i + 1]
  if (args[i] == "--output-dir" && i < length(args))  out_dir <- args[i + 1]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Source theme + palettes
source(file.path(script_dir, "theme_research.R"))
source(file.path(script_dir, "palettes.R"))

message("Reading data from: ", data_dir)
message("Writing results to: ", out_dir)

# -- Load data ----------------------------------------------------------------

long <- read_csv(file.path(data_dir, "dac_long.csv"), show_col_types = FALSE) |>
  mutate(
    visit = factor(visit, levels = c("V1", "V2", "V3", "V4")),
    visit_date = as.Date(visit_date),
    blood_date = as.Date(blood_date)
  )

baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("Female", "Male")))

message(sprintf("Long: %d rows, %d patients", nrow(long), n_distinct(long$record_id)))
message(sprintf("Baseline: %d patients", nrow(baseline)))

# Identify longitudinal subcohort (>=2 visits)
visit_counts <- long |>
  group_by(record_id) |>
  summarise(n_visits = n_distinct(visit), .groups = "drop")

longitudinal_ids <- visit_counts |>
  filter(n_visits >= 2) |>
  pull(record_id)

baseline <- baseline |>
  mutate(longitudinal = record_id %in% longitudinal_ids)

message(sprintf("Longitudinal subcohort: %d patients", length(longitudinal_ids)))

# =============================================================================
# A. TABLE 1 -- Cohort characterization
# =============================================================================
message("A. Table 1...")

summarise_group <- function(df, group_label) {
  n <- nrow(df)
  tibble(
    group = group_label,
    N = n,
    # Age
    age_mean = sprintf("%.1f (%.1f)", mean(df$age_at_enrol, na.rm = TRUE),
                       sd(df$age_at_enrol, na.rm = TRUE)),
    age_median = sprintf("%.0f [%.0f, %.0f]",
                         median(df$age_at_enrol, na.rm = TRUE),
                         quantile(df$age_at_enrol, 0.25, na.rm = TRUE),
                         quantile(df$age_at_enrol, 0.75, na.rm = TRUE)),
    age_n = sum(!is.na(df$age_at_enrol)),
    # Sex
    female_n = sum(df$sex == "Female", na.rm = TRUE),
    female_pct = sprintf("%d (%.1f%%)", sum(df$sex == "Female", na.rm = TRUE),
                         100 * mean(df$sex == "Female", na.rm = TRUE)),
    male_n = sum(df$sex == "Male", na.rm = TRUE),
    # Race (top categories)
    race_white = sum(str_detect(df$race %||% "", "White"), na.rm = TRUE),
    race_black = sum(str_detect(df$race %||% "", "Black"), na.rm = TRUE),
    race_asian = sum(str_detect(df$race %||% "", "Asian"), na.rm = TRUE),
    race_other = sum(str_detect(df$race %||% "", "Other"), na.rm = TRUE),
    race_missing = sum(is.na(df$race)),
    # Ethnicity
    hispanic = sum(str_detect(df$ethnicity %||% "", "Hispanic") &
                     !str_detect(df$ethnicity %||% "", "Non"), na.rm = TRUE),
    non_hispanic = sum(str_detect(df$ethnicity %||% "", "Non-Hispanic"), na.rm = TRUE),
    # Education
    edu_hs_or_less = sum(df$education %in% c("No high school", "Some high school",
                                              "High school diploma"), na.rm = TRUE),
    edu_college = sum(df$education %in% c("2-year degree", "4-year degree"), na.rm = TRUE),
    edu_graduate = sum(df$education == "Graduate+", na.rm = TRUE),
    edu_missing = sum(is.na(df$education)),
    # Comorbidities
    htn = sum(df$htn == "Yes", na.rm = TRUE),
    hld = sum(df$hld == "Yes", na.rm = TRUE),
    dm2 = sum(df$dm2 == "Yes", na.rm = TRUE),
    heartdisease = sum(df$heartdisease == "Yes", na.rm = TRUE),
    stroke = sum(df$stroke == "Yes", na.rm = TRUE),
    # Cognitive
    moca_mean = sprintf("%.1f (%.1f)", mean(df$moca_raw, na.rm = TRUE),
                        sd(df$moca_raw, na.rm = TRUE)),
    moca_n = sum(!is.na(df$moca_raw)),
    mmse_mean = sprintf("%.1f (%.1f)", mean(df$mmse_raw, na.rm = TRUE),
                        sd(df$mmse_raw, na.rm = TRUE)),
    mmse_n = sum(!is.na(df$mmse_raw)),
    # Diagnosis (where assigned)
    dx_mci = sum(str_detect(df$diagnosis %||% "", "MCI"), na.rm = TRUE),
    dx_ad = sum(str_detect(df$diagnosis %||% "", "AD"), na.rm = TRUE),
    dx_assigned = sum(!is.na(df$diagnosis)),
    # APOE (where available)
    apoe_e4_carrier = sum(df$apoe_genotype %in% c("E2/E4", "E3/E4", "E4/E4"),
                          na.rm = TRUE),
    apoe_n = sum(!is.na(df$apoe_genotype) & df$apoe_genotype != "Unknown"),
    # Study status
    dropout = sum(df$dropout == TRUE, na.rm = TRUE),
    deceased = sum(df$deceased == TRUE, na.rm = TRUE),
    ltfu = sum(df$ltfu == TRUE, na.rm = TRUE),
  )
}

t1_all <- summarise_group(baseline, "All V1")
t1_long <- summarise_group(baseline |> filter(longitudinal), "Longitudinal (>=2 visits)")

table1 <- bind_rows(t1_all, t1_long)

write_csv(table1, file.path(out_dir, "table1.csv"))

# Formatted Table 1 via gt
fmt_row <- function(label, all_val, long_val) {
  tibble(Characteristic = label, `All V1` = all_val, `Longitudinal` = long_val)
}

t1_fmt <- bind_rows(
  fmt_row("N", as.character(t1_all$N), as.character(t1_long$N)),
  fmt_row("Age, mean (SD)", t1_all$age_mean, t1_long$age_mean),
  fmt_row("  N with age", as.character(t1_all$age_n), as.character(t1_long$age_n)),
  fmt_row("Female, n (%)", t1_all$female_pct, t1_long$female_pct),
  fmt_row("Race: White", as.character(t1_all$race_white), as.character(t1_long$race_white)),
  fmt_row("Race: Black/AA", as.character(t1_all$race_black), as.character(t1_long$race_black)),
  fmt_row("Race: Asian/PI", as.character(t1_all$race_asian), as.character(t1_long$race_asian)),
  fmt_row("Race: Other", as.character(t1_all$race_other), as.character(t1_long$race_other)),
  fmt_row("Hispanic", as.character(t1_all$hispanic), as.character(t1_long$hispanic)),
  fmt_row("Education: HS or less", as.character(t1_all$edu_hs_or_less), as.character(t1_long$edu_hs_or_less)),
  fmt_row("Education: College", as.character(t1_all$edu_college), as.character(t1_long$edu_college)),
  fmt_row("Education: Graduate+", as.character(t1_all$edu_graduate), as.character(t1_long$edu_graduate)),
  fmt_row("HTN", as.character(t1_all$htn), as.character(t1_long$htn)),
  fmt_row("HLD", as.character(t1_all$hld), as.character(t1_long$hld)),
  fmt_row("DM2", as.character(t1_all$dm2), as.character(t1_long$dm2)),
  fmt_row("Heart disease", as.character(t1_all$heartdisease), as.character(t1_long$heartdisease)),
  fmt_row("Stroke", as.character(t1_all$stroke), as.character(t1_long$stroke)),
  fmt_row("MoCA, mean (SD)", t1_all$moca_mean, t1_long$moca_mean),
  fmt_row("  N with MoCA", as.character(t1_all$moca_n), as.character(t1_long$moca_n)),
  fmt_row("MMSE, mean (SD)", t1_all$mmse_mean, t1_long$mmse_mean),
  fmt_row("  N with MMSE", as.character(t1_all$mmse_n), as.character(t1_long$mmse_n)),
  fmt_row("Dx: MCI", as.character(t1_all$dx_mci), as.character(t1_long$dx_mci)),
  fmt_row("Dx: AD", as.character(t1_all$dx_ad), as.character(t1_long$dx_ad)),
  fmt_row("Dx assigned (any)", as.character(t1_all$dx_assigned), as.character(t1_long$dx_assigned)),
  fmt_row("APOE E4 carrier", as.character(t1_all$apoe_e4_carrier), as.character(t1_long$apoe_e4_carrier)),
  fmt_row("  N with APOE", as.character(t1_all$apoe_n), as.character(t1_long$apoe_n)),
  fmt_row("Dropout", as.character(t1_all$dropout), as.character(t1_long$dropout)),
  fmt_row("Deceased", as.character(t1_all$deceased), as.character(t1_long$deceased)),
  fmt_row("LTFU", as.character(t1_all$ltfu), as.character(t1_long$ltfu)),
)

gt_table <- t1_fmt |>
  gt() |>
  tab_header(title = "Table 1. DAC CRF Cohort Characteristics") |>
  tab_source_note(sprintf("Data pulled %s. Longitudinal = patients with >=2 visits.",
                          Sys.Date()))
gtsave(gt_table, file.path(out_dir, "table1.html"))
message("  -> table1.csv, table1.html")


# =============================================================================
# B. ENROLLMENT FUNNEL
# =============================================================================
message("B. Enrollment funnel...")

funnel_data <- long |>
  group_by(record_id) |>
  summarise(
    has_v1 = any(visit == "V1"),
    has_v2 = any(visit == "V2"),
    has_v3 = any(visit == "V3"),
    has_v4 = any(visit == "V4"),
    .groups = "drop"
  )

funnel_counts <- tibble(
  stage = factor(c("V1", "V2", "V3", "V4"), levels = c("V1", "V2", "V3", "V4")),
  n = c(
    sum(funnel_data$has_v1),
    sum(funnel_data$has_v2),
    sum(funnel_data$has_v3),
    sum(funnel_data$has_v4)
  )
)

p_funnel <- ggplot(funnel_counts, aes(x = stage, y = n)) +
  geom_col(fill = "#377EB8", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, fontface = "bold", size = 5) +
  labs(
    title = "Enrollment Funnel: Patients per Visit Stage",
    x = "Visit", y = "N patients"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_research()

ggsave(file.path(out_dir, "enrollment_funnel.png"), p_funnel,
       width = 8, height = 6, dpi = 300)
message("  -> enrollment_funnel.png")


# =============================================================================
# C. TEMPORAL ALIGNMENT
# =============================================================================
message("C. Temporal alignment...")

# Inter-visit intervals
intervals <- long |>
  filter(!is.na(visit_interval), visit != "V1") |>
  mutate(transition = case_when(
    visit == "V2" ~ "V1 -> V2",
    visit == "V3" ~ "V2 -> V3",
    visit == "V4" ~ "V3 -> V4",
  )) |>
  filter(!is.na(transition))

p_timing <- ggplot(intervals, aes(x = transition, y = visit_interval)) +
  geom_violin(fill = "#377EB8", alpha = 0.3, color = "#377EB8") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, color = "#377EB8") +
  geom_hline(yintercept = c(14, 30, 90), linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  annotate("text", x = 0.5, y = c(14, 30, 90), label = c("14d", "30d", "90d"),
           hjust = 0, vjust = -0.5, size = 3, color = "grey40") +
  labs(
    title = "Inter-Visit Intervals",
    subtitle = "Dashed lines: 14, 30, 90 day reference windows",
    x = "Transition", y = "Days between visits"
  ) +
  theme_research()

ggsave(file.path(out_dir, "visit_timing.png"), p_timing,
       width = 8, height = 6, dpi = 300)

# Temporal alignment windows table
windows <- intervals |>
  group_by(transition) |>
  summarise(
    N = n(),
    median_days = median(visit_interval, na.rm = TRUE),
    iqr_lo = quantile(visit_interval, 0.25, na.rm = TRUE),
    iqr_hi = quantile(visit_interval, 0.75, na.rm = TRUE),
    min_days = min(visit_interval, na.rm = TRUE),
    max_days = max(visit_interval, na.rm = TRUE),
    within_14d = sum(visit_interval <= 14),
    within_30d = sum(visit_interval <= 30),
    within_90d = sum(visit_interval <= 90),
    outlier_6mo = sum(visit_interval > 180),
    .groups = "drop"
  )

write_csv(windows, file.path(out_dir, "temporal_windows.csv"))
message("  -> visit_timing.png, temporal_windows.csv")


# =============================================================================
# D. BIOMARKER TRAJECTORIES (spaghetti plots, days from V1)
# =============================================================================
message("D. Biomarker trajectories...")

biomarkers <- c("ptau217", "nfl", "gfap", "ab42_40", "lucent_ad", "ptau217_ab42_ratio")
biomarker_labels <- c(
  ptau217 = "pTau217", nfl = "NfL (pg/mL)", gfap = "GFAP (pg/mL)",
  ab42_40 = "Ab42/40 Ratio", lucent_ad = "Lucent AD Score",
  ptau217_ab42_ratio = "pTau217/Ab42-40"
)

long_bio <- long |>
  filter(record_id %in% longitudinal_ids) |>
  select(record_id, visit, days_from_v1, all_of(biomarkers)) |>
  pivot_longer(cols = all_of(biomarkers), names_to = "biomarker", values_to = "value") |>
  filter(!is.na(value), !is.na(days_from_v1)) |>
  mutate(biomarker = factor(biomarker, levels = biomarkers,
                            labels = biomarker_labels[biomarkers]))

# Group means
bio_means <- long_bio |>
  group_by(biomarker, visit, days_from_v1) |>
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") |>
  group_by(biomarker, visit) |>
  summarise(
    days_from_v1 = median(days_from_v1),
    mean_val = mean(mean_val),
    .groups = "drop"
  )

# Use visit-level means for cleaner display
bio_summary <- long_bio |>
  group_by(biomarker, visit) |>
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    se_val = sd(value, na.rm = TRUE) / sqrt(n()),
    median_days = median(days_from_v1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_traj <- ggplot(long_bio, aes(x = days_from_v1, y = value)) +
  geom_line(aes(group = record_id), alpha = 0.25, color = "grey60") +
  geom_point(aes(group = record_id), alpha = 0.4, size = 1.5, color = "grey60") +
  geom_pointrange(
    data = bio_summary,
    aes(x = median_days, y = mean_val,
        ymin = mean_val - se_val, ymax = mean_val + se_val),
    color = "#E41A1C", size = 0.8, linewidth = 0.8
  ) +
  facet_wrap(~biomarker, scales = "free_y", ncol = 3) +
  labs(
    title = "Biomarker Trajectories (Longitudinal Subcohort)",
    subtitle = "Grey: individual patients. Red: group mean +/- SE",
    x = "Days from V1", y = "Value"
  ) +
  theme_research()

ggsave(file.path(out_dir, "biomarker_trajectories.png"), p_traj,
       width = 14, height = 8, dpi = 300)
message("  -> biomarker_trajectories.png")


# =============================================================================
# E. BIOMARKER DISTRIBUTIONS per visit (all patients at each visit)
# =============================================================================
message("E. Biomarker distributions per visit...")

long_bio_all <- long |>
  select(record_id, visit, all_of(biomarkers)) |>
  pivot_longer(cols = all_of(biomarkers), names_to = "biomarker", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(biomarker = factor(biomarker, levels = biomarkers,
                            labels = biomarker_labels[biomarkers]))

p_dist <- ggplot(long_bio_all, aes(x = visit, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "#377EB8", alpha = 0.2, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1, color = "#377EB8") +
  facet_wrap(~biomarker, scales = "free_y", ncol = 3) +
  labs(
    title = "Biomarker Distributions by Visit",
    x = "Visit", y = "Value"
  ) +
  theme_research()

ggsave(file.path(out_dir, "biomarker_distributions.png"), p_dist,
       width = 14, height = 8, dpi = 300)
message("  -> biomarker_distributions.png")


# =============================================================================
# F. CLINICAL SCORE DISTRIBUTIONS (baseline only -- single form)
# =============================================================================
message("F. Clinical score distributions...")

clin_long <- baseline |>
  select(record_id, sex, moca_raw, mmse_raw, linus_total, badl) |>
  pivot_longer(cols = c(moca_raw, mmse_raw, linus_total, badl),
               names_to = "score", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(score = factor(score,
    levels = c("moca_raw", "mmse_raw", "linus_total", "badl"),
    labels = c("MoCA", "MMSE", "Linus Health Total", "bADL (/6)")
  ))

if (nrow(clin_long) > 0) {
  p_clin <- ggplot(clin_long, aes(x = score, y = value)) +
    geom_boxplot(outlier.shape = NA, fill = "#4DAF4A", alpha = 0.2, width = 0.5) +
    geom_jitter(aes(color = sex), width = 0.15, alpha = 0.5, size = 1.5) +
    facet_wrap(~score, scales = "free", ncol = 2) +
    labs(
      title = "Baseline Clinical Scores",
      subtitle = "Colored by sex",
      x = NULL, y = "Score"
    ) +
    scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"),
                       na.value = "grey60") +
    theme_research()

  ggsave(file.path(out_dir, "clinical_scores.png"), p_clin,
         width = 10, height = 8, dpi = 300)
  message("  -> clinical_scores.png")
} else {
  message("  -> skipped (no clinical score data)")
}


# =============================================================================
# G. DELTA ANALYSIS
# =============================================================================
message("G. Delta analysis...")

# Compute deltas relative to V1 for each biomarker
v1_vals <- long |>
  filter(visit == "V1") |>
  select(record_id, all_of(paste0(biomarkers))) |>
  rename_with(~paste0(.x, "_v1"), all_of(biomarkers))

deltas <- long |>
  filter(visit != "V1", record_id %in% longitudinal_ids) |>
  select(record_id, visit, days_from_v1, all_of(biomarkers)) |>
  left_join(v1_vals, by = "record_id")

# Compute deltas
for (bm in biomarkers) {
  deltas[[paste0("delta_", bm)]] <- deltas[[bm]] - deltas[[paste0(bm, "_v1")]]
}

delta_long <- deltas |>
  select(record_id, visit, days_from_v1,
         starts_with("delta_"), ends_with("_v1")) |>
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "biomarker",
    values_to = "delta",
    names_prefix = "delta_"
  ) |>
  filter(!is.na(delta)) |>
  mutate(biomarker = factor(biomarker, levels = biomarkers,
                            labels = biomarker_labels[biomarkers]))

# Add baseline value for regression-to-mean check
delta_with_base <- deltas |>
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "biomarker",
    values_to = "delta",
    names_prefix = "delta_"
  ) |>
  filter(!is.na(delta))

# Pivot baseline values
base_vals <- deltas |>
  select(record_id, visit, ends_with("_v1")) |>
  pivot_longer(
    cols = ends_with("_v1"),
    names_to = "biomarker",
    values_to = "baseline_val",
    names_pattern = "(.+)_v1"
  )

delta_merged <- delta_with_base |>
  select(record_id, visit, days_from_v1, biomarker, delta) |>
  left_join(base_vals |> select(record_id, biomarker, baseline_val),
            by = c("record_id", "biomarker"),
            relationship = "many-to-many") |>
  distinct() |>
  filter(!is.na(baseline_val)) |>
  mutate(biomarker = factor(biomarker, levels = biomarkers,
                            labels = biomarker_labels[biomarkers]))

  # -- TEa reference: 20% total allowable error (Mayfield et al. stability criterion) --
  # Compute percent change and flag TEa exceedance
  delta_merged <- delta_merged |>
    mutate(
      pct_change = ifelse(baseline_val != 0, 100 * delta / baseline_val, NA_real_),
      exceeds_tea = abs(pct_change) > 20
    )

  # TEa bounds for delta-vs-baseline plot: +/-20% of baseline_val
  # These appear as diagonal lines through origin on the delta vs baseline scatter
  tea_pct <- 0.20

if (nrow(delta_long) > 0) {
  # Delta distributions (with percent change overlay)
  p_delta_hist <- ggplot(delta_long, aes(x = delta)) +
    geom_histogram(fill = "#377EB8", alpha = 0.6, bins = 20) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~biomarker, scales = "free", ncol = 3) +
    labs(
      title = "Delta from V1: Distribution of Changes",
      x = "Change from V1", y = "Count"
    ) +
    theme_research()

  # Delta vs baseline with +/-20% TEa bands
  p_delta_base <- ggplot(delta_merged, aes(x = baseline_val, y = delta)) +
    geom_abline(slope = tea_pct, intercept = 0,
                linetype = "dotted", color = "#E41A1C", linewidth = 0.6) +
    geom_abline(slope = -tea_pct, intercept = 0,
                linetype = "dotted", color = "#E41A1C", linewidth = 0.6) +
    geom_point(aes(shape = exceeds_tea), alpha = 0.5, color = "#377EB8", size = 2) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4),
                       labels = c("Within TEa", "Exceeds 20% TEa"),
                       name = NULL) +
    geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~biomarker, scales = "free", ncol = 3) +
    labs(
      title = "Delta vs Baseline (Regression to Mean + 20% TEa Bands)",
      subtitle = "Red dotted lines: +/-20% total allowable error. X = exceeds TEa.",
      x = "V1 Baseline Value", y = "Delta from V1"
    ) +
    theme_research()

  # Percent change distribution with TEa threshold
  p_pct_change <- ggplot(delta_merged |> filter(!is.na(pct_change)),
                         aes(x = pct_change)) +
    geom_histogram(fill = "#377EB8", alpha = 0.6, bins = 25) +
    geom_vline(xintercept = c(-20, 20), linetype = "dotted",
               color = "#E41A1C", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~biomarker, scales = "free_y", ncol = 3) +
    labs(
      title = "Percent Change from V1",
      subtitle = "Red dotted lines: +/-20% TEa (Mayfield et al. stability threshold)",
      x = "% Change from V1", y = "Count"
    ) +
    theme_research()

  # Delta vs time interval
  p_delta_time <- ggplot(delta_long |> filter(!is.na(days_from_v1)),
                         aes(x = days_from_v1, y = delta)) +
    geom_point(alpha = 0.5, color = "#377EB8") +
    geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~biomarker, scales = "free", ncol = 3) +
    labs(
      title = "Delta vs Time from V1 (Time-Dependency Check)",
      x = "Days from V1", y = "Delta from V1"
    ) +
    theme_research()

  p_delta_combined <- p_delta_hist / p_pct_change / p_delta_base / p_delta_time +
    plot_annotation(
      title = "Delta Analysis: Biomarker Changes from V1",
      subtitle = "TEa = 20% total allowable error (Mayfield et al. pTau217/Ab1-42 stability)",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )

  ggsave(file.path(out_dir, "delta_analysis.png"), p_delta_combined,
         width = 14, height = 22, dpi = 300)

  # Delta summary table with TEa exceedance counts
  delta_summary <- delta_merged |>
    group_by(biomarker, visit) |>
    summarise(
      n = n(),
      mean_delta = mean(delta, na.rm = TRUE),
      sd_delta = sd(delta, na.rm = TRUE),
      median_delta = median(delta, na.rm = TRUE),
      iqr_lo = quantile(delta, 0.25, na.rm = TRUE),
      iqr_hi = quantile(delta, 0.75, na.rm = TRUE),
      mean_pct_change = mean(pct_change, na.rm = TRUE),
      median_pct_change = median(pct_change, na.rm = TRUE),
      n_exceeds_tea_20pct = sum(exceeds_tea, na.rm = TRUE),
      pct_exceeds_tea = 100 * mean(exceeds_tea, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(delta_summary, file.path(out_dir, "delta_summary.csv"))
  message("  -> delta_analysis.png, delta_summary.csv")
} else {
  message("  -> skipped (no delta data available)")
}


# =============================================================================
# H. VITALS TRAJECTORIES
# =============================================================================
message("H. Vitals trajectories...")

vitals <- c("bmi", "systolic", "diastolic")
vitals_labels <- c(bmi = "BMI", systolic = "Systolic BP (mmHg)",
                   diastolic = "Diastolic BP (mmHg)")

long_vitals <- long |>
  filter(record_id %in% longitudinal_ids) |>
  select(record_id, visit, days_from_v1, all_of(vitals)) |>
  pivot_longer(cols = all_of(vitals), names_to = "vital", values_to = "value") |>
  filter(!is.na(value), !is.na(days_from_v1)) |>
  mutate(vital = factor(vital, levels = vitals, labels = vitals_labels[vitals]))

if (nrow(long_vitals) > 0) {
  vital_summary <- long_vitals |>
    group_by(vital, visit) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      se_val = sd(value, na.rm = TRUE) / sqrt(n()),
      median_days = median(days_from_v1, na.rm = TRUE),
      .groups = "drop"
    )

  p_vitals <- ggplot(long_vitals, aes(x = days_from_v1, y = value)) +
    geom_line(aes(group = record_id), alpha = 0.25, color = "grey60") +
    geom_point(aes(group = record_id), alpha = 0.4, size = 1.5, color = "grey60") +
    geom_pointrange(
      data = vital_summary,
      aes(x = median_days, y = mean_val,
          ymin = mean_val - se_val, ymax = mean_val + se_val),
      color = "#E41A1C", size = 0.8, linewidth = 0.8
    ) +
    facet_wrap(~vital, scales = "free_y", ncol = 3) +
    labs(
      title = "Vitals Trajectories (Longitudinal Subcohort)",
      subtitle = "Grey: individual patients. Red: group mean +/- SE",
      x = "Days from V1", y = "Value"
    ) +
    theme_research()

  ggsave(file.path(out_dir, "vitals_trajectories.png"), p_vitals,
         width = 14, height = 5, dpi = 300)
  message("  -> vitals_trajectories.png")
} else {
  message("  -> skipped (no longitudinal vitals)")
}


# =============================================================================
# I. BASELINE CORRELATION HEATMAP
# =============================================================================
message("I. Baseline correlation heatmap...")

# Get V1 biomarker values
v1_data <- long |>
  filter(visit == "V1") |>
  select(record_id, ptau217, nfl, gfap, ab42_40, lucent_ad, ptau217_ab42_ratio,
         bmi, systolic, diastolic)

corr_data <- baseline |>
  select(record_id, age_at_enrol, moca_raw, mmse_raw) |>
  left_join(v1_data, by = "record_id") |>
  select(-record_id)

# Rename for display
names(corr_data) <- c("Age", "MoCA", "MMSE", "pTau217", "NfL", "GFAP",
                       "Ab42/40", "Lucent AD", "pTau217/Ab42-40",
                       "BMI", "Systolic", "Diastolic")

# Pairwise complete correlations
corr_mat <- cor(corr_data, use = "pairwise.complete.obs")

# Convert to long for ggplot
corr_long <- as.data.frame(as.table(corr_mat)) |>
  rename(var1 = Var1, var2 = Var2, r = Freq) |>
  mutate(
    var1 = factor(var1, levels = colnames(corr_mat)),
    var2 = factor(var2, levels = rev(colnames(corr_mat)))
  )

p_corr <- ggplot(corr_long, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))),
            size = 3) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(
    title = "Baseline Correlation Matrix",
    subtitle = "Pairwise complete observations",
    x = NULL, y = NULL
  ) +
  theme_research() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "baseline_correlations.png"), p_corr,
       width = 10, height = 8, dpi = 300)
message("  -> baseline_correlations.png")


# =============================================================================
# SUMMARY
# =============================================================================
message("\n=== COMPLETE ===")
message(sprintf("All outputs in: %s", out_dir))
message(sprintf("  table1.csv / table1.html"))
message(sprintf("  enrollment_funnel.png"))
message(sprintf("  visit_timing.png / temporal_windows.csv"))
message(sprintf("  biomarker_trajectories.png"))
message(sprintf("  biomarker_distributions.png"))
message(sprintf("  clinical_scores.png"))
message(sprintf("  delta_analysis.png / delta_summary.csv"))
message(sprintf("  vitals_trajectories.png"))
message(sprintf("  baseline_correlations.png"))
