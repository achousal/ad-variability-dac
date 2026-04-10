#!/usr/bin/env Rscript
# dac_ratio.R -- pTau217/Ab42-40 ratio variability (Analysis 7)
#
# NOTE: Only 6 patients have paired pTau217 (Lumipulse) + Ab42/40 (Lucent) data
# at both V1 and V2. This is insufficient for inferential statistics. All output
# is descriptive/exploratory only.
#
# The ratio here is ptau217 / ab42_40. Note the Ab42/40 column is the Lucent
# composite ratio (not raw Ab1-42 pg/mL from Lumipulse), so this ratio is
# not directly comparable to the Mayfield et al. pTau217/Ab1-42 ratio.
#
# Usage:
#   Rscript analysis/dac_ratio.R [--data-dir PATH] [--output-dir PATH]
#
# Requires: tidyverse, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(patchwork)
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
out_dir  <- file.path(project_dir, "results", "dac_ratio")

for (i in seq_along(args)) {
  if (args[i] == "--data-dir"   && i < length(args)) data_dir <- args[i + 1]
  if (args[i] == "--output-dir" && i < length(args)) out_dir  <- args[i + 1]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "theme_research.R"))
source(file.path(script_dir, "palettes.R"))

message("Reading data from: ", data_dir)
message("Writing results to: ", out_dir)
message("NOTE: N=6 paired patients -- descriptive only, no inference.")

# -- Load data -----------------------------------------------------------------

long <- read_csv(file.path(data_dir, "dac_long.csv"), show_col_types = FALSE) |>
  mutate(
    visit      = factor(visit, levels = c("V1", "V2", "V3", "V4")),
    blood_date = as.Date(blood_date)
  )

# Patients with >=2 measurements of BOTH ptau217 and ab42_40
paired_ids <- long |>
  filter(!is_redraw) |>
  group_by(record_id) |>
  summarise(
    n_ptau = sum(!is.na(ptau217)),
    n_ab   = sum(!is.na(ab42_40)),
    .groups = "drop"
  ) |>
  filter(n_ptau >= 2, n_ab >= 2) |>
  pull(record_id)

message(sprintf("Patients with >=2 paired ptau217+ab42_40 measurements: %d", length(paired_ids)))

if (length(paired_ids) < 2) {
  message("Insufficient N for analysis. Exiting.")
  quit(status = 0)
}

sub <- long |>
  filter(record_id %in% paired_ids, !is_redraw, !is.na(ptau217), !is.na(ab42_40)) |>
  mutate(ratio = ptau217 / ab42_40)

message(sprintf("Observations: %d rows, %d patients", nrow(sub), n_distinct(sub$record_id)))

# -- Per-patient CV for each component + ratio --------------------------------

cv_df <- sub |>
  group_by(record_id) |>
  summarise(
    n_obs      = n(),
    cv_ptau217 = sd(ptau217)  / mean(ptau217)  * 100,
    cv_ab42_40 = sd(ab42_40)  / mean(ab42_40)  * 100,
    cv_ratio   = sd(ratio)    / mean(ratio)     * 100,
    .groups    = "drop"
  )

message("\nCV summary (N=", nrow(cv_df), " patients):")
cv_summary <- cv_df |>
  summarise(across(starts_with("cv_"), list(
    median = ~ median(.x, na.rm = TRUE),
    q25    = ~ quantile(.x, 0.25, na.rm = TRUE),
    q75    = ~ quantile(.x, 0.75, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))
print(cv_summary)

# Does the ratio reduce CV relative to individual components?
cv_long <- cv_df |>
  pivot_longer(starts_with("cv_"), names_to = "analyte", values_to = "cv_pct") |>
  mutate(
    analyte = recode(analyte,
      "cv_ptau217" = "pTau217",
      "cv_ab42_40" = "Ab42/40",
      "cv_ratio"   = "pTau217/Ab42-40 ratio"
    ),
    analyte = factor(analyte, levels = c("pTau217", "Ab42/40", "pTau217/Ab42-40 ratio"))
  )

message("\nCV by analyte:")
print(cv_long |> group_by(analyte) |>
  summarise(median_cv = median(cv_pct), .groups = "drop"))

# -- Figure 1: CV comparison (paired lines per patient) -----------------------

fig_cv <- ggplot(cv_long, aes(x = analyte, y = cv_pct, group = record_id)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey50") +
  geom_line(alpha = 0.5, color = "#4575b4", linewidth = 0.6) +
  geom_point(aes(color = analyte), size = 2.5) +
  stat_summary(aes(group = 1), fun = median, geom = "point",
               size = 5, shape = 18, color = "black") +
  scale_color_manual(
    values = c("pTau217" = "#E41A1C", "Ab42/40" = "#4DAF4A",
               "pTau217/Ab42-40 ratio" = "#377EB8"),
    guide = "none"
  ) +
  annotate("text", x = 3.4, y = 21, label = "20% TEa",
           color = "grey40", size = 3.2, hjust = 0) +
  labs(
    title    = "Within-Patient CV: Individual Analytes vs. Ratio",
    subtitle = sprintf("N=%d patients with >=2 paired measurements; diamond = median; lines connect same patient", nrow(cv_df)),
    x        = NULL,
    y        = "Within-patient CV (%)"
  ) +
  theme_research() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

fsz <- get_figure_size("box_single")
ggsave(file.path(out_dir, "fig_ratio_cv_comparison.pdf"),
       fig_cv, width = fsz["width"], height = fsz["height"])
ggsave(file.path(out_dir, "fig_ratio_cv_comparison.png"),
       fig_cv, width = fsz["width"], height = fsz["height"], dpi = 300)
message("Saved: fig_ratio_cv_comparison")

# -- Figure 2: Individual trajectories for all 3 analytes ---------------------

traj_df <- sub |>
  select(record_id, visit, days_from_v1, ptau217, ab42_40, ratio) |>
  pivot_longer(c(ptau217, ab42_40, ratio),
               names_to = "analyte", values_to = "value") |>
  mutate(
    analyte = recode(analyte,
      "ptau217" = "pTau217", "ab42_40" = "Ab42/40",
      "ratio"   = "pTau217/Ab42-40 ratio"
    ),
    analyte = factor(analyte, levels = c("pTau217", "Ab42/40", "pTau217/Ab42-40 ratio"))
  )

fig_traj <- ggplot(traj_df, aes(x = days_from_v1, y = value,
                                 group = record_id, color = factor(record_id))) +
  geom_line(alpha = 0.7, linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ analyte, scales = "free_y", ncol = 3) +
  labs(
    title    = "pTau217, Ab42/40, and Ratio Trajectories",
    subtitle = sprintf("N=%d paired patients; each line = one patient", nrow(cv_df)),
    x        = "Days from V1",
    y        = "Value"
  ) +
  theme_research() +
  theme(legend.position = "none")

fsz_m <- get_figure_size("scatter_multi")
ggsave(file.path(out_dir, "fig_ratio_trajectories.pdf"),
       fig_traj, width = fsz_m["width"] * 0.8, height = fsz_m["height"] * 0.45)
ggsave(file.path(out_dir, "fig_ratio_trajectories.png"),
       fig_traj, width = fsz_m["width"] * 0.8, height = fsz_m["height"] * 0.45, dpi = 300)
message("Saved: fig_ratio_trajectories")

# -- Summary table -------------------------------------------------------------

gt_ratio <- cv_long |>
  group_by(analyte) |>
  summarise(
    n          = n(),
    median_cv  = median(cv_pct),
    q25        = quantile(cv_pct, 0.25),
    q75        = quantile(cv_pct, 0.75),
    pct_above_20 = mean(cv_pct > 20) * 100,
    .groups    = "drop"
  ) |>
  mutate(cv_fmt = sprintf("%.1f (%.1f-%.1f)", median_cv, q25, q75)) |>
  select(Analyte = analyte, `N patients` = n,
         `Median CV % (IQR)` = cv_fmt,
         `% > 20% TEa` = pct_above_20) |>
  gt() |>
  tab_header(
    title    = "Within-Patient CV: Ratio vs Individual Analytes",
    subtitle = "Exploratory -- N=6 patients with paired Lumipulse + Lucent measurements"
  ) |>
  tab_footnote(
    footnote = "CAUTION: N=6. All values are descriptive only. Ratio = pTau217 (Lumipulse) / Ab42/40 (Lucent composite). Not directly comparable to Mayfield et al. pTau217/Ab1-42 ratio."
  ) |>
  fmt_number(columns = `% > 20% TEa`, decimals = 0, pattern = "{x}%") |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_ratio, file.path(out_dir, "table_ratio_cv.html"))
gtsave(gt_ratio, file.path(out_dir, "table_ratio_cv.png"), expand = 20)
message("Saved: table_ratio_cv")

write_csv(cv_df, file.path(out_dir, "ratio_cv_per_patient.csv"))

message("\nDone. Outputs written to: ", out_dir)
