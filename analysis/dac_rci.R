#!/usr/bin/env Rscript
# dac_rci.R -- Reliable Change Index analysis for DAC CRF longitudinal biomarker data
#
# Reads dac_long.csv, dac_baseline.csv, and (optionally) ICC results from
# results/dac_icc_cv/icc_lmm_results.csv.
#
# Outputs tables and figures to results/dac_rci/.
#
# Usage:
#   Rscript analysis/dac_rci.R [--data-dir PATH] [--output-dir PATH] [--icc-dir PATH]
#
# Methods:
#   Empirical RCI  -- RCI = (Vx - V1) / SD_diff
#                     SD_diff = empirical SD of Vx-V1 differences for that visit pair.
#                     Computed for V1->V2, V1->V3, and V1->V4 separately.
#                     Threshold: |RCI| > 1.96 => reliable change (p < 0.05, two-tailed).
#
#   Jacobson-Truax RCI (sensitivity) -- RCI_JT = (Vx - V1) / SE_diff
#                     SE_diff = sqrt(2) * SEM,  SEM = SD_V1 * sqrt(1 - ICC)
#                     Requires ICC from dac_icc_cv.R.
#
#   Log-scale RCI (sensitivity for skewed analytes) -- applied to pTau217 and NfL.
#                     RCI_log = log(Vx/V1) / SD_log_diff
#
# Requires: tidyverse, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
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
out_dir  <- file.path(project_dir, "results", "dac_rci")
icc_dir  <- file.path(project_dir, "results", "dac_icc_cv")

for (i in seq_along(args)) {
  if (args[i] == "--data-dir"   && i < length(args)) data_dir <- args[i + 1]
  if (args[i] == "--output-dir" && i < length(args)) out_dir  <- args[i + 1]
  if (args[i] == "--icc-dir"    && i < length(args)) icc_dir  <- args[i + 1]
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
    visit_date = as.Date(visit_date),
    blood_date = as.Date(blood_date)
  )

baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("Female", "Male")))

icc_file <- file.path(icc_dir, "icc_lmm_results.csv")
icc_results <- if (file.exists(icc_file)) {
  message("Loading ICC results from: ", icc_file)
  read_csv(icc_file, show_col_types = FALSE)
} else {
  message("ICC results not found. Jacobson-Truax method will be skipped.")
  NULL
}

# -- Biomarker metadata --------------------------------------------------------

biomarkers <- list(
  ptau217   = list(label = "pTau217",   unit = "pg/mL", log_rci = TRUE),
  nfl       = list(label = "NfL",       unit = "pg/mL", log_rci = TRUE),
  gfap      = list(label = "GFAP",      unit = "pg/mL", log_rci = FALSE),
  ab42_40   = list(label = "Ab42/40",   unit = "ratio",  log_rci = FALSE),
  lucent_ad = list(label = "Lucent AD", unit = "score",  log_rci = FALSE)
)

RCI_THRESHOLD <- 1.96
VISIT_PAIRS   <- c("V2", "V3", "V4")  # each paired against V1

# -- Build V1 anchor -----------------------------------------------------------

v1 <- long |>
  filter(visit == "V1", !is_redraw) |>
  select(record_id, all_of(names(biomarkers)), v1_date = blood_date) |>
  distinct(record_id, .keep_all = TRUE)

# -- Helper: build V1 -> Vx pair table -----------------------------------------
# Deduplicate Vx by keeping non-redraw first, then earliest blood_date.

build_pairs <- function(long_df, v1_df, target_visit) {
  vx <- long_df |>
    filter(visit == target_visit) |>
    arrange(record_id, is_redraw, blood_date) |>
    distinct(record_id, .keep_all = TRUE) |>
    select(record_id, all_of(names(biomarkers)),
           vx_date = blood_date, interval_days = days_from_v1)

  v1_df |>
    inner_join(vx, by = "record_id", suffix = c("_v1", "_vx")) |>
    mutate(target_visit = target_visit)
}

# -- Helper: compute empirical RCI ---------------------------------------------

compute_rci <- function(pairs_df, bm, icc_df = NULL) {
  bm_v1 <- paste0(bm, "_v1")
  bm_vx <- paste0(bm, "_vx")
  meta   <- biomarkers[[bm]]

  sub <- pairs_df |>
    select(record_id, target_visit, interval_days,
           v1 = all_of(bm_v1), vx = all_of(bm_vx)) |>
    filter(!is.na(v1), !is.na(vx))

  if (nrow(sub) < 3) return(NULL)

  diff    <- sub$vx - sub$v1
  sd_diff <- sd(diff)

  sub <- sub |>
    mutate(
      diff_raw   = vx - v1,
      pct_change = diff_raw / v1 * 100,
      rci        = diff_raw / sd_diff,
      rci_class  = case_when(
        rci >  RCI_THRESHOLD ~ "Reliable increase",
        rci < -RCI_THRESHOLD ~ "Reliable decrease",
        TRUE                 ~ "No reliable change"
      )
    )

  # Jacobson-Truax (uses ICC-based SEM)
  if (!is.null(icc_df)) {
    icc_val <- icc_df |> filter(biomarker == bm) |> pull(icc)
    if (length(icc_val) == 1 && !is.na(icc_val)) {
      sem_jt  <- sd(sub$v1) * sqrt(1 - icc_val)
      se_diff_jt <- sqrt(2) * sem_jt
      sub <- sub |>
        mutate(
          rci_jt       = diff_raw / se_diff_jt,
          rci_jt_class = case_when(
            rci_jt >  RCI_THRESHOLD ~ "Reliable increase",
            rci_jt < -RCI_THRESHOLD ~ "Reliable decrease",
            TRUE                    ~ "No reliable change"
          )
        )
    }
  }

  # Log-scale sensitivity
  if (meta$log_rci && all(sub$v1 > 0, sub$vx > 0)) {
    log_diff    <- log(sub$vx) - log(sub$v1)
    sd_log_diff <- sd(log_diff)
    sub <- sub |>
      mutate(
        rci_log       = log(vx / v1) / sd_log_diff,
        rci_log_class = case_when(
          rci_log >  RCI_THRESHOLD ~ "Reliable increase",
          rci_log < -RCI_THRESHOLD ~ "Reliable decrease",
          TRUE                     ~ "No reliable change"
        )
      )
  }

  list(
    data    = sub |> mutate(biomarker = bm, label = meta$label),
    sd_diff = sd_diff,
    n_pairs = nrow(sub)
  )
}

# -- Run for all visit pairs x biomarkers --------------------------------------

message("\nComputing RCI for V1->V2, V1->V3, V1->V4...")

all_pairs <- map(VISIT_PAIRS, ~ build_pairs(long, v1, .x)) |>
  set_names(VISIT_PAIRS)

rci_all <- map_dfr(VISIT_PAIRS, function(vx_label) {
  pairs_vx <- all_pairs[[vx_label]]
  message(sprintf("\n  --- V1 -> %s (N=%d patients) ---", vx_label, nrow(pairs_vx)))
  map_dfr(names(biomarkers), function(bm) {
    res <- compute_rci(pairs_vx, bm, icc_results)
    if (is.null(res)) return(NULL)
    message(sprintf("    %s: %d pairs, SD_diff=%.3g, reliable changers: %d",
                    biomarkers[[bm]]$label, res$n_pairs, res$sd_diff,
                    sum(res$data$rci_class != "No reliable change")))
    res$data |>
      mutate(visit_pair = paste0("V1->", vx_label),
             sd_diff    = res$sd_diff,
             n_pair_bm  = res$n_pairs)
  })
})

# -- Classification summary (per biomarker x visit pair) -----------------------

rci_summary <- rci_all |>
  group_by(visit_pair, biomarker, label, sd_diff, n_pair_bm) |>
  summarise(
    n_increase   = sum(rci_class == "Reliable increase"),
    n_decrease   = sum(rci_class == "Reliable decrease"),
    n_reliable   = n_increase + n_decrease,
    pct_reliable = n_reliable / n() * 100,
    .groups      = "drop"
  )

message("\nRCI summary (all visit pairs):")
print(rci_summary |> select(visit_pair, label, n_pair_bm, sd_diff, n_reliable, pct_reliable))

# -- Table: RCI summary across visit pairs ------------------------------------

message("\nBuilding summary table...")

tbl_wide <- rci_summary |>
  mutate(
    cell = sprintf("%d/%d (%.0f%%)\n+%d/-%d",
                   n_reliable, n_pair_bm, pct_reliable, n_increase, n_decrease)
  ) |>
  select(label, visit_pair, cell) |>
  pivot_wider(names_from = visit_pair, values_from = cell, values_fill = "N/A")

gt_rci <- tbl_wide |>
  gt() |>
  tab_header(
    title    = "Reliable Change Index: Proportion of Reliable Changers by Visit Pair",
    subtitle = "DAC CRF cohort -- empirical RCI (|RCI| > 1.96, two-tailed p < 0.05)"
  ) |>
  tab_footnote(
    footnote = "Each cell: reliable changers / N pairs (%), increases / decreases. SD_diff computed separately per visit pair from all available difference scores. |RCI| > 1.96 threshold."
  ) |>
  cols_label(label = "Biomarker") |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_rci, file.path(out_dir, "table_rci_summary.html"))
gtsave(gt_rci, file.path(out_dir, "table_rci_summary.png"), expand = 20)
message("Saved: table_rci_summary")

# -- Figure 1: % reliable changers by visit pair (line plot) ------------------

# Only plot biomarkers with >=3 pairs in at least one visit comparison
bm_with_data <- rci_summary |>
  filter(n_pair_bm >= 3) |>
  pull(label) |> unique()

fig_trend <- rci_summary |>
  filter(label %in% bm_with_data) |>
  mutate(
    visit_pair = factor(visit_pair, levels = paste0("V1->", VISIT_PAIRS)),
    label      = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label))
  ) |>
  ggplot(aes(x = visit_pair, y = pct_reliable, group = label, color = label)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_pair_bm), alpha = 0.85) +
  scale_size_continuous(range = c(2, 6), name = "N pairs") +
  scale_color_manual(
    values = c(
      "pTau217"   = "#E41A1C",
      "NfL"       = "#377EB8",
      "GFAP"      = "#4DAF4A",
      "Ab42/40"   = "#984EA3",
      "Lucent AD" = "#FF7F00"
    ),
    name = "Biomarker"
  ) +
  labs(
    title    = "Reliable Changers by Follow-up Interval",
    subtitle = "% patients with |RCI| > 1.96; point size = N pairs available",
    x        = "Visit comparison",
    y        = "% reliable changers"
  ) +
  theme_research()

fsz <- get_figure_size("scatter_bivar")
ggsave(file.path(out_dir, "fig_rci_by_visit_pair.pdf"),
       fig_trend, width = fsz["width"], height = fsz["height"])
ggsave(file.path(out_dir, "fig_rci_by_visit_pair.png"),
       fig_trend, width = fsz["width"], height = fsz["height"], dpi = 300)
message("Saved: fig_rci_by_visit_pair")

# -- Figure 2: RCI distributions per biomarker, faceted by visit pair ---------

rci_colors <- c(
  "Reliable increase"  = "#E41A1C",
  "No reliable change" = "#999999",
  "Reliable decrease"  = "#377EB8"
)

bm_order <- map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)

dist_df <- rci_all |>
  filter(label %in% bm_with_data) |>
  mutate(
    label      = factor(label, levels = bm_order),
    visit_pair = factor(visit_pair, levels = paste0("V1->", VISIT_PAIRS)),
    rci_class  = factor(rci_class, levels = names(rci_colors))
  )

fig_dist <- ggplot(dist_df, aes(x = rci, fill = rci_class, color = rci_class)) +
  geom_vline(xintercept = c(-RCI_THRESHOLD, RCI_THRESHOLD),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_histogram(bins = 12, alpha = 0.7, position = "identity") +
  scale_fill_manual(values  = rci_colors, name = NULL) +
  scale_color_manual(values = rci_colors, name = NULL) +
  facet_grid(label ~ visit_pair, scales = "free") +
  labs(
    title    = "RCI Distributions by Biomarker and Visit Pair",
    subtitle = "Dashed lines: RCI = +/-1.96 threshold; rows = biomarker, columns = visit comparison",
    x        = "RCI",
    y        = "N patients"
  ) +
  theme_research(base_size = 11) +
  theme(legend.position = "bottom")

fsz_g <- get_figure_size("heatmap")
n_bm  <- length(unique(dist_df$label))
n_vp  <- length(VISIT_PAIRS)
ggsave(file.path(out_dir, "fig_rci_distributions.pdf"),
       fig_dist, width = fsz_g["width"] * (n_vp / 2), height = fsz_g["height"] * (n_bm / 3))
ggsave(file.path(out_dir, "fig_rci_distributions.png"),
       fig_dist, width = fsz_g["width"] * (n_vp / 2), height = fsz_g["height"] * (n_bm / 3),
       dpi = 300)
message("Saved: fig_rci_distributions")

# -- Figure 3: SD_diff by visit pair (shows whether noise floor grows) ---------

sd_trend_df <- rci_summary |>
  filter(label %in% bm_with_data) |>
  mutate(
    visit_pair = factor(visit_pair, levels = paste0("V1->", VISIT_PAIRS)),
    label      = factor(label, levels = bm_order)
  )

fig_sd <- ggplot(sd_trend_df,
                 aes(x = visit_pair, y = sd_diff, group = label, color = label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "pTau217"   = "#E41A1C",
      "NfL"       = "#377EB8",
      "GFAP"      = "#4DAF4A",
      "Ab42/40"   = "#984EA3",
      "Lucent AD" = "#FF7F00"
    ),
    name = "Biomarker"
  ) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    title    = "SD of Difference Scores by Follow-up Interval",
    subtitle = "Growing SD_diff means noise accumulates; stable SD_diff means oscillatory pattern persists",
    x        = "Visit comparison",
    y        = "SD_diff (raw units)"
  ) +
  theme_research() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fig_sd_diff_by_visit_pair.pdf"),
       fig_sd, width = fsz_g["width"] * 1.2, height = fsz_g["height"] * 0.8)
ggsave(file.path(out_dir, "fig_sd_diff_by_visit_pair.png"),
       fig_sd, width = fsz_g["width"] * 1.2, height = fsz_g["height"] * 0.8, dpi = 300)
message("Saved: fig_sd_diff_by_visit_pair")

# -- Figure 4: V1->V4 spaghetti (full longitudinal, reliable changers flagged) -

make_spaghetti_full <- function(bm) {
  meta <- biomarkers[[bm]]

  # Get RCI class from V1->V4 (or V1->V3 if V4 not available)
  rci_class_df <- rci_all |>
    filter(biomarker == bm) |>
    arrange(desc(visit_pair)) |>          # V1->V4 first
    distinct(record_id, .keep_all = TRUE) |>
    select(record_id, rci_class_final = rci_class)

  if (nrow(rci_class_df) == 0) return(NULL)

  plot_data <- long |>
    filter(
      record_id %in% rci_class_df$record_id,
      !is.na(.data[[bm]]),
      !is_redraw
    ) |>
    left_join(rci_class_df, by = "record_id", relationship = "many-to-one") |>
    mutate(
      rci_class_final = factor(rci_class_final, levels = names(rci_colors)),
      visit           = droplevels(visit)
    )

  n_reliable <- sum(rci_class_df$rci_class_final != "No reliable change")

  ggplot(
    plot_data |> arrange(rci_class_final == "No reliable change"),
    aes(x = visit, y = .data[[bm]], group = record_id,
        color = rci_class_final, alpha = rci_class_final)
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    scale_color_manual(values = rci_colors, name = "RCI class (latest pair)") +
    scale_alpha_manual(
      values = c("Reliable increase" = 1, "No reliable change" = 0.25,
                 "Reliable decrease" = 1),
      guide  = "none"
    ) +
    labs(
      title    = meta$label,
      subtitle = sprintf("N=%d; %d reliable changers (latest pair)", nrow(rci_class_df), n_reliable),
      x        = "Visit",
      y        = sprintf("%s (%s)", meta$label, meta$unit)
    ) +
    theme_research()
}

spag_bms <- intersect(c("ptau217", "nfl", "lucent_ad"), names(biomarkers))
spag_plots <- compact(map(spag_bms, make_spaghetti_full))

if (length(spag_plots) > 0) {
  fig_spag <- wrap_plots(spag_plots, ncol = length(spag_plots)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  fsz_s <- get_figure_size("scatter_multi")
  ggsave(file.path(out_dir, "fig_rci_spaghetti_full.pdf"),
         fig_spag, width = fsz_s["width"], height = fsz_s["height"] * 0.55)
  ggsave(file.path(out_dir, "fig_rci_spaghetti_full.png"),
         fig_spag, width = fsz_s["width"], height = fsz_s["height"] * 0.55, dpi = 300)
  message("Saved: fig_rci_spaghetti_full")
}

# -- Write CSVs ----------------------------------------------------------------

write_csv(rci_all,     file.path(out_dir, "rci_per_patient_all_pairs.csv"))
write_csv(rci_summary, file.path(out_dir, "rci_summary_all_pairs.csv"))

message("\nDone. Outputs written to: ", out_dir)
