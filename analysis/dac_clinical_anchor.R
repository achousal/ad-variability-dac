#!/usr/bin/env Rscript
# dac_clinical_anchor.R -- Clinical anchoring of biomarker variability
#
# NOTE: Longitudinal cognitive data are not available in this REDCap export.
# MoCA/MMSE are baseline-only. Delta-delta analysis (biomarker change vs
# cognitive change) is therefore not possible. This script performs:
#
#   1. Cross-sectional: baseline MoCA vs. baseline biomarker level
#   2. Baseline MoCA vs. per-patient CV (does impairment predict more noise?)
#   3. Baseline MoCA vs. RCI class (V1->V2)
#
# All analyses treat MoCA as a proxy for cognitive status at enrollment, not
# as a longitudinal outcome.
#
# Usage:
#   Rscript analysis/dac_clinical_anchor.R [--data-dir PATH] [--output-dir PATH]
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

data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "results", "dac_clinical_anchor")
icc_cv_dir  <- file.path(project_dir, "results", "dac_icc_cv")
rci_dir     <- file.path(project_dir, "results", "dac_rci")

for (i in seq_along(args)) {
  if (args[i] == "--data-dir"   && i < length(args)) data_dir <- args[i + 1]
  if (args[i] == "--output-dir" && i < length(args)) out_dir  <- args[i + 1]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "theme_research.R"))
source(file.path(script_dir, "palettes.R"))

message("Reading data from: ", data_dir)
message("Writing results to: ", out_dir)
message("NOTE: Longitudinal cognitive data unavailable. Baseline MoCA only.")

# -- Load data -----------------------------------------------------------------

long <- read_csv(file.path(data_dir, "dac_long.csv"), show_col_types = FALSE) |>
  mutate(visit = factor(visit, levels = c("V1", "V2", "V3", "V4")))

baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE)

long_ids <- long |>
  filter(!is_redraw) |>
  group_by(record_id) |>
  filter(n_distinct(visit) >= 2) |>
  pull(record_id) |> unique()

# V1 biomarker values
v1 <- long |>
  filter(visit == "V1", !is_redraw, record_id %in% long_ids) |>
  select(record_id, ptau217, nfl, lucent_ad, ptau217_risk) |>
  distinct(record_id, .keep_all = TRUE)

# Per-patient CV (from dac_icc_cv output)
cv_file <- file.path(icc_cv_dir, "cv_per_patient.csv")
cv_wide <- if (file.exists(cv_file)) {
  read_csv(cv_file, show_col_types = FALSE) |>
    select(record_id, biomarker, cv_pct) |>
    pivot_wider(names_from = biomarker, values_from = cv_pct,
                names_prefix = "cv_")
} else NULL

# RCI class (V1->V2) from dac_rci output
rci_file <- file.path(rci_dir, "rci_per_patient_all_pairs.csv")
rci_v1v2 <- if (file.exists(rci_file)) {
  read_csv(rci_file, show_col_types = FALSE) |>
    filter(visit_pair == "V1->V2") |>
    select(record_id, biomarker, rci_class) |>
    mutate(reliable = rci_class != "No reliable change") |>
    pivot_wider(names_from = biomarker, values_from = c(rci_class, reliable),
                names_sep = "_")
} else NULL

# Join with baseline MoCA
anchor <- baseline |>
  filter(record_id %in% long_ids) |>
  select(record_id, moca_raw, mmse_raw, diagnosis) |>
  left_join(v1,      by = "record_id") |>
  left_join(cv_wide, by = "record_id") |>
  left_join(rci_v1v2, by = "record_id")

message(sprintf("\nAnchor dataset: %d patients", nrow(anchor)))
message(sprintf("  MoCA available: %d", sum(!is.na(anchor$moca_raw))))
message(sprintf("  MMSE available: %d", sum(!is.na(anchor$mmse_raw))))

biomarkers <- list(
  ptau217   = list(label = "pTau217",   cv_col = "cv_ptau217"),
  nfl       = list(label = "NfL",       cv_col = "cv_nfl"),
  lucent_ad = list(label = "Lucent AD", cv_col = "cv_lucent_ad")
)

# -- Analysis 1: Baseline MoCA vs baseline biomarker level --------------------

message("\nBaseline MoCA vs biomarker level (Spearman):")
cor_level <- map_dfr(names(biomarkers), function(bm) {
  sub <- anchor |> filter(!is.na(moca_raw), !is.na(.data[[bm]]))
  if (nrow(sub) < 5) return(tibble(biomarker=bm, n=nrow(sub), rho=NA, p=NA))
  ct <- cor.test(sub$moca_raw, sub[[bm]], method = "spearman", exact = FALSE)
  tibble(biomarker = bm, label = biomarkers[[bm]]$label,
         n = nrow(sub), rho = ct$estimate, p = ct$p.value)
})
print(cor_level)

# -- Analysis 2: Baseline MoCA vs per-patient CV -------------------------------

message("\nBaseline MoCA vs CV (Spearman):")
cor_cv <- map_dfr(names(biomarkers), function(bm) {
  cv_col <- biomarkers[[bm]]$cv_col
  if (is.null(cv_wide) || !cv_col %in% names(anchor)) {
    return(tibble(biomarker=bm, n=0, rho=NA, p=NA))
  }
  sub <- anchor |> filter(!is.na(moca_raw), !is.na(.data[[cv_col]]))
  if (nrow(sub) < 5) return(tibble(biomarker=bm, n=nrow(sub), rho=NA, p=NA))
  ct <- cor.test(sub$moca_raw, sub[[cv_col]], method = "spearman", exact = FALSE)
  tibble(biomarker = bm, label = biomarkers[[bm]]$label,
         n = nrow(sub), rho = ct$estimate, p = ct$p.value)
})
print(cor_cv)

# -- Analysis 3: Baseline MoCA by RCI class ------------------------------------

message("\nBaseline MoCA by RCI class:")
rci_moca <- map_dfr(names(biomarkers), function(bm) {
  cls_col <- paste0("rci_class_", bm)
  if (is.null(rci_v1v2) || !cls_col %in% names(anchor)) return(NULL)
  anchor |>
    filter(!is.na(moca_raw), !is.na(.data[[cls_col]])) |>
    group_by(class = .data[[cls_col]]) |>
    summarise(n = n(), median_moca = median(moca_raw), .groups = "drop") |>
    mutate(biomarker = bm, label = biomarkers[[bm]]$label)
})
if (nrow(rci_moca) > 0) print(rci_moca)

# -- Figure 1: Scatter MoCA vs baseline biomarker level -----------------------

moca_df <- anchor |>
  filter(!is.na(moca_raw)) |>
  pivot_longer(cols = all_of(names(biomarkers)),
               names_to = "biomarker", values_to = "level") |>
  filter(!is.na(level)) |>
  mutate(label = map_chr(biomarker, ~ biomarkers[[.x]]$label))

fig_level <- ggplot(moca_df, aes(x = moca_raw, y = level)) +
  geom_point(alpha = 0.7, size = 2, color = "#4575b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A1C",
              linewidth = 0.8, alpha = 0.15) +
  geom_text(
    data = cor_level |> filter(!is.na(rho)),
    aes(x = -Inf, y = Inf,
        label = sprintf("rho=%.2f, p=%.3f, n=%d", rho, p, n)),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.5,
    size = 3.2, color = "grey30"
  ) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    title    = "Baseline MoCA vs. Baseline Biomarker Level",
    subtitle = "Cross-sectional only -- no longitudinal cognitive data available",
    x        = "MoCA (baseline)",
    y        = "Biomarker level (V1)"
  ) +
  theme_research()

fsz <- get_figure_size("scatter_multi")
ggsave(file.path(out_dir, "fig_moca_vs_level.pdf"),
       fig_level, width = fsz["width"] * 0.8, height = fsz["height"] * 0.5)
ggsave(file.path(out_dir, "fig_moca_vs_level.png"),
       fig_level, width = fsz["width"] * 0.8, height = fsz["height"] * 0.5, dpi = 300)
message("Saved: fig_moca_vs_level")

# -- Figure 2: Baseline MoCA vs per-patient CV ---------------------------------

if (!is.null(cv_wide)) {
  cv_moca_df <- anchor |>
    filter(!is.na(moca_raw)) |>
    pivot_longer(
      cols      = starts_with("cv_"),
      names_to  = "biomarker",
      values_to = "cv_pct",
      names_prefix = "cv_"
    ) |>
    filter(!is.na(cv_pct), biomarker %in% names(biomarkers)) |>
    mutate(label = map_chr(biomarker, ~ biomarkers[[.x]]$label))

  fig_cv_moca <- ggplot(cv_moca_df, aes(x = moca_raw, y = cv_pct)) +
    geom_hline(yintercept = 20, linetype = "dashed", color = "grey50",
               linewidth = 0.5) +
    geom_point(alpha = 0.7, size = 2, color = "#4575b4") +
    geom_smooth(method = "lm", se = TRUE, color = "#E41A1C",
                linewidth = 0.8, alpha = 0.15) +
    geom_text(
      data = cor_cv |> filter(!is.na(rho)),
      aes(x = -Inf, y = Inf,
          label = sprintf("rho=%.2f, p=%.3f, n=%d", rho, p, n)),
      inherit.aes = FALSE, hjust = -0.05, vjust = 1.5,
      size = 3.2, color = "grey30"
    ) +
    facet_wrap(~ label, scales = "free_y", ncol = 3) +
    labs(
      title    = "Baseline MoCA vs. Within-Patient CV",
      subtitle = "Dashed: 20% TEa. Does cognitive impairment predict higher variability?",
      x        = "MoCA (baseline)",
      y        = "Within-patient CV (%)"
    ) +
    theme_research()

  ggsave(file.path(out_dir, "fig_moca_vs_cv.pdf"),
         fig_cv_moca, width = fsz["width"] * 0.8, height = fsz["height"] * 0.5)
  ggsave(file.path(out_dir, "fig_moca_vs_cv.png"),
         fig_cv_moca, width = fsz["width"] * 0.8, height = fsz["height"] * 0.5,
         dpi = 300)
  message("Saved: fig_moca_vs_cv")
}

# -- Summary table -------------------------------------------------------------

tbl_df <- cor_level |>
  rename(rho_level = rho, p_level = p, n_level = n) |>
  left_join(
    cor_cv |> rename(rho_cv = rho, p_cv = p, n_cv = n),
    by = c("biomarker", "label")
  ) |>
  mutate(
    level_fmt = if_else(!is.na(rho_level),
                        sprintf("%.2f (p=%.3f, n=%d)", rho_level, p_level, n_level),
                        "insufficient N"),
    cv_fmt    = if_else(!is.na(rho_cv),
                        sprintf("%.2f (p=%.3f, n=%d)", rho_cv, p_cv, n_cv),
                        "insufficient N")
  )

gt_anchor <- tbl_df |>
  select(
    Biomarker                       = label,
    `Spearman rho: MoCA vs level`   = level_fmt,
    `Spearman rho: MoCA vs CV`      = cv_fmt
  ) |>
  gt() |>
  tab_header(
    title    = "Clinical Anchoring: Baseline MoCA Correlations",
    subtitle = "Longitudinal subcohort (N<=47)"
  ) |>
  tab_footnote(
    footnote = "LIMITATION: No longitudinal cognitive data in current REDCap export. MoCA is baseline-only. Delta-delta analysis (biomarker change vs cognitive change) is not possible."
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_anchor, file.path(out_dir, "table_clinical_anchor.html"))
gtsave(gt_anchor, file.path(out_dir, "table_clinical_anchor.png"), expand = 20)
message("Saved: table_clinical_anchor")

# -- Write CSV -----------------------------------------------------------------

write_csv(cor_level, file.path(out_dir, "cor_moca_vs_level.csv"))
write_csv(cor_cv,    file.path(out_dir, "cor_moca_vs_cv.csv"))

message("\nDone. Outputs written to: ", out_dir)
