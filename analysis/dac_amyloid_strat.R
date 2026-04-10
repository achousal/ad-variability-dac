#!/usr/bin/env Rscript
# dac_amyloid_strat.R -- Amyloid-stratified variability analysis, DAC CRF cohort
#
# Stratifies by baseline pTau217 risk category (primary) and Lucent AD
# interpretation (secondary). Compares within-patient CV, RCI exceedance
# rates, and trajectory monotonicity between amyloid risk groups.
#
# Note: PET is not usable in the longitudinal subcohort (50/51 missing).
#
# Usage:
#   Rscript analysis/dac_amyloid_strat.R [--data-dir PATH] [--output-dir PATH]
#
# Requires: tidyverse, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(forcats)
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
out_dir  <- file.path(project_dir, "results", "dac_amyloid_strat")

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

baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE)

# Longitudinal subcohort (>=2 visits, non-redraw)
long_ids <- long |>
  filter(!is_redraw) |>
  group_by(record_id) |>
  filter(n_distinct(visit) >= 2) |>
  pull(record_id) |>
  unique()

long_sub <- long |> filter(record_id %in% long_ids, !is_redraw)

# Baseline amyloid classifiers from V1
v1_amyloid <- long |>
  filter(visit == "V1", !is_redraw, record_id %in% long_ids) |>
  select(record_id, ptau217_risk, lucent_interpretation) |>
  distinct(record_id, .keep_all = TRUE)

message(sprintf("Longitudinal subcohort: %d patients", length(long_ids)))
message("pTau217 risk distribution:")
print(table(v1_amyloid$ptau217_risk, useNA = "ifany"))
message("Lucent interpretation distribution:")
print(table(v1_amyloid$lucent_interpretation, useNA = "ifany"))

# -- Biomarker metadata --------------------------------------------------------

biomarkers <- list(
  ptau217   = list(label = "pTau217",   unit = "pg/mL"),
  nfl       = list(label = "NfL",       unit = "pg/mL"),
  lucent_ad = list(label = "Lucent AD", unit = "score")
)

# Focus on pTau217, NfL, Lucent AD (adequate N in longitudinal subcohort)

RCI_THRESHOLD <- 1.96

# -- Compute per-patient CV (reuse same logic as dac_icc_cv.R) ----------------

cv_df <- long_sub |>
  pivot_longer(cols = names(biomarkers), names_to = "biomarker", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(record_id, biomarker) |>
  filter(n() >= 2) |>
  summarise(
    n_obs    = n(),
    mean_val = mean(value),
    cv_pct   = sd(value) / mean(value) * 100,
    .groups  = "drop"
  ) |>
  left_join(v1_amyloid, by = "record_id") |>
  mutate(label = map_chr(biomarker, ~ biomarkers[[.x]]$label))

# -- Compute V1->V2 RCI (reuse logic from dac_rci.R) -------------------------

v1 <- long_sub |>
  filter(visit == "V1") |>
  select(record_id, all_of(names(biomarkers)), v1_date = blood_date) |>
  distinct(record_id, .keep_all = TRUE)

v2 <- long_sub |>
  filter(visit == "V2") |>
  arrange(record_id, blood_date) |>
  distinct(record_id, .keep_all = TRUE) |>
  select(record_id, all_of(names(biomarkers)), v2_date = blood_date)

pairs <- v1 |>
  inner_join(v2, by = "record_id", suffix = c("_v1", "_v2")) |>
  left_join(v1_amyloid, by = "record_id")

rci_df <- map_dfr(names(biomarkers), function(bm) {
  sub <- pairs |>
    select(record_id, ptau217_risk, lucent_interpretation,
           v1 = paste0(bm, "_v1"), vx = paste0(bm, "_v2")) |>
    filter(!is.na(v1), !is.na(vx))
  if (nrow(sub) < 3) return(NULL)
  sd_diff <- sd(sub$vx - sub$v1)
  sub |>
    mutate(
      biomarker  = bm,
      label      = biomarkers[[bm]]$label,
      rci        = (vx - v1) / sd_diff,
      rci_class  = case_when(
        rci >  RCI_THRESHOLD ~ "Reliable increase",
        rci < -RCI_THRESHOLD ~ "Reliable decrease",
        TRUE                 ~ "No reliable change"
      ),
      reliable   = rci_class != "No reliable change"
    )
})

# -- Monotonicity per patient --------------------------------------------------
# Classify each patient's pTau217 and NfL trajectory as:
#   monotone_up, monotone_down, non-monotone (bouncing)
# Requires >=3 visits.

classify_monotone <- function(values) {
  v <- values[!is.na(values)]
  if (length(v) < 3) return(NA_character_)
  diffs <- diff(v)
  if (all(diffs >= 0)) "Monotone up"
  else if (all(diffs <= 0)) "Monotone down"
  else "Non-monotone"
}

mono_df <- long_sub |>
  arrange(record_id, visit) |>
  group_by(record_id) |>
  filter(n_distinct(visit) >= 3) |>
  summarise(
    mono_ptau217   = classify_monotone(ptau217),
    mono_nfl       = classify_monotone(nfl),
    mono_lucent    = classify_monotone(lucent_ad),
    .groups        = "drop"
  ) |>
  left_join(v1_amyloid, by = "record_id")

message(sprintf("\nPatients with >=3 visits for monotonicity: %d", nrow(mono_df)))

# -- Stratification groups -----------------------------------------------------
# Primary:  pTau217 High vs Low (drop Intermediate for cleaner contrast)
# Secondary: Lucent High vs Low

risk_colors <- c("High" = "#E41A1C", "Intermediate" = "#FF7F00", "Low" = "#4DAF4A")

# -- Analysis 1: CV by pTau217 risk group -------------------------------------

cv_risk <- cv_df |>
  filter(!is.na(ptau217_risk)) |>
  mutate(ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High")))

# Kruskal-Wallis + pairwise Wilcoxon (High vs Low)
kw_cv <- cv_risk |>
  group_by(label) |>
  summarise(
    kw_p = tryCatch(
      kruskal.test(cv_pct ~ ptau217_risk)$p.value,
      error = function(e) NA_real_
    ),
    n_high = sum(ptau217_risk == "High"),
    n_int  = sum(ptau217_risk == "Intermediate"),
    n_low  = sum(ptau217_risk == "Low"),
    med_high = median(cv_pct[ptau217_risk == "High"], na.rm = TRUE),
    med_low  = median(cv_pct[ptau217_risk == "Low"],  na.rm = TRUE),
    wil_p_hl = tryCatch(
      wilcox.test(cv_pct[ptau217_risk == "High"],
                  cv_pct[ptau217_risk == "Low"],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

message("\nCV by pTau217 risk (Kruskal-Wallis):")
print(kw_cv)

# -- Analysis 2: RCI exceedance by pTau217 risk --------------------------------

rci_risk <- rci_df |>
  filter(!is.na(ptau217_risk)) |>
  mutate(ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High")))

rci_exceedance <- rci_risk |>
  group_by(label, ptau217_risk) |>
  summarise(
    n          = n(),
    n_reliable = sum(reliable),
    pct        = mean(reliable) * 100,
    .groups    = "drop"
  )

message("\nRCI exceedance by pTau217 risk:")
print(rci_exceedance)

# Fisher exact: High vs Low
fisher_rci <- rci_risk |>
  filter(ptau217_risk %in% c("High", "Low")) |>
  group_by(label) |>
  summarise(
    fisher_p = tryCatch({
      tab <- table(ptau217_risk, reliable)
      if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
      fisher.test(tab)$p.value
    }, error = function(e) NA_real_),
    .groups = "drop"
  )

message("\nFisher exact (RCI High vs Low):")
print(fisher_rci)

# -- Analysis 3: Monotonicity by pTau217 risk ----------------------------------

mono_risk <- mono_df |>
  filter(!is.na(ptau217_risk)) |>
  mutate(ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High")))

mono_summary <- bind_rows(
  mono_risk |>
    filter(!is.na(mono_ptau217)) |>
    count(ptau217_risk, class = mono_ptau217) |>
    mutate(label = "pTau217"),
  mono_risk |>
    filter(!is.na(mono_nfl)) |>
    count(ptau217_risk, class = mono_nfl) |>
    mutate(label = "NfL")
) |>
  group_by(label, ptau217_risk) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup()

message("\nMonotonicity by pTau217 risk:")
print(mono_summary)

# -- Figure 1: CV by pTau217 risk (violin + jitter) ----------------------------

fig_cv_risk <- cv_risk |>
  ggplot(aes(x = ptau217_risk, y = cv_pct, fill = ptau217_risk, color = ptau217_risk)) +
  geom_violin(alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.6) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  scale_fill_manual(values  = risk_colors, guide = "none") +
  scale_color_manual(values = risk_colors, guide = "none") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  # Annotate KW p-values
  geom_text(
    data = kw_cv |>
      mutate(
        x     = 2,
        y     = Inf,
        ptau217_risk = factor("Intermediate", levels = c("Low", "Intermediate", "High")),
        label_txt = sprintf("KW p=%.3f", kw_p)
      ),
    aes(x = x, y = y, label = label_txt),
    inherit.aes = FALSE, vjust = 1.5, size = 3.2, color = "grey30"
  ) +
  labs(
    title    = "Within-Patient CV by Baseline pTau217 Risk",
    subtitle = "Dashed: 20% TEa; KW = Kruskal-Wallis p-value",
    x        = "pTau217 risk category",
    y        = "Within-patient CV (%)"
  ) +
  theme_research()

fsz <- get_figure_size("violin_grouped")
ggsave(file.path(out_dir, "fig_cv_by_risk.pdf"),
       fig_cv_risk, width = fsz["width"] * 1.3, height = fsz["height"])
ggsave(file.path(out_dir, "fig_cv_by_risk.png"),
       fig_cv_risk, width = fsz["width"] * 1.3, height = fsz["height"], dpi = 300)
message("Saved: fig_cv_by_risk")

# -- Figure 2: RCI distribution by pTau217 risk --------------------------------

rci_colors_cls <- c(
  "Reliable increase"  = "#E41A1C",
  "No reliable change" = "#999999",
  "Reliable decrease"  = "#377EB8"
)

fig_rci_risk <- rci_risk |>
  mutate(rci_class = factor(rci_class, levels = names(rci_colors_cls))) |>
  ggplot(aes(x = rci, fill = rci_class, color = rci_class)) +
  geom_vline(xintercept = c(-RCI_THRESHOLD, RCI_THRESHOLD),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_histogram(bins = 10, alpha = 0.7, position = "identity") +
  scale_fill_manual(values  = rci_colors_cls, name = NULL) +
  scale_color_manual(values = rci_colors_cls, name = NULL) +
  facet_grid(ptau217_risk ~ label, scales = "free_x") +
  labs(
    title    = "RCI Distribution by pTau217 Risk Category",
    subtitle = "Rows = risk group (Low/Intermediate/High), columns = biomarker",
    x        = "RCI (V1->V2)",
    y        = "N patients"
  ) +
  theme_research(base_size = 11) +
  theme(legend.position = "bottom")

fsz_g <- get_figure_size("heatmap")
ggsave(file.path(out_dir, "fig_rci_by_risk.pdf"),
       fig_rci_risk, width = fsz_g["width"] * 1.2, height = fsz_g["height"])
ggsave(file.path(out_dir, "fig_rci_by_risk.png"),
       fig_rci_risk, width = fsz_g["width"] * 1.2, height = fsz_g["height"], dpi = 300)
message("Saved: fig_rci_by_risk")

# -- Figure 3: Monotonicity stacked bar by pTau217 risk ------------------------

mono_colors <- c(
  "Monotone up"   = "#E41A1C",
  "Non-monotone"  = "#999999",
  "Monotone down" = "#377EB8"
)

if (nrow(mono_summary) > 0) {
  fig_mono <- mono_summary |>
    mutate(class = factor(class, levels = names(mono_colors))) |>
    ggplot(aes(x = ptau217_risk, y = pct, fill = class)) +
    geom_col(position = "stack", width = 0.6) +
    geom_text(aes(label = sprintf("n=%d", n)),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = mono_colors, name = "Trajectory") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~ label, ncol = 2) +
    labs(
      title    = "Trajectory Monotonicity by pTau217 Risk",
      subtitle = "Patients with >=3 visits; proportions within risk group",
      x        = "pTau217 risk category",
      y        = "% patients"
    ) +
    theme_research()

  fsz_b <- get_figure_size("bar")
  ggsave(file.path(out_dir, "fig_monotonicity_by_risk.pdf"),
         fig_mono, width = fsz_b["width"], height = fsz_b["height"])
  ggsave(file.path(out_dir, "fig_monotonicity_by_risk.png"),
         fig_mono, width = fsz_b["width"], height = fsz_b["height"], dpi = 300)
  message("Saved: fig_monotonicity_by_risk")
}

# -- Figure 4: pTau217 spaghetti stratified by risk ----------------------------

spag_data <- long_sub |>
  filter(!is.na(ptau217)) |>
  select(-ptau217_risk) |>
  left_join(v1_amyloid |> select(record_id, ptau217_risk), by = "record_id") |>
  filter(!is.na(ptau217_risk)) |>
  mutate(ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High")))

fig_spag <- ggplot(spag_data,
                   aes(x = visit, y = ptau217, group = record_id,
                       color = ptau217_risk)) +
  geom_line(alpha = 0.45, linewidth = 0.5) +
  geom_point(alpha = 0.6, size = 1.5) +
  stat_summary(aes(group = ptau217_risk), fun = median,
               geom = "line", linewidth = 1.5, linetype = "solid") +
  scale_color_manual(values = risk_colors, name = "pTau217 risk") +
  facet_wrap(~ ptau217_risk, ncol = 3) +
  labs(
    title    = "pTau217 Trajectories by Baseline Risk Category",
    subtitle = "Individual lines (faded); thick line = group median",
    x        = "Visit",
    y        = "pTau217 (pg/mL)"
  ) +
  theme_research() +
  theme(legend.position = "none")

fsz_s <- get_figure_size("scatter_multi")
ggsave(file.path(out_dir, "fig_ptau217_spaghetti_by_risk.pdf"),
       fig_spag, width = fsz_s["width"] * 0.8, height = fsz_s["height"] * 0.5)
ggsave(file.path(out_dir, "fig_ptau217_spaghetti_by_risk.png"),
       fig_spag, width = fsz_s["width"] * 0.8, height = fsz_s["height"] * 0.5, dpi = 300)
message("Saved: fig_ptau217_spaghetti_by_risk")

# -- Table: summary across analyses -------------------------------------------

tbl_summary <- kw_cv |>
  select(label, n_high, n_int, n_low, med_high, med_low, kw_p, wil_p_hl) |>
  left_join(
    rci_exceedance |>
      select(label, ptau217_risk, pct) |>
      pivot_wider(names_from = ptau217_risk, values_from = pct,
                  names_prefix = "rci_pct_"),
    by = "label"
  ) |>
  left_join(fisher_rci, by = "label") |>
  mutate(
    med_fmt    = sprintf("%.1f vs %.1f", med_high, med_low),
    kw_p_fmt   = sprintf("%.3f", kw_p),
    wil_p_fmt  = sprintf("%.3f", wil_p_hl),
    fisher_fmt = sprintf("%.3f", fisher_p),
    rci_high   = sprintf("%.0f%%", coalesce(rci_pct_High, 0)),
    rci_low    = sprintf("%.0f%%", coalesce(rci_pct_Low, 0))
  )

gt_summary <- tbl_summary |>
  select(
    Biomarker     = label,
    `N High/Int/Low` = n_high,
    `Median CV% (High vs Low)` = med_fmt,
    `KW p` = kw_p_fmt,
    `Wilcox p (H vs L)` = wil_p_fmt,
    `% RCI+ High` = rci_high,
    `% RCI+ Low`  = rci_low,
    `Fisher p (H vs L)` = fisher_fmt
  ) |>
  gt() |>
  tab_header(
    title    = "Biomarker Variability by Baseline pTau217 Risk",
    subtitle = "Within-patient CV and RCI exceedance stratified by amyloid risk group"
  ) |>
  tab_footnote(
    footnote = "KW = Kruskal-Wallis across Low/Intermediate/High. Wilcox and Fisher = High vs Low only. RCI threshold: |RCI| > 1.96 (V1->V2)."
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_summary, file.path(out_dir, "table_stratification_summary.html"))
gtsave(gt_summary, file.path(out_dir, "table_stratification_summary.png"), expand = 20)
message("Saved: table_stratification_summary")

# -- Write CSVs ----------------------------------------------------------------

write_csv(cv_risk,        file.path(out_dir, "cv_by_risk.csv"))
write_csv(rci_risk,       file.path(out_dir, "rci_by_risk.csv"))
write_csv(rci_exceedance, file.path(out_dir, "rci_exceedance_by_risk.csv"))
write_csv(mono_summary,   file.path(out_dir, "monotonicity_by_risk.csv"))
write_csv(kw_cv,          file.path(out_dir, "cv_stats_by_risk.csv"))

message("\nDone. Outputs written to: ", out_dir)
