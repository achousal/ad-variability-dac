#!/usr/bin/env Rscript
# dac_varstudy.R -- Variability sub-study cohort analysis
#
# Isolates the 40 patients who consented to the tighter-cadence variability
# sub-study. Answers two questions:
#
#   1. Interval-stratified variability: does CV/RCI compress at very short
#      intervals (<=30d), or is the biological noise floor interval-independent?
#      Interval groups: Short (<=30d), Medium (31-90d), Long (>90d) by V1->V2.
#
#   2. 4-visit trajectory analysis: the 23 patients with complete 4-visit data
#      provide the richest individual trajectories. Monotonicity, CV, and
#      fitted trajectories for this subgroup.
#
# Usage:
#   Rscript analysis/dac_varstudy.R [--data-dir PATH] [--output-dir PATH]
#
# Requires: tidyverse, lme4, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(forcats)
  library(patchwork)
  library(lme4)
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
out_dir  <- file.path(project_dir, "results", "dac_varstudy")

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

var_ids <- baseline |> filter(variability_consented == TRUE) |> pull(record_id)
message(sprintf("Variability sub-study patients: %d", length(var_ids)))

sub <- long |>
  filter(record_id %in% var_ids, !is_redraw) |>
  mutate(years_from_v1 = days_from_v1 / 365.25)

# V1->V2 interval group per patient
v1v2_interval <- sub |>
  filter(visit == "V2") |>
  arrange(record_id, blood_date) |>
  distinct(record_id, .keep_all = TRUE) |>
  transmute(
    record_id,
    v1v2_days = days_from_v1,
    interval_group = cut(
      days_from_v1,
      breaks = c(0, 30, 90, Inf),
      labels = c("Short (<=30d)", "Medium (31-90d)", "Long (>90d)"),
      right  = TRUE
    )
  )

sub <- sub |> left_join(v1v2_interval, by = "record_id")

# 4-visit subgroup
ids_4v <- sub |>
  group_by(record_id) |>
  filter(n_distinct(visit) == 4) |>
  pull(record_id) |> unique()

sub4 <- sub |> filter(record_id %in% ids_4v)

message(sprintf("  4-visit subgroup: %d patients", length(ids_4v)))
message("\nInterval group distribution:")
print(v1v2_interval |> count(interval_group))

# V1 amyloid risk
v1_risk <- sub |>
  filter(visit == "V1") |>
  select(record_id, ptau217_risk) |>
  distinct(record_id, .keep_all = TRUE)

# -- Biomarkers ----------------------------------------------------------------

biomarkers <- list(
  ptau217   = list(label = "pTau217",   unit = "pg/mL"),
  nfl       = list(label = "NfL",       unit = "pg/mL"),
  lucent_ad = list(label = "Lucent AD", unit = "score")
)

RCI_THRESHOLD <- 1.96
interval_colors <- c(
  "Short (<=30d)"  = "#377EB8",
  "Medium (31-90d)"= "#FF7F00",
  "Long (>90d)"    = "#E41A1C"
)

# -- Part 1: CV by interval group ---------------------------------------------

cv_df <- sub |>
  pivot_longer(cols = names(biomarkers), names_to = "biomarker", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(record_id, biomarker, interval_group) |>
  filter(n() >= 2) |>
  summarise(
    n_obs    = n(),
    mean_val = mean(value),
    cv_pct   = sd(value) / mean(value) * 100,
    .groups  = "drop"
  ) |>
  mutate(label = map_chr(biomarker, ~ biomarkers[[.x]]$label))

cv_summary <- cv_df |>
  filter(!is.na(interval_group)) |>
  group_by(label, interval_group) |>
  summarise(
    n         = n(),
    median_cv = median(cv_pct),
    q25       = quantile(cv_pct, 0.25),
    q75       = quantile(cv_pct, 0.75),
    pct_above_20 = mean(cv_pct > 20) * 100,
    .groups   = "drop"
  )

message("\nCV by interval group:")
print(cv_summary |> select(label, interval_group, n, median_cv, pct_above_20))

# Kruskal-Wallis across interval groups
kw_interval <- cv_df |>
  filter(!is.na(interval_group)) |>
  group_by(label) |>
  summarise(
    kw_p = tryCatch(
      kruskal.test(cv_pct ~ interval_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )
message("\nKruskal-Wallis CV ~ interval group:")
print(kw_interval)

# -- Part 2: RCI by interval group --------------------------------------------

v1 <- sub |>
  filter(visit == "V1") |>
  select(record_id, all_of(names(biomarkers)), v1_date = blood_date) |>
  distinct(record_id, .keep_all = TRUE)

v2 <- sub |>
  filter(visit == "V2") |>
  arrange(record_id, blood_date) |>
  distinct(record_id, .keep_all = TRUE) |>
  select(record_id, all_of(names(biomarkers)), v2_date = blood_date)

pairs <- v1 |>
  inner_join(v2, by = "record_id", suffix = c("_v1", "_v2")) |>
  left_join(v1v2_interval, by = "record_id")

rci_df <- map_dfr(names(biomarkers), function(bm) {
  sub_p <- pairs |>
    select(record_id, interval_group,
           v1 = paste0(bm, "_v1"), vx = paste0(bm, "_v2")) |>
    filter(!is.na(v1), !is.na(vx))
  if (nrow(sub_p) < 3) return(NULL)
  sd_diff <- sd(sub_p$vx - sub_p$v1)
  sub_p |>
    mutate(
      biomarker = bm,
      label     = biomarkers[[bm]]$label,
      rci       = (vx - v1) / sd_diff,
      sd_diff   = sd_diff,
      rci_class = case_when(
        rci >  RCI_THRESHOLD ~ "Reliable increase",
        rci < -RCI_THRESHOLD ~ "Reliable decrease",
        TRUE                 ~ "No reliable change"
      )
    )
})

rci_summary <- rci_df |>
  filter(!is.na(interval_group)) |>
  group_by(label, interval_group) |>
  summarise(
    n          = n(),
    n_reliable = sum(rci_class != "No reliable change"),
    pct        = mean(rci_class != "No reliable change") * 100,
    sd_diff    = first(sd_diff),
    .groups    = "drop"
  )

message("\nRCI exceedance by interval group:")
print(rci_summary)

# -- Part 3: 4-visit trajectory analysis --------------------------------------

# Monotonicity
classify_mono <- function(values) {
  v <- values[!is.na(values)]
  if (length(v) < 3) return(NA_character_)
  diffs <- diff(v)
  if (all(diffs >= 0)) "Monotone up"
  else if (all(diffs <= 0)) "Monotone down"
  else "Non-monotone"
}

mono_4v <- sub4 |>
  arrange(record_id, visit) |>
  group_by(record_id) |>
  summarise(
    mono_ptau217   = classify_mono(ptau217),
    mono_nfl       = classify_mono(nfl),
    mono_lucent    = classify_mono(lucent_ad),
    max_follow_days = max(days_from_v1),
    .groups        = "drop"
  ) |>
  left_join(v1_risk, by = "record_id") |>
  left_join(v1v2_interval |> select(record_id, interval_group), by = "record_id")

mono_summary_4v <- bind_rows(
  mono_4v |> filter(!is.na(mono_ptau217)) |>
    count(class = mono_ptau217) |> mutate(label = "pTau217"),
  mono_4v |> filter(!is.na(mono_nfl)) |>
    count(class = mono_nfl) |> mutate(label = "NfL"),
  mono_4v |> filter(!is.na(mono_lucent)) |>
    count(class = mono_lucent) |> mutate(label = "Lucent AD")
) |>
  group_by(label) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup()

message("\nMonotonicity (4-visit subgroup):")
print(mono_summary_4v)

# -- Figure 1: CV by interval group (violin + jitter) -------------------------

fig_cv_interval <- cv_df |>
  filter(!is.na(interval_group)) |>
  mutate(
    label          = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)),
    interval_group = factor(interval_group, levels = names(interval_colors))
  ) |>
  ggplot(aes(x = interval_group, y = cv_pct,
             fill = interval_group, color = interval_group)) +
  geom_violin(alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.6) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(
    data = kw_interval |>
      mutate(label = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label))),
    aes(x = 2, y = Inf, label = sprintf("KW p=%.3f", kw_p)),
    inherit.aes = FALSE, vjust = 1.5, size = 3.2, color = "grey30"
  ) +
  scale_fill_manual(values  = interval_colors, guide = "none") +
  scale_color_manual(values = interval_colors, guide = "none") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    title    = "Within-Patient CV by V1->V2 Interval (Variability Sub-study, N=40)",
    subtitle = "Dashed: 20% TEa. Does shorter interval compress the noise floor?",
    x        = "V1->V2 interval",
    y        = "Within-patient CV (%)"
  ) +
  theme_research() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

fsz <- get_figure_size("violin_grouped")
ggsave(file.path(out_dir, "fig_cv_by_interval.pdf"),
       fig_cv_interval, width = fsz["width"] * 1.4, height = fsz["height"])
ggsave(file.path(out_dir, "fig_cv_by_interval.png"),
       fig_cv_interval, width = fsz["width"] * 1.4, height = fsz["height"], dpi = 300)
message("Saved: fig_cv_by_interval")

# -- Figure 2: RCI distribution by interval group ----------------------------

rci_colors <- c(
  "Reliable increase"  = "#E41A1C",
  "No reliable change" = "#999999",
  "Reliable decrease"  = "#377EB8"
)

fig_rci_interval <- rci_df |>
  filter(!is.na(interval_group)) |>
  mutate(
    label          = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)),
    interval_group = factor(interval_group, levels = names(interval_colors)),
    rci_class      = factor(rci_class, levels = names(rci_colors))
  ) |>
  ggplot(aes(x = rci, fill = rci_class, color = rci_class)) +
  geom_vline(xintercept = c(-RCI_THRESHOLD, RCI_THRESHOLD),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_histogram(bins = 10, alpha = 0.7, position = "identity") +
  scale_fill_manual(values  = rci_colors, name = NULL) +
  scale_color_manual(values = rci_colors, name = NULL) +
  facet_grid(interval_group ~ label, scales = "free_x") +
  labs(
    title    = "RCI Distribution by Interval Group (Variability Sub-study)",
    subtitle = "Rows = interval, columns = biomarker; dashed: |RCI| = 1.96",
    x        = "RCI (V1->V2)",
    y        = "N patients"
  ) +
  theme_research(base_size = 11) +
  theme(legend.position = "bottom")

fsz_g <- get_figure_size("heatmap")
ggsave(file.path(out_dir, "fig_rci_by_interval.pdf"),
       fig_rci_interval, width = fsz_g["width"] * 1.2, height = fsz_g["height"])
ggsave(file.path(out_dir, "fig_rci_by_interval.png"),
       fig_rci_interval, width = fsz_g["width"] * 1.2, height = fsz_g["height"], dpi = 300)
message("Saved: fig_rci_by_interval")

# -- Figure 3: 4-visit spaghetti by biomarker (colored by interval group) -----

spag4_df <- sub4 |>
  pivot_longer(cols = names(biomarkers), names_to = "biomarker", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    label          = map_chr(biomarker, ~ biomarkers[[.x]]$label),
    label          = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)),
    interval_group = factor(interval_group, levels = names(interval_colors))
  )

# Group median per visit
group_median <- spag4_df |>
  group_by(label, visit) |>
  summarise(median_val = median(value, na.rm = TRUE), .groups = "drop")

fig_spag4 <- ggplot() +
  geom_line(
    data    = spag4_df,
    aes(x = visit, y = value, group = record_id, color = interval_group),
    alpha   = 0.35, linewidth = 0.5
  ) +
  geom_point(
    data  = spag4_df,
    aes(x = visit, y = value, color = interval_group),
    alpha = 0.5, size = 1.2
  ) +
  geom_line(
    data      = group_median,
    aes(x = visit, y = median_val, group = 1),
    color     = "black", linewidth = 1.3, linetype = "solid"
  ) +
  scale_color_manual(values = interval_colors, name = "V1->V2 interval") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    title    = "4-Visit Trajectories: Variability Sub-study",
    subtitle = "N=23 patients with complete 4-visit data; black = cohort median; color = V1->V2 interval",
    x        = "Visit",
    y        = "Biomarker value"
  ) +
  theme_research()

fsz_s <- get_figure_size("scatter_multi")
ggsave(file.path(out_dir, "fig_4visit_trajectories.pdf"),
       fig_spag4, width = fsz_s["width"] * 0.85, height = fsz_s["height"] * 0.5)
ggsave(file.path(out_dir, "fig_4visit_trajectories.png"),
       fig_spag4, width = fsz_s["width"] * 0.85, height = fsz_s["height"] * 0.5, dpi = 300)
message("Saved: fig_4visit_trajectories")

# -- Figure 4: Monotonicity stacked bar (4-visit) -----------------------------

mono_colors <- c(
  "Monotone up"   = "#E41A1C",
  "Non-monotone"  = "#999999",
  "Monotone down" = "#377EB8"
)

fig_mono4 <- mono_summary_4v |>
  mutate(
    class = factor(class, levels = names(mono_colors)),
    label = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label))
  ) |>
  ggplot(aes(x = label, y = pct, fill = class)) +
  geom_col(position = "stack", width = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = mono_colors, name = "Trajectory") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Trajectory Monotonicity: 4-Visit Sub-study Patients",
    subtitle = "N=23 patients with complete 4-visit data",
    x        = NULL,
    y        = "% patients"
  ) +
  theme_research()

fsz_b <- get_figure_size("bar")
ggsave(file.path(out_dir, "fig_4visit_monotonicity.pdf"),
       fig_mono4, width = fsz_b["width"] * 0.7, height = fsz_b["height"])
ggsave(file.path(out_dir, "fig_4visit_monotonicity.png"),
       fig_mono4, width = fsz_b["width"] * 0.7, height = fsz_b["height"], dpi = 300)
message("Saved: fig_4visit_monotonicity")

# -- Figure 5: CV vs follow-up duration (continuous, all sub-study) -----------

cv_follow <- cv_df |>
  left_join(
    sub |> group_by(record_id) |>
      summarise(max_days = max(days_from_v1), .groups = "drop"),
    by = "record_id"
  ) |>
  filter(!is.na(interval_group)) |>
  mutate(
    label          = factor(label, levels = map_chr(names(biomarkers), ~ biomarkers[[.x]]$label)),
    interval_group = factor(interval_group, levels = names(interval_colors))
  )

fig_cv_follow <- ggplot(cv_follow,
                        aes(x = max_days, y = cv_pct, color = interval_group)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "black", linewidth = 0.8, alpha = 0.12) +
  scale_color_manual(values = interval_colors, name = "V1->V2 interval") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    title    = "CV vs. Total Follow-up Duration",
    subtitle = "Does longer total follow-up reveal more variability?",
    x        = "Total follow-up (days from V1 to last visit)",
    y        = "Within-patient CV (%)"
  ) +
  theme_research()

ggsave(file.path(out_dir, "fig_cv_vs_followup.pdf"),
       fig_cv_follow, width = fsz_s["width"] * 0.85, height = fsz_s["height"] * 0.5)
ggsave(file.path(out_dir, "fig_cv_vs_followup.png"),
       fig_cv_follow, width = fsz_s["width"] * 0.85, height = fsz_s["height"] * 0.5,
       dpi = 300)
message("Saved: fig_cv_vs_followup")

# -- Summary table -------------------------------------------------------------

gt_cv_int <- cv_summary |>
  mutate(
    cv_fmt       = sprintf("%.1f (%.1f-%.1f)", median_cv, q25, q75),
    above_fmt    = sprintf("%.0f%%", pct_above_20),
    interval_group = factor(interval_group, levels = names(interval_colors))
  ) |>
  select(Biomarker = label, `Interval group` = interval_group,
         N = n, `Median CV% (IQR)` = cv_fmt, `% >20% TEa` = above_fmt) |>
  arrange(Biomarker, `Interval group`) |>
  gt(groupname_col = "Biomarker") |>
  tab_header(
    title    = "Within-Patient CV by V1->V2 Interval: Variability Sub-study",
    subtitle = "N=40 consented patients; all-visit CV pooled across available visits"
  ) |>
  tab_footnote(
    footnote = "CV computed across all available visits (>=2). Interval group assigned by V1->V2 gap. KW p-values: pTau217, NfL, Lucent AD tested separately."
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_cv_int, file.path(out_dir, "table_cv_by_interval.html"))
gtsave(gt_cv_int, file.path(out_dir, "table_cv_by_interval.png"), expand = 20)
message("Saved: table_cv_by_interval")

# -- Write CSVs ----------------------------------------------------------------

write_csv(cv_df,            file.path(out_dir, "cv_varstudy.csv"))
write_csv(cv_summary,       file.path(out_dir, "cv_summary_by_interval.csv"))
write_csv(rci_df,           file.path(out_dir, "rci_varstudy.csv"))
write_csv(rci_summary,      file.path(out_dir, "rci_summary_by_interval.csv"))
write_csv(mono_summary_4v,  file.path(out_dir, "monotonicity_4visit.csv"))

message("\nDone. Outputs written to: ", out_dir)
