#!/usr/bin/env Rscript
# dac_narrative_figures.R -- Two composite narrative figures
#
# Figure 1 (4 panels): "The noise floor is fixed"
#   A. pTau217 CV by V1->V2 interval -- shortening intervals doesn't compress noise
#   B. ICC vs CV per biomarker -- stable ranking != stable values
#   C. Biological CV vs 20% preanalytical TEa -- floor sits at the boundary
#   D. SD_diff stability across V1->V2, V1->V3, V1->V4 -- more time doesn't help
#
# Figure 2 (2 panels): "Signal within the noise"
#   A. NfL monotonicity asymmetry (4-visit subgroup) -- directional decline in subset
#   B. MoCA vs Lucent AD within-patient CV -- cognitive impairment predicts instability
#
# Usage:
#   Rscript analysis/dac_narrative_figures.R [--data-dir PATH] [--results-dir PATH] [--output-dir PATH]
#
# Requires: tidyverse, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(purrr)
  library(forcats)
  library(patchwork)
})

# -- Paths ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) normalizePath(file.path(getwd(), "analysis"), mustWork = FALSE)
)
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
out_dir     <- file.path(project_dir, "results", "dac_narrative")

for (i in seq_along(args)) {
  if (args[i] == "--data-dir"    && i < length(args)) data_dir    <- args[i + 1]
  if (args[i] == "--results-dir" && i < length(args)) results_dir <- args[i + 1]
  if (args[i] == "--output-dir"  && i < length(args)) out_dir     <- args[i + 1]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "theme_research.R"))
source(file.path(script_dir, "palettes.R"))

message("Writing narrative figures to: ", out_dir)

# -- Load pre-computed results -------------------------------------------------

cv_summary    <- read_csv(file.path(results_dir, "dac_icc_cv",    "cv_summary.csv"),           show_col_types = FALSE)
icc_results   <- read_csv(file.path(results_dir, "dac_icc_cv",    "icc_lmm_results.csv"),       show_col_types = FALSE)
cv_per_pt     <- read_csv(file.path(results_dir, "dac_icc_cv",    "cv_per_patient.csv"),         show_col_types = FALSE)
rci_pairs     <- read_csv(file.path(results_dir, "dac_rci",       "rci_summary_all_pairs.csv"),  show_col_types = FALSE)
cv_varstudy   <- read_csv(file.path(results_dir, "dac_varstudy",  "cv_varstudy.csv"),            show_col_types = FALSE)
mono_4v       <- read_csv(file.path(results_dir, "dac_varstudy",  "monotonicity_4visit.csv"),    show_col_types = FALSE)
cor_moca_cv   <- read_csv(file.path(results_dir, "dac_clinical_anchor", "cor_moca_vs_cv.csv"),   show_col_types = FALSE)

# Raw data needed for panels A and B (Fig. 2)
long     <- read_csv(file.path(data_dir, "dac_long.csv"),     show_col_types = FALSE)
baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE)

var_ids <- baseline |> filter(variability_consented == TRUE) |> pull(record_id)

long_ids <- long |>
  filter(!is_redraw) |>
  group_by(record_id) |> filter(n_distinct(visit) >= 2) |>
  pull(record_id) |> unique()

lucent_cv <- cv_per_pt |>
  filter(biomarker == "lucent_ad") |>
  select(record_id, cv_pct)

moca_cv <- baseline |>
  filter(record_id %in% long_ids) |>
  select(record_id, moca_raw) |>
  left_join(lucent_cv, by = "record_id") |>
  filter(!is.na(moca_raw), !is.na(cv_pct))

# -- Shared theme layer --------------------------------------------------------

THEME <- theme_research(base_size = 12) +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    strip.text       = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.key.size  = unit(0.4, "cm"),
    plot.tag         = element_text(size = 13, face = "bold")
  )

TEA <- 20   # Mayfield preanalytical TEa (%)

# Elahi lab palette (from palettes.yaml)
elahi <- lab_palette("elahi")

# Reserve elahi[1]=red, elahi[2]=blue for M/F comparisons.
#
# Core trio: purple/orange/green used for interval groups and the 3 main biomarkers.
# GFAP/Ab42/40 (panel B only) get pink and grey as secondary slots.

# Interval group colors (panel A)
interval_colors <- c(
  "Short\n(<=30d)"   = elahi[4],  # purple
  "Medium\n(31-90d)" = elahi[5],  # orange
  "Long\n(>90d)"     = elahi[3]   # green
)

# Biomarker colors (panels B/C/D and Fig. 2 scatter)
bm_colors <- c(
  "pTau217"   = elahi[4],         # purple
  "NfL"       = elahi[5],         # orange
  "GFAP"      = elahi[7],         # pink
  "Ab42/40"   = elahi[8],         # grey
  "Lucent AD" = elahi[3]          # green
)

mono_colors <- c(
  "Monotone up"   = elahi[5],     # orange
  "Non-monotone"  = elahi[8],     # grey
  "Monotone down" = elahi[4]      # purple
)

bm_order <- c("pTau217", "NfL", "GFAP", "Ab42/40", "Lucent AD")

# =============================================================================
# FIGURE 1 — THE NOISE FLOOR IS FIXED
# =============================================================================

# -- Panel A: |% change| vs interval days (all consecutive pairs) -------------
# Use every Vn->Vn+1 pair, not just V1->V2, to test interval independence.

# Build consecutive-pair dataset from variability sub-study patients
consec_pairs <- local({
  sub_raw <- long |>
    filter(record_id %in% var_ids, !is_redraw) |>
    arrange(record_id, days_from_v1)

  map_dfr(split(sub_raw, sub_raw$record_id), function(grp) {
    grp <- grp[order(grp$days_from_v1), ]
    if (nrow(grp) < 2) return(NULL)
    map_dfr(seq_len(nrow(grp) - 1), function(i) {
      r1 <- grp[i, ]; r2 <- grp[i + 1, ]
      interval <- r2$days_from_v1 - r1$days_from_v1
      pair_lbl <- paste0(as.character(r1$visit), "->", as.character(r2$visit))
      map_dfr(c("ptau217", "nfl", "lucent_ad"), function(bm) {
        v1_val <- r1[[bm]]; v2_val <- r2[[bm]]
        if (is.na(v1_val) || is.na(v2_val) || v1_val <= 0) return(NULL)
        tibble(
          record_id      = r1$record_id,
          pair           = pair_lbl,
          interval_days  = interval,
          biomarker      = bm,
          abs_pct_change = abs(v2_val - v1_val) / v1_val * 100
        )
      })
    })
  }) |>
    mutate(
      label = c(ptau217 = "pTau217", nfl = "NfL", lucent_ad = "Lucent AD")[biomarker],
      label = factor(label, levels = c("pTau217", "NfL", "Lucent AD")),
      pair  = factor(pair, levels = c("V1->V2", "V2->V3", "V3->V4"))
    )
})

# Correlation of |% change| with interval (Spearman, per biomarker)
cor_interval <- consec_pairs |>
  group_by(label) |>
  summarise(
    rho = cor(interval_days, abs_pct_change, method = "spearman"),
    p   = cor.test(interval_days, abs_pct_change,
                   method = "spearman", exact = FALSE)$p.value,
    n   = n(),
    .groups = "drop"
  )
message("\nSpearman rho (|%change| vs interval):")
print(cor_interval)

consec_binned <- consec_pairs |>
  mutate(
    interval_group = cut(
      interval_days,
      breaks = c(0, 30, 90, Inf),
      labels = c("Short\n(<=30d)", "Medium\n(31-90d)", "Long\n(>90d)"),
      right  = TRUE
    )
  ) |>
  filter(!is.na(interval_group))

# KW p per biomarker across interval bins
kw_consec <- consec_binned |>
  group_by(label) |>
  summarise(
    kw_p = tryCatch(
      kruskal.test(abs_pct_change ~ interval_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

pa <- ggplot(consec_binned,
             aes(x = interval_group, y = abs_pct_change,
                 fill = interval_group, color = interval_group)) +
  geom_hline(yintercept = TEA, linetype = "dashed",
             color = "grey30", linewidth = 0.6) +
  geom_violin(alpha = 0.18, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.14, outlier.shape = NA,
               fill = "white", linewidth = 0.55) +
  geom_jitter(width = 0.09, size = 1.6, alpha = 0.65) +
  scale_fill_manual(values = interval_colors, guide = "none") +
  scale_color_manual(values = interval_colors, guide = "none") +
  geom_text(
    data = kw_consec,
    aes(x = 2, y = Inf,
        label = sprintf("KW p=%.2f", kw_p)),
    inherit.aes = FALSE, vjust = 1.5, size = 3, color = "grey30"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    tag      = "A",
    title    = "Noise floor is interval-independent",
    subtitle = "All consecutive Vn->Vn+1 pairs; dashed = 20% TEa",
    x        = "Inter-visit interval",
    y        = "|% change| from prior visit"
  ) +
  THEME

# -- Panel B: ICC vs median CV scatter ----------------------------------------

icc_cv_plot <- icc_results |>
  left_join(cv_summary |> select(label, median_cv, pct_above_20),
            by = "label") |>
  filter(!is.na(icc)) |>
  mutate(
    label = factor(label, levels = bm_order),
    # Quadrant annotation: top-right = bad (high CV + low ICC), top-left = paradox zone
    quadrant = case_when(
      icc >= 0.75 & median_cv >= 15 ~ "High reliability\nHigh noise",
      icc >= 0.75 & median_cv < 15  ~ "High reliability\nLow noise",
      TRUE                           ~ "Low reliability"
    )
  )

pb <- ggplot(icc_cv_plot, aes(x = median_cv, y = icc, color = label)) +
  # Quadrant shading
  annotate("rect", xmin = 15, xmax = Inf, ymin = 0.75, ymax = Inf,
           fill = "#fff3cd", alpha = 0.6) +
  annotate("text", x = 15.5, y = 0.77, label = "High reliability,\nhigh noise",
           hjust = 0, size = 2.9, color = "#856404", fontface = "italic") +
  geom_vline(xintercept = TEA, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "grey50", linewidth = 0.5) +
  geom_point(size = 4) +
  scale_color_manual(values = bm_colors, guide = "none") +
  ggrepel::geom_text_repel(aes(label = label), size = 3,
    box.padding = 0.4, segment.size = 0.3, segment.color = "grey60",
    color = "grey20") +
  scale_y_continuous(limits = c(0.75, 1.0), breaks = seq(0.75, 1.0, 0.05)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(
    tag      = "B",
    title    = "Reliable ranking, noisy values",
    subtitle = "High ICC does not mean low variability",
    x        = "Median within-patient CV (%)",
    y        = "ICC (LMM)"
  ) +
  THEME

# -- Panel C: Biological CV vs preanalytical TEa floor -------------------------

# Show biological CV distribution per biomarker with 20% TEa as reference.
# Annotate combined minimum detectable change (MDC) for a serial pair:
#   MDC = 1.96 * sqrt(2) * sqrt(TEa^2 + bio_CV^2)
cv_noise <- cv_summary |>
  filter(label %in% c("pTau217", "NfL", "Lucent AD")) |>
  mutate(
    label = factor(label, levels = bm_order),
    mdc_combined = 1.96 * sqrt(2) * sqrt((TEA / 100)^2 + (median_cv / 100)^2) * 100
  )

# Pull per-patient CV for these biomarkers for the distribution
cv_dist <- cv_per_pt |>
  filter(label %in% c("pTau217", "NfL", "Lucent AD")) |>
  mutate(label = factor(label, levels = bm_order))

# Legend labels carry the per-biomarker MDC value; avoids in-figure text clutter.
mdc_labels <- setNames(
  sprintf("%s  (MDC = %.0f%%)", as.character(cv_noise$label), cv_noise$mdc_combined),
  as.character(cv_noise$label)
)

pc <- ggplot() +
  geom_hline(yintercept = TEA, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_violin(data = cv_dist,
              aes(x = label, y = cv_pct, fill = label),
              alpha = 0.2, color = NA, trim = TRUE) +
  geom_boxplot(data = cv_dist,
               aes(x = label, y = cv_pct, color = label),
               width = 0.12, outlier.shape = NA,
               fill = "white", linewidth = 0.55) +
  # MDC threshold lines (labeled via legend, not in-figure text)
  geom_segment(data = cv_noise,
               aes(x = as.integer(label) - 0.35,
                   xend = as.integer(label) + 0.35,
                   y = mdc_combined, yend = mdc_combined,
                   color = label),
               linewidth = 1.1, linetype = "solid") +
  scale_fill_manual(values = bm_colors, labels = mdc_labels,
                    name = "Combined MDC", guide = "none") +
  scale_color_manual(values = bm_colors, labels = mdc_labels,
                     name = "Combined MDC",
                     guide = guide_legend(override.aes = list(linewidth = 1.5))) +
  annotate("text", x = 3.48, y = TEA + 1.2, label = "20% TEa (preanalytical)",
           size = 2.9, hjust = 1, color = "grey30") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
  labs(
    tag      = "C",
    title    = "Noise budget",
    subtitle = "Biological CV sits at the 20% preanalytical floor;\ncombined MDC = 1.96*sqrt(2*(TEa^2+bio_CV^2))",
    x        = NULL,
    y        = "Within-patient CV (%)"
  ) +
  THEME

# -- Panel D: SD_diff stability across visit pairs ----------------------------
# Normalize SD_diff to % of V1 cohort mean so pTau217 (raw ~0.02) and NfL
# (raw ~9) are on the same relative scale.

v1_means <- long |>
  filter(visit == "V1", !is_redraw) |>
  summarise(
    pTau217   = mean(ptau217,   na.rm = TRUE),
    NfL       = mean(nfl,       na.rm = TRUE),
    `Lucent AD` = mean(lucent_ad, na.rm = TRUE)
  ) |>
  pivot_longer(everything(), names_to = "label", values_to = "v1_mean")

sd_trend <- rci_pairs |>
  filter(label %in% c("pTau217", "NfL", "Lucent AD"),
         n_pair_bm >= 5) |>
  left_join(v1_means, by = "label") |>
  mutate(
    sd_diff_pct = sd_diff / v1_mean * 100,
    label       = factor(label, levels = bm_order),
    visit_pair  = factor(visit_pair, levels = c("V1->V2", "V1->V3", "V1->V4"))
  )

pd <- ggplot(sd_trend,
             aes(x = visit_pair, y = sd_diff_pct,
                 group = label, color = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n_pair_bm)) +
  scale_color_manual(values = bm_colors, name = "Biomarker") +
  scale_size_continuous(range = c(2, 5), name = "N pairs") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(
    tag      = "D",
    title    = "Noise floor doesn't shrink with time",
    subtitle = "SD_diff (% of V1 mean) stable from V2 to V4",
    x        = "Visit comparison",
    y        = "SD_diff (% of V1 cohort mean)"
  ) +
  THEME +
  theme(legend.position = "right")

# -- Assemble Figure 1 ---------------------------------------------------------

fig1 <- pa / (pb | pc | pd) +
  plot_annotation(
    title   = "Figure 1. The biological noise floor is fixed",
    caption = "TEa = total allowable error (Mayfield et al.); MDC = minimum detectable change; ICC = intraclass correlation coefficient (LMM); KW = Kruskal-Wallis.",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, color = "grey50")
    )
  )

ggsave(file.path(out_dir, "fig1_noise_floor_fixed.pdf"),
       fig1, width = 14, height = 12)
ggsave(file.path(out_dir, "fig1_noise_floor_fixed.png"),
       fig1, width = 14, height = 12, dpi = 300)
message("Saved: fig1_noise_floor_fixed")

# =============================================================================
# FIGURE 2 — SIGNAL WITHIN THE NOISE
# =============================================================================

# -- Panel A: NfL monotonicity asymmetry (4-visit) ----------------------------

nfl_mono <- mono_4v |>
  filter(label == "NfL") |>
  mutate(class = factor(class, levels = names(mono_colors)))

ptau_mono <- mono_4v |>
  filter(label == "pTau217") |>
  mutate(class = factor(class, levels = names(mono_colors)))

mono_plot_df <- bind_rows(nfl_mono, ptau_mono) |>
  mutate(label = factor(label, levels = c("pTau217", "NfL")))

pe <- ggplot(mono_plot_df,
             aes(x = label, y = pct, fill = class)) +
  geom_col(position = "stack", width = 0.45) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = mono_colors, name = "Trajectory") +
  # Arrow annotation for NfL asymmetry; elahi[2] = blue (reserved for directional annotation)
  annotate("segment",
           x = 2.3, xend = 2.3, y = 60, yend = 90,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           color = elahi[2], linewidth = 0.8) +
  annotate("text", x = 2.45, y = 75,
           label = "Asymmetric\ndecline", color = elahi[2],
           size = 3, hjust = 0, fontface = "italic") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    tag      = "A",
    title    = "NfL shows directional decline in 1-in-5",
    subtitle = "4-visit patients (N=23); pTau217 shows symmetric noise",
    x        = NULL,
    y        = "% patients"
  ) +
  THEME

# -- Panel B: MoCA vs Lucent AD CV scatter ------------------------------------

rho_val <- cor_moca_cv |> filter(label == "Lucent AD") |> pull(rho)
p_val   <- cor_moca_cv |> filter(label == "Lucent AD") |> pull(p)

pf <- ggplot(moca_cv, aes(x = moca_raw, y = cv_pct)) +
  geom_hline(yintercept = TEA, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.75, color = bm_colors[["Lucent AD"]]) +
  geom_smooth(method = "lm", se = TRUE, color = bm_colors[["Lucent AD"]],
              fill = bm_colors[["Lucent AD"]], alpha = 0.12, linewidth = 0.9) +
  annotate("text", x = max(moca_cv$moca_raw, na.rm = TRUE),
           y = max(moca_cv$cv_pct, na.rm = TRUE) * 0.97,
           label = sprintf("rho = %.2f\np = %.3f\nn = %d", rho_val, p_val, nrow(moca_cv)),
           hjust = 1, vjust = 1, size = 3.4, color = "grey20") +
  annotate("text", x = max(moca_cv$moca_raw, na.rm = TRUE),
           y = TEA + 1.5, label = "20% TEa",
           hjust = 1, size = 3, color = "grey40") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    tag      = "B",
    title    = "Cognitive impairment predicts Lucent AD instability",
    subtitle = "Lower MoCA -> higher within-patient CV (Spearman rho = -0.53)",
    x        = "MoCA (baseline)",
    y        = "Lucent AD within-patient CV (%)"
  ) +
  THEME

# -- Assemble Figure 2 ---------------------------------------------------------

fig2 <- (pe | pf) +
  plot_annotation(
    title   = "Figure 2. Signal within the noise",
    caption = "4-visit subgroup: N=23 patients with complete V1-V4 data. MoCA analysis: N=29 longitudinal patients with baseline MoCA.",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, color = "grey50")
    )
  )

ggsave(file.path(out_dir, "fig2_signal_within_noise.pdf"),
       fig2, width = 14, height = 6)
ggsave(file.path(out_dir, "fig2_signal_within_noise.png"),
       fig2, width = 14, height = 6, dpi = 300)
message("Saved: fig2_signal_within_noise")

# =============================================================================
# FIGURE 3 — MONOTONE PATIENT TRAJECTORIES
# =============================================================================
# mono_4v is a summary table; re-derive patient-level NfL class from long data.

pts_4v <- long |>
  filter(!is_redraw) |>
  group_by(record_id) |>
  filter(n_distinct(visit) == 4) |>
  pull(record_id) |>
  unique()

# Relaxed monotonicity: allow 1 step in the wrong direction provided the
# V1->V4 overall direction is conserved.  Strict (0 violations) gave N=2
# monotone cases; relaxed gives N=14, making stratification tractable.
classify_mono <- function(x) {
  d            <- diff(x)
  overall_up   <- x[length(x)] > x[1]
  overall_down <- x[length(x)] < x[1]
  if      (sum(d < 0) <= 1 && overall_up)   "Monotone up"
  else if (sum(d > 0) <= 1 && overall_down) "Monotone down"
  else "Non-monotone"
}

nfl_class <- long |>
  filter(record_id %in% pts_4v, !is_redraw, !is.na(nfl)) |>
  arrange(record_id, visit) |>
  group_by(record_id) |>
  filter(n() == 4) |>
  summarise(class = classify_mono(nfl), .groups = "drop") |>
  mutate(class = factor(class, levels = names(mono_colors)))

# Longitudinal biomarkers for these patients, normalized to % change from V1
traj <- long |>
  filter(record_id %in% nfl_class$record_id, !is_redraw) |>
  select(record_id, visit, ptau217, nfl, lucent_ad) |>
  pivot_longer(c(ptau217, nfl, lucent_ad),
               names_to = "biomarker", values_to = "value") |>
  mutate(
    label = c(ptau217 = "pTau217", nfl = "NfL", lucent_ad = "Lucent AD")[biomarker],
    label = factor(label, levels = c("pTau217", "NfL", "Lucent AD")),
    visit = factor(visit, levels = c("V1", "V2", "V3", "V4"))
  ) |>
  filter(!is.na(value)) |>
  left_join(nfl_class, by = "record_id") |>
  group_by(record_id, label) |>
  mutate(
    v1_val     = value[visit == "V1"][1],
    pct_change = (value - v1_val) / v1_val * 100
  ) |>
  filter(!is.na(v1_val)) |>
  ungroup()

# Baseline MoCA by NfL class
moca_traj <- baseline |>
  filter(record_id %in% nfl_class$record_id) |>
  select(record_id, moca_raw) |>
  left_join(nfl_class, by = "record_id") |>
  filter(!is.na(moca_raw), !is.na(class))

# -- Panel A: spaghetti trajectories, faceted by biomarker --------------------

p3a <- ggplot(traj,
              aes(x = visit, y = pct_change,
                  group = record_id, color = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_line(alpha = 0.55, linewidth = 0.8) +
  geom_point(size = 1.8, alpha = 0.75) +
  stat_summary(aes(group = class),
               fun = median, geom = "line",
               linewidth = 1.6, linetype = "solid", alpha = 0.9) +
  scale_color_manual(values = mono_colors, name = "NfL class") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    tag      = "A",
    title    = "Biomarker trajectories by NfL monotonicity class",
    subtitle = "Thin lines = individual patients; thick = class median; % change from V1",
    x        = "Visit",
    y        = "% change from V1"
  ) +
  THEME +
  theme(legend.position = "right")

# -- Panel B: MoCA at baseline by NfL class -----------------------------------

p3b <- ggplot(moca_traj, aes(x = class, y = moca_raw, color = class)) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "white", linewidth = 0.6) +
  geom_jitter(width = 0.1, size = 2.2, alpha = 0.7) +
  scale_color_manual(values = mono_colors, guide = "none") +
  labs(
    tag      = "B",
    title    = "Baseline MoCA by NfL class",
    subtitle = sprintf("N = %d patients with complete V1-V4 NfL", nrow(nfl_class)),
    x        = NULL,
    y        = "MoCA (baseline)"
  ) +
  THEME

# -- Assemble Figure 3 ---------------------------------------------------------

fig3 <- p3a / p3b +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title   = "Figure 3. Monotone patient trajectories",
    caption = "NfL class derived from V1-V4 consecutive differences. % change = (value - V1) / V1 * 100.",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, color = "grey50")
    )
  )

ggsave(file.path(out_dir, "fig3_monotone_trajectories.pdf"),
       fig3, width = 14, height = 10)
ggsave(file.path(out_dir, "fig3_monotone_trajectories.png"),
       fig3, width = 14, height = 10, dpi = 300)
message("Saved: fig3_monotone_trajectories")

# =============================================================================
# FIGURE 4 — STRATIFIED MONOTONE ANALYSIS
# =============================================================================
# Restrict to patients with V1 Lucent AD for consistent N across both stratifiers.

v1_strat <- long |>
  filter(record_id %in% nfl_class$record_id, visit == "V1", !is_redraw) |>
  select(record_id, nfl_v1 = nfl, lucent_ad_v1 = lucent_ad) |>
  filter(!is.na(lucent_ad_v1))   # n = 19

nfl_mean_cut    <- mean(v1_strat$nfl_v1,       na.rm = TRUE)
lucent_med_cut  <- median(v1_strat$lucent_ad_v1, na.rm = TRUE)

strat <- v1_strat |>
  mutate(
    nfl_group    = factor(if_else(nfl_v1       <= nfl_mean_cut,   "Low", "High"),
                          levels = c("Low", "High")),
    lucent_group = factor(if_else(lucent_ad_v1 <= lucent_med_cut, "Low", "High"),
                          levels = c("Low", "High"))
  ) |>
  left_join(nfl_class, by = "record_id")

message(sprintf("Strat N=%d | NfL cut=%.1f | Lucent AD cut=%.0f",
                nrow(strat), nfl_mean_cut, lucent_med_cut))
message("NfL group x class:")
print(table(strat$nfl_group, strat$class))
message("Lucent AD group x class:")
print(table(strat$lucent_group, strat$class))

# Trajectory data restricted to same 19 patients
traj_strat <- traj |>
  filter(record_id %in% strat$record_id) |>
  left_join(strat |> select(record_id, nfl_group, lucent_group),
            by = "record_id")

# -- Panel A: spaghetti, color = class, linetype = NfL group ------------------

p4a <- ggplot(traj_strat,
              aes(x = visit, y = pct_change,
                  group = record_id, color = class, linetype = nfl_group)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.3) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_color_manual(values = mono_colors, name = "NfL class") +
  scale_linetype_manual(values = c("Low" = "dashed", "High" = "solid"),
                        name = "Baseline NfL") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(
    tag      = "A",
    title    = "Trajectories by class and baseline NfL level",
    subtitle = sprintf(
      "N=%d (Lucent AD-complete); solid = High NfL, dashed = Low NfL (cut = %.1f pg/mL)",
      nrow(strat), nfl_mean_cut),
    x = "Visit", y = "% change from V1"
  ) +
  THEME +
  theme(legend.position = "right")

# -- Panels B & C: stacked bars with Fisher's exact ---------------------------

make_class_bar <- function(df, group_var, group_label, cut_note, tag_lbl, title_lbl) {
  counts <- df |>
    count(group = .data[[group_var]], class, .drop = FALSE) |>
    group_by(group) |>
    mutate(pct = n / sum(n) * 100) |>
    ungroup()

  ct <- table(df[[group_var]], df$class)
  fp <- fisher.test(ct, simulate.p.value = TRUE, B = 9999)$p.value

  ggplot(counts, aes(x = group, y = pct, fill = class)) +
    geom_col(position = "stack", width = 0.5) +
    geom_text(
      aes(label = if_else(n > 0, sprintf("%d\n(%.0f%%)", n, pct), "")),
      position = position_stack(vjust = 0.5),
      size = 3, color = "white", fontface = "bold"
    ) +
    annotate("text", x = 1.5, y = 107,
             label = sprintf("Fisher p = %.3f", fp),
             size = 3, color = "grey30") +
    scale_fill_manual(values = mono_colors, name = "NfL class") +
    scale_y_continuous(limits = c(0, 115),
                       labels = function(x) paste0(x, "%")) +
    labs(tag = tag_lbl, title = title_lbl,
         subtitle = cut_note, x = group_label, y = "% patients") +
    THEME +
    theme(legend.position = "none")
}

p4b <- make_class_bar(
  strat, "nfl_group", "Baseline NfL",
  sprintf("Mean cut = %.1f pg/mL", nfl_mean_cut),
  "B", "NfL class by baseline NfL"
)

p4c <- make_class_bar(
  strat, "lucent_group", "Baseline Lucent AD",
  sprintf("Median cut = %.0f  (no validated threshold)", lucent_med_cut),
  "C", "NfL class by Lucent AD level"
)

# -- Assemble Figure 4 --------------------------------------------------------

fig4 <- p4a / (p4b | p4c) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title   = "Figure 4. Stratified monotone analysis",
    caption = sprintf(
      paste0("NfL Low/High: mean split (%.1f pg/mL, N=%d). ",
             "Lucent AD Low/High: median split (%.0f, no validated threshold). ",
             "Fisher's exact p (B=9999). Exploratory."),
      nfl_mean_cut, nrow(strat), lucent_med_cut
    ),
    theme = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, color = "grey50")
    )
  )

ggsave(file.path(out_dir, "fig4_stratified_monotone.pdf"),
       fig4, width = 14, height = 10)
ggsave(file.path(out_dir, "fig4_stratified_monotone.png"),
       fig4, width = 14, height = 10, dpi = 300)
message("Saved: fig4_stratified_monotone")

message("\nDone. Outputs written to: ", out_dir)
