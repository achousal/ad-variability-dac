#!/usr/bin/env Rscript
# dac_lcmm.R -- Latent class trajectory analysis for DAC CRF biomarker data
#
# Fits 1-3 latent class linear trajectory models per biomarker using lcmm::hlme.
# Selects best model by BIC. Characterizes classes by trajectory shape and
# baseline clinical features.
#
# Usage:
#   Rscript analysis/dac_lcmm.R [--data-dir PATH] [--output-dir PATH]
#
# Requires: tidyverse, lcmm, gt, patchwork

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(forcats)
  library(patchwork)
  library(lcmm)
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
out_dir  <- file.path(project_dir, "results", "dac_lcmm")

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

v1_amyloid <- long |>
  filter(visit == "V1", !is_redraw) |>
  select(record_id, ptau217_risk, lucent_interpretation) |>
  distinct(record_id, .keep_all = TRUE)

baseline <- read_csv(file.path(data_dir, "dac_baseline.csv"), show_col_types = FALSE)

long_sub <- long |>
  filter(!is_redraw, !is.na(days_from_v1)) |>
  group_by(record_id) |>
  filter(n_distinct(visit) >= 2) |>
  ungroup() |>
  mutate(years_from_v1 = days_from_v1 / 365.25)

message(sprintf("Longitudinal subcohort: %d patients", n_distinct(long_sub$record_id)))

# -- Biomarkers to model -------------------------------------------------------

# pTau217 and NfL: log-transformed (right-skewed). Lucent AD: raw.
biomarkers <- list(
  ptau217   = list(label = "pTau217",   unit = "pg/mL", log_t = TRUE),
  nfl       = list(label = "NfL",       unit = "pg/mL", log_t = TRUE),
  lucent_ad = list(label = "Lucent AD", unit = "score",  log_t = FALSE)
)

MAX_CLASSES <- 3
N_RANDOM_STARTS <- 30   # random restarts to avoid local optima
set.seed(42)

# Build initial B vector for ng-class hlme (random intercepts only, shared variance)
# Layout: (ng-1) class probs | ng*(intercept+slope) | rand_int_var | resid_stderr
make_b_init <- function(m1_fit, ng, jitter_sd = 0.3) {
  fe <- m1_fit$best[1:2]             # population intercept, slope
  se <- abs(tail(m1_fit$best, 1))    # residual std error
  ri <- abs(m1_fit$best[5])          # approx random intercept var
  class_probs <- rep(-0.5, ng - 1)
  fixed_by_class <- as.vector(
    t(matrix(fe, nrow = 1)[rep(1, ng), ]) +
    rnorm(ng * 2, 0, jitter_sd * abs(fe + 0.01))
  )
  c(class_probs, fixed_by_class, max(ri, 1e-4), max(se, 1e-4))
}

# -- Helper: fit 1-k class models and return selection table ------------------

fit_lcmm_models <- function(dat, bm, max_ng = MAX_CLASSES) {
  meta <- biomarkers[[bm]]

  mod_dat <- dat |>
    select(record_id, years_from_v1, value = all_of(bm)) |>
    filter(!is.na(value)) |>
    mutate(
      value_fit = if (meta$log_t) log(value) else value,
      record_id = as.integer(record_id)
    ) |>
    arrange(record_id, years_from_v1)

  n_pts <- n_distinct(mod_dat$record_id)
  n_obs <- nrow(mod_dat)
  message(sprintf("  %s: %d patients, %d obs", meta$label, n_pts, n_obs))

  if (n_pts < 10) {
    message("  Skipping: N < 10")
    return(NULL)
  }

  # Convert to data.frame (lcmm doesn't handle tibbles)
  mod_dat <- as.data.frame(mod_dat)

  # -- Fit ng=1 as reference (needed for B= warm start in ng>1) --
  m1 <- tryCatch(
    hlme(value_fit ~ years_from_v1,
         random    = ~ years_from_v1,
         subject   = "record_id",
         ng        = 1,
         data      = mod_dat),
    error = function(e) {
      message("  ng=1 failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(m1)) return(NULL)

  models <- list("1" = m1)

  # -- Fit ng=2 and ng=3 with multiple random starts ------------------------
  for (ng in 2:max_ng) {
    if (n_pts < ng * 5) {
      message(sprintf("  Skipping ng=%d: N=%d too small (need %d)", ng, n_pts, ng * 5))
      next
    }

    message(sprintf("  Fitting ng=%d (%d random starts)...", ng, N_RANDOM_STARTS))

    # gridsearch() needs a model object already configured for ng classes.
    # Multi-start: fit hlme ng>1 with N_RANDOM_STARTS different B vectors,
    # keep the fit with lowest loglikelihood (best convergence).
    candidates <- list()
    for (s in seq_len(N_RANDOM_STARTS)) {
      b_try <- make_b_init(m1, ng, jitter_sd = 0.5)
      fit_s <- tryCatch(
        hlme(value_fit ~ years_from_v1,
             mixture  = ~ years_from_v1,
             random   = ~ 1,
             subject  = "record_id",
             ng       = ng,
             data     = mod_dat,
             B        = b_try,
             nwg      = FALSE,
             maxiter  = 100),
        error   = function(e) NULL,
        warning = function(w) suppressWarnings(
          hlme(value_fit ~ years_from_v1, mixture = ~years_from_v1, random = ~1,
               subject = "record_id", ng = ng, data = mod_dat,
               B = b_try, nwg = FALSE, maxiter = 100)
        )
      )
      if (!is.null(fit_s) && fit_s$conv == 1) candidates <- c(candidates, list(fit_s))
    }

    if (length(candidates) == 0) {
      message(sprintf("  ng=%d: no converged solutions", ng))
      next
    }

    # Select solution with lowest -2*loglik (best fit)
    best <- candidates[[which.min(map_dbl(candidates, ~ -.x$loglik))]]
    message(sprintf("  ng=%d: %d/%d converged, best BIC=%.1f",
                    ng, length(candidates), N_RANDOM_STARTS, best$BIC))
    models[[as.character(ng)]] <- best
  }

  # -- Model selection table --
  sel <- map_dfr(names(models), function(k) {
    m <- models[[k]]
    ng_val  <- as.integer(k)
    # Smallest class size
    if (ng_val > 1) {
      pp   <- m$pprob
      cls_n <- table(pp$class)
      min_class_n <- min(cls_n)
      min_class_pct <- min(cls_n) / sum(cls_n) * 100
    } else {
      min_class_n   <- n_pts
      min_class_pct <- 100
    }
    # Entropy (separation quality): ranges 0-1; higher = better separation
    entropy <- if (ng_val > 1) {
      pp <- m$pprob
      prob_cols <- grep("^prob", names(pp), value = TRUE)
      probs <- as.matrix(pp[, prob_cols])
      # Normalised entropy
      ent_raw <- -sum(probs * log(probs + 1e-10)) / nrow(probs)
      1 - ent_raw / log(ng_val)
    } else NA_real_

    tibble(
      ng            = ng_val,
      loglik        = m$loglik,
      AIC           = m$AIC,
      BIC           = m$BIC,
      n_patients    = n_pts,
      min_class_n   = min_class_n,
      min_class_pct = min_class_pct,
      entropy       = entropy
    )
  })

  list(
    models = models,
    sel    = sel,
    data   = mod_dat,
    meta   = meta,
    bm     = bm
  )
}

# -- Run for each biomarker ----------------------------------------------------

message("\nFitting latent class models (this may take ~1-2 min)...")
lcmm_results <- map(names(biomarkers), function(bm) {
  message(sprintf("\n=== %s ===", biomarkers[[bm]]$label))
  fit_lcmm_models(long_sub, bm)
}) |>
  set_names(names(biomarkers)) |>
  compact()

# -- Print selection tables ----------------------------------------------------

for (bm in names(lcmm_results)) {
  message(sprintf("\nModel selection: %s", biomarkers[[bm]]$label))
  print(lcmm_results[[bm]]$sel)
}

# -- Identify best model per biomarker (lowest BIC, min class >=5) ------------

get_best_model <- function(res) {
  sel <- res$sel |>
    filter(min_class_n >= 5) |>    # exclude solutions with tiny classes
    arrange(BIC)
  if (nrow(sel) == 0) sel <- res$sel |> arrange(BIC)
  best_ng <- sel$ng[1]
  list(model = res$models[[as.character(best_ng)]], ng = best_ng)
}

best_models <- map(lcmm_results, get_best_model)

for (bm in names(best_models)) {
  message(sprintf("%s: best model ng=%d", biomarkers[[bm]]$label, best_models[[bm]]$ng))
}

# -- Extract class assignments + posterior probabilities ----------------------

class_df <- map_dfr(names(best_models), function(bm) {
  bst  <- best_models[[bm]]
  mod  <- bst$model
  meta <- biomarkers[[bm]]
  ng   <- bst$ng

  if (ng == 1) {
    # Everyone in class 1
    ids <- unique(lcmm_results[[bm]]$data$record_id)
    return(tibble(
      record_id    = ids,
      class        = 1L,
      class_label  = "Class 1",
      max_prob     = 1,
      biomarker    = bm,
      label        = meta$label
    ))
  }

  pp <- mod$pprob |>
    as_tibble() |>
    rename(record_id = record_id, class = class) |>
    mutate(record_id = as.integer(record_id))

  # Maximum posterior probability (class certainty)
  prob_cols <- grep("^prob", names(pp), value = TRUE)
  pp <- pp |>
    mutate(max_prob = apply(as.matrix(pp[, prob_cols]), 1, max))

  pp |>
    select(record_id, class, max_prob) |>
    mutate(
      biomarker   = bm,
      label       = meta$label,
      class_label = paste0("Class ", class)
    )
})

message("\nClass assignments:")
print(class_df |> count(label, class_label))

# -- GT: model selection table -------------------------------------------------

sel_all <- map_dfr(names(lcmm_results), function(bm) {
  lcmm_results[[bm]]$sel |>
    mutate(
      biomarker = bm,
      label     = biomarkers[[bm]]$label,
      best      = ng == best_models[[bm]]$ng
    )
})

gt_sel <- sel_all |>
  mutate(
    bic_fmt     = sprintf("%.1f", BIC),
    aic_fmt     = sprintf("%.1f", AIC),
    ent_fmt     = if_else(!is.na(entropy), sprintf("%.2f", entropy), "—"),
    cls_fmt     = sprintf("%d (%.0f%%)", min_class_n, min_class_pct),
    best_marker = if_else(best, "*", "")
  ) |>
  select(
    Biomarker       = label,
    `N classes`     = ng,
    BIC             = bic_fmt,
    AIC             = aic_fmt,
    Entropy         = ent_fmt,
    `Min class (n, %)` = cls_fmt,
    `Best*`         = best_marker
  ) |>
  gt() |>
  tab_header(
    title    = "Latent Class Trajectory Model Selection",
    subtitle = "DAC CRF longitudinal subcohort -- pTau217, NfL, Lucent AD"
  ) |>
  tab_footnote(
    footnote = "Best model = lowest BIC with minimum class size >=5. Entropy (0-1): higher = better class separation. * = selected model."
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = `Best*`, rows = `Best*` == "*")
  ) |>
  opt_table_font(font = "Arial") |>
  opt_stylize(style = 1)

gtsave(gt_sel, file.path(out_dir, "table_model_selection.html"))
gtsave(gt_sel, file.path(out_dir, "table_model_selection.png"), expand = 20)
message("Saved: table_model_selection")

# -- Figure 1: Predicted class trajectories ------------------------------------

class_colors <- c(
  "Class 1" = "#377EB8",
  "Class 2" = "#E41A1C",
  "Class 3" = "#4DAF4A"
)

make_class_traj_plot <- function(bm) {
  bst  <- best_models[[bm]]
  mod  <- bst$model
  meta <- biomarkers[[bm]]
  ng   <- bst$ng
  dat  <- lcmm_results[[bm]]$data

  t_seq <- seq(0, max(dat$years_from_v1), length.out = 100)

  if (ng == 1) {
    # Single class: plot fitted line from M1
    fe     <- mod$best[grep("^intercept|^years", names(mod$best), ignore.case = TRUE)]
    # Use predictY for single class too
    nd <- data.frame(years_from_v1 = t_seq)
    pred <- tryCatch(
      predictY(mod, newdata = nd, var.time = "years_from_v1", draws = FALSE),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    pred_df <- tibble(
      years_from_v1 = t_seq,
      pred          = if (meta$log_t) exp(pred$pred[, 1]) else pred$pred[, 1],
      class_label   = "Class 1"
    )
  } else {
    nd <- data.frame(years_from_v1 = t_seq)
    pred <- tryCatch(
      predictY(mod, newdata = nd, var.time = "years_from_v1", draws = FALSE),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    pred_df <- map_dfr(seq_len(ng), function(k) {
      col_nm <- paste0("Ypred_class", k)
      if (!col_nm %in% colnames(pred$pred)) return(NULL)
      tibble(
        years_from_v1 = t_seq,
        pred          = if (meta$log_t) exp(pred$pred[, col_nm]) else pred$pred[, col_nm],
        class_label   = paste0("Class ", k)
      )
    })
  }

  # Class sizes for labels
  cls_sizes <- class_df |>
    filter(biomarker == bm) |>
    count(class_label) |>
    mutate(lbl = paste0(class_label, " (n=", n, ")"))
  pred_df <- pred_df |>
    left_join(cls_sizes |> select(class_label, lbl), by = "class_label")

  # Raw observations (spaghetti, faded)
  obs_df <- dat |>
    as_tibble() |>
    mutate(
      value_orig    = if (meta$log_t) exp(value_fit) else value_fit,
      record_id     = as.integer(record_id)
    ) |>
    left_join(class_df |> filter(biomarker == bm) |> select(record_id, class_label),
              by = "record_id")

  ggplot() +
    geom_line(
      data    = obs_df,
      aes(x = years_from_v1, y = value_orig, group = record_id, color = class_label),
      alpha   = 0.2, linewidth = 0.4
    ) +
    geom_line(
      data    = pred_df,
      aes(x = years_from_v1, y = pred, color = class_label, group = class_label),
      linewidth = 1.5
    ) +
    scale_color_manual(
      values = setNames(class_colors[seq_len(ng)], paste0("Class ", seq_len(ng))),
      labels = if (nrow(cls_sizes) > 0) setNames(cls_sizes$lbl, cls_sizes$class_label) else waiver(),
      name   = NULL
    ) +
    labs(
      title    = meta$label,
      subtitle = sprintf("ng=%d classes; bold = class mean trajectory", ng),
      x        = "Years from V1",
      y        = sprintf("%s (%s)", meta$label, meta$unit)
    ) +
    theme_research()
}

traj_plots <- compact(map(names(biomarkers), make_class_traj_plot))

if (length(traj_plots) > 0) {
  fig_traj <- wrap_plots(traj_plots, ncol = length(traj_plots)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  fsz_s <- get_figure_size("scatter_multi")
  ggsave(file.path(out_dir, "fig_class_trajectories.pdf"),
         fig_traj, width = fsz_s["width"], height = fsz_s["height"] * 0.55)
  ggsave(file.path(out_dir, "fig_class_trajectories.png"),
         fig_traj, width = fsz_s["width"], height = fsz_s["height"] * 0.55, dpi = 300)
  message("Saved: fig_class_trajectories")
}

# -- Figure 2: Class membership posterior probabilities (quality check) --------

multi_class_labels <- map_chr(
  names(best_models)[map_int(best_models, "ng") > 1],
  ~ biomarkers[[.x]]$label
)

fig_probs_df <- class_df |>
  filter(label %in% multi_class_labels)

fig_probs <- fig_probs_df |>
  ggplot(aes(x = max_prob, fill = class_label)) +
  geom_histogram(bins = 15, alpha = 0.75, color = "white") +
  scale_fill_manual(values = class_colors, name = "Class") +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey40") +
  facet_wrap(~ label, ncol = 2) +
  labs(
    title    = "Posterior Class Membership Probabilities",
    subtitle = "Dashed: 0.7 threshold for high-confidence assignment",
    x        = "Maximum posterior probability",
    y        = "N patients"
  ) +
  theme_research()

fsz_b <- get_figure_size("bar")
if (nrow(fig_probs_df) > 0) {
  ggsave(file.path(out_dir, "fig_membership_probs.pdf"),
         fig_probs, width = fsz_b["width"], height = fsz_b["height"])
  ggsave(file.path(out_dir, "fig_membership_probs.png"),
         fig_probs, width = fsz_b["width"], height = fsz_b["height"], dpi = 300)
  message("Saved: fig_membership_probs")
} else {
  message("Skipping fig_membership_probs: all biomarkers have ng=1")
}

# -- Figure 3: Class characterisation by amyloid risk -------------------------

char_df <- class_df |>
  left_join(v1_amyloid, by = "record_id") |>
  filter(!is.na(ptau217_risk)) |>
  mutate(
    ptau217_risk = factor(ptau217_risk, levels = c("Low", "Intermediate", "High")),
    class_label  = factor(class_label, levels = paste0("Class ", 1:3))
  )

if (nrow(char_df) > 0 && any(map_int(best_models, "ng") > 1)) {
  fig_char <- char_df |>
    filter(label %in% map_chr(
      names(best_models)[map_int(best_models, "ng") > 1],
      ~ biomarkers[[.x]]$label
    )) |>
    count(label, class_label, ptau217_risk) |>
    group_by(label, class_label) |>
    mutate(pct = n / sum(n) * 100) |>
    ungroup() |>
    ggplot(aes(x = class_label, y = pct, fill = ptau217_risk)) +
    geom_col(position = "stack", width = 0.6) +
    geom_text(aes(label = if_else(n >= 2, sprintf("n=%d", n), "")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("High" = "#E41A1C", "Intermediate" = "#FF7F00", "Low" = "#4DAF4A"),
      name   = "pTau217 risk"
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~ label, ncol = 2) +
    labs(
      title    = "Amyloid Risk Composition by Trajectory Class",
      subtitle = "pTau217 risk category within each latent class",
      x        = "Trajectory class",
      y        = "% patients"
    ) +
    theme_research()

  ggsave(file.path(out_dir, "fig_class_characterisation.pdf"),
         fig_char, width = fsz_b["width"], height = fsz_b["height"])
  ggsave(file.path(out_dir, "fig_class_characterisation.png"),
         fig_char, width = fsz_b["width"], height = fsz_b["height"], dpi = 300)
  message("Saved: fig_class_characterisation")
}

# -- Write CSVs ----------------------------------------------------------------

write_csv(sel_all,   file.path(out_dir, "model_selection.csv"))
write_csv(class_df,  file.path(out_dir, "class_assignments.csv"))

message("\nDone. Outputs written to: ", out_dir)
