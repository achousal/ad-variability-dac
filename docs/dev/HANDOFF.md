---
updated: "2026-04-08T00:00"
project: "ad-variability-dac"
---

## What I Was Doing

Built the complete longitudinal variability analysis pipeline across five analysis modules:
EDA -> ICC/CV -> RCI -> Amyloid stratification -> Random slopes mixed model -> Latent class trajectory analysis.

## Current State

All scripts running end-to-end. No commits made -- all work on disk.

| Script | Output dir | Status |
|--------|-----------|--------|
| `dac_pull.py` | `data/` | Done |
| `dac_eda.R` | `results/dac_eda/` | Done |
| `dac_icc_cv.R` | `results/dac_icc_cv/` | Done |
| `dac_rci.R` | `results/dac_rci/` | Done |
| `dac_amyloid_strat.R` | `results/dac_amyloid_strat/` | Done |
| `dac_mixed_model.R` | `results/dac_mixed_model/` | Done |
| `dac_lcmm.R` | `results/dac_lcmm/` | Done |

## Key Findings (complete summary)

### EDA
- No systematic directional trend. pTau217: 40% up / 40% down / 19% flat.
- 74% of patients with >=3 visits are non-monotone (oscillatory trajectories).
- ~25% mean absolute % change for pTau217 regardless of interval length.
- NfL hint of drift at >90d (|%change| jumps from ~16% to 46%).

### ICC + CV
- pTau217 ICC=0.893 (Good), NfL=0.832 (Good). GFAP/Ab42-40/Lucent AD: Excellent (0.94-0.96).
- High ICC = stable ranking; pTau217 median CV 16.8%, 36% above 20% TEa.
- No significant CV predictors at this N (sex trending, rho=0.27, p=0.09).

### RCI (V1->V2, V1->V3, V1->V4)
- Very few reliable changers: 0-10% across all visit pairs.
- SD_diff stable across visit pairs -- noise floor does not shrink with longer follow-up.
- NfL 10% reliable changers at V1->V4 (N=20); fragile.
- pTau217 0% at V1->V4. Multi-visit comparison was key finding.

### Amyloid stratification (pTau217 risk: High vs Low)
- No difference in pTau217 or NfL CV by risk group (KW p>0.5).
- Lucent AD CV: High-risk patients have *lower* CV (2.8% vs 14.3%, p=0.007) -- likely ceiling/floor compression.
- pTau217 RCI exceedance: 20% in High vs 0% in Low (Fisher p=0.16, underpowered N=10).
- Monotonicity: High-risk patients are *more* non-monotone (89%), not less.

### Random slopes mixed model
- No significant population-level time trend (all LRT p>0.5).
- pTau217 and NfL: M2 (random slopes) converges without singularity -- patients vary in trajectory shape but distribution centered near zero.
- Lucent AD M2 singular -- slope variance collapses to zero; variability is intercept-driven.

### Latent class trajectory analysis (lcmm::hlme, 1-3 classes)
- pTau217 and NfL: BIC marginally favors 2 classes but minimum class sizes are 4 and 2 patients -- unstable, not interpretable.
- Lucent AD: 3-class solution preferred (BIC 1095->1053, delta=42), entropy=0.96, min class n=6. Only biomarker with meaningful trajectory subgroups.
- Conclusion: oscillatory pattern is cohort-wide for pTau217/NfL, not a mixture of hidden progressors.

## Coherent narrative across all analyses

At clinical visit timescales (weeks to ~1 year), plasma AD biomarker variability is dominated by biological oscillation, not progressive disease signal:
- Rankings are stable (high ICC), absolute values are noisy (high CV)
- Very few patients show reliable change by any method
- The noise floor does not shrink with longer follow-up
- No distinct "progressor" subgroup separable from the noise
- Amyloid burden does not appear to modulate biological variability at these timescales

Clinical implication: a single pTau217 measurement reliably classifies patients (good ICC) but serial measurements cannot be interpreted as directional without a large change (well exceeding the 20% TEa + biological noise). The pTau217/Ab42-40 ratio (Mayfield's primary finding) may reduce preanalytical noise but this cannot be evaluated here (N=6 with paired data).

### Clinical anchoring (Analysis 6, 2026-04-08)
- No longitudinal cognitive data -- MoCA baseline-only; delta-delta impossible.
- MoCA vs baseline Lucent AD: rho=-0.40, p=0.04 (only significant biomarker-cognition correlation).
- MoCA vs Lucent AD CV: rho=+0.53, p=0.003 -- higher cognitive impairment -> higher Lucent AD variability.
- pTau217/NfL show no significant MoCA correlations.
- RCI changers (n=2-3 per biomarker) too few for MoCA comparisons.

### Ratio variability (Analysis 7, 2026-04-08)
- N=10 patients with >=2 paired pTau217+Ab42-40 measurements (more than expected).
- Ratio median CV 8.49% vs pTau217 alone 8.36% -- ratio does NOT reduce biological variability at N=10.
- Ab42/40 remains the most stable component (median CV 2.85%).
- CAVEAT: Ab42-40 here is Lucent composite ratio, not raw Lumipulse Ab1-42 -- not directly comparable to Mayfield.

### Variability sub-study cohort (2026-04-08)
- 40 consented patients, all longitudinal. Interval groups: Short <=30d (n=22), Medium 31-90d (n=11), Long >90d (n=7).
- **pTau217 CV is completely interval-independent**: 16.8% / 16.3% / 19.3% across Short/Medium/Long (KW p=0.92). Biological oscillation operates at timescales shorter than the shortest measured interval (~2 weeks).
- NfL and Lucent AD similarly non-significant (p=0.49, p=0.91).
- **4-visit subgroup (N=23)**: pTau217 80% non-monotone, NfL 78% non-monotone. NfL has 22% monotone-down (no monotone-up) -- possible mild directional decline signal in a subset.
- RCI reliable changers cluster in Long group (14-17%) vs Short (5%) but n=1 patient each.

## Narrative Figures Session (2026-04-08)

Built `analysis/dac_narrative_figures.R` — single entry point producing 4 figures to `results/dac_narrative/`:

| Figure | File | Content |
|--------|------|---------|
| Fig 1 | `fig1_noise_floor_fixed` | 4-panel: CV vs interval (A), ICC vs CV (B), noise budget with MDC legend (C), SD_diff stability (D) |
| Fig 2 | `fig2_signal_within_noise` | NfL monotonicity stacked bar (A) + MoCA vs Lucent AD CV scatter (B) |
| Fig 3 | `fig3_monotone_trajectories` | Spaghetti % change from V1 by biomarker + MoCA boxplot by NfL class |
| Fig 4 | `fig4_stratified_monotone` | Spaghetti with NfL Low/High linetype + stacked bars (NfL × class, Lucent AD × class) |

**Design decisions:**
- Core color trio: purple=pTau217, orange=NfL, green=Lucent AD (elahi[4/5/3]); GFAP=pink, Ab42/40=grey for panel B only
- Panel 1C MDC values moved from in-figure text to legend
- **Monotonicity rule relaxed**: 1 step in wrong direction allowed if V1→V4 direction conserved. Strict: N=2 monotone; relaxed: N=14 (8 down, 6 up, 5 non-monotone)
- Fig 4 restricted to N=16 with V1 Lucent AD (consistent N across both stratifiers)
- NfL Low/High: mean split (26.8 pg/mL) — matches ad-prot-adrc convention
- Lucent AD Low/High: median split (cut=36) — no validated threshold; labeled exploratory
- No PET or diagnosis data exists for the 4-visit subgroup (all NA) — Lucent AD is only amyloid proxy

**Fig 4 contingency tables (N=16):**
- NfL Low: 1 up / 4 non-mono / 4 down; NfL High: 3 up / 1 non-mono / 3 down
- Lucent AD Low: 3 up / 2 non-mono / 3 down; High: 1 up / 3 non-mono / 4 down

## Next Steps

1. **Manuscript**: all analyses complete. Core story: interval-independent biological noise floor (sub-study result is the strongest statement), ICC, RCI, LCMM negative.
2. **Fig 4 interpretation**: High-NfL patients skew monotone-up (worsening, biologically expected). Lucent AD High skews monotone-down — may reflect ceiling/compression artifact (same pattern seen in amyloid strat analysis). Flag this in presentation.
3. **Potential Fig 5**: cross-biomarker concordance — of 8 monotone-down NfL patients, how many are also monotone in pTau217 / Lucent AD?
4. **Consider /research** — seed hypotheses from oscillatory variability and NfL monotone-down subset
5. **Commit everything** — no git history yet

## Key Decisions / Technical Notes

- Redraws (is_redraw=True) excluded from all pair analyses; 2 patients (2494, 2534) had duplicate V2 rows
- RCI computed for V1->V2, V1->V3, V1->V4 (not just V1->V2 as originally planned)
- Log-scale RCI as sensitivity for pTau217/NfL; 90-100% agreement with raw scale
- LCMM: gridsearch() was broken; replaced with manual 30-start loop using make_b_init()
- LCMM min class size floor: 5 patients (ng solutions with smaller classes excluded from selection)
- lmerTest loaded after lme4 for Satterthwaite p-values on fixed effects
