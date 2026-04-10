# Analysis Plan: Individual-Level Biomarker Variability

## Motivation

Group-level descriptive analysis (EDA, 2026-04-07) shows no systematic directional trend in plasma AD biomarkers across serial visits. Mean deltas are near zero, but individual-level variability is high (~25% for pTau217) and 74% of trajectories are non-monotone. Group averages wash out because patients go in opposite directions. We need methods that respect individual trajectories and separate "true changers" from "oscillators."

## Analytical Strategies

### 1. Within-Patient Coefficient of Variation (CV)

**Question**: Which patients are genuinely variable vs. stable?

Compute each patient's CV (SD/mean) across visits for each biomarker. Cluster patients into low-CV (stable set point, measurement noise only) vs. high-CV (genuinely variable). Then model: what predicts high CV? Candidate predictors: amyloid status, diagnosis, age, sex, comorbidities (eGFR, HTN, DM2), inter-visit interval.

- **Data requirement**: >=2 visits per patient (N~42-51)
- **Output**: CV distributions per biomarker, predictor analysis

### 2. Reliable Change Index (RCI)

**Question**: Which individual patients show changes exceeding measurement noise?

Standard in neuropsychology for distinguishing true change from error:

```
RCI = (V2 - V1) / (SD_diff * sqrt(2))
```

`SD_diff` estimated from V2-V1 pairs in the full cohort. |RCI| > 1.96 = reliable change at p < 0.05. Mayfield et al. provides the preanalytical component of SD_diff (TEa = 20%); our data contributes the biological component.

- **Data requirement**: V1-V2 pairs (N~42-53)
- **Output**: Count and proportion of "reliable changers" per biomarker, direction of reliable changes

### 3. Stratify by Amyloid Status

**Question**: Does biological variability differ between amyloid-positive and amyloid-negative patients?

Use baseline amyloid classification from:
- PET results (28 positive, 23 negative)
- pTau217 risk category (Low/Intermediate/High, N~695)
- Lucent AD interpretation (High/Intermediate/Low, N~305)

Compare CV, RCI exceedance rates, and trajectory patterns between groups. Amyloid-positive patients on a degenerative trajectory should show directional NfL/pTau217 drift; amyloid-negative patients should oscillate. If both groups bounce equally, biological variability dominates even in active disease at these timescales.

- **Data requirement**: Amyloid status + >=2 visits
- **Output**: Stratified CV, RCI, and trajectory comparisons

### 4. Intraclass Correlation Coefficient (ICC)

**Question**: How reliably do serial measurements rank patients?

ICC(3,1) from a two-way mixed model quantifies the proportion of total variance attributable to between-patient differences vs. within-patient fluctuation.

| ICC | Interpretation |
|-----|---------------|
| >0.90 | Excellent reliability -- single measurement sufficient |
| 0.75-0.90 | Good -- single measurement usually adequate |
| 0.50-0.75 | Moderate -- repeat testing recommended |
| <0.50 | Poor -- single measurement unreliable for classification |

This is the number clinicians need for repeat-testing guidelines.

- **Data requirement**: >=2 visits per patient
- **Output**: ICC per biomarker, with 95% CI

### 5. Patient-Level Slopes with Shrinkage (Random Effects Model)

**Question**: What is the distribution of individual rates of change after accounting for noise?

Linear mixed model with random intercepts + random slopes:

```r
lmer(ptau217 ~ days_from_v1 + (1 + days_from_v1 | record_id))
```

Random slopes give each patient's estimated rate of change, shrunk toward the group mean. Better than raw per-patient regression because it borrows strength across patients and regularizes noisy 2-3 point trajectories. Plot the distribution of patient-level slopes: centered at zero with wide tails, or distinct clusters?

- **Data requirement**: >=2 visits per patient with days_from_v1
- **Output**: Distribution of patient-level slopes, variance components, slope vs. baseline scatter

### 6. Latent Trajectory Classes (Group-Based Trajectory Modeling)

**Question**: Are there distinct subgroups following different patterns?

Fit 2-3 latent trajectory classes (e.g., "stable," "increasing," "decreasing") using `lcmm` or `flexmix` in R. Even with N~42, a 2-class model (stable vs. changing) is feasible. Compare classes on clinical characteristics.

- **Data requirement**: >=2 visits, ideally 3+ (N~27-30 with 3+ visits)
- **Output**: Class assignments, class-specific trajectories, clinical predictors of class membership

### 7. Anchor to Clinical Change

**Question**: Do biomarker changes track clinical changes?

Link biomarker deltas to:
- MoCA/MMSE deltas (limited -- mostly single-timepoint, but `most_recent_moca` field may capture updates)
- Clinical milestones: diagnosis change, treatment initiation (anti-amyloid vs. symptomatic), PET conversion
- Functional status (bADL, N=33)

A biomarker change that tracks cognitive decline is signal; one that doesn't is noise. Limited by single-timepoint cognitive data in current REDCap structure.

- **Data requirement**: Paired biomarker + clinical score data
- **Output**: Correlation of deltas, concordance of direction

## Priority Order

| Priority | Analysis | Answers | Feasibility |
|----------|----------|---------|-------------|
| 1 | ICC + CV | "How reliable? Who is variable?" | Immediate |
| 2 | RCI | "Which patients truly changed?" | Immediate |
| 3 | Amyloid stratification | "Does disease state modulate variability?" | Immediate |
| 4 | Random slopes model | "What's the distribution of individual trajectories?" | Short-term |
| 5 | Latent trajectory classes | "Are there distinct subgroups?" | Short-term, needs >=3 visits |
| 6 | Clinical anchoring | "Does biomarker change = clinical change?" | Limited by data |
| 7 | Ratio variability | "Does pTau217/Ab reduce biological noise?" | Blocked -- need more Lucent panel overlap |

## Dependencies

- R: `lme4` (mixed models), `irr` or `psych` (ICC), `lcmm` (trajectory classes)
- Existing: `dac_long.csv`, `dac_baseline.csv` from current pipeline

## Relation to Companion Manuscript

Mayfield et al. establishes the preanalytical noise floor (20% TEa, 48h fresh stability). This analysis plan builds the biological variability layer on top. Together they answer: **of the total variability in a serial pTau217 measurement, how much is preanalytical vs. biological, and can we identify patients whose changes exceed both noise sources?**
