# Research Notes: Plasma AD Biomarker Longitudinal Variability

## Research Question

What is the magnitude and character of short-term (weeks to months) biological variability in plasma AD biomarkers, and how does it compare to the preanalytical noise floor established by Mayfield et al. (20% TEa for pTau217/Ab1-42)?

## Cohort

- **Source**: DAC CRF, Mount Sinai (REDCap project 21855, PI: Fanny Elahi)
- **N**: 738 enrolled, 51 with >=2 visits, 23 with all 4 visits
- **Variability sub-study**: 40 consented
- **Visit timing**: V1->V2 median 32d (range 0-407d); V2->V3 and V3->V4 ~14d median
- **Demographics**: Mean age 77.2 (SD 10.2), 62% female, 59% White, 15% Black/AA

## Key Findings (Exploratory, 2026-04-07)

### 1. No systematic directional trend

| Biomarker | N | Increasing | Decreasing | Flat (<5%) | Median % Change |
|-----------|---|------------|------------|------------|-----------------|
| pTau217 | 42 | 40% | 40% | 19% | -0.3% |
| NfL | 42 | 33% | 40% | 26% | -0.7% |
| GFAP | 10 | 40% | 30% | 30% | +2.4% |
| Ab42/40 | 10 | 20% | 20% | 60% | +2.9% |
| Lucent AD | 39 | 28% | 38% | 33% | 0.0% |
| pTau217/Ab42-40 | 10 | 40% | 40% | 20% | -2.7% |

Biomarkers are not systematically increasing or decreasing. Changes scatter symmetrically around zero.

### 2. Oscillatory, not progressive

Among patients with 3+ visits:
- **pTau217**: 74% non-monotone (bounce), 15% monotone up, 11% monotone down
- **NfL**: 77% non-monotone, 0% monotone up, 23% monotone down
- **Lucent AD**: 61% non-monotone, 21% monotone up, 18% monotone down

Most patients fluctuate around a set point rather than drifting in one direction.

### 3. ~25% biological variability floor

pTau217 mean |% change| is ~25% regardless of inter-visit interval:
- <=30d: 25.1%
- 31-90d: 26.0%
- >90d: 23.6%

This sits at the 20% preanalytical TEa boundary (Mayfield et al.), meaning biological noise is roughly equal to or slightly exceeds the preanalytical floor.

### 4. TEa exceedance is common

| Biomarker | Visit | N | % Exceeding 20% TEa |
|-----------|-------|---|---------------------|
| pTau217 | V2 | 41 | 46% |
| pTau217 | V3 | 24 | 71% |
| pTau217 | V4 | 15 | 60% |
| NfL | V2 | 41 | 32% |
| Lucent AD | V2 | 33 | 42% |

Nearly half of serial pTau217 measurements exceed the preanalytical total allowable error threshold purely from biological variability.

### 5. NfL may show genuine drift at longer intervals

NfL mean |% change| jumps from ~16% (<=90d) to 46% (>90d). This subset may reflect progressive neurodegeneration signal, consistent with NfL's role as a neuronal injury rather than amyloid-specific marker.

### 6. Ab42/40 is the most stable individual analyte

60% of patients remain within 5% change. Mean |% change| <10%. This supports its use as the denominator in the pTau217/Ab ratio for variance reduction.

## Implications

1. **Clinical interpretation**: A single pTau217 measurement carries ~25% biological noise. Repeat-testing protocols should account for this -- only changes exceeding ~25-30% can be confidently attributed to biology rather than fluctuation.

2. **Clinical trial enrichment**: Cutoffs based on single measurements will misclassify ~46% of patients into different risk categories on repeat testing due to biological variability alone.

3. **Treatment monitoring**: Anti-amyloid therapy response monitoring with serial pTau217 requires changes substantially exceeding the 25% biological noise band to be interpretable.

4. **Ratio advantage**: The pTau217/Ab1-42 ratio reduces preanalytical noise (Mayfield et al.), but we have insufficient data (N=6 paired) to determine whether it also reduces biological variability. This is a critical gap.

## Companion Manuscript

`Stability manuscript draft 1.12.26.docx` -- Mayfield, John, Pennington, Kemna, Kroeger, Dasgupta, Burns, Elahi, Nelson, Morris (KU ADRC). Establishes preanalytical stability of pTau217/Ab1-42 ratio at 48h fresh / 24h frozen with 20% TEa criterion. Our biological variability analysis builds directly on this preanalytical foundation.

## Open Questions

- Does the pTau217/Ab1-42 ratio reduce biological variability as it does preanalytical noise? (Need more Lucent panel overlap data)
- What patient characteristics predict high vs. low intra-individual variability? (Age, diagnosis, amyloid status, comorbidities)
- Is the ~25% floor consistent across amyloid-positive vs. amyloid-negative subgroups?
- Does the variability sub-study (N=40, tighter cadence) show the same oscillatory pattern?
- How does this biological variability compare to published test-retest reliability (ICC) for Lumipulse pTau217?

## Pipeline

```bash
# Pull from REDCap and reshape
python analysis/dac_pull.py --output-dir data/

# Run EDA with TEa-annotated delta analysis
Rscript analysis/dac_eda.R --data-dir data/ --output-dir results/
```

## Outputs

| File | Content |
|------|---------|
| results/table1.csv, table1.html | Cohort characteristics (all V1 vs longitudinal) |
| results/enrollment_funnel.png | 735 -> 51 -> 37 -> 23 visit attrition |
| results/visit_timing.png | Inter-visit interval distributions |
| results/temporal_windows.csv | Alignment window counts |
| results/biomarker_trajectories.png | Spaghetti plots, 6 biomarkers |
| results/biomarker_distributions.png | Per-visit boxplots |
| results/clinical_scores.png | Baseline MoCA, MMSE, Linus, bADL |
| results/delta_analysis.png | Delta distributions + % change + TEa bands + time dependency |
| results/delta_summary.csv | Delta statistics with TEa exceedance counts |
| results/vitals_trajectories.png | BMI, BP longitudinal |
| results/baseline_correlations.png | 12x12 pairwise correlations |
