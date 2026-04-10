# ad-variability-dac

Blood biomarker variability study using the DAC CRF cohort (REDCap project 21855, Mount Sinai). Longitudinal descriptive analysis of plasma AD biomarkers across up to 4 serial visits.

## Study design

- **Cohort**: 738 patients enrolled through Dean/Alzheimer's Center clinics
- **Visits**: V1 (baseline), V2, V3, V4 at variable intervals; redraws treated as V2
- **Variability sub-study**: 40 consented patients with tighter visit cadence
- **PI**: Fanny Elahi

## Data source

REDCap API at `redcap.mountsinai.org`. Credentials in vault `_code/.env` (`REDCAP_URL`, `REDCAP_TOKEN`). The project is classic (non-longitudinal) -- visits are encoded as separate instruments with suffixed field names.

## Key biomarkers

| Analyte | V1 N | Platform |
|---------|------|----------|
| pTau217 | 695 | C2N/Lumipulse |
| NfL | 612 | C2N/Lumipulse |
| GFAP | 64 | Lucent panel |
| Ab42/40 | 64 | Lucent panel |
| Lucent AD | 305 | Lucent composite |

## Project structure

```
analysis/     Scripts (Python pull, R EDA, theme/palette helpers)
data/         Extracted CSVs (long, baseline, data dictionary)
results/      Figures and tables
```

## Running the pipeline

```bash
# 1. Pull from REDCap and reshape
python analysis/dac_pull.py --output-dir data/

# 2. Run EDA
Rscript analysis/dac_eda.R --data-dir data/ --output-dir results/
```

## Data notes

- **Cognitive scores** (MoCA, MMSE) are baseline-only (single `mocammseneuroexamadl` form, not per-visit)
- **GFAP and Ab42/40** are sparse (Lucent panel only, ~64 at V1)
- **Diagnosis** assigned for ~96/738; rest pending
- **Inter-visit timing** is highly variable, especially V1->V2 (median 32d, range 0-407d)
- **No PII** in exported CSVs -- MRN, DOB, contact info excluded by dac_pull.py

## Dependencies

- Python: `requests`
- R: `tidyverse`, `gt`, `patchwork`
