#!/usr/bin/env python3
"""dac_pull.py -- Pull DAC CRF data from REDCap and reshape for longitudinal analysis.

Connects to REDCap API (Mount Sinai DAC CRF project), pulls all relevant forms,
reshapes the classic (non-longitudinal) flat record structure into:
  - dac_long.csv:     one row per patient-visit (longitudinal biomarkers + vitals)
  - dac_baseline.csv: one row per patient (demographics + clinical scores + labs)
  - dac_data_dictionary.csv: field mapping for reproducibility

Usage:
    python dac_pull.py [--output-dir PATH] [--no-pull]

Requires REDCAP_URL and REDCAP_TOKEN in environment or _code/.env.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import requests
except ImportError:
    sys.exit("requests library required: pip install requests")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# REDCap field mapping: canonical name -> {visit: redcap_field}
# ---------------------------------------------------------------------------

VISIT_FIELDS: dict[str, dict[str, str]] = {
    "visit_date": {
        "V1": "pat_enroll",
        "V2": "pat_enroll_b5f0e3",
        "V3": "pat_enroll_1d8ed6",
        "V4": "pat_enroll_v4",
        "Redraw": "date_of_initial_visit_294463",
    },
    "blood_date": {
        "V1": "date_of_initial_visit",
        "V2": "date_of_initial_visit_08f636",
        "V3": "date_of_initial_visit_46f612",
        "V4": "date_of_initial_visit_v4",
        "Redraw": "date_of_initial_visit_294463",
    },
    "blood_time": {
        "V1": "time_blood_drawn",
        "V2": "time_blood_drawn_679688",
        "V3": "time_blood_drawn_a2f04e",
        "V4": "time_blood_drawn_v4",
        "Redraw": "time_blood_drawn_666c60",
    },
    "bmi": {
        "V1": "enc_bmi",
        "V2": "enc_bmi_36bd90",
        "V3": "enc_bmi_be0f87",
        "V4": "enc_bmi_v4",
        "Redraw": "enc_bmi_3f2871",
    },
    "systolic": {
        "V1": "enc_systolic",
        "V2": "enc_systolic_7c06fe",
        "V3": "enc_systolic_addb5d",
        "V4": "enc_systolic_v4",
        "Redraw": "enc_systolic_414821",
    },
    "diastolic": {
        "V1": "enc_diastolic",
        "V2": "enc_diastolic_5b05b4",
        "V3": "enc_diastolic_3cd377",
        "V4": "enc_diastolic_v4",
        "Redraw": "enc_diastolic_90ec36",
    },
    "ptau217": {
        "V1": "ptau217_numerical_results_a4ce1e",
        "V2": "ptau217_numerical_results_a4ce1e_0b2a4e",
        "V3": "ptau217_numerical_results_a4ce1e_0b2a4e_522059",
        "V4": "ptau217_numerical_results_a4ce1e_0b2a4e_522059_365877",
    },
    "ptau217_risk": {
        "V1": "ptau_217",
        "V2": "ptau_217_a07651",
        "V3": "ptau_217_a07651_13538c",
        "V4": "ptau_217_a07651_13538c_bbc4c2",
    },
    "nfl": {
        "V1": "nfl_pdf_numerical_results_d30f3f",
        "V2": "nfl_pdf_numerical_results_d30f3f_143f97",
        "V3": "nfl_pdf_numerical_results_d30f3f_143f97_38088c",
        "V4": "nfl_pdf_numerical_results_d30f3f_143f97_38088c_7b5015",
    },
    "gfap": {
        "V1": "gfap_numerical_results_6bbd30",
        "V2": "gfap_numerical_results_6bbd30_d1485c",
        "V3": "gfap_numerical_results_6bbd30_d1485c_f10748",
        "V4": "gfap_numerical_results_6bbd30_d1485c_f10748_3e7f1f",
    },
    "ab42_40": {
        "V1": "numerical_results_6bbd31",
        "V2": "gfap_numerical_results_6bbd30_d1485c_2",
        "V3": "gfap_numerical_results_6bbd30_d1485c_f10749",
        "V4": "gfap_numerical_results_6bbd30_d1485c_f10748_3e7f1f_2",
    },
    "lucent_ad": {
        "V1": "lucent_ad",
        "V2": "lucent_ad_92539d",
        "V3": "lucent_ad_92539d_be8973",
        "V4": "lucent_ad_92539d_be8973_e34ab0",
    },
    "lucent_interpretation": {
        "V1": "lucent_ad_interpretation",
        "V2": "lucent_ad_interpretation_1424db",
        "V3": "lucent_ad_interpretation_1424db_8f57f6",
        "V4": "lucent_ad_interpretation_1424db_8f57f6_830e67",
    },
    "physician": {
        "V1": "physician",
        "V2": "physician_cfb80b",
        "V3": "physician_cb1f55",
        "V4": "physician_v4",
    },
    "clinic": {
        "V1": "clinic",
        "V2": "clinic_d4f9ef",
        "V3": "clinic_1c54c4",
        "V4": "",
    },
}

# Demographic fields (baseline only, wide)
DEMO_FIELDS = [
    "record_id", "age_at_enrol", "sex", "pat_edu", "symptomonset",
    "htn", "hld", "afib", "heartdisease", "stroke", "parkinsons",
    "ms", "dm1", "dm2", "sleepapnea", "rbd", "sleep", "cancer",
    "ckd", "liver", "autoimmune", "majordepression", "bipolar",
    "schizophrenia", "anxiety", "ptsd", "personality",
    "hearing_impairment", "hearing_impairment_2",
]

# Clinical score fields (baseline only)
CLINICAL_FIELDS = [
    "moca_raw_score", "raw_mmse_score", "moca_score", "mmse_total_score",
    "most_recent_moca", "moca_score_date", "recent_mmse", "mmse_score_date_2",
    "badl_summary_of_6", "clinical_visit_date",
    "totalscore", "clockdrawing", "delayedrecall", "dementia_risk_forecast",
]

# Follow-up fields (baseline only -- on V1 study_follow_up form)
FOLLOWUP_FIELDS = [
    "apoe4_genotype", "apoe4_status",
    "results_pet", "centiloid_score_pet", "avg_suvr_score", "pet_scan_date",
    "pat_dropout", "pat_dropout_date", "deceased", "ltfu",
    "case_management",
    "mri", "date_of_mri",
    "date_of_results",
]

# Lab fields
LAB_FIELDS = [
    "cmb", "cmbdate", "egfr", "creatinine", "fasting_glucose",
    "lipid", "lipiddate", "ldl", "hdl",
    "gluca1c", "gluca1cdate", "glua1cval",
    "vitaminb12", "b12value", "vitamind", "vitamindvalue",
    "lipoprotein_date", "lipoproteina",
]

# V2-specific fields
V2_EXTRA_FIELDS = [
    "varibaility_study", "cohort_stratification", "referral",
]

# Categorical code maps
SEX_MAP = {"1": "Male", "2": "Female"}
EDUCATION_MAP = {
    "1": "No high school", "2": "Some high school",
    "3": "High school diploma", "4": "2-year degree",
    "5": "4-year degree", "6": "Graduate+",
}
RACE_CODES = {1: "White", 2: "Black/AA", 3: "Asian/PI", 4: "AI/AN", 5: "Other"}
ETHNICITY_CODES = {1: "Hispanic", 2: "Non-Hispanic"}
PTAU_RISK_MAP = {"1": "Low", "2": "Intermediate", "3": "High"}
LUCENT_INTERP_MAP = {"1": "High", "2": "Intermediate", "3": "Low"}
APOE_MAP = {
    "1": "E2/E2", "2": "E2/E3", "3": "E2/E4", "9": "E3/E2",
    "4": "E3/E3", "5": "E3/E4", "6": "E4/E4", "7": "Unknown", "8": "Other",
}
DIAGNOSIS_CODES = {
    1: "MCI", 2: "AD", 3: "Aging", 4: "Psych", 5: "Other neurodegen",
    6: "LBD", 7: "FTD", 8: "Parkinsons", 9: "No dx/other priority", 10: "Mixed",
}
PET_MAP = {"1": "Positive", "2": "Negative", "3": "Indeterminate", "4": "Intermediate"}
SYMPT_MAP = {"1": "Symptomatic", "2": "Asymptomatic", "3": "Not sure"}
YESNO_MAP = {"1": "Yes", "0": "No"}
COMORBIDITY_MAP = {"1": "Yes", "2": "No"}


# ---------------------------------------------------------------------------
# REDCap API helpers
# ---------------------------------------------------------------------------

def load_env(env_path: Path) -> None:
    """Source a .env file into os.environ (simple key=value, no export)."""
    if not env_path.exists():
        return
    with open(env_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, val = line.partition("=")
            os.environ.setdefault(key.strip(), val.strip())


def redcap_export(url: str, token: str, **params: Any) -> list[dict]:
    """POST to REDCap API and return parsed JSON."""
    payload = {"token": token, "format": "json", "returnFormat": "json", **params}
    resp = requests.post(url, data=payload, timeout=120)
    resp.raise_for_status()
    data = resp.json()
    if isinstance(data, dict) and "error" in data:
        raise RuntimeError(f"REDCap API error: {data['error']}")
    return data


def pull_by_forms(
    url: str, token: str, forms: list[str], extra_fields: list[str] | None = None,
) -> list[dict]:
    """Pull records for specific forms, always including record_id."""
    payload: dict[str, str] = {"content": "record"}
    for i, form in enumerate(forms):
        payload[f"forms[{i}]"] = form
    # record_id lives on the first form; force-include it so every export has it
    fi = 0
    for field in (extra_fields or ["record_id"]):
        payload[f"fields[{fi}]"] = field
        fi += 1
    if "record_id" not in (extra_fields or []):
        payload[f"fields[{fi}]"] = "record_id"
    return redcap_export(url, token, **payload)


# ---------------------------------------------------------------------------
# Data reshaping
# ---------------------------------------------------------------------------

def _val(record: dict, field: str) -> str | None:
    """Get a non-empty string value or None."""
    v = record.get(field, "")
    if v is None or str(v).strip() in ("", "NaN", "nan"):
        return None
    return str(v).strip()


def _float(record: dict, field: str) -> float | None:
    """Get a float value or None."""
    v = _val(record, field)
    if v is None:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _date(record: dict, field: str) -> str | None:
    """Get a date string (YYYY-MM-DD) or None."""
    v = _val(record, field)
    if v is None:
        return None
    # Validate date format
    try:
        datetime.strptime(v, "%Y-%m-%d")
        return v
    except ValueError:
        return None


def build_long_rows(record: dict) -> list[dict]:
    """Reshape one flat REDCap record into visit-level rows."""
    rid = record.get("record_id", "")
    rows = []

    for visit in ("V1", "V2", "V3", "V4", "Redraw"):
        # Check if this visit has a date
        date_field = VISIT_FIELDS["visit_date"].get(visit, "")
        if not date_field or not _val(record, date_field):
            continue

        row: dict[str, Any] = {"record_id": rid, "visit": visit}

        # Extract all visit-level fields
        for canonical, field_map in VISIT_FIELDS.items():
            rc_field = field_map.get(visit, "")
            if not rc_field:
                row[canonical] = None
                continue

            if canonical in ("visit_date", "blood_date"):
                row[canonical] = _date(record, rc_field)
            elif canonical in (
                "bmi", "systolic", "diastolic",
                "ptau217", "nfl", "gfap", "ab42_40", "lucent_ad",
            ):
                row[canonical] = _float(record, rc_field)
            elif canonical == "ptau217_risk":
                row[canonical] = PTAU_RISK_MAP.get(_val(record, rc_field) or "", None)
            elif canonical == "lucent_interpretation":
                row[canonical] = LUCENT_INTERP_MAP.get(
                    _val(record, rc_field) or "", None
                )
            elif canonical in ("blood_time",):
                row[canonical] = _val(record, rc_field)
            else:
                row[canonical] = _val(record, rc_field)

        # Derived: pTau217/Ab1-42 ratio (the stability-validated composite)
        if row.get("ptau217") is not None and row.get("ab42_40") is not None:
            # ab42_40 is already a ratio (0.02-0.09 range); ptau217 is pg/mL
            # The Lumipulse ratio uses raw pg/mL for both components
            # Here ab42_40 is the Ab42/40 ratio from Lucent, not raw Ab1-42
            # So this composite is only meaningful when both come from the same platform
            row["ptau217_ab42_ratio"] = row["ptau217"] / row["ab42_40"] if row["ab42_40"] != 0 else None
        else:
            row["ptau217_ab42_ratio"] = None

        # Mark redraws
        row["is_redraw"] = visit == "Redraw"

        rows.append(row)

    return rows


def merge_redraws(rows: list[dict]) -> list[dict]:
    """Relabel Redraw visits as V2 with a flag."""
    for row in rows:
        if row["visit"] == "Redraw":
            row["visit"] = "V2"
            row["is_redraw"] = True
    return rows


def compute_derived(rows: list[dict]) -> list[dict]:
    """Add days_from_v1, visit_interval, and age_at_visit."""
    # Group by record_id
    by_patient: dict[str, list[dict]] = {}
    for row in rows:
        by_patient.setdefault(row["record_id"], []).append(row)

    for rid, patient_rows in by_patient.items():
        # Sort by visit date
        patient_rows.sort(key=lambda r: r.get("visit_date") or "9999")

        v1_date = None
        for r in patient_rows:
            if r["visit"] == "V1" and r.get("visit_date"):
                v1_date = datetime.strptime(r["visit_date"], "%Y-%m-%d")
                break

        prev_date = None
        for r in patient_rows:
            vd = r.get("visit_date")
            if vd and v1_date:
                current = datetime.strptime(vd, "%Y-%m-%d")
                r["days_from_v1"] = (current - v1_date).days
            else:
                r["days_from_v1"] = None

            if vd and prev_date:
                current = datetime.strptime(vd, "%Y-%m-%d")
                r["visit_interval"] = (current - prev_date).days
            else:
                r["visit_interval"] = 0 if r["visit"] == "V1" else None

            if vd:
                prev_date = datetime.strptime(vd, "%Y-%m-%d")

    return rows


def build_baseline(
    records: list[dict],
    demo_data: list[dict],
    followup_data: list[dict],
    lab_data: list[dict],
    merged_records: dict[str, dict] | None = None,
) -> list[dict]:
    """Build one-row-per-patient baseline table."""
    # Index supplemental data by record_id
    demo_idx = {r["record_id"]: r for r in demo_data}
    fu_idx = {r["record_id"]: r for r in followup_data}
    lab_idx = {r["record_id"]: r for r in lab_data}
    # Merged records have all fields from all forms
    merged_idx = merged_records or {}

    baseline_rows = []

    all_ids = set()
    for source in (demo_idx, fu_idx, lab_idx):
        all_ids.update(source.keys())
    # Also include records from the main pull
    for r in records:
        all_ids.add(r["record_id"])

    for rid in sorted(all_ids, key=lambda x: int(x) if x.isdigit() else x):
        demo = demo_idx.get(rid, {})
        fu = fu_idx.get(rid, {})
        lab = lab_idx.get(rid, {})

        row: dict[str, Any] = {"record_id": rid}

        # Demographics
        row["age_at_enrol"] = _float(demo, "age_at_enrol")
        row["sex"] = SEX_MAP.get(_val(demo, "sex") or "", None)
        row["education"] = EDUCATION_MAP.get(_val(demo, "pat_edu") or "", None)
        row["symptom_onset_age"] = _float(demo, "symptomonset")

        # Race (checkbox -- multi-select)
        races = []
        for code, label in RACE_CODES.items():
            if demo.get(f"race___{code}") == "1":
                races.append(label)
        row["race"] = "; ".join(races) if races else None

        # Ethnicity (checkbox)
        ethnicities = []
        for code, label in ETHNICITY_CODES.items():
            if demo.get(f"ethnicity___{code}") == "1":
                ethnicities.append(label)
        row["ethnicity"] = "; ".join(ethnicities) if ethnicities else None

        # Comorbidities
        for field in [
            "htn", "hld", "afib", "heartdisease", "stroke", "parkinsons",
            "ms", "dm1", "dm2", "sleepapnea", "rbd", "sleep", "cancer",
            "ckd", "liver", "autoimmune", "majordepression", "bipolar",
            "schizophrenia", "anxiety", "ptsd", "personality",
            "hearing_impairment", "hearing_impairment_2",
        ]:
            val = _val(demo, field)
            row[field] = COMORBIDITY_MAP.get(val or "", None)

        # Clinical scores (from mocammseneuroexamadl + linus_health forms, in merged)
        m = merged_idx.get(rid, {})
        row["moca_raw"] = _float(m, "moca_raw_score")
        row["mmse_raw"] = _float(m, "raw_mmse_score")
        row["moca_calc"] = _float(m, "moca_score")
        row["mmse_calc"] = _float(m, "mmse_total_score")
        row["badl"] = _float(m, "badl_summary_of_6")
        row["linus_total"] = _float(m, "totalscore")
        row["linus_clock"] = _float(m, "clockdrawing")
        row["linus_delayed_recall"] = _float(m, "delayedrecall")

        # Diagnosis (checkbox on V1 follow-up)
        diagnoses = []
        for code, label in DIAGNOSIS_CODES.items():
            if fu.get(f"diagnosis___{code}") == "1":
                diagnoses.append(label)
        row["diagnosis"] = "; ".join(diagnoses) if diagnoses else None

        # APOE
        row["apoe_genotype"] = APOE_MAP.get(_val(fu, "apoe4_genotype") or "", None)

        # PET
        row["pet_result"] = PET_MAP.get(_val(fu, "results_pet") or "", None)
        row["centiloid"] = _float(fu, "centiloid_score_pet")
        row["avg_suvr"] = _float(fu, "avg_suvr_score")
        row["pet_date"] = _date(fu, "pet_scan_date")

        # Study status
        row["dropout"] = _val(fu, "pat_dropout") == "1"
        row["deceased"] = _val(fu, "deceased") == "1"
        row["ltfu"] = _val(fu, "ltfu") == "1"
        row["case_management"] = _val(fu, "case_management")

        # V2 extras (from V2 visit form, in merged)
        row["variability_consented"] = _val(m, "varibaility_study") == "1"
        row["cohort_stratification"] = SYMPT_MAP.get(
            _val(m, "cohort_stratification") or "", None
        )

        # Lab values
        row["egfr"] = _float(lab, "egfr")
        row["creatinine"] = _float(lab, "creatinine")
        row["fasting_glucose"] = _float(lab, "fasting_glucose")
        row["ldl"] = _float(lab, "ldl")
        row["hdl"] = _float(lab, "hdl")
        row["a1c"] = _float(lab, "glua1cval")
        row["vitamin_b12"] = _float(lab, "b12value")
        row["vitamin_d"] = _float(lab, "vitamindvalue")
        row["lipoprotein_a"] = _float(lab, "lipoproteina")

        baseline_rows.append(row)

    return baseline_rows


def export_data_dictionary() -> list[dict]:
    """Generate field mapping documentation."""
    rows = []
    for canonical, field_map in VISIT_FIELDS.items():
        for visit, rc_field in field_map.items():
            if rc_field:
                rows.append({
                    "canonical_name": canonical,
                    "visit": visit,
                    "redcap_field": rc_field,
                    "table": "long",
                })
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: _code/data/ relative to script)",
    )
    parser.add_argument(
        "--no-pull",
        action="store_true",
        help="Skip API pull, use cached /tmp/dac_data.json",
    )
    args = parser.parse_args()

    # Resolve paths
    script_dir = Path(__file__).resolve().parent
    code_dir = script_dir.parent
    if args.output_dir:
        out_dir = args.output_dir
    else:
        out_dir = code_dir / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load env
    env_path = code_dir / ".env"
    load_env(env_path)

    url = os.environ.get("REDCAP_URL", "")
    token = os.environ.get("REDCAP_TOKEN", "")
    if not url or not token:
        sys.exit("REDCAP_URL and REDCAP_TOKEN must be set in environment or _code/.env")

    # --- Pull data ---
    log.info("Pulling data from REDCap...")

    # Pull visit + follow-up forms (contains visit-level biomarkers)
    visit_forms = [
        "study_visit", "study_follow_up",
        "study_visit_d45fec", "study_follow_up_v2",
        "study_visit_9901e8", "study_follow_up_v3",
        "study_visit_v4", "study_follow_up_v4",
        "study_visit_1bfb53",
    ]
    log.info("  Pulling visit/follow-up forms...")
    visit_data = pull_by_forms(url, token, visit_forms)
    log.info(f"  Got {len(visit_data)} records from visit forms")

    # Pull demographics (separate to get checkbox fields)
    log.info("  Pulling demographics...")
    demo_data = pull_by_forms(url, token, ["demographics"])
    log.info(f"  Got {len(demo_data)} demographic records")

    # Pull clinical scores
    log.info("  Pulling clinical scores...")
    clinical_data = pull_by_forms(
        url, token, ["mocammseneuroexamadl", "linus_health"]
    )
    log.info(f"  Got {len(clinical_data)} clinical records")

    # Pull follow-up (for diagnosis checkboxes, APOE, PET, status)
    log.info("  Pulling V1 follow-up for diagnosis/APOE/PET...")
    followup_data = pull_by_forms(url, token, ["study_follow_up"])
    log.info(f"  Got {len(followup_data)} follow-up records")

    # Pull V2 follow-up for variability/stratification
    log.info("  Pulling V2 visit for variability/stratification...")
    v2_data = pull_by_forms(url, token, ["study_visit_d45fec"])
    log.info(f"  Got {len(v2_data)} V2 visit records")

    # Pull labs
    log.info("  Pulling laboratory information...")
    lab_data = pull_by_forms(url, token, ["laboratory_information"])
    log.info(f"  Got {len(lab_data)} lab records")

    # --- Merge all record data into a single dict per record_id ---
    log.info("Merging record data...")
    merged: dict[str, dict] = {}
    for source in (visit_data, demo_data, clinical_data, followup_data, v2_data, lab_data):
        for rec in source:
            rid = rec.get("record_id", "")
            if rid not in merged:
                merged[rid] = {}
            # Merge non-empty values
            for k, v in rec.items():
                if v not in ("", None) and (k not in merged[rid] or merged[rid][k] in ("", None)):
                    merged[rid][k] = v

    all_records = list(merged.values())
    log.info(f"Merged {len(all_records)} unique records")

    # --- Build long format ---
    log.info("Building longitudinal (long) format...")
    long_rows = []
    for record in all_records:
        long_rows.extend(build_long_rows(record))
    long_rows = merge_redraws(long_rows)
    long_rows = compute_derived(long_rows)
    log.info(f"Long format: {len(long_rows)} patient-visit rows")

    # --- Build baseline ---
    log.info("Building baseline (wide) format...")
    baseline_rows = build_baseline(
        all_records, demo_data, followup_data, lab_data, merged_records=merged,
    )
    log.info(f"Baseline format: {len(baseline_rows)} patients")

    # --- Export ---
    long_path = out_dir / "dac_long.csv"
    baseline_path = out_dir / "dac_baseline.csv"
    dict_path = out_dir / "dac_data_dictionary.csv"

    # Long
    if long_rows:
        long_keys = list(long_rows[0].keys())
        with open(long_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=long_keys)
            writer.writeheader()
            writer.writerows(long_rows)
        log.info(f"Wrote {long_path} ({len(long_rows)} rows)")

    # Baseline
    if baseline_rows:
        baseline_keys = list(baseline_rows[0].keys())
        with open(baseline_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=baseline_keys)
            writer.writeheader()
            writer.writerows(baseline_rows)
        log.info(f"Wrote {baseline_path} ({len(baseline_rows)} rows)")

    # Data dictionary
    dd_rows = export_data_dictionary()
    with open(dict_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["canonical_name", "visit", "redcap_field", "table"])
        writer.writeheader()
        writer.writerows(dd_rows)
    log.info(f"Wrote {dict_path} ({len(dd_rows)} mappings)")

    # --- Summary ---
    visits_count = {}
    for r in long_rows:
        v = r["visit"]
        visits_count[v] = visits_count.get(v, 0) + 1
    log.info("Visit counts: %s", visits_count)
    log.info("Done.")


if __name__ == "__main__":
    main()
