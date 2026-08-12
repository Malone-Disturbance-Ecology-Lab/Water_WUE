# -*- coding: utf-8 -*-
"""
Console-only updated EDI result extractor.

Purpose:
Print updated values for the Results section without writing new CSV files.

It checks:
1. Growing-season months
2. SPEI-3 hydroclimate summary
3. Spatial WUE_T response summary
4. Q3 fitted curve range
5. Q3 smooth support/significance
6. Site counts by coast
7. Sentence-ready values for manuscript update
"""

from pathlib import Path
import pandas as pd
import numpy as np

# ============================================================
# PATHS
# ============================================================
BASE = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6")

Q3_DIR = BASE / "Q3" / "Q3_WUE_T_SPEI_sensitivity_outputs"
Q4_DIR = BASE / "Q4" / "Q4_Ecological_Impacts_outputs"
SPATIAL_DIR = BASE / "upscaling"

monthly_spatial_file = SPATIAL_DIR / "EDI_response_v2_upland_monthly_coast_summary.csv"
annual_impact_file = SPATIAL_DIR / "EDI_response_communication_annual_impacted_area_by_coast_2000_2025.csv"
validation_file = SPATIAL_DIR / "EDI_response_v2_regional_validation_Q3_Q4.csv"

#q4_curve_file = Q4_DIR / "EDI_Q3_aligned_prediction_scores_upland.csv"

Q3_FINAL_DIR = BASE / "Q3" / "Q3_august_update"

q4_curve_file = (
    Q3_FINAL_DIR
    / "Q3_reviewer_coast_threshold_predictions_all_ecosystems_all_months.csv"
)

q4_monthly_scores_file = Q4_DIR / "EDI_Q3_aligned_monthly_scores.csv"

q3_smooth_file = Q3_DIR / "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv"
q3_site_sens_file = Q3_DIR / "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv"
q3_model_comp_file = Q3_DIR / "Q3_WUE_T_SPEI_model_comparison.csv"

COAST_ORDER = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]

DISPLAY = {
    "AK Coast": "Alaska Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast",
}

pd.set_option("display.max_columns", 100)
pd.set_option("display.width", 220)

# ============================================================
# HELPERS
# ============================================================
def read_csv_safe(path, required=True):
    if not path.exists():
        msg = f"Missing file: {path}"
        if required:
            raise FileNotFoundError(msg)
        print("WARNING:", msg)
        return pd.DataFrame()
    return pd.read_csv(path)

def fmt(x, digits=2):
    if pd.isna(x) or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"

def first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def order_coasts(df):
    if "coast_region" in df.columns:
        df = df.copy()
        df["coast_region"] = pd.Categorical(df["coast_region"], COAST_ORDER, ordered=True)
        df = df.sort_values("coast_region")
    return df

def print_section(title):
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)

# ============================================================
# READ FILES
# ============================================================
monthly = read_csv_safe(monthly_spatial_file)
annual_impact = read_csv_safe(annual_impact_file)
validation = read_csv_safe(validation_file, required=False)

q4_curve = read_csv_safe(q4_curve_file, required=False)
q4_monthly = read_csv_safe(q4_monthly_scores_file, required=False)

q3_smooth = read_csv_safe(q3_smooth_file, required=False)
q3_site_sens = read_csv_safe(q3_site_sens_file, required=False)
q3_model_comp = read_csv_safe(q3_model_comp_file, required=False)

# ============================================================
# 0. GROWING SEASON CHECK
# ============================================================
print_section("0. GROWING-SEASON CHECK")

months_present = sorted(monthly["month"].dropna().astype(int).unique())
years_present = sorted(monthly["year"].dropna().astype(int).unique())

print("Months present:", months_present)
print("Year range:", min(years_present), "-", max(years_present))
print("Rows:", len(monthly))
print("Expected rows for 26 years × 7 months × 4 coasts:", 26 * 7 * 4)

if months_present == [4, 5, 6, 7, 8, 9, 10]:
    print("PASS: Spatial results are April-October growing season only.")
else:
    print("WARNING: Spatial results are not April-October only.")

# ============================================================
# 1. SPEI-3 HYDROCLIMATE SUMMARY
# ============================================================
print_section("1. SPEI-3 HYDROCLIMATE SUMMARY BY COAST")

spei_summary = (
    monthly
    .groupby("coast_region")
    .agg(
        n_months=("mean_SPEI3", "count"),
        mean_SPEI3=("mean_SPEI3", "mean"),
        min_monthly_mean_SPEI3=("mean_SPEI3", "min"),
        max_monthly_mean_SPEI3=("mean_SPEI3", "max"),
        cumulative_SPEI3=("mean_SPEI3", "sum"),
        dry_month_pct=("mean_SPEI3", lambda x: (x <= -1).mean() * 100),
        wet_month_pct=("mean_SPEI3", lambda x: (x >= 1).mean() * 100),
    )
    .reset_index()
)

spei_summary = order_coasts(spei_summary)
print(spei_summary.round(3).to_string(index=False))

# ============================================================
# 2. SPATIAL RESPONSE SUMMARY
# ============================================================
print_section("2. SPATIAL WUE_T RESPONSE SUMMARY BY COAST")

impact_summary = (
    annual_impact
    .groupby("coast_region")
    .agg(
        mean_annual_negative_event_area_pct=("annual_mean_negative_impacted_area_pct", "mean"),
        mean_annual_positive_event_area_pct=("annual_mean_positive_impacted_area_pct", "mean"),
        mean_annual_any_strong_response_area_pct=("annual_mean_any_impacted_area_pct", "mean"),
        mean_annual_net_positive_minus_negative_pct=("annual_mean_net_positive_minus_negative_pct", "mean"),
        mean_annual_dry_area_pct=("annual_mean_dry_area_pct", "mean"),
        max_annual_negative_event_area_pct=("annual_mean_negative_impacted_area_pct", "max"),
        max_annual_positive_event_area_pct=("annual_mean_positive_impacted_area_pct", "max"),
        max_annual_any_strong_response_area_pct=("annual_mean_any_impacted_area_pct", "max"),
        max_annual_dry_area_pct=("annual_mean_dry_area_pct", "max"),
    )
    .reset_index()
)

if not validation.empty:
    keep_cols = [
        "coast_region",
        "spatial_period_mean_WUE_T_pct_change",
        "spatial_period_mean_abs_WUE_T_pct_change",
        "spatial_pct_pixels_strong_decrease",
        "spatial_pct_pixels_strong_increase",
        "spatial_dry_pct_pixels",
        "spatial_dry_mean_WUE_T_pct_change",
        "spatial_dry_mean_abs_WUE_T_pct_change",
        "spatial_dry_pct_pixels_strong_decrease",
        "spatial_dry_pct_pixels_strong_increase",
        "Q3_n_sites",
        "Q3_SPEI3_mean_slope",
        "Q3_pct_negative_site_slopes",
    ]
    keep_cols = [c for c in keep_cols if c in validation.columns]
    impact_summary = impact_summary.merge(validation[keep_cols], on="coast_region", how="left")

impact_summary = order_coasts(impact_summary)
print(impact_summary.round(3).to_string(index=False))

# ============================================================
# 3. Q3 FITTED CURVE RANGE
# ============================================================
print_section("3. Q3-ALIGNED SPEI-3 FITTED CURVE RANGE BY COAST")

curve_summary = pd.DataFrame()

if not q4_curve.empty:
    curve_sp3 = q4_curve[
        (q4_curve["SPEI_timescale"] == "SPEI_3") &
        (q4_curve["water_class"] == "Upland")
    ].copy()

    curve_summary = (
        curve_sp3
        .groupby("coast_region")
        .agg(
            curve_min_SPEI=("SPEI_value", "min"),
            curve_max_SPEI=("SPEI_value", "max"),
            min_predicted_pct_change=("predicted_pct_change", "min"),
            max_predicted_pct_change=("predicted_pct_change", "max"),
            max_abs_predicted_pct_change=("predicted_pct_change", lambda x: np.nanmax(np.abs(x))),
            mean_abs_predicted_pct_change=("predicted_pct_change", lambda x: np.nanmean(np.abs(x))),
        )
        .reset_index()
    )

    curve_summary["crosses_5pct_threshold"] = curve_summary["max_abs_predicted_pct_change"] >= 5
    curve_summary = order_coasts(curve_summary)

    print(curve_summary.round(4).to_string(index=False))
else:
    print("Skipped: Q4 curve file missing.")

# ============================================================
# 4. Q3 SMOOTH SUPPORT / SIGNIFICANCE
# ============================================================
print_section("4. Q3 SPEI-3 SMOOTH TERM SUPPORT")

smooth_summary = pd.DataFrame()

if not q3_smooth.empty:
    smooth = q3_smooth.copy()

    term_col = first_existing_col(smooth, ["smooth_term", "term"])
    p_col = first_existing_col(smooth, ["p-value", "p.value", "p_value", "p"])
    edf_col = first_existing_col(smooth, ["edf", "EDF"])
    f_col = first_existing_col(smooth, ["F", "F_value", "F.value"])

    if term_col is None:
        print("Could not identify smooth term column.")
    else:
        smooth_sp3 = smooth[smooth[term_col].astype(str).str.contains("SPEI_3", case=False, na=False)].copy()

        def coast_from_term(term):
            term = str(term)
            for coast in COAST_ORDER:
                if coast in term:
                    return coast
                if coast.replace(" ", ".") in term:
                    return coast
            return np.nan

        smooth_sp3["coast_region"] = smooth_sp3[term_col].apply(coast_from_term)

        cols = ["coast_region", term_col]
        for c in [edf_col, f_col, p_col]:
            if c is not None:
                cols.append(c)

        smooth_summary = smooth_sp3[cols].copy()

        if p_col is not None:
            smooth_summary["significant_p_lt_0.05"] = smooth_summary[p_col] < 0.05
        if edf_col is not None:
            smooth_summary["smooth_shrunk_near_zero"] = smooth_summary[edf_col] < 0.1

        smooth_summary = order_coasts(smooth_summary)
        print(smooth_summary.round(5).to_string(index=False))
else:
    print("Skipped: Q3 smooth file missing.")

# ============================================================
# 5. SITE COUNTS
# ============================================================
print_section("5. SITE COUNTS USED IN Q4 MONTHLY UPLAND EDI SCORES")

site_count_summary = pd.DataFrame()

if not q4_monthly.empty:
    site_count_summary = (
        q4_monthly
        .groupby("coast_region")
        .agg(
            n_sites=("site_name", "nunique"),
            n_site_months=("site_name", "count"),
            first_year=("Year", "min"),
            last_year=("Year", "max"),
            mean_SPEI3_site_month=("SPEI_3", "mean"),
            mean_abs_EDI_Q3_pct_change=("EDI_Q3_magnitude", "mean"),
            median_abs_EDI_Q3_pct_change=("EDI_Q3_magnitude", "median"),
            pct_site_months_5pct_or_more=("EDI_Q3_magnitude", lambda x: (x >= 5).mean() * 100),
            pct_site_months_10pct_or_more=("EDI_Q3_magnitude", lambda x: (x >= 10).mean() * 100),
        )
        .reset_index()
    )

    site_count_summary = order_coasts(site_count_summary)
    print(site_count_summary.round(3).to_string(index=False))

    print("\nTotal unique upland sites in Q4 monthly scores:", q4_monthly["site_name"].nunique())
else:
    print("Skipped: Q4 monthly scores file missing.")

# ============================================================
# 6. Q3 SITE SENSITIVITY SUMMARY
# ============================================================
print_section("6. Q3 SITE SENSITIVITY SUMMARY FOR SPEI-3")

if not q3_site_sens.empty:
    sens = q3_site_sens.copy()

    if "SPEI_timescale" in sens.columns:
        sens = sens[sens["SPEI_timescale"] == "SPEI_3"].copy()

    sens = order_coasts(sens)
    print(sens.round(3).to_string(index=False))
else:
    print("Skipped: Q3 site sensitivity file missing.")

# ============================================================
# 7. MODEL COMPARISON
# ============================================================
print_section("7. Q3 MODEL COMPARISON")

if not q3_model_comp.empty:
    print(q3_model_comp.round(4).to_string(index=False))
else:
    print("Skipped: Q3 model comparison file missing.")

# ============================================================
# 8. SENTENCE-READY VALUES
# ============================================================
print_section("8. SENTENCE-READY VALUES FOR RESULTS SECTION")

combined = spei_summary[[
    "coast_region",
    "mean_SPEI3",
    "min_monthly_mean_SPEI3",
    "max_monthly_mean_SPEI3",
    "cumulative_SPEI3",
    "dry_month_pct",
    "wet_month_pct",
]].merge(
    impact_summary,
    on="coast_region",
    how="left"
)

if not curve_summary.empty:
    combined = combined.merge(
        curve_summary[[
            "coast_region",
            "min_predicted_pct_change",
            "max_predicted_pct_change",
            "max_abs_predicted_pct_change",
            "crosses_5pct_threshold",
        ]],
        on="coast_region",
        how="left"
    )

combined = order_coasts(combined)

for _, r in combined.iterrows():
    coast = r["coast_region"]
    label = DISPLAY.get(coast, coast)

    print(f"\n{label}")
    print("-" * len(label))

    print(
        f"SPEI-3 exposure: mean SPEI-3 = {fmt(r.get('mean_SPEI3'))}; "
        f"dry monthly regional means <= -1 = {fmt(r.get('dry_month_pct'))}%; "
        f"wet monthly regional means >= 1 = {fmt(r.get('wet_month_pct'))}%; "
        f"cumulative SPEI-3 = {fmt(r.get('cumulative_SPEI3'))}; "
        f"monthly mean range = {fmt(r.get('min_monthly_mean_SPEI3'))} to {fmt(r.get('max_monthly_mean_SPEI3'))}."
    )

    print(
        f"Spatial WUE_T response: mean absolute predicted WUE_T change = "
        f"{fmt(r.get('spatial_period_mean_abs_WUE_T_pct_change'))}%; "
        f"negative event area = {fmt(r.get('mean_annual_negative_event_area_pct'))}%; "
        f"positive event area = {fmt(r.get('mean_annual_positive_event_area_pct'))}%; "
        f"any strong response area = {fmt(r.get('mean_annual_any_strong_response_area_pct'))}%."
    )

    print(
        f"Dry pixels/months: mean dry area = {fmt(r.get('mean_annual_dry_area_pct'))}%; "
        f"dry mean WUE_T change = {fmt(r.get('spatial_dry_mean_WUE_T_pct_change'))}%; "
        f"dry strong decrease = {fmt(r.get('spatial_dry_pct_pixels_strong_decrease'))}%; "
        f"dry strong increase = {fmt(r.get('spatial_dry_pct_pixels_strong_increase'))}%."
    )

    if "max_abs_predicted_pct_change" in r.index:
        print(
            f"Q3 fitted SPEI-3 curve: predicted WUE_T change range = "
            f"{fmt(r.get('min_predicted_pct_change'))}% to {fmt(r.get('max_predicted_pct_change'))}%; "
            f"max absolute fitted change = {fmt(r.get('max_abs_predicted_pct_change'))}%; "
            f"crosses 5% threshold = {r.get('crosses_5pct_threshold')}."
        )

# ============================================================
# 9. SIMPLE INTERPRETATION FLAGS
# ============================================================
print_section("9. SIMPLE INTERPRETATION FLAGS")

for _, r in combined.iterrows():
    coast = r["coast_region"]
    label = DISPLAY.get(coast, coast)

    mean_spei = r.get("mean_SPEI3", np.nan)
    any_resp = r.get("mean_annual_any_strong_response_area_pct", np.nan)
    dry_change = r.get("spatial_dry_mean_WUE_T_pct_change", np.nan)

    if pd.isna(mean_spei):
        hydro = "unknown hydroclimate"
    elif mean_spei <= -0.10:
        hydro = "dry tendency"
    elif mean_spei >= 0.10:
        hydro = "wet tendency"
    else:
        hydro = "near-normal mean SPEI-3"

    if pd.isna(any_resp):
        response = "unknown response"
    elif any_resp >= 25:
        response = "strong threshold-level spatial response"
    elif any_resp >= 5:
        response = "moderate threshold-level spatial response"
    else:
        response = "weak/no threshold-level spatial response"

    if pd.isna(dry_change):
        dryflag = "unknown dry-month response"
    elif dry_change >= 5:
        dryflag = "dry-month WUE_T increase"
    elif dry_change <= -5:
        dryflag = "dry-month WUE_T decrease"
    else:
        dryflag = "dry-month response below 5% threshold"

    print(f"{label}: {hydro}; {response}; {dryflag}")

print("\nDONE. Nothing was written to the output directory.")