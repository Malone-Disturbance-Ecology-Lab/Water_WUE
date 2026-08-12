import os
import pandas as pd
import numpy as np

# ============================================================
# PATHS
# ============================================================

Q3 = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs"

wide_file = os.path.join(
    Q3,
    "Q3_WUE_T_SPEI_model_base_wide.csv"
)

long_file = os.path.join(
    Q3,
    "Q3_WUE_T_SPEI_model_data_long.csv"
)

# ============================================================
# LOAD
# ============================================================

wide = pd.read_csv(wide_file)
long = pd.read_csv(long_file)

print("=" * 100)
print("COASTAL SITE REPRESENTATION AND TEMPORAL COVERAGE")
print("=" * 100)

required = [
    "site_name", "Year", "month",
    "water_class", "coast_region"
]

missing = [c for c in required if c not in wide.columns]

if missing:
    raise ValueError(f"Missing columns: {missing}")

# One unique site-month only
d = (
    wide[required]
    .drop_duplicates()
    .copy()
)

d["date"] = pd.to_datetime(
    dict(
        year=d["Year"].astype(int),
        month=d["month"].astype(int),
        day=15
    )
)

# ============================================================
# 1. SITE-LEVEL RECORD LENGTH
# ============================================================

site_summary = (
    d.groupby(
        ["coast_region", "water_class", "site_name"],
        observed=True
    )
    .agg(
        n_months=("date", "nunique"),
        n_years=("Year", "nunique"),
        first_date=("date", "min"),
        last_date=("date", "max")
    )
    .reset_index()
)

site_summary["record_span_years"] = (
    (
        site_summary["last_date"] -
        site_summary["first_date"]
    ).dt.days / 365.25
)

# ============================================================
# 2. NUMBER OF SITES BY COAST AND ECOSYSTEM
# ============================================================

print("\n" + "=" * 100)
print("1. SITE COUNTS BY COAST AND ECOSYSTEM")
print("=" * 100)

counts = (
    site_summary.groupby(
        ["coast_region", "water_class"],
        observed=True
    )
    .agg(
        n_sites=("site_name", "nunique")
    )
    .reset_index()
)

count_pivot = (
    counts.pivot(
        index="coast_region",
        columns="water_class",
        values="n_sites"
    )
    .fillna(0)
)

for eco in ["Upland", "Freshwater", "Saline"]:
    if eco not in count_pivot.columns:
        count_pivot[eco] = 0

count_pivot["Wetland_total"] = (
    count_pivot["Freshwater"] +
    count_pivot["Saline"]
)

count_pivot["All_sites"] = (
    count_pivot["Upland"] +
    count_pivot["Freshwater"] +
    count_pivot["Saline"]
)

print(
    count_pivot[
        [
            "Upland",
            "Freshwater",
            "Saline",
            "Wetland_total",
            "All_sites"
        ]
    ]
    .astype(int)
    .sort_values("All_sites", ascending=False)
    .to_string()
)

# ============================================================
# 3. TEMPORAL COVERAGE BY COAST
# ============================================================

print("\n" + "=" * 100)
print("2. TEMPORAL COVERAGE BY COAST")
print("=" * 100)

coast_time = (
    site_summary.groupby(
        "coast_region",
        observed=True
    )
    .agg(
        n_sites=("site_name", "nunique"),
        mean_months=("n_months", "mean"),
        median_months=("n_months", "median"),
        mean_years_observed=("n_years", "mean"),
        median_years_observed=("n_years", "median"),
        mean_record_span_years=("record_span_years", "mean"),
        median_record_span_years=("record_span_years", "median")
    )
    .reset_index()
)

print(
    coast_time
    .round(2)
    .sort_values(
        "median_months",
        ascending=False
    )
    .to_string(index=False)
)

# ============================================================
# 4. WETLAND TEMPORAL COVERAGE BY COAST
# ============================================================

print("\n" + "=" * 100)
print("3. WETLAND-ONLY COVERAGE BY COAST")
print("=" * 100)

wet = site_summary[
    site_summary["water_class"].isin(
        ["Freshwater", "Saline"]
    )
].copy()

wet_time = (
    wet.groupby(
        "coast_region",
        observed=True
    )
    .agg(
        n_wetland_sites=("site_name", "nunique"),
        mean_months=("n_months", "mean"),
        median_months=("n_months", "median"),
        mean_years_observed=("n_years", "mean"),
        median_years_observed=("n_years", "median"),
        median_record_span_years=(
            "record_span_years",
            "median"
        )
    )
    .reset_index()
)

print(
    wet_time
    .round(2)
    .sort_values(
        "n_wetland_sites",
        ascending=False
    )
    .to_string(index=False)
)

# ============================================================
# 5. SPEI-48 SUPPORT
#
# Directly checks data supporting the longest-timescale model.
# ============================================================

print("\n" + "=" * 100)
print("4. SPEI-48 SUPPORT BY COAST")
print("=" * 100)

sp48 = long[
    long["SPEI_timescale"] == "SPEI_48"
].copy()

sp48_summary = (
    sp48.groupby(
        "coast_region",
        observed=True
    )
    .agg(
        n_sites=("site_name", "nunique"),
        n_site_months=("WUE_T", "size"),
        first_year=("Year", "min"),
        last_year=("Year", "max")
    )
    .reset_index()
)

print(
    sp48_summary
    .sort_values(
        "n_site_months",
        ascending=False
    )
    .to_string(index=False)
)

print("\n" + "=" * 100)
print("5. SPEI-48 SUPPORT BY COAST × ECOSYSTEM")
print("=" * 100)

sp48_eco = (
    sp48.groupby(
        ["coast_region", "water_class"],
        observed=True
    )
    .agg(
        n_sites=("site_name", "nunique"),
        n_site_months=("WUE_T", "size")
    )
    .reset_index()
)

print(
    sp48_eco
    .sort_values(
        ["coast_region", "water_class"]
    )
    .to_string(index=False)
)

# ============================================================
# 6. SHORTEST RECORDS
# ============================================================

print("\n" + "=" * 100)
print("6. SHORTEST SITE RECORDS")
print("=" * 100)

print(
    site_summary[
        [
            "coast_region",
            "water_class",
            "site_name",
            "n_months",
            "n_years",
            "first_date",
            "last_date"
        ]
    ]
    .sort_values(
        ["n_months", "n_years"]
    )
    .head(20)
    .to_string(index=False)
)

print("\nDONE. Nothing was saved.")