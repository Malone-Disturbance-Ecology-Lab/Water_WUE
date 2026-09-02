# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:48:26 2026

@author: ammar
"""

import os
import pandas as pd
import numpy as np

# ============================================================
# PATHS
# ============================================================
MONTHLY_DATA = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

GAPFILLED_DIR = r"M:\Research\WUE_CUE\ameri_data\ameri_fill"

# ============================================================
# LOAD FINAL MONTHLY DATA
# Apply EXACTLY the same filters used in Study_area_august.py
# ============================================================
df = pd.read_csv(MONTHLY_DATA)

for col in df.select_dtypes(include="object").columns:
    df[col] = df[col].where(
        df[col].isna(),
        df[col].astype(str).str.strip()
    )
    df[col] = df[col].replace({
        "": np.nan,
        "nan": np.nan,
        "NaN": np.nan
    })

required = [
    'site_name',
    'Year',
    'month',
    'water_class',
    'lat',
    'long',
    'Trans_ratio',
    'WUE_tra'
]

study = df.dropna(subset=required).copy()

study = study[
    study['water_class'].isin(
        ['Upland', 'Freshwater', 'Saline']
    )
].copy()

study = study[
    np.isfinite(study['lat']) &
    np.isfinite(study['long']) &
    np.isfinite(study['Trans_ratio']) &
    (study['Trans_ratio'] >= 0) &
    (study['Trans_ratio'] <= 1) &
    np.isfinite(study['WUE_tra'])
].copy()

# ============================================================
# FINAL STUDY SITE LIST
# ============================================================
study_sites = sorted(study['site_name'].unique())

print("=" * 100)
print("FINAL ANALYTICAL DATASET USED BY STUDY-AREA CODE")
print("=" * 100)

print(f"\nFinal study sites: {len(study_sites)}")
print(", ".join(study_sites))

# ============================================================
# SITE-YEAR / SITE-MONTH COVERAGE IN ACTUAL ANALYSIS
# ============================================================
coverage = (
    study.groupby('site_name')
    .agg(
        first_year=('Year', 'min'),
        last_year=('Year', 'max'),
        n_site_months=('Year', 'size'),
        n_years=('Year', 'nunique'),
        ecosystem=('water_class', 'first')
    )
    .reset_index()
)

# Number of actual contributing months by year
site_year = (
    study.groupby(['site_name', 'Year'])
    .size()
    .reset_index(name='n_months')
)

# identify whether years are continuous
def summarize_years(x):
    years = sorted(x.astype(int).unique())
    return ", ".join(map(str, years))

year_lists = (
    study.groupby('site_name')['Year']
    .apply(summarize_years)
    .reset_index(name='years_present')
)

coverage = coverage.merge(
    year_lists,
    on='site_name',
    how='left'
)

print("\n" + "=" * 100)
print("ACTUAL YEARS AND MONTHS CONTRIBUTED TO ANALYSIS")
print("=" * 100)

print(
    coverage.sort_values(
        ['first_year', 'site_name']
    ).to_string(index=False)
)

# ============================================================
# DISTRIBUTION OF SITE CONTRIBUTION
# ============================================================
print("\n" + "=" * 100)
print("SITE CONTRIBUTION SUMMARY")
print("=" * 100)

print(f"Study period represented: "
      f"{int(study['Year'].min())}–{int(study['Year'].max())}")

print(f"Total study site-months: {len(study):,}")

print("\nYears contributed per site:")
print(
    coverage['n_years']
    .describe()
    .round(2)
    .to_string()
)

print("\nSite-months contributed per site:")
print(
    coverage['n_site_months']
    .describe()
    .round(2)
    .to_string()
)

# ============================================================
# NUMBER OF SITES CONTRIBUTING EACH YEAR
# ============================================================
year_summary = (
    study.groupby('Year')
    .agg(
        n_sites=('site_name', 'nunique'),
        n_site_months=('site_name', 'size')
    )
    .reset_index()
)

print("\n" + "=" * 100)
print("NUMBER OF STUDY SITES CONTRIBUTING BY YEAR")
print("=" * 100)

print(year_summary.to_string(index=False))

# ============================================================
# COMPARE FINAL STUDY SITES AGAINST GAP-FILLED DIRECTORY
# ============================================================
gapfill_sites = []

for f in os.listdir(GAPFILLED_DIR):
    if f.endswith("_fill.csv"):
        site = f.replace("_fill.csv", "")
        gapfill_sites.append(site)

gapfill_sites = sorted(set(gapfill_sites))

study_set = set(study_sites)
gap_set = set(gapfill_sites)

print("\n" + "=" * 100)
print("STUDY SITES VS ALL GAP-FILLED SITES")
print("=" * 100)

print(f"Sites processed through gap filling: {len(gapfill_sites)}")
print(f"Sites retained in final study dataset: {len(study_sites)}")

excluded = sorted(gap_set - study_set)

print(f"\nGap-filled sites NOT represented in final study dataset: {len(excluded)}")
print(", ".join(excluded) if excluded else "None")

# ============================================================
# OPTIONAL: YEAR-BY-YEAR MONTH CONTRIBUTION PER SITE
# ============================================================
print("\n" + "=" * 100)
print("SITE-YEAR MONTH COUNTS")
print("Each number = months from that site/year surviving into analysis")
print("=" * 100)

pivot = site_year.pivot(
    index='site_name',
    columns='Year',
    values='n_months'
).fillna(0).astype(int)

print(pivot.to_string())

print("\nDone. Nothing was saved.")