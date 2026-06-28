# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:46:15 2026

@author: ammar
"""

# =============================================================================
# CLIMATE–BIOME Q2 DATA PREPARATION
# Creates climate–biome grouped dataset for Q2-style analyses
# =============================================================================

import pandas as pd
import numpy as np
import os

# =============================================================================
# FILE PATHS
# =============================================================================

BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"

# Q2 site-level percent-change file
Q2_FILE = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\results\SPEI_analysis_results"
    r"\SPEI_site_level_details_FINAL.csv"
)

# Climate–biome grouped dataset from previous workflow
CLIMATE_FILE = (
    BASE_DIR +
    r"\ClimateBiome_WUE_WUET_FinalDataset.csv"
)

# Output
OUTPUT_FILE = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\results\SPEI_analysis_results"
    r"\Q2_ClimateBiome_Dataset.csv"
)

# =============================================================================
# LOAD DATA
# =============================================================================

print("="*80)
print("Q2 CLIMATE–BIOME DATASET PREPARATION")
print("="*80)

q2_df = pd.read_csv(Q2_FILE)
climate_df = pd.read_csv(CLIMATE_FILE)

print(f"\nLoaded Q2 dataset: {len(q2_df)} rows")
print(f"Loaded climate–biome dataset: {len(climate_df)} rows")

# =============================================================================
# KEEP ONLY WUE + WUE_T
# =============================================================================

metrics_keep = ['WUE', 'WUE_tra']

q2_df = q2_df[
    q2_df['WUE_Metric'].isin(metrics_keep)
].copy()

print(f"\nRows after metric filter: {len(q2_df)}")

# =============================================================================
# EXTRACT SITE → FINAL_GROUP MAPPING
# =============================================================================

group_map = (
    climate_df[
        ['site_name', 'Final_Group']
    ]
    .drop_duplicates()
    .rename(columns={'site_name': 'Site'})
)

print(f"\nClimate–biome mapped sites: {len(group_map)}")

# =============================================================================
# MERGE FINAL GROUPS INTO Q2 DATA
# =============================================================================

q2_df = q2_df.merge(
    group_map,
    on='Site',
    how='left'
)

# =============================================================================
# CHECK MISSING GROUPS
# =============================================================================

missing_sites = (
    q2_df[
        q2_df['Final_Group'].isna()
    ]['Site']
    .drop_duplicates()
)

print("\n" + "="*80)
print("CHECKING UNASSIGNED SITES")
print("="*80)

if len(missing_sites) > 0:

    print(f"\n⚠️ Sites without climate–biome group: {len(missing_sites)}")

    for s in sorted(missing_sites):
        print(f"  - {s}")

else:
    print("\n✓ All sites assigned successfully")

# =============================================================================
# REMOVE UNASSIGNED SITES
# =============================================================================

q2_df = q2_df[
    q2_df['Final_Group'].notna()
].copy()

# =============================================================================
# STRICT DOUBLE INTERSECTION
# SAME LOGIC AS Q2 FIGURES
# =============================================================================

print("\n" + "="*80)
print("APPLYING STRICT DOUBLE INTERSECTION")
print("="*80)

filtered_results = []

timescales = ['SPEI_6', 'SPEI_48']

conditions = [
    ('PASS C', 'Dry (all)'),
    ('PASS C', 'Wet (all)'),
    ('PASS B', 'Severe Dry'),
    ('PASS B', 'Severe Wet')
]

for pass_name, condition in conditions:

    for timescale in timescales:

        subset = q2_df[
            (q2_df['Pass'] == pass_name) &
            (q2_df['Condition'] == condition) &
            (q2_df['SPEI_Timescale'] == timescale)
        ].copy()

        # strict shared-site intersection
        wue_sites = set(
            subset[
                subset['WUE_Metric'] == 'WUE'
            ]['Site'].unique()
        )

        tra_sites = set(
            subset[
                subset['WUE_Metric'] == 'WUE_tra'
            ]['Site'].unique()
        )

        shared_sites = wue_sites.intersection(tra_sites)

        sub_filtered = subset[
            subset['Site'].isin(shared_sites)
        ].copy()

        filtered_results.append(sub_filtered)

        print(f"\n{pass_name} | {condition} | {timescale}")

        print(f"  Shared sites: {len(shared_sites)}")

# =============================================================================
# COMBINE FILTERED RESULTS
# =============================================================================

final_df = pd.concat(
    filtered_results,
    ignore_index=True
)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("\n" + "="*80)
print("SUPPORT SUMMARY BY CLIMATE–BIOME GROUP")
print("="*80)

summary = (
    final_df[
        ['Site', 'Final_Group', 'Condition', 'SPEI_Timescale']
    ]
    .drop_duplicates()
    .groupby(
        ['Final_Group', 'Condition', 'SPEI_Timescale']
    )
    .size()
    .reset_index(name='N_sites')
)

print(summary.to_string(index=False))

# =============================================================================
# METRIC BALANCE CHECK
# =============================================================================

print("\n" + "="*80)
print("METRIC BALANCE CHECK")
print("="*80)

for metric in metrics_keep:

    sub = final_df[
        final_df['WUE_Metric'] == metric
    ]

    n_sites = len(sub['Site'].unique())

    print(f"\n{metric}: {n_sites} sites")

# =============================================================================
# SAVE DATASET
# =============================================================================

final_df.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\n" + "="*80)
print("DATASET SAVED")
print("="*80)

print(f"\nSaved file:")
print(OUTPUT_FILE)

print("\n✓ Q2 climate–biome dataset preparation complete")