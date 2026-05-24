# -*- coding: utf-8 -*-
"""
Climate × Broad Biome Group Analysis
WUE_ET vs WUE_T only

IMPORTANT:
Filtering is STILL based on ALL THREE metrics:
    - WUE
    - WUE_eva
    - WUE_tra

This preserves consistency with Figures 1–2.

Final analysis uses only:
    - WUE
    - WUE_tra
"""

import pandas as pd
import numpy as np
from scipy.stats import (
    kruskal,
    mannwhitneyu,
    wilcoxon,
    friedmanchisquare
)
from statsmodels.stats.multitest import multipletests

# =============================================================================
# FILE PATHS
# =============================================================================

BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"

FILE = BASE_DIR + r"\wue_site_level_summary_SPEI_1.csv"
MONTHLY_FILE = BASE_DIR + r"\monthly_data_after_outlier_removal.csv"

# =============================================================================
# LOAD DATA
# =============================================================================

print("="*80)
print("CLIMATE × BROAD BIOME ANALYSIS")
print("WUE_ET and WUE_T only")
print("="*80)

df = pd.read_csv(FILE)
monthly_df = pd.read_csv(MONTHLY_FILE)

print(f"Loaded site-level file: {len(df)} rows")
print(f"Loaded monthly file: {len(monthly_df)} rows")

# =============================================================================
# STEP 1: NN FILTER + STRICT TRIPLE INTERSECTION
# =============================================================================

print("\n" + "="*80)
print("STEP 1: NN FILTER + STRICT TRIPLE INTERSECTION")
print("="*80)

# Site-level
nn_df = df[df['SPEI_Class'] == 'NN'].copy()

# Monthly
monthly_nn = monthly_df[
    monthly_df['SPEI_1_Cat'] == 'NN'
].copy()

print(f"\nNN site-level rows: {len(nn_df)}")
print(f"NN monthly rows: {len(monthly_nn)}")

print("\nApplying strict triple-intersection filtering...")

sites_wue = set(
    nn_df[
        nn_df['WUE_Metric'] == 'WUE'
    ]['site_name'].unique()
)

sites_eva = set(
    nn_df[
        nn_df['WUE_Metric'] == 'WUE_eva'
    ]['site_name'].unique()
)

sites_tra = set(
    nn_df[
        nn_df['WUE_Metric'] == 'WUE_tra'
    ]['site_name'].unique()
)

shared_sites = (
    sites_wue
    .intersection(sites_eva)
    .intersection(sites_tra)
)

print(f"\nSites with WUE: {len(sites_wue)}")
print(f"Sites with WUE_E: {len(sites_eva)}")
print(f"Sites with WUE_T: {len(sites_tra)}")

print(f"\nSTRICT triple intersection: {len(shared_sites)} sites")

nn_df = nn_df[
    nn_df['site_name'].isin(shared_sites)
].copy()

# Keep ONLY WUE + WUE_T for analysis
nn_df = nn_df[
    nn_df['WUE_Metric'].isin(['WUE', 'WUE_tra'])
].copy()

print("\nFinal metrics retained:")
print(sorted(nn_df['WUE_Metric'].unique()))

print(f"\nRows after filtering: {len(nn_df)}")
print(f"Unique shared sites retained: {nn_df['site_name'].nunique()}")

print("\nVerification by metric:")

for metric in ['WUE', 'WUE_tra']:

    n_sites = (
        nn_df[
            nn_df['WUE_Metric'] == metric
        ]['site_name']
        .nunique()
    )

    print(f"  {metric}: {n_sites} sites")

print("\n✓ STEP 1 COMPLETE")

# =============================================================================
# STEP 2: GLOBAL NN FILTER + VALID T:ET FILTER
# =============================================================================

print("\n" + "="*80)
print("STEP 2: GLOBAL NN FILTER + VALID T:ET FILTER")
print("="*80)

print("\nApplying global NN observation filter (>=3 observations)...")

nn_obs_counts = (
    monthly_nn
    .groupby('site_name')
    .size()
    .reset_index(name='NN_obs')
)

valid_nn_sites = set(
    nn_obs_counts[
        nn_obs_counts['NN_obs'] >= 3
    ]['site_name']
)

print(f"Sites with >=3 NN observations: {len(valid_nn_sites)}")

nn_df = nn_df[
    nn_df['site_name'].isin(valid_nn_sites)
].copy()

print(f"Sites retained after NN filter: {nn_df['site_name'].nunique()}")

print("\nApplying valid T:ET filtering...")

monthly_valid_tet = monthly_nn.copy()

monthly_valid_tet = monthly_valid_tet[
    monthly_valid_tet['Trans_ratio'].notna()
]

monthly_valid_tet = monthly_valid_tet[
    monthly_valid_tet['Trans_ratio'] > 0
]

valid_tet_sites = set(
    monthly_valid_tet['site_name'].unique()
)

print(f"Sites with valid T:ET values: {len(valid_tet_sites)}")

nn_df = nn_df[
    nn_df['site_name'].isin(valid_tet_sites)
].copy()

print("\nFinal dataset after all filtering:")

print(f"  Total rows: {len(nn_df)}")
print(f"  Unique sites: {nn_df['site_name'].nunique()}")

print("\nVerification by metric:")

for metric in ['WUE', 'WUE_tra']:

    sub = nn_df[
        nn_df['WUE_Metric'] == metric
    ]

    print(
        f"  {metric}: "
        f"{sub['site_name'].nunique()} sites, "
        f"{len(sub)} rows"
    )

site_support = (
    monthly_valid_tet[
        monthly_valid_tet['site_name'].isin(
            nn_df['site_name'].unique()
        )
    ]
    .groupby('site_name')
    .size()
)

print("\nNN observation support:")

print(f"  Total retained NN observations: {site_support.sum()}")
print(f"  Median observations per site: {site_support.median():.1f}")
print(f"  Min observations per site: {site_support.min()}")
print(f"  Max observations per site: {site_support.max()}")

print("\n✓ STEP 2 COMPLETE")

# =============================================================================
# STEP 3: CREATE BROAD BIOME GROUPS
# =============================================================================

print("\n" + "="*80)
print("STEP 3: CREATE BROAD BIOME GROUPS")
print("="*80)

biome_map = {

    'WET': 'Wetlands',

    'ENF': 'Forests',
    'DBF': 'Forests',
    'MF': 'Forests',

    'GRA': 'Drylands',
    'CSH': 'Drylands',
    'OSH': 'Drylands',

    'CRO': 'Croplands',

    'BSV': 'Sparse/Barren'
}

print("\nOriginal biome classes:\n")

orig_counts = (
    nn_df[
        ['site_name', 'biome']
    ]
    .drop_duplicates()['biome']
    .value_counts()
)

for biome, n in orig_counts.items():

    print(f"{biome}: {n} sites")

# APPLY BIOME MAPPING
nn_df['Broad_Biome'] = (
    nn_df['biome']
    .map(biome_map)
)

print("\nChecking for unmapped biome classes...")

unmapped = (
    nn_df[
        nn_df['Broad_Biome'].isna()
    ]['biome']
    .dropna()
    .unique()
)

if len(unmapped) > 0:

    print("\n⚠️ Unmapped biome classes detected:")

    for b in sorted(unmapped):

        print(f"  - {b}")

else:

    print("\n✓ All biome classes mapped successfully")

print("\n" + "-"*70)
print("BROAD BIOME SUMMARY")
print("-"*70)

broad_counts = (
    nn_df[
        ['site_name', 'Broad_Biome']
    ]
    .drop_duplicates()['Broad_Biome']
    .value_counts()
)

for biome, n in broad_counts.items():

    print(f"{biome}: {n} sites")

print("\n" + "-"*70)
print("CLIMATE SUMMARY")
print("-"*70)

climate_counts = (
    nn_df[
        ['site_name', 'climate']
    ]
    .drop_duplicates()['climate']
    .value_counts()
)

for climate, n in climate_counts.items():

    print(f"{climate}: {n} sites")

print("\n" + "-"*70)
print("CLIMATE × BROAD BIOME CROSS-TAB")
print("-"*70)

cross_tab = pd.crosstab(

    nn_df[
        ['site_name', 'climate', 'Broad_Biome']
    ].drop_duplicates()['climate'],

    nn_df[
        ['site_name', 'climate', 'Broad_Biome']
    ].drop_duplicates()['Broad_Biome']
)

print(cross_tab)

print("\n✓ STEP 3 COMPLETE")

# =============================================================================
# STEP 4+: KEEP EVERYTHING BELOW EXACTLY THE SAME
# =============================================================================

# DO NOT CHANGE:
# - STEP 4
# - STEP 5
# - STEP 6
# - OUTPUT FILE NAME
# - SAVED CSV NAME
# - FIGURE WORKFLOW

# Continue using your existing downstream code unchanged.



# =============================================================================
# STEP 4: CREATE FINAL CLIMATE × BIOME GROUPS
# =============================================================================

print("\n" + "="*80)
print("STEP 4: CREATE FINAL CLIMATE × BIOME GROUPS")
print("="*80)

# ---------------------------------------------------------------------
# CLEAN CLIMATE STRINGS
# Removes accidental trailing parentheses/spaces
# ---------------------------------------------------------------------

nn_df['climate_clean'] = (
    nn_df['climate']
    .astype(str)
    .str.strip()
    .str.replace(")", "", regex=False)
)

# ---------------------------------------------------------------------
# CREATE FINAL ECOLOGICAL GROUPS
# ---------------------------------------------------------------------

def assign_final_group(row):

    climate = row['climate_clean']
    biome = row['Broad_Biome']

    # ================================================================
    # HUMID COASTAL WETLANDS
    # ================================================================
    if (
        (
            'Humid Subtropical' in climate
            or 'Humid Continental' in climate
        )
        and biome == 'Wetlands'
    ):
        return 'Humid Coastal Wetlands'

    # ================================================================
    # HUMID COASTAL FORESTS
    # ================================================================
    if (
        (
            'Humid Subtropical' in climate
            or 'Humid Continental' in climate
        )
        and biome == 'Forests'
    ):
        return 'Humid Coastal Forests'

    # ================================================================
    # MEDITERRANEAN WETLANDS
    # ================================================================
    if (
        'Mediterranean' in climate
        and biome == 'Wetlands'
    ):
        return 'Mediterranean Wetlands'

    # ================================================================
    # MEDITERRANEAN DRYLANDS
    # ================================================================
    if (
        'Mediterranean' in climate
        and biome == 'Drylands'
    ):
        return 'Mediterranean Drylands'

    # ================================================================
    # COLD COASTAL SYSTEMS
    # ================================================================
    if (
        (
            'Tundra' in climate
            or 'Marine West Coast' in climate
            or 'Dry Continental' in climate
        )
    ):
        return 'Cold Coastal Systems'

    # ================================================================
    # OTHERWISE
    # ================================================================
    return np.nan

# Apply grouping
nn_df['Final_Group'] = nn_df.apply(
    assign_final_group,
    axis=1
)

# ---------------------------------------------------------------------
# CHECK UNASSIGNED SITES
# ---------------------------------------------------------------------

print("\nChecking for unassigned sites...")

unassigned = (
    nn_df[
        nn_df['Final_Group'].isna()
    ][
        ['site_name', 'climate_clean', 'Broad_Biome']
    ]
    .drop_duplicates()
)

if len(unassigned) > 0:

    print(f"\n⚠️ Unassigned sites: {len(unassigned)}\n")

    print(unassigned.sort_values('climate_clean'))

else:
    print("\n✓ All sites assigned successfully")

# ---------------------------------------------------------------------
# REMOVE UNASSIGNED
# ---------------------------------------------------------------------

nn_df = nn_df[
    nn_df['Final_Group'].notna()
].copy()

# ---------------------------------------------------------------------
# FINAL GROUP SUMMARY
# ---------------------------------------------------------------------

print("\n" + "-"*70)
print("FINAL GROUP SUMMARY")
print("-"*70)

group_summary = (
    nn_df[
        ['site_name', 'Final_Group']
    ]
    .drop_duplicates()['Final_Group']
    .value_counts()
)

for grp, n in group_summary.items():

    print(f"{grp}: {n} sites")

# ---------------------------------------------------------------------
# GROUP × BIOME TABLE
# ---------------------------------------------------------------------

print("\n" + "-"*70)
print("FINAL GROUP × BROAD BIOME")
print("-"*70)

group_biome = pd.crosstab(

    nn_df[
        ['site_name', 'Final_Group', 'Broad_Biome']
    ].drop_duplicates()['Final_Group'],

    nn_df[
        ['site_name', 'Final_Group', 'Broad_Biome']
    ].drop_duplicates()['Broad_Biome']
)

print(group_biome)

# ---------------------------------------------------------------------
# GROUP × CLIMATE TABLE
# ---------------------------------------------------------------------

print("\n" + "-"*70)
print("FINAL GROUP × CLIMATE")
print("-"*70)

group_climate = pd.crosstab(

    nn_df[
        ['site_name', 'Final_Group', 'climate_clean']
    ].drop_duplicates()['Final_Group'],

    nn_df[
        ['site_name', 'Final_Group', 'climate_clean']
    ].drop_duplicates()['climate_clean']
)

print(group_climate)

# ---------------------------------------------------------------------
# VERIFY METRIC BALANCE
# ---------------------------------------------------------------------

print("\n" + "-"*70)
print("METRIC BALANCE CHECK")
print("-"*70)

for metric in ['WUE', 'WUE_tra']:

    sub = nn_df[
        nn_df['WUE_Metric'] == metric
    ]

    counts = (
        sub[
            ['site_name', 'Final_Group']
        ]
        .drop_duplicates()['Final_Group']
        .value_counts()
    )

    print(f"\n{metric}")

    for grp, n in counts.items():

        print(f"  {grp}: {n} sites")

print("\n✓ STEP 4 COMPLETE")



# =============================================================================
# STEP 5: DESCRIPTIVE STATISTICS BY FINAL GROUP
# =============================================================================

print("\n" + "="*80)
print("STEP 5: DESCRIPTIVE STATISTICS BY FINAL GROUP")
print("="*80)

for metric in ['WUE', 'WUE_tra']:

    print("\n" + "="*70)
    print(f"METRIC: {metric}")
    print("="*70)

    metric_df = nn_df[
        nn_df['WUE_Metric'] == metric
    ].copy()

    for grp in sorted(metric_df['Final_Group'].unique()):

        vals = metric_df[
            metric_df['Final_Group'] == grp
        ]['Median'].dropna().values

        n = len(vals)

        median = np.median(vals)
        q25 = np.percentile(vals, 25)
        q75 = np.percentile(vals, 75)

        mean = np.mean(vals)
        sd = np.std(vals, ddof=1)
        se = sd / np.sqrt(n)

        min_val = np.min(vals)
        max_val = np.max(vals)

        print(f"\n{grp}")

        print(f"  n = {n}")
        print(f"  Median = {median:.3f}")
        print(f"  IQR = {q25:.3f} - {q75:.3f}")
        print(f"  Mean ± SE = {mean:.3f} ± {se:.3f}")
        print(f"  SD = {sd:.3f}")
        print(f"  Range = {min_val:.3f} - {max_val:.3f}")

print("\n✓ STEP 5 COMPLETE")


# =============================================================================
# STEP 6: KRUSKAL-WALLIS + PAIRWISE TESTS
# =============================================================================

print("\n" + "="*80)
print("STEP 6: KRUSKAL-WALLIS + PAIRWISE TESTS")
print("="*80)

for metric in ['WUE', 'WUE_tra']:

    print("\n" + "="*70)
    print(f"METRIC: {metric}")
    print("="*70)

    metric_df = nn_df[
        nn_df['WUE_Metric'] == metric
    ].copy()

    groups = []
    group_names = []

    # ---------------------------------------------------------------
    # BUILD GROUPS
    # ---------------------------------------------------------------

    for grp in sorted(metric_df['Final_Group'].unique()):

        vals = metric_df[
            metric_df['Final_Group'] == grp
        ]['Median'].dropna().values

        groups.append(vals)
        group_names.append(grp)

        print(f"\n{grp}")
        print(f"  n = {len(vals)}")
        print(f"  Median = {np.median(vals):.3f}")

    # ---------------------------------------------------------------
    # KRUSKAL-WALLIS
    # ---------------------------------------------------------------

    h_stat, p_val = kruskal(*groups)

    print("\n" + "-"*60)
    print("KRUSKAL-WALLIS TEST")
    print("-"*60)

    print(f"H-statistic = {h_stat:.3f}")
    print(f"p-value = {p_val:.6f}")

    if p_val < 0.05:

        print("\n✓ Significant group differences detected")

        # -----------------------------------------------------------
        # PAIRWISE MANN-WHITNEY TESTS
        # -----------------------------------------------------------

        pairwise_results = []

        for i in range(len(group_names)):
            for j in range(i + 1, len(group_names)):

                g1 = group_names[i]
                g2 = group_names[j]

                vals1 = groups[i]
                vals2 = groups[j]

                u_stat, p_pair = mannwhitneyu(
                    vals1,
                    vals2,
                    alternative='two-sided'
                )

                pairwise_results.append({

                    'comparison': f"{g1} vs {g2}",
                    'u_stat': u_stat,
                    'p_raw': p_pair,
                    'median1': np.median(vals1),
                    'median2': np.median(vals2)
                })

        # -----------------------------------------------------------
        # BONFERRONI CORRECTION
        # -----------------------------------------------------------

        raw_p = [x['p_raw'] for x in pairwise_results]

        reject, p_corr, _, _ = multipletests(
            raw_p,
            method='bonferroni'
        )

        print("\n" + "-"*60)
        print("PAIRWISE MANN-WHITNEY TESTS")
        print("-"*60)

        for i, res in enumerate(pairwise_results):

            sig = (
                '***' if p_corr[i] < 0.001 else
                '**' if p_corr[i] < 0.01 else
                '*' if p_corr[i] < 0.05 else
                'ns'
            )

            print(
                f"\n{res['comparison']}"
            )

            print(
                f"  U = {res['u_stat']:.1f}"
            )

            print(
                f"  p_raw = {res['p_raw']:.4f}"
            )

            print(
                f"  p_corr = {p_corr[i]:.4f} ({sig})"
            )

            print(
                f"  Median difference = "
                f"{res['median1']:.3f} vs "
                f"{res['median2']:.3f}"
            )

    else:

        print("\nNo significant differences among groups")

print("\n✓ STEP 6 COMPLETE")

# =============================================================================
# SAVE FINAL FIGURE DATASET
# =============================================================================

OUTPUT_FILE = (
    BASE_DIR +
    r"\ClimateBiome_WUE_WUET_FinalDataset.csv"
)

nn_df.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nSaved final figure dataset:")
print(OUTPUT_FILE)


