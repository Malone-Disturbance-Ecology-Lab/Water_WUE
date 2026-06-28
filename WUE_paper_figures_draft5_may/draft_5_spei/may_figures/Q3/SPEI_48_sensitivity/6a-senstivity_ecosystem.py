# -*- coding: utf-8 -*-
"""
CHUNK 1: WUE_T Sensitivity by Ecosystem/Salinity Class - SPEI-48
Purpose: Load data, run statistical tests, prepare data for visualization

Includes:
- Full 4-group analysis (Upland, Freshwater, Brackish, Saline)
- ALL pairwise comparisons (Freshwater-Upland, Freshwater-Saline, Upland-Saline, etc.)
- Freshwater vs Saline (excluding Brackish)
- Freshwater vs Salt-affected (Brackish + Saline combined)
- Analysis WITH and WITHOUT Brackish included
- Sensitivity magnitude tests for all comparisons
"""

import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================
# CONFIGURATION
# ============================================

# Input files
SENSITIVITY_FILE = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\Q3_analysis\Q3_WUET_SPEI48_slopes.csv")
TET_FILE = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\T_ET_ratio\TET_gradient_site_NN_medians.csv")
OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april\diagnostics_SPEI48")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Ecosystem colors
ECOSYSTEM_COLORS = {
    'Upland': '#800080',      # Purple
    'Freshwater': '#0000FF',  # Blue
    'Brackish': '#008080',    # Teal
    'Saline': '#FFA500'       # Orange
}

ECOSYSTEM_ORDER = ['Upland', 'Freshwater', 'Brackish', 'Saline']

print("="*80)
print("CHUNK 1: Data Preparation and Statistical Analysis (SPEI-48)")
print("="*80)

# ============================================
# STEP 1: LOAD AND MERGE DATA
# ============================================

print("\n📂 Loading SPEI-48 sensitivity data...")
df_slopes = pd.read_csv(SENSITIVITY_FILE)
print(f"   Slopes file shape: {df_slopes.shape}")

print("\n📂 Loading TET file for Salinity_Category...")
df_tet = pd.read_csv(TET_FILE)
print(f"   TET file shape: {df_tet.shape}")

# Merge to add Salinity_Category
df = df_slopes.merge(df_tet[['site_name', 'Salinity_Category']], on='site_name', how='left')
print(f"\n   Merged shape: {df.shape}")

# Check for missing salinity
missing_salinity = df['Salinity_Category'].isna().sum()
if missing_salinity > 0:
    print(f"   ⚠️ {missing_salinity} sites missing Salinity_Category - dropping")
    df = df.dropna(subset=['Salinity_Category']).copy()

# Get slope column
if 'slope_theilsen' in df.columns:
    slope_col = 'slope_theilsen'
elif 'slope' in df.columns:
    slope_col = 'slope'
else:
    print(f"❌ No slope column found")
    exit()

df = df.rename(columns={slope_col: 'slope'})
df['abs_slope'] = df['slope'].abs()

print(f"\n✅ Final dataset: {len(df)} sites")
print(f"   Slope column: slope (SPEI-48)")

# ============================================
# STEP 2: GROUP DISTRIBUTION
# ============================================

print("\n" + "="*60)
print("SALINITY CATEGORY DISTRIBUTION")
print("="*60)

sal_counts = df['Salinity_Category'].value_counts()
for cat in ECOSYSTEM_ORDER:
    count = sal_counts.get(cat, 0)
    print(f"   {cat}: {count} sites")

# Keep only the 4 main categories for analysis
df_clean = df[df['Salinity_Category'].isin(ECOSYSTEM_ORDER)].copy()
print(f"\n   Total for analysis (4 categories): {len(df_clean)} sites")

# ============================================
# STEP 3: GROUP STATISTICS
# ============================================

print("\n" + "="*60)
print("GROUP STATISTICS (SPEI-48)")
print("="*60)

# Calculate statistics for each group
stats_results = []
for group in ECOSYSTEM_ORDER:
    group_df = df_clean[df_clean['Salinity_Category'] == group]
    n = len(group_df)
    
    if n > 0:
        slope_median = group_df['slope'].median()
        slope_q25 = group_df['slope'].quantile(0.25)
        slope_q75 = group_df['slope'].quantile(0.75)
        slope_min = group_df['slope'].min()
        slope_max = group_df['slope'].max()
        n_pos = (group_df['slope'] > 0).sum()
        n_neg = (group_df['slope'] < 0).sum()
        pct_pos = (n_pos / n * 100)
        
        n_sig = 0
        if 'is_significant' in df.columns:
            n_sig = (group_df['is_significant'] == True).sum()
        
        stats_results.append({
            'Group': group,
            'N': n,
            'Median_Slope': slope_median,
            'Q1_Slope': slope_q25,
            'Q3_Slope': slope_q75,
            'Min_Slope': slope_min,
            'Max_Slope': slope_max,
            'N_Positive': n_pos,
            'N_Negative': n_neg,
            'Percent_Positive': pct_pos,
            'N_Significant': n_sig if n_sig > 0 else None
        })
        
        print(f"\n{group} (n={n}):")
        print(f"   Median = {slope_median:.3f} [{slope_q25:.3f}, {slope_q75:.3f}]")
        print(f"   Range = {slope_min:.2f} to {slope_max:.2f}")
        print(f"   Positive: {n_pos} ({pct_pos:.0f}%)")
        print(f"   Negative: {n_neg} ({100-pct_pos:.0f}%)")

# Save summary
summary_df = pd.DataFrame(stats_results)
summary_csv = OUTPUT_DIR / 'ecosystem_sensitivity_stats_SPEI48.csv'
summary_df.to_csv(summary_csv, index=False)
print(f"\n✅ Saved statistics: {summary_csv}")

# ============================================
# STEP 4: FULL 4-GROUP KRUSKAL-WALLIS TEST
# ============================================

print("\n" + "="*60)
print("FULL 4-GROUP ANALYSIS (Upland, Freshwater, Brackish, Saline)")
print("="*60)

groups_present = [g for g in ECOSYSTEM_ORDER if g in df_clean['Salinity_Category'].values]
group_data = [df_clean[df_clean['Salinity_Category'] == g]['slope'].values for g in groups_present]

if len(group_data) >= 2:
    kw_stat, kw_p = kruskal(*group_data)
    n_total = len(df_clean)
    epsilon_sq = kw_stat / (n_total - 1)
    
    print(f"\n📊 Kruskal-Wallis Test (Slope to SPEI-48):")
    print(f"   H-statistic = {kw_stat:.3f}")
    print(f"   p-value = {kw_p:.4f}")
    print(f"   Effect size (ε²) = {epsilon_sq:.3f}")
    
    if kw_p < 0.05:
        print(f"   ✓ Significant difference detected among groups")
    else:
        print(f"   ✗ No significant difference detected")
else:
    kw_stat, kw_p, epsilon_sq = None, None, None

# ============================================
# STEP 5: ALL PAIRWISE COMPARISONS (FDR-CORRECTED)
# ============================================

print("\n" + "="*60)
print("ALL PAIRWISE COMPARISONS (FDR-Corrected)")
print("="*60)

all_pairwise_results = []

for i, group1 in enumerate(groups_present):
    for group2 in groups_present[i+1:]:
        data1 = df_clean[df_clean['Salinity_Category'] == group1]['slope'].values
        data2 = df_clean[df_clean['Salinity_Category'] == group2]['slope'].values
        
        if len(data1) > 0 and len(data2) > 0:
            u_stat, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
            
            # Calculate rank-biserial correlation
            n1, n2 = len(data1), len(data2)
            r_biserial = 1 - (2 * u_stat) / (n1 * n2)
            
            # Also test absolute slopes
            abs1 = np.abs(data1)
            abs2 = np.abs(data2)
            u_abs, p_abs = mannwhitneyu(abs1, abs2, alternative='two-sided')
            r_biserial_abs = 1 - (2 * u_abs) / (n1 * n2)
            
            all_pairwise_results.append({
                'Group1': group1,
                'Group2': group2,
                'n1': n1,
                'n2': n2,
                'median1': np.median(data1),
                'median2': np.median(data2),
                'U_statistic': u_stat,
                'p_value': p_val,
                'effect_size_r': r_biserial,
                'abs_U_statistic': u_abs,
                'abs_p_value': p_abs,
                'abs_effect_size_r': r_biserial_abs
            })

# Apply FDR correction
if all_pairwise_results:
    p_values = [r['p_value'] for r in all_pairwise_results]
    reject, p_corrected = fdrcorrection(p_values, alpha=0.05)
    
    for i, result in enumerate(all_pairwise_results):
        result['p_corrected'] = p_corrected[i]
        result['significant'] = reject[i]

print("\n📊 Pairwise Comparison Results (Slope):")
for result in all_pairwise_results:
    sig_marker = "✓ SIGNIFICANT" if result['significant'] else "✗ NOT significant"
    print(f"\n   {result['Group1']} vs {result['Group2']}:")
    print(f"      n = {result['n1']} vs {result['n2']}")
    print(f"      median = {result['median1']:.3f} vs {result['median2']:.3f}")
    print(f"      U = {result['U_statistic']:.1f}, p_corrected = {result['p_corrected']:.4f}")
    print(f"      effect size r = {result['effect_size_r']:.3f}")
    print(f"      {sig_marker}")

# ============================================
# STEP 6: ANALYSIS WITHOUT BRACKISH (3 GROUPS)
# ============================================

print("\n" + "="*60)
print("ANALYSIS WITHOUT BRACKISH (Upland, Freshwater, Saline)")
print("="*60)

df_no_brackish = df_clean[df_clean['Salinity_Category'] != 'Brackish'].copy()
groups_no_brackish = [g for g in ['Upland', 'Freshwater', 'Saline'] if g in df_no_brackish['Salinity_Category'].values]
group_data_no_brackish = [df_no_brackish[df_no_brackish['Salinity_Category'] == g]['slope'].values for g in groups_no_brackish]

if len(group_data_no_brackish) >= 2:
    kw_no_brackish_stat, kw_no_brackish_p = kruskal(*group_data_no_brackish)
    n_no_brackish = len(df_no_brackish)
    epsilon_sq_no_brackish = kw_no_brackish_stat / (n_no_brackish - 1)
    
    print(f"\n📊 Kruskal-Wallis Test (Without Brackish):")
    print(f"   Groups: {groups_no_brackish}")
    print(f"   H-statistic = {kw_no_brackish_stat:.3f}")
    print(f"   p-value = {kw_no_brackish_p:.4f}")
    print(f"   Effect size (ε²) = {epsilon_sq_no_brackish:.3f}")
    
    if kw_no_brackish_p < 0.05:
        print(f"   ✓ Significant difference detected among groups (without Brackish)")
    else:
        print(f"   ✗ No significant difference detected (without Brackish)")
    
    # Pairwise comparisons without Brackish
    print(f"\n📊 Pairwise Comparisons (Without Brackish):")
    for i, group1 in enumerate(groups_no_brackish):
        for group2 in groups_no_brackish[i+1:]:
            data1 = df_no_brackish[df_no_brackish['Salinity_Category'] == group1]['slope'].values
            data2 = df_no_brackish[df_no_brackish['Salinity_Category'] == group2]['slope'].values
            
            if len(data1) > 0 and len(data2) > 0:
                u_stat, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
                n1, n2 = len(data1), len(data2)
                r_biserial = 1 - (2 * u_stat) / (n1 * n2)
                
                sig_marker = "✓ SIGNIFICANT" if p_val < 0.05 else "✗ NOT significant"
                print(f"\n      {group1} vs {group2}:")
                print(f"         n = {n1} vs {n2}")
                print(f"         U = {u_stat:.1f}, p = {p_val:.4f}")
                print(f"         effect size r = {r_biserial:.3f}")
                print(f"         {sig_marker}")

# ============================================
# STEP 7: SPECIFIC COMPARISONS (Freshwater vs Specific Groups)
# ============================================

print("\n" + "="*60)
print("SPECIFIC COMPARISONS")
print("="*60)

# Comparison 1: Freshwater vs Saline (excluding Brackish)
print("\n📊 Comparison 1: Freshwater vs Saline (Excluding Brackish)")
fw_saline_df = df_clean[df_clean['Salinity_Category'].isin(['Freshwater', 'Saline'])].copy()
fw_slopes = fw_saline_df[fw_saline_df['Salinity_Category'] == 'Freshwater']['slope'].values
saline_slopes = fw_saline_df[fw_saline_df['Salinity_Category'] == 'Saline']['slope'].values

if len(fw_slopes) > 0 and len(saline_slopes) > 0:
    u_stat, p_val = mannwhitneyu(fw_slopes, saline_slopes, alternative='two-sided')
    n1, n2 = len(fw_slopes), len(saline_slopes)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    print(f"   Freshwater: n={n1}, median={np.median(fw_slopes):.3f}")
    print(f"   Saline: n={n2}, median={np.median(saline_slopes):.3f}")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {p_val:.4f}")
    print(f"   Effect size r = {r_biserial:.3f}")
    
    # Absolute slopes
    fw_abs = np.abs(fw_slopes)
    saline_abs = np.abs(saline_slopes)
    u_abs, p_abs = mannwhitneyu(fw_abs, saline_abs, alternative='two-sided')
    r_biserial_abs = 1 - (2 * u_abs) / (n1 * n2)
    print(f"\n   |Slope| comparison:")
    print(f"   U = {u_abs:.1f}, p = {p_abs:.4f}, r = {r_biserial_abs:.3f}")

# Comparison 2: Freshwater vs Upland
print("\n📊 Comparison 2: Freshwater vs Upland")
fw_upland_df = df_clean[df_clean['Salinity_Category'].isin(['Freshwater', 'Upland'])].copy()
fw_slopes = fw_upland_df[fw_upland_df['Salinity_Category'] == 'Freshwater']['slope'].values
upland_slopes = fw_upland_df[fw_upland_df['Salinity_Category'] == 'Upland']['slope'].values

if len(fw_slopes) > 0 and len(upland_slopes) > 0:
    u_stat, p_val = mannwhitneyu(fw_slopes, upland_slopes, alternative='two-sided')
    n1, n2 = len(fw_slopes), len(upland_slopes)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    print(f"   Freshwater: n={n1}, median={np.median(fw_slopes):.3f}")
    print(f"   Upland: n={n2}, median={np.median(upland_slopes):.3f}")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {p_val:.4f}")
    print(f"   Effect size r = {r_biserial:.3f}")

# Comparison 3: Upland vs Saline
print("\n📊 Comparison 3: Upland vs Saline")
upland_saline_df = df_clean[df_clean['Salinity_Category'].isin(['Upland', 'Saline'])].copy()
upland_slopes = upland_saline_df[upland_saline_df['Salinity_Category'] == 'Upland']['slope'].values
saline_slopes = upland_saline_df[upland_saline_df['Salinity_Category'] == 'Saline']['slope'].values

if len(upland_slopes) > 0 and len(saline_slopes) > 0:
    u_stat, p_val = mannwhitneyu(upland_slopes, saline_slopes, alternative='two-sided')
    n1, n2 = len(upland_slopes), len(saline_slopes)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    print(f"   Upland: n={n1}, median={np.median(upland_slopes):.3f}")
    print(f"   Saline: n={n2}, median={np.median(saline_slopes):.3f}")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {p_val:.4f}")
    print(f"   Effect size r = {r_biserial:.3f}")

# Comparison 4: Freshwater vs Salt-affected (Brackish + Saline combined)
print("\n📊 Comparison 4: Freshwater vs Salt-affected (Brackish + Saline)")
df_clean['two_group'] = df_clean['Salinity_Category'].copy()
df_clean.loc[df_clean['Salinity_Category'].isin(['Brackish', 'Saline']), 'two_group'] = 'Salt-affected'
df_clean.loc[df_clean['Salinity_Category'] == 'Freshwater', 'two_group'] = 'Freshwater'

df_two_group = df_clean[df_clean['two_group'].isin(['Freshwater', 'Salt-affected'])].copy()

fw_slopes = df_two_group[df_two_group['two_group'] == 'Freshwater']['slope'].values
sa_slopes = df_two_group[df_two_group['two_group'] == 'Salt-affected']['slope'].values

if len(fw_slopes) > 0 and len(sa_slopes) > 0:
    u_stat, mw_p = mannwhitneyu(fw_slopes, sa_slopes, alternative='two-sided')
    n1, n2 = len(fw_slopes), len(sa_slopes)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    print(f"   Freshwater: n={n1}, median={np.median(fw_slopes):.3f}")
    print(f"   Salt-affected: n={n2}, median={np.median(sa_slopes):.3f}")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {mw_p:.4f}")
    print(f"   Effect size r = {r_biserial:.3f}")
    
    # Test for |slope|
    fw_abs = np.abs(fw_slopes)
    sa_abs = np.abs(sa_slopes)
    u_abs, mw_abs_p = mannwhitneyu(fw_abs, sa_abs, alternative='two-sided')
    r_biserial_abs = 1 - (2 * u_abs) / (n1 * n2)
    
    print(f"\n   |Slope| comparison (sensitivity magnitude):")
    print(f"   Mann-Whitney U = {u_abs:.1f}, p = {mw_abs_p:.4f}")
    print(f"   Effect size r = {r_biserial_abs:.3f}")
    
    two_group_results = {
        'Freshwater_n': n1,
        'Salt_affected_n': n2,
        'Freshwater_median': np.median(fw_slopes),
        'Salt_affected_median': np.median(sa_slopes),
        'Mann_Whitney_U': u_stat,
        'p_value': mw_p,
        'effect_size_r': r_biserial,
        'abs_p_value': mw_abs_p,
        'abs_effect_size_r': r_biserial_abs
    }
else:
    two_group_results = None

# Comparison 5: Brackish vs Saline
print("\n📊 Comparison 5: Brackish vs Saline")
brackish_saline_df = df_clean[df_clean['Salinity_Category'].isin(['Brackish', 'Saline'])].copy()
brackish_slopes = brackish_saline_df[brackish_saline_df['Salinity_Category'] == 'Brackish']['slope'].values
saline_slopes = brackish_saline_df[brackish_saline_df['Salinity_Category'] == 'Saline']['slope'].values

if len(brackish_slopes) > 0 and len(saline_slopes) > 0:
    u_stat, p_val = mannwhitneyu(brackish_slopes, saline_slopes, alternative='two-sided')
    n1, n2 = len(brackish_slopes), len(saline_slopes)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    print(f"   Brackish: n={n1}, median={np.median(brackish_slopes):.3f}")
    print(f"   Saline: n={n2}, median={np.median(saline_slopes):.3f}")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {p_val:.4f}")
    print(f"   Effect size r = {r_biserial:.3f}")

# ============================================
# STEP 8: SAVE RESULTS FOR FIGURE
# ============================================

print("\n" + "="*60)
print("SAVING PREPARED DATA FOR FIGURE (SPEI-48)")
print("="*60)

# Prepare data for visualization
figure_data = {
    'df_clean': df_clean,
    'ECOSYSTEM_COLORS': ECOSYSTEM_COLORS,
    'ECOSYSTEM_ORDER': ECOSYSTEM_ORDER,
    'groups_present': groups_present,
    'spei_version': 'SPEI-48',
    'slope_column': 'slope',
    'statistics': {
        'kw_stat': kw_stat,
        'kw_p': kw_p,
        'epsilon_sq': epsilon_sq,
        'kw_no_brackish_stat': kw_no_brackish_stat if 'kw_no_brackish_stat' in dir() else None,
        'kw_no_brackish_p': kw_no_brackish_p if 'kw_no_brackish_p' in dir() else None,
        'pairwise_results': all_pairwise_results,
        'two_group_results': two_group_results,
        'freshwater_vs_saline': {
            'u_stat': u_stat if 'u_stat' in dir() else None,
            'p_value': p_val if 'p_val' in dir() else None,
            'r_biserial': r_biserial if 'r_biserial' in dir() else None
        } if 'u_stat' in dir() else None,
        'freshwater_vs_upland': {
            'u_stat': u_stat_fw_upland if 'u_stat_fw_upland' in dir() else None,
            'p_value': p_val_fw_upland if 'p_val_fw_upland' in dir() else None
        } if 'u_stat_fw_upland' in dir() else None
    }
}

# Save to pickle
import pickle
pickle_file = OUTPUT_DIR / 'ecosystem_diagnostic_data_SPEI48.pkl'
with open(pickle_file, 'wb') as f:
    pickle.dump(figure_data, f)
print(f"✅ Saved prepared data: {pickle_file}")

# ============================================
# FINAL SUMMARY FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT RESULTS SUMMARY (SPEI-48)")
print("="*80)

print(f"\n📊 SAMPLE SIZES:")
for group in ECOSYSTEM_ORDER:
    n = len(df_clean[df_clean['Salinity_Category'] == group])
    if n > 0:
        print(f"   {group}: n={n}")

print(f"\n📊 FULL 4-GROUP COMPARISON:")
if kw_p is not None:
    if kw_p < 0.05:
        print(f"   ✓ Significant difference (Kruskal-Wallis H={kw_stat:.2f}, p={kw_p:.4f})")
    else:
        print(f"   ✗ No significant difference (Kruskal-Wallis H={kw_stat:.2f}, p={kw_p:.4f})")

print(f"\n📊 SIGNIFICANT PAIRWISE DIFFERENCES (FDR-corrected, α=0.05):")
sig_pairs = [r for r in all_pairwise_results if r['significant']]
if sig_pairs:
    for r in sig_pairs:
        print(f"   ✓ {r['Group1']} vs {r['Group2']}: p_corr={r['p_corrected']:.4f}, r={r['effect_size_r']:.3f}")
else:
    print("   ✗ No significant pairwise differences")

print("\n" + "="*80)
print("CHUNK 1 COMPLETE (SPEI-48)")
print("="*80)
print(f"\n📁 Output saved to: {OUTPUT_DIR}")
print(f"   - Statistics CSV: ecosystem_sensitivity_stats_SPEI48.csv")
print(f"   - Pickle data: ecosystem_diagnostic_data_SPEI48.pkl")
print("="*80)