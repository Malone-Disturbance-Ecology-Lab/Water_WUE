# -*- coding: utf-8 -*-
"""
Created on Mon May 11 23:36:01 2026

@author: ammar
"""

"""
CHUNK 1: WUE_T Sensitivity by Ecosystem/Salinity Class - Data Prep & Statistics
Purpose: Load data, run statistical tests, prepare data for visualization
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
SENSITIVITY_FILE = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\Q3_analysis\Q3_WUET_SPEI6_slopes.csv")
OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april\diagnostics")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Ecosystem colors (matching supplementary figure)
ECOSYSTEM_COLORS = {
    'Upland': '#800080',      # Purple
    'Freshwater': '#0000FF',  # Blue
    'Brackish': '#008080',    # Teal
    'Saline': '#FFA500'       # Orange
}

ECOSYSTEM_ORDER = ['Upland', 'Freshwater', 'Brackish', 'Saline']

print("="*80)
print("CHUNK 1: Data Preparation and Statistical Analysis")
print("="*80)

# ============================================
# STEP 1: LOAD DATA
# ============================================

print("\n📂 Loading sensitivity data...")
df = pd.read_csv(SENSITIVITY_FILE)
print(f"   Shape: {df.shape}")
print(f"   Columns: {df.columns.tolist()}")

# ============================================
# STEP 2: DATA CHECKS AND CLEANING
# ============================================

print("\n" + "="*60)
print("DATA CHECKS")
print("="*60)

# Check required columns
required_cols = ['site_name', 'Salinity_Category', 'slope']
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    print(f"❌ Missing columns: {missing_cols}")
    exit()

print(f"\n✅ Using:")
print(f"   - Slope column: slope")
print(f"   - Group column: Salinity_Category")
print(f"   - Site column: site_name")

# Check distribution
print(f"\n📊 Salinity_Category distribution:")
sal_counts = df['Salinity_Category'].value_counts()
for cat in ECOSYSTEM_ORDER:
    count = sal_counts.get(cat, 0)
    print(f"   {cat}: {count} sites")

# Remove missing values
df_clean = df.dropna(subset=['slope', 'Salinity_Category']).copy()
print(f"\n✅ Clean data: {len(df_clean)} sites")

# Add absolute slope
df_clean['abs_slope'] = df_clean['slope'].abs()

# ============================================
# STEP 3: GROUP STATISTICS
# ============================================

print("\n" + "="*60)
print("GROUP STATISTICS")
print("="*60)

# Calculate median slope for ordering (for figure only - keep statistical order)
group_medians = df_clean.groupby('Salinity_Category')['slope'].median().sort_values()
print(f"\n📊 Median slopes for ordering:")
for group, median in group_medians.items():
    print(f"   {group}: {median:.3f}")

# Store results
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
        
        # Check for significance if available
        n_sig = 0
        if 'spearman_p' in df.columns:
            n_sig = (group_df['spearman_p'] < 0.05).sum()
        elif 'pearson_p' in df.columns:
            n_sig = (group_df['pearson_p'] < 0.05).sum()
        
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
summary_csv = OUTPUT_DIR / 'ecosystem_sensitivity_stats.csv'
summary_df.to_csv(summary_csv, index=False)
print(f"\n✅ Saved statistics: {summary_csv}")

# ============================================
# STEP 4: STATISTICAL TESTS
# ============================================

print("\n" + "="*60)
print("STATISTICAL TESTS")
print("="*60)

# Get groups for testing (only those with data)
groups_present = [g for g in ECOSYSTEM_ORDER if g in df_clean['Salinity_Category'].values]
group_data = [df_clean[df_clean['Salinity_Category'] == g]['slope'].values for g in groups_present]

# Test 1: Kruskal-Wallis for slopes
if len(group_data) >= 2:
    kw_stat, kw_p = kruskal(*group_data)
    n_total = len(df_clean)
    epsilon_sq = kw_stat / (n_total - 1)
    
    print(f"\n📊 Kruskal-Wallis Test (Slope):")
    print(f"   H-statistic = {kw_stat:.3f}")
    print(f"   p-value = {kw_p:.4f}")
    print(f"   Effect size (ε²) = {epsilon_sq:.3f}")
    
    if kw_p < 0.05:
        print(f"   ✓ Significant difference detected")
    else:
        print(f"   ✗ No significant difference detected")
else:
    kw_stat, kw_p, epsilon_sq = None, None, None
    print("\n⚠️ Insufficient groups for Kruskal-Wallis test")

# Test 2: Kruskal-Wallis for |slope|
abs_group_data = [df_clean[df_clean['Salinity_Category'] == g]['abs_slope'].values for g in groups_present]
if len(abs_group_data) >= 2:
    kw_abs_stat, kw_abs_p = kruskal(*abs_group_data)
    epsilon_sq_abs = kw_abs_stat / (n_total - 1)
    
    print(f"\n📊 Kruskal-Wallis Test (|Slope| - Sensitivity Magnitude):")
    print(f"   H-statistic = {kw_abs_stat:.3f}")
    print(f"   p-value = {kw_abs_p:.4f}")
    print(f"   Effect size (ε²) = {epsilon_sq_abs:.3f}")
    
    if kw_abs_p < 0.05:
        print(f"   ✓ Significant difference detected")
    else:
        print(f"   ✗ No significant difference detected")
else:
    kw_abs_stat, kw_abs_p, epsilon_sq_abs = None, None, None

# Test 3: Pairwise Mann-Whitney with FDR (if KW significant)
pairwise_results = []
if kw_p is not None and kw_p < 0.05:
    print(f"\n📊 Pairwise Mann-Whitney U Tests (FDR-corrected):")
    
    for i, group1 in enumerate(groups_present):
        for group2 in groups_present[i+1:]:
            data1 = df_clean[df_clean['Salinity_Category'] == group1]['slope'].values
            data2 = df_clean[df_clean['Salinity_Category'] == group2]['slope'].values
            
            if len(data1) > 0 and len(data2) > 0:
                u_stat, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
                
                # Calculate rank-biserial correlation
                n1, n2 = len(data1), len(data2)
                r_biserial = 1 - (2 * u_stat) / (n1 * n2)
                
                pairwise_results.append({
                    'Group1': group1,
                    'Group2': group2,
                    'U_statistic': u_stat,
                    'p_value': p_val,
                    'effect_size_r': r_biserial
                })
    
    # Apply FDR correction
    if pairwise_results:
        p_values = [r['p_value'] for r in pairwise_results]
        reject, p_corrected = fdrcorrection(p_values, alpha=0.05)
        
        for i, result in enumerate(pairwise_results):
            result['p_corrected'] = p_corrected[i]
            result['significant'] = reject[i]
            
            if reject[i]:
                print(f"\n   ✓ {result['Group1']} vs {result['Group2']}:")
                print(f"      p_corrected = {result['p_corrected']:.4f}")
                print(f"      effect size r = {result['effect_size_r']:.3f}")

# Test 4: Two-group comparison (Freshwater vs Salt-affected)
print(f"\n📊 Two-Group Comparison: Freshwater vs Salt-affected")
print(f"   Salt-affected = Brackish + Saline")

# Create Salt-affected group
df_clean['two_group'] = df_clean['Salinity_Category'].copy()
df_clean.loc[df_clean['Salinity_Category'].isin(['Brackish', 'Saline']), 'two_group'] = 'Salt-affected'
df_clean.loc[df_clean['Salinity_Category'] == 'Freshwater', 'two_group'] = 'Freshwater'

# Exclude Upland from this comparison
df_two_group = df_clean[df_clean['two_group'].isin(['Freshwater', 'Salt-affected'])].copy()

fw_slopes = df_two_group[df_two_group['two_group'] == 'Freshwater']['slope'].values
sa_slopes = df_two_group[df_two_group['two_group'] == 'Salt-affected']['slope'].values

if len(fw_slopes) > 0 and len(sa_slopes) > 0:
    u_stat, mw_p = mannwhitneyu(fw_slopes, sa_slopes, alternative='two-sided')
    n1, n2 = len(fw_slopes), len(sa_slopes)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    print(f"\n   Freshwater: n = {len(fw_slopes)}, median = {np.median(fw_slopes):.3f}")
    print(f"   Salt-affected: n = {len(sa_slopes)}, median = {np.median(sa_slopes):.3f}")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {mw_p:.4f}")
    print(f"   Effect size r = {r_biserial:.3f}")
    
    # Test for |slope|
    fw_abs = df_two_group[df_two_group['two_group'] == 'Freshwater']['abs_slope'].values
    sa_abs = df_two_group[df_two_group['two_group'] == 'Salt-affected']['abs_slope'].values
    u_abs, mw_abs_p = mannwhitneyu(fw_abs, sa_abs, alternative='two-sided')
    r_biserial_abs = 1 - (2 * u_abs) / (n1 * n2)
    
    print(f"\n   |Slope| comparison:")
    print(f"   Mann-Whitney U = {u_abs:.1f}, p = {mw_abs_p:.4f}")
    print(f"   Effect size r = {r_biserial_abs:.3f}")
    
    two_group_results = {
        'Freshwater_n': len(fw_slopes),
        'Salt_affected_n': len(sa_slopes),
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

# ============================================
# STEP 5: SAVE RESULTS FOR FIGURE
# ============================================

print("\n" + "="*60)
print("SAVING PREPARED DATA FOR FIGURE")
print("="*60)

# Prepare data for visualization
figure_data = {
    'df_clean': df_clean,
    'ECOSYSTEM_COLORS': ECOSYSTEM_COLORS,
    'ECOSYSTEM_ORDER': ECOSYSTEM_ORDER,
    'groups_present': groups_present,
    'statistics': {
        'kw_stat': kw_stat,
        'kw_p': kw_p,
        'epsilon_sq': epsilon_sq,
        'kw_abs_stat': kw_abs_stat,
        'kw_abs_p': kw_abs_p,
        'epsilon_sq_abs': epsilon_sq_abs,
        'pairwise_results': pairwise_results,
        'two_group_results': two_group_results
    }
}

# Save to pickle
import pickle
pickle_file = OUTPUT_DIR / 'ecosystem_diagnostic_data.pkl'
with open(pickle_file, 'wb') as f:
    pickle.dump(figure_data, f)
print(f"✅ Saved prepared data: {pickle_file}")

# ============================================
# PRINT SUMMARY FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT RESULTS SUMMARY")
print("="*80)

if kw_p is not None:
    if kw_p < 0.05:
        sig_text = f"differed significantly (H = {kw_stat:.2f}, p = {kw_p:.4f}, ε² = {epsilon_sq:.3f})"
    else:
        sig_text = f"did not differ significantly (H = {kw_stat:.2f}, p = {kw_p:.4f}, ε² = {epsilon_sq:.3f})"
    
    print(f"\nSite-level WUE_T sensitivity slopes {sig_text} among ecosystem/salinity classes.")
    
    if kw_abs_p is not None:
        if kw_abs_p < 0.05:
            abs_text = f"differed significantly (H = {kw_abs_stat:.2f}, p = {kw_abs_p:.4f})"
        else:
            abs_text = f"did not differ significantly (H = {kw_abs_stat:.2f}, p = {kw_abs_p:.4f})"
        print(f"Sensitivity magnitude {abs_text} among groups.")
    
    if pairwise_results and any(r['significant'] for r in pairwise_results):
        print("\nSignificant pairwise differences (FDR-corrected):")
        for r in pairwise_results:
            if r['significant']:
                print(f"  {r['Group1']} vs {r['Group2']}: p = {r['p_corrected']:.4f}, r = {r['effect_size_r']:.3f}")

print("\n" + "="*80)
print("CHUNK 1 COMPLETE")
print("="*80)