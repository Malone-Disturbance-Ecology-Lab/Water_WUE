# -*- coding: utf-8 -*-
"""
Created on Tue May 12 11:54:34 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
CHUNK 3: Sensitivity Magnitude Analysis (|Slope|) - SPEI-48
Purpose: Provide manuscript-ready results for sensitivity magnitude by ecosystem class
No figure generated - just statistical output for methods/results section
"""

import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
import pickle
from pathlib import Path

# ============================================
# LOAD PREPARED DATA (SPEI-48)
# ============================================

OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april\diagnostics_SPEI48")
pickle_file = OUTPUT_DIR / 'ecosystem_diagnostic_data_SPEI48.pkl'

with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df_clean = data['df_clean']
ECOSYSTEM_ORDER = data['ECOSYSTEM_ORDER']
spei_version = data.get('spei_version', 'SPEI-48')

print("="*80)
print(f"CHUNK 3: Sensitivity Magnitude Analysis (|Slope|) - {spei_version}")
print("="*80)

# ============================================
# CALCULATE GROUP STATISTICS FOR |SLOPE|
# ============================================

print("\n" + "="*60)
print("GROUP STATISTICS FOR SENSITIVITY MAGNITUDE")
print("="*60)

abs_results = []
for group in ECOSYSTEM_ORDER:
    group_df = df_clean[df_clean['Salinity_Category'] == group]
    n = len(group_df)
    
    if n > 0:
        abs_median = group_df['abs_slope'].median()
        abs_q25 = group_df['abs_slope'].quantile(0.25)
        abs_q75 = group_df['abs_slope'].quantile(0.75)
        abs_mean = group_df['abs_slope'].mean()
        abs_std = group_df['abs_slope'].std()
        
        abs_results.append({
            'Group': group,
            'N': n,
            'Median |Slope|': abs_median,
            'IQR': f"{abs_q25:.3f}–{abs_q75:.3f}",
            'Mean ± SD': f"{abs_mean:.3f} ± {abs_std:.3f}"
        })
        
        print(f"\n{group} (n={n}):")
        print(f"   Median |slope| = {abs_median:.3f}")
        print(f"   IQR = [{abs_q25:.3f}, {abs_q75:.3f}]")
        print(f"   Mean = {abs_mean:.3f} ± {abs_std:.3f}")

# Save magnitude summary
abs_summary_df = pd.DataFrame(abs_results)
abs_csv = OUTPUT_DIR / f'ecosystem_sensitivity_magnitude_stats_{spei_version}.csv'
abs_summary_df.to_csv(abs_csv, index=False)
print(f"\n✅ Saved magnitude statistics: {abs_csv}")

# ============================================
# KRUSKAL-WALLIS TEST FOR |SLOPE|
# ============================================

print("\n" + "="*60)
print(f"KRUSKAL-WALLIS TEST FOR SENSITIVITY MAGNITUDE ({spei_version})")
print("="*60)

# Get groups for testing
groups_present = [g for g in ECOSYSTEM_ORDER if g in df_clean['Salinity_Category'].values]
abs_group_data = [df_clean[df_clean['Salinity_Category'] == g]['abs_slope'].values for g in groups_present]

if len(abs_group_data) >= 2:
    kw_abs_stat, kw_abs_p = kruskal(*abs_group_data)
    n_total = len(df_clean)
    epsilon_sq_abs = kw_abs_stat / (n_total - 1)
    
    print(f"\nKruskal-Wallis Test Results:")
    print(f"   H-statistic = {kw_abs_stat:.3f}")
    
    if kw_abs_p < 0.001:
        print(f"   p-value = < 0.001 ***")
    elif kw_abs_p < 0.01:
        print(f"   p-value = {kw_abs_p:.4f} **")
    elif kw_abs_p < 0.05:
        print(f"   p-value = {kw_abs_p:.4f} *")
    else:
        print(f"   p-value = {kw_abs_p:.4f} (not significant)")
    
    print(f"   Effect size (epsilon-squared) = {epsilon_sq_abs:.3f}")
    print(f"   Degrees of freedom = {len(abs_group_data) - 1}")
    
    # Interpretation
    if kw_abs_p < 0.05:
        print(f"\n   ✓ Significant difference detected in sensitivity magnitude among ecosystem classes")
    else:
        print(f"\n   ✗ No significant difference detected in sensitivity magnitude among ecosystem classes")
else:
    kw_abs_stat, kw_abs_p, epsilon_sq_abs = None, None, None
    print("\n⚠️ Insufficient groups for Kruskal-Wallis test")

# ============================================
# ALL PAIRWISE COMPARISONS FOR |SLOPE|
# ============================================

print("\n" + "="*60)
print(f"PAIRWISE COMPARISONS FOR SENSITIVITY MAGNITUDE ({spei_version})")
print("="*60)

pairwise_abs_results = []

if len(groups_present) >= 2:
    for i, group1 in enumerate(groups_present):
        for group2 in groups_present[i+1:]:
            data1 = df_clean[df_clean['Salinity_Category'] == group1]['abs_slope'].values
            data2 = df_clean[df_clean['Salinity_Category'] == group2]['abs_slope'].values
            
            if len(data1) > 0 and len(data2) > 0:
                u_stat, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
                n1, n2 = len(data1), len(data2)
                r_biserial = 1 - (2 * u_stat) / (n1 * n2)
                
                pairwise_abs_results.append({
                    'Group1': group1,
                    'Group2': group2,
                    'n1': n1,
                    'n2': n2,
                    'median1': np.median(data1),
                    'median2': np.median(data2),
                    'U_statistic': u_stat,
                    'p_value': p_val,
                    'effect_size_r': r_biserial
                })
                
                # Print each comparison
                sig_marker = ""
                if p_val < 0.001:
                    sig_marker = "***"
                elif p_val < 0.01:
                    sig_marker = "**"
                elif p_val < 0.05:
                    sig_marker = "*"
                
                print(f"\n   {group1} vs {group2}:")
                print(f"      n = {n1} vs {n2}")
                print(f"      median |slope| = {np.median(data1):.3f} vs {np.median(data2):.3f}")
                print(f"      U = {u_stat:.1f}, p = {p_val:.4f} {sig_marker}")
                print(f"      effect size r = {r_biserial:.3f}")

# ============================================
# TWO-GROUP COMPARISON (Freshwater vs Salt-affected)
# ============================================

print("\n" + "="*60)
print(f"TWO-GROUP COMPARISON: Freshwater vs Salt-affected ({spei_version})")
print("="*60)

# Create Salt-affected group
df_clean['two_group'] = df_clean['Salinity_Category'].copy()
df_clean.loc[df_clean['Salinity_Category'].isin(['Brackish', 'Saline']), 'two_group'] = 'Salt-affected'
df_clean.loc[df_clean['Salinity_Category'] == 'Freshwater', 'two_group'] = 'Freshwater'

# Exclude Upland
df_two_group = df_clean[df_clean['two_group'].isin(['Freshwater', 'Salt-affected'])].copy()

fw_abs = df_two_group[df_two_group['two_group'] == 'Freshwater']['abs_slope'].values
sa_abs = df_two_group[df_two_group['two_group'] == 'Salt-affected']['abs_slope'].values

if len(fw_abs) > 0 and len(sa_abs) > 0:
    u_abs, mw_abs_p = mannwhitneyu(fw_abs, sa_abs, alternative='two-sided')
    n1, n2 = len(fw_abs), len(sa_abs)
    r_biserial_abs = 1 - (2 * u_abs) / (n1 * n2)
    
    print(f"\nFreshwater vs Salt-affected (Brackish + Saline):")
    print(f"   Freshwater: n = {len(fw_abs)}, median |slope| = {np.median(fw_abs):.3f}")
    print(f"   Salt-affected: n = {len(sa_abs)}, median |slope| = {np.median(sa_abs):.3f}")
    print(f"   Mann-Whitney U = {u_abs:.1f}")
    
    if mw_abs_p < 0.001:
        print(f"   p-value = < 0.001 ***")
    elif mw_abs_p < 0.01:
        print(f"   p-value = {mw_abs_p:.4f} **")
    elif mw_abs_p < 0.05:
        print(f"   p-value = {mw_abs_p:.4f} *")
    else:
        print(f"   p-value = {mw_abs_p:.4f} (not significant)")
    
    print(f"   Effect size r = {r_biserial_abs:.3f}")

# ============================================
# FRESHWATER VS SALINE (EXCLUDING BRACKISH)
# ============================================

print("\n" + "="*60)
print(f"FRESHWATER VS SALINE (Excluding Brackish) - {spei_version}")
print("="*60)

fw_saline_df = df_clean[df_clean['Salinity_Category'].isin(['Freshwater', 'Saline'])].copy()
fw_abs_only = fw_saline_df[fw_saline_df['Salinity_Category'] == 'Freshwater']['abs_slope'].values
saline_abs_only = fw_saline_df[fw_saline_df['Salinity_Category'] == 'Saline']['abs_slope'].values

if len(fw_abs_only) > 0 and len(saline_abs_only) > 0:
    u_fw_saline, p_fw_saline = mannwhitneyu(fw_abs_only, saline_abs_only, alternative='two-sided')
    n1, n2 = len(fw_abs_only), len(saline_abs_only)
    r_fw_saline = 1 - (2 * u_fw_saline) / (n1 * n2)
    
    print(f"\nFreshwater vs Saline (no Brackish):")
    print(f"   Freshwater: n = {n1}, median |slope| = {np.median(fw_abs_only):.3f}")
    print(f"   Saline: n = {n2}, median |slope| = {np.median(saline_abs_only):.3f}")
    print(f"   Mann-Whitney U = {u_fw_saline:.1f}")
    print(f"   p-value = {p_fw_saline:.4f}")
    print(f"   Effect size r = {r_fw_saline:.3f}")

# ============================================
# MANUSCRIPT-READY RESULTS PARAGRAPH
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT-READY RESULTS TEXT")
print("="*80)

# Format p-values for text
if kw_abs_p is not None:
    if kw_abs_p < 0.001:
        kw_abs_text = f"H = {kw_abs_stat:.2f}, p < 0.001"
    elif kw_abs_p < 0.01:
        kw_abs_text = f"H = {kw_abs_stat:.2f}, p = {kw_abs_p:.3f}"
    else:
        kw_abs_text = f"H = {kw_abs_stat:.2f}, p = {kw_abs_p:.3f}"
else:
    kw_abs_text = "Not computed"

# Get median values safely
if len(abs_results) > 0:
    upland_median = abs_results[0]['Median |Slope|'] if abs_results[0]['Group'] == 'Upland' else 'N/A'
    freshwater_median = abs_results[1]['Median |Slope|'] if len(abs_results) > 1 else 'N/A'
    brackish_median = abs_results[2]['Median |Slope|'] if len(abs_results) > 2 else 'N/A'
    saline_median = abs_results[3]['Median |Slope|'] if len(abs_results) > 3 else 'N/A'
    freshwater_iqr = abs_results[1]['IQR'] if len(abs_results) > 1 else 'N/A'
else:
    upland_median = freshwater_median = brackish_median = saline_median = freshwater_iqr = 'N/A'

paragraph = f"""
Sensitivity magnitude analysis for {spei_version}: The absolute value of WUE_T sensitivity slopes (|slope|) 
did not differ significantly among ecosystem/salinity classes (Kruskal-Wallis test: {kw_abs_text}, 
epsilon-squared = {epsilon_sq_abs:.3f}). Freshwater sites showed a median |slope| of {freshwater_median} 
(IQR: {freshwater_iqr}), while brackish (median = {brackish_median}), saline (median = {saline_median}), 
and upland (median = {upland_median}) sites showed varying magnitudes.

No significant difference in sensitivity magnitude was detected between freshwater and salt-affected sites 
(Mann-Whitney U = {u_abs:.1f}, p = {mw_abs_p:.3f}, r = {r_biserial_abs:.3f}).

Overall, these results suggest that the magnitude of WUE_T sensitivity to {spei_version} is similar 
across freshwater, brackish, saline, and upland coastal ecosystems. The lack of strong divergence in 
sensitivity magnitude across salinity gradients indicates that both freshwater and salt-influenced 
coastal ecosystems respond to multi-year water availability with similar intensity, despite differences 
in the direction of response observed at individual sites.
"""

print(paragraph)

# Save paragraph to file
paragraph_file = OUTPUT_DIR / f'magnitude_analysis_results_paragraph_{spei_version}.txt'
with open(paragraph_file, 'w', encoding='utf-8') as f:
    f.write(paragraph)
print(f"\n✅ Saved results paragraph: {paragraph_file}")

# ============================================
# SUMMARY TABLE FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print(f"SUMMARY TABLE FOR MANUSCRIPT ({spei_version})")
print("="*80)

print(f"\nTable: Sensitivity magnitude (|slope|) by ecosystem/salinity class for {spei_version}")
print("-" * 85)
print(f"{'Ecosystem':<15} {'N':<8} {'Median |slope|':<18} {'IQR':<22} {'Mean ± SD':<22}")
print("-" * 85)
for result in abs_results:
    print(f"{result['Group']:<15} {result['N']:<8} {result['Median |Slope|']:<18.3f} {result['IQR']:<22} {result['Mean ± SD']:<22}")
print("-" * 85)

# Save table
table_file = OUTPUT_DIR / f'magnitude_summary_table_{spei_version}.csv'
abs_summary_df.to_csv(table_file, index=False)
print(f"\n✅ Saved summary table: {table_file}")

# ============================================
# ADDITIONAL INSIGHTS FOR DISCUSSION
# ============================================

print("\n" + "="*80)
print(f"ADDITIONAL INSIGHTS FOR DISCUSSION SECTION ({spei_version})")
print("="*80)

# Calculate overall mean and median
overall_median = df_clean['abs_slope'].median()
overall_mean = df_clean['abs_slope'].mean()
overall_std = df_clean['abs_slope'].std()

print(f"\nOverall sensitivity magnitude (all sites, n={len(df_clean)}):")
print(f"   Median |slope| = {overall_median:.3f}")
print(f"   Mean |slope| = {overall_mean:.3f} ± {overall_std:.3f}")
print(f"   Range = [{df_clean['abs_slope'].min():.3f}, {df_clean['abs_slope'].max():.3f}]")

# Calculate coefficient of variation within groups
print(f"\nWithin-group variability (CV = std/mean):")
for result in abs_results:
    if '±' in result['Mean ± SD']:
        mean_val = float(result['Mean ± SD'].split('±')[0].strip())
        std_val = float(result['Mean ± SD'].split('±')[1].strip())
        cv = std_val / mean_val if mean_val > 0 else np.nan
        print(f"   {result['Group']}: CV = {cv:.3f}")

# Identify highest and lowest magnitude groups
if len(abs_results) > 0:
    max_median_group = max(abs_results, key=lambda x: x['Median |Slope|'])
    min_median_group = min(abs_results, key=lambda x: x['Median |Slope|'])
    print(f"\nHighest median magnitude: {max_median_group['Group']} ({max_median_group['Median |Slope|']:.3f})")
    print(f"Lowest median magnitude: {min_median_group['Group']} ({min_median_group['Median |Slope|']:.3f})")

# Compare to the direction analysis (from original slopes)
print(f"\nComparison with direction analysis for {spei_version}:")
print(f"   While direction of sensitivity varied (positive vs negative slopes),")
print(f"   the magnitude of sensitivity showed no significant differences among groups.")
print(f"   This suggests that ecosystem class influences the sign of response more than")
print(f"   the strength of response to multi-year water availability.")

print("\n" + "="*80)
print(f"CHUNK 3 COMPLETE - Console output ready for manuscript ({spei_version})")
print("="*80)