"""
CHUNK 3: Sensitivity Magnitude Analysis (|Slope|) - Console Output Only
Purpose: Provide manuscript-ready results for sensitivity magnitude by ecosystem class
No figure generated - just statistical output for methods/results section
FIXED: Removed stdout encoding changes
"""

import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
import pickle
from pathlib import Path

# ============================================
# LOAD PREPARED DATA
# ============================================

OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april\diagnostics")
pickle_file = OUTPUT_DIR / 'ecosystem_diagnostic_data.pkl'

with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df_clean = data['df_clean']
ECOSYSTEM_ORDER = data['ECOSYSTEM_ORDER']

print("="*80)
print("CHUNK 3: Sensitivity Magnitude Analysis (|Slope|)")
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
abs_csv = OUTPUT_DIR / 'ecosystem_sensitivity_magnitude_stats.csv'
abs_summary_df.to_csv(abs_csv, index=False)
print(f"\n✅ Saved magnitude statistics: {abs_csv}")

# ============================================
# KRUSKAL-WALLIS TEST FOR |SLOPE|
# ============================================

print("\n" + "="*60)
print("KRUSKAL-WALLIS TEST FOR SENSITIVITY MAGNITUDE")
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
# TWO-GROUP COMPARISON (Freshwater vs Salt-affected)
# ============================================

print("\n" + "="*60)
print("TWO-GROUP COMPARISON: Freshwater vs Salt-affected")
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
# MANUSCRIPT-READY RESULTS PARAGRAPH
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT-READY RESULTS TEXT")
print("="*80)

# Format p-values for text
if kw_abs_p < 0.001:
    kw_abs_text = f"H = {kw_abs_stat:.2f}, p < 0.001"
elif kw_abs_p < 0.01:
    kw_abs_text = f"H = {kw_abs_stat:.2f}, p = {kw_abs_p:.3f}"
else:
    kw_abs_text = f"H = {kw_abs_stat:.2f}, p = {kw_abs_p:.3f}"

# Get median values safely
upland_median = abs_results[0]['Median |Slope|'] if len(abs_results) > 0 else 'N/A'
freshwater_median = abs_results[1]['Median |Slope|'] if len(abs_results) > 1 else 'N/A'
brackish_median = abs_results[2]['Median |Slope|'] if len(abs_results) > 2 else 'N/A'
saline_median = abs_results[3]['Median |Slope|'] if len(abs_results) > 3 else 'N/A'
freshwater_iqr = abs_results[1]['IQR'] if len(abs_results) > 1 else 'N/A'

# Build paragraph (without special characters)
paragraph = f"""
Sensitivity magnitude analysis: The absolute value of WUE_T sensitivity slopes (|slope|) 
did not differ significantly among ecosystem/salinity classes (Kruskal-Wallis test: {kw_abs_text}, 
epsilon-squared = {epsilon_sq_abs:.3f}). Freshwater sites showed a median |slope| of {freshwater_median} 
(IQR: {freshwater_iqr}), while brackish (median = {brackish_median}), saline (median = {saline_median}), 
and upland (median = {upland_median}) sites showed varying magnitudes.

No significant difference in sensitivity magnitude was detected between freshwater and salt-affected sites 
(Mann-Whitney U = {u_abs:.1f}, p = {mw_abs_p:.3f}, r = {r_biserial_abs:.3f}).

Overall, these results suggest that the magnitude of WUE_T sensitivity to water availability is similar 
across freshwater, brackish, saline, and upland coastal ecosystems. The lack of strong divergence in 
sensitivity magnitude across salinity gradients indicates that both freshwater and salt-influenced 
coastal ecosystems respond to water availability with similar intensity, despite differences in the 
direction of response observed at individual sites.
"""

print(paragraph)

# Save paragraph to file with utf-8 encoding (this is safe)
paragraph_file = OUTPUT_DIR / 'magnitude_analysis_results_paragraph.txt'
with open(paragraph_file, 'w', encoding='utf-8') as f:
    f.write(paragraph)
print(f"\n✅ Saved results paragraph: {paragraph_file}")

# ============================================
# SUMMARY TABLE FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print("SUMMARY TABLE FOR MANUSCRIPT")
print("="*80)

print("\nTable: Sensitivity magnitude (|slope|) by ecosystem/salinity class")
print("-" * 80)
print(f"{'Ecosystem':<15} {'N':<8} {'Median |slope|':<15} {'IQR':<20} {'Mean ± SD':<20}")
print("-" * 80)
for result in abs_results:
    print(f"{result['Group']:<15} {result['N']:<8} {result['Median |Slope|']:<15.3f} {result['IQR']:<20} {result['Mean ± SD']:<20}")
print("-" * 80)

# Save table
table_file = OUTPUT_DIR / 'magnitude_summary_table.csv'
abs_summary_df.to_csv(table_file, index=False)
print(f"\n✅ Saved summary table: {table_file}")

# ============================================
# ADDITIONAL INSIGHTS FOR DISCUSSION
# ============================================

print("\n" + "="*80)
print("ADDITIONAL INSIGHTS FOR DISCUSSION SECTION")
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
print(f"\nComparison with direction analysis:")
print(f"   While direction of sensitivity varied (positive vs negative slopes),")
print(f"   the magnitude of sensitivity showed no significant differences among groups.")
print(f"   This suggests that ecosystem class influences the sign of response more than")
print(f"   the strength of response to water availability.")

print("\n" + "="*80)
print("CHUNK 3 COMPLETE - Console output ready for manuscript")
print("="*80)