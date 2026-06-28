"""
DIAGNOSTIC: Compare WUE_T Sensitivity Across SPEI Timescales
Tests which SPEI timescale yields the most significant relationship with WUE_T
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from sklearn.linear_model import TheilSenRegressor
import warnings
warnings.filterwarnings('ignore')

# ============================================
# SETUP
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
diagnostic_dir = output_dir / "spei_timescale_diagnostic"
diagnostic_dir.mkdir(parents=True, exist_ok=True)

# Input files
monthly_filtered_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
nn_medians_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"

print("="*80)
print("DIAGNOSTIC: WUE_T Sensitivity Across SPEI Timescales")
print("="*80)

# ============================================
# STEP 1: LOAD DATA
# ============================================

print("\n📂 STEP 1: Loading filtered data...")

df_monthly = pd.read_csv(monthly_filtered_path)
print(f"  Loaded monthly data: {len(df_monthly)} rows")

# Get strict triple intersection sites (58 sites)
nn_data = pd.read_csv(nn_medians_path)
wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

print(f"  Strict triple intersection sites: {len(shared_sites)}")

# Filter to shared sites
df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
print(f"  Filtered monthly data: {len(df_monthly)} rows")

# ============================================
# STEP 2: SPEI TIMESCALES TO TEST
# ============================================

spei_timescales = {
    'SPEI_1': 1,
    'SPEI_3': 3,
    'SPEI_6': 6,
    'SPEI_12': 12,
    'SPEI_24': 24,
    'SPEI_36': 36,
    'SPEI_48': 48
}

print(f"\n📊 STEP 2: Testing SPEI timescales: {list(spei_timescales.keys())}")

# ============================================
# STEP 3: COMPUTE SITE-LEVEL SIGNIFICANCE FOR EACH TIMESCALE
# ============================================

print("\n📈 STEP 3: Computing site-level slopes and p-values for each SPEI timescale...")

all_results = []

for spei_col, spei_value in spei_timescales.items():
    print(f"\n  Processing {spei_col} (SPEI-{spei_value})...")
    
    site_results = []
    
    for site in shared_sites:
        site_data = df_monthly[df_monthly['site_name'] == site].sort_values(spei_col).copy()
        
        # Drop rows with missing values for this SPEI column
        site_data_clean = site_data.dropna(subset=[spei_col, 'WUE_tra'])
        n_obs = len(site_data_clean)
        
        if n_obs >= 5:
            X = site_data_clean[spei_col].values
            y = site_data_clean['WUE_tra'].values
            
            # Theil-Sen slope (robust)
            try:
                theil = TheilSenRegressor(random_state=42)
                theil.fit(X.reshape(-1, 1), y)
                slope_theilsen = theil.coef_[0]
            except:
                slope_theilsen = np.nan
            
            # Linear regression for p-value
            try:
                slope_ols, intercept, r_value, p_value, std_err = stats.linregress(X, y)
            except:
                slope_ols = np.nan
                p_value = np.nan
            
            site_results.append({
                'site_name': site,
                'spei_timescale': spei_col,
                'spei_value': spei_value,
                'n_obs': n_obs,
                'slope_theilsen': slope_theilsen,
                'slope_ols': slope_ols,
                'p_value': p_value,
                'significant_05': p_value < 0.05 if not np.isnan(p_value) else False,
                'significant_01': p_value < 0.01 if not np.isnan(p_value) else False,
                'r_squared': r_value**2 if not np.isnan(r_value) else np.nan
            })
    
    df_results = pd.DataFrame(site_results)
    df_results = df_results.dropna(subset=['p_value'])
    
    # Calculate summary statistics for this timescale
    n_sites_analyzed = len(df_results)
    n_sig_05 = df_results['significant_05'].sum()
    n_sig_01 = df_results['significant_01'].sum()
    pct_sig_05 = (n_sig_05 / n_sites_analyzed * 100) if n_sites_analyzed > 0 else 0
    pct_sig_01 = (n_sig_01 / n_sites_analyzed * 100) if n_sites_analyzed > 0 else 0
    median_p = df_results['p_value'].median()
    mean_p = df_results['p_value'].mean()
    median_slope = df_results['slope_theilsen'].median()
    
    all_results.append({
        'spei_timescale': spei_col,
        'spei_value': spei_value,
        'n_sites_analyzed': n_sites_analyzed,
        'n_sig_p05': n_sig_05,
        'pct_sig_p05': pct_sig_05,
        'n_sig_p01': n_sig_01,
        'pct_sig_p01': pct_sig_01,
        'median_p_value': median_p,
        'mean_p_value': mean_p,
        'median_slope': median_slope,
        'total_sites_possible': len(shared_sites)
    })
    
    print(f"    Sites analyzed: {n_sites_analyzed}/{len(shared_sites)}")
    print(f"    Significant (p<0.05): {n_sig_05} ({pct_sig_05:.1f}%)")
    print(f"    Significant (p<0.01): {n_sig_01} ({pct_sig_01:.1f}%)")
    print(f"    Median p-value: {median_p:.4f}")
    print(f"    Median slope: {median_slope:+.4f}")

# ============================================
# STEP 4: COMPARE TIMESCALES
# ============================================

print("\n" + "="*80)
print("STEP 4: COMPARISON ACROSS SPEI TIMESCALES")
print("="*80)

df_comparison = pd.DataFrame(all_results)
df_comparison = df_comparison.sort_values('spei_value')

print("\n📊 SUMMARY TABLE:")
print("-"*80)
print(f"{'SPEI':<8} {'Sites':<8} {'p<0.05':<12} {'p<0.01':<12} {'Median p':<12} {'Median Slope':<15}")
print("-"*80)

for _, row in df_comparison.iterrows():
    print(f"SPEI-{row['spei_value']:<3} {row['n_sites_analyzed']:<8} "
          f"{row['n_sig_p05']} ({row['pct_sig_p05']:.0f}%){'':<4} "
          f"{row['n_sig_p01']} ({row['pct_sig_p01']:.0f}%){'':<4} "
          f"{row['median_p_value']:.4f}      "
          f"{row['median_slope']:+.4f}")

print("-"*80)

# ============================================
# STEP 5: IDENTIFY BEST TIMESCALE
# ============================================

print("\n" + "="*80)
print("STEP 5: BEST TIMESCALE IDENTIFICATION")
print("="*80)

# Find timescale with most significant sites (p<0.05)
best_by_sig = df_comparison.loc[df_comparison['n_sig_p05'].idxmax()]
# Find timescale with lowest median p-value
best_by_median_p = df_comparison.loc[df_comparison['median_p_value'].idxmin()]
# Find timescale with best balance (high significance + low p-value)
df_comparison['score'] = df_comparison['n_sig_p05'] / df_comparison['median_p_value']
best_by_score = df_comparison.loc[df_comparison['score'].idxmax()]

print(f"\n🏆 MOST SIGNIFICANT (most sites with p<0.05):")
print(f"   SPEI-{best_by_sig['spei_value']}: {best_by_sig['n_sig_p05']}/{best_by_sig['n_sites_analyzed']} sites ({best_by_sig['pct_sig_p05']:.1f}%)")

print(f"\n📉 LOWEST MEDIAN P-VALUE:")
print(f"   SPEI-{best_by_median_p['spei_value']}: median p = {best_by_median_p['median_p_value']:.4f}")

print(f"\n⭐ BEST BALANCE (most significant + lowest p-value):")
print(f"   SPEI-{best_by_score['spei_value']}: {best_by_score['n_sig_p05']} sig sites, median p={best_by_score['median_p_value']:.4f}")

print(f"\n📊 BEST COVERAGE (most sites with ≥5 obs):")
print(f"   SPEI-{best_by_coverage['spei_value']}: {best_by_coverage['n_sites_analyzed']}/{best_by_coverage['total_sites_possible']} sites")

# ============================================
# STEP 6: SAVE RESULTS
# ============================================

print("\n💾 STEP 6: Saving results...")

# Save comparison summary
comparison_file = diagnostic_dir / 'SPEI_timescale_comparison.csv'
df_comparison.to_csv(comparison_file, index=False)
print(f"  ✓ Saved: {comparison_file}")

# Save detailed site-level results for each timescale
for spei_col in spei_timescales.keys():
    # Get results for this timescale
    site_details = []
    for site in shared_sites:
        site_data = df_monthly[df_monthly['site_name'] == site].dropna(subset=[spei_col, 'WUE_tra'])
        if len(site_data) >= 5:
            X = site_data[spei_col].values
            y = site_data['WUE_tra'].values
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(X, y)
                site_details.append({
                    'site_name': site,
                    'n_obs': len(site_data),
                    'slope': slope,
                    'p_value': p_value,
                    'r_squared': r_value**2,
                    'significant_05': p_value < 0.05
                })
            except:
                continue
    
    if len(site_details) > 0:
        df_detail = pd.DataFrame(site_details)
        detail_file = diagnostic_dir / f'{spei_col}_site_results.csv'
        df_detail.to_csv(detail_file, index=False)
        print(f"  ✓ Saved: {detail_file}")

# ============================================
# STEP 7: RECOMMENDATION
# ============================================

print("\n" + "="*80)
print("STEP 7: RECOMMENDATION")
print("="*80)

# Get current SPEI-6 stats
spei6_stats = df_comparison[df_comparison['spei_timescale'] == 'SPEI_6'].iloc[0]

print(f"\nCurrent Q3 analysis uses: SPEI-6")
print(f"  • Significant sites (p<0.05): {spei6_stats['n_sig_p05']}/{spei6_stats['n_sites_analyzed']} ({spei6_stats['pct_sig_p05']:.1f}%)")
print(f"  • Median p-value: {spei6_stats['median_p_value']:.4f}")

# Compare with best timescale
if best_by_sig['spei_timescale'] != 'SPEI_6':
    print(f"\n⚠️ SPEI-{best_by_sig['spei_value']} shows MORE significant sites ({best_by_sig['n_sig_p05']} vs {spei6_stats['n_sig_p05']})")
    print(f"   Consider using SPEI-{best_by_sig['spei_value']} for future sensitivity analysis")
    print(f"\n   RECOMMENDED: Switch from SPEI-6 to SPEI-{best_by_sig['spei_value']}")
else:
    print(f"\n✅ SPEI-6 already has the most significant sites. Keep using SPEI-6.")

# Print summary of best timescale
print("\n" + "="*80)
print("SUMMARY OF BEST TIMESCALE")
print("="*80)
print(f"""
Based on the analysis of {len(shared_sites)} sites from strict triple intersection:

BEST TIMESCALE: SPEI-{best_by_sig['spei_value']}
  • {best_by_sig['n_sig_p05']}/{best_by_sig['n_sites_analyzed']} sites significant at p<0.05 ({best_by_sig['pct_sig_p05']:.1f}%)
  • {best_by_sig['n_sig_p01']}/{best_by_sig['n_sites_analyzed']} sites significant at p<0.01 ({best_by_sig['pct_sig_p01']:.1f}%)
  • Median p-value: {best_by_sig['median_p_value']:.4f}
  • Median slope: {best_by_sig['median_slope']:+.4f}

RECOMMENDATION: {'Switch to SPEI-' + str(best_by_sig['spei_value']) if best_by_sig['spei_timescale'] != 'SPEI_6' else 'Keep using SPEI-6'}
""")

print("="*80)
print("DIAGNOSTIC COMPLETE")
print("="*80)
print(f"\n📁 Results saved to: {diagnostic_dir}")