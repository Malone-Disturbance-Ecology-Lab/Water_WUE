# -*- coding: utf-8 -*-
"""
DIAGNOSTIC ANALYSIS: WUE_T Response to Extreme Drought/Wet Conditions by Ecosystem Class
=====================================================================================
Test whether WUE_T responses under severe drought and wet periods differ among 
ecosystem/salinity classes using SPEI-48.

Hypotheses (H3):
1. Saline sites show smaller WUE_T declines than upland/freshwater sites under severe drought
2. Saline sites show positive WUE_T changes under wet conditions
3. Upland sites show more negative WUE_T changes under severe drought than saline sites

Filtering criteria:
- Triple intersection sites (WUE, WUE_eva, WUE_tra all available)
- Trans_ratio > 0 (no NaN, no zero)
- Minimum 3 near-normal months for baseline
====================================================================================
"""

import pandas as pd
import numpy as np
from scipy.stats import wilcoxon, kruskal, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# File paths
BASE_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results")
MONTHLY_FILE = BASE_DIR / "monthly_data_after_outlier_removal.csv"
NN_MEDIANS_FILE = BASE_DIR / "wue_site_level_NN_medians_SPEI_1.csv"
TET_FILE = BASE_DIR / "T_ET_ratio" / "TET_gradient_site_NN_medians.csv"
SLOPES_FILE = BASE_DIR / "Q3_analysis" / "Q3_WUET_SPEI48_slopes.csv"

# Output directory
OUTPUT_DIR = BASE_DIR / "diagnostic_extreme_conditions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Minimum samples for analysis
MIN_MONTHS_BASELINE = 3
MIN_SITES_PER_GROUP = 3
MIN_SITES_FOR_TEST = 5

# SPEI-48 thresholds
SPEI_CLASSES = {
    'ExtremeDry': (-np.inf, -2.0),
    'SevereDry': (-2.0, -1.5),
    'ModerateDry': (-1.5, -1.0),
    'NearNormal': (-0.5, 0.5),
    'ModerateWet': (1.0, 1.5),
    'SevereWet': (1.5, 2.0),
    'ExtremeWet': (2.0, np.inf)
}

# Broader classes for statistical power
BROAD_CLASSES = {
    'Dry_all': (-np.inf, -1.0),
    'SevereExtremeDry': (-np.inf, -1.5),
    'Wet_all': (1.0, np.inf),
    'SevereExtremeWet': (1.5, np.inf)
}

# Ecosystem order for plotting
ECOSYSTEM_ORDER = ['Upland', 'Freshwater', 'Brackish', 'Saline']
ECOSYSTEM_COLORS = {
    'Upland': '#800080',
    'Freshwater': '#0000FF',
    'Brackish': '#008080',
    'Saline': '#FFA500'
}

print("="*80)
print("DIAGNOSTIC: WUE_T Response to Extreme Drought/Wet Conditions")
print("="*80)

# ============================================================================
# STEP 1: LOAD ALL DATA
# ============================================================================

print("\n[STEP 1] Loading data...")

# Load monthly data
df_monthly = pd.read_csv(MONTHLY_FILE)
print(f"  Monthly data: {len(df_monthly):,} rows")

# Load NN medians for triple intersection
nn_data = pd.read_csv(NN_MEDIANS_FILE)
print(f"  NN medians data: {len(nn_data)} rows")

# Load TET file for salinity categories and Trans_ratio
df_tet = pd.read_csv(TET_FILE)
print(f"  TET data: {len(df_tet)} rows")

# Load SPEI-48 slopes file (as backup for salinity)
df_slopes = pd.read_csv(SLOPES_FILE)
print(f"  SPEI-48 slopes data: {len(df_slopes)} rows")

# ============================================================================
# STEP 2: GET TRIPLE INTERSECTION SITES (WUE, WUE_eva, WUE_tra all available)
# ============================================================================

print("\n[STEP 2] Getting triple intersection sites...")

wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
print(f"  Triple intersection sites: {len(shared_sites)}")

# ============================================================================
# STEP 3: FILTER MONTHLY DATA
# ============================================================================

print("\n[STEP 3] Filtering monthly data...")

# Filter to triple intersection sites
df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
print(f"  After triple intersection: {len(df_monthly):,} rows")

# Apply strict month filter (ALL three WUE metrics must have data)
strict_mask = (
    df_monthly['WUE'].notna() &
    df_monthly['WUE_eva'].notna() &
    df_monthly['WUE_tra'].notna()
)
df_monthly = df_monthly[strict_mask].copy()
print(f"  After strict month filter: {len(df_monthly):,} rows")

# Apply Trans_ratio > 0 filter (consistent with SPEI-48 workflow)
if 'Trans_ratio' in df_monthly.columns:
    initial_rows = len(df_monthly)
    df_monthly = df_monthly[df_monthly['Trans_ratio'] > 0].copy()
    print(f"  After Trans_ratio > 0 filter: {len(df_monthly):,} rows (removed {initial_rows - len(df_monthly)})")
else:
    print(f"  Warning: Trans_ratio column not found, skipping ratio filter")

# Keep only required columns
required_cols = ['site_name', 'year', 'month', 'WUE_tra', 'SPEI_48']
available_cols = [col for col in required_cols if col in df_monthly.columns]
df_monthly = df_monthly[available_cols].copy()
df_monthly = df_monthly.dropna(subset=['WUE_tra', 'SPEI_48'])
print(f"  After dropping missing WUE_tra/SPEI_48: {len(df_monthly):,} rows")

# ============================================================================
# STEP 4: GET SALINITY CATEGORIES (from TET file or slopes file)
# ============================================================================

print("\n[STEP 4] Getting Salinity Categories...")

# Try to get salinity from TET file first
if 'Salinity_Category' in df_tet.columns:
    salinity_df = df_tet[['site_name', 'Salinity_Category']].drop_duplicates()
    print(f"  Salinity from TET file: {len(salinity_df)} sites")
else:
    # Try from slopes file
    if 'Salinity_Category' in df_slopes.columns:
        salinity_df = df_slopes[['site_name', 'Salinity_Category']].drop_duplicates()
        print(f"  Salinity from slopes file: {len(salinity_df)} sites")
    else:
        print(f"  ERROR: Salinity_Category not found in any file")
        print(f"    TET columns: {df_tet.columns.tolist()}")
        print(f"    Slopes columns: {df_slopes.columns.tolist()}")
        exit()

# Merge salinity to monthly data
df_monthly = df_monthly.merge(salinity_df, on='site_name', how='left')
initial_sites = df_monthly['site_name'].nunique()
df_monthly = df_monthly.dropna(subset=['Salinity_Category'])
final_sites = df_monthly['site_name'].nunique()
print(f"  Sites with salinity data: {final_sites} (dropped {initial_sites - final_sites})")

# Keep only the 4 main ecosystem classes
df_monthly = df_monthly[df_monthly['Salinity_Category'].isin(ECOSYSTEM_ORDER)].copy()
print(f"  Sites after filtering to main classes: {df_monthly['site_name'].nunique()}")

# ============================================================================
# STEP 5: ASSIGN SPEI-48 CLASSES
# ============================================================================

print("\n[STEP 5] Assigning SPEI-48 classes...")

def assign_spei_class(spei_value):
    """Assign SPEI value to a class"""
    if spei_value <= -2.0:
        return 'ExtremeDry'
    elif spei_value <= -1.5:
        return 'SevereDry'
    elif spei_value <= -1.0:
        return 'ModerateDry'
    elif -0.5 <= spei_value <= 0.5:
        return 'NearNormal'
    elif 1.0 <= spei_value < 1.5:
        return 'ModerateWet'
    elif 1.5 <= spei_value < 2.0:
        return 'SevereWet'
    elif spei_value >= 2.0:
        return 'ExtremeWet'
    else:
        return 'Other'

df_monthly['spei_class'] = df_monthly['SPEI_48'].apply(assign_spei_class)

# Assign broad classes
df_monthly['dry_class'] = 'None'
df_monthly.loc[df_monthly['SPEI_48'] <= -1.0, 'dry_class'] = 'Dry_all'
df_monthly.loc[df_monthly['SPEI_48'] <= -1.5, 'dry_class'] = 'SevereExtremeDry'

df_monthly['wet_class'] = 'None'
df_monthly.loc[df_monthly['SPEI_48'] >= 1.0, 'wet_class'] = 'Wet_all'
df_monthly.loc[df_monthly['SPEI_48'] >= 1.5, 'wet_class'] = 'SevereExtremeWet'

print(f"\n  Class distribution:")
for class_name in ['ExtremeDry', 'SevereDry', 'ModerateDry', 'NearNormal', 
                   'ModerateWet', 'SevereWet', 'ExtremeWet']:
    count = len(df_monthly[df_monthly['spei_class'] == class_name])
    if count > 0:
        print(f"    {class_name}: {count} months")

# ============================================================================
# STEP 6: CALCULATE SITE-LEVEL NEAR-NORMAL BASELINES
# ============================================================================

print("\n[STEP 6] Calculating site-level near-normal WUE_T baselines...")

site_baselines = []
for site in df_monthly['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site]
    near_normal_data = site_data[site_data['spei_class'] == 'NearNormal']
    
    if len(near_normal_data) >= MIN_MONTHS_BASELINE:
        baseline_median = near_normal_data['WUE_tra'].median()
        site_baselines.append({
            'site_name': site,
            'baseline_wue_tra': baseline_median,
            'n_near_normal': len(near_normal_data)
        })
    else:
        print(f"  WARNING: {site} has insufficient near-normal months ({len(near_normal_data)})")

df_baseline = pd.DataFrame(site_baselines)
print(f"  Sites with valid baselines: {len(df_baseline)}")

# Merge baselines back to monthly data
df_monthly = df_monthly.merge(df_baseline, on='site_name', how='inner')
print(f"  Monthly data after baseline filter: {len(df_monthly):,} rows")

# Calculate anomalies
df_monthly['wue_anomaly'] = df_monthly['WUE_tra'] - df_monthly['baseline_wue_tra']
df_monthly['wue_pct_change'] = 100 * (df_monthly['WUE_tra'] - df_monthly['baseline_wue_tra']) / df_monthly['baseline_wue_tra']

# ============================================================================
# STEP 7: SITE-LEVEL SUMMARIES BY HYDROCLIMATE CLASS
# ============================================================================

print("\n[STEP 7] Calculating site-level summaries...")

def get_site_summary(df, condition_col, condition_value):
    """Calculate site-level median percent change for a given condition"""
    results = []
    for site in df['site_name'].unique():
        site_data = df[df['site_name'] == site]
        condition_data = site_data[site_data[condition_col] == condition_value]
        
        if len(condition_data) >= 2:  # At least 2 months for median
            eco = site_data['Salinity_Category'].iloc[0]
            median_pct = condition_data['wue_pct_change'].median()
            q25 = condition_data['wue_pct_change'].quantile(0.25)
            q75 = condition_data['wue_pct_change'].quantile(0.75)
            n_months = len(condition_data)
            
            results.append({
                'site_name': site,
                'Salinity_Category': eco,
                'condition': condition_value,
                'median_pct_change': median_pct,
                'q25': q25,
                'q75': q75,
                'n_months': n_months
            })
    return pd.DataFrame(results)

# Summarize for each condition
dry_summary = get_site_summary(df_monthly, 'dry_class', 'SevereExtremeDry')
dry_all_summary = get_site_summary(df_monthly, 'dry_class', 'Dry_all')
wet_summary = get_site_summary(df_monthly, 'wet_class', 'SevereExtremeWet')
wet_all_summary = get_site_summary(df_monthly, 'wet_class', 'Wet_all')

print(f"  SevereExtremeDry: {len(dry_summary)} sites")
print(f"  Dry_all: {len(dry_all_summary)} sites")
print(f"  SevereExtremeWet: {len(wet_summary)} sites")
print(f"  Wet_all: {len(wet_all_summary)} sites")

# ============================================================================
# STEP 8: STATISTICAL TESTS - WITHIN GROUP (DEVIATION FROM ZERO)
# ============================================================================

print("\n" + "="*80)
print("STATISTICAL TESTS: Within-group deviation from zero")
print("="*80)

def test_within_group(df_summary, condition_name, group_col='Salinity_Category'):
    """Test whether median percent change differs from zero for each group"""
    results = []
    for group in df_summary[group_col].unique():
        group_data = df_summary[df_summary[group_col] == group]
        
        if len(group_data) >= MIN_SITES_FOR_TEST:
            try:
                stat, p_val = wilcoxon(group_data['median_pct_change'])
                n_sites = len(group_data)
                median_val = group_data['median_pct_change'].median()
                q25_val = group_data['median_pct_change'].quantile(0.25)
                q75_val = group_data['median_pct_change'].quantile(0.75)
                n_positive = (group_data['median_pct_change'] > 0).sum()
                n_negative = (group_data['median_pct_change'] < 0).sum()
                
                results.append({
                    'Group': group,
                    'Condition': condition_name,
                    'N_sites': n_sites,
                    'Median_pct': median_val,
                    'IQR': f"{q25_val:.1f}–{q75_val:.1f}",
                    'N_positive': n_positive,
                    'N_negative': n_negative,
                    'p_value': p_val,
                    'significant': p_val < 0.05
                })
                
                sig_marker = "✓ SIGNIFICANT" if p_val < 0.05 else "✗ NOT significant"
                print(f"\n  {group} under {condition_name} (n={n_sites}):")
                print(f"    Median % change = {median_val:.1f}% [{q25_val:.1f}, {q75_val:.1f}]")
                print(f"    Positive: {n_positive}, Negative: {n_negative}")
                print(f"    Wilcoxon p = {p_val:.4f} → {sig_marker}")
                
            except Exception as e:
                print(f"  WARNING: Could not test {group} under {condition_name}: {e}")
        else:
            print(f"\n  {group} under {condition_name}: insufficient sites (n={len(group_data)} < {MIN_SITES_FOR_TEST})")
    
    return pd.DataFrame(results)

print("\n📊 SevereExtremeDry (SPEI-48 ≤ -1.5):")
dry_within_results = test_within_group(dry_summary, 'SevereExtremeDry')

print("\n📊 SevereExtremeWet (SPEI-48 ≥ 1.5):")
wet_within_results = test_within_group(wet_summary, 'SevereExtremeWet')

# ============================================================================
# STEP 9: STATISTICAL TESTS - SPECIFIC HYPOTHESES (H3)
# ============================================================================

print("\n" + "="*80)
print("HYPOTHESIS TESTS (H3) - Specific Comparisons")
print("="*80)

# H3.1: Saline vs Freshwater under severe drought
print("\n📊 H3.1: Saline vs Freshwater under SevereExtremeDry")

saline_dry = dry_summary[dry_summary['Salinity_Category'] == 'Saline']['median_pct_change'].values
freshwater_dry = dry_summary[dry_summary['Salinity_Category'] == 'Freshwater']['median_pct_change'].values
upland_dry = dry_summary[dry_summary['Salinity_Category'] == 'Upland']['median_pct_change'].values
brackish_dry = dry_summary[dry_summary['Salinity_Category'] == 'Brackish']['median_pct_change'].values

h3_results = {}

if len(saline_dry) >= MIN_SITES_FOR_TEST and len(freshwater_dry) >= MIN_SITES_FOR_TEST:
    u_stat, p_val = mannwhitneyu(saline_dry, freshwater_dry, alternative='two-sided')
    n1, n2 = len(saline_dry), len(freshwater_dry)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    saline_median = np.median(saline_dry)
    freshwater_median = np.median(freshwater_dry)
    
    h3_results['saline_vs_freshwater'] = {
        'saline_n': n1, 'saline_median': saline_median,
        'freshwater_n': n2, 'freshwater_median': freshwater_median,
        'u_stat': u_stat, 'p_value': p_val, 'r': r_biserial
    }
    
    print(f"   Saline: n={n1}, median % change = {saline_median:.1f}%")
    print(f"   Freshwater: n={n2}, median % change = {freshwater_median:.1f}%")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {p_val:.4f}, r = {r_biserial:.3f}")
    
    if p_val < 0.05:
        if saline_median > freshwater_median:
            print(f"   ✓ SUPPORTS: Saline sites show SMALLER declines (less negative) than Freshwater")
        else:
            print(f"   ✗ DOES NOT support: Saline declines are not smaller")
    else:
        print(f"   ✗ No significant difference detected")

# Saline vs Upland
if len(saline_dry) >= MIN_SITES_FOR_TEST and len(upland_dry) >= MIN_SITES_FOR_TEST:
    u_stat, p_val = mannwhitneyu(saline_dry, upland_dry, alternative='two-sided')
    n1, n2 = len(saline_dry), len(upland_dry)
    r_biserial = 1 - (2 * u_stat) / (n1 * n2)
    
    upland_median = np.median(upland_dry)
    
    h3_results['saline_vs_upland'] = {
        'saline_n': n1, 'saline_median': saline_median,
        'upland_n': n2, 'upland_median': upland_median,
        'u_stat': u_stat, 'p_value': p_val, 'r': r_biserial
    }
    
    print(f"\n   Saline vs Upland:")
    print(f"   Saline: n={n1}, median = {saline_median:.1f}%")
    print(f"   Upland: n={n2}, median = {upland_median:.1f}%")
    print(f"   Mann-Whitney U = {u_stat:.1f}, p = {p_val:.4f}, r = {r_biserial:.3f}")
    
    if p_val < 0.05 and saline_median > upland_median:
        print(f"   ✓ SUPPORTS: Saline shows SMALLER declines than Upland")
    elif p_val < 0.05:
        print(f"   ⚠️ Significant but opposite direction")
    else:
        print(f"   ✗ No significant difference")

# H3.2: Saline sites under wet conditions
print("\n📊 H3.2: Saline sites under wet conditions")

saline_wet = wet_all_summary[wet_all_summary['Salinity_Category'] == 'Saline']['median_pct_change'].values
saline_severe_wet = wet_summary[wet_summary['Salinity_Category'] == 'Saline']['median_pct_change'].values

if len(saline_wet) >= MIN_SITES_FOR_TEST:
    median_wet = np.median(saline_wet)
    n_pos = (saline_wet > 0).sum()
    n_neg = (saline_wet < 0).sum()
    
    h3_results['saline_wet'] = {
        'n': len(saline_wet),
        'median': median_wet,
        'n_positive': n_pos,
        'n_negative': n_neg
    }
    
    print(f"   Wet_all (SPEI-48 ≥ 1.0): n={len(saline_wet)}")
    print(f"     Median % change = {median_wet:.1f}%")
    print(f"     Positive: {n_pos}, Negative: {n_neg}")
    
    if median_wet > 0:
        print(f"     ✓ Median POSITIVE - consistent with hypothesis")
    else:
        print(f"     ✗ Median NOT positive")
    
    # Test if median > 0
    if len(saline_wet) >= MIN_SITES_FOR_TEST:
        stat, p_val = wilcoxon(saline_wet)
        print(f"     Wilcoxon test (H0: median=0): p = {p_val:.4f}")
        h3_results['saline_wet_p'] = p_val

if len(saline_severe_wet) >= MIN_SITES_FOR_TEST:
    median_severe = np.median(saline_severe_wet)
    n_pos = (saline_severe_wet > 0).sum()
    n_neg = (saline_severe_wet < 0).sum()
    
    h3_results['saline_severe_wet'] = {
        'n': len(saline_severe_wet),
        'median': median_severe,
        'n_positive': n_pos,
        'n_negative': n_neg
    }
    
    print(f"\n   SevereExtremeWet (SPEI-48 ≥ 1.5): n={len(saline_severe_wet)}")
    print(f"     Median % change = {median_severe:.1f}%")
    print(f"     Positive: {n_pos}, Negative: {n_neg}")

# ============================================================================
# STEP 10: CREATE FIGURES
# ============================================================================

print("\n[STEP 10] Creating figures...")

# Figure 1: Site-level median percent change by hydroclimate class
fig1, ax1 = plt.subplots(1, 1, figsize=(14, 8))

# Combine data for plotting
all_conditions = []
for cond in ['SevereExtremeDry', 'Dry_all', 'Wet_all', 'SevereExtremeWet']:
    if cond == 'SevereExtremeDry':
        df_plot = dry_summary.copy()
    elif cond == 'Dry_all':
        df_plot = dry_all_summary.copy()
    elif cond == 'Wet_all':
        df_plot = wet_all_summary.copy()
    elif cond == 'SevereExtremeWet':
        df_plot = wet_summary.copy()
    else:
        continue
    
    if len(df_plot) > 0:
        df_plot['condition'] = cond
        all_conditions.append(df_plot)

if all_conditions:
    df_plot_all = pd.concat(all_conditions, ignore_index=True)
    
    # Order conditions from dry to wet
    condition_order = ['SevereExtremeDry', 'Dry_all', 'Wet_all', 'SevereExtremeWet']
    
    # Create point plot
    for i, eco in enumerate(ECOSYSTEM_ORDER):
        eco_data = df_plot_all[df_plot_all['Salinity_Category'] == eco]
        if len(eco_data) > 0:
            x_positions = []
            for cond in condition_order:
                cond_data = eco_data[eco_data['condition'] == cond]
                if len(cond_data) > 0:
                    x_base = condition_order.index(cond)
                    jitter = np.random.normal(0, 0.1, len(cond_data))
                    x_positions.extend([x_base + j for j in jitter])
                    ax1.scatter([x_base + j for j in jitter], cond_data['median_pct_change'],
                              color=ECOSYSTEM_COLORS[eco], s=80, alpha=0.7, 
                              edgecolors='black', linewidth=0.5, label=eco if i==0 else "")
    
    ax1.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax1.set_xticks(range(len(condition_order)))
    ax1.set_xticklabels(condition_order, fontsize=12, rotation=45)
    ax1.set_ylabel('WUE$_T$ Percent Change from Baseline (%)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Hydroclimate Condition (SPEI-48)', fontsize=14, fontweight='bold')
    ax1.set_title('WUE$_T$ Response to Extreme Conditions by Ecosystem Class', fontsize=16, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_axisbelow(True)
    
    plt.tight_layout()
    fig1.savefig(OUTPUT_DIR / 'Figure_H3_extreme_conditions_response.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved Figure 1")

# Figure 2: Focused comparison
fig2, axes = plt.subplots(1, 2, figsize=(16, 8))

for idx, (condition, df_summ, title) in enumerate([
    ('SevereExtremeDry', dry_summary, 'Severe Drought (SPEI-48 ≤ -1.5)'),
    ('Wet_all', wet_all_summary, 'Wet Conditions (SPEI-48 ≥ 1.0)')
]):
    if len(df_summ) > 0:
        ax = axes[idx]
        
        eco_order = [e for e in ECOSYSTEM_ORDER if e in df_summ['Salinity_Category'].unique()]
        box_data = [df_summ[df_summ['Salinity_Category'] == eco]['median_pct_change'].values for eco in eco_order]
        
        bp = ax.boxplot(box_data, positions=range(len(eco_order)), widths=0.6,
                        patch_artist=True, showfliers=False)
        
        for i, patch in enumerate(bp['boxes']):
            patch.set_facecolor(ECOSYSTEM_COLORS[eco_order[i]])
            patch.set_alpha(0.7)
            patch.set_edgecolor('black')
        
        for i, eco in enumerate(eco_order):
            eco_data = df_summ[df_summ['Salinity_Category'] == eco]
            x_jitter = np.random.normal(i, 0.08, len(eco_data))
            ax.scatter(x_jitter, eco_data['median_pct_change'], 
                      color=ECOSYSTEM_COLORS[eco], s=60, alpha=0.6, 
                      edgecolors='black', linewidth=0.5, zorder=3)
        
        ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_xticks(range(len(eco_order)))
        ax.set_xticklabels(eco_order, fontsize=12, rotation=45)
        ax.set_ylabel('WUE$_T$ Percent Change (%)', fontsize=14, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

plt.tight_layout()
fig2.savefig(OUTPUT_DIR / 'Figure_H3_focused_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved Figure 2")

# ============================================================================
# STEP 11: MANUSCRIPT-READY INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("MANUSCRIPT-READY INTERPRETATION")
print("="*80)

# Generate interpretation
interpretation = f"""
DIAGNOSTIC ANALYSIS: WUE_T Response to Extreme Conditions

DATA SUMMARY:
  • Sites analyzed: {len(df_baseline)}
  • Near-normal baseline months per site: ≥{MIN_MONTHS_BASELINE}
  • SPEI-48 thresholds: 
    - Severe drought: ≤ -1.5
    - Wet conditions: ≥ 1.0

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

H3.1: Saline vs Freshwater/Upland under severe drought
"""

if 'saline_vs_freshwater' in h3_results:
    r = h3_results['saline_vs_freshwater']
    interpretation += f"""
   Saline (n={r['saline_n']}): median = {r['saline_median']:.1f}%
   Freshwater (n={r['freshwater_n']}): median = {r['freshwater_median']:.1f}%
   Mann-Whitney U = {r['u_stat']:.1f}, p = {r['p_value']:.4f}, r = {r['r']:.3f}
   
   {'✓ Sites with higher salinity showed less negative WUE_T responses during severe drought, consistent with the hypothesis that salt-adapted coastal wetlands maintain greater functional resilience to water limitation.' if r['p_value'] < 0.05 and r['saline_median'] > r['freshwater_median'] else 'No significant difference in WUE_T response was detected between saline and freshwater sites under severe drought.'}
"""

if 'saline_vs_upland' in h3_results:
    r = h3_results['saline_vs_upland']
    interpretation += f"""
   
   Saline (n={r['saline_n']}): median = {r['saline_median']:.1f}%
   Upland (n={r['upland_n']}): median = {r['upland_median']:.1f}%
   Mann-Whitney U = {r['u_stat']:.1f}, p = {r['p_value']:.4f}, r = {r['r']:.3f}
   
   {'✓ Upland sites showed more negative WUE_T declines than saline sites, consistent with expectations that upland ecosystems are more vulnerable to prolonged drought.' if r['p_value'] < 0.05 and r['upland_median'] < r['saline_median'] else 'No significant difference in WUE_T response was detected between upland and saline sites.'}
"""

interpretation += f"""

H3.2: Saline sites under wet conditions
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

if 'saline_wet' in h3_results:
    r = h3_results['saline_wet']
    interpretation += f"""
   Wet_all (SPEI-48 ≥ 1.0): n={r['n']}, median = {r['median']:.1f}%
   Positive: {r['n_positive']}, Negative: {r['n_negative']}
   {'✓ The positive median response is consistent with the hypothesis that saline sites increase WUE_T during pluvial periods.' if r['median'] > 0 else 'The median response was not positive, which does not support the hypothesis.'}
"""

if 'saline_wet_p' in h3_results:
    interpretation += f"""
   Wilcoxon test (H0: median = 0): p = {h3_results['saline_wet_p']:.4f}
"""

interpretation += f"""

LIMITATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   • Sample sizes varied by ecosystem class (Brackish had few sites)
   • Results are correlational and do not imply causation
   • WUE_T responses integrate multiple physiological processes
   • Near-normal baseline definition (±0.5 SPEI) may influence anomaly calculations

SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   This analysis provides {'support for' if any(r.get('p_value', 1) < 0.05 for r in h3_results.values() if isinstance(r, dict) and 'p_value' in r) else 'limited support for'} 
   the hypothesis that ecosystem salinity class modulates WUE_T responses to extreme 
   hydroclimatic conditions. Differences in drought response between saline and 
   freshwater/upland sites were {'detected' if h3_results.get('saline_vs_freshwater', {}).get('p_value', 1) < 0.05 else 'not detected'}.
"""

print(interpretation)

# Save interpretation
with open(OUTPUT_DIR / 'H3_interpretation.txt', 'w') as f:
    f.write(interpretation)
print(f"\n✅ Interpretation saved to: {OUTPUT_DIR / 'H3_interpretation.txt'}")

# Save summary CSVs
dry_summary.to_csv(OUTPUT_DIR / 'summary_SevereExtremeDry.csv', index=False)
wet_summary.to_csv(OUTPUT_DIR / 'summary_SevereExtremeWet.csv', index=False)
dry_all_summary.to_csv(OUTPUT_DIR / 'summary_Dry_all.csv', index=False)
wet_all_summary.to_csv(OUTPUT_DIR / 'summary_Wet_all.csv', index=False)
print(f"  ✓ Saved summary CSV files")

print("\n" + "="*80)
print("DIAGNOSTIC COMPLETE")
print("="*80)
print(f"\n📁 Output directory: {OUTPUT_DIR}")
print("   - Figure_H3_extreme_conditions_response.png")
print("   - Figure_H3_focused_comparison.png")
print("   - summary_*.csv files")
print("   - H3_interpretation.txt")
print("="*80)