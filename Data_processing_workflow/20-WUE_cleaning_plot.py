# -*- coding: utf-8 -*-
"""
Outlier Removal for WUE Metrics Only - CONDITIONAL FILTERING
Applies per-site IQR filtering ONLY for sites with sufficient data (n>=8)
Sites with low data (n<8) receive ONLY hard bounds (no IQR)

CORRECT ORDER:
1. Low ET filter (remove unrealistic values when ET is tiny)
2. IQR per-site outlier detection (ONLY if n >= 8)
3. Hard bounds after IQR (applied to ALL sites)
4. Validation

Outputs:
- Cleaned CSV files (SAME NAMES as original)
- IQR removal summary by site
- Hard bounds removal summary by site  
- Combined audit summary
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input paths
DATA_DIR = r"M:\Research\WUE_CUE\data_products"
MONTHLY_INPUT = os.path.join(DATA_DIR, "WUE_CUE_monthly_merged_indices.csv")
YEARLY_INPUT = os.path.join(DATA_DIR, "WUE_CUE_yearly_merged_indices.csv")

# Output paths (SAME NAMES - no downstream impact)
MONTHLY_OUTPUT = os.path.join(DATA_DIR, "WUE_CUE_monthly_merged_indices_clean.csv")
YEARLY_OUTPUT = os.path.join(DATA_DIR, "WUE_CUE_yearly_merged_indices_clean.csv")

# Summary outputs (SAME NAMES)
IQR_SUMMARY = os.path.join(DATA_DIR, "outlier_summary_iqr_by_site_clean.csv")
HARD_BOUNDS_SUMMARY = os.path.join(DATA_DIR, "outlier_summary_hard_bounds_after_iqr_by_site_clean.csv")
AUDIT_SUMMARY = os.path.join(DATA_DIR, "outlier_filter_audit_by_site_clean.csv")

# Plots directory
PLOTS_DIR = os.path.join(DATA_DIR, "plots_outliers")
os.makedirs(PLOTS_DIR, exist_ok=True)

# Variables to clean
METRICS_TO_CLEAN = [
    'WUE', 'WUE_tra', 'WUE_eva', 'CUE', 'Evap_ratio', 'Trans_ratio'
]

METRICS_FOR_PLOTS = ['WUE', 'WUE_tra', 'WUE_eva', 'Evap_ratio', 'Trans_ratio']

# Hard bounds (applied to ALL sites, regardless of sample size)
HARD_BOUNDS = {
    'WUE': {'min': 0, 'max': 15},
    'WUE_tra': {'min': 0, 'max': 50},
    'WUE_eva': {'min': 0, 'max': 50},
    'CUE': {'min': -10, 'max': 1.5},
    'Evap_ratio': {'min': 0, 'max': 1},
    'Trans_ratio': {'min': 0, 'max': 1}
}

# Climate label shortener
CLIMATE_SHORTEN = {
    'Humid Subtropical': 'Humid Subtrop.',
    'Humid Continental': 'Humid Cont.',
    'Marine West Coast': 'Marine WC',
    'Mediterranean': 'Medit.',
    'Subarctic': 'Subarctic',
    'Arid': 'Arid',
    'Semi-Arid': 'Semi-Arid',
    'Tropical': 'Tropical'
}

# IQR threshold: only apply IQR if site has at least this many records
MIN_RECORDS_FOR_IQR = 8
IQR_MULTIPLIER = 1.5

ET_THRESHOLDS = {
    'monthly': {'ET': 5},
    'yearly': {'ET': 50}
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def shorten_climate_label(label):
    if pd.isna(label):
        return label
    label = str(label).strip()
    return CLIMATE_SHORTEN.get(label, label)


def load_data(filepath, timescale, log_messages):
    """Load and prepare data"""
    try:
        df = pd.read_csv(filepath)
        log_messages.append(f"  Loaded {timescale}: {len(df)} rows")
    except Exception as e:
        log_messages.append(f"  ERROR loading {timescale}: {e}")
        return None, log_messages
    
    df.columns = df.columns.str.strip()
    
    for col in ['site_name', 'IGBP', 'climate', 'water_class']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    
    for col in METRICS_TO_CLEAN:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    df = df.replace([np.inf, -np.inf], np.nan)
    
    return df, log_messages


def apply_low_et_filter(df, timescale, log_messages, audit_dict, initial_counts):
    """Step 1: Remove values when ET is very low"""
    et_threshold = ET_THRESHOLDS[timescale]['ET']
    
    if 'ET' not in df.columns:
        return df, log_messages, audit_dict, initial_counts
    
    wue_metrics = ['WUE', 'WUE_tra', 'WUE_eva']
    existing_wue = [m for m in wue_metrics if m in df.columns]
    
    if not existing_wue:
        return df, log_messages, audit_dict, initial_counts
    
    low_et_mask = df['ET'] < et_threshold
    n_low_et = low_et_mask.sum()
    
    if n_low_et > 0:
        for metric in existing_wue:
            before = df[metric].notna().sum()
            df.loc[low_et_mask, metric] = np.nan
            after = df[metric].notna().sum()
            removed = before - after
            
            # For low ET, we need to track by site too
            for site in df['site_name'].unique():
                site_mask = df['site_name'] == site
                site_before = df.loc[site_mask & ~low_et_mask, metric].notna().sum() + df.loc[site_mask & low_et_mask, metric].notna().sum()
                site_after = df.loc[site_mask & ~low_et_mask, metric].notna().sum()
                site_removed = site_before - site_after
                
                if site_removed > 0:
                    key = (timescale, 'low_ET', metric, site)
                    audit_dict[key] = site_removed
                    
                    # Update initial counts for this site
                    init_key = (timescale, metric, site)
                    if init_key not in initial_counts:
                        initial_counts[init_key] = site_before
        
        log_messages.append(f"  Low ET filter ({et_threshold}): {n_low_et} rows affected")
    
    return df, log_messages, audit_dict, initial_counts


def apply_iqr_by_site_conditional(df, timescale, log_messages, audit_dict, iqr_records):
    """
    Step 2: Apply IQR outlier detection PER SITE
    ONLY if site has >= MIN_RECORDS_FOR_IQR records for that metric
    Sites with < MIN_RECORDS_FOR_IQR skip IQR (only hard bounds later)
    """
    
    for metric in METRICS_TO_CLEAN:
        if metric not in df.columns:
            continue
        
        sites = df['site_name'].unique()
        iqr_applied_count = 0
        iqr_skipped_count = 0
        
        for site in sites:
            site_mask = df['site_name'] == site
            site_data = df.loc[site_mask, metric].dropna()
            n_records = len(site_data)
            
            # Check if site has enough records for IQR
            if n_records < MIN_RECORDS_FOR_IQR:
                # Skip IQR for this site (log once per metric)
                if n_records > 0:
                    iqr_skipped_count += 1
                continue
            
            # Site has enough records - apply IQR
            iqr_applied_count += 1
            
            Q1 = site_data.quantile(0.25)
            Q3 = site_data.quantile(0.75)
            IQR = Q3 - Q1
            
            if IQR == 0:
                continue
            
            lower_bound = Q1 - IQR_MULTIPLIER * IQR
            upper_bound = Q3 + IQR_MULTIPLIER * IQR
            
            before = df.loc[site_mask, metric].notna().sum()
            
            outliers_mask = (df.loc[site_mask, metric] < lower_bound) | (df.loc[site_mask, metric] > upper_bound)
            n_outliers = outliers_mask.sum()
            
            if n_outliers > 0:
                df.loc[site_mask & outliers_mask, metric] = np.nan
                after = df.loc[site_mask, metric].notna().sum()
                removed = before - after
                
                # Update audit dictionary
                key = (timescale, 'IQR', metric, site)
                audit_dict[key] = removed
                
                # Record IQR summary
                iqr_records.append({
                    'timescale': timescale,
                    'site_name': site,
                    'variable': metric,
                    'n_before_iqr': n_records,
                    'n_removed_iqr': n_outliers,
                    'pct_removed_iqr': (n_outliers / n_records) * 100,
                    'Q1': Q1,
                    'Q3': Q3,
                    'IQR': IQR,
                    'iqr_lower_bound': lower_bound,
                    'iqr_upper_bound': upper_bound
                })
        
        # Log summary for this metric
        total_removed = sum(v for k, v in audit_dict.items() if len(k) == 4 and k[1] == 'IQR' and k[2] == metric)
        if total_removed > 0 and metric in METRICS_FOR_PLOTS:
            log_messages.append(f"  IQR {metric}: {total_removed} outliers removed across {iqr_applied_count} sites (skipped {iqr_skipped_count} sites with <{MIN_RECORDS_FOR_IQR} records)")
        elif iqr_skipped_count > 0:
            log_messages.append(f"  IQR {metric}: Skipped {iqr_skipped_count} sites with <{MIN_RECORDS_FOR_IQR} records (only hard bounds will apply)")
    
    return df, log_messages, audit_dict, iqr_records


def apply_hard_bounds_after_iqr(df, timescale, log_messages, audit_dict, hard_bounds_records):
    """
    Step 3: Apply hard bounds AFTER IQR as final safety check
    Applied to ALL sites regardless of sample size
    """
    
    for metric, bounds in HARD_BOUNDS.items():
        if metric not in df.columns:
            continue
        
        sites = df['site_name'].unique()
        
        for site in sites:
            site_mask = df['site_name'] == site
            
            # Count before hard bounds
            before = df.loc[site_mask, metric].notna().sum()
            
            if before == 0:
                continue
            
            # Count values outside bounds
            if bounds['min'] is not None:
                below_min = ((df.loc[site_mask, metric] < bounds['min']) & (df.loc[site_mask, metric].notna())).sum()
            else:
                below_min = 0
                
            if bounds['max'] is not None and bounds['max'] != np.inf:
                above_max = ((df.loc[site_mask, metric] > bounds['max']) & (df.loc[site_mask, metric].notna())).sum()
            else:
                above_max = 0
            
            n_removed = below_min + above_max
            
            if n_removed > 0:
                # Apply hard bounds
                if bounds['min'] is not None:
                    df.loc[site_mask & (df[metric] < bounds['min']), metric] = np.nan
                if bounds['max'] is not None and bounds['max'] != np.inf:
                    df.loc[site_mask & (df[metric] > bounds['max']), metric] = np.nan
                
                after = df.loc[site_mask, metric].notna().sum()
                removed = before - after
                
                # Update audit dictionary
                key = (timescale, 'hard_bounds', metric, site)
                audit_dict[key] = removed
                
                # Record hard bounds summary
                hard_bounds_records.append({
                    'timescale': timescale,
                    'site_name': site,
                    'variable': metric,
                    'n_before_hard_bounds_after_iqr': before,
                    'n_removed_below_min_after_iqr': below_min,
                    'n_removed_above_max_after_iqr': above_max,
                    'n_removed_total_hard_bounds_after_iqr': n_removed,
                    'pct_removed_hard_bounds_after_iqr': (n_removed / before) * 100,
                    'hard_min': bounds['min'],
                    'hard_max': bounds['max']
                })
        
        total_removed = sum(v for k, v in audit_dict.items() if len(k) == 4 and k[1] == 'hard_bounds' and k[2] == metric)
        if total_removed > 0 and metric in METRICS_FOR_PLOTS:
            log_messages.append(f"  Hard bounds after IQR {metric}: {total_removed} additional values removed")
    
    return df, log_messages, audit_dict, hard_bounds_records


def validate_final_bounds(df, timescale, log_messages):
    """Step 4: Final validation - ensure no WUE > 15 remains"""
    if 'WUE' in df.columns:
        wue_values = df['WUE'].dropna()
        if len(wue_values) > 0:
            max_wue = wue_values.max()
            if max_wue > 15:
                log_messages.append(f"  ⚠️ FINAL VALIDATION: Max WUE = {max_wue:.2f} exceeds 15!")
                offending = df[df['WUE'] > 15][['site_name', 'Year'] + (['month'] if timescale == 'monthly' else []) + ['WUE']]
                log_messages.append(f"  Offending rows:")
                for _, row in offending.iterrows():
                    if timescale == 'monthly':
                        log_messages.append(f"    {row['site_name']} {row['Year']}-{row['month']:02d}: {row['WUE']:.2f}")
                    else:
                        log_messages.append(f"    {row['site_name']} {row['Year']}: {row['WUE']:.2f}")
                df.loc[df['WUE'] > 15, 'WUE'] = np.nan
                log_messages.append(f"  ✅ Forced {len(offending)} values exceeding 15 to NaN")
            else:
                log_messages.append(f"  ✅ Final validation passed: max WUE = {max_wue:.2f} (≤15)")
    
    return df, log_messages


def create_plots(df_before, df_after, timescale, log_messages):
    """Create diagnostic plots"""
    
    try:
        # Per-site boxplots for each metric
        for metric in ['WUE', 'WUE_tra', 'WUE_eva']:
            if metric not in df_before.columns:
                continue
            
            fig, axes = plt.subplots(1, 2, figsize=(18, 8))
            all_sites = sorted(set(df_before['site_name'].unique()) | set(df_after['site_name'].unique()))
            
            before_data = [df_before[(df_before['site_name'] == site) & (df_before[metric].notna())][metric].values for site in all_sites]
            after_data = [df_after[(df_after['site_name'] == site) & (df_after[metric].notna())][metric].values for site in all_sites]
            
            # Before
            bp1 = axes[0].boxplot(before_data, labels=all_sites, patch_artist=True)
            for patch in bp1['boxes']:
                patch.set_facecolor('lightblue')
            axes[0].set_title(f'{metric} - BEFORE', fontsize=14, fontweight='bold')
            axes[0].tick_params(axis='x', rotation=90, labelsize=8)
            if metric == 'WUE':
                axes[0].axhline(y=15, color='r', linestyle='--', label='Upper bound (15)')
                axes[0].legend()
                axes[0].set_ylim(0, 20)
            
            # After
            bp2 = axes[1].boxplot(after_data, labels=all_sites, patch_artist=True)
            for patch in bp2['boxes']:
                patch.set_facecolor('lightgreen')
            axes[1].set_title(f'{metric} - AFTER (Hard Bounds + IQR for n≥{MIN_RECORDS_FOR_IQR})', fontsize=14, fontweight='bold')
            axes[1].tick_params(axis='x', rotation=90, labelsize=8)
            if metric == 'WUE':
                axes[1].axhline(y=15, color='r', linestyle='--', label='Upper bound (15)')
                axes[1].legend()
                axes[1].set_ylim(0, 20)
            
            plt.suptitle(f'{timescale.upper()} Data - {metric} Cleaning', fontsize=16, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(PLOTS_DIR, f'{timescale}_{metric}_boxplot_comparison.png'), dpi=150, bbox_inches='tight')
            plt.close()
            log_messages.append(f"  Created: {timescale}_{metric}_boxplot_comparison.png")
        
        # Histograms
        for metric in METRICS_FOR_PLOTS:
            if metric not in df_before.columns:
                continue
            
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            
            before_vals = df_before[metric].dropna()
            after_vals = df_after[metric].dropna()
            
            axes[0].hist(before_vals, bins=50, alpha=0.7, color='blue', edgecolor='black')
            axes[0].set_xlabel(metric)
            axes[0].set_ylabel('Frequency')
            axes[0].set_title(f'{metric} - BEFORE (n={len(before_vals)})')
            if metric == 'WUE':
                axes[0].axvline(x=15, color='r', linestyle='--', label='Upper bound (15)')
                axes[0].legend()
            
            axes[1].hist(after_vals, bins=50, alpha=0.7, color='green', edgecolor='black')
            axes[1].set_xlabel(metric)
            axes[1].set_ylabel('Frequency')
            axes[1].set_title(f'{metric} - AFTER (n={len(after_vals)})')
            if metric == 'WUE':
                axes[1].axvline(x=15, color='r', linestyle='--', label='Upper bound (15)')
                axes[1].legend()
            
            plt.suptitle(f'{timescale.upper()} Data - {metric} Distribution', fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(PLOTS_DIR, f'{timescale}_{metric}_histogram_comparison.png'), dpi=150, bbox_inches='tight')
            plt.close()
            log_messages.append(f"  Created: {timescale}_{metric}_histogram_comparison.png")
        
        # WUE vs ET scatter
        if 'WUE' in df_before.columns and 'ET' in df_before.columns:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            
            axes[0].scatter(df_before['ET'], df_before['WUE'], alpha=0.5, s=10, c='blue')
            axes[0].set_xlabel('ET (mm)')
            axes[0].set_ylabel('WUE')
            axes[0].set_title(f'{timescale.upper()} - BEFORE')
            axes[0].axhline(y=15, color='r', linestyle='--')
            axes[0].set_ylim(0, 25)
            
            axes[1].scatter(df_after['ET'], df_after['WUE'], alpha=0.5, s=10, c='green')
            axes[1].set_xlabel('ET (mm)')
            axes[1].set_ylabel('WUE')
            axes[1].set_title(f'{timescale.upper()} - AFTER')
            axes[1].axhline(y=15, color='r', linestyle='--')
            axes[1].set_ylim(0, 25)
            
            plt.tight_layout()
            plt.savefig(os.path.join(PLOTS_DIR, f'{timescale}_WUE_vs_ET.png'), dpi=150, bbox_inches='tight')
            plt.close()
            log_messages.append(f"  Created: {timescale}_WUE_vs_ET.png")
        
    except Exception as e:
        log_messages.append(f"  ERROR creating plots: {e}")
    
    return log_messages


def create_audit_summary(audit_dict, initial_counts, log_messages):
    """Create final audit summary combining all stages"""
    
    audit_records = []
    
    # Get all unique site-metric combinations
    all_keys = set()
    for (timescale, stage, metric, site) in audit_dict.keys():
        all_keys.add((timescale, metric, site))
    
    for (timescale, metric, site) in all_keys:
        initial = initial_counts.get((timescale, metric, site), 0)
        
        # Get removals by stage
        low_et = audit_dict.get((timescale, 'low_ET', metric, site), 0)
        iqr = audit_dict.get((timescale, 'IQR', metric, site), 0)
        hard = audit_dict.get((timescale, 'hard_bounds', metric, site), 0)
        
        final = initial - low_et - iqr - hard
        
        if initial > 0:
            audit_records.append({
                'timescale': timescale,
                'site_name': site,
                'variable': metric,
                'n_initial': initial,
                'n_removed_low_ET': low_et,
                'n_removed_iqr': iqr,
                'n_removed_hard_bounds_after_iqr': hard,
                'n_final': final,
                'pct_removed_total': ((initial - final) / initial * 100)
            })
    
    if audit_records:
        audit_df = pd.DataFrame(audit_records)
        audit_df = audit_df.sort_values(['timescale', 'site_name', 'variable'])
        audit_df.to_csv(AUDIT_SUMMARY, index=False)
        log_messages.append(f"\n✅ Saved audit summary: {AUDIT_SUMMARY}")
        log_messages.append(f"   Total records: {len(audit_df)}")
    else:
        log_messages.append(f"\n⚠️ No audit records generated")
        audit_df = pd.DataFrame()
    
    return audit_df, log_messages


def print_summary_stats(df, timescale, log_messages, audit_dict):
    """Print summary statistics after cleaning"""
    
    log_messages.append(f"\n{'='*60}")
    log_messages.append(f"SUMMARY AFTER CLEANING ({timescale.upper()})")
    log_messages.append(f"{'='*60}")
    
    for metric in METRICS_FOR_PLOTS:
        if metric in df.columns:
            values = df[metric].dropna()
            if len(values) > 0:
                log_messages.append(f"\n{metric}:")
                log_messages.append(f"  Mean: {values.mean():.3f}")
                log_messages.append(f"  Median: {values.median():.3f}")
                log_messages.append(f"  Max: {values.max():.3f}")
                log_messages.append(f"  N: {len(values)}")
    
    # Print sites that still needed hard bounds after IQR
    hard_bounds_sites = set()
    for (timescale_key, stage, metric, site), removed in audit_dict.items():
        if stage == 'hard_bounds' and removed > 0:
            hard_bounds_sites.add((site, metric))
    
    if hard_bounds_sites:
        log_messages.append(f"\n  Sites that needed hard bounds after IQR:")
        for site, metric in sorted(hard_bounds_sites):
            log_messages.append(f"    {site} - {metric}")


def process_file(filepath, output_path, timescale, log_messages):
    """Process a single file: low ET → IQR (conditional) → hard bounds → validation"""
    
    log_messages.append(f"\n{'='*60}")
    log_messages.append(f"Processing {timescale.upper()} file")
    log_messages.append(f"  IQR applied ONLY for sites with ≥{MIN_RECORDS_FOR_IQR} records")
    log_messages.append(f"  Hard bounds applied to ALL sites")
    log_messages.append(f"{'='*60}")
    
    # Load data
    df, log_messages = load_data(filepath, timescale, log_messages)
    if df is None:
        return None, log_messages
    
    # Save copy for plots
    df_before = df.copy()
    
    # Initialize tracking
    audit_dict = {}
    iqr_records = []
    hard_bounds_records = []
    initial_counts = {}
    
    # Record initial counts
    for metric in METRICS_TO_CLEAN:
        if metric in df.columns:
            for site in df['site_name'].unique():
                site_mask = df['site_name'] == site
                count = df.loc[site_mask, metric].notna().sum()
                if count > 0:
                    initial_counts[(timescale, metric, site)] = count
    
    # STEP 1: Low ET filter
    df, log_messages, audit_dict, initial_counts = apply_low_et_filter(df, timescale, log_messages, audit_dict, initial_counts)
    
    # STEP 2: IQR by site (CONDITIONAL - only if n >= MIN_RECORDS_FOR_IQR)
    df, log_messages, audit_dict, iqr_records = apply_iqr_by_site_conditional(df, timescale, log_messages, audit_dict, iqr_records)
    
    # STEP 3: Hard bounds after IQR (applied to ALL sites)
    df, log_messages, audit_dict, hard_bounds_records = apply_hard_bounds_after_iqr(df, timescale, log_messages, audit_dict, hard_bounds_records)
    
    # STEP 4: Final validation
    df, log_messages = validate_final_bounds(df, timescale, log_messages)
    
    # Save IQR summary
    if iqr_records:
        iqr_df = pd.DataFrame(iqr_records)
        iqr_df.to_csv(IQR_SUMMARY, index=False)
        log_messages.append(f"\n✅ Saved IQR summary: {IQR_SUMMARY}")
        log_messages.append(f"   Total records: {len(iqr_df)}")
    
    # Save hard bounds summary
    if hard_bounds_records:
        hard_df = pd.DataFrame(hard_bounds_records)
        hard_df.to_csv(HARD_BOUNDS_SUMMARY, index=False)
        log_messages.append(f"✅ Saved hard bounds summary: {HARD_BOUNDS_SUMMARY}")
        log_messages.append(f"   Total records: {len(hard_df)}")
    
    # Create audit summary
    audit_df, log_messages = create_audit_summary(audit_dict, initial_counts, log_messages)
    
    # Create plots
    log_messages = create_plots(df_before, df, timescale, log_messages)
    
    # Print summary statistics
    print_summary_stats(df, timescale, log_messages, audit_dict)
    
    # Save cleaned file (SAME NAME - no downstream impact)
    try:
        df.to_csv(output_path, index=False)
        log_messages.append(f"\n✅ Saved cleaned {timescale} file: {os.path.basename(output_path)}")
        log_messages.append(f"   Final rows: {len(df)}")
    except Exception as e:
        log_messages.append(f"❌ ERROR saving: {e}")
    
    return df, log_messages


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 80)
    print("WUE METRICS OUTLIER REMOVAL - CONDITIONAL FILTERING")
    print(f"IQR applied ONLY for sites with ≥{MIN_RECORDS_FOR_IQR} records")
    print("Hard bounds applied to ALL sites")
    print("Order: Low ET → IQR (conditional) → Hard Bounds → Validation")
    print("=" * 80)
    print(f"WUE hard bound: 0 < WUE ≤ 15 (applied to ALL sites)")
    print()
    
    log_messages = []
    log_messages.append(f"Outlier removal started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_messages.append(f"Processing order: Low ET → IQR (conditional) → Hard Bounds → Validation")
    log_messages.append(f"IQR threshold: n ≥ {MIN_RECORDS_FOR_IQR} records")
    log_messages.append(f"Hard bounds: WUE ≤ 15, WUE_tra ≤ 50, WUE_eva ≤ 50 (applied to ALL sites)")
    
    # Process monthly
    if os.path.exists(MONTHLY_INPUT):
        df_monthly, log_messages = process_file(MONTHLY_INPUT, MONTHLY_OUTPUT, 'monthly', log_messages)
    else:
        log_messages.append(f"❌ Monthly file not found")
    
    # Process yearly
    if os.path.exists(YEARLY_INPUT):
        df_yearly, log_messages = process_file(YEARLY_INPUT, YEARLY_OUTPUT, 'yearly', log_messages)
    else:
        log_messages.append(f"❌ Yearly file not found")
    
    # Final summary
    log_messages.append(f"\n{'='*60}")
    log_messages.append("PROCESSING COMPLETE")
    log_messages.append(f"{'='*60}")
    log_messages.append(f"Output files (SAME NAMES - no downstream impact):")
    log_messages.append(f"  Monthly: {MONTHLY_OUTPUT}")
    log_messages.append(f"  Yearly: {YEARLY_OUTPUT}")
    log_messages.append(f"  IQR Summary: {IQR_SUMMARY}")
    log_messages.append(f"  Hard Bounds Summary: {HARD_BOUNDS_SUMMARY}")
    log_messages.append(f"  Audit Summary: {AUDIT_SUMMARY}")
    
    # Save log
    log_file = os.path.join(PLOTS_DIR, "outlier_cleaning_log.txt")
    with open(log_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(log_messages))
    print(f"\n✅ Log saved to: {log_file}")
    
    print("\n" + "\n".join(log_messages[-50:]))


if __name__ == "__main__":
    main()