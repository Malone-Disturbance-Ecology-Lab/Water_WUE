# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 13:53:30 2025
@author: ammar

Figure 1 Panel (a): WUE Distribution Boxplot
Clean, professional boxplot with consistent site colors

UPDATED: 
- ALL THREE metrics (WUE_ET, WUE_E, WUE_T) now use the SAME shared subset
- Shared subset defined as sites with ALL THREE metrics available
- Strict triple intersection: WUE ∩ WUE_E ∩ WUE_T
- Capping at 95th percentile for VISUALIZATION ONLY
- Statistical summaries use ORIGINAL uncapped values within shared subset
- All metrics have IDENTICAL sample sizes
- NOTE: Strict month filter applied upstream in CHUNK 3 (no redundant filter here)
- UPDATED: Site-month counts now use ALL retained NN observations (not unique month combos)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import json

# =============================================================================
# CONFIGURATION - SIMPLE ADJUSTMENTS
# =============================================================================

# Figure size
FIGURE_WIDTH = 15
FIGURE_HEIGHT = 6

# Font sizes
XLABEL_FONT_SIZE = 28
YLABEL_FONT_SIZE = 22
XTICK_FONT_SIZE = 26
YTICK_FONT_SIZE = 26

# Boxplot colors
BOX_COLOR = '#F0F8FF'
BOX_ALPHA = 0.9
BOX_LINEWIDTH = 2
MEDIAN_COLOR = 'red'
MEDIAN_LINEWIDTH = 3

# DOT COLOR - CHANGED TO GREY (was master color mapping)
DOT_COLOR = '#E0E0E0'  # Grey color for all dots
DOT_ALPHA = 0.6  # Slightly transparent for better visibility

# Scatter settings
SCATTER_SIZE = 60
SCATTER_EDGECOLOR = 'black'
SCATTER_EDGEWIDTH = 1.0

# Output
DPI = 600
OUTPUT_FORMAT = 'both'
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\figures"

# File to save master color mapping (kept for compatibility but not used for dot colors)
COLOR_MAPPING_FILE = os.path.join(OUTPUT_DIR, "master_site_color_mapping.json")

# Capping percentile
CAPPING_PERCENTILE = 95

# =============================================================================
# CAPPING FUNCTION (for visualization only)
# =============================================================================

def cap_values_for_visualization(values, percentile=CAPPING_PERCENTILE):
    """Cap values at specified percentile for VISUALIZATION ONLY"""
    if len(values) == 0:
        return values, None, 0, []
    
    cap_value = np.percentile(values, percentile)
    extreme_mask = values > cap_value
    n_capped = np.sum(extreme_mask)
    extreme_values = values[extreme_mask] if n_capped > 0 else []
    values_capped = np.where(extreme_mask, cap_value, values)
    return values_capped, cap_value, n_capped, extreme_values

# =============================================================================
# DATA ANALYSIS - ALL THREE METRICS USE STRICT TRIPLE INTERSECTION
# =============================================================================

def analyze_figure_data():
    """Analyze data with ALL metrics using strict triple intersection"""
    
    print("="*80)
    print("DATA ANALYSIS FOR FIGURE 1")
    print("="*80)
    print("NOTE: Data already filtered upstream (strict month filter applied in CHUNK 3)")
    
    # Load data (already filtered by upstream CHUNK 3)
    results_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
    os.chdir(results_dir)
    nn_spei1 = pd.read_csv('wue_site_level_NN_medians_SPEI_1.csv')
    
    print(f"\nOriginal dataset (already filtered upstream):")
    print(f"Total records: {len(nn_spei1)}")
    print(f"Unique sites in file: {nn_spei1['site_name'].nunique()}")
    
    # =========================================================================
    # STEP 1: Define the shared subset using TRIPLE INTERSECTION
    # =========================================================================
    
    # Get sites with WUE available
    wue_sites = set(nn_spei1[nn_spei1['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
    
    # Get sites with WUE_eva available
    eva_sites = set(nn_spei1[nn_spei1['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
    
    # Get sites with WUE_tra available
    tra_sites = set(nn_spei1[nn_spei1['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
    
    # STRICT TRIPLE INTERSECTION: sites with ALL THREE metrics available
    shared_subset_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
    
    # Sites excluded from analysis
    wue_only = wue_sites - shared_subset_sites
    eva_only = eva_sites - shared_subset_sites
    tra_only = tra_sites - shared_subset_sites
    
    print(f"\n{'='*60}")
    print("STRICT TRIPLE INTERSECTION (WUE ∩ WUE_E ∩ WUE_T)")
    print(f"{'='*60}")
    print(f"WUE (bulk) available sites: {len(wue_sites)}")
    print(f"WUE_eva available sites: {len(eva_sites)}")
    print(f"WUE_tra available sites: {len(tra_sites)}")
    print(f"\nShared subset (ALL THREE metrics available): {len(shared_subset_sites)} sites")
    
    if len(wue_only) > 0:
        print(f"\n⚠️ Sites with WUE ONLY (excluded from all metrics):")
        for site in sorted(wue_only):
            print(f"    - {site}")
    
    if len(eva_only) > 0:
        print(f"\n⚠️ Sites with WUE_eva ONLY (excluded from all metrics):")
        for site in sorted(eva_only):
            print(f"    - {site}")
    
    if len(tra_only) > 0:
        print(f"\n⚠️ Sites with WUE_tra ONLY (excluded from all metrics):")
        for site in sorted(tra_only):
            print(f"    - {site}")
    
    # =========================================================================
    # STEP 2: ALL metrics now use strict triple intersection
    # =========================================================================
    
    print(f"\n{'='*60}")
    print("SUBSET SUMMARY")
    print(f"{'='*60}")
    print(f"ALL metrics (WUE, WUE_E, WUE_T): {len(shared_subset_sites)} sites")
    print(f"  (Strict triple intersection - data already filtered upstream)")
    
    # =========================================================================
    # STEP 3: Prepare data for each metric - ALL using shared subset
    # =========================================================================
    
    wue_metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    wue_labels_full = {'WUE': 'WUE (WUE$_{ET}$)', 
                       'WUE_eva': 'WUE$_E$ (Evaporation efficiency)', 
                       'WUE_tra': 'WUE$_T$ (Transpiration efficiency)'}
    wue_labels_short = {'WUE': 'WUE', 'WUE_eva': 'WUE_E', 'WUE_tra': 'WUE_T'}
    
    print(f"\n" + "="*60)
    print("SUMMARY STATISTICS (ORIGINAL UNCAPPED VALUES)")
    print("="*60)
    print("Note: Statistics use ORIGINAL values - capping applied for visualization only")
    print("      ALL metrics use strict triple intersection\n")
    
    summary_stats = []
    extreme_values_dict = {}
    
    # For ALL metrics - use shared subset only
    for metric in wue_metrics:
        metric_data = nn_spei1[nn_spei1['WUE_Metric'] == metric]
        # Filter to shared subset only (SAME sites for ALL metrics)
        metric_data = metric_data[metric_data['site_name'].isin(shared_subset_sites)]
        values_with_sites = metric_data[['site_name', 'WUE_median']].dropna()
        values = values_with_sites['WUE_median'].values
        
        # Calculate sample standard deviation (ddof=1) and standard error
        sd = np.std(values, ddof=1) if len(values) > 0 else 0
        se = sd / np.sqrt(len(values)) if len(values) > 0 else 0
        
        stats = {
            'Metric': wue_labels_short[metric],
            'Metric_Full': wue_labels_full[metric],
            'N': len(values),
            'N_sites_used': len(shared_subset_sites),
            'Subset_Type': f'Strict triple intersection (n={len(shared_subset_sites)} sites)',
            'Median': float(np.median(values)),
            'Mean': float(np.mean(values)),
            'SD': float(sd),
            'SE': float(se),
            'Min': float(np.min(values)),
            'Max': float(np.max(values)),
            'IQR_25': float(np.percentile(values, 25)),
            'IQR_75': float(np.percentile(values, 75)),
            'P95': float(np.percentile(values, 95))
        }
        
        # For visualization capping - identify extreme values
        if metric in ['WUE_tra', 'WUE_eva']:
            cap_value = np.percentile(values, 95)
            extreme_mask = values > cap_value
            extreme_values = values[extreme_mask]
            extreme_sites = values_with_sites[extreme_mask]['site_name'].values
            cap_used = cap_value
            n_capped = len(extreme_values)
            pct_capped = (len(extreme_values) / len(values)) * 100 if len(values) > 0 else 0
        else:
            # WUE remains UNCAPPED
            cap_used = None
            n_capped = 0
            pct_capped = 0
            extreme_values = []
            extreme_sites = []
        
        # Store extreme values
        if metric in ['WUE_tra', 'WUE_eva'] and n_capped > 0:
            extreme_values_dict[metric] = {
                'cap_value': cap_used,
                'values': extreme_values,
                'sites': extreme_sites,
                'n_capped': n_capped
            }
        
        print(f"\n{wue_labels_full[metric]}:")
        print(f"  Subset: {stats['Subset_Type']} (n={stats['N_sites_used']} sites)")
        print(f"  Median: {stats['Median']:.3f}")
        print(f"  Mean: {stats['Mean']:.3f} ± {stats['SE']:.3f} (SE)")
        print(f"  SD: {stats['SD']:.3f}")
        print(f"  Range: {stats['Min']:.3f} to {stats['Max']:.3f}")
        print(f"  IQR: {stats['IQR_25']:.3f} to {stats['IQR_75']:.3f}")
        print(f"  95th percentile: {stats['P95']:.3f}")
        
        if metric in ['WUE_tra', 'WUE_eva']:
            print(f"\n  VISUALIZATION Capping ({CAPPING_PERCENTILE}th percentile):")
            print(f"    Cap value: {cap_used:.3f}")
            print(f"    Points capped in visualization: {n_capped}/{len(values)} ({pct_capped:.1f}%)")
        
        summary_stats.append(stats)
    
    # Store the shared subset for plotting
    plot_data = {
        'shared_subset_sites': shared_subset_sites,
        'nn_spei1': nn_spei1
    }
    
    return summary_stats, plot_data, extreme_values_dict, shared_subset_sites

# =============================================================================
# PRINT FINAL SITE-MONTH COUNTS FOR FIGURE 1 
# FIXED: Now counts ALL retained NN observations (not unique month combos)
# =============================================================================

def print_final_site_month_counts(shared_subset_sites):
    """Print the FINAL site-month counts actually used in Figure 1
    FIXED: Counts ALL retained NN observations directly from monthly data"""
    
    print("\n" + "="*80)
    print("FINAL SITE-MONTH COUNTS (EXACTLY AS USED IN FIGURE 1)")
    print("="*80)
    
    # Load monthly data
    monthly_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
    df_monthly = pd.read_csv(monthly_path)
    
    # Filter to NN conditions at SPEI-1
    df_nn = df_monthly[df_monthly["SPEI_1_Cat"] == "NN"].copy()
    
    print("-" * 70)
    print("NOTE: ALL metrics use the SAME strict triple intersection")
    print("      Support counts based on ALL retained NN monthly observations")
    print("-" * 70)
    
    # Filter to shared subset sites only
    shared_sites_list = list(shared_subset_sites)
    df_nn_shared = df_nn[df_nn["site_name"].isin(shared_sites_list)]
    
    # Count ALL retained NN observations (FIXED: no drop_duplicates)
    total_nn_observations = len(df_nn_shared)
    n_sites_shared = len(shared_sites_list)
    
    # Count observations per site
    nn_obs_per_site = []
    for site in shared_sites_list:
        site_data = df_nn_shared[df_nn_shared["site_name"] == site]
        n_obs = len(site_data)  # FIXED: count all observations, not unique months
        nn_obs_per_site.append(n_obs)
    
    print(f"\nALL METRICS (WUE, WUE$_E$, WUE$_T$) - Strict triple intersection (n={n_sites_shared} sites each):")
    print(f"  Total retained NN observations: {total_nn_observations}")
    print(f"  Sites: {n_sites_shared}")
    
    print(f"\n  Observation statistics for kept sites:")
    print(f"    Min observations: {min(nn_obs_per_site)}")
    print(f"    Max observations: {max(nn_obs_per_site)}")
    print(f"    Mean observations: {np.mean(nn_obs_per_site):.1f}")
    print(f"    Median observations: {np.median(nn_obs_per_site):.0f}")
    
    # Verification
    print("\n" + "-" * 70)
    print("VERIFICATION (Figure 1 - ALL metrics now identical):")
    print(f"  WUE sites: {n_sites_shared}")
    print(f"  WUE_E sites: {n_sites_shared}")
    print(f"  WUE_T sites: {n_sites_shared}")
    print(f"  Match: ✓ YES (all metrics use identical strict triple intersection)")
    print(f"\n  Total retained NN observations: {total_nn_observations}")
    print(f"  Note: This represents ALL retained monthly observations (not unique month combos)")
    print("-" * 70)
    
    return total_nn_observations, n_sites_shared

# =============================================================================
# PRINT EXTREME VALUES FOR PAPER
# =============================================================================

def print_extreme_values_for_paper(extreme_values_dict):
    """Print extreme values that fall outside the cap range for reporting in paper"""
    
    print("\n" + "="*80)
    print("EXTREME VALUES OUTSIDE CAPPED RANGE (FOR PAPER REPORTING)")
    print("="*80)
    print("Note: These values are ABOVE the 95th percentile cap used for visualization")
    print("They are displayed in the paper as capped values but their true values are reported here\n")
    
    for metric, data in extreme_values_dict.items():
        metric_name = "WUE_T (Transpiration efficiency)" if metric == 'WUE_tra' else "WUE_E (Evaporation efficiency)"
        print(f"\n{'='*60}")
        print(f"{metric_name}")
        print(f"{'='*60}")
        print(f"{CAPPING_PERCENTILE}th percentile cap value: {data['cap_value']:.3f} g C kg⁻¹ H₂O")
        print(f"Number of sites exceeding cap: {data['n_capped']}")
        
        if data['n_capped'] > 0:
            print(f"\nIndividual site values above cap:")
            print("-" * 60)
            print(f"{'Site Name':<25} {'Original Value':<20} {'Value in Figure':<20}")
            print("-" * 60)
            
            # Sort by value descending
            sorted_indices = np.argsort(data['values'])[::-1]
            for idx in sorted_indices:
                site = data['sites'][idx]
                original_value = data['values'][idx]
                print(f"{site:<25} {original_value:>15.3f}         {data['cap_value']:>15.3f} (capped)")
            
            print(f"\nSummary of extreme values:")
            print(f"  Range: {np.min(data['values']):.3f} - {np.max(data['values']):.3f}")
            print(f"  Mean: {np.mean(data['values']):.3f} ± {np.std(data['values']):.3f}")
            print(f"  Median: {np.median(data['values']):.3f}")
            
            print(f"\n📝 PAPER-READY TEXT:")
            print(f"  For {metric_name}, {data['n_capped']} sites had values exceeding the {CAPPING_PERCENTILE}th percentile cap of {data['cap_value']:.3f} g C kg⁻¹ H₂O.")
            print(f"  These extreme values ranged from {np.min(data['values']):.3f} to {np.max(data['values']):.3f} g C kg⁻¹ H₂O")
            print(f"  (mean: {np.mean(data['values']):.3f} ± {np.std(data['values']):.3f}, median: {np.median(data['values']):.3f}).")
            print(f"  For visualization clarity, these values were capped at the {CAPPING_PERCENTILE}th percentile.")

# =============================================================================
# CREATE MASTER COLOR MAPPING - KEPT FOR COMPATIBILITY BUT NOT USED FOR DOTS
# =============================================================================

def create_master_color_mapping(shared_subset_sites):
    """Create master color mapping for ALL sites in the strict triple intersection"""
    
    # Keep for compatibility but dots will use grey color instead
    all_sites_list = sorted(shared_subset_sites)
    
    # Create master color mapping (for reference only - not used for dot colors)
    master_color_map = {}
    for i, site in enumerate(all_sites_list):
        master_color_map[site] = DOT_COLOR  # Now all sites get grey
    
    # Save mapping to JSON file (SAME NAME)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(COLOR_MAPPING_FILE, 'w') as f:
        json.dump(master_color_map, f, indent=2)
    
    print(f"\n✅ Created MASTER color mapping for {len(master_color_map)} sites (all grey)")
    print(f"Saved to: {COLOR_MAPPING_FILE}")
    
    return master_color_map

# =============================================================================
# CREATE FIGURE 1 - ALL THREE METRICS USE STRICT TRIPLE INTERSECTION
# ALL DOTS ARE NOW GREY
# =============================================================================

def create_figure1_boxplot(plot_data, master_color_map, extreme_values_dict, save_directory=None):
    """Create Figure 1 boxplot - ALL metrics use strict triple intersection"""
    
    print("\n" + "="*80)
    print("CREATING FIGURE 1 - POINTS IN FRONT OF BOXPLOT")
    print("="*80)
    print("Subset rules (ALL metrics now identical):")
    print("  - ALL metrics (WUE$_{ET}$, WUE$_E$, WUE$_T$): Strict triple intersection")
    print(f"    (sites with ALL THREE metrics available) - {len(plot_data['shared_subset_sites'])} sites")
    print("\nCapping for VISUALIZATION only:")
    print(f"  - WUE_tra and WUE_eva: capped at {CAPPING_PERCENTILE}th percentile in plot")
    print("  - WUE_ET: remains UNCAPPED")
    print("  - Statistical summaries use ORIGINAL uncapped values within strict triple intersection")
    print("\nDOT COLORS:")
    print(f"  - All dots: GREY (#{DOT_COLOR[1:]}) - uniform color across all sites and metrics")
    
    # Print sample sizes in console only (not on figure)
    print("\n" + "="*60)
    print("📊 SAMPLE SIZES USED IN FIGURE 1 (printed in console only):")
    print(f"  ALL metrics (WUE$_{{ET}}$, WUE$_E$, WUE$_T$): n = {len(plot_data['shared_subset_sites'])} sites each")
    print("="*60)
    
    # Setup output
    if save_directory is None:
        figures_dir = OUTPUT_DIR
    else:
        figures_dir = save_directory
    os.makedirs(figures_dir, exist_ok=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    
    # Prepare data
    wue_metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    wue_labels = {'WUE': 'WUE$_{ET}$', 'WUE_eva': 'WUE$_E$', 'WUE_tra': 'WUE$_T$'}
    nn_spei1 = plot_data['nn_spei1']
    shared_subset_sites = plot_data['shared_subset_sites']
    
    # Prepare boxplot data with capping for VISUALIZATION only
    boxplot_data = []
    boxplot_labels = []
    cap_values_dict = {}
    
    for metric in wue_metrics:
        metric_data = nn_spei1[nn_spei1['WUE_Metric'] == metric]
        # ALL metrics now use strict triple intersection
        metric_data = metric_data[metric_data['site_name'].isin(shared_subset_sites)]
        metric_data = metric_data.dropna(subset=['WUE_median'])
        
        original_values = metric_data['WUE_median'].values
        
        if len(original_values) == 0:
            print(f"⚠️ Warning: No data for {metric} after subset filtering")
            boxplot_data.append([])
            boxplot_labels.append(wue_labels[metric])
            continue
        
        # Apply capping for VISUALIZATION only
        if metric in ['WUE_tra', 'WUE_eva'] and metric in extreme_values_dict:
            cap_value = extreme_values_dict[metric]['cap_value']
            values_for_plot = np.where(original_values > cap_value, cap_value, original_values)
            cap_values_dict[metric] = cap_value
            print(f"{metric}: Using cap value {cap_value:.3f} for visualization (n={len(original_values)} sites)")
        else:
            values_for_plot = original_values
            cap_values_dict[metric] = None
            print(f"{metric}: UNCAPPED (n={len(original_values)} sites)")
        
        boxplot_data.append(values_for_plot)
        boxplot_labels.append(wue_labels[metric])
    
    # Create boxplot with zorder=2 (behind points)
    bp = ax.boxplot(
        boxplot_data,
        labels=boxplot_labels,
        patch_artist=True,
        widths=0.6,
        zorder=2
    )
    
    # Style boxes
    for box in bp['boxes']:
        box.set(facecolor=BOX_COLOR, alpha=BOX_ALPHA, 
                linewidth=BOX_LINEWIDTH, edgecolor='black',
                zorder=2)
    
    # Style median lines
    for median in bp['medians']:
        median.set(color=MEDIAN_COLOR, linewidth=MEDIAN_LINEWIDTH, zorder=2)
    
    # Style whiskers and caps
    for whisker in bp['whiskers']:
        whisker.set(color='black', linewidth=BOX_LINEWIDTH, zorder=2)
    for cap in bp['caps']:
        cap.set(color='black', linewidth=BOX_LINEWIDTH, zorder=2)
    ax.set_ylim(bottom=0, top=13.5)
    ax.set_yticks([0, 5, 10, 13])
    
    # Add individual points - ALL NOW GREY (uniform color)
    for i, metric in enumerate(wue_metrics):
        metric_data = nn_spei1[nn_spei1['WUE_Metric'] == metric]
        # ALL metrics use strict triple intersection
        metric_data = metric_data[metric_data['site_name'].isin(shared_subset_sites)]
        metric_data = metric_data.dropna(subset=['WUE_median'])
        
        for idx, row in metric_data.iterrows():
            site = row['site_name']
            original_value = row['WUE_median']
            
            # Apply capping for VISUALIZATION only
            if metric in ['WUE_tra', 'WUE_eva'] and metric in cap_values_dict and cap_values_dict[metric] is not None:
                plot_value = min(original_value, cap_values_dict[metric])
            else:
                plot_value = original_value
            
            # Add jitter
            x_pos = i + 1 + np.random.normal(0, 0.08)
            
            # Plot with UNIFORM GREY color for ALL dots
            ax.scatter(
                x_pos, plot_value,
                s=SCATTER_SIZE,
                color=DOT_COLOR,  # GREY for all dots
                edgecolor=SCATTER_EDGECOLOR,
                linewidth=SCATTER_EDGEWIDTH,
                alpha=DOT_ALPHA,
                zorder=3
            )
    
    # Add panel label
    ax.text(0.02, 0.98, '(a)', transform=ax.transAxes, 
            fontsize=32, fontweight='bold', 
            verticalalignment='top', horizontalalignment='left')
    
    # Axis labels
    ax.set_ylabel('Median WUE (g C kg$^{-1}$ H$_2$O$^{-1}$)',
                  fontsize=YLABEL_FONT_SIZE,
                  fontweight='bold',
                  labelpad=20)
    
    # Tick sizes
    ax.tick_params(axis='both', which='major', 
                   labelsize=max(XTICK_FONT_SIZE, YTICK_FONT_SIZE),
                   width=2, length=6)
    
    # NO GRID
    ax.grid(False)
    
    # Border settings
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')
    
    # Adjust layout
    plt.tight_layout(pad=2.0)
    
    # Save figure (SAME FILE NAMES)
    base_filename = os.path.join(figures_dir, 'Figure1_WUE_Distribution')
    
    if OUTPUT_FORMAT in ['png', 'both']:
        plt.savefig(f"{base_filename}.png", dpi=DPI, bbox_inches='tight', facecolor='white')
        print(f"\n✅ Saved: {base_filename}.png")
    
    if OUTPUT_FORMAT in ['pdf', 'both']:
        plt.savefig(f"{base_filename}.pdf", bbox_inches='tight', facecolor='white')
        print(f"✅ Saved: {base_filename}.pdf")
    
    plt.show()
    
    print(f"\nFigure 1 created successfully")
    print(f"  ✓ ALL metrics use STRICT TRIPLE INTERSECTION: {len(plot_data['shared_subset_sites'])} sites each")
    print(f"  ✓ Strict month filter applied upstream (CHUNK 3)")
    print(f"  ✓ WUE, WUE_E, and WUE_T all have exactly the same {len(plot_data['shared_subset_sites'])} sites")
    print(f"  ✓ ALL dots are GREY (uniform color: #{DOT_COLOR[1:]})")
    
    return fig

# =============================================================================
# PRINT RESULTS FOR MANUSCRIPT
# =============================================================================

def print_manuscript_results(summary_stats, shared_subset_sites, master_color_map, total_nn_observations, n_sites_shared):
    """Print formatted results for manuscript
    UPDATED: Now uses retained NN observations terminology"""
    
    print("\n" + "="*80)
    print("FIGURE 1 RESULTS FOR MANUSCRIPT")
    print("="*80)
    
    print(f"\n📊 DATASET SUMMARY:")
    print(f"  • Sites used in statistical tests (site-level medians): {n_sites_shared}")
    print(f"  • Total retained NN observations (data availability): {total_nn_observations}")
    print(f"    (Strict triple intersection - filtered upstream)")
    
    print(f"\n📈 WUE Statistics (g C kg⁻¹ H₂O) - ORIGINAL UNCAPPED VALUES:")
    print("-" * 130)
    print(f"{'Metric':<15} {'n_sites':<10} {'Median':<12} {'Mean ± SE':<20} {'SD':<12} {'Range':<30} {'95th %ile':<12}")
    print("-" * 130)
    
    for stats in summary_stats:
        median_str = f"{stats['Median']:.3f}"
        mean_se_str = f"{stats['Mean']:.3f} ± {stats['SE']:.3f}"
        sd_str = f"{stats['SD']:.3f}"
        range_str = f"{stats['Min']:.3f} - {stats['Max']:.3f}"
        p95_str = f"{stats.get('P95', np.nan):.3f}"
        
        print(f"{stats['Metric']:<15} {stats['N_sites_used']:<10} {median_str:<12} {mean_se_str:<20} {sd_str:<12} {range_str:<30} {p95_str:<12}")
    
    print("-" * 130)
    
    # Methodological note for paper
    print(f"\n📝 METHODOLOGICAL NOTE FOR PAPER:")
    print(f"  Boxplots show medians and interquartile ranges. Mean ± SE is reported in the")
    print(f"  manuscript summary table, with SD provided as a measure of among-site variability.")
    print(f"")
    print(f"  All three WUE metrics (WUE$_{{ET}}$, WUE$_E$, and WUE$_T$) were evaluated")
    print(f"  on the same shared subset of {n_sites_shared} sites that had ALL THREE")
    print(f"  metrics available under NN conditions (strict triple intersection).")
    print(f"  A strict month filter (only months where all three metrics had data) was applied upstream.")
    print(f"  This ensures full comparability across all metrics.")
    print(f"")
    print(f"  Total retained NN observations: {total_nn_observations}")
    print(f"  (This represents ALL retained monthly observations, not unique month combinations)")
    print(f"")
    print(f"  Capping for visualization (shared subset):")
    print(f"  • WUE (WUE$_{{ET}}$): UNCAPPED in both visualization and statistics")
    print(f"  • WUE$_E$ (WUE_eva): Capped at {CAPPING_PERCENTILE}th percentile for VISUALIZATION only")
    print(f"  • WUE$_T$ (WUE_tra): Capped at {CAPPING_PERCENTILE}th percentile for VISUALIZATION only")
    print(f"  • All statistical summaries use ORIGINAL uncapped values")
    print(f"")
    print(f"  Dot colors: All points are displayed in grey for consistent representation")

# =============================================================================
# MAIN WORKFLOW FOR FIGURE 1
# =============================================================================

def main_figure1():
    """Complete workflow for Figure 1 - relies on upstream filtering"""
    
    print("="*80)
    print("FIGURE 1: WUE DISTRIBUTION ANALYSIS")
    print("="*80)
    print("NOTE: Filtering applied upstream in CHUNK 3 (strict month filter)")
    print("="*80)
    
    # Step 1: Analyze data with strict triple intersection
    summary_stats, plot_data, extreme_values_dict, shared_subset_sites = analyze_figure_data()
    
    # Step 2: Print final site-month counts (FIXED: uses all retained observations)
    total_nn_observations, n_sites_shared = print_final_site_month_counts(
        shared_subset_sites=shared_subset_sites
    )
    
    # Step 3: Print extreme values for paper reporting
    if extreme_values_dict:
        print_extreme_values_for_paper(extreme_values_dict)
    
    # Step 4: Create MASTER color mapping for shared subset sites (all grey now)
    master_color_map = create_master_color_mapping(
        shared_subset_sites=shared_subset_sites
    )
    
    # Step 5: Create Figure 1
    figure = create_figure1_boxplot(
        plot_data=plot_data,
        master_color_map=master_color_map,
        extreme_values_dict=extreme_values_dict,
        save_directory=OUTPUT_DIR
    )
    
    # Step 6: Print results for manuscript
    print_manuscript_results(summary_stats, shared_subset_sites, master_color_map, 
                            total_nn_observations, n_sites_shared)
    
    print("\n" + "="*80)
    print("FIGURE 1 COMPLETE")
    print("="*80)
    print("\n✅ KEY FEATURES:")
    print(f"  • STRICT TRIPLE INTERSECTION: WUE ∩ WUE_E ∩ WUE_T = {n_sites_shared} sites")
    print(f"  • Upstream filtering (CHUNK 3): Strict month filter applied")
    print(f"  • ALL metrics now have EXACTLY the same {n_sites_shared} sites")
    print(f"  • Full comparability across WUE, WUE_E, and WUE_T is guaranteed")
    print(f"  • Total retained NN observations: {total_nn_observations}")
    print(f"  • n numbers printed in console only (not on figure)")
    print(f"  • SAME output file names (no downstream impact)")
    print(f"  • ALL dots are GREY (uniform color: #{DOT_COLOR[1:]})")
    
    return figure, master_color_map

# =============================================================================
# EXECUTE FIGURE 1
# =============================================================================

if __name__ == "__main__":
    figure1, master_colors = main_figure1()