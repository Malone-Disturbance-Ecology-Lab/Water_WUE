# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 14:05:35 2025
Updated: Simplified version - filtering already done by CHUNK 3
- ALL THREE metrics use SAME strict triple intersection (WUE ∩ WUE_E ∩ WUE_T)
- No redundant filters (CHUNK 3 already applied per-metric ≥3 months)
- NO Y-LIMIT CONSTRAINTS - matplotlib auto-scales
- Now reports Mean ± SE and SD for group-level summary statistics
- UPDATED: Points now colored by ECOSYSTEM TYPE (same color as boxplot face)
@author: ammar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import json
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch

# =============================================================================
# CONFIGURATION - ADJUST THESE PARAMETERS ONLY
# =============================================================================

# Figure dimensions
FIGURE_WIDTH = 15  # Width for 3 categories
FIGURE_HEIGHT = 6  # Same as Figure 1 height

# Font sizes - ADJUSTED PROPORTIONALLY
X_TICK_LABEL_SIZE = 28    # X-axis tick labels (WUE, WUE_E, WUE_T)
Y_TICK_LABEL_SIZE = 24    # Y-axis tick labels (numbers)
Y_AXIS_LABEL_SIZE = 22    # Y-axis title
LEGEND_FONT_SIZE = 25     # Legend text size

# Boxplot styling
BOX_WIDTH = 0.6           # Width of boxes
BOX_BORDER_WIDTH = 2.5    # Thickness of box borders
WHISKER_WIDTH = 2.5       # Thickness of whiskers
MEDIAN_WIDTH = 3.5        # Thickness of median lines
POINT_SIZE = 80           # Size of individual points

# Spacing - for 3 categories per metric
METRIC_SPACING = 2.5      # Space between different WUE metrics
WITHIN_METRIC_SPACING = 0.8  # Space between Freshwater, Saline, and Upland

# Colors for BOXES and POINTS (3 categories)
SALINITY_BOX_COLORS = {
    'Freshwater': '#0000FF',  # Blue for FRESHWATER BOXES and POINTS
    'Saline': '#FFA500',       # Orange for SALINE BOXES and POINTS
    'Upland': '#800080'        # Purple for UPLAND BOXES and POINTS
}

# Capping for visualization (SAME as Figure 1)
USE_CAPPING = True
CAPPING_PERCENTILE = 95  # Same as Figure 1

# Y-axis limits - SET TO None for auto-scaling (no constraints)
Y_AXIS_MIN = None  # Auto-scale minimum
Y_AXIS_MAX = None  # Auto-scale maximum

# File with MASTER color mapping (created by Figure 1) - NOT USED for points anymore
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\figures"
COLOR_MAPPING_FILE = os.path.join(OUTPUT_DIR, "master_site_color_mapping.json")

# Output settings
DPI = 600
OUTPUT_FORMAT = 'both'  # 'png', 'pdf', or 'both'

# Allowed salinity categories for Figure 2 (Brackish excluded)
ALLOWED_CATEGORIES = ['Freshwater', 'Saline', 'Upland']

# =============================================================================
# CAPPING FUNCTION (for visualization only - SAME as Figure 1)
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
# LOAD MASTER COLOR MAPPING FROM FIGURE 1 (For reference only - not used for coloring)
# =============================================================================

def load_master_color_mapping():
    """Load master color mapping created by Figure 1 (for reference only)"""
    
    if not os.path.exists(COLOR_MAPPING_FILE):
        print(f"❌ ERROR: Master color mapping not found!")
        print(f"Please run Figure 1 workflow first to create: {COLOR_MAPPING_FILE}")
        raise FileNotFoundError(f"Master color mapping not found: {COLOR_MAPPING_FILE}")
    
    with open(COLOR_MAPPING_FILE, 'r') as f:
        master_color_map = json.load(f)
    
    print(f"✅ Loaded MASTER color mapping with {len(master_color_map)} sites")
    print(f"File: {COLOR_MAPPING_FILE}")
    print(f"NOTE: Points will be colored by ECOSYSTEM TYPE, NOT by site")
    
    return master_color_map

# =============================================================================
# DATA LOADING FUNCTIONS - SIMPLIFIED (filtering already done by CHUNK 3)
# =============================================================================

def load_and_prepare_data(master_color_map):
    """
    Load data - filtering already applied by CHUNK 3 (per-metric ≥3 months)
    Figure 2 only applies salinity filter and triple intersection
    """
    BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
    
    # Load data (already filtered by CHUNK 3: each metric ≥3 months NN)
    file_path = os.path.join(BASE_DIR, "wue_site_level_summary_SPEI_1.csv")
    data = pd.read_csv(file_path)
    
    print(f"\n[Data Loading] Loaded {len(data)} rows from {file_path}")
    print(f"  Note: CHUNK 3 already applied per-metric ≥3 months NN filter")
    
    # Step 1: Filter to near-normal conditions (SPEI_Class == 'NN')
    nn_data = data[data['SPEI_Class'] == 'NN'].copy()
    print(f"  After SPEI_Class == 'NN': {len(nn_data)} rows")
    
    # Step 2: Filter to allowed salinity categories (Freshwater, Saline, Upland)
    coastal_data = nn_data[nn_data['Salinity_Category'].isin(ALLOWED_CATEGORIES)].copy()
    print(f"  After filtering to {ALLOWED_CATEGORIES}: {len(coastal_data)} rows")
    
    # Step 3: Get unique sites for each metric (for triple intersection)
    sites_with_wue = set(coastal_data[coastal_data['WUE_Metric'] == 'WUE']['site_name'].unique())
    sites_with_eva = set(coastal_data[coastal_data['WUE_Metric'] == 'WUE_eva']['site_name'].unique())
    sites_with_tra = set(coastal_data[coastal_data['WUE_Metric'] == 'WUE_tra']['site_name'].unique())
    
    # STRICT TRIPLE INTERSECTION: sites with ALL THREE metrics
    shared_sites = sites_with_wue.intersection(sites_with_eva).intersection(sites_with_tra)
    
    print(f"\n[Strict Triple Intersection - ALL METRICS]")
    print(f"  Sites with WUE available: {len(sites_with_wue)}")
    print(f"  Sites with WUE_E available: {len(sites_with_eva)}")
    print(f"  Sites with WUE_T available: {len(sites_with_tra)}")
    print(f"  Strict triple intersection (ALL metrics): {len(shared_sites)} sites")
    
    # Filter ALL metrics to strict triple intersection only
    final_data = coastal_data[coastal_data['site_name'].isin(shared_sites)].copy()
    
    print(f"\n[Final Dataset Summary - ALL metrics use strict triple intersection]")
    for metric in ['WUE', 'WUE_eva', 'WUE_tra']:
        metric_data = final_data[final_data['WUE_Metric'] == metric]
        print(f"  {metric}: {len(metric_data)} rows, {metric_data['site_name'].nunique()} sites")
    
    # Verify all metrics have same number of sites
    wue_sites_final = set(final_data[final_data['WUE_Metric'] == 'WUE']['site_name'].unique())
    eva_sites_final = set(final_data[final_data['WUE_Metric'] == 'WUE_eva']['site_name'].unique())
    tra_sites_final = set(final_data[final_data['WUE_Metric'] == 'WUE_tra']['site_name'].unique())
    
    if len(wue_sites_final) == len(eva_sites_final) == len(tra_sites_final) == len(shared_sites):
        print(f"\n✓ VERIFIED: ALL three metrics have identical sample size ({len(shared_sites)} sites each)")
    else:
        print(f"\n⚠️ WARNING: Sample sizes differ!")
        print(f"  WUE: {len(wue_sites_final)} sites")
        print(f"  WUE_E: {len(eva_sites_final)} sites")
        print(f"  WUE_T: {len(tra_sites_final)} sites")
    
    # Print sample sizes by category for each metric
    print(f"\n[Sample Sizes by Category - Strict triple intersection]")
    for metric in ['WUE', 'WUE_eva', 'WUE_tra']:
        print(f"\n  {metric}:")
        metric_data = final_data[final_data['WUE_Metric'] == metric]
        for salinity in ALLOWED_CATEGORIES:
            n_sites = len(metric_data[metric_data['Salinity_Category'] == salinity]['site_name'].unique())
            print(f"    {salinity}: {n_sites} sites")
    
    return final_data, shared_sites

# =============================================================================
# STATISTICAL TESTING FUNCTIONS - ALL METRICS USE STRICT TRIPLE INTERSECTION
# =============================================================================

def perform_statistical_tests(data, shared_sites):
    """Perform Kruskal-Wallis tests with post-hoc pairwise comparisons using ORIGINAL values
       ALL metrics now use the same strict triple intersection"""
    
    print("\n" + "="*80)
    print("STATISTICAL ANALYSIS (3 Groups: Freshwater, Saline, Upland)")
    print("NOTE: Using ORIGINAL uncapped values (capping applied for visualization only)")
    print("      ALL metrics use the same strict triple intersection")
    print("="*80)
    
    metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    salinities = ['Upland', 'Freshwater', 'Saline']
    results = {}
    
    for metric in metrics:
        print(f"\n{'='*60}")
        print(f"METRIC: {metric}")
        print(f"{'='*60}")
        
        # Get data for this metric - ALL now use strict triple intersection
        metric_data = data[data['WUE_Metric'] == metric].copy()
        print(f"  Using strict triple intersection ({len(shared_sites)} sites)")
        print(f"  Statistical test uses ORIGINAL uncapped values")
        
        # Extract ORIGINAL values for each salinity category (no capping for stats)
        groups = {}
        group_names = []
        group_values = []
        
        for salinity in salinities:
            vals = metric_data[metric_data['Salinity_Category'] == salinity]['Median'].dropna()
            groups[salinity] = vals
            group_names.append(salinity)
            group_values.append(vals)
            print(f"  {salinity}: n={len(vals)} sites")
        
        # Kruskal-Wallis test (only include groups with data)
        non_empty_groups = [(name, vals) for name, vals in groups.items() if len(vals) > 0]
        if len(non_empty_groups) >= 2:
            names, values_lists = zip(*non_empty_groups)
            h_stat, p_value = kruskal(*values_lists)
            print(f"\n  Kruskal-Wallis test (using ORIGINAL uncapped values):")
            print(f"    H-statistic = {h_stat:.3f}")
            print(f"    p-value = {p_value:.6f}")
            
            if p_value < 0.05:
                print(f"    ✓ Significant differences detected (p < 0.05)")
                
                # Post-hoc pairwise Mann-Whitney U tests with Bonferroni correction
                pairwise_results = []
                
                for i in range(len(names)):
                    for j in range(i+1, len(names)):
                        group1 = names[i]
                        group2 = names[j]
                        vals1 = groups[group1]
                        vals2 = groups[group2]
                        
                        if len(vals1) > 0 and len(vals2) > 0:
                            u_stat, p_pair = mannwhitneyu(vals1, vals2, alternative='two-sided')
                            pairwise_results.append({
                                'comparison': f"{group1} vs {group2}",
                                'group1': group1,
                                'group2': group2,
                                'u_stat': u_stat,
                                'p_value': p_pair,
                                'n1': len(vals1),
                                'n2': len(vals2),
                                'median1': np.median(vals1),
                                'median2': np.median(vals2)
                            })
                
                # Apply Bonferroni correction
                if pairwise_results:
                    p_values = [res['p_value'] for res in pairwise_results]
                    reject, p_corrected, _, _ = multipletests(p_values, method='bonferroni')
                    
                    print(f"\n  Pairwise comparisons (Bonferroni-corrected):")
                    for idx, res in enumerate(pairwise_results):
                        sig_symbol = '***' if p_corrected[idx] < 0.001 else '**' if p_corrected[idx] < 0.01 else '*' if p_corrected[idx] < 0.05 else 'ns'
                        print(f"    {res['comparison']}: U={res['u_stat']:.1f}, p_raw={res['p_value']:.4f}, p_corr={p_corrected[idx]:.4f} {sig_symbol}")
                        print(f"       Median {res['group1']}: {res['median1']:.3f} (n={res['n1']})")
                        print(f"       Median {res['group2']}: {res['median2']:.3f} (n={res['n2']})")
            
            else:
                print(f"    No significant differences among groups (p >= 0.05)")
        else:
            print(f"\n  Insufficient data for Kruskal-Wallis test (need at least 2 groups with data)")
        
        results[metric] = groups
    
    return results

# =============================================================================
# FIGURE CREATION FUNCTIONS - ALL METRICS USE STRICT TRIPLE INTERSECTION
# =============================================================================

def create_figure_formatting(data, master_color_map, shared_sites):
    """Create SINGLE PANEL Figure 2 with capping for visualization - NO Y-LIMIT CONSTRAINTS
       ALL metrics now use the same strict triple intersection
       POINTS are colored by ECOSYSTEM TYPE (same color as boxplot face)"""
    
    print(f"\nCreating SINGLE PANEL FIGURE 2:")
    print(f"  • BOXES colored by SALINITY:")
    for salinity, color in SALINITY_BOX_COLORS.items():
        print(f"    - {salinity} boxes: {color}")
    print(f"  • POINTS colored by ECOSYSTEM TYPE (matching box colors)")
    print(f"  • NOT using site-specific colors from Figure 1")
    print(f"  • Strict triple intersection (ALL metrics): {len(shared_sites)} sites")
    print(f"  • Capping for visualization: {CAPPING_PERCENTILE}th percentile (same as Figure 1)")
    print(f"  • Y-axis: AUTO-SCALED (no constraints)")
    print(f"  • Statistical tests use ORIGINAL uncapped values")
    print(f"  Figure size: {FIGURE_WIDTH}x{FIGURE_HEIGHT} inches")
    
    # Create SINGLE figure
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    
    # Metrics to plot
    metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    salinities = ['Upland', 'Freshwater', 'Saline']
    
    # Store cap values for reporting
    cap_values_dict = {}
    
    # Calculate positions for boxplots (3 categories per metric)
    positions = []
    box_data = []
    box_labels = []
    
    for i, metric in enumerate(metrics):
        base_pos = i * (len(salinities) + METRIC_SPACING)
        
        for j, salinity in enumerate(salinities):
            pos = base_pos + j * WITHIN_METRIC_SPACING
            positions.append(pos)
            
            # Get data for this box - ALL metrics use strict triple intersection
            subset = data[(data['WUE_Metric'] == metric) & 
                          (data['Salinity_Category'] == salinity)]
            
            original_values = subset['Median'].dropna().values
            
            if len(original_values) == 0:
                box_data.append([])
                box_labels.append(f"{metric}_{salinity}")
                continue
            
            # Apply capping for VISUALIZATION ONLY (same as Figure 1)
            if USE_CAPPING and metric in ['WUE_eva', 'WUE_tra']:
                values_capped, cap_value, n_capped, extreme_values = cap_values_for_visualization(original_values, CAPPING_PERCENTILE)
                box_data.append(values_capped)
                cap_key = f"{metric}_{salinity}"
                cap_values_dict[cap_key] = {
                    'cap_value': cap_value,
                    'n_capped': n_capped,
                    'total_n': len(original_values),
                    'extreme_values': extreme_values
                }
                if n_capped > 0:
                    print(f"  {metric} - {salinity}: Capped {n_capped}/{len(original_values)} values at {cap_value:.3f}")
            else:
                # WUE remains UNCAPPED (same as Figure 1)
                box_data.append(original_values)
                cap_key = f"{metric}_{salinity}"
                cap_values_dict[cap_key] = None
    
    # Create boxplot
    boxplot = ax.boxplot(
        box_data,
        positions=positions,
        patch_artist=True,
        widths=BOX_WIDTH,
        showfliers=False
    )
    
    # Color BOXES by SALINITY
    for i, (patch, pos) in enumerate(zip(boxplot['boxes'], positions)):
        salinity_idx = i % len(salinities)
        patch.set_facecolor(SALINITY_BOX_COLORS[salinities[salinity_idx]])
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
        patch.set_linewidth(BOX_BORDER_WIDTH)
    
    # Customize whiskers and caps
    for element in ['whiskers', 'caps']:
        for line in boxplot[element]:
            line.set_color('black')
            line.set_linewidth(WHISKER_WIDTH)
    
    # Customize median lines (red)
    for median in boxplot['medians']:
        median.set_color('red')
        median.set_linewidth(MEDIAN_WIDTH)
    
    # POINTS: Color by ECOSYSTEM TYPE (matching box colors)
    for i, metric in enumerate(metrics):
        base_pos = i * (len(salinities) + METRIC_SPACING)
        
        for j, salinity in enumerate(salinities):
            pos = base_pos + j * WITHIN_METRIC_SPACING
            subset = data[(data['WUE_Metric'] == metric) & 
                          (data['Salinity_Category'] == salinity)]
            
            # Get the color for this ecosystem type
            point_color = SALINITY_BOX_COLORS[salinity]
            
            for idx, row in subset.iterrows():
                site = row['site_name']
                original_value = row['Median']
                
                # Apply capping for VISUALIZATION ONLY
                if USE_CAPPING and metric in ['WUE_eva', 'WUE_tra']:
                    cap_key = f"{metric}_{salinity}"
                    if cap_key in cap_values_dict and cap_values_dict[cap_key] is not None:
                        plot_value = min(original_value, cap_values_dict[cap_key]['cap_value'])
                    else:
                        plot_value = original_value
                else:
                    plot_value = original_value
                
                # Add jitter
                x_jitter = np.random.normal(pos, 0.08)
                
                # Color points by ECOSYSTEM TYPE (matching box color)
                ax.scatter(
                    x_jitter, plot_value,
                    s=POINT_SIZE,
                    color=point_color,  # Same color as boxplot for this ecosystem
                    edgecolor='black',
                    linewidth=0.5,
                    alpha=0.6,
                    zorder=3
                )
    
    # Set y-axis to AUTO-SCALE (no constraints)
    if Y_AXIS_MIN is not None:
        ax.set_ylim(bottom=Y_AXIS_MIN)
    if Y_AXIS_MAX is not None:
        ax.set_ylim(top=Y_AXIS_MAX)
    
    # If both are None, matplotlib auto-scales
    if Y_AXIS_MIN is None and Y_AXIS_MAX is None:
        print(f"\n  Y-axis: AUTO-SCALED (no constraints applied)")
    elif Y_AXIS_MIN is not None and Y_AXIS_MAX is not None:
        print(f"\n  Y-axis: Manually set to [{Y_AXIS_MIN}, {Y_AXIS_MAX}]")
    elif Y_AXIS_MIN is not None:
        print(f"\n  Y-axis: Min = {Y_AXIS_MIN}, Max = auto-scaled")
    else:
        print(f"\n  Y-axis: Min = auto-scaled, Max = {Y_AXIS_MAX}")
    
    # Set x-axis ticks and labels
    x_tick_positions = []
    for i in range(len(metrics)):
        base_pos = i * (len(salinities) + METRIC_SPACING)
        center_pos = base_pos + (len(salinities) - 1) * WITHIN_METRIC_SPACING / 2
        x_tick_positions.append(center_pos)
    
    ax.set_ylim(bottom=0, top=13.5)
    ax.set_yticks([0, 5, 10, 13])
    ax.set_xticks(x_tick_positions)
    ax.set_xticklabels(['WUE$_{ET}$', 'WUE$_E$', 'WUE$_T$'], 
                       fontsize=X_TICK_LABEL_SIZE)
    
    # Set y-axis label
    ax.set_ylabel('Median WUE (g C kg$^{-1}$ H$_2$O$^{-1}$)',
                  fontsize=Y_AXIS_LABEL_SIZE, 
                  fontweight='bold', 
                  labelpad=20)
    
    # Set tick parameters with VISIBLE TICKS
    ax.tick_params(axis='x', which='major', 
                   labelsize=X_TICK_LABEL_SIZE,
                   length=10,
                   width=3,
                   pad=15)
    
    ax.tick_params(axis='y', which='major',
                   labelsize=Y_TICK_LABEL_SIZE,
                   length=10,
                   width=3,
                   pad=15)
    
    # Add vertical separators between metrics
    separator_positions = []
    for i in range(len(metrics) - 1):
        base_pos_current = i * (len(salinities) + METRIC_SPACING)
        base_pos_next = (i + 1) * (len(salinities) + METRIC_SPACING)
        sep_pos = (base_pos_current + len(salinities) * WITHIN_METRIC_SPACING + base_pos_next) / 2
        separator_positions.append(sep_pos)
    
    for sep_pos in separator_positions:
        ax.axvline(x=sep_pos, color='gray', linestyle='--', alpha=0.3, linewidth=2)
    
    # Remove grid
    ax.grid(False)
    
    # Add thick border
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')
    
    # LEGEND: Show SALINITY BOX colors (3 categories)
    legend_elements = []
    for salinity in salinities:
        if salinity in SALINITY_BOX_COLORS:
            legend_elements.append(
                Patch(facecolor=SALINITY_BOX_COLORS[salinity], alpha=0.7, 
                      edgecolor='black', linewidth=BOX_BORDER_WIDTH, label=salinity)
            )
    
    # Position legend
    ax.legend(handles=legend_elements, loc='upper left', 
             fontsize=LEGEND_FONT_SIZE, 
             framealpha=0.95, 
             fancybox=True,
             borderpad=1,
             labelspacing=0.5,
             handletextpad=0.5,
             handlelength=1.0)
    
    # Adjust layout
    plt.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    
    # Add panel label (b) in top-left corner
    ax.text(0.02, 0.98, '(b)', transform=ax.transAxes, 
            fontsize=32, fontweight='bold', zorder=10,
            verticalalignment='top', horizontalalignment='left')
    
    return fig, cap_values_dict

def save_figure(fig, filename="Figure2_WUE_Salinity_ThreeCategories"):
    """Save figure with current settings"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    if OUTPUT_FORMAT in ['png', 'both']:
        png_path = os.path.join(OUTPUT_DIR, f"{filename}.png")
        fig.savefig(png_path, dpi=DPI, bbox_inches='tight', facecolor='white')
        print(f"✅ Saved PNG: {png_path}")
    
    if OUTPUT_FORMAT in ['pdf', 'both']:
        pdf_path = os.path.join(OUTPUT_DIR, f"{filename}.pdf")
        fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
        print(f"✅ Saved PDF: {pdf_path}")
    
    plt.show()

# =============================================================================
# PRINT SUMMARY STATISTICS - ALL METRICS USE STRICT TRIPLE INTERSECTION
# =============================================================================

def print_summary_statistics(data, shared_sites):
    """Print comprehensive summary statistics for manuscript using ORIGINAL values
       ALL metrics now use the same strict triple intersection"""
    
    print("\n" + "="*80)
    print("FIGURE 2 SUMMARY STATISTICS (ORIGINAL UNCAPPED VALUES)")
    print("NOTE: Statistics use original values - capping applied for visualization only")
    print("      ALL metrics use the same strict triple intersection")
    print("="*80)
    
    metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    salinities = ['Freshwater', 'Saline', 'Upland']
    
    for metric in metrics:
        print(f"\n{'='*60}")
        print(f"METRIC: {metric}")
        print(f"{'='*60}")
        
        metric_data = data[data['WUE_Metric'] == metric].copy()
        
        print(f"  Data source: wue_site_level_summary_SPEI_1.csv (filtered by CHUNK 3)")
        used_sites = metric_data['site_name'].unique()
        print(f"  Total sites used: {len(used_sites)}")
        
        print(f"\n  Sample sizes by salinity category (ORIGINAL uncapped values):")
        
        for salinity in salinities:
            subset = metric_data[metric_data['Salinity_Category'] == salinity]
            values = subset['Median'].dropna()
            n_sites = len(values)
            
            if n_sites > 0:
                # Calculate statistics
                median_val = np.median(values)
                mean_val = np.mean(values)
                sd_val = np.std(values, ddof=1)  # sample SD
                se_val = sd_val / np.sqrt(n_sites)  # standard error
                iqr_25 = np.percentile(values, 25)
                iqr_75 = np.percentile(values, 75)
                min_val = np.min(values)
                max_val = np.max(values)
                
                print(f"\n  {salinity}:")
                print(f"    n = {n_sites} sites")
                print(f"    Median = {median_val:.3f}")
                print(f"    Mean ± SE = {mean_val:.3f} ± {se_val:.3f}")
                print(f"    SD = {sd_val:.3f}")
                print(f"    IQR = {iqr_25:.3f} - {iqr_75:.3f}")
                print(f"    Range = {min_val:.3f} - {max_val:.3f}")
            else:
                print(f"\n  {salinity}: No data available")

# =============================================================================
# PRINT CAPPING SUMMARY
# =============================================================================

def print_capping_summary(cap_values_dict):
    """Print summary of values capped for visualization"""
    
    print("\n" + "="*80)
    print("VISUALIZATION CAPPING SUMMARY (Same as Figure 1)")
    print("="*80)
    print(f"Capping percentile: {CAPPING_PERCENTILE}th")
    print("Note: Capping applied ONLY for visualization - statistics use original values\n")
    
    for key, cap_info in cap_values_dict.items():
        if cap_info is not None and cap_info['n_capped'] > 0:
            # Handle key format: "WUE_eva_Freshwater" or "WUE_tra_Upland"
            parts = key.split('_')
            if len(parts) >= 3:
                metric = f"{parts[0]}_{parts[1]}"
                salinity = '_'.join(parts[2:])
            else:
                metric = parts[0]
                salinity = parts[1] if len(parts) > 1 else "Unknown"
            
            if metric == 'WUE_eva':
                metric_display = 'WUE_E'
            elif metric == 'WUE_tra':
                metric_display = 'WUE_T'
            else:
                metric_display = metric
            
            print(f"{metric_display} - {salinity}:")
            print(f"  Cap value: {cap_info['cap_value']:.3f} g C kg⁻¹ H₂O")
            print(f"  Capped {cap_info['n_capped']}/{cap_info['total_n']} ({cap_info['n_capped']/cap_info['total_n']*100:.1f}%) values")
            if len(cap_info['extreme_values']) > 0:
                print(f"  Extreme values range: {np.min(cap_info['extreme_values']):.3f} - {np.max(cap_info['extreme_values']):.3f}")

# =============================================================================
# MAIN WORKFLOW FOR FIGURE 2
# =============================================================================

def main():
    """Main workflow - creates SINGLE PANEL Figure 2 with capping for visualization (NO Y-LIMIT)
       ALL metrics now use the same strict triple intersection
       Points are colored by ECOSYSTEM TYPE (matching box colors)"""
    
    print("="*80)
    print("FIGURE 2: WUE BY SALINITY - SINGLE PANEL (UPDATED)")
    print("="*80)
    print("BOXES: Colored by SALINITY (Blue=Freshwater, Orange=Saline, Purple=Upland)")
    print("POINTS: Colored by ECOSYSTEM TYPE (matching box colors)")
    print("STRICT TRIPLE INTERSECTION: ALL metrics now use the same shared subset")
    print("NOTE: Filtering already applied by CHUNK 3 (each metric ≥3 months NN)")
    print(f"CAPPING: {CAPPING_PERCENTILE}th percentile for visualization (same as Figure 1)")
    print("Y-AXIS: AUTO-SCALED (no constraints applied)")
    print("SINGLE PANEL: All WUE metrics in one figure")
    print("="*80)
    
    # Step 1: Load MASTER color mapping from Figure 1 (for reference only)
    print("\n[Step 1] Loading MASTER color mapping from Figure 1 (for reference only)...")
    master_color_map = load_master_color_mapping()
    
    # Step 2: Load and prepare data (filtering already done by CHUNK 3)
    print("\n[Step 2] Loading and preparing Figure 2 data...")
    data, shared_sites = load_and_prepare_data(master_color_map)
    
    # Step 3: Print summary statistics (using ORIGINAL uncapped values)
    print("\n[Step 3] Calculating statistics (using ORIGINAL values)...")
    print_summary_statistics(data, shared_sites)
    
    # Step 4: Perform statistical tests (using ORIGINAL uncapped values)
    print("\n[Step 4] Performing statistical tests (using ORIGINAL values)...")
    stat_results = perform_statistical_tests(data, shared_sites)
    
    # Step 5: Create SINGLE PANEL Figure 2 with capping for visualization (NO Y-LIMIT)
    print("\n[Step 5] Creating SINGLE PANEL Figure 2 with capping (auto-scaled y-axis)...")
    fig, cap_values_dict = create_figure_formatting(data, master_color_map, shared_sites)
    
    # Step 6: Print capping summary
    print_capping_summary(cap_values_dict)
    
    # Step 7: Save figure
    print("\n[Step 6] Saving Figure 2...")
    save_figure(fig, "Figure2_WUE_By_Salinity_ThreeCategories_Capped_AutoScale_EcosystemPoints")
    
    print("\n" + "="*80)
    print("FIGURE 2 WORKFLOW COMPLETE")
    print("="*80)
    print(f"\nFigure Specifications:")
    print(f"  Dimensions: {FIGURE_WIDTH} × {FIGURE_HEIGHT} inches")
    print(f"  Box colors: Blue=Freshwater, Orange=Saline, Purple=Upland")
    print(f"  Point colors: SAME as box colors (by ecosystem type)")
    print(f"  Capping: {CAPPING_PERCENTILE}th percentile for visualization only")
    print(f"  Y-axis: AUTO-SCALED (no constraints)")
    print(f"  DPI: {DPI}")
    print(f"  Output format: {OUTPUT_FORMAT}")
    print(f"\nData Specifications:")
    print(f"  Filtering applied by CHUNK 3 (each metric ≥3 months NN)")
    print(f"  ALL metrics (WUE, WUE_E, WUE_T): Strict triple intersection ({len(shared_sites)} sites)")
    print(f"  ✓ WUE, WUE_E, and WUE_T all have IDENTICAL sample sizes")
    print(f"  WUE_E and WUE_T: Capped at {CAPPING_PERCENTILE}th percentile for visualization only")
    print(f"  WUE: UNCAPPED in visualization")
    print(f"  Salinity categories: Freshwater, Saline, Upland")
    print(f"\nStatistical tests: ALL use ORIGINAL uncapped values")
    print(f"\nSummary statistics include: Mean ± SE and SD")
    
    return data, shared_sites, master_color_map, fig

# =============================================================================
# EXECUTE FIGURE 2
# =============================================================================

if __name__ == "__main__":
    # Run Figure 2 workflow
    data, shared_sites, master_colors, figure = main()