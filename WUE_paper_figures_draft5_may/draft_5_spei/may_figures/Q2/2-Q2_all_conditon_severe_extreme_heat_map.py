# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 11:09:49 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
COMBINED FIGURE: Two-panel figure with Panel A (All conditions) and Panel B (Severe/Extreme conditions)
- Shared calculations using strict double intersection (sites with both metrics)
- Single colorbar for both panels
- Panel (a) on left: All conditions (SPEI-6/48 × Dry/Wet)
- Panel (b) on right: Severe/Extreme conditions (SPEI-6/48 × Severe Dry/Severe Wet)
- Both panels show median % change for WUE_ET and WUE_T
- Includes N values in y-axis labels (larger font)
- Fixed: Proper spacing between panels, colorbar below panels, no overlap
- UPDATED: Uses NN bootstrap CI site-level classifications from SPEI_site_level_details_FINAL.csv
- FIXED: Colormap properly centered at 0 (negative=gray, positive=green)
- FIXED: ALL positive values show green, ALL negative values show gray
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Patch
from scipy import stats
from scipy.stats import wilcoxon, ttest_1samp
import os

# =============================================================================
# GLOBAL SETTINGS - LARGER FONTS FOR BETTER LEGIBILITY
# =============================================================================
plt.rcParams.update({
    'font.size': 20,           # Base font size
    'axes.titlesize': 22,      
    'axes.labelsize': 20,      
    'xtick.labelsize': 18,     
    'ytick.labelsize': 18,     
    'legend.fontsize': 18,     
})

# =============================================================================
# NOTE:
# The old bootstrap_median_ci() function is no longer used for direction classification.
# Direction now comes from the updated site-level NN bootstrap CI classification
# already stored in SPEI_site_level_details_FINAL.csv.
# =============================================================================

# =============================================================================
# LOAD AND RECOMPUTE STATISTICS FOR BOTH WORKFLOWS (SHARED)
# UPDATED: Uses site-level NN bootstrap CI classifications
# =============================================================================
def load_and_recompute_all_stats(site_level_path):
    """
    Load site-level data and recompute statistics for both:
    - ALL conditions (Dry, Wet) - for Panel A
    - Severe/Extreme conditions (Severe Dry, Severe Wet) - for Panel B
    Both use strict double intersection (sites with both WUE and WUE_tra)
    
    UPDATED: Direction comes from site-level NN bootstrap CI classifications,
    not from recomputed bootstrap around group medians.
    """
    print("=" * 100)
    print("LOADING DATA AND RECOMPUTING STATISTICS FOR BOTH PANELS")
    print("Using strict double intersection (sites with both WUE_ET and WUE_T)")
    print("Direction classification uses site-level NN bootstrap CI results")
    print("=" * 100)
    
    df = pd.read_csv(site_level_path)
    print(f"\n[Loaded] {len(df)} rows from site-level file")
    
    # Define conditions for each panel
    panel_a_conditions = ['Dry (all)', 'Wet (all)']
    panel_b_conditions = ['Severe Dry', 'Severe Wet']
    timescales = ['SPEI_6', 'SPEI_48']
    metrics = ['WUE', 'WUE_tra']  # Only WUE_ET and WUE_T
    
    results_panel_a = []
    results_panel_b = []
    
    # Process each timescale and condition
    for timescale in timescales:
        for condition in panel_a_conditions + panel_b_conditions:
            print(f"\n{'='*50}")
            print(f"Processing: {timescale} - {condition}")
            
            # CRITICAL FIX: Filter to ONLY 'All' salinity to avoid duplicate rows
            condition_data = df[
                (df['Salinity'] == 'All') &  # <-- THIS IS THE FIX
                (df['SPEI_Timescale'] == timescale) & 
                (df['Condition'] == condition) &
                (df['WUE_Metric'].isin(metrics))
            ].copy()
            
            if len(condition_data) == 0:
                print(f"  No data found")
                continue
            
            # Get sites for each metric
            wue_sites = set(condition_data[condition_data['WUE_Metric'] == 'WUE']['Site'].unique())
            tra_sites = set(condition_data[condition_data['WUE_Metric'] == 'WUE_tra']['Site'].unique())
            
            print(f"  Sites with WUE_ET: {len(wue_sites)}")
            print(f"  Sites with WUE_T: {len(tra_sites)}")
            
            # STRICT DOUBLE INTERSECTION
            shared_sites = wue_sites.intersection(tra_sites)
            print(f"  Shared sites (both metrics): {len(shared_sites)}")
            
            # Process each metric
            for metric in metrics:
                metric_data = condition_data[
                    (condition_data['WUE_Metric'] == metric) & 
                    (condition_data['Site'].isin(shared_sites))
                ].copy()
                
                if len(metric_data) == 0:
                    continue
                
                # Aggregate to site level
                site_medians = metric_data.groupby('Site')['Median_%_Change'].median()
                pct_changes = site_medians.values
                n_sites = len(pct_changes)
                
                # Calculate heatmap value (median of site medians)
                median_val = np.median(pct_changes) if n_sites > 0 else np.nan
                
                # Calculate mean, SD, SE (keep for reporting)
                if n_sites >= 2:
                    mean_val = np.mean(pct_changes)
                    sd_val = np.std(pct_changes, ddof=1)
                    se_val = sd_val / np.sqrt(n_sites)
                else:
                    mean_val, sd_val, se_val = np.nan, np.nan, np.nan
                
                # ========== UPDATED: Use site-level NN bootstrap CI classifications ==========
                # Count site-level directions from the updated columns
                # Now these counts will be correct because we filtered Salinity == 'All'
                direction_counts = metric_data['Direction'].value_counts()
                increase_sites = direction_counts.get('Increase', 0)
                decrease_sites = direction_counts.get('Decrease', 0)
                no_change_sites = direction_counts.get('No change', 0)
                insufficient_data_sites = direction_counts.get('Insufficient data', 0)
                total_valid = increase_sites + decrease_sites
                
                # Determine group-level direction based on site-level counts
                if increase_sites > decrease_sites and increase_sites > 0:
                    direction = "More sites increased beyond NN bootstrap CI"
                elif decrease_sites > increase_sites and decrease_sites > 0:
                    direction = "More sites decreased beyond NN bootstrap CI"
                elif increase_sites == 0 and decrease_sites == 0:
                    direction = "No sites exceeded NN bootstrap CI"
                else:
                    direction = "Mixed site-level responses"
                
                # Get Direction_Method and Threshold_Method from first non-null row
                direction_method = metric_data['Direction_Method'].iloc[0] if 'Direction_Method' in metric_data.columns else 'NN_Bootstrap_Median_CI_95'
                threshold_method = metric_data['Threshold_Method'].iloc[0] if 'Threshold_Method' in metric_data.columns else 'Site_specific_NN_bootstrap_95CI_median'
                
                # Calculate median of site-level CI bounds (for reporting only)
                if 'CI_Lower' in metric_data.columns:
                    ci_lower = metric_data['CI_Lower'].median()
                else:
                    ci_lower = np.nan
                    
                if 'CI_Upper' in metric_data.columns:
                    ci_upper = metric_data['CI_Upper'].median()
                else:
                    ci_upper = np.nan
                    
                if 'CI_Width' in metric_data.columns:
                    ci_width = metric_data['CI_Width'].median()
                else:
                    ci_width = np.nan
                
                # Set t_statistic and p_value to NaN (one-sided t-test removed)
                t_stat = np.nan
                p_val = np.nan
                
                result = {
                    'SPEI_Timescale': timescale,
                    'Condition': condition,
                    'WUE_Metric': metric,
                    'n_sites': n_sites,
                    'Median_%_Change': median_val,
                    'CI_Lower': ci_lower,
                    'CI_Upper': ci_upper,
                    'CI_Width': ci_width,
                    'Mean_%_Change': mean_val,
                    'SD_%_Change': sd_val,
                    'SE_%_Change': se_val,
                    'Direction': direction,
                    'Direction_Method': direction_method,
                    'Threshold_Method': threshold_method,
                    'Increase_sites': increase_sites,
                    'Decrease_sites': decrease_sites,
                    'No_change_sites': no_change_sites,
                    'Insufficient_data_sites': insufficient_data_sites,
                    't_statistic': t_stat,
                    'p_value': p_val
                }
                
                if condition in panel_a_conditions:
                    results_panel_a.append(result)
                else:
                    results_panel_b.append(result)
                
                # Print summary
                metric_label = 'WUE_ET' if metric == 'WUE' else 'WUE_T'
                print(f"\n  {metric_label} (n={n_sites}):")
                print(f"    Median = {median_val:+.1f}%" if not np.isnan(median_val) else "    Median = NA")
                if not np.isnan(median_val):
                    print(f"    Mean ± SE = {mean_val:+.1f}% ± {se_val:.1f}%" if not np.isnan(se_val) else f"    Mean = {mean_val:+.1f}%")
                    print(f"    Direction: {direction}")
                    print(f"    Site-level counts: Increase={increase_sites}, Decrease={decrease_sites}, No change={no_change_sites}, Insufficient={insufficient_data_sites}")
                    # Verify counts sum correctly
                    total_recorded = increase_sites + decrease_sites + no_change_sites + insufficient_data_sites
                    if total_recorded != n_sites:
                        print(f"    WARNING: Count mismatch! N_sites={n_sites}, Total recorded={total_recorded}")
    
    df_panel_a = pd.DataFrame(results_panel_a)
    df_panel_b = pd.DataFrame(results_panel_b)
    
    print(f"\n{'='*60}")
    print(f"RECOMPUTED SUMMARY")
    print(f"Panel A (All conditions): {len(df_panel_a)} records")
    print(f"Panel B (Severe/Extreme): {len(df_panel_b)} records")
    
    return df_panel_a, df_panel_b

# =============================================================================
# CREATE CUSTOM COLORMAP - FIXED: POSITIVE = GREEN, NEGATIVE = GRAY
# =============================================================================
def create_diverging_colormap(vmin, vmax, vcenter=0, n_colors=256):
    """
    Create a diverging colormap where:
    - Values below vcenter (negative) use gray shades (light to dark gray)
    - Values above vcenter (positive) use green shades (light to dark green)
    - ALL positive values are green, ALL negative values are gray
    """
    # Calculate proportion for negative and positive sides
    if vmax > vcenter and vmin < vcenter:
        # Negative side proportion
        neg_ratio = (vcenter - vmin) / (vmax - vmin)
        # Positive side proportion
        pos_ratio = (vmax - vcenter) / (vmax - vmin)
        
        n_neg = int(n_colors * neg_ratio)
        n_pos = int(n_colors * pos_ratio)
        
        # Ensure at least some colors on each side
        if n_neg < 2:
            n_neg = 2
        if n_pos < 2:
            n_pos = 2
        
        # Gray colors for negative side (from light gray to medium gray)
        # Using reversed grays so more negative = darker gray
        gray_colors = plt.cm.Greys_r(np.linspace(0.2, 0.7, n_neg))
        
        # Green colors for positive side (from light green to dark green)
        green_colors = plt.cm.Greens(np.linspace(0.4, 0.9, n_pos))
        
        # Combine colors
        colors = np.vstack((gray_colors, green_colors))
    else:
        # Fallback if no negative or no positive values
        if vmax > vcenter:
            # Only positive values
            colors = plt.cm.Greens(np.linspace(0.4, 0.9, n_colors))
        else:
            # Only negative values
            colors = plt.cm.Greys_r(np.linspace(0.2, 0.7, n_colors))
    
    return ListedColormap(colors)

# =============================================================================
# CREATE TWO-PANEL FIGURE - FIXED VERSION
# =============================================================================
def create_two_panel_figure(df_panel_a, df_panel_b, save_path):
    """
    Create side-by-side heatmaps for Panel A and Panel B
    with shared colorbar, larger labels, and clear panel distinctions
    FIXED: Proper spacing, colorbar below panels, no overlap
    FIXED: Colormap properly centered at 0 (negative=gray, positive=green)
    """
    
    print("\n" + "=" * 100)
    print("CREATING TWO-PANEL FIGURE")
    print("Panel A (left): All conditions")
    print("Panel B (right): Severe/Extreme conditions")
    print("=" * 100)
    
    # Define row orders for each panel
    rows_panel_a = [
        ('SPEI_6', 'Dry (all)'),
        ('SPEI_48', 'Dry (all)'),
        ('SPEI_6', 'Wet (all)'),
        ('SPEI_48', 'Wet (all)')
    ]
    
    rows_panel_b = [
        ('SPEI_6', 'Severe Dry'),
        ('SPEI_48', 'Severe Dry'),
        ('SPEI_6', 'Severe Wet'),
        ('SPEI_48', 'Severe Wet')
    ]
    
    metric_order = ['WUE', 'WUE_tra']
    metric_labels_display = ['WUE$_{ET}$', 'WUE$_T$']
    
    # Extract data matrices
    def extract_data_matrix(df, rows):
        """Extract median values into matrix"""
        matrix = []
        n_sites_list = []
        row_labels = []
        
        for timescale, condition in rows:
            row_vals = []
            n_sites_val = None
            
            for metric in metric_order:
                row_data = df[
                    (df['SPEI_Timescale'] == timescale) &
                    (df['Condition'] == condition) &
                    (df['WUE_Metric'] == metric)
                ]
                
                if len(row_data) > 0:
                    val = row_data.iloc[0]['Median_%_Change']
                    if n_sites_val is None:
                        n_sites_val = int(row_data.iloc[0]['n_sites'])
                else:
                    val = np.nan
                
                row_vals.append(val)
            
            matrix.append(row_vals)
            n_sites_list.append(n_sites_val if n_sites_val else 0)
            
            # Create row label
            ts_display = timescale.replace('_', '-')
            cond_display = condition.replace(' (all)', '').replace('Severe ', 'Ext. ')
            row_labels.append(f"{ts_display}\n{cond_display}")
        
        return np.array(matrix), row_labels, n_sites_list
    
    data_a, labels_a, n_a = extract_data_matrix(df_panel_a, rows_panel_a)
    data_b, labels_b, n_b = extract_data_matrix(df_panel_b, rows_panel_b)
    
    # Determine common color scale for both panels - MANUALLY FIXED
    vmax = 8
    vmin = -15
    vcenter = 0
    
    # Print verification
    all_valid = np.concatenate([data_a[~np.isnan(data_a)], data_b[~np.isnan(data_b)]])
    print(f"Data range: min={np.nanmin(all_valid):.1f}, max={np.nanmax(all_valid):.1f}")
    print(f"Color scale: vmin={vmin:.1f}, vcenter={vcenter:.1f}, vmax={vmax:.1f}")
    
    # Create custom diverging colormap (negative=gray, positive=green)
    custom_cmap = create_diverging_colormap(vmin, vmax, vcenter)
    
    # Create figure with GridSpec for better control
    fig = plt.figure(figsize=(20, 11))
    
    # Use GridSpec: [top_margin, main_plot_area, bottom_margin_for_colorbar, bottom_margin_for_legend]
    gs = fig.add_gridspec(3, 2, 
                          height_ratios=[0.08, 1, 0.06],
                          width_ratios=[0.48, 0.48],
                          hspace=0.2,
                          wspace=0.5)
    
    # Create axes for the two panels
    ax_left = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1])
    
    # Use custom normalization to ensure correct mapping
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    # ========== PANEL A (LEFT) ==========
    im_left = ax_left.imshow(data_a, cmap=custom_cmap, aspect='auto', norm=norm)
    
    # X-axis labels
    ax_left.set_xticks(range(len(metric_labels_display)))
    ax_left.set_xticklabels(metric_labels_display, fontweight='bold', fontsize=22)
    
    # Y-axis labels with N values (larger font)
    y_labels_a = []
    for i, (label, n) in enumerate(zip(labels_a, n_a)):
        lines = label.split('\n')
        y_labels_a.append(f"{lines[0]}\n{lines[1]} (N={n})")
    
    ax_left.set_yticks(range(len(y_labels_a)))
    ax_left.set_yticklabels(y_labels_a, fontweight='bold', fontsize=20)
    
    # Add value annotations
    for i in range(data_a.shape[0]):
        for j in range(data_a.shape[1]):
            val = data_a[i, j]
            if not np.isnan(val):
                # Determine text color based on value
                if val > 0:
                    # Positive values - green background, use white text for dark green
                    if val > 4:
                        text_color = 'white'
                    else:
                        text_color = 'black'
                else:
                    # Negative values - gray background
                    text_color = 'black'
                val_str = f"{val:+.1f}%" if val != 0 else "0.0%"
                ax_left.text(j, i, val_str, ha='center', va='center', 
                            fontsize=20, fontweight='bold', color=text_color)
    
    # Panel label (a) - LARGE and bold, positioned above the plot area
    ax_left.text(-0.25, 1.08, '(a)', transform=ax_left.transAxes,
                fontsize=36, fontweight='bold', va='bottom', ha='left')
    
    # Add title for Panel A
    ax_left.set_title('All conditions (Dry/Wet)', fontsize=22, fontweight='bold', pad=20)
    
    # Add separation line between Dry and Wet
    ax_left.axhline(y=1.5, color='black', linewidth=2.5, linestyle='-', alpha=0.8)
    
    # Invert y-axis
    ax_left.invert_yaxis()
    
    # Remove top and right spines
    ax_left.spines['top'].set_visible(False)
    ax_left.spines['right'].set_visible(False)
    ax_left.spines['left'].set_linewidth(1.5)
    ax_left.spines['bottom'].set_linewidth(1.5)
    
    # ========== PANEL B (RIGHT) ==========
    ax_right.imshow(data_b, cmap=custom_cmap, aspect='auto', norm=norm)
    
    # X-axis labels
    ax_right.set_xticks(range(len(metric_labels_display)))
    ax_right.set_xticklabels(metric_labels_display, fontweight='bold', fontsize=22)
    
    # Y-axis labels with N values
    y_labels_b = []
    for i, (label, n) in enumerate(zip(labels_b, n_b)):
        lines = label.split('\n')
        y_labels_b.append(f"{lines[0]}\n{lines[1]} (N={n})")
    
    ax_right.set_yticks(range(len(y_labels_b)))
    ax_right.set_yticklabels(y_labels_b, fontweight='bold', fontsize=20)
    
    # Add value annotations
    for i in range(data_b.shape[0]):
        for j in range(data_b.shape[1]):
            val = data_b[i, j]
            if not np.isnan(val):
                # Determine text color based on value
                if val > 0:
                    if val > 4:
                        text_color = 'white'
                    else:
                        text_color = 'black'
                else:
                    text_color = 'black'
                val_str = f"{val:+.1f}%" if val != 0 else "0.0%"
                ax_right.text(j, i, val_str, ha='center', va='center', 
                             fontsize=20, fontweight='bold', color=text_color)
    
    # Panel label (b)
    ax_right.text(-0.25, 1.08, '(b)', transform=ax_right.transAxes,
                 fontsize=36, fontweight='bold', va='bottom', ha='left')
    
    # Add title for Panel B
    ax_right.set_title('Severe/Extreme conditions', fontsize=22, fontweight='bold', pad=20)
    
    # Add separation line between Dry and Wet
    ax_right.axhline(y=1.5, color='black', linewidth=2.5, linestyle='-', alpha=0.8)
    
    # Invert y-axis
    ax_right.invert_yaxis()
    
    # Remove spines - keep left spine for consistency but hide ticks
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(True)
    ax_right.spines['left'].set_linewidth(1.5)
    ax_right.spines['bottom'].set_linewidth(1.5)
    
    # ========== SHARED COLORBAR (Below panels, centered) ==========
    cbar_ax = fig.add_subplot(gs[2, :])
    cbar = fig.colorbar(im_left, cax=cbar_ax, orientation='horizontal', fraction=0.7)
    cbar.set_label('Median % Change', fontweight='bold', fontsize=20, labelpad=10)
    
    # MANUALLY SET COLORBAR TICKS (based on vmin=-15, vmax=8)
    tick_values = [-15, -10, -5, 0, 5, 8]
    cbar.set_ticks(tick_values)
    cbar.set_ticklabels([str(tick) for tick in tick_values])
    cbar.ax.tick_params(labelsize=20)
    
    # ========== SHARED LEGEND (Below colorbar) ==========
    legend_elements = [
        Patch(facecolor='#27ae60', edgecolor='black', label='Increase (+ve)'),
        Patch(facecolor='#7b7d7d', edgecolor='black', label='Decrease (-ve)'),
    ]
    
    # Add legend below the figure
    fig.legend(handles=legend_elements, loc='lower center', 
               bbox_to_anchor=(0.5, -0.07), ncol=2, fontsize=22,
               frameon=True, framealpha=0.95, borderpad=1)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure with high DPI
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Two-panel figure saved to: {save_path}")
    
    plt.show()
    return fig

# =============================================================================
# PRINT RESULTS SUMMARY
# =============================================================================
def print_results_summary(df_panel_a, df_panel_b):
    """Print formatted results summary for both panels"""
    
    print("\n" + "=" * 100)
    print("RESULTS SUMMARY")
    print("=" * 100)
    print("NOTE: Direction based on site-level NN bootstrap CI classifications")
    print("=" * 100)
    
    print("\n" + "-" * 50)
    print("PANEL A: ALL CONDITIONS (Dry/Wet)")
    print("-" * 50)
    
    for timescale in ['SPEI_6', 'SPEI_48']:
        ts_display = timescale.replace('_', '-')
        print(f"\n{ts_display}:")
        
        for condition in ['Dry (all)', 'Wet (all)']:
            cond_display = condition.replace(' (all)', '')
            print(f"\n  {cond_display}:")
            
            for metric in ['WUE', 'WUE_tra']:
                row = df_panel_a[
                    (df_panel_a['SPEI_Timescale'] == timescale) &
                    (df_panel_a['Condition'] == condition) &
                    (df_panel_a['WUE_Metric'] == metric)
                ]
                
                if len(row) > 0:
                    row = row.iloc[0]
                    metric_label = 'WUE_ET' if metric == 'WUE' else 'WUE_T'
                    print(f"    {metric_label} (N={int(row['n_sites'])}): Median = {row['Median_%_Change']:+.1f}%, Mean ± SE = {row['Mean_%_Change']:+.1f}% ± {row['SE_%_Change']:.1f}%")
                    print(f"      Direction: {row['Direction']}")
                    print(f"      Site counts: ↑{row['Increase_sites']} ↓{row['Decrease_sites']} ={row['No_change_sites']} ?{row['Insufficient_data_sites']}")
    
    print("\n" + "-" * 50)
    print("PANEL B: SEVERE/EXTREME CONDITIONS")
    print("-" * 50)
    
    for timescale in ['SPEI_6', 'SPEI_48']:
        ts_display = timescale.replace('_', '-')
        print(f"\n{ts_display}:")
        
        for condition in ['Severe Dry', 'Severe Wet']:
            cond_display = condition.replace('Severe ', 'Ext. ')
            print(f"\n  {cond_display}:")
            
            for metric in ['WUE', 'WUE_tra']:
                row = df_panel_b[
                    (df_panel_b['SPEI_Timescale'] == timescale) &
                    (df_panel_b['Condition'] == condition) &
                    (df_panel_b['WUE_Metric'] == metric)
                ]
                
                if len(row) > 0:
                    row = row.iloc[0]
                    metric_label = 'WUE_ET' if metric == 'WUE' else 'WUE_T'
                    print(f"    {metric_label} (N={int(row['n_sites'])}): Median = {row['Median_%_Change']:+.1f}%, Mean ± SE = {row['Mean_%_Change']:+.1f}% ± {row['SE_%_Change']:.1f}%")
                    print(f"      Direction: {row['Direction']}")
                    print(f"      Site counts: ↑{row['Increase_sites']} ↓{row['Decrease_sites']} ={row['No_change_sites']} ?{row['Insufficient_data_sites']}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def main():
    # File paths
    site_level_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\figures_SPEI6_SPEI48"
    output_fig_path = os.path.join(output_dir, 'Combined_Figure_TwoPanel.png')
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 100)
    print("COMBINED TWO-PANEL FIGURE WORKFLOW")
    print("Panel (a): All conditions (Dry/Wet) - LEFT")
    print("Panel (b): Severe/Extreme conditions - RIGHT")
    print("Shared calculations using strict double intersection")
    print("Single colorbar for both panels")
    print("=" * 100)
    print("\nThis downstream figure uses updated NN bootstrap CI site-level classifications")
    print("from SPEI_site_level_details_FINAL.csv. The heatmap still shows median percent")
    print("change; direction summaries are based on site-level Increase/Decrease/No change")
    print("counts relative to NN bootstrap CI.")
    print("=" * 100)
    
    # Load and recompute statistics for both panels
    print("\n1. LOADING AND RECOMPUTING STATISTICS...")
    df_panel_a, df_panel_b = load_and_recompute_all_stats(site_level_path)
    
    if len(df_panel_a) == 0 and len(df_panel_b) == 0:
        print("❌ ERROR: No data after recomputation.")
        return None
    
    # Print results summary
    print("\n2. PRINTING RESULTS SUMMARY...")
    print_results_summary(df_panel_a, df_panel_b)
    
    # Create two-panel figure
    print("\n3. CREATING TWO-PANEL FIGURE...")
    fig = create_two_panel_figure(df_panel_a, df_panel_b, output_fig_path)
    
    # Save statistics to CSV
    stats_path_a = os.path.join(output_dir, 'PanelA_AllConditions_Statistics.csv')
    stats_path_b = os.path.join(output_dir, 'PanelB_SevereExtreme_Statistics.csv')
    df_panel_a.to_csv(stats_path_a, index=False)
    df_panel_b.to_csv(stats_path_b, index=False)
    print(f"\n✅ Saved Panel A statistics to: {stats_path_a}")
    print(f"✅ Saved Panel B statistics to: {stats_path_b}")
    
    print("\n" + "=" * 100)
    print("WORKFLOW COMPLETE ✓")
    print("=" * 100)
    
    return fig, df_panel_a, df_panel_b

if __name__ == "__main__":
    fig, stats_a, stats_b = main()


# =============================================================================
# Q2 SUMMARY (UPDATED - uses site-level NN bootstrap CI data)
# =============================================================================

import pandas as pd
import numpy as np

file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"

df = pd.read_csv(file)

metrics = ["WUE", "WUE_tra"]
timescales = ["SPEI_6", "SPEI_48"]

conditions = [
    ("PASS C", "Dry (all)"),
    ("PASS C", "Wet (all)"),
    ("PASS B", "Severe Dry"),
    ("PASS B", "Severe Wet")
]

rows = []

for pass_name, condition in conditions:
    for timescale in timescales:
        subset_condition = df[
            (df["Pass"] == pass_name) &
            (df["Salinity"] == "All") &
            (df["SPEI_Timescale"] == timescale) &
            (df["Condition"] == condition) &
            (df["WUE_Metric"].isin(metrics))
        ].copy()

        # strict double intersection, same as your figure workflow
        wue_sites = set(subset_condition[subset_condition["WUE_Metric"] == "WUE"]["Site"].unique())
        tra_sites = set(subset_condition[subset_condition["WUE_Metric"] == "WUE_tra"]["Site"].unique())
        shared_sites = wue_sites.intersection(tra_sites)

        for metric in metrics:
            sub = subset_condition[
                (subset_condition["WUE_Metric"] == metric) &
                (subset_condition["Site"].isin(shared_sites))
            ].copy()

            site_vals = sub.groupby("Site")["Median_%_Change"].median().dropna().values

            if len(site_vals) == 0:
                continue

            q25 = np.percentile(site_vals, 25)
            q75 = np.percentile(site_vals, 75)
            median = np.median(site_vals)
            mean = np.mean(site_vals)
            sd = np.std(site_vals, ddof=1) if len(site_vals) > 1 else np.nan
            se = sd / np.sqrt(len(site_vals)) if len(site_vals) > 1 else np.nan

            rows.append({
                "Pass": pass_name,
                "Condition": condition,
                "Timescale": timescale,
                "Metric": "WUE_ET" if metric == "WUE" else "WUE_T",
                "N_sites": len(site_vals),
                "Median_%": median,
                "IQR_low_%": q25,
                "IQR_high_%": q75,
                "IQR_width_%": q75 - q25,
                "Mean_%": mean,
                "SD_%": sd,
                "SE_%": se
            })

out = pd.DataFrame(rows)

print("\n" + "=" * 100)
print("Q2 SUMMARY: MEDIAN, IQR, MEAN, SD, SE")
print("UPDATED: Uses site-level NN bootstrap CI data from SPEI_site_level_details_FINAL.csv")
print("=" * 100)
print(out.to_string(index=False))

out_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\Q2_median_IQR_mean_SE_summary.csv"
out.to_csv(out_file, index=False)
print("\nSaved:", out_file)