# -*- coding: utf-8 -*-
"""
FIGURE 4 (SUPPLEMENTARY) WORKFLOW - Freshwater vs Saline vs Upland Median % Change Heatmap
UPDATED: Split into SPEI-6 and SPEI-48 panels for better readability
Q2 methodology - WUE_ET and WUE_T only, strict double intersection
USES: Existing Direction column from Step 1 NN bootstrap CI classifications
STATISTICS: Pairwise Mann-Whitney U tests between ecosystem types
COLORBAR: Zero-gap positioning using make_axes_locatable

PANEL A: SPEI-6 (Short-term) - Dry and Wet conditions (DRY first, then WET)
PANEL B: SPEI-48 (Long-term) - Dry and Wet conditions (DRY first, then WET)
ORDER: Within each condition: Upland, Freshwater, Saline
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
import os
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# =============================================================================
# STYLING - OPTIMIZED FOR READABILITY
# =============================================================================
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 30,
    'axes.titlesize': 0,
    'axes.labelsize': 34,
    'xtick.labelsize': 32,
    'ytick.labelsize': 30,
    'legend.fontsize': 30,
})

# =============================================================================
# 1. DATA LOADING - USING SITE-LEVEL FILE
# =============================================================================
def load_data(filepath):
    """Load and filter site-level data."""
    df = pd.read_csv(filepath)
    
    # NOTE: Do NOT filter Salinity == "All" - this is salinity/ecosystem comparison
    filters = (
        (df['Pass'] == 'PASS C') &
        (df['Salinity'].isin(['Freshwater', 'Saline', 'Upland'])) &
        (df['SPEI_Timescale'].isin(['SPEI_6', 'SPEI_48'])) &
        (df['Condition'].isin(['Dry (all)', 'Wet (all)'])) &
        (df['WUE_Metric'].isin(['WUE', 'WUE_tra']))
    )
    
    filtered_df = df[filters].copy()
    
    print(f"\n[LOAD] Loaded {len(filtered_df)} records after filtering")
    print(f"  Unique sites: {filtered_df['Site'].nunique()}")
    
    return filtered_df

# =============================================================================
# 2. STATISTICAL TESTING - PAIRWISE MANN-WHITNEY U TESTS AT SITE LEVEL
# =============================================================================
def perform_statistical_tests_site_level(filtered_df):
    """
    Perform pairwise Mann-Whitney U tests between ecosystem types at site level.
    Tests: Upland vs Freshwater, Upland vs Saline, Freshwater vs Saline
    Applies strict double intersection within each salinity type separately.
    """
    print("\n" + "=" * 100)
    print("PAIRWISE STATISTICAL TESTS: Upland vs Freshwater vs Saline")
    print("Mann-Whitney U tests comparing site-level Median_%_Change values")
    print("=" * 100)
    
    all_results = []
    
    # Define test combinations
    for timescale in ['SPEI_6', 'SPEI_48']:
        for condition in ['Dry (all)', 'Wet (all)']:
            for metric in ['WUE', 'WUE_tra']:  # Original metric names in data
                
                print(f"\n{'='*70}")
                print(f"Testing: {timescale} | {condition} | {metric}")
                print(f"{'='*70}")
                
                # Dictionary to store site-level values for each salinity type
                salinity_values = {}
                salinity_sites = {}
                
                # For each salinity type, apply strict double intersection
                for salinity in ['Upland', 'Freshwater', 'Saline']:
                    # Get data for this salinity
                    salinity_data = filtered_df[
                        (filtered_df['SPEI_Timescale'] == timescale) &
                        (filtered_df['Condition'] == condition) &
                        (filtered_df['Salinity'] == salinity)
                    ].copy()
                    
                    if len(salinity_data) == 0:
                        print(f"  {salinity}: No data available")
                        continue
                    
                    # Get sites for each metric (WUE and WUE_tra)
                    wue_sites = set(salinity_data[salinity_data['WUE_Metric'] == 'WUE']['Site'].unique())
                    tra_sites = set(salinity_data[salinity_data['WUE_Metric'] == 'WUE_tra']['Site'].unique())
                    
                    # Strict double intersection: sites with BOTH metrics
                    shared_sites = wue_sites.intersection(tra_sites)
                    
                    if len(shared_sites) == 0:
                        print(f"  {salinity}: No sites with both WUE and WUE_tra")
                        continue
                    
                    # Extract site-level Median_%_Change for the specific metric
                    metric_data = salinity_data[
                        (salinity_data['WUE_Metric'] == metric) &
                        (salinity_data['Site'].isin(shared_sites))
                    ].copy()
                    
                    # Collect values
                    values = []
                    for site in shared_sites:
                        site_rows = metric_data[metric_data['Site'] == site]
                        if len(site_rows) > 0:
                            values.append(site_rows['Median_%_Change'].iloc[0])
                    
                    if len(values) > 0:
                        salinity_values[salinity] = values
                        salinity_sites[salinity] = shared_sites
                        print(f"  {salinity}: n={len(values)} sites, median={np.median(values):+.1f}%")
                    else:
                        print(f"  {salinity}: No valid data for metric {metric}")
                
                # Perform pairwise tests if we have at least 2 groups with data
                salinity_list = list(salinity_values.keys())
                
                if len(salinity_list) >= 2:
                    # Define pairwise combinations in new order
                    pairs = []
                    if 'Upland' in salinity_list and 'Freshwater' in salinity_list:
                        pairs.append(('Upland', 'Freshwater'))
                    if 'Upland' in salinity_list and 'Saline' in salinity_list:
                        pairs.append(('Upland', 'Saline'))
                    if 'Freshwater' in salinity_list and 'Saline' in salinity_list:
                        pairs.append(('Freshwater', 'Saline'))
                    
                    for group1, group2 in pairs:
                        vals1 = salinity_values[group1]
                        vals2 = salinity_values[group2]
                        
                        # Perform Mann-Whitney U test
                        stat, p_value = mannwhitneyu(vals1, vals2, alternative='two-sided')
                        
                        # Calculate medians and difference
                        median1 = np.median(vals1)
                        median2 = np.median(vals2)
                        diff = median1 - median2
                        
                        all_results.append({
                            'SPEI_Timescale': timescale,
                            'Condition': condition.replace(' (all)', ''),
                            'WUE_Metric': 'WUE_ET' if metric == 'WUE' else 'WUE_T',
                            'Group1': group1,
                            'Group2': group2,
                            'Group1_n': len(vals1),
                            'Group2_n': len(vals2),
                            'Group1_median': median1,
                            'Group2_median': median2,
                            'Difference_(G1-G2)': diff,
                            'U_statistic': stat,
                            'p_value': p_value
                        })
                        
                        print(f"\n  {group1} vs {group2}:")
                        print(f"    {group1}: n={len(vals1)}, median={median1:+.1f}%")
                        print(f"    {group2}: n={len(vals2)}, median={median2:+.1f}%")
                        print(f"    Difference = {diff:+.1f}%")
                        print(f"    U = {stat:.2f}, p = {p_value:.4f}")
                else:
                    print(f"\n  WARNING: Insufficient groups for pairwise comparison")
    
    # Apply FDR correction across all tests
    if len(all_results) > 0:
        results_df = pd.DataFrame(all_results)
        p_values = results_df['p_value'].values
        rejected, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        
        results_df['q_value'] = q_values
        results_df['Significant'] = ['significant' if rej else 'not significant' for rej in rejected]
        
        print("\n" + "=" * 100)
        print("FDR CORRECTION APPLIED (Benjamini-Hochberg)")
        print(f"  Total tests: {len(results_df)}")
        print(f"  Significant at q<0.05: {sum(rejected)}")
        print("=" * 100)
    else:
        results_df = pd.DataFrame()
        print("\nWARNING: No statistical tests could be performed")
    
    return results_df

# =============================================================================
# 3. APPLY STRICT DOUBLE INTERSECTION AND RECOMPUTE STATISTICS
# =============================================================================
def recompute_with_double_intersection(filtered_df):
    """
    Apply strict double intersection and recompute statistics.
    UPDATED: Direction from existing 'Direction' column (Step 1 NN bootstrap CI)
    """
    print("\n" + "=" * 100)
    print("APPLYING STRICT DOUBLE INTERSECTION (WUE_ET and WUE_T only)")
    print("Direction classification uses existing 'Direction' column from CSV (Step 1 NN bootstrap CI)")
    print("=" * 100)
    
    groups = []
    for condition in ['Dry (all)', 'Wet (all)']:  # DRY first, then WET
        for timescale in ['SPEI_6', 'SPEI_48']:
            for salinity in ['Upland', 'Freshwater', 'Saline']:  # Upland first, then Freshwater, then Saline
                groups.append((condition, timescale, salinity))
    
    recomputed_results = []
    
    for condition, timescale, salinity in groups:
        combo_data = filtered_df[
            (filtered_df['Condition'] == condition) &
            (filtered_df['SPEI_Timescale'] == timescale) &
            (filtered_df['Salinity'] == salinity)
        ].copy()
        
        if len(combo_data) == 0:
            continue
        
        # Get sites for each metric
        wue_sites = set(combo_data[combo_data['WUE_Metric'] == 'WUE']['Site'].unique())
        tra_sites = set(combo_data[combo_data['WUE_Metric'] == 'WUE_tra']['Site'].unique())
        
        # Strict double intersection
        shared_sites = wue_sites.intersection(tra_sites)
        
        if len(shared_sites) == 0:
            continue
        
        # Process each metric
        for metric in ['WUE', 'WUE_tra']:
            metric_data = combo_data[
                (combo_data['WUE_Metric'] == metric) &
                (combo_data['Site'].isin(shared_sites))
            ].copy()
            
            if len(metric_data) == 0:
                continue
            
            site_medians = []
            site_increase = 0
            site_decrease = 0
            site_no_change = 0
            site_insufficient = 0
            
            for site in shared_sites:
                site_rows = metric_data[metric_data['Site'] == site]
                if len(site_rows) == 0:
                    continue
                
                median_val = site_rows['Median_%_Change'].median()
                site_medians.append(median_val)
                
                # Use existing Direction column from CSV
                direction = site_rows['Direction'].iloc[0] if 'Direction' in site_rows.columns else 'Unknown'
                
                if direction == 'Increase':
                    site_increase += 1
                elif direction == 'Decrease':
                    site_decrease += 1
                elif direction == 'No change':
                    site_no_change += 1
                elif direction == 'Insufficient data':
                    site_insufficient += 1
            
            n_sites = len(site_medians)
            valid_sites = site_increase + site_decrease + site_no_change
            
            if n_sites > 0:
                median_pct = np.median(site_medians)
                
                # Calculate percentages based on valid sites
                pct_increase = 100 * site_increase / valid_sites if valid_sites > 0 else 0
                pct_decrease = 100 * site_decrease / valid_sites if valid_sites > 0 else 0
                pct_no_change = 100 * site_no_change / valid_sites if valid_sites > 0 else 0
                pct_insufficient = 100 * site_insufficient / n_sites if n_sites > 0 else 0
                
                recomputed_results.append({
                    'Condition': condition.replace(' (all)', ''),
                    'SPEI_Timescale': timescale,
                    'Salinity': salinity,
                    'WUE_Metric': 'WUE_ET' if metric == 'WUE' else 'WUE_T',
                    'n_sites': n_sites,
                    'Valid_sites': valid_sites,
                    'Median_%_Change': median_pct,
                    'Increase_n': site_increase,
                    'Increase_pct': pct_increase,
                    'Decrease_n': site_decrease,
                    'Decrease_pct': pct_decrease,
                    'No_change_n': site_no_change,
                    'No_change_pct': pct_no_change,
                    'Insufficient_n': site_insufficient,
                    'Insufficient_pct': pct_insufficient
                })
    
    results_df = pd.DataFrame(recomputed_results)
    print(f"\n[RECOMPUTED] {len(results_df)} records")
    
    return results_df

# =============================================================================
# 4. CREATE CUSTOM COLORMAP
# =============================================================================
def create_figure1_colormap():
    """Create the stitched gray/green colormap."""
    n_colors = 256
    gray_colors = plt.cm.Greys_r(np.linspace(0.2, 0.8, n_colors // 2))
    green_colors = plt.cm.Greens(np.linspace(0.3, 0.9, n_colors // 2))
    return ListedColormap(np.vstack((gray_colors, green_colors)))

# =============================================================================
# 5. BUILD DATA MATRIX FOR A SINGLE TIMESCALE
# =============================================================================
def build_data_matrix_for_timescale(recomputed_df, timescale):
    """Build data matrix for a specific timescale."""
    metrics = ['WUE_ET', 'WUE_T']
    
    # Row order (top → bottom): DRY then WET, each with Upland→Freshwater→Saline
    row_order = [
        ('Dry', 'Upland'),
        ('Dry', 'Freshwater'),
        ('Dry', 'Saline'),
        ('Wet', 'Upland'),
        ('Wet', 'Freshwater'),
        ('Wet', 'Saline'),
    ]
    
    row_data = []
    row_labels = []
    n_sites_list = []
    salinity_list = []
    
    for condition, salinity in row_order:
        row_vals = []
        n_sites_val = None
        
        for metric in metrics:
            subset = recomputed_df[
                (recomputed_df['Condition'] == condition) &
                (recomputed_df['SPEI_Timescale'] == timescale) &
                (recomputed_df['Salinity'] == salinity) &
                (recomputed_df['WUE_Metric'] == metric)
            ]
            
            if len(subset) > 0:
                row_vals.append(subset.iloc[0]['Median_%_Change'])
                if n_sites_val is None:
                    n_sites_val = subset.iloc[0]['n_sites']
            else:
                row_vals.append(np.nan)
        
        row_data.append(row_vals)
        n_sites_list.append(n_sites_val if n_sites_val is not None else 0)
        salinity_list.append(salinity)
        row_labels.append(f"{condition}")
    
    return np.array(row_data, dtype=float), row_labels, n_sites_list, salinity_list

# =============================================================================
# 6. CREATE INDIVIDUAL HEATMAP PANEL
# =============================================================================
def create_heatmap_panel(recomputed_df, timescale, salinity_colors, ax, panel_label, vmin, vmax):
    """Create a single heatmap panel."""
    
    # Build data matrix for this timescale
    data_matrix, row_labels, n_sites, salinity_list = build_data_matrix_for_timescale(recomputed_df, timescale)
    n_rows = len(data_matrix)
    
    # Plot heatmap with fixed color scale
    custom_cmap = create_figure1_colormap()
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    im = ax.imshow(data_matrix, cmap=custom_cmap, norm=norm, aspect='auto', origin='upper')
    
    # Add cell annotations
    fontsize_cell = 32
    for i in range(n_rows):
        for j in range(2):
            val = data_matrix[i, j]
            if np.isnan(val):
                continue
            val_str = f"{val:+.0f}%"
            # Determine text color based on value
            if abs(val) >= abs(vmax) * 0.6:
                txt_color = 'white'
            else:
                txt_color = 'black'
            ax.text(j, i, val_str, ha='center', va='center',
                   fontsize=fontsize_cell, fontweight='bold', color=txt_color)
    
    # Salinity left strips - Upland (purple), Freshwater (blue), Saline (orange)
    strip_x0 = -0.45
    strip_w = 0.22
    ax.set_xlim(strip_x0, 2.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    
    for i, sal in enumerate(salinity_list):
        ax.add_patch(Rectangle((strip_x0, i-0.5), strip_w, 1,
                              facecolor=salinity_colors[sal],
                              edgecolor='white', linewidth=2))
    
    # Axis labels
    metric_labels = ['WUE$_{ET}$', 'WUE$_T$']
    ax.set_xticks(range(2))
    ax.set_xticklabels(metric_labels, fontweight='bold', fontsize=25)
    
    # Y-axis labels with sample sizes
    y_labels = []
    for label, n in zip(row_labels, n_sites):
        n_str = str(int(n)) if n > 0 else 'NA'
        y_labels.append(f"{label} (N={n_str})")
    
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(y_labels, fontweight='bold', fontsize=25)
    
    # Separator between Dry and Wet (after row 3, index 2.5)
    ax.hlines(y=2.5, xmin=-0.45, xmax=1.5, color='black', linewidth=3, alpha=0.8)
    
    # Vertical separators between metrics
    for x in [0.5, 1.5]:
        ax.axvline(x=x, color='gray', linewidth=1.5, alpha=0.5, linestyle=':')
    
    # Clean up spines
    for spine in ['top', 'right', 'bottom']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)
    
    # Title
    ts_display = 'SPEI-6 (Short-term)' if timescale == 'SPEI_6' else 'SPEI-48 (Long-term)'
    ax.set_title(ts_display, fontsize=30, fontweight='bold', pad=15, x=0.3, loc='center')
    
    # Panel label
    ax.text(-0.35, 1.08, panel_label, transform=ax.transAxes,
            fontsize=48, fontweight='bold',
            verticalalignment='top', horizontalalignment='left',
            zorder=100)
    
    return im

# =============================================================================
# 7. CREATE TWO-PANEL FIGURE
# =============================================================================
def create_two_panel_figure(recomputed_df, save_path):
    """Create a figure with two separate heatmap panels (SPEI-6 and SPEI-48)."""
    
    if len(recomputed_df) == 0:
        print("ERROR: No data for Figure 4!")
        return None
    
    # Find global color scale limits across both timescales
    all_values = recomputed_df['Median_%_Change'].values
    absmax = np.max(np.abs(all_values))
    vmax = absmax * 1.1
    vmin = -vmax
    
    print(f"Color scale: vmin={vmin:.1f}, vmax={vmax:.1f}")
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # Salinity colors - Updated order for legend: Upland, Freshwater, Saline
    salinity_colors = {
        'Upland': '#800080',      # Purple
        'Freshwater': '#0000FF',  # Blue
        'Saline': '#FFA500'       # Orange
    }
    
    # Create individual panels with FIXED color scale
    im1 = create_heatmap_panel(recomputed_df, 'SPEI_6', salinity_colors, ax1, '(a)', vmin, vmax)
    im2 = create_heatmap_panel(recomputed_df, 'SPEI_48', salinity_colors, ax2, '(b)', vmin, vmax)
    
    # Create colorbar
    
    pos = ax2.get_position()
    colorbar_offset = -0.08  # Change this value to move left/right (negative = left, positive = right)
    cax = fig.add_axes([pos.x1 + colorbar_offset, pos.y0, 0.015, pos.height])
    cbar = plt.colorbar(im2, cax=cax)
    cbar.set_label('Median % Change (%)', fontweight='bold', fontsize=25)
    cbar.ax.tick_params(labelsize=28)
    
    # Legend for salinity - ORDER: Upland, Freshwater, Saline
    legend_elements = [
        Patch(facecolor='#800080', edgecolor='black', label='Upland', linewidth=2),
        Patch(facecolor='#0000FF', edgecolor='black', label='Freshwater', linewidth=2),
        Patch(facecolor='#FFA500', edgecolor='black', label='Saline', linewidth=2),
    ]
    
    fig.legend(handles=legend_elements, loc='upper center',
               bbox_to_anchor=(0.5, 1.04), ncol=3, fontsize=30,
               frameon=True, framealpha=0.95, borderpad=0.8,
               handlelength=1.5, handleheight=1.5)
    
    # Adjust layout
    fig.subplots_adjust(left=0.10, right=0.92, top=0.78, bottom=0.12, wspace=0.04)
    
    # Save figure
    png_path = f"{save_path}.png"
    pdf_path = f"{save_path}.pdf"
    
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSUCCESS: Figure saved to:")
    print(f"   PNG: {png_path} (300 DPI)")
    print(f"   PDF: {pdf_path}")
    
    plt.show()
    return fig

# =============================================================================
# 8. PRINT RESULTS
# =============================================================================
def print_results(recomputed_df, stats_results):
    """Print formatted results for manuscript including statistical tests."""
    print("\n" + "=" * 100)
    print("RESULTS - Upland vs Freshwater vs Saline Comparison (PASS C)")
    print("Q2 methodology: WUE_ET and WUE_T only, strict double intersection")
    print("Direction classification uses Step 1 NN bootstrap CI (Direction column)")
    print("=" * 100)
    
    if len(recomputed_df) == 0:
        return
    
    print("\nMEDIAN PERCENT CHANGE BY ECOSYSTEM TYPE:")
    print("-" * 70)
    
    # Sort for display
    for timescale in ['SPEI_6', 'SPEI_48']:
        ts_display = timescale.replace('_', '-')
        print(f"\n{ts_display}:")
        
        for condition in ['Dry', 'Wet']:  # DRY first, then WET
            print(f"  {condition}:")
            for salinity in ['Upland', 'Freshwater', 'Saline']:
                subset = recomputed_df[
                    (recomputed_df['SPEI_Timescale'] == timescale) &
                    (recomputed_df['Condition'] == condition) &
                    (recomputed_df['Salinity'] == salinity)
                ]
                if len(subset) > 0:
                    wue_et = subset[subset['WUE_Metric'] == 'WUE_ET']['Median_%_Change'].values[0]
                    wue_t = subset[subset['WUE_Metric'] == 'WUE_T']['Median_%_Change'].values[0]
                    n_sites = subset.iloc[0]['n_sites']
                    print(f"    {salinity:12s} (n={n_sites}): WUE_ET={wue_et:+.1f}%, WUE_T={wue_t:+.1f}%")
    
    print("\nRESPONSE DIRECTION (% Increase, % Decrease, % No Change):")
    print("-" * 70)
    
    for _, row in recomputed_df.iterrows():
        ts_display = row['SPEI_Timescale'].replace('_', '-')
        print(f"  {ts_display} {row['Condition']} {row['Salinity']:12s} {row['WUE_Metric']:6s}: "
              f"Inc={row['Increase_pct']:.0f}%, Dec={row['Decrease_pct']:.0f}%, NC={row['No_change_pct']:.0f}% (n={row['n_sites']})")
    
    # Print statistical test results
    if len(stats_results) > 0:
        print("\n" + "=" * 100)
        print("PAIRWISE STATISTICAL TEST RESULTS (Mann-Whitney U with FDR correction):")
        print("-" * 70)
        print("Tests compare site-level median percent change distributions")
        print("between ecosystem types under the same conditions.")
        print("-" * 70)
        
        for _, row in stats_results.iterrows():
            sig_marker = "✓" if row['Significant'] == 'significant' else "✗"
            print(f"\n  {sig_marker} {row['SPEI_Timescale']} {row['Condition']} {row['WUE_Metric']}: {row['Group1']} vs {row['Group2']}")
            print(f"      {row['Group1']}: n={row['Group1_n']:.0f}, median={row['Group1_median']:+.1f}%")
            print(f"      {row['Group2']}: n={row['Group2_n']:.0f}, median={row['Group2_median']:+.1f}%")
            print(f"      Difference = {row['Difference_(G1-G2)']:+.1f}%")
            print(f"      U = {row['U_statistic']:.2f}, p = {row['p_value']:.4f}, q = {row['q_value']:.4f}")
            print(f"      → {row['Significant'].upper()} at q<0.05")
    
    print("\n" + "=" * 100)
    print("INTERPRETATION:")
    print("-" * 70)
    print("  • Strict double intersection applied (sites with both WUE_ET and WUE_T)")
    print("  • Direction classification uses existing 'Direction' column from Step 1")
    print("  • Direction reflects site-level NN bootstrap 95% CI threshold test")
    print("  • Split into SPEI-6 and SPEI-48 panels for improved readability")
    print("  • Row order: DRY first (top 3 rows), then WET (bottom 3 rows)")
    print("  • Within each condition: Upland, Freshwater, Saline")
    print("  • Statistical tests: Two-sided Mann-Whitney U tests with FDR correction (q<0.05)")
    print("=" * 100)

# =============================================================================
# 9. SAVE OUTPUTS
# =============================================================================
def save_outputs(recomputed_df, stats_results, output_dir):
    """Save CSV outputs."""
    summary_path = os.path.join(output_dir, 'Table_Upland_FW_Saline_ResponseSummary_WUEET_WUET.csv')
    recomputed_df.to_csv(summary_path, index=False)
    
    if len(stats_results) > 0:
        stats_path = os.path.join(output_dir, 'Table_Upland_FW_Saline_PairwiseStats_MannWhitney.csv')
        stats_results.to_csv(stats_path, index=False)
        print(f"\n✅ Saved response summary to: {summary_path}")
        print(f"✅ Saved statistical test results to: {stats_path}")
    else:
        print(f"\n✅ Saved response summary to: {summary_path}")

# =============================================================================
# 10. MAIN WORKFLOW
# =============================================================================
def main():
    """Main execution function."""
    # File paths
    input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\figures_SPEI6_SPEI48"
    output_base = os.path.join(output_dir, 'FigSx_Upland_FW_Saline_MedianPctChange_WUEET_WUET')
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 100)
    print("FIGURE 4 (SUPPLEMENTARY): UPLAND VS FRESHWATER VS SALINE HEATMAP")
    print("UPDATED: Split into SPEI-6 and SPEI-48 panels for better readability")
    print("Q2 methodology: WUE_ET and WUE_T only")
    print("=" * 100)
    print("\nPANEL A: SPEI-6 (Short-term)")
    print("  Row order: DRY (Upland → Freshwater → Saline) then WET (Upland → Freshwater → Saline)")
    print("\nPANEL B: SPEI-48 (Long-term)")
    print("  Row order: DRY (Upland → Freshwater → Saline) then WET (Upland → Freshwater → Saline)")
    print("\nCOLORS: Figure 1 gray/green for % change, Purple/Blue/Orange for salinity")
    print("METHOD: Strict double intersection (WUE_ET and WUE_T only)")
    print("DIRECTION: Uses existing 'Direction' column from Step 1 NN bootstrap CI")
    print("STATISTICS: Pairwise Mann-Whitney U tests with FDR correction")
    print("=" * 100)
    
    # Step 1: Load data
    print("\n[STEP 1] Loading site-level data...")
    filtered_data = load_data(input_file)
    
    if len(filtered_data) == 0:
        print("ERROR: No data found!")
        return
    
    # Step 2: Perform statistical tests BEFORE aggregation
    print("\n[STEP 2] Performing pairwise statistical tests...")
    stats_results = perform_statistical_tests_site_level(filtered_data)
    
    # Step 3: Apply double intersection and recompute statistics for heatmap
    print("\n[STEP 3] Applying strict double intersection...")
    recomputed_df = recompute_with_double_intersection(filtered_data)
    
    if len(recomputed_df) == 0:
        print("ERROR: No data after recomputation!")
        return
    
    # Step 4: Create two-panel figure
    print("\n[STEP 4] Creating two-panel heatmap...")
    fig = create_two_panel_figure(recomputed_df, output_base)
    
    # Step 5: Print results with statistics
    print_results(recomputed_df, stats_results)
    
    # Step 6: Save outputs
    print("\n[STEP 5] Saving outputs...")
    save_outputs(recomputed_df, stats_results, output_dir)
    
    print("\n" + "=" * 100)
    print("WORKFLOW COMPLETE")
    print("=" * 100)
    print("\n✅ KEY UPDATES:")
    print("  1. Split into two panels: SPEI-6 (a) and SPEI-48 (b)")
    print("  2. WUE_E completely removed")
    print("  3. Strict double intersection (WUE_ET and WUE_T only)")
    print("  4. Direction uses existing 'Direction' column from Step 1")
    print("  5. Fixed black separator line - now stays within heatmap")
    print("  6. Reduced space between panels")
    print("  7. Colorbar properly positioned next to Panel B")
    print("  8. Added pairwise Mann-Whitney U tests with FDR correction")
    print("  9. NEW ORDER: DRY first, then WET")
    print("  10. NEW ORDER within each condition: Upland, Freshwater, Saline")
    print("  11. Legend order: Upland, Freshwater, Saline")
    
    return fig, recomputed_df, stats_results

# =============================================================================
# EXECUTION
# =============================================================================
if __name__ == "__main__":
    main()