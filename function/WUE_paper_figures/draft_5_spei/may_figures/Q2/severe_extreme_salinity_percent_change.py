# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:27:56 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
DIAGNOSTIC WORKFLOW: SEVERE/EXTREME CONDITIONS ONLY (PASS B)
- Uses PASS B (Severe Dry, Severe Wet) - NOT PASS C
- Freshwater vs Saline vs Upland comparison
- Strict double intersection (WUE_ET and WUE_T)
- Direction from Step 1 NN bootstrap CI
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
import os
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# =============================================================================
# STYLING
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
# 1. DATA LOADING - PASS B ONLY (SEVERE/EXTREME CONDITIONS)
# =============================================================================
def load_data_severe(filepath):
    """Load and filter site-level data for SEVERE/EXTREME conditions (PASS B)."""
    df = pd.read_csv(filepath)
    
    # CRITICAL: Use PASS B for severe/extreme conditions
    filters = (
        (df['Pass'] == 'PASS B') &  # <- CHANGED from PASS C to PASS B
        (df['Salinity'].isin(['Freshwater', 'Saline', 'Upland'])) &
        (df['SPEI_Timescale'].isin(['SPEI_6', 'SPEI_48'])) &
        (df['Condition'].isin(['Severe Dry', 'Severe Wet'])) &  # Severe only
        (df['WUE_Metric'].isin(['WUE', 'WUE_tra']))
    )
    
    filtered_df = df[filters].copy()
    
    print(f"\n[LOAD SEVERE] Loaded {len(filtered_df)} records after filtering")
    print(f"  Unique sites: {filtered_df['Site'].nunique()}")
    print(f"  Conditions found: {filtered_df['Condition'].unique().tolist()}")
    print(f"  Salinities found: {filtered_df['Salinity'].unique().tolist()}")
    
    return filtered_df

# =============================================================================
# 2. CREATE CUSTOM COLORMAP
# =============================================================================
def create_colormap():
    """Create the stitched gray/green colormap."""
    n_colors = 256
    gray_colors = plt.cm.Greys_r(np.linspace(0.2, 0.8, n_colors // 2))
    green_colors = plt.cm.Greens(np.linspace(0.3, 0.9, n_colors // 2))
    return ListedColormap(np.vstack((gray_colors, green_colors)))

# =============================================================================
# 3. STATISTICAL TESTS - PAIRWISE MANN-WHITNEY U (SEVERE ONLY)
# =============================================================================
def perform_statistical_tests_severe(filtered_df):
    """
    Perform pairwise Mann-Whitney U tests between ecosystem types.
    USES PASS B: Severe Dry, Severe Wet only
    """
    print("\n" + "=" * 100)
    print("DIAGNOSTIC: SEVERE/EXTREME CONDITIONS ONLY (PASS B)")
    print("PAIRWISE STATISTICAL TESTS: Freshwater vs Saline vs Upland")
    print("Mann-Whitney U tests comparing site-level Median_%_Change values")
    print("Conditions: Severe Dry, Severe Wet")
    print("=" * 100)
    
    all_results = []
    
    for timescale in ['SPEI_6', 'SPEI_48']:
        for condition in ['Severe Wet', 'Severe Dry']:
            for metric in ['WUE', 'WUE_tra']:
                
                print(f"\n{'='*70}")
                print(f"Testing: {timescale} | {condition} | {metric}")
                print(f"{'='*70}")
                
                salinity_values = {}
                
                for salinity in ['Freshwater', 'Saline', 'Upland']:
                    salinity_data = filtered_df[
                        (filtered_df['SPEI_Timescale'] == timescale) &
                        (filtered_df['Condition'] == condition) &
                        (filtered_df['Salinity'] == salinity)
                    ].copy()
                    
                    if len(salinity_data) == 0:
                        print(f"  {salinity}: No data available")
                        continue
                    
                    # Strict double intersection
                    wue_sites = set(salinity_data[salinity_data['WUE_Metric'] == 'WUE']['Site'].unique())
                    tra_sites = set(salinity_data[salinity_data['WUE_Metric'] == 'WUE_tra']['Site'].unique())
                    shared_sites = wue_sites.intersection(tra_sites)
                    
                    if len(shared_sites) == 0:
                        print(f"  {salinity}: No shared sites")
                        continue
                    
                    # Extract values
                    metric_data = salinity_data[
                        (salinity_data['WUE_Metric'] == metric) &
                        (salinity_data['Site'].isin(shared_sites))
                    ].copy()
                    
                    values = []
                    for site in shared_sites:
                        site_rows = metric_data[metric_data['Site'] == site]
                        if len(site_rows) > 0:
                            values.append(site_rows['Median_%_Change'].iloc[0])
                    
                    if len(values) > 0:
                        salinity_values[salinity] = values
                        print(f"  {salinity}: n={len(values)} sites, median={np.median(values):+.1f}%")
                
                # Pairwise tests
                salinity_list = list(salinity_values.keys())
                if len(salinity_list) >= 2:
                    pairs = []
                    if 'Freshwater' in salinity_list and 'Saline' in salinity_list:
                        pairs.append(('Freshwater', 'Saline'))
                    if 'Freshwater' in salinity_list and 'Upland' in salinity_list:
                        pairs.append(('Freshwater', 'Upland'))
                    if 'Saline' in salinity_list and 'Upland' in salinity_list:
                        pairs.append(('Saline', 'Upland'))
                    
                    for group1, group2 in pairs:
                        vals1 = salinity_values[group1]
                        vals2 = salinity_values[group2]
                        
                        stat, p_value = mannwhitneyu(vals1, vals2, alternative='two-sided')
                        median1 = np.median(vals1)
                        median2 = np.median(vals2)
                        diff = median1 - median2
                        
                        all_results.append({
                            'Test_Type': 'SEVERE_ONLY_PASS_B',
                            'SPEI_Timescale': timescale,
                            'Condition': condition.replace('Severe ', ''),
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
                        print(f"    U = {stat:.2f}, p = {p_value:.4f}")
                        print(f"    Median {group1}: {median1:+.1f}%, Median {group2}: {median2:+.1f}%")
    
    # FDR correction
    if len(all_results) > 0:
        results_df = pd.DataFrame(all_results)
        p_values = results_df['p_value'].values
        rejected, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        
        results_df['q_value'] = q_values
        results_df['Significant'] = ['significant' if rej else 'not significant' for rej in rejected]
        
        print("\n" + "=" * 100)
        print(f"FDR CORRECTION: {len(results_df)} tests, {sum(rejected)} significant at q<0.05")
        print("=" * 100)
    
    return results_df

# =============================================================================
# 4. RECOMPUTE WITH DOUBLE INTERSECTION (SEVERE ONLY)
# =============================================================================
def recompute_severe(filtered_df):
    """Recompute statistics for severe/extreme conditions using Direction column."""
    print("\n" + "=" * 100)
    print("RECOMPUTING FOR SEVERE/EXTREME CONDITIONS (PASS B)")
    print("Using Direction column from Step 1 NN bootstrap CI")
    print("=" * 100)
    
    results = []
    
    for condition in ['Severe Wet', 'Severe Dry']:
        for timescale in ['SPEI_6', 'SPEI_48']:
            for salinity in ['Freshwater', 'Saline', 'Upland']:
                
                combo_data = filtered_df[
                    (filtered_df['Condition'] == condition) &
                    (filtered_df['SPEI_Timescale'] == timescale) &
                    (filtered_df['Salinity'] == salinity)
                ].copy()
                
                if len(combo_data) == 0:
                    continue
                
                # Strict double intersection
                wue_sites = set(combo_data[combo_data['WUE_Metric'] == 'WUE']['Site'].unique())
                tra_sites = set(combo_data[combo_data['WUE_Metric'] == 'WUE_tra']['Site'].unique())
                shared_sites = wue_sites.intersection(tra_sites)
                
                if len(shared_sites) == 0:
                    continue
                
                for metric in ['WUE', 'WUE_tra']:
                    metric_data = combo_data[
                        (combo_data['WUE_Metric'] == metric) &
                        (combo_data['Site'].isin(shared_sites))
                    ].copy()
                    
                    site_medians = []
                    site_inc = site_dec = site_nc = site_insuff = 0
                    
                    for site in shared_sites:
                        site_rows = metric_data[metric_data['Site'] == site]
                        if len(site_rows) == 0:
                            continue
                        
                        site_medians.append(site_rows['Median_%_Change'].median())
                        
                        direction = site_rows['Direction'].iloc[0] if 'Direction' in site_rows.columns else 'Unknown'
                        if direction == 'Increase':
                            site_inc += 1
                        elif direction == 'Decrease':
                            site_dec += 1
                        elif direction == 'No change':
                            site_nc += 1
                        elif direction == 'Insufficient data':
                            site_insuff += 1
                    
                    n_sites = len(site_medians)
                    valid = site_inc + site_dec + site_nc
                    
                    if n_sites > 0:
                        results.append({
                            'Condition': condition.replace('Severe ', ''),
                            'SPEI_Timescale': timescale,
                            'Salinity': salinity,
                            'WUE_Metric': 'WUE_ET' if metric == 'WUE' else 'WUE_T',
                            'n_sites': n_sites,
                            'Median_%_Change': np.median(site_medians),
                            'Increase_pct': 100 * site_inc / valid if valid > 0 else 0,
                            'Decrease_pct': 100 * site_dec / valid if valid > 0 else 0,
                            'No_change_pct': 100 * site_nc / valid if valid > 0 else 0,
                        })
    
    df_results = pd.DataFrame(results)
    print(f"\n[SEVERE RECOMPUTED] {len(df_results)} records")
    return df_results

# =============================================================================
# 5. CREATE HEATMAP (SEVERE ONLY)
# =============================================================================
def create_heatmap_severe(recomputed_df, save_path):
    """Create two-panel heatmap for severe/extreme conditions."""
    
    if len(recomputed_df) == 0:
        print("ERROR: No data for severe heatmap")
        return None
    
    # Color limits
    all_vals = recomputed_df['Median_%_Change'].values
    absmax = np.max(np.abs(all_vals))
    vmax = absmax * 1.1
    vmin = -vmax
    
    print(f"Color scale: vmin={vmin:.1f}, vmax={vmax:.1f}")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 10))
    
    salinity_colors = {'Freshwater': '#0000FF', 'Saline': '#FFA500', 'Upland': '#800080'}
    row_order = [('Wet', 'Freshwater'), ('Wet', 'Saline'), ('Wet', 'Upland'),
                 ('Dry', 'Freshwater'), ('Dry', 'Saline'), ('Dry', 'Upland')]
    
    for ax_idx, (ax, timescale) in enumerate([(ax1, 'SPEI_6'), (ax2, 'SPEI_48')]):
        data_matrix = []
        n_sites_list = []
        salinity_list = []
        
        for condition, salinity in row_order:
            row_vals = []
            n_val = None
            for metric in ['WUE_ET', 'WUE_T']:
                subset = recomputed_df[
                    (recomputed_df['Condition'] == condition) &
                    (recomputed_df['SPEI_Timescale'] == timescale) &
                    (recomputed_df['Salinity'] == salinity) &
                    (recomputed_df['WUE_Metric'] == metric)
                ]
                if len(subset) > 0:
                    row_vals.append(subset.iloc[0]['Median_%_Change'])
                    if n_val is None:
                        n_val = subset.iloc[0]['n_sites']
                else:
                    row_vals.append(np.nan)
            data_matrix.append(row_vals)
            n_sites_list.append(n_val if n_val else 0)
            salinity_list.append(salinity)
        
        data_matrix = np.array(data_matrix)
        
        # Plot
        cmap = create_colormap()
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        im = ax.imshow(data_matrix, cmap=cmap, norm=norm, aspect='auto', origin='upper')
        
        # Annotations
        for i in range(6):
            for j in range(2):
                if not np.isnan(data_matrix[i, j]):
                    color = 'white' if abs(data_matrix[i, j]) > vmax * 0.6 else 'black'
                    ax.text(j, i, f"{data_matrix[i, j]:+.0f}%", ha='center', va='center',
                           fontsize=32, fontweight='bold', color=color)
        
        # Salinity strips
        strip_x0 = -0.45
        for i, sal in enumerate(salinity_list):
            ax.add_patch(Rectangle((strip_x0, i-0.5), 0.22, 1,
                          facecolor=salinity_colors[sal], edgecolor='white', linewidth=2))
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['WUE$_{ET}$', 'WUE$_T$'], fontweight='bold', fontsize=25)
        
        y_labels = [f"{cond}\n(n={n})" if n > 0 else f"{cond}\n(NA)" 
                    for (cond, _), n in zip(row_order, n_sites_list)]
        ax.set_yticks(range(6))
        ax.set_yticklabels(y_labels, fontweight='bold', fontsize=22)
        ax.set_xlim(strip_x0, 2.5)
        
        # Separator
        ax.hlines(y=2.5, xmin=strip_x0, xmax=1.5, color='black', linewidth=3)
        
        # Title
        ts_label = 'SPEI-6 (Short-term)' if timescale == 'SPEI_6' else 'SPEI-48 (Long-term)'
        ax.set_title(f"{ts_label}\nSEVERE/EXTREME CONDITIONS", fontsize=28, fontweight='bold', pad=15)
        ax.text(-0.35, 1.08, '(a)' if ax_idx == 0 else '(b)', transform=ax.transAxes,
               fontsize=48, fontweight='bold', va='top', ha='left')
        
        for spine in ['top', 'right', 'bottom']:
            ax.spines[spine].set_visible(False)
    
    # Colorbar
    pos = ax2.get_position()
    cax = fig.add_axes([pos.x1 + 0, pos.y0, 0.015, pos.height])
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('Median % Change (%)', fontweight='bold', fontsize=25)
    cbar.ax.tick_params(labelsize=28)
    
    # Legend
    legend_elements = [
        Patch(facecolor='#0000FF', edgecolor='black', label='Freshwater'),
        Patch(facecolor='#FFA500', edgecolor='black', label='Saline'),
        Patch(facecolor='#800080', edgecolor='black', label='Upland'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.04),
               ncol=3, fontsize=30, frameon=True, framealpha=0.95)
    
    fig.subplots_adjust(left=0.10, right=0.92, top=0.78, bottom=0.12, wspace=0.04)
    
    # Save
    fig.savefig(f"{save_path}_SEVERE_ONLY.png", dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(f"{save_path}_SEVERE_ONLY.pdf", bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {save_path}_SEVERE_ONLY.png")
    
    plt.show()
    return fig

# =============================================================================
# 6. PRINT RESULTS
# =============================================================================
def print_results_severe(recomputed_df, stats_df):
    """Print results for severe/extreme conditions."""
    print("\n" + "=" * 100)
    print("SEVERE/EXTREME CONDITIONS ONLY (PASS B) - RESULTS")
    print("=" * 100)
    
    print("\nMEDIAN PERCENT CHANGE BY ECOSYSTEM TYPE:")
    for timescale in ['SPEI_6', 'SPEI_48']:
        print(f"\n{timescale.replace('_', '-')}:")
        for condition in ['Wet', 'Dry']:
            print(f"  {condition} (Severe):")
            for salinity in ['Freshwater', 'Saline', 'Upland']:
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
    
    if len(stats_df) > 0:
        print("\n" + "=" * 100)
        print("STATISTICAL TESTS (Mann-Whitney U with FDR correction):")
        for _, row in stats_df.iterrows():
            sig = "✓ SIGNIFICANT" if row['Significant'] == 'significant' else "✗ NOT significant"
            print(f"\n  {sig} | {row['SPEI_Timescale']} {row['Condition']} {row['WUE_Metric']}")
            print(f"    {row['Group1']} vs {row['Group2']}: p={row['p_value']:.4f}, q={row['q_value']:.4f}")

# =============================================================================
# 7. MAIN DIAGNOSTIC WORKFLOW (SEVERE ONLY)
# =============================================================================
def main_severe_diagnostic():
    """Main diagnostic workflow for severe/extreme conditions (PASS B)."""
    
    input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\figures_SPEI6_SPEI48"
    output_base = os.path.join(output_dir, 'DIAGNOSTIC_SEVERE_ONLY_PASS_B')
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 100)
    print("DIAGNOSTIC WORKFLOW: SEVERE/EXTREME CONDITIONS ONLY")
    print("Using PASS B (Severe Dry, Severe Wet)")
    print("=" * 100)
    
    # Step 1: Load data with PASS B
    print("\n[STEP 1] Loading data (PASS B only)...")
    filtered_data = load_data_severe(input_file)
    
    if len(filtered_data) == 0:
        print("ERROR: No severe/extreme data found in PASS B")
        return None, None, None
    
    # Step 2: Statistical tests
    print("\n[STEP 2] Performing statistical tests...")
    stats_results = perform_statistical_tests_severe(filtered_data)
    
    # Step 3: Recompute for heatmap
    print("\n[STEP 3] Recomputing statistics...")
    recomputed_df = recompute_severe(filtered_data)
    
    if len(recomputed_df) == 0:
        print("ERROR: No data after recomputation")
        return None, None, None
    
    # Step 4: Create heatmap
    print("\n[STEP 4] Creating heatmap...")
    fig = create_heatmap_severe(recomputed_df, output_base)
    
    # Step 5: Print results
    print_results_severe(recomputed_df, stats_results)
    
    # Step 6: Save CSV outputs
    recomputed_df.to_csv(os.path.join(output_dir, 'DIAGNOSTIC_SEVERE_ONLY_ResponseSummary.csv'), index=False)
    if len(stats_results) > 0:
        stats_results.to_csv(os.path.join(output_dir, 'DIAGNOSTIC_SEVERE_ONLY_Stats.csv'), index=False)
    
    print("\n" + "=" * 100)
    print("DIAGNOSTIC WORKFLOW COMPLETE")
    print("=" * 100)
    
    return fig, recomputed_df, stats_results


# =============================================================================
# RUN
# =============================================================================
if __name__ == "__main__":
    fig, df_results, stats = main_severe_diagnostic()