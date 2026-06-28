# -*- coding: utf-8 -*-
"""
COMPLETE 4-PANEL BAR CHART FIGURE - UPDATED
- TOP ROW: Panel (a) and (b) - All conditions (PASS C) for SPEI-6 and SPEI-48
- BOTTOM ROW: Panel (c) and (d) - Severe/Extreme conditions (PASS B) for SPEI-6 and SPEI-48
- Shared calculations using strict double intersection (sites with both WUE_ET and WUE_T)
- Shows percent of sites responding (Increase vs Decrease vs No change)
- Saves Google Docs compatible JPG version
- UPDATED: Uses existing site-level NN bootstrap CI classifications from CSV
- FIXED: Y-axis limit increased to 50, text labels no longer overlap horizontal lines
- ADJUSTED: Figure dimensions increased for better proportion
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch

# ============================================================================
# SETUP: UPDATED FIGURE SIZE FOR BETTER PROPORTIONS
# ============================================================================
plt.rcParams['figure.figsize'] = [52, 58]  # Increased from [45, 50] for better proportions
plt.rcParams['font.size'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 84
plt.rcParams['axes.labelsize'] = 78
plt.rcParams['xtick.labelsize'] = 72
plt.rcParams['ytick.labelsize'] = 72
plt.rcParams['legend.fontsize'] = 80

# ============================================================================
# FILE PATHS
# ============================================================================
SITE_LEVEL_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"
output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\figures_SPEI6_SPEI48"
os.makedirs(output_folder, exist_ok=True)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def compute_percent_responses(df, pass_filter, conditions, timescales, metrics):
    """
    Compute percent responses for a given pass filter.
    Uses existing Direction column from site-level CSV (NN bootstrap CI classifications).
    Returns: (summary_df, site_counts_dict)
    """
    print(f"\n{'='*60}")
    print(f"COMPUTING FOR {pass_filter}")
    print(f"{'='*60}")
    print("NOTE: Direction classification uses updated site-level NN bootstrap CI classifications from SPEI_site_level_details_FINAL.csv.")
    
    # Initial filter
    filtered = df[
        (df['Pass'] == pass_filter) &
        (df['Salinity'] == "All") &
        (df['SPEI_Timescale'].isin(timescales)) &
        (df['Condition'].isin(conditions)) &
        (df['WUE_Metric'].isin(metrics))
    ].copy()
    
    print(f"After filtering: {len(filtered):,} records, {filtered['Site'].nunique()} unique sites")
    
    # Apply strict double intersection per combination
    shared_sites_dict = {}
    for timescale in timescales:
        for condition in conditions:
            subset = filtered[
                (filtered['SPEI_Timescale'] == timescale) &
                (filtered['Condition'] == condition)
            ]
            
            wue_sites = set(subset[subset['WUE_Metric'] == 'WUE']['Site'].unique())
            tra_sites = set(subset[subset['WUE_Metric'] == 'WUE_tra']['Site'].unique())
            shared_sites = wue_sites.intersection(tra_sites)
            key = f"{timescale}_{condition}"
            shared_sites_dict[key] = shared_sites
            print(f"  {timescale} | {condition}: {len(shared_sites)} shared sites")
    
    # Filter to shared sites only
    filtered_shared = filtered[
        filtered.apply(lambda row: 
            row['Site'] in shared_sites_dict.get(f"{row['SPEI_Timescale']}_{row['Condition']}", set()), 
            axis=1)
    ].copy()
    
    # Compute summary statistics using existing Direction column
    summary_list = []
    for timescale in timescales:
        for condition in conditions:
            for metric in metrics:
                subset = filtered_shared[
                    (filtered_shared['SPEI_Timescale'] == timescale) &
                    (filtered_shared['Condition'] == condition) &
                    (filtered_shared['WUE_Metric'] == metric)
                ]
                
                if len(subset) > 0:
                    total = len(subset)
                    increase = len(subset[subset['Direction'] == 'Increase'])
                    decrease = len(subset[subset['Direction'] == 'Decrease'])
                    no_change = len(subset[subset['Direction'] == 'No change'])
                    insufficient = len(subset[subset['Direction'] == 'Insufficient data'])
                    
                    # Use ONLY valid classified sites in denominator
                    valid_total = increase + decrease + no_change
                    
                    if valid_total > 0:
                        pct_inc = 100 * increase / valid_total
                        pct_dec = 100 * decrease / valid_total
                        pct_nc = 100 * no_change / valid_total
                    else:
                        pct_inc = pct_dec = pct_nc = 0
                    
                    summary_list.append({
                        'SPEI_Timescale': timescale,
                        'Condition': condition,
                        'WUE_Metric': metric,
                        'Total_sites': total,
                        'Valid_sites': valid_total,
                        'Insufficient': insufficient,
                        'Increase': increase,
                        'Decrease': decrease,
                        'No_change': no_change,
                        'Pct_Increase': pct_inc,
                        'Pct_Decrease': pct_dec,
                        'Pct_No_Change': pct_nc
                    })
                    
                    print(f"  {timescale} | {condition} | {metric}: "
                          f"Inc={increase} ({pct_inc:.0f}%), Dec={decrease} ({pct_dec:.0f}%), "
                          f"NC={no_change} ({pct_nc:.0f}%), Insuff={insufficient}")
    
    site_counts = {k: len(v) for k, v in shared_sites_dict.items()}
    
    return pd.DataFrame(summary_list), site_counts

# ============================================================================
# CREATE 4-PANEL FIGURE
# ============================================================================

def create_4panel_figure(df_all, site_counts_all, df_severe, site_counts_severe, output_path):
    """
    Create 4-panel figure with updated y-axis limit (50) and adjusted text positioning
    """
    
    print("\n" + "=" * 80)
    print("CREATING 4-PANEL FIGURE")
    print("=" * 80)
    
    # UPDATED: Y-axis limit increased to 50 to prevent label overlap
    Y_AXIS_MAX = 50
    Y_AXIS_TICKS = np.arange(0, 51, 10)
    bar_width = 0.45
    
    colors = {
        'Efficient': '#66BB6A',      # Green
        'Less Efficient': '#CCCCCC'  # Grey
    }
    
    metric_order = ["WUE$_{ET}$", "WUE$_T$"]
    
    # Create figure with 2 rows, 2 columns - increased size
    fig, axes = plt.subplots(2, 2, figsize=(52, 58))
    
    # Define panels configuration
    panels = [
        (axes[0, 0], df_all, site_counts_all, ["Dry (all)", "Wet (all)"], 
         "SPEI-6 (Short-term)", "Dry", "Wet", "(a)", "SPEI_6"),
        
        (axes[0, 1], df_all, site_counts_all, ["Dry (all)", "Wet (all)"], 
         "SPEI-48 (Long-term)", "Dry", "Wet", "(b)", "SPEI_48"),
        
        (axes[1, 0], df_severe, site_counts_severe, ["Severe Dry", "Severe Wet"], 
         "SPEI-6 (Short-term)", "Ext. Dry", "Ext. Wet", "(c)", "SPEI_6"),
        
        (axes[1, 1], df_severe, site_counts_severe, ["Severe Dry", "Severe Wet"], 
         "SPEI-48 (Long-term)", "Ext. Dry", "Ext. Wet", "(d)", "SPEI_48")
    ]
    
    # Function to create a single panel
    def create_panel(ax, data, site_counts, conditions, title, dry_label, wet_label, panel_label, timescale):
        x_base = np.arange(len(metric_order))
        
        for i, condition in enumerate(conditions):
            cond_data = data[data['Condition'] == condition]
            x_positions = x_base + i * bar_width
            
            # Get values for each metric
            eff_values = []
            less_values = []
            
            for metric_display in metric_order:
                if metric_display == "WUE$_T$":
                    orig_metric = "WUE_tra"
                else:
                    orig_metric = "WUE"
                
                metric_row = cond_data[cond_data['WUE_Metric'] == orig_metric]
                if len(metric_row) > 0:
                    row = metric_row.iloc[0]
                    eff_values.append(row['Pct_Increase'])
                    less_values.append(row['Pct_Decrease'])
                else:
                    eff_values.append(0)
                    less_values.append(0)
            
            # Draw bars
            for j in range(len(eff_values)):
                x_pos = x_positions[j]
                eff = eff_values[j]
                less = less_values[j]
                
                # Efficient (Increase) - Green
                ax.bar(x_pos, eff, bar_width, 
                      color=colors['Efficient'], 
                      edgecolor='black', linewidth=3.5,
                      alpha=0.95, zorder=5)
                
                # Less Efficient (Decrease) - Grey
                ax.bar(x_pos, less, bar_width, 
                      bottom=eff,
                      color=colors['Less Efficient'],
                      edgecolor='black', linewidth=3.5,
                      alpha=0.95, zorder=5)
            
            # Add percentage labels with adjusted positioning to avoid overlap
            for j, (eff, less) in enumerate(zip(eff_values, less_values)):
                x_pos = x_positions[j]
                bbox_props = dict(boxstyle="round,pad=0.3", facecolor="white", 
                                edgecolor="none", alpha=0.9)
                
                # Position label for Increase (green) section
                if eff > 2:  # Only label if > 2%
                    # Place in middle of green bar, but avoid conflict with horizontal lines
                    eff_y = eff / 2
                    # Ensure label stays within plot area with margin
                    if eff_y > Y_AXIS_MAX - 3:
                        eff_y = Y_AXIS_MAX - 4
                    ax.text(x_pos, eff_y, f'{eff:.0f}%', 
                           ha='center', va='center', 
                           fontsize=50, fontweight='bold',
                           color='black', bbox=bbox_props, zorder=25)
                
                # Position label for Decrease (grey) section
                if less > 2:  # Only label if > 2%
                    less_y = eff + (less / 2)
                    # Ensure label stays within plot area with margin
                    if less_y > Y_AXIS_MAX - 3:
                        less_y = Y_AXIS_MAX - 4
                    ax.text(x_pos, less_y, f'{less:.0f}%', 
                           ha='center', va='center', 
                           fontsize=50, fontweight='bold',
                           color='black', bbox=bbox_props, zorder=25)
            
            # Add horizontal lines for Wet/Ext. Wet bars (positioned carefully to avoid text)
            if "Wet" in condition:
                for j in range(len(eff_values)):
                    x_pos = x_positions[j]
                    eff = eff_values[j]
                    less = less_values[j]
                    line_length = bar_width * 0.7
                    
                    # Only draw lines if there's enough space and not conflicting with labels
                    if eff > 8:  # Only draw if bar is tall enough
                        line_y1 = eff * 0.3  # Place at 30% height of green bar
                        # Ensure line not too close to top
                        if line_y1 < Y_AXIS_MAX - 5:
                            ax.hlines(line_y1, x_pos - line_length/2, x_pos + line_length/2,
                                     colors='#777777', linewidth=5, zorder=15)
                    
                    if less > 8:  # Only draw if bar is tall enough
                        line_y2 = eff + (less * 0.3)  # Place at 30% height of grey bar
                        if line_y2 < Y_AXIS_MAX - 5:
                            ax.hlines(line_y2, x_pos - line_length/2, x_pos + line_length/2,
                                     colors='#777777', linewidth=5, zorder=15)
        
        # Formatting with updated y-axis
        ax.set_xticks(x_base + bar_width/2)
        ax.set_xticklabels(metric_order, fontsize=78, fontweight='bold')
        ax.set_ylim(0, Y_AXIS_MAX)
        ax.set_ylabel('Percent of sites (%)', fontsize=84, fontweight='bold', labelpad=40)
        ax.grid(True, alpha=0.08, axis='y', linestyle=':', linewidth=2.0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(5.0)
        ax.spines['bottom'].set_linewidth(5.0)
        ax.set_yticks(Y_AXIS_TICKS)
        ax.tick_params(axis='y', which='major', length=25, width=4.0)
        ax.tick_params(axis='x', which='major', length=20, width=4.0)
        ax.set_title(title, fontsize=90, fontweight='bold', pad=28)
        ax.axhline(y=0, color='black', linewidth=3.0, alpha=0.3)
        
        # Add n= labels below axes (adjusted position for new y-limit)
        dry_x = x_base[0] - 0.15
        wet_x = x_base[0] + bar_width + 0.15
        
        dry_key = f"{timescale}_{conditions[0]}"
        wet_key = f"{timescale}_{conditions[1]}"
        
        dry_count = site_counts.get(dry_key, 0)
        wet_count = site_counts.get(wet_key, 0)
        
        # Position N labels below the x-axis
        ax.text(dry_x, -6, f'{dry_label}: N={dry_count}', 
               ha='center', va='top', fontsize=52, fontweight='bold', fontstyle='italic')
        ax.text(wet_x, -6, f'{wet_label}: N={wet_count}', 
               ha='center', va='top', fontsize=52, fontweight='bold', fontstyle='italic')
        
        # Add panel label
        ax.text(-0.18, 1.05, panel_label, transform=ax.transAxes, 
               fontsize=96, fontweight='bold',
               verticalalignment='top', horizontalalignment='left', zorder=100)
    
    # Create all panels
    for ax, data, counts, conditions, title, dry_label, wet_label, panel_label, timescale in panels:
        panel_data = data[data['SPEI_Timescale'] == timescale]
        create_panel(ax, panel_data, counts, conditions, title, dry_label, wet_label, panel_label, timescale)
    
    # Adjust layout with more space
    plt.tight_layout(rect=[0, 0.04, 1, 0.58])
    fig.subplots_adjust(top=0.90, hspace=0.32, wspace=0.28)
    
    # ========== LEGEND ==========
    class HandlerWetMarker(HandlerPatch):
        def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
            y1 = ydescent + height * 0.33
            y2 = ydescent + height * 0.67
            x0 = xdescent + width * 0.1
            x1 = xdescent + width * 0.9
            line1 = plt.Line2D([x0, x1], [y1, y1], color='#777777', lw=6, transform=trans)
            line2 = plt.Line2D([x0, x1], [y2, y2], color='#777777', lw=6, transform=trans)
            return [line1, line2]
    
    wet_marker_handle = mpatches.Rectangle((0, 0), 1, 1, facecolor='none', edgecolor='none')
    
    legend_elements = [
        mpatches.Patch(facecolor=colors['Efficient'], edgecolor='black', linewidth=3.0, alpha=0.9),
        mpatches.Patch(facecolor=colors['Less Efficient'], edgecolor='black', linewidth=3.0, alpha=0.9),
        mpatches.Patch(facecolor='white', edgecolor='black', linewidth=3.0, alpha=0.9),
        wet_marker_handle
    ]
    
    legend_labels = ['Efficient', 'Less Efficient', 'No change', 'Wet / Ext. Wet']
    
    legend = fig.legend(
        handles=legend_elements, labels=legend_labels, loc='upper center',
        bbox_to_anchor=(0.5, 1.06), borderaxespad=0.4, ncol=4, fontsize=60,
        frameon=True, framealpha=0.95, edgecolor='black', borderpad=1.5,
        handlelength=2.0, handleheight=1.0, handletextpad=1.0,
        columnspacing=2.0, labelspacing=1.0,
        handler_map={wet_marker_handle: HandlerWetMarker()}
    )
    legend.get_frame().set_facecolor('#f8f9fa')
    legend.get_frame().set_linewidth(2.0)
    legend.get_frame().set_boxstyle('round', pad=0.4, rounding_size=0.3)
    
    # Ensure text is in front
    for ax_row in axes:
        for ax in ax_row:
            for text in ax.texts:
                text.set_zorder(100)
            ax.title.set_zorder(100)
            ax.xaxis.label.set_zorder(100)
            ax.yaxis.label.set_zorder(100)
    legend.set_zorder(100)
    
    # Save figure
    print("\n" + "=" * 80)
    print("SAVING FIGURES")
    print("=" * 80)
    
    # Save high-res PNG
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ High-res PNG: {output_path}")
    
    # Save JPG for Google Docs
    jpg_path = output_path.replace('.png', '.jpg')
    plt.savefig(jpg_path, dpi=150, bbox_inches='tight', facecolor='white', format='jpg')
    print(f"✓ Google Docs JPG: {jpg_path}")
    
    # Save PDF for Word
    pdf_path = output_path.replace('.png', '.pdf')
    plt.savefig(pdf_path, dpi=150, bbox_inches='tight', facecolor='white', format='pdf')
    print(f"✓ Word PDF: {pdf_path}")
    
    plt.show()
    return fig

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("=" * 100)
    print("4-PANEL BAR CHART FIGURE - UPDATED VERSION")
    print("Top row (a,b): All conditions (PASS C)")
    print("Bottom row (c,d): Severe/Extreme conditions (PASS B)")
    print("=" * 100)
    print("\n✓ Y-axis limit increased to 50 for better label spacing")
    print("✓ Figure dimensions increased (52x58) for improved proportions")
    print("✓ Text labels repositioned to avoid overlap with horizontal lines")
    print("✓ Using updated NN bootstrap CI classifications from site-level CSV")
    
    # Load data
    print("\n1. LOADING SITE-LEVEL DATA...")
    df = pd.read_csv(SITE_LEVEL_FILE)
    print(f"   Loaded {len(df):,} rows")
    print(f"   Direction column unique values: {df['Direction'].unique()}")
    
    timescales = ['SPEI_6', 'SPEI_48']
    metrics = ['WUE', 'WUE_tra']
    
    # Compute Panel A (All conditions - PASS C)
    print("\n2. COMPUTING ALL CONDITIONS (PASS C)...")
    df_all, site_counts_all = compute_percent_responses(
        df, 'PASS C', ['Dry (all)', 'Wet (all)'], timescales, metrics)
    
    # Compute Panel B (Severe/Extreme - PASS B)
    print("\n3. COMPUTING SEVERE/EXTREME CONDITIONS (PASS B)...")
    df_severe, site_counts_severe = compute_percent_responses(
        df, 'PASS B', ['Severe Dry', 'Severe Wet'], timescales, metrics)
    
    # Print results
    print("\n" + "=" * 100)
    print("RESULTS SUMMARY")
    print("=" * 100)
    
    print("\nPANEL A (All conditions - PASS C):")
    for timescale in timescales:
        print(f"\n  {timescale.replace('_', '-')}:")
        for condition in ['Dry (all)', 'Wet (all)']:
            for metric in ['WUE', 'WUE_tra']:
                row = df_all[
                    (df_all['SPEI_Timescale'] == timescale) &
                    (df_all['Condition'] == condition) &
                    (df_all['WUE_Metric'] == metric)
                ]
                if len(row) > 0:
                    row = row.iloc[0]
                    metric_label = 'WUE_ET' if metric == 'WUE' else 'WUE_T'
                    print(f"    {condition} | {metric_label}: "
                          f"Inc={row['Pct_Increase']:.0f}%, Dec={row['Pct_Decrease']:.0f}%, NC={row['Pct_No_Change']:.0f}% "
                          f"(n={row['Valid_sites']} valid, {row['Insufficient']} insufficient)")
    
    print("\nPANEL B (Severe/Extreme - PASS B):")
    for timescale in timescales:
        print(f"\n  {timescale.replace('_', '-')}:")
        for condition in ['Severe Dry', 'Severe Wet']:
            for metric in ['WUE', 'WUE_tra']:
                row = df_severe[
                    (df_severe['SPEI_Timescale'] == timescale) &
                    (df_severe['Condition'] == condition) &
                    (df_severe['WUE_Metric'] == metric)
                ]
                if len(row) > 0:
                    row = row.iloc[0]
                    metric_label = 'WUE_ET' if metric == 'WUE' else 'WUE_T'
                    print(f"    {condition} | {metric_label}: "
                          f"Inc={row['Pct_Increase']:.0f}%, Dec={row['Pct_Decrease']:.0f}%, NC={row['Pct_No_Change']:.0f}% "
                          f"(n={row['Valid_sites']} valid, {row['Insufficient']} insufficient)")
    
    # Create 4-panel figure
    print("\n4. CREATING 4-PANEL FIGURE...")
    output_path = os.path.join(output_folder, "Combined_4Panel_BarChart.png")
    fig = create_4panel_figure(df_all, site_counts_all, df_severe, site_counts_severe, output_path)
    
    # Save data
    df_all.to_csv(os.path.join(output_folder, 'Panel_AllConditions_Data.csv'), index=False)
    df_severe.to_csv(os.path.join(output_folder, 'Panel_SevereConditions_Data.csv'), index=False)
    
    print("\n" + "=" * 100)
    print("WORKFLOW COMPLETE ✓")
    print("=" * 100)
    print("\n📊 OUTPUT FILES:")
    print(f"   Google Docs: {output_path.replace('.png', '.jpg')}")
    print(f"   Microsoft Word: {output_path.replace('.png', '.pdf')}")
    print(f"   High-res backup: {output_path}")
    print(f"   Data (All conditions): {os.path.join(output_folder, 'Panel_AllConditions_Data.csv')}")
    print(f"   Data (Severe conditions): {os.path.join(output_folder, 'Panel_SevereConditions_Data.csv')}")
    print("\n✅ Improvements applied:")
    print("   • Y-axis limit increased from 40 to 50")
    print("   • Figure size increased from 45x50 to 52x58")
    print("   • Text labels repositioned to avoid overlap")
    print("   • Horizontal lines repositioned (30% bar height)")
    print("   • Using NN bootstrap CI classifications from upstream CSV")

if __name__ == "__main__":
    main()