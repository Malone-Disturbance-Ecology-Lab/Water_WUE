# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 20:01:50 2025

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 20:22:57 2025

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
FIGURE 2 WORKFLOW - Percent Responses for Freshwater vs Saline (PASS C)

Creates grouped bar charts showing Efficient/Less Efficient responses
for Freshwater vs Saline sites across SPEI-6 and SPEI-48.

Data filtering:
- Pass = "PASS C"
- Salinity in ["Freshwater", "Saline"]
- SPEI_Timescale in ["SPEI_6", "SPEI_48"]
- Condition in ["Dry (all)", "Wet (all)"]
- WUE_Metric in ["WUE", "WUE_eva", "WUE_tra"]

Plot structure:
- Left panel: SPEI-6
- Right panel: SPEI-48
- Each panel: Four grouped bars per metric (Dry-FW, Dry-Saline, Wet-FW, Wet-Saline)
- Colors: Green (Efficient), Grey (Less Efficient)
- Salinity encoding: Blue top (Freshwater), Orange top (Saline)
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import matplotlib
import os

# =============================================================================
# STYLING - CLEAN AND READABLE
# =============================================================================
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 18,
    'axes.titlesize': 22,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 16,
    'figure.constrained_layout.use': True,
})

# =============================================================================
# 1. DATA LOADING AND FILTERING
# =============================================================================
def load_figure2_data(filepath):
    """Load and filter data for Figure 2."""
    df = pd.read_csv(filepath)
    
    # Apply filters - SAME AS HEATMAP BUT WITH SALINITY SPLIT
    filters = (
        (df['Pass'] == 'PASS C') &
        (df['Salinity'].isin(['Freshwater', 'Saline'])) &
        (df['SPEI_Timescale'].isin(['SPEI_6', 'SPEI_48'])) &
        (df['Condition'].isin(['Dry (all)', 'Wet (all)'])) &
        (df['WUE_Metric'].isin(['WUE', 'WUE_eva', 'WUE_tra']))
    )
    
    filtered_df = df[filters].copy()
    
    # Remove Brackish if present
    filtered_df = filtered_df[filtered_df['Salinity'] != 'Brackish']
    
    # Ensure consistent ordering
    for col, categories in [
        ('Condition', ['Dry (all)', 'Wet (all)']),
        ('SPEI_Timescale', ['SPEI_6', 'SPEI_48']),
        ('WUE_Metric', ['WUE', 'WUE_eva', 'WUE_tra']),
        ('Salinity', ['Freshwater', 'Saline'])
    ]:
        filtered_df[col] = pd.Categorical(filtered_df[col], categories=categories, ordered=True)
    
    # Sort for consistent plotting
    filtered_df = filtered_df.sort_values(['SPEI_Timescale', 'Condition', 'Salinity', 'WUE_Metric'])
    
    return filtered_df

# =============================================================================
# 2. CALCULATE PERCENTAGES
# =============================================================================
def calculate_response_percentages(df):
    """Calculate Efficient and Less Efficient percentages for each group."""
    results = []
    
    # Group by all relevant dimensions
    for (timescale, condition, salinity, metric), group in df.groupby(
        ['SPEI_Timescale', 'Condition', 'Salinity', 'WUE_Metric']
    ):
        if len(group) == 0:
            continue
            
        row = group.iloc[0]
        total_sites = row['Total_sites']
        
        # Calculate percentages (same as original Figure 2)
        pct_efficient = (row['Increase'] / total_sites * 100) if total_sites > 0 else 0
        pct_less_efficient = (row['Decrease'] / total_sites * 100) if total_sites > 0 else 0
        pct_no_change = (row['No_change'] / total_sites * 100) if total_sites > 0 else 0
        
        # Clean up condition name
        condition_clean = condition.replace(' (all)', '')
        
        results.append({
            'SPEI_Timescale': timescale,
            'Condition': condition_clean,
            'Salinity': salinity,
            'WUE_Metric': metric,
            'Total_sites': total_sites,
            'Pct_Efficient': pct_efficient,
            'Pct_Less_Efficient': pct_less_efficient,
            'Pct_No_Change': pct_no_change
        })
    
    return pd.DataFrame(results)

# =============================================================================
# 3. CREATE PLOT
# =============================================================================

def create_figure2_plot(df_percent, save_path):
    """Create Figure 2 grouped bar charts."""
    if len(df_percent) == 0:
        print("ERROR: No data for Figure 2!")
        return None
    
    # Set up figure with increased size
    fig, axes = plt.subplots(1, 2, figsize=(55, 28), constrained_layout=True)
    
    panel_fontsize = 96
    # Left panel (b)
    axes[0].text(-0.12, 1.05, '(b)', transform=axes[0].transAxes,
            fontsize=panel_fontsize, fontweight='bold',
            verticalalignment='top', horizontalalignment='left',
            zorder=100)

    # Right panel (c)
    axes[1].text(-0.12, 1.05, '(c)', transform=axes[1].transAxes,
            fontsize=panel_fontsize, fontweight='bold',
            verticalalignment='top', horizontalalignment='left',
            zorder=100) 
    
    # Define colors and styles
    colors = {
        'Efficient': '#27ae60',
        'Less_Efficient': '#CCCCCC',
    }
    
    salinity_colors = {
        'Freshwater': '#0000FF',
        'Saline': '#FFA500'
    }
    
    bar_width = 0.45
    group_spacing = 0.60
    cap_height = 6
    
    n_values = {}
    timescales = ['SPEI_6', 'SPEI_48']
    timescale_labels = ['SPEI-6', 'SPEI-48']
    
    for ax_idx, (timescale, timescale_label) in enumerate(zip(timescales, timescale_labels)):
        ax = axes[ax_idx]
        ts_data = df_percent[df_percent['SPEI_Timescale'] == timescale]
        
        metrics = ['WUE', 'WUE_eva', 'WUE_tra']
        metric_labels = ['WUE$_{ET}$', 'WUE$_E$', 'WUE$_T$']
        
        n_metrics = len(metrics)
        x_positions = np.arange(n_metrics) * (4 * bar_width + group_spacing)
        
        for metric_idx, (metric, metric_label) in enumerate(zip(metrics, metric_labels)):
            metric_data = ts_data[ts_data['WUE_Metric'] == metric]
            base_x = x_positions[metric_idx]
            
            conditions_order = [
                ('Dry', 'Freshwater'),
                ('Dry', 'Saline'),
                ('Wet', 'Freshwater'),
                ('Wet', 'Saline')
            ]
            
            for bar_idx, (condition, salinity) in enumerate(conditions_order):
                match = metric_data[
                    (metric_data['Condition'] == condition) &
                    (metric_data['Salinity'] == salinity)
                ]
                
                if len(match) == 0:
                    continue
                
                row = match.iloc[0]
                x_pos = base_x + bar_idx * bar_width
                total_height = row['Pct_Efficient'] + row['Pct_Less_Efficient']
                
                n_key = (timescale, condition, salinity)
                if n_key not in n_values:
                    n_values[n_key] = row['Total_sites']
                
                cap_color = salinity_colors[salinity]
                
                # Plot Efficient (green) bars
                ax.bar(x_pos, row['Pct_Efficient'], 
                       width=bar_width, color=colors['Efficient'],
                       edgecolor='black', linewidth=4.0, zorder=2)
                
                # Plot Less Efficient (grey) bars
                ax.bar(x_pos, row['Pct_Less_Efficient'], 
                       width=bar_width, color=colors['Less_Efficient'],
                       edgecolor='black', linewidth=4.0,
                       bottom=row['Pct_Efficient'], zorder=2)
                
                # Colored cap at top
                cap_y = total_height
                ax.add_patch(Rectangle(
                    (x_pos - bar_width/2, cap_y - 2.5),
                    bar_width, cap_height,
                    facecolor=cap_color, edgecolor='none', linewidth=0, zorder=3))
                
                # Percentage labels
                if row['Pct_Efficient'] >= 5:
                    label_y = row['Pct_Efficient'] / 2
                    ax.text(x_pos, label_y, f"{row['Pct_Efficient']:.0f}%",
                           ha='center', va='center', fontsize=44, fontweight='bold',
                           bbox=dict(boxstyle="round,pad=0.2", facecolor='white',
                                    edgecolor='none', alpha=0.9, linewidth=0), zorder=4)
                
                if row['Pct_Less_Efficient'] >= 5:
                    label_y = row['Pct_Efficient'] + (row['Pct_Less_Efficient'] / 2)
                    ax.text(x_pos, label_y, f"{row['Pct_Less_Efficient']:.0f}%",
                           ha='center', va='center', fontsize=44, fontweight='bold',
                           bbox=dict(boxstyle="round,pad=0.2", facecolor='white',
                                    edgecolor='none', alpha=0.9, linewidth=0), zorder=4)
                
                # Dry/Wet labels
                if bar_idx == 0:
                    ax.text(x_pos, -2.5, 'Dry', ha='center', va='top', fontsize=36,
                           fontweight='bold', bbox=dict(boxstyle="round,pad=0.3",
                           facecolor='white', edgecolor='gray', alpha=0.95, linewidth=1.5), zorder=4)
                elif bar_idx == 2:
                    ax.text(x_pos, -2.5, 'Wet', ha='center', va='top', fontsize=36,
                           fontweight='bold', bbox=dict(boxstyle="round,pad=0.3",
                           facecolor='white', edgecolor='gray', alpha=0.95, linewidth=1.5), zorder=4)
        
        ax.set_xticks(x_positions)
        ax.set_xticklabels(metric_labels, fontweight='bold', fontsize=78)
        ax.set_ylabel('Percentage of Sites (%)', fontweight='bold', fontsize=70, labelpad=40)
        ax.set_ylim(-10, 110)
        ax.set_yticks(range(0, 101, 20))
        ax.tick_params(axis='y', labelsize=66)
        ax.tick_params(axis='x', labelsize=66)
        ax.grid(True, axis='y', alpha=0.15, linestyle=':', linewidth=3.0, zorder=0)
        
        if ax_idx == 0:
            ax.set_title(f'{timescale_label}', fontweight='bold', fontsize=90, 
                        pad=40, loc='left', x=0.02)
        else:
            ax.set_title(f'{timescale_label}', fontweight='bold', fontsize=90, 
                        pad=40, loc='right', x=0.98)
        
        for spine in ax.spines.values():
            spine.set_linewidth(5)
    
    # Add n-values
    for ax_idx, (timescale, timescale_label) in enumerate(zip(timescales, timescale_labels)):
        ax = axes[ax_idx]
        first_metric_x = 0
        conditions_order = [
            ('Dry', 'Freshwater'), ('Dry', 'Saline'),
            ('Wet', 'Freshwater'), ('Wet', 'Saline')
        ]
        for bar_idx, (condition, salinity) in enumerate(conditions_order):
            n_key = (timescale, condition, salinity)
            if n_key in n_values:
                n_val = int(n_values[n_key])
                x_pos = first_metric_x + bar_idx * bar_width
                ax.text(x_pos, -6, f"n={n_val}", ha='center', va='top',
                       fontsize=38, fontweight='bold', zorder=4)
    
    # Legend
    legend_elements = [
        Patch(facecolor=colors['Efficient'], edgecolor='black', linewidth=4, label='Efficient'),
        Patch(facecolor=colors['Less_Efficient'], edgecolor='black', linewidth=4, label='Less Efficient'),
        Patch(facecolor=salinity_colors['Freshwater'], edgecolor='black', linewidth=4, label='Freshwater'),
        Patch(facecolor=salinity_colors['Saline'], edgecolor='black', linewidth=4, label='Saline'),
    ]
    
    legend = fig.legend(handles=legend_elements, loc='upper center',
                       bbox_to_anchor=(0.5, 1.02), ncol=4, fontsize=60,
                       frameon=True, framealpha=0.95, borderpad=1.2,
                       handlelength=2.5, handleheight=1.5, handletextpad=0.8,
                       columnspacing=1.0, labelspacing=0.5)
    legend.get_frame().set_linewidth(4.0)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_facecolor('#f8f9fa')
    
    fig.subplots_adjust(top=0.86)
    
    # Save
    png_path = f"{save_path}.png"
    pdf_path = f"{save_path}.pdf"
    fig.savefig(png_path, dpi=500, bbox_inches='tight', facecolor='white')
    print(f"   Google Docs PNG: {png_path}")
    
    plt.show()
    return fig

# THEN DELETE ALL OTHER def create_figure2_plot definitions
# Keep print_figure2_results and main unchanged
# =============================================================================
# 4. PRINT RESULTS SUMMARY
# =============================================================================
def print_figure2_results(df_percent):
    """Print formatted results for manuscript."""
    print("\n" + "="*120)
    print("FIGURE 2 RESULTS - Freshwater vs Saline Percent Responses (PASS C)")
    print("="*120)
    
    if len(df_percent) == 0:
        return
    
    # Group by timescale and metric
    for timescale in ['SPEI_6', 'SPEI_48']:
        ts_display = timescale.replace('_', '-')
        print(f"\n{ts_display}:")
        print("-"*60)
        
        ts_data = df_percent[df_percent['SPEI_Timescale'] == timescale]
        
        for metric, metric_name in [('WUE', 'WUE'),
                                   ('WUE_eva', 'WUE_E'),
                                   ('WUE_tra', 'WUE_T')]:
            print(f"\n  {metric_name}:")
            
            metric_data = ts_data[ts_data['WUE_Metric'] == metric]
            
            for condition in ['Dry', 'Wet']:
                for salinity in ['Freshwater', 'Saline']:
                    match = metric_data[
                        (metric_data['Condition'] == condition) &
                        (metric_data['Salinity'] == salinity)
                    ]
                    
                    if len(match) > 0:
                        row = match.iloc[0]
                        print(f"    {condition}-{salinity:12s} (n={int(row['Total_sites']):3d}): "
                              f"Efficient={row['Pct_Efficient']:5.1f}%, "
                              f"Less Efficient={row['Pct_Less_Efficient']:5.1f}%")

# =============================================================================
# 5. MAIN WORKFLOW
# =============================================================================
def main():
    """Main execution function."""
    # File paths - SAME AS FIGURE 4
    input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_analysis_summary_FINAL.csv"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\figures_SPEI6_SPEI48"
    output_base = os.path.join(output_dir, 'Fig2_FWvsSaline_PercentResponses_PASSC_SPEI6_SPEI48')
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*120)
    print("FIGURE 2: FRESHWATER VS SALINE PERCENT RESPONSES")
    print("="*120)
    print("DATA FILTERS:")
    print("  - Pass = PASS C")
    print("  - Salinity = Freshwater, Saline")
    print("  - SPEI Timescale = SPEI-6, SPEI-48")
    print("  - Condition = Dry (all), Wet (all)")
    print("  - WUE Metrics = WUE, WUE_E, WUE_T")
    print("\nPLOT STRUCTURE:")
    print("  - Left panel: SPEI-6")
    print("  - Right panel: SPEI-48")
    print("  - Each panel: 4 bars per metric (Dry-FW, Dry-Saline, Wet-FW, Wet-Saline)")
    print("  - Colors: Green (Efficient), Grey (Less Efficient)")
    print("  - Salinity encoding: Blue top (Freshwater), Orange top (Saline)")
    print("  - Percentage labels on bars (>5%)")
    print("  - n-values shown once below SPEI-6 panel")
    print("="*120)
    
    # Load and filter data
    print("\nLoading data...")
    data = load_figure2_data(input_file)
    
    if len(data) == 0:
        print("ERROR: No data found! Check file path and filters.")
        return
    
    print(f"SUCCESS: Loaded {len(data)} data points")
    print(f"   Conditions: {data['Condition'].unique().tolist()}")
    print(f"   Salinities: {data['Salinity'].unique().tolist()}")
    print(f"   Timescales: {data['SPEI_Timescale'].unique().tolist()}")
    print(f"   Metrics: {data['WUE_Metric'].unique().tolist()}")
    
    # Calculate percentages
    print("\nCalculating response percentages...")
    df_percent = calculate_response_percentages(data)
    print(f"Calculated percentages for {len(df_percent)} groups")
    
    # Create plot
    print("\nCreating plot...")
    fig = create_figure2_plot(df_percent, output_base)
    
    # Print results summary
    print("\n" + "="*120)
    print("RESULTS SUMMARY")
    print("="*120)
    print_figure2_results(df_percent)
    
    print("\n" + "="*120)
    print("WORKFLOW COMPLETE")
    print("="*120)
    
    return fig

# =============================================================================
# EXECUTION
# =============================================================================
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    