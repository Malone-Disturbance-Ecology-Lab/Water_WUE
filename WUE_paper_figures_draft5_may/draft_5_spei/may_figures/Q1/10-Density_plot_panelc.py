# -*- coding: utf-8 -*-
"""
Created for Panel C: Distribution of site sample sizes (number of NN observations)
SINGLE CONTINUOUS PANEL with three ecosystem distributions side-by-side
ORDER: Upland (left) | Freshwater (middle) | Saline (right)
Colors: EXACT same as Figure 2
SHOWS: Histogram (site counts) + thin KDE overlay scaled to counts
FIXED: Now counts ALL retained NN observations (not unique month combos)
FIXED: Removed duplicate function blocks
FIXED: Updated terminology from "months" to "observations"
FIXED: Added strict metric filter (WUE, WUE_eva, WUE_tra ALL must have data)
@author: ammar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import gaussian_kde

# =============================================================================
# CONFIGURATION - EXACT SAME AS PANEL B (FIGURE 2)
# =============================================================================

# Figure dimensions (SAME as Figure 2)
FIGURE_WIDTH = 11
FIGURE_HEIGHT = 5

# Font sizes (SAME as Figure 2)
X_TICK_LABEL_SIZE = 20
Y_TICK_LABEL_SIZE = 24
X_AXIS_LABEL_SIZE = 22
Y_AXIS_LABEL_SIZE = 22
ECOSYSTEM_LABEL_SIZE = 24
ANNOTATION_FONT_SIZE = 14

# Histogram styling
HISTOGRAM_ALPHA = 0.5
HISTOGRAM_EDGE_WIDTH = 1

# KDE overlay styling (thin, semi-transparent)
KDE_LINE_WIDTH = 1.5
KDE_ALPHA = 0.6

# Spacing between ecosystem sections
SECTION_SPACING = 0.15

# Colors for ECOSYSTEM TYPES - EXACT SAME RGB as Figure 2
SALINITY_COLORS = {
    'Upland': '#800080',     # Purple - matches Figure 2
    'Freshwater': '#0000FF',  # Blue - matches Figure 2  
    'Saline': '#FFA500'       # Orange - matches Figure 2
}

# Order for display (left to right)
DISPLAY_ORDER = ['Upland', 'Freshwater', 'Saline']

# File paths
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\figures"

# Input data files
BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
MONTHLY_DATA_FILE = os.path.join(BASE_DIR, "monthly_data_after_outlier_removal.csv")
SUMMARY_DATA_FILE = os.path.join(BASE_DIR, "wue_site_level_summary_SPEI_1.csv")

# Output settings
DPI = 600
OUTPUT_FORMAT = 'both'

# Near-normal SPEI range (SAME as Figure 2)
SPEI_MIN = -0.5
SPEI_MAX = 0.5

# =============================================================================
# DATA LOADING - EXACT SAME FILTERING AS FIGURE 2 (FIXED WITH STRICT FILTER)
# =============================================================================

def load_site_sample_sizes():
    """Load monthly data and count retained NN observations per site for each salinity category"""
    print(f"\n[Data Loading] Loading monthly data from {MONTHLY_DATA_FILE}")
    
    if not os.path.exists(MONTHLY_DATA_FILE):
        print(f"❌ ERROR: Monthly data file not found!")
        raise FileNotFoundError(f"Monthly data not found: {MONTHLY_DATA_FILE}")
    
    df_monthly = pd.read_csv(MONTHLY_DATA_FILE)
    print(f"  Loaded {len(df_monthly):,} rows from monthly data")
    
    # Filter to NN conditions - use SPEI_1_Cat if available
    if "SPEI_1_Cat" in df_monthly.columns:
        df_nn = df_monthly[df_monthly["SPEI_1_Cat"] == "NN"].copy()
    else:
        nn_mask = (df_monthly["SPEI_1"] >= SPEI_MIN) & (df_monthly["SPEI_1"] <= SPEI_MAX)
        df_nn = df_monthly[nn_mask].copy()
    print(f"  After NN filter: {len(df_nn):,} rows")
    
    # =========================================================================
    # CRITICAL FIX: APPLY STRICT METRIC FILTER (SAME AS CHUNK 3)
    # Only keep rows where WUE, WUE_eva, AND WUE_tra ALL have data
    # =========================================================================
    print("\n  Applying strict metric filter (WUE, WUE_eva, WUE_tra ALL must have data)...")
    strict_mask = (
        df_nn['WUE'].notna() &
        df_nn['WUE_eva'].notna() &
        df_nn['WUE_tra'].notna()
    )
    df_nn = df_nn[strict_mask].copy()
    print(f"  After strict metric filter: {len(df_nn):,} rows")
    
    # Get strict triple intersection sites
    print(f"\n[Data Loading] Loading summary data from {SUMMARY_DATA_FILE}")
    df_summary = pd.read_csv(SUMMARY_DATA_FILE)
    nn_summary = df_summary[df_summary['SPEI_Class'] == 'NN'].copy()
    
    sites_with_wue = set(nn_summary[nn_summary['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
    sites_with_eva = set(nn_summary[nn_summary['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
    sites_with_tra = set(nn_summary[nn_summary['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
    
    shared_sites = sites_with_wue.intersection(sites_with_eva).intersection(sites_with_tra)
    print(f"  Strict triple intersection sites: {len(shared_sites)}")
    
    # Filter to triple intersection
    df_nn = df_nn[df_nn['site_name'].isin(shared_sites)].copy()
    print(f"  After site filter: {len(df_nn):,} rows")
    print(f"  Unique sites: {df_nn['site_name'].nunique()}")
    
    # Count retained NN OBSERVATIONS per site per salinity category
    salinity_dict = {}
    
    for salinity in DISPLAY_ORDER:
        mask = (df_nn['Salinity_Category'] == salinity)
        site_data = df_nn[mask]  # All rows for this salinity category
        
        # Count ALL retained observations (simple row count, not unique month combos)
        obs_counts = site_data.groupby('site_name').size()
        values = obs_counts.values
        
        salinity_dict[salinity] = values
        
        if len(values) > 0:
            total_obs = np.sum(values)
            print(f"\n  {salinity}: n_sites={len(values)}, total_obs={total_obs}, range={np.min(values)}-{np.max(values)}, median={np.median(values):.0f}")
    
    return salinity_dict, shared_sites

# =============================================================================
# CREATE SINGLE CONTINUOUS PANEL WITH HISTOGRAMS
# =============================================================================

def create_continuous_panel(sample_size_data, shared_sites):
    """Create SINGLE continuous panel with histograms + thin KDE overlay"""
    
    print(f"\nCreating PANEL C: Histograms + KDE overlay")
    print(f"  • Order: Upland (left) | Freshwater (middle) | Saline (right)")
    print(f"  • Colors: EXACT same as Figure 2")
    print(f"  • Y-axis: Number of sites (counts)")
    
    # Create figure with one axis
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    
    # x-axis range for each section: 0-60 observations (to capture Saline max=60)
    ACTUAL_MIN = 0
    ACTUAL_MAX = 60
    SECTION_WIDTH = 40
    
    # Define section boundaries with gaps
    GAP_WIDTH = SECTION_WIDTH * SECTION_SPACING
    
    section_x_starts = []
    section_x_ends = []
    
    # Section 1: Upland (left)
    section1_start = 0
    section1_end = section1_start + SECTION_WIDTH
    section_x_starts.append(section1_start)
    section_x_ends.append(section1_end)
    
    # Section 2: Freshwater (middle)
    section2_start = section1_end + GAP_WIDTH
    section2_end = section2_start + SECTION_WIDTH
    section_x_starts.append(section2_start)
    section_x_ends.append(section2_end)
    
    # Section 3: Saline (right)
    section3_start = section2_end + GAP_WIDTH
    section3_end = section3_start + SECTION_WIDTH
    section_x_starts.append(section3_start)
    section_x_ends.append(section3_end)
    
    # Number of bins for histogram
    N_BINS = 12
    
    # Calculate global maximum count for y-axis scaling
    max_count = 0
    
    for salinity in DISPLAY_ORDER:
        vals = sample_size_data[salinity]
        if len(vals) > 0:
            vals_filtered = vals[(vals >= ACTUAL_MIN) & (vals <= ACTUAL_MAX)]
            if len(vals_filtered) > 0:
                counts, _ = np.histogram(vals_filtered, bins=N_BINS, range=(ACTUAL_MIN, ACTUAL_MAX))
                max_count = max(max_count, np.max(counts))
    
    max_count = max_count * 1.2 if max_count > 0 else 5
    
    # X-axis tick positions for each section (0-60 observations range)
    X_TICK_POSITIONS = [0, 10, 20, 30, 40]
    
    # Plot each ecosystem distribution
    section_centers = []
    
    for idx, salinity in enumerate(DISPLAY_ORDER):
        vals = sample_size_data[salinity]
        x_start = section_x_starts[idx]
        x_end = section_x_ends[idx]
        section_center = (x_start + x_end) / 2
        section_centers.append(section_center)
        color = SALINITY_COLORS[salinity]
        
        # Add vertical separator lines between sections
        if idx > 0:
            ax.axvline(x=x_start - 2, color='gray', linestyle='--', alpha=0.4, linewidth=2)
        
        if len(vals) > 0:
            vals_filtered = vals[(vals >= ACTUAL_MIN) & (vals <= ACTUAL_MAX)]
            
            if len(vals_filtered) > 0:
                # Create histogram (counts, not density)
                counts, bin_edges = np.histogram(vals_filtered, bins=N_BINS, range=(ACTUAL_MIN, ACTUAL_MAX))
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                bin_width = bin_edges[1] - bin_edges[0]
                
                # Transform x-coordinates to section coordinates
                bin_centers_transformed = x_start + bin_centers
                
                # Plot histogram bars (actual counts)
                ax.bar(bin_centers_transformed, counts, width=bin_width * 0.9, 
                      color=color, alpha=HISTOGRAM_ALPHA, 
                      edgecolor=color, linewidth=HISTOGRAM_EDGE_WIDTH)
                
                # Add KDE overlay scaled to histogram count scale
                if len(vals_filtered) >= 2:
                    try:
                        kde = gaussian_kde(vals_filtered, bw_method='scott')
                        x_actual_grid = np.linspace(ACTUAL_MIN, ACTUAL_MAX, 300)
                        x_transformed_grid = x_start + x_actual_grid
                        
                        # CORRECT SCALING: KDE * n * bin_width
                        kde_scaled = kde(x_actual_grid) * len(vals_filtered) * bin_width
                        
                        # Plot KDE as thin line overlay
                        ax.plot(x_transformed_grid, kde_scaled, 
                               color=color, linewidth=KDE_LINE_WIDTH, alpha=KDE_ALPHA)
                        
                    except Exception as e:
                        print(f"  Warning: KDE failed for {salinity}: {e}")
        
        # Add x-axis ticks (ONLY black ticks - NO colored ticks)
        for tick_val in X_TICK_POSITIONS:
            tick_x = x_start + tick_val
            ax.plot([tick_x, tick_x], [0, -1.5], color='black', linewidth=1.5)
            
            # Add tick labels with proper spacing
            ax.text(tick_x, -2.5, f'{tick_val}',
                   ha='center', va='top', fontsize=X_TICK_LABEL_SIZE)
        
        # Add annotation (n sites and total observations)
        n_sites = len(vals)
        total_obs = np.sum(vals)
        
        # Position annotation at right side
        x_annotation = x_start + (SECTION_WIDTH * 0.83)
        y_annotation = max_count * 0.70
        
        ax.text(x_annotation, y_annotation,
               f'N = {n_sites}\n∑ = {total_obs:.0f} months',
               ha='center', va='center',
               fontsize=ANNOTATION_FONT_SIZE,
               color=color)
    
    # Add ecosystem labels above each section
    y_label_pos = max_count * 0.84
    for idx, salinity in enumerate(DISPLAY_ORDER):
        center = section_centers[idx] + 6
        ax.text(center, y_label_pos, salinity,
               ha='center', va='bottom',
               fontsize=ECOSYSTEM_LABEL_SIZE,
               fontweight='bold',
               color=SALINITY_COLORS[salinity])
    
    # Set y-axis
    ax.set_ylim(-4, max_count)
    ax.set_ylabel('Number of sites', fontsize=Y_AXIS_LABEL_SIZE, fontweight='bold', labelpad=15)
    
    # Set y-axis ticks (integer values for counts)
    if max_count <= 10:
        y_ticks = np.arange(0, max_count + 2, 2)
    elif max_count <= 20:
        y_ticks = np.arange(0, max_count + 5, 5)
    else:
        y_ticks = np.arange(0, max_count + 5, 5)
    
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{int(tick)}' for tick in y_ticks], fontsize=Y_TICK_LABEL_SIZE)
    
    # Set x-axis limits
    ax.set_xlim(-8, section_x_ends[-1] + 8)
    
    # Add x-axis label (centered under middle section)
    x_label_pos = (section_x_starts[1] + section_x_ends[1]) / 2
    ax.text(x_label_pos, -5.5, 
           'Near-normal (NN) observations per site',
           ha='center', va='top', fontsize=X_AXIS_LABEL_SIZE)
    
    # Style axes - thick border on all sides
    for spine in ax.spines.values():
        spine.set_linewidth(2.5)
        spine.set_color('black')
    
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Tick parameters
    ax.tick_params(axis='y', which='major', labelsize=Y_TICK_LABEL_SIZE, length=8, width=2)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
    # Add light grid on y-axis only
    ax.grid(True, linestyle=':', alpha=0.25, axis='y', linewidth=0.8)
    
    # Add panel label (c)
    fig.text(0.09, 0.92, '(c)', fontsize=32, fontweight='bold',
            transform=fig.transFigure,
            verticalalignment='top', horizontalalignment='left')
    
    # Adjust layout - increased bottom margin to prevent overlap
    plt.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.22)
    
    return fig

# =============================================================================
# PRINT STATISTICS (UPDATED TERMINOLOGY)
# =============================================================================

def print_sample_size_statistics(sample_size_data, shared_sites):
    """Print comprehensive statistics about site sample sizes"""
    
    print("\n" + "="*80)
    print("PANEL C: SAMPLE SIZE DISTRIBUTION STATISTICS")
    print("="*80)
    print(f"Strict triple intersection sites: {len(shared_sites)}")
    print("Order: Upland (left), Freshwater (middle), Saline (right)")
    print("NOTE: Strict metric filter applied (WUE, WUE_E, WUE_T ALL have data)")
    print("="*80)
    
    summary_data = []
    
    for salinity in DISPLAY_ORDER:
        values = sample_size_data[salinity]
        n_sites = len(values)
        
        if n_sites > 0:
            median_val = np.median(values)
            mean_val = np.mean(values)
            std_val = np.std(values, ddof=1) if n_sites > 1 else 0
            se_val = std_val / np.sqrt(n_sites) if n_sites > 1 else 0
            min_val = np.min(values)
            max_val = np.max(values)
            q25 = np.percentile(values, 25)
            q75 = np.percentile(values, 75)
            total_obs = np.sum(values)
            
            print(f"\n  {salinity}:")
            print(f"    n_sites = {n_sites}")
            print(f"    Total observations = {total_obs:.0f}")
            print(f"    Median = {median_val:.1f} observations")
            print(f"    Mean ± SE = {mean_val:.1f} ± {se_val:.1f} observations")
            print(f"    SD = {std_val:.1f} observations")
            print(f"    Q1-Q3 = {q25:.1f} - {q75:.1f} observations")
            print(f"    Range = {min_val:.0f} - {max_val:.0f} observations")
            
            summary_data.append({
                'Salinity': salinity,
                'n_sites': n_sites,
                'Total_observations': total_obs,
                'Median_observations': median_val,
                'Mean_observations': mean_val,
                'SE_observations': se_val,
                'SD_observations': std_val,
                'Q1_observations': q25,
                'Q3_observations': q75,
                'Min_observations': min_val,
                'Max_observations': max_val
            })
        else:
            print(f"\n  {salinity}: No data available")
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = os.path.join(OUTPUT_DIR, "PanelC_sample_size_distribution_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n✅ Saved summary statistics to: {summary_csv}")
    
    return summary_df

# =============================================================================
# SAVE FIGURE
# =============================================================================

def save_figure(fig, filename="PanelC_SampleSize_Distributions"):
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
# MAIN WORKFLOW
# =============================================================================

def main():
    """Main workflow for Panel C"""
    
    print("="*80)
    print("PANEL C: SAMPLE SIZE DISTRIBUTIONS")
    print("="*80)
    print("Order: Upland (left) | Freshwater (middle) | Saline (right)")
    print("Colors: EXACT same as Figure 2")
    print("Y-axis: Number of sites (counts)")
    print("="*80)
    
    # Step 1: Load data
    print("\n[Step 1] Loading site sample size data...")
    sample_size_data, shared_sites = load_site_sample_sizes()
    
    # Step 2: Print statistics
    print("\n[Step 2] Calculating sample size statistics...")
    summary_df = print_sample_size_statistics(sample_size_data, shared_sites)
    
    # Step 3: Create continuous panel
    print("\n[Step 3] Creating continuous panel with histograms...")
    fig = create_continuous_panel(sample_size_data, shared_sites)
    
    # Step 4: Save figure
    print("\n[Step 4] Saving Panel C figure...")
    save_figure(fig, "PanelC_SampleSize_Distributions_Histogram")
    
    print("\n" + "="*80)
    print("PANEL C WORKFLOW COMPLETE")
    print("="*80)
    
    return sample_size_data, shared_sites, fig

# =============================================================================
# EXECUTE
# =============================================================================

if __name__ == "__main__":
    sample_size_data, shared_sites, figure = main()