# -*- coding: utf-8 -*-
"""
Panel D (Supplementary): T:ET Ratio Density Distribution by Ecosystem Type
FOUR SUBPLOTS sharing y-axis (cleaner than single continuous panel)
ORDER: Upland | Freshwater | Brackish | Saline
Each subplot has x-axis 0-1 with ticks at 0, 0.5, 1.0
FIXED: Added strict month filter (same as CHUNK 3) to match Panel C
FIXED: Now counts ALL retained observations consistent with other panels
FIXED: Updated verification block to match corrected Panel C totals
@author: ammar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import gaussian_kde

# =============================================================================
# CONFIGURATION
# =============================================================================

# Figure dimensions
FIGURE_WIDTH = 15
FIGURE_HEIGHT = 6

# Font sizes
X_TICK_LABEL_SIZE = 18
Y_TICK_LABEL_SIZE = 20
X_AXIS_LABEL_SIZE = 20
Y_AXIS_LABEL_SIZE = 22
ECOSYSTEM_TITLE_SIZE = 24
ANNOTATION_FONT_SIZE = 18

# Histogram styling
HISTOGRAM_ALPHA = 0.5
HISTOGRAM_EDGE_WIDTH = 1

# KDE overlay styling
KDE_LINE_WIDTH = 1.5
KDE_ALPHA = 0.6

# Colors for ECOSYSTEM TYPES
ECOSYSTEM_COLORS = {
    'Upland': '#800080',     # Purple
    'Freshwater': '#0000FF',  # Blue
    'Brackish': '#008080',    # Teal
    'Saline': '#FFA500'       # Orange
}

# Order for display (left to right)
DISPLAY_ORDER = ['Upland', 'Freshwater', 'Brackish', 'Saline']

# File paths
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\figures\supplementary"

# Input data files
BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
MONTHLY_FILE = os.path.join(BASE_DIR, "monthly_data_after_outlier_removal.csv")
NN_MEDIANS_FILE = os.path.join(BASE_DIR, "wue_site_level_NN_medians_SPEI_1.csv")

# Output settings
DPI = 600
OUTPUT_FORMAT = 'both'

# T:ET ratio range (0-1)
TET_MIN = 0
TET_MAX = 1

# Number of bins for histogram
N_BINS = 10

# =============================================================================
# STEP 1: GET STRICT TRIPLE INTERSECTION SITES
# =============================================================================

print("="*80)
print("PANEL D: T:ET RATIO DISTRIBUTIONS (4 SUBPLOTS)")
print("="*80)

# Load the NN medians file
nn_medians = pd.read_csv(NN_MEDIANS_FILE)
print(f"Loaded NN medians file: {len(nn_medians)} rows")

# Get sites with ALL THREE metrics (strict triple intersection)
wue_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())

shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
print(f"Strict triple intersection sites: {len(shared_sites)}")

# =============================================================================
# STEP 2: LOAD MONTHLY DATA AND APPLY STRICT MONTH FILTER (SAME AS CHUNK 3)
# =============================================================================

print("\nLoading monthly data...")
df_monthly = pd.read_csv(MONTHLY_FILE)
print(f"  Loaded {len(df_monthly):,} rows")

# Filter to NN conditions
if 'SPEI_1_Cat' in df_monthly.columns:
    df_nn = df_monthly[df_monthly['SPEI_1_Cat'] == 'NN'].copy()
else:
    nn_mask = (df_monthly["SPEI_1"] >= -0.5) & (df_monthly["SPEI_1"] <= 0.5)
    df_nn = df_monthly[nn_mask].copy()

print(f"  After NN filter: {len(df_nn):,} rows")

# =============================================================================
# CRITICAL FIX: APPLY STRICT MONTH FILTER (SAME AS CHUNK 3)
# Only keep rows where WUE, WUE_eva, AND WUE_tra ALL have data
# =============================================================================
print("\nApplying strict month filter (same as CHUNK 3)...")
strict_mask = (
    df_nn['WUE'].notna() &
    df_nn['WUE_eva'].notna() &
    df_nn['WUE_tra'].notna()
)
df_nn_strict = df_nn[strict_mask].copy()
print(f"  Rows after strict month filter: {len(df_nn_strict):,}")
print(f"  Rows removed: {len(df_nn) - len(df_nn_strict)} ({((len(df_nn) - len(df_nn_strict))/len(df_nn)*100):.1f}%)")

# Filter to strict triple intersection sites
df_nn_strict = df_nn_strict[df_nn_strict['site_name'].isin(shared_sites)].copy()
print(f"\nAfter site filter: {len(df_nn_strict):,} rows")
print(f"Unique sites: {df_nn_strict['site_name'].nunique()} (expected: {len(shared_sites)})")

# Check Salinity_Category column
if 'Salinity_Category' not in df_nn_strict.columns:
    print("ERROR: Salinity_Category column not found!")
    raise KeyError("Salinity_Category column missing")

print(f"Ecosystem types: {df_nn_strict['Salinity_Category'].unique()}")

# =============================================================================
# STEP 3: EXTRACT T:ET RATIO VALUES BY ECOSYSTEM
# =============================================================================

print("\nExtracting T:ET ratio values...")
tet_data = {}

for ecosystem in DISPLAY_ORDER:
    mask = (df_nn_strict['Salinity_Category'] == ecosystem)
    values = df_nn_strict[mask]['Trans_ratio'].dropna().values
    tet_data[ecosystem] = values
    
    if len(values) > 0:
        print(f"  {ecosystem}: n_obs={len(values)}, range={np.min(values):.3f}-{np.max(values):.3f}, median={np.median(values):.3f}")
    else:
        print(f"  {ecosystem}: No data")

# =============================================================================
# STEP 4: CREATE 4 SUBPLOTS SHARING Y-AXIS
# =============================================================================

print("\nCreating 4 subplots with shared y-axis...")

fig, axes = plt.subplots(1, 4, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), sharey=True)

# Calculate global y-axis max for consistent scaling
global_ymax = 0
for ecosystem in DISPLAY_ORDER:
    vals = tet_data[ecosystem]
    if len(vals) >= 2:
        counts, _ = np.histogram(vals, bins=N_BINS, range=(TET_MIN, TET_MAX))
        global_ymax = max(global_ymax, np.max(counts))
global_ymax = global_ymax * 1.25 if global_ymax > 0 else 15

# Plot each ecosystem in its own subplot
for idx, ecosystem in enumerate(DISPLAY_ORDER):
    ax = axes[idx]
    vals = tet_data[ecosystem]
    color = ECOSYSTEM_COLORS[ecosystem]
    
    if len(vals) >= 2:
        # Histogram
        counts, bin_edges = np.histogram(vals, bins=N_BINS, range=(TET_MIN, TET_MAX))
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = bin_edges[1] - bin_edges[0]
        
        # Plot histogram bars
        ax.bar(bin_centers, counts, width=bin_width * 0.9,
              color=color, alpha=HISTOGRAM_ALPHA,
              edgecolor=color, linewidth=HISTOGRAM_EDGE_WIDTH)
        
        # KDE overlay
        try:
            kde = gaussian_kde(vals, bw_method='scott')
            x_grid = np.linspace(TET_MIN, TET_MAX, 300)
            density = kde(x_grid)
            kde_scaled = density * len(vals) * bin_width
            ax.plot(x_grid, kde_scaled,
                   color=color, linewidth=KDE_LINE_WIDTH, alpha=KDE_ALPHA)
        except Exception as e:
            print(f"  KDE failed for {ecosystem}: {e}")
    
    # Set x-axis limits and ticks
    ax.set_xlim(TET_MIN, TET_MAX)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(['0', '0.5', '1.0'], fontsize=X_TICK_LABEL_SIZE)
    
    # Set y-axis limits (only for first subplot, others share)
    if idx == 0:
        ax.set_ylim(0, global_ymax)
        ax.set_ylabel('Number of observations', fontsize=Y_AXIS_LABEL_SIZE, fontweight='bold')
    else:
        ax.tick_params(axis='y', labelleft=False)
    
    # Set x-axis label only for bottom row
    ax.set_xlabel('T:ET ratio', fontsize=X_AXIS_LABEL_SIZE)
    
    # Set title (ecosystem name with color)
    ax.set_title(ecosystem, fontsize=ECOSYSTEM_TITLE_SIZE, fontweight='bold', color=color)
    
    # Add annotation with n_obs
    n_obs = len(vals)
    median_val = np.median(vals) if len(vals) > 0 else 0
    ax.text(0.95, 0.93, f'n = {n_obs}\nmedian = {median_val:.2f}',
           transform=ax.transAxes, ha='right', va='top',
           fontsize=ANNOTATION_FONT_SIZE, color=color,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none'))
    
    # Style axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # Tick parameters
    ax.tick_params(axis='x', which='major', labelsize=X_TICK_LABEL_SIZE, length=6, width=1.5)
    if idx == 0:
        ax.tick_params(axis='y', which='major', labelsize=Y_TICK_LABEL_SIZE, length=6, width=1.5)
    
    # Add light grid on y-axis only
    ax.grid(True, linestyle=':', alpha=0.25, axis='y', linewidth=0.8)

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.12, wspace=0.25)

# Save
os.makedirs(OUTPUT_DIR, exist_ok=True)
png_path = os.path.join(OUTPUT_DIR, "PanelD_TET_Ratio_Distributions_Subplots.png")
pdf_path = os.path.join(OUTPUT_DIR, "PanelD_TET_Ratio_Distributions_Subplots.pdf")
plt.savefig(png_path, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {png_path}")
print(f"✅ Saved: {pdf_path}")

plt.show()

# =============================================================================
# PRINT STATISTICS
# =============================================================================

print("\n" + "="*80)
print("PANEL D: T:ET RATIO DISTRIBUTION STATISTICS")
print("="*80)
print("NOTE: Strict month filter applied (same as CHUNK 3)")
print("      Only rows where WUE, WUE_E, WUE_T ALL have data")
print("="*80)

for ecosystem in DISPLAY_ORDER:
    values = tet_data[ecosystem]
    n_obs = len(values)
    
    if n_obs > 0:
        median_val = np.median(values)
        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1) if n_obs > 1 else 0
        se_val = std_val / np.sqrt(n_obs) if n_obs > 1 else 0
        min_val = np.min(values)
        max_val = np.max(values)
        q25 = np.percentile(values, 25)
        q75 = np.percentile(values, 75)
        
        print(f"\n  {ecosystem}:")
        print(f"    n_observations = {n_obs}")
        print(f"    Median = {median_val:.3f}")
        print(f"    Mean ± SE = {mean_val:.3f} ± {se_val:.3f}")
        print(f"    SD = {std_val:.3f}")
        print(f"    Q1-Q3 = {q25:.3f} - {q75:.3f}")
        print(f"    Range = {min_val:.3f} - {max_val:.3f}")

# =============================================================================
# VERIFICATION: Compare with Panel C totals (UPDATED)
# =============================================================================
print("\n" + "="*80)
print("VERIFICATION: Panel C vs Panel D Consistency")
print("="*80)

# Calculate totals from Panel D
total_obs_panel_d = sum(len(vals) for vals in tet_data.values())
print(f"Panel D total T:ET observations (strict filter): {total_obs_panel_d}")

# Panel C totals after strict metric filter (Brackish excluded from Panel C)
# Upland=270, Freshwater=245, Saline=127, Total=642
total_obs_panel_c = 642

# Get Brackish observations count
brackish_obs = len(tet_data.get('Brackish', []))
print(f"Panel C total NN observations (Upland+Freshwater+Saline): {total_obs_panel_c}")
print(f"Panel D Brackish observations: {brackish_obs}")

# Verify match
if total_obs_panel_d - brackish_obs == total_obs_panel_c:
    print(f"\n✓ MATCH: Panel D totals match Panel C after excluding Brackish observations")
    print(f"  {total_obs_panel_d} - {brackish_obs} = {total_obs_panel_d - brackish_obs} = {total_obs_panel_c}")
else:
    print(f"\n⚠️ MISMATCH: Panel D without Brackish ({total_obs_panel_d - brackish_obs}) != Panel C ({total_obs_panel_c})")
    print("   Check filtering consistency between panels")

# Print ecosystem breakdown for reference
print("\n" + "-"*50)
print("Ecosystem breakdown (Panel D):")
for ecosystem in DISPLAY_ORDER:
    n_obs = len(tet_data.get(ecosystem, []))
    print(f"  {ecosystem}: {n_obs} observations")
print("-"*50)

print("\n" + "="*80)
print("PANEL D WORKFLOW COMPLETE")
print("="*80)