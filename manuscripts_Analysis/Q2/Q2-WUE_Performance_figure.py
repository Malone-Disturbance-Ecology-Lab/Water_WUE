"""
Q2_diagnostic_figure.py - Visual inspection of Q2 performance metrics by coast

This script creates a 3x3 diagnostic figure using already saved Q2 output CSV files.
It does NOT recalculate Q2 metrics or save new CSV files.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("Q2 Diagnostic Figure - Visual Inspection Only")
print("="*60)

# ============================================================================
# PATHS
# ============================================================================

input_file = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q_performace metric_T_ET_relationship\Q2_performance_by_coast_site_level_data.csv"

print(f"\nLoading data from: {input_file}")

# ============================================================================
# LOAD DATA
# ============================================================================

df = pd.read_csv(input_file)

# Clean column names if needed
df.columns = df.columns.str.strip()

print(f"\nData loaded:")
print(f"  Rows: {len(df):,}")
print(f"  Columns: {df.columns.tolist()}")
print(f"  Sites: {df['site_name'].nunique():,}")
print(f"  Ecosystem classes: {df['water_class'].unique().tolist()}")
print(f"  Coast regions: {df['coast_region'].unique().tolist()}")

# ============================================================================
# CONSTANTS
# ============================================================================

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
ECOSYSTEM_COLORS = {
    "Upland": "#800080",     # Purple
    "Freshwater": "#0000FF", # Blue
    "Saline": "#FFA500"      # Orange
}

# Coast region order
COAST_ORDER = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
COAST_LABELS_SHORT = ["Atlantic", "Pacific", "Gulf", "AK"]

# Distinct colors for each coast (consistent across panels)
COAST_COLORS = {
    "Atlantic Coast": "#1f77b4",  # Blue
    "Pacific Coast": "#ff7f0e",   # Orange
    "Gulf Coast": "#2ca02c",      # Green
    "AK Coast": "#d62728"         # Red
}

# Lighter versions for violin fills
COAST_COLORS_LIGHT = {
    "Atlantic Coast": "#4a9bd9",
    "Pacific Coast": "#ffa64d",
    "Gulf Coast": "#5cb85c",
    "AK Coast": "#e64a4a"
}

# ============================================================================
# THEME
# ============================================================================

def theme_wue(base_size=12):
    """Malone-style theme"""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': base_size,
        'axes.labelsize': base_size,
        'axes.titlesize': base_size + 2,
        'xtick.labelsize': base_size - 1,
        'ytick.labelsize': base_size - 1,
        'legend.fontsize': base_size - 1,
        'figure.titlesize': base_size + 4,
        'axes.titleweight': 'bold',
    })

theme_wue()

# ============================================================================
# CREATE FIGURE - 3x3 LAYOUT
# ============================================================================

print("\n" + "="*60)
print("Creating 3x3 diagnostic figure")
print("="*60)

fig = plt.figure(figsize=(15, 12))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1],
                      hspace=0.35, wspace=0.3)

# -------------------------------------------------------------------------
# PANEL A: Stability (Top Left) - UNCHANGED
# -------------------------------------------------------------------------
ax_a = fig.add_subplot(gs[0, 0])

# Plot points by ecosystem with colors
for ecosystem in ECOSYSTEM_CLASSES:
    subset = df[df['water_class'] == ecosystem]
    if len(subset) > 0:
        ax_a.scatter(subset['mean_TET'], subset['stability'], 
                    color=ECOSYSTEM_COLORS[ecosystem], alpha=0.75, s=50, 
                    label=ecosystem)
        
        # Add trend line per ecosystem
        if len(subset) >= 2:
            x = subset['mean_TET'].values
            y = subset['stability'].values
            # Fit polynomial in log space
            log_y = np.log10(y)
            coeffs = np.polyfit(x, log_y, 1)
            poly = np.poly1d(coeffs)
            x_smooth = np.linspace(x.min(), x.max(), 100)
            y_smooth = 10 ** poly(x_smooth)
            ax_a.plot(x_smooth, y_smooth, color=ECOSYSTEM_COLORS[ecosystem], 
                     linewidth=2, alpha=0.7)

ax_a.set_xlabel('Mean T:ET ratio')
ax_a.set_ylabel('Stability index (log scale)')
ax_a.set_xlim(0, 1)
ax_a.set_yscale('log')
ax_a.set_title('A: Stability', fontweight='bold')
ax_a.legend(loc='best', title='Ecosystem class')
ax_a.grid(True, alpha=0.3)

# -------------------------------------------------------------------------
# PANEL B: Stability by Coast (Top Middle) - UNCHANGED
# -------------------------------------------------------------------------
ax_b = fig.add_subplot(gs[0, 1])

# Prepare data by coast
coast_data_b = []
for coast in COAST_ORDER:
    subset = df[df['coast_region'] == coast]['stability'].dropna().values
    coast_data_b.append(subset)

# Create violin + boxplot with coast colors
positions = range(len(COAST_ORDER))
for i, (coast, data) in enumerate(zip(COAST_ORDER, coast_data_b)):
    if len(data) > 0:
        # Violin plot
        parts = ax_b.violinplot(data, positions=[i], widths=0.6, showmeans=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(COAST_COLORS_LIGHT[coast])
            pc.set_alpha(0.3)
        
        # Boxplot with coast color
        bp = ax_b.boxplot(data, positions=[i], widths=0.18, patch_artist=True,
                         boxprops=dict(facecolor=COAST_COLORS[coast], alpha=0.7),
                         whiskerprops=dict(color=COAST_COLORS[coast]),
                         capprops=dict(color=COAST_COLORS[coast]),
                         medianprops=dict(color='black', linewidth=1.5))

ax_b.set_xticks(positions)
ax_b.set_xticklabels(COAST_LABELS_SHORT)
ax_b.set_ylabel('Stability index')
ax_b.set_title('B: Stability by coast', fontweight='bold')
ax_b.grid(True, alpha=0.3)

# Add legend for coast colors
legend_elements = [
    Patch(facecolor=COAST_COLORS['Atlantic Coast'], alpha=0.7, label='Atlantic Coast'),
    Patch(facecolor=COAST_COLORS['Pacific Coast'], alpha=0.7, label='Pacific Coast'),
    Patch(facecolor=COAST_COLORS['Gulf Coast'], alpha=0.7, label='Gulf Coast'),
    Patch(facecolor=COAST_COLORS['AK Coast'], alpha=0.7, label='AK Coast')
]
ax_b.legend(handles=legend_elements, loc='best', fontsize=9, title='Coast region')

# -------------------------------------------------------------------------
# PANEL C: Plasticity Slope - DOT PLOT WITH ERROR BARS (NEW STYLE)
# -------------------------------------------------------------------------
ax_c = fig.add_subplot(gs[0, 2])

# Calculate means and standard deviations for each coast
coast_stats_c = {}
for coast in COAST_ORDER:
    data = df[df['coast_region'] == coast]['plasticity_slope_SPEI3'].dropna().values
    if len(data) > 0:
        coast_stats_c[coast] = {
            'mean': np.mean(data),
            'std': np.std(data),
            'n': len(data),
            'data': data
        }

positions_c = range(len(COAST_ORDER))
for i, coast in enumerate(COAST_ORDER):
    if coast in coast_stats_c:
        stats = coast_stats_c[coast]
        # Add jittered points (show all values)
        x_jitter = np.random.normal(positions_c[i], 0.08, size=len(stats['data']))
        ax_c.scatter(x_jitter, stats['data'], color=COAST_COLORS[coast], 
                    alpha=0.4, s=30, zorder=1)
        
        # Add mean with error bar
        ax_c.errorbar(positions_c[i], stats['mean'], yerr=stats['std'], 
                     fmt='o', color=COAST_COLORS[coast], 
                     markersize=12, markeredgecolor='black', markeredgewidth=2,
                     capsize=6, capthick=2, elinewidth=2, zorder=10)
        
        # Add n label
        ax_c.text(positions_c[i], ax_c.get_ylim()[0] + 0.05, f'n={stats["n"]}', 
                 ha='center', fontsize=8, color=COAST_COLORS[coast])

ax_c.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
ax_c.set_xticks(positions_c)
ax_c.set_xticklabels(COAST_LABELS_SHORT)
ax_c.set_ylabel('SPEI-3 slope')
ax_c.set_title('C: Plasticity slope (Mean ± SD)', fontweight='bold')
ax_c.grid(True, alpha=0.3)

# -------------------------------------------------------------------------
# PANEL D: Plasticity Range - SWARM + VIOLIN PLOT (CHANGED STYLE)
# -------------------------------------------------------------------------
ax_d = fig.add_subplot(gs[1, 0])

# Prepare data by coast
coast_data_d = {}
for coast in COAST_ORDER:
    coast_data_d[coast] = df[df['coast_region'] == coast]['plasticity_p95_p05'].dropna().values

positions_d = range(len(COAST_ORDER))

# Create swarm plot with violin overlay
for i, coast in enumerate(COAST_ORDER):
    data = coast_data_d[coast]
    if len(data) > 0:
        # Violin plot (transparent)
        parts = ax_d.violinplot(data, positions=[i], widths=0.6, showmeans=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(COAST_COLORS_LIGHT[coast])
            pc.set_alpha(0.2)
            pc.set_edgecolor(COAST_COLORS[coast])
            pc.set_linewidth(1)
        
        # Swarm plot (jittered points with controlled positions)
        # Use deterministic jitter for consistent visualization
        n_points = len(data)
        if n_points > 0:
            # Sort data for better visualization
            sorted_data = np.sort(data)
            # Create jitter positions
            jitter_width = 0.2
            positions_jitter = np.linspace(-jitter_width/2, jitter_width/2, n_points)
            # Add small random jitter
            positions_jitter = positions_jitter + np.random.normal(0, 0.03, n_points)
            # Clip to avoid overlap
            positions_jitter = np.clip(positions_jitter, -jitter_width/2, jitter_width/2)
            
            ax_d.scatter(positions_d[i] + positions_jitter, sorted_data, 
                       color=COAST_COLORS[coast], alpha=0.7, s=40, 
                       edgecolors='black', linewidth=0.5, zorder=5)

ax_d.set_xticks(positions_d)
ax_d.set_xticklabels(COAST_LABELS_SHORT)
ax_d.set_ylabel('WUE_T 95th / 5th')
ax_d.set_title('D: Plasticity range (Violin + Swarm)', fontweight='bold')
ax_d.grid(True, alpha=0.3)

# -------------------------------------------------------------------------
# PANEL E: Drought Resistance - BAR PLOT WITH POINTS (NEW STYLE)
# -------------------------------------------------------------------------
ax_e = fig.add_subplot(gs[1, 1])

# Prepare data by coast for resistance
coast_stats_e = {}
for coast in COAST_ORDER:
    data = df[df['coast_region'] == coast]['resistance'].dropna().values
    if len(data) > 0:
        coast_stats_e[coast] = {
            'mean': np.mean(data),
            'std': np.std(data),
            'n': len(data),
            'data': data
        }

positions_e = range(len(COAST_ORDER))
bar_width = 0.6

# Create bar plot with error bars
for i, coast in enumerate(COAST_ORDER):
    if coast in coast_stats_e:
        stats = coast_stats_e[coast]
        
        # Add bar
        bar = ax_e.bar(positions_e[i], stats['mean'], width=bar_width,
                      color=COAST_COLORS[coast], alpha=0.5, 
                      edgecolor=COAST_COLORS[coast], linewidth=2)
        
        # Add error bar
        ax_e.errorbar(positions_e[i], stats['mean'], yerr=stats['std'], 
                     fmt='none', color='black', 
                     capsize=6, capthick=2, elinewidth=2, zorder=10)
        
        # Add jittered points (show all values)
        x_jitter = np.random.normal(positions_e[i], 0.1, size=len(stats['data']))
        ax_e.scatter(x_jitter, stats['data'], color=COAST_COLORS[coast], 
                    alpha=0.6, s=25, zorder=5, edgecolors='black', linewidth=0.5)
        
        # Add n label
        ax_e.text(positions_e[i], 0.05, f'n={stats["n"]}', 
                 ha='center', fontsize=8, color='black')

ax_e.axhline(y=1, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
ax_e.text(0.02, 1.02, 'Reference: 1.0', color='red', fontsize=9, transform=ax_e.transAxes)
ax_e.set_xticks(positions_e)
ax_e.set_xticklabels(COAST_LABELS_SHORT)
ax_e.set_ylabel('Resistance')
ax_e.set_title('E: Drought resistance (Bar ± SD)', fontweight='bold')
ax_e.grid(True, alpha=0.3)

# -------------------------------------------------------------------------
# PANEL F: Drought Recovery - STRIPPLOT WITH MEAN LINE (NEW STYLE)
# -------------------------------------------------------------------------
ax_f = fig.add_subplot(gs[1, 2])

# Prepare data by coast for recovery
coast_stats_f = {}
for coast in COAST_ORDER:
    data = df[df['coast_region'] == coast]['mean_recovery'].dropna().values
    if len(data) > 0:
        coast_stats_f[coast] = {
            'mean': np.mean(data),
            'std': np.std(data),
            'n': len(data),
            'data': data
        }

positions_f = range(len(COAST_ORDER))

# Create stripplot with mean lines
for i, coast in enumerate(COAST_ORDER):
    if coast in coast_stats_f:
        stats = coast_stats_f[coast]
        
        # Add jittered points (show all values)
        x_jitter = np.random.normal(positions_f[i], 0.08, size=len(stats['data']))
        ax_f.scatter(x_jitter, stats['data'], color=COAST_COLORS[coast], 
                    alpha=0.6, s=40, zorder=5, edgecolors='black', linewidth=0.5)
        
        # Add mean line
        ax_f.hlines(stats['mean'], positions_f[i] - 0.3, positions_f[i] + 0.3,
                   color=COAST_COLORS[coast], linewidth=3, zorder=10)
        
        # Add mean value text
        ax_f.text(positions_f[i], stats['mean'] + 0.02, f'μ={stats["mean"]:.2f}', 
                 ha='center', fontsize=8, color=COAST_COLORS[coast], fontweight='bold')
        
        # Add n label at bottom
        ax_f.text(positions_f[i], ax_f.get_ylim()[0] + 0.02, f'n={stats["n"]}', 
                 ha='center', fontsize=8, color='black')

ax_f.axhline(y=1, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
ax_f.text(0.02, 1.02, 'Reference: 1.0', color='red', fontsize=9, transform=ax_f.transAxes)
ax_f.set_xticks(positions_f)
ax_f.set_xticklabels(COAST_LABELS_SHORT)
ax_f.set_ylabel('Recovery')
ax_f.set_title('F: Drought recovery (Stripplot + Mean)', fontweight='bold')
ax_f.grid(True, alpha=0.3)

# -------------------------------------------------------------------------
# EMPTY PANELS (Row 2, Col 3 and Row 3)
# -------------------------------------------------------------------------
ax_empty1 = fig.add_subplot(gs[2, 0])
ax_empty1.axis('off')
ax_empty2 = fig.add_subplot(gs[2, 1])
ax_empty2.axis('off')
ax_empty3 = fig.add_subplot(gs[2, 2])
ax_empty3.axis('off')

# -------------------------------------------------------------------------
# FINISH
# -------------------------------------------------------------------------

fig.suptitle('Q2 WUE_T Performance Metrics Diagnostic', fontsize=16, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])

print("\nDisplaying figure...")
print("  Layout: 3x3 panels (bottom row intentionally blank)")
print("  - Panel A: Stability with ecosystem colors and trend lines (UNCHANGED)")
print("  - Panel B: Stability by coast (Violin + Boxplot) (UNCHANGED)")
print("  - Panel C: Plasticity slope - Dot plot with error bars (NEW STYLE)")
print("  - Panel D: Plasticity range - Violin + Swarm plot (CHANGED STYLE)")
print("  - Panel E: Drought resistance - Bar plot with error bars and points (NEW STYLE)")
print("  - Panel F: Drought recovery - Stripplot with mean lines (NEW STYLE)")
print("\nKey features:")
print("  - All panels show individual site values clearly")
print("  - Panels C, E, F use different plot types (dot plot, bar plot, stripplot)")
print("  - Panel D changed from boxplot to violin + swarm")
print("  - Consistent coast colors across all panels")
print("  - Summary statistics (mean, SD, n) are displayed")
print("  - Reference lines at 0 or 1 where appropriate")
print("\nThis figure is for visual inspection only.")
print("Close the figure window to end the script.")

plt.show()

print("\n" + "="*60)
print("Diagnostic figure displayed successfully!")
print("="*60)