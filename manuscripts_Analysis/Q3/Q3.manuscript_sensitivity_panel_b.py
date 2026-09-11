# -*- coding: utf-8 -*-
"""
generate_panel_B_Q3.py – FINAL WITH COMMON X‑AXIS SCALE

Panel B: Site sensitivity grouped by coastline.
- Exact Panel A coast colors with alpha=0.35 for box fills.
- Box outlines, whiskers, and caps now use the same coast color.
- Y‑axis coast labels only on left panel, colored black.
- Panel labels a), b), c) outside the plotting box.
- X‑axis tick marks visible.
- Common x‑axis range (0 to max, ticks every 25) across all panels.
- Diagnostic print for data check.
- Coast name font size reduced to 22 (only y‑tick labels on left panel).
- Added vertical dotted reference lines at 5%, 10%, and 20%.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_rgba
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)

# ============================================================================
# CONSTANTS – EXACT COAST COLORS FROM PANEL A
# ============================================================================
COAST_COLORS = {
    "Atlantic Coast": "#A50F15",
    "Pacific Coast":  "#0072B2",
    "Gulf Coast":     "#4D4D4D",
    "AK Coast":       "#009E73"
}


COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]

GRAY_POINT = "#595959"

# ============================================================================
# STYLE AND FONT SIZES – MATCH PANEL A
# ============================================================================
sns.set_style("white")
base_font = 30
tick_font = 28
plt.rcParams.update({
    'font.size': base_font,
    'axes.labelsize': base_font,
    'axes.titlesize': base_font + 4,
    'xtick.labelsize': tick_font,
    'ytick.labelsize': tick_font,
    'axes.titleweight': 'bold',
    'axes.grid': False,
})

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data for Panel B...")
site_sensitivity_rank = pd.read_csv(
    os.path.join(output_dir, "Q3_WUE_T_SPEI_site_sensitivity_rank.csv")
)

# ============================================================================
# CREATE PANEL B – MANUAL Y‑AXIS (sharey=False)
# ============================================================================
print("Creating Panel B...")

all_system_sensitivity = site_sensitivity_rank[
    site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES) &
    (site_sensitivity_rank['n_months'] >= 6)
]

# ----------------------------------------------------------------------------
# DIAGNOSTIC PRINT – check data by group
# ----------------------------------------------------------------------------
print("\nPanel B data check (count, min, median, max by timescale & coast):")
print(
    all_system_sensitivity
    .groupby(['SPEI_timescale', 'coast_region'])['mean_abs_pct_change']
    .agg(['count', 'min', 'median', 'max'])
)

# ----------------------------------------------------------------------------
# COMMON X‑AXIS SCALE FOR ALL PANELS
# ----------------------------------------------------------------------------
xmax = np.nanmax(all_system_sensitivity['mean_abs_pct_change'])
xmax = max(75, np.ceil(xmax / 25) * 25)   # at least 75, rounded up to multiple of 25
common_xticks = np.arange(0, xmax + 1, 25)

# ----------------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharey=False)
panel_labels = ['a)', 'b)', 'c)']

for idx, timescale in enumerate(SELECTED_TIMESCALES):
    ax = axes[idx]
    subset = all_system_sensitivity[all_system_sensitivity['SPEI_timescale'] == timescale]

    if len(subset) > 0:
        data_for_box = []
        for coast in COAST_REGION_LEVELS:
            coast_subset = subset[subset['coast_region'] == coast]
            if len(coast_subset) > 0:
                data_for_box.append(coast_subset['mean_abs_pct_change'].values)
            else:
                data_for_box.append([])

        positions = range(len(COAST_REGION_LEVELS))

        bp = ax.boxplot(
            data_for_box,
            positions=positions,
            widths=0.6,
            patch_artist=True,
            vert=False,
            showfliers=False,
            whiskerprops=dict(color='black', linewidth=1.2),   # temporary, will be overridden
            capprops=dict(color='black', linewidth=1.2),       # temporary
            medianprops=dict(color='black', linewidth=2)
        )

        # Now set colours for each box, whiskers, and caps
        for i, box in enumerate(bp['boxes']):
            coast = COAST_REGION_LEVELS[i]
            color = COAST_COLORS[coast]
            box.set_facecolor(to_rgba(color, 0.35))
            box.set_edgecolor(color)          # box outline in coast colour
            box.set_linewidth(1.2)

        # Set whisker and cap colours to coast colour
        whiskers = bp['whiskers']
        caps = bp['caps']
        for i, coast in enumerate(COAST_REGION_LEVELS):
            color = COAST_COLORS[coast]
            # Each box has two whiskers (lower and upper) -> indices 2*i and 2*i+1
            whiskers[2*i].set_color(color)
            whiskers[2*i+1].set_color(color)
            caps[2*i].set_color(color)
            caps[2*i+1].set_color(color)

        all_points = subset[['mean_abs_pct_change', 'coast_region']].dropna()
        for i, coast in enumerate(COAST_REGION_LEVELS):
            coast_points = all_points[all_points['coast_region'] == coast]['mean_abs_pct_change'].values
            if len(coast_points) > 0:
                y_jitter = np.random.normal(i, 0.08, len(coast_points))
                ax.scatter(
                    coast_points,
                    y_jitter,
                    color=GRAY_POINT,
                    s=38,
                    alpha=0.40,
                    edgecolors='none',
                    zorder=3
                )

        # Y‑axis: manual ticks and limits (Atlantic at top)
        ax.set_yticks(list(positions))
        ax.set_ylim(len(COAST_REGION_LEVELS) - 0.5, -0.5)

        if idx == 0:
            ax.set_yticklabels(COAST_REGION_LEVELS)
            for tick_label in ax.get_yticklabels():
                tick_label.set_color("black")
                tick_label.set_fontweight("bold")
                tick_label.set_fontsize(22)          # ← REDUCED COAST NAME FONT SIZE
            ax.tick_params(axis='y', length=6, pad=8, left=True)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis='y', left=False, labelleft=False)

        ax.set_ylabel('')
        ax.set_xlabel('')

        ax.set_title(f'SPEI-{timescale.split("_")[1]}', fontweight='bold', fontsize=25)

        # Use fixed x‑axis range and ticks (common for all panels)
        ax.set_xlim(0, xmax)
        ax.set_xticks(common_xticks)

        # Make tick marks visible
        ax.tick_params(
            axis='x',
            which='major',
            bottom=True,
            top=False,
            labelbottom=True,
            length=7,
            width=1.2,
            direction='out',
            colors='black'
        )
        ax.spines['bottom'].set_visible(True)
        ax.spines['bottom'].set_linewidth(1.0)

        # ===== ADD VERTICAL REFERENCE LINES =====
        # Dotted vertical lines at 5%, 10%, 20% (matching Malone's R workflow)
        for xref in [5, 10, 20]:
            ax.axvline(
                x=xref,
                color="0.65",
                linestyle=":",
                linewidth=1.0,
                alpha=1.0,
                zorder=1
            )
        # ========================================

        ax.grid(False)

        # Panel label – outside the plotting box, above top‑left
        ax.text(
            -0.10, 1.06, panel_labels[idx],
            transform=ax.transAxes,
            fontsize=tick_font + 2,
            fontweight='bold',
            va='bottom',
            ha='left',
            clip_on=False
        )

# ----------------------------------------------------------------------------
# SHARED X‑AXIS LABEL
# ----------------------------------------------------------------------------
fig.text(
    0.6, 0.055,
    r'Mean absolute WUE$_T$ change from near-normal (%)',
    ha='center',
    va='center',
    fontsize=24
)

# Layout with more bottom space
fig.subplots_adjust(left=0.30, bottom=0.24, top=0.84, wspace=0.22)

# ----------------------------------------------------------------------------
# SAVE AND DISPLAY
# ----------------------------------------------------------------------------
output_path = os.path.join(figure_dir, "Q3_panel_B_site_sensitivity_coastline.png")
fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.show(block=True)
print(f"Panel B saved to: {output_path}")