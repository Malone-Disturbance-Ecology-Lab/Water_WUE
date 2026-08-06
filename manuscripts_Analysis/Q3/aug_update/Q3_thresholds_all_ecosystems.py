# -*- coding: utf-8 -*-
"""
Q3_figure_panel_A_updated_all_ecosystems.py

Updated Panel A figure using reviewer outputs:
- All ecosystems (Upland, Freshwater, Saline)
- All months (month‑averaged curves)
- Thresholds plotted: 5%, 10%, 20% + top two highest detected above 20%
- No shaded ribbons (curves only)
- Panel labels moved to left margin (a–d for rows)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
input_dir = os.path.join(base_dir, "Q3_august_update")
figure_dir = os.path.join(base_dir, "Q3_WUE_T_SPEI_sensitivity_figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)

# ============================================================================
# CONSTANTS
# ============================================================================
COAST_COLORS = {
    'Atlantic Coast': '#A50F15',
    'Pacific Coast':  '#0072B2',
    'Gulf Coast':     '#4D4D4D',
    'AK Coast':       '#009E73'
}

ECOSYSTEM_COLORS = {
    'Upland':      '#800080',
    'Freshwater':  '#0000FF',
    'Saline':      '#FFA500'
}

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]

COAST_REGION_LEVELS = [
    "Atlantic Coast",
    "Pacific Coast",
    "Gulf Coast",
    "AK Coast"
]

COAST_DISPLAY_LABELS = {
    "Atlantic Coast": "Atlantic Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "AK Coast": "Alaska Coast"
}

# Row labels (for left margin)
ROW_LABELS = ['a', 'b', 'c', 'd']

SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]
TIMESCALE_LABELS = {"SPEI_1": "SPEI-1", "SPEI_3": "SPEI-3", "SPEI_48": "SPEI-48"}

# ---- full marker map for all possible threshold levels ----
THRESHOLD_MARKER_MAP = {
    "5%":  "o",
    "10%": "s",
    "15%": "v",
    "20%": "^",
    "25%": "D",
    "30%": "P",
    "35%": "X",
    "40%": "*",
    "50%": "h",
    "75%": "8"
}

REF_GRAY = "#9A9A9A"
ZERO_GRAY = "#707070"

def get_threshold_marker(threshold_pct):
    label = str(threshold_pct).strip()
    return THRESHOLD_MARKER_MAP.get(label, "o")

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading reviewer data (all ecosystems, all months, averaged)...")
curve_file = os.path.join(input_dir, "Q3_reviewer_coast_threshold_prediction_curves_month_averaged_all_ecosystems.csv")
marker_file = os.path.join(input_dir, "Q3_reviewer_threshold_markers_month_averaged_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv")

curve_data = pd.read_csv(curve_file)
marker_data = pd.read_csv(marker_file)

# Ensure threshold_pct is string with '%' for marker mapping
if 'threshold_pct' in marker_data.columns:
    if marker_data['threshold_pct'].dtype in ['int64', 'float64']:
        marker_data['threshold_pct'] = marker_data['threshold_pct'].astype(int).astype(str) + '%'
    elif not marker_data['threshold_pct'].astype(str).str.endswith('%').all():
        marker_data['threshold_pct'] = marker_data['threshold_pct'].astype(int).astype(str) + '%'

# Filter to selected timescales
curve_data = curve_data[curve_data['SPEI_timescale'].isin(SELECTED_TIMESCALES)]
marker_data = marker_data[marker_data['SPEI_timescale'].isin(SELECTED_TIMESCALES)]

print(f"Curve data: {len(curve_data)} rows")
print(f"Marker data (before threshold selection): {len(marker_data)} rows")

# ============================================================================
# DETERMINE WHICH THRESHOLDS TO PLOT
# ============================================================================
# Extract unique detected threshold percentages (as integers)
detected_thresholds = marker_data['threshold_pct'].str.replace('%', '').astype(int).unique()
detected_thresholds = sorted(detected_thresholds)
print(f"\nAll detected threshold levels in marker data: {detected_thresholds}")

# Base thresholds always plotted
base_thresholds = [5, 10, 20]

# Thresholds above 20%
above_20 = [t for t in detected_thresholds if t > 20]
top_two_above_20 = sorted(above_20, reverse=True)[:2]  # highest two

# Final selection
selected_thresholds = sorted(set(base_thresholds + top_two_above_20))
print(f"Selected thresholds for plotting: {selected_thresholds}")

# Convert to string labels for filtering and legends
selected_labels = [f"{t}%" for t in selected_thresholds]

# Filter marker_data to only selected thresholds
marker_data = marker_data[marker_data['threshold_pct'].isin(selected_labels)]

print(f"Marker data after threshold selection: {len(marker_data)} rows")

# Count markers per selected threshold for console summary
marker_counts = marker_data.groupby('threshold_pct').size()
print("\nNumber of markers per selected threshold:")
for lbl in selected_labels:
    count = marker_counts.get(lbl, 0)
    print(f"  {lbl}: {count}")

# ============================================================================
# STYLE – NO GRID, WITH VISIBLE TICKS
# ============================================================================
sns.set_style("white")
base_font = 30
tick_font = 28
legend_font = 26
plt.rcParams.update({
    'font.size': base_font,
    'axes.labelsize': base_font,
    'axes.titlesize': base_font + 4,
    'xtick.labelsize': tick_font,
    'ytick.labelsize': tick_font,
    'legend.fontsize': legend_font,
    'legend.title_fontsize': legend_font,
    'axes.titleweight': 'bold',
    'axes.grid': False,
    'xtick.bottom': True,
    'ytick.left': True,
    'xtick.top': False,
    'ytick.right': False,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.edgecolor': 'black',
    'axes.linewidth': 1.5,
})

# ============================================================================
# CREATE PANEL A – 4 ROWS x 3 COLUMNS
# ============================================================================
fig, axes = plt.subplots(4, 3, figsize=(18, 16))

# Adjust subplots to make room for row labels
plt.subplots_adjust(left=0.15, right=0.78, top=0.95, bottom=0.08, wspace=0.3, hspace=0.4)

# Handles for legends
eco_handles = []
eco_labels = []
marker_handles = []
marker_labels = []

for idx_c, coast in enumerate(COAST_REGION_LEVELS):
    for idx_t, timescale in enumerate(SELECTED_TIMESCALES):
        ax = axes[idx_c, idx_t]

        ax.grid(False)
        ax.minorticks_off()
        ax.set_facecolor("white")

        # ------------------------------------------------------------------
        # REFERENCE LINES (SPEI anomalies, zero, and threshold lines)
        # ------------------------------------------------------------------
        ax.axvline(x=-1, color=REF_GRAY, linestyle='--', linewidth=1.0, alpha=0.85, zorder=2)
        ax.axvline(x=1,  color=REF_GRAY, linestyle='--', linewidth=1.0, alpha=0.85, zorder=2)
        ax.axhline(y=0, color=ZERO_GRAY, linestyle='-', linewidth=1.1, alpha=0.9, zorder=2)

        # ---- horizontal reference lines for the selected thresholds ----
        for thr in selected_thresholds:
            for yref in [-thr, thr]:
                ax.axhline(
                    y=yref,
                    color="#595959",
                    linestyle="--",
                    linewidth=0.9,
                    alpha=0.75,
                    zorder=2.4
                )

        # ------------------------------------------------------------------
        # CURVES FOR EACH ECOSYSTEM
        # ------------------------------------------------------------------
        for ecosystem in ECOSYSTEM_CLASSES:
            sub = curve_data[
                (curve_data['coast_region'] == coast) &
                (curve_data['SPEI_timescale'] == timescale) &
                (curve_data['water_class'] == ecosystem)
            ]
            if len(sub) == 0:
                continue
            sub = sub.sort_values("SPEI_value")
            color = ECOSYSTEM_COLORS[ecosystem]

            ax.plot(
                sub['SPEI_value'],
                sub['mean_pct'],
                color=color,
                linewidth=1.8,
                label=ecosystem,
                zorder=3
            )

            # ------------------------------------------------------------------
            # THRESHOLD MARKERS (only for selected thresholds)
            # ------------------------------------------------------------------
            markers_sub = marker_data[
                (marker_data['coast_region'] == coast) &
                (marker_data['SPEI_timescale'] == timescale) &
                (marker_data['water_class'] == ecosystem)
            ]
            for _, row in markers_sub.iterrows():
                thr_label = row['threshold_pct']
                marker = get_threshold_marker(thr_label)
                x_pos = row['SPEI_threshold']
                y_pos = row['pct_change_threshold']
                x_low = row['SPEI_threshold_lower']
                x_high = row['SPEI_threshold_upper']

                # Horizontal bar indicating uncertainty range
                if not np.isnan(x_low) and not np.isnan(x_high):
                    ax.hlines(
                        y=y_pos,
                        xmin=x_low,
                        xmax=x_high,
                        color=color,
                        linewidth=1.5,
                        alpha=0.7,
                        zorder=5
                    )

                # Marker point
                ax.scatter(
                    x=x_pos,
                    y=y_pos,
                    color=color,
                    s=120,
                    marker=marker,
                    edgecolors='white',
                    linewidths=0.7,
                    zorder=6
                )

                # Collect legend handles for each unique threshold label
                if thr_label not in marker_labels:
                    marker_labels.append(thr_label)
                    marker_handles.append(
                        plt.Line2D(
                            [0], [0],
                            marker=marker,
                            color='w',
                            markerfacecolor='gray',
                            markersize=14,
                            label=thr_label
                        )
                    )

        # Collect ecosystem legend handles
        for eco in ECOSYSTEM_CLASSES:
            if eco not in eco_labels:
                eco_labels.append(eco)
                eco_handles.append(
                    plt.Line2D([0], [0], color=ECOSYSTEM_COLORS[eco], lw=3, label=eco)
                )

        # ------------------------------------------------------------------
        # AXIS LABELS, TICKS
        # ------------------------------------------------------------------
        if idx_t == 0:
            ax.set_ylabel(COAST_DISPLAY_LABELS[coast], fontweight='bold', fontsize=base_font)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.xaxis.set_major_locator(plt.MaxNLocator(5, integer=True))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5, integer=True))

        ax.tick_params(
            axis='both',
            which='major',
            direction='out',
            length=7,
            width=1.5,
            colors='black',
            bottom=True,
            left=True,
            top=False,
            right=False,
            pad=6
        )

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.5)

        if idx_c == 0:
            ax.set_title(TIMESCALE_LABELS[timescale], fontweight='bold', fontsize=base_font+4)

# ----------------------------------------------------------------------------
# ADD ROW LABELS (a, b, c, d) IN LEFT MARGIN
# ----------------------------------------------------------------------------
for i, coast in enumerate(COAST_REGION_LEVELS):
    ax0 = axes[i, 0]
    bbox = ax0.get_position()
    x_pos = bbox.x0 - 0.025   # adjust as needed
    y_pos = bbox.y0 + bbox.height / 2.0 + 0.1   # shifted up slightly
    fig.text(
        x_pos, y_pos,
        ROW_LABELS[i] + ')',
        fontsize=tick_font + 4,
        fontweight='bold',
        va='center',
        ha='right',
        transform=fig.transFigure
    )

# ----------------------------------------------------------------------------
# GLOBAL Y-AXIS LABEL
# ----------------------------------------------------------------------------
fig.text(
    0.01, 0.5,
    r'WUE$_{T}$ Change (%)',
    ha='center',
    va='center',
    rotation=90,
    fontsize=base_font
)

# ----------------------------------------------------------------------------
# LEGENDS (ecosystem and threshold markers)
# ----------------------------------------------------------------------------
if eco_handles:
    fig.legend(
        handles=eco_handles,
        labels=eco_labels,
        title='Ecosystem',
        loc='center left',
        bbox_to_anchor=(0.82, 0.65),
        framealpha=0.9,
        edgecolor='black',
        fontsize=legend_font,
        title_fontsize=legend_font
    )

if marker_handles:
    fig.legend(
        handles=marker_handles,
        labels=marker_labels,
        title='Threshold Change (%)',
        loc='center left',
        bbox_to_anchor=(0.82, 0.35),
        framealpha=0.9,
        edgecolor='black',
        fontsize=legend_font,
        title_fontsize=legend_font
    )

plt.suptitle('')

# ----------------------------------------------------------------------------
# SAVE
# ----------------------------------------------------------------------------
output_path = os.path.join(figure_dir, "Q3_panel_A_coast_GAM_thresholds_all_ecosystems.png")
fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.show()
print(f"Updated Panel A saved to: {output_path}")