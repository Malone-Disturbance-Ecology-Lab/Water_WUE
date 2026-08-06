# -*- coding: utf-8 -*-
"""
generate_panel_A_only.py – FINAL CLEAN VERSION
Only zero line, SPEI anomaly lines, GAM curve, ribbon, and threshold markers.
No gray dotted threshold lines.
ADDED: global y‑axis label centered on the left, no overlap with coast labels.
FIXED: visible x/y tick marks and axis borders.
FIXED: Alaska coast name consistency (AK Coast vs Alaska Coast).
ADDED: horizontal WUE_T threshold reference lines at ±5%, ±10%, ±20% (darker).
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
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)

# ============================================================================
# CONSTANTS
# ============================================================================
# Use consistent coast naming: "AK Coast" in code, display as "Alaska Coast"
COAST_COLORS = {
    'Atlantic Coast': '#A50F15',   # dark red
    'Pacific Coast':  '#0072B2',   # blue
    'Gulf Coast':     '#4D4D4D',   # dark charcoal
    'AK Coast':       '#009E73'    # bluish green
}

COAST_REGION_LEVELS = [
    "Atlantic Coast",
    "Pacific Coast",
    "Gulf Coast",
    "AK Coast"
]

# Display labels for y-axis (row titles)
COAST_DISPLAY_LABELS = {
    "Atlantic Coast": "Atlantic Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "AK Coast": "Alaska Coast"
}

SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]
TIMESCALE_LABELS = {"SPEI_1": "SPEI-1", "SPEI_3": "SPEI-3", "SPEI_48": "SPEI-48"}

THRESHOLD_MARKER_MAP = {"5%": "o", "10%": "s", "20%": "^"}

REF_GRAY = "#9A9A9A"   # vertical dashed
ZERO_GRAY = "#707070"  # zero line

def get_threshold_marker(threshold_pct):
    threshold_label = str(threshold_pct)
    if threshold_label in THRESHOLD_MARKER_MAP:
        return THRESHOLD_MARKER_MAP[threshold_label]
    try:
        pct_val = int(float(threshold_label))
        label_map = {5: "5%", 10: "10%", 20: "20%"}
        return THRESHOLD_MARKER_MAP.get(label_map.get(pct_val, "10%"), "s")
    except:
        return "s"

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data for Panel A...")
coast_prediction_impact = pd.read_csv(
    os.path.join(output_dir, "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv")
)
coast_threshold_markers = pd.read_csv(
    os.path.join(output_dir, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv")
)

if 'threshold_pct' in coast_threshold_markers.columns:
    if coast_threshold_markers['threshold_pct'].dtype in ['int64', 'float64']:
        coast_threshold_markers['threshold_pct'] = coast_threshold_markers['threshold_pct'].astype(int).astype(str) + '%'

# --------------------------------------------------------------------------
# Diagnostic checks for coast names
# --------------------------------------------------------------------------
print("\nUnique coast names in prediction file:")
print(sorted(coast_prediction_impact["coast_region"].dropna().unique()))

print("\nUnique coast names in threshold marker file:")
print(sorted(coast_threshold_markers["coast_region"].dropna().unique()))

# Validate that all coasts in COAST_REGION_LEVELS have colors defined
missing_colors = [c for c in COAST_REGION_LEVELS if c not in COAST_COLORS]
if missing_colors:
    raise ValueError(f"Missing colors for coasts: {missing_colors}")

# Check if any coast in COAST_REGION_LEVELS is missing from the data
prediction_coasts = set(coast_prediction_impact["coast_region"].dropna().unique())
marker_coasts = set(coast_threshold_markers["coast_region"].dropna().unique())
missing_prediction = [c for c in COAST_REGION_LEVELS if c not in prediction_coasts]
if missing_prediction:
    print(f"WARNING: These coasts are missing from prediction data: {missing_prediction}")
missing_marker = [c for c in COAST_REGION_LEVELS if c not in marker_coasts]
if missing_marker:
    print(f"WARNING: These coasts are missing from threshold marker data: {missing_marker}")

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
    # ----- ensure tick marks are visible -----
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
panel_labels = [chr(97 + i) + ')' for i in range(12)]

marker_handles = []
marker_labels = []

for idx_c, coast in enumerate(COAST_REGION_LEVELS):
    for idx_t, timescale in enumerate(SELECTED_TIMESCALES):
        ax = axes[idx_c, idx_t]
        panel_idx = idx_c * 3 + idx_t

        ax.grid(False)
        ax.minorticks_off()
        ax.set_facecolor("white")

        # ------------------------------------------------------------------
        # 1. CONFIDENCE RIBBON
        # ------------------------------------------------------------------
        subset = coast_prediction_impact[
            (coast_prediction_impact['coast_region'] == coast) &
            (coast_prediction_impact['SPEI_timescale'] == timescale)
        ]

        if len(subset) > 0:
            ax.fill_between(
                subset['SPEI_value'],
                subset['predicted_lower_pct'],
                subset['predicted_upper_pct'],
                color=COAST_COLORS[coast],
                alpha=0.10,
                zorder=1
            )

        # ------------------------------------------------------------------
        # 2. REFERENCE LINES (SPEI anomalies and zero)
        # ------------------------------------------------------------------
        # Vertical lines at SPEI = -1 and +1
        ax.axvline(x=-1, color=REF_GRAY, linestyle='--', linewidth=1.0, alpha=0.85, zorder=2)
        ax.axvline(x=1,  color=REF_GRAY, linestyle='--', linewidth=1.0, alpha=0.85, zorder=2)
        # Zero line (stronger)
        ax.axhline(y=0, color=ZERO_GRAY, linestyle='-', linewidth=1.1, alpha=0.9, zorder=2)

        # ------------------------------------------------------------------
        # 2b. WUE_T threshold reference lines, matching Malone logic
        #     Fixed horizontal dotted lines at ±5%, ±10%, ±20%
        #     Updated to be darker: color="#595959", alpha=0.75, lw=0.9
        # ------------------------------------------------------------------
        for yref in [-20, -10, -5, 5, 10, 20]:
            ax.axhline(
                y=yref,
                color="#595959",   # close to ggplot gray35
                linestyle="--",
                linewidth=0.9,
                alpha=0.75,
                zorder=2.4
            )

        # ------------------------------------------------------------------
        # 3. GAM CURVE
        # ------------------------------------------------------------------
        if len(subset) > 0:
            ax.plot(
                subset['SPEI_value'],
                subset['predicted_pct_change'],
                color=COAST_COLORS[coast],
                linewidth=1.8,
                label=coast,
                zorder=3
            )

            # ------------------------------------------------------------------
            # 4. THRESHOLD MARKERS ONLY
            # ------------------------------------------------------------------
            markers_subset = coast_threshold_markers[
                (coast_threshold_markers['coast_region'] == coast) &
                (coast_threshold_markers['SPEI_timescale'] == timescale)
            ]

            for _, row in markers_subset.iterrows():
                marker = get_threshold_marker(row['threshold_pct'])

                ax.scatter(
                    x=row['SPEI_threshold'],
                    y=row['pct_change_threshold'],
                    color=COAST_COLORS[coast],
                    s=150,
                    marker=marker,
                    edgecolors='white',
                    linewidths=0.7,
                    zorder=6
                )

                if row['threshold_pct'] not in marker_labels:
                    marker_labels.append(row['threshold_pct'])
                    marker_handles.append(
                        plt.Line2D(
                            [0], [0],
                            marker=marker,
                            color='w',
                            markerfacecolor='gray',
                            markersize=14,
                            label=row['threshold_pct']
                        )
                    )

        # ------------------------------------------------------------------
        # 5. AXIS LABELS, PANEL LABELS, TICKS
        # ------------------------------------------------------------------
        ax.set_title('')
        ax.text(0.02, 0.95, panel_labels[panel_idx],
                transform=ax.transAxes, fontsize=tick_font+4,
                fontweight='bold', va='top', ha='left')

        # Use the display label for the row title
        if idx_t == 0:
            ax.set_ylabel(COAST_DISPLAY_LABELS[coast], fontweight='bold', fontsize=base_font)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.xaxis.set_major_locator(plt.MaxNLocator(5, integer=True))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5, integer=True))

        # ----- make x and y tick marks visible -----
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
            pad=6          # extra space between tick marks and labels
        )

        # make sure axis spines are visible and black
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.5)

        if idx_c == 0:
            ax.set_title(TIMESCALE_LABELS[timescale], fontweight='bold', fontsize=base_font+4)

# ----------------------------------------------------------------------------
# LAYOUT AND LEGEND
# ----------------------------------------------------------------------------
plt.subplots_adjust(left=0.12, right=0.78, top=0.95, bottom=0.08, wspace=0.3, hspace=0.4)

# ---- ADD GLOBAL Y‑AXIS LABEL ----
fig.text(
    0.01, 0.5,
    r'Upland WUE$_{T}$ Change (%)',
    ha='center',
    va='center',
    rotation=90,
    fontsize=base_font
)

if marker_handles:
    fig.legend(
        handles=marker_handles,
        labels=marker_labels,
        title='Threshold Change (%)',
        loc='center left',
        bbox_to_anchor=(0.82, 0.5),
        framealpha=0.9,
        edgecolor='black',
        fontsize=legend_font,
        title_fontsize=legend_font
    )

plt.suptitle('')

# ----------------------------------------------------------------------------
# SAVE AND SHOW
# ----------------------------------------------------------------------------
output_path = os.path.join(figure_dir, "Q3_panel_A_coast_GAM_thresholds.png")
fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.show()
print(f"Panel A saved to: {output_path}")