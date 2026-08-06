# -*- coding: utf-8 -*-
"""
Q3_threshold_summary_spei_gradient_dry_only.py

Dry‑side only summary of threshold crossings for significant coast × SPEI combinations.
- Only dry‑side thresholds (SPEI < -1).
- Only significant coast-SPEI pairs.
- Y-axis: SPEI-X ↑ or ↓ (no "D" since all are dry).
- Ecosystem represented by color (no black outline).
- Threshold percentage represented by marker shape.
- Top two highest thresholds per group retained.
- PNG only.
- Legends: Ecosystem (lower right, top), Threshold (lower right, below).
- No title or subtitle.
- Large font sizes for publication.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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

input_file = os.path.join(input_dir,
    "Q3_reviewer_threshold_markers_month_averaged_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv")

# ============================================================================
# CONSTANTS
# ============================================================================
ECOSYSTEM_COLORS = {
    'Upland':      '#800080',
    'Freshwater':  '#0000FF',
    'Saline':      '#FFA500'
}
ECOSYSTEM_ORDER = ['Upland', 'Freshwater', 'Saline']

SIGNIFICANT_COAST_SPEI = [
    ("AK Coast", "SPEI_3"),
    ("AK Coast", "SPEI_48"),
    ("Pacific Coast", "SPEI_1"),
    ("Pacific Coast", "SPEI_3"),
    ("Atlantic Coast", "SPEI_3"),
    ("Atlantic Coast", "SPEI_48"),
]

COAST_ORDER = ["AK Coast", "Pacific Coast", "Atlantic Coast"]
COAST_DISPLAY = {
    "AK Coast":       "Alaska Coast",
    "Pacific Coast":  "Pacific Coast",
    "Atlantic Coast": "Atlantic Coast"
}

THRESHOLD_LEVELS = [5, 10, 20, 30, 35]
THRESHOLD_MARKER_MAP = {
    5:  "o",
    10: "s",
    20: "^",
    30: "P",
    35: "X"
}

DODGE = {
    "Upland": -0.18,
    "Freshwater": 0.00,
    "Saline": 0.18
}

# ============================================================================
# LOAD AND FILTER DATA (DRY-ONLY)
# ============================================================================
print("Loading threshold data...")
df = pd.read_csv(input_file)

# Ensure threshold_pct is numeric
if df['threshold_pct'].dtype == 'object':
    df['threshold_pct_num'] = df['threshold_pct'].str.replace('%', '').astype(int)
else:
    df['threshold_pct_num'] = df['threshold_pct']

df = df[df['threshold_pct_num'].isin(THRESHOLD_LEVELS)]

# Clean side/direction
df['anomaly_side_clean'] = df['anomaly_side'].astype(str).str.lower()
df['impact_direction_clean'] = df['impact_direction'].astype(str).str.lower()

# Dry-only
df = df[
    (df['anomaly_side_clean'] == 'dry') &
    (df['SPEI_threshold'] < -1)
]

# Filter to significant coast-SPEI pairs
df['coast_spei'] = list(zip(df['coast_region'], df['SPEI_timescale']))
df = df[df['coast_spei'].isin(SIGNIFICANT_COAST_SPEI)]
df = df[df['coast_region'].isin(COAST_ORDER)]

print(f"Data rows after dry-side filters: {len(df)}")
print("Coasts:", sorted(df['coast_region'].unique()))
print("Timescales:", sorted(df['SPEI_timescale'].unique()))

if df.empty:
    raise ValueError("No dry-side data remaining after filtering.")

# ---- Keep top two highest thresholds per group ----
group_cols = [
    'coast_region',
    'water_class',
    'SPEI_timescale',
    'impact_direction_clean'
]
df = df.sort_values(group_cols + ['threshold_pct_num'], ascending=[True, True, True, True, False])
df = df.groupby(group_cols, observed=True).head(2).copy()
print(f"Rows after top-2 per group: {len(df)}")

# ---- Create concise y-label: SPEI-X ↑ or ↓ (no "D") ----
df['dir_symbol'] = df['impact_direction_clean'].map({'increase': '↑', 'decrease': '↓'})
df['SPEI_label'] = df['SPEI_timescale'].astype(str).str.replace('SPEI_', 'SPEI-', regex=False)
df['y_label'] = df['SPEI_label'] + ' ' + df['dir_symbol']  # e.g., "SPEI-1 ↑"

# ---- Define global y-label order ----
timescale_order = ['SPEI_1', 'SPEI_3', 'SPEI_48']
dir_order = ['↓', '↑']  # decrease before increase

global_y_labels = []
for ts in timescale_order:
    if ts not in df['SPEI_timescale'].unique():
        continue
    for direction in dir_order:
        label = ts.replace('SPEI_', 'SPEI-') + ' ' + direction
        if label in df['y_label'].unique():
            global_y_labels.append(label)

global_y_pos = {label: i for i, label in enumerate(global_y_labels)}

# ---- Prepare per-coast data ----
coast_panels = {}
for coast in COAST_ORDER:
    sub = df[df['coast_region'] == coast].copy()
    if sub.empty:
        coast_panels[coast] = None
        continue
    coast_labels = [lbl for lbl in global_y_labels if lbl in sub['y_label'].unique()]
    sub['y_pos'] = sub['y_label'].map(global_y_pos)
    sub['y_pos'] = pd.to_numeric(sub['y_pos'], errors='coerce')
    sub = sub.dropna(subset=['y_pos'])
    coast_panels[coast] = {'data': sub, 'y_labels': coast_labels}

# ============================================================================
# PLOTTING – SINGLE PANEL PER COAST (STACKED)
# ============================================================================
sns.set_style("white")
# ---- Large fonts for publication ----
plt.rcParams.update({
    'font.size': 26,
    'axes.labelsize': 28,
    'axes.titlesize': 30,
    'xtick.labelsize': 26,
    'ytick.labelsize': 26,
    'legend.fontsize': 26,
    'legend.title_fontsize': 28,
    'figure.titlesize': 32,
})

n_coasts = len(COAST_ORDER)
fig, axes = plt.subplots(n_coasts, 1, figsize=(16, 4 * n_coasts), sharex=True)
if n_coasts == 1:
    axes = [axes]

# ---- Adjust subplot to make room for right-side legends ----
plt.subplots_adjust(right=0.78, left=0.12, top=0.95, bottom=0.08)

color_map = ECOSYSTEM_COLORS
marker_map = THRESHOLD_MARKER_MAP

for idx, coast in enumerate(COAST_ORDER):
    ax = axes[idx]
    panel = coast_panels.get(coast)

    if panel is None or panel['data'].empty:
        ax.text(0.5, 0.5, f"No dry thresholds for {COAST_DISPLAY[coast]}",
                transform=ax.transAxes, ha='center', va='center', fontsize=26, color='gray')
        ax.set_ylabel(COAST_DISPLAY[coast], fontweight='bold', fontsize=28)
        ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        continue

    data = panel['data'].copy()
    coast_labels = panel['y_labels']

    # Compute dodge
    data['y_dodge'] = data['water_class'].astype(str).map(DODGE)
    data['y_dodge'] = pd.to_numeric(data['y_dodge'], errors='coerce')
    data['y_plot'] = data['y_pos'] + data['y_dodge']

    # Plot points – NO BLACK OUTLINE
    for _, row in data.iterrows():
        x = row['SPEI_threshold']
        y = row['y_plot']
        ecosystem = row['water_class']
        color = color_map[ecosystem]
        thr = row['threshold_pct_num']
        marker = marker_map.get(thr, 'o')

        ax.scatter(x, y, s=120, color=color, marker=marker,
                   edgecolors='none', linewidth=0, zorder=3)

    # ---- Set y-ticks only for labels that exist in this coast ----
    y_positions = [global_y_pos[lbl] for lbl in coast_labels]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(coast_labels, fontsize=26)
    ax.set_ylabel(COAST_DISPLAY[coast], fontweight='bold', fontsize=28)

    if y_positions:
        ax.set_ylim(min(y_positions) - 0.8, max(y_positions) + 0.8)
    else:
        ax.set_ylim(-0.5, 0.5)

    # X-axis range (dry side only)
    ax.set_xlim(-5.0, 0.2)
    ax.invert_yaxis()

    # Reference lines
    ax.axvline(x=-1, color='gray', linestyle='--', linewidth=1.5, alpha=0.6, zorder=1)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=2, alpha=0.8, zorder=1)
    ax.axvspan(-1, 0, facecolor='lightgray', alpha=0.15, zorder=0)

    # Horizontal grid lines
    for pos in y_positions:
        ax.axhline(y=pos - 0.5, color='lightgray', linestyle='-', linewidth=0.8, alpha=0.3, zorder=0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # X ticks at 0.5 intervals
    ax.set_xticks(np.arange(-5, 0.5, 0.5))

# X-axis label for bottom panel
axes[-1].set_xlabel('SPEI value at threshold crossing', fontsize=30, fontweight='bold')

# ---- NO TITLE OR SUBTITLE ----

# ---- LEGENDS placed lower right, stacked vertically ----
# Ecosystem (top of the two)
eco_handles = [mpatches.Patch(color=color_map[eco], label=eco) for eco in ECOSYSTEM_ORDER]
leg_eco = fig.legend(handles=eco_handles, title='Ecosystem',
                     loc='center left', bbox_to_anchor=(0.95, 0.75),
                     frameon=True, edgecolor='black', fontsize=26, title_fontsize=28)

# Threshold (below ecosystem)
marker_handles = []
for thr in sorted(marker_map.keys()):
    if thr in df['threshold_pct_num'].unique():
        marker_handles.append(
            mlines.Line2D([], [], color='gray', marker=marker_map[thr], linestyle='None',
                          markersize=14, label=f"{thr}%")
        )
if marker_handles:
    leg_marker = fig.legend(handles=marker_handles, title='Threshold',
                            loc='center left', bbox_to_anchor=(0.95, 0.45),
                            frameon=True, edgecolor='black', fontsize=26, title_fontsize=28)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# ============================================================================
# SAVE PNG ONLY
# ============================================================================
png_path = os.path.join(figure_dir, "Q3_threshold_summary_spei_gradient_dry_only.png")
plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white')
print(f"PNG saved to: {png_path}")

plt.show()