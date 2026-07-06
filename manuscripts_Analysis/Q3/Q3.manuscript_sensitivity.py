# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 14:11:39 2026

@author: ammar
"""

"""
generate_figure_04_only.py

Creates ONLY the multipanel composite Figure 4 from the Q3 WUE analysis,
matching the exact appearance of the original workflow.

Assumes all required CSV files are present in:
    M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs

Saves output to:
    M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_figures\talib_manuscript\
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS (match original)
# ============================================================================
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures", "talib_manuscript")

os.makedirs(figure_dir, exist_ok=True)

# ============================================================================
# CONSTANTS (identical to original)
# ============================================================================
SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
ECOSYSTEM_COLORS = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500"
}
COAST_COLORS = {
    "Atlantic Coast": "#008B8B",
    "Pacific Coast": "#C44E52",
    "Gulf Coast": "#8C564B",
    "AK Coast": "#4D4D4D"
}
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]
GRAY_COLOR = "#595959"

THRESHOLD_MARKER_MAP = {"5%": "o", "10%": "s", "20%": "^"}

def get_threshold_marker(threshold_pct):
    """Get marker symbol for threshold percentage label (same as original)"""
    threshold_label = str(threshold_pct)
    if threshold_label in THRESHOLD_MARKER_MAP:
        return THRESHOLD_MARKER_MAP[threshold_label]
    try:
        pct_val = int(float(threshold_label))
        label_map = {5: "5%", 10: "10%", 20: "20%"}
        return THRESHOLD_MARKER_MAP.get(label_map.get(pct_val, "10%"), "s")
    except:
        return "s"

def theme_wue(base_size=13):
    """Identical to Malone's theme_wue()"""
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
        'legend.title_fontsize': base_size - 1,
        'grid.alpha': 0.3,
    })

# ============================================================================
# LOAD REQUIRED DATA
# ============================================================================
print("Loading required data files...")
files_needed = {
    "coast_prediction_impact": "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv",
    "coast_threshold_markers": "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv",
    "site_sensitivity_rank": "Q3_WUE_T_SPEI_site_sensitivity_rank.csv",
}

data = {}
for name, fname in files_needed.items():
    path = os.path.join(output_dir, fname)
    if os.path.exists(path):
        data[name] = pd.read_csv(path)
        print(f"  Loaded: {fname}")
    else:
        raise FileNotFoundError(f"Required file not found: {path}")

coast_prediction_impact = data["coast_prediction_impact"]
coast_threshold_markers = data["coast_threshold_markers"]
site_sensitivity_rank = data["site_sensitivity_rank"]

# Ensure threshold_pct column is string with % (as in original)
for df in [coast_threshold_markers]:
    if 'threshold_pct' in df.columns:
        if df['threshold_pct'].dtype in ['int64', 'float64']:
            df['threshold_pct'] = df['threshold_pct'].astype(int).astype(str) + '%'

# ============================================================================
# CREATE FIGURE 4 – MULTIPANEL COMPOSITE
# ============================================================================
print("\nCreating Figure 4 composite...")
theme_wue()

# ----------------------------------------------------------------------------
# Panel A: Coast‑specific GAM thresholds (top)
# ----------------------------------------------------------------------------
fig_a, axes_a = plt.subplots(4, 3, figsize=(12, 10))

for idx_c, coast in enumerate(COAST_REGION_LEVELS):
    for idx_t, timescale in enumerate(SELECTED_TIMESCALES):
        ax = axes_a[idx_c, idx_t]

        subset = coast_prediction_impact[
            (coast_prediction_impact['coast_region'] == coast) &
            (coast_prediction_impact['SPEI_timescale'] == timescale)
        ]

        if len(subset) > 0:
            ax.fill_between(
                subset['SPEI_value'],
                subset['predicted_lower_pct'],
                subset['predicted_upper_pct'],
                color=COAST_COLORS.get(coast, '#808080'),
                alpha=0.1
            )
            ax.plot(
                subset['SPEI_value'],
                subset['predicted_pct_change'],
                color=COAST_COLORS.get(coast, '#808080'),
                linewidth=0.9,
                label=coast
            )

            markers_subset = coast_threshold_markers[
                (coast_threshold_markers['coast_region'] == coast) &
                (coast_threshold_markers['SPEI_timescale'] == timescale)
            ]
            for _, row in markers_subset.iterrows():
                marker = get_threshold_marker(row['threshold_pct'])
                ax.hlines(
                    y=row['pct_change_threshold'],
                    xmin=row['SPEI_threshold_lower'],
                    xmax=row['SPEI_threshold_upper'],
                    color=COAST_COLORS.get(coast, '#808080'),
                    linewidth=1.5,
                    alpha=0.7
                )
                ax.scatter(
                    x=row['SPEI_threshold'],
                    y=row['pct_change_threshold'],
                    color=COAST_COLORS.get(coast, '#808080'),
                    s=60,
                    marker=marker,
                    edgecolor='black',
                    linewidth=1,
                    zorder=5
                )

        ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.5)
        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.5)
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
        ax.axhline(y=-20, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
        ax.axhline(y=-10, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
        ax.axhline(y=-5, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
        ax.axhline(y=5, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
        ax.axhline(y=10, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
        ax.axhline(y=20, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)

        if idx_c == 0:
            ax.set_title(timescale, fontweight='bold')
        if idx_t == 0:
            ax.set_ylabel(coast, fontweight='bold')
        if idx_c == 3:
            ax.set_xlabel('SPEI')

plt.suptitle('Coast-Specific GAM Thresholds for Upland WUE_T Response', fontsize=14, fontweight='bold')
plt.tight_layout()

temp_a_path = os.path.join(figure_dir, "temp_panel_A.png")
fig_a.savefig(temp_a_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig_a)

# ----------------------------------------------------------------------------
# Panel B: Coastline grouped sensitivity (bottom) – FLIPPED
# ----------------------------------------------------------------------------
all_system_sensitivity = site_sensitivity_rank[
    site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES) &
    (site_sensitivity_rank['n_months'] >= 6)
]

fig_b, axes_b = plt.subplots(1, 3, figsize=(14, 6))

if len(all_system_sensitivity) > 0:
    for idx, timescale in enumerate(SELECTED_TIMESCALES):
        ax = axes_b[idx]
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
                boxprops=dict(facecolor='lightgray', alpha=0.5),
                whiskerprops=dict(color=GRAY_COLOR),
                capprops=dict(color=GRAY_COLOR),
                medianprops=dict(color='black', linewidth=1.5)
            )

            for ecosystem in ECOSYSTEM_CLASSES:
                eco_subset = subset[subset['water_class'] == ecosystem]
                if len(eco_subset) > 0:
                    for i, coast in enumerate(COAST_REGION_LEVELS):
                        coast_eco_subset = eco_subset[eco_subset['coast_region'] == coast]
                        if len(coast_eco_subset) > 0:
                            y_jitter = np.random.normal(i, 0.08, len(coast_eco_subset))
                            ax.scatter(
                                coast_eco_subset['mean_abs_pct_change'],
                                y_jitter,
                                color=ECOSYSTEM_COLORS[ecosystem],
                                s=60,
                                alpha=0.7,
                                edgecolor='black',
                                linewidth=0.5,
                                label=ecosystem if i == 0 else ""
                            )

            ax.axvline(x=5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
            ax.axvline(x=10, color='gray', linestyle=':', linewidth=1, alpha=0.5)
            ax.axvline(x=20, color='gray', linestyle=':', linewidth=1, alpha=0.5)

            ax.set_yticks(range(len(COAST_REGION_LEVELS)))
            ax.set_yticklabels(COAST_REGION_LEVELS)
            ax.set_title(timescale, fontweight='bold')
            ax.set_xlabel('Mean absolute WUE_T change from near-normal (%)')
            if idx == 0:
                ax.set_ylabel('Coast region')
            ax.grid(True, alpha=0.3, axis='x')

# Add legend (handles from ECOSYSTEM_CLASSES)
handles = [Patch(facecolor=ECOSYSTEM_COLORS[eco], alpha=0.7, label=eco) for eco in ECOSYSTEM_CLASSES]
fig_b.legend(handles=handles, loc='lower center', ncol=3, title='Ecosystem class', fontsize=11)

plt.suptitle('Site Sensitivity Grouped by Coastline', fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

temp_b_path = os.path.join(figure_dir, "temp_panel_B.png")
fig_b.savefig(temp_b_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig_b)

# ----------------------------------------------------------------------------
# Composite using PIL (exactly as original)
# ----------------------------------------------------------------------------
try:
    from PIL import Image, ImageDraw, ImageFont

    if os.path.exists(temp_a_path) and os.path.exists(temp_b_path):
        img_a = Image.open(temp_a_path)
        img_b = Image.open(temp_b_path)

        composite_width = max(img_a.width, img_b.width)
        img_a_resized = img_a.resize((composite_width, img_a.height))
        img_b_resized = img_b.resize((composite_width, img_b.height))

        composite_height = img_a_resized.height + img_b_resized.height + 50
        composite = Image.new('RGB', (composite_width, composite_height), 'white')
        composite.paste(img_a_resized, (0, 0))
        composite.paste(img_b_resized, (0, img_a_resized.height + 20))

        draw = ImageDraw.Draw(composite)
        try:
            font = ImageFont.truetype("arial.ttf", 40)
        except:
            font = ImageFont.load_default()

        draw.text((20, 20), "A", fill='black', font=font)
        draw.text((20, img_a_resized.height + 30), "B", fill='black', font=font)

        # Final save
        final_path = os.path.join(figure_dir, "Q3_plot_04_WUE_T_SPEI_sensitivity_multipanel.png")
        composite.save(final_path, dpi=(300, 300))
        print(f"Composite figure saved to:\n{final_path}")

        # Clean up temp files
        os.remove(temp_a_path)
        os.remove(temp_b_path)
        print("Temporary panel files removed.")

    else:
        raise RuntimeError("Temporary panel files not found after saving.")

except ImportError:
    print("ERROR: PIL (Pillow) is not installed. Cannot create composite figure.")
    print("Please install Pillow: pip install Pillow")
    raise
except Exception as e:
    print(f"Error during composite creation: {e}")
    raise

print("\nDone. Only Figure 4 was generated.")