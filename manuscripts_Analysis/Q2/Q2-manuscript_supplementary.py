# -*- coding: utf-8 -*-
"""
Q2_supplementary_figure.py – 2×2 supplementary figure for Q2 site‑level
WUE_T performance metrics by ecosystem class.

Panels:
a) Plasticity slope (SPEI-3)
b) Plasticity range
c) Drought resistance
d) Post-drought recovery

Boxplots have black outlines. Panel a p‑value moved slightly higher.
Saves to Q2_WUE_performance_figures/Q2_supplementary_figure.png
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2"
input_dir = os.path.join(base_dir, "Q2_WUE_performance_outputs")
output_dir = os.path.join(base_dir, "Q2_WUE_performance_figures")
os.makedirs(output_dir, exist_ok=True)   # create if not exists

slope_file = os.path.join(input_dir, "Q2_site_plasticity_slope.csv")
range_file = os.path.join(input_dir, "Q2_site_plasticity_range.csv")
resistance_file = os.path.join(input_dir, "Q2_site_resistance.csv")
recovery_file = os.path.join(input_dir, "Q2_site_recovery.csv")
stats_file = os.path.join(input_dir, "Q2_ecosystem_class_statistical_tests.csv")

# ============================================================================
# ECOSYSTEM CONFIGURATION
# ============================================================================
ECOSYSTEM_ORDER = ["Upland", "Freshwater", "Saline"]
ECOSYSTEM_COLORS = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500"
}

# ============================================================================
# FONT SIZES
# ============================================================================
TITLE_SIZE = 26
AXIS_LABEL_SIZE = 26
TICK_LABEL_SIZE = 22
P_VALUE_SIZE = 24

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def standardize_ecosystem_column(df):
    df = df.copy()
    if "ecosystem_class" in df.columns:
        eco_col = "ecosystem_class"
    elif "water_class" in df.columns:
        eco_col = "water_class"
    else:
        raise KeyError(f"No ecosystem column. Columns: {list(df.columns)}")
    df["ecosystem_class"] = df[eco_col].astype(str).str.strip()
    df = df[df["ecosystem_class"].isin(ECOSYSTEM_ORDER)].copy()
    df["ecosystem_class"] = pd.Categorical(df["ecosystem_class"],
                                            categories=ECOSYSTEM_ORDER,
                                            ordered=True)
    return df

def get_n_per_ecosystem(df, eco_col="ecosystem_class"):
    return df.groupby(eco_col, observed=True).size().to_dict()

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading Q2 performance data...")

# Panel A: slope (SPEI-3 only)
df_slope = pd.read_csv(slope_file)
if 'SPEI_timescale' in df_slope.columns:
    df_slope = df_slope[df_slope['SPEI_timescale'] == 'SPEI_3'].copy()
df_slope = standardize_ecosystem_column(df_slope)
slope_col = 'slope' if 'slope' in df_slope.columns else 'plasticity_slope'

# Panel B: range – use plasticity_p95_p05
df_range = pd.read_csv(range_file)
df_range = standardize_ecosystem_column(df_range)
if 'plasticity_p95_p05' in df_range.columns:
    range_col = 'plasticity_p95_p05'
elif 'plasticity_range' in df_range.columns:
    range_col = 'plasticity_range'
else:
    raise KeyError(f"No plasticity range column found. Columns: {list(df_range.columns)}")

# Panel C: resistance
df_resist = pd.read_csv(resistance_file)
df_resist = standardize_ecosystem_column(df_resist)
resist_col = 'resistance'

# Panel D: recovery
df_recov = pd.read_csv(recovery_file)
df_recov = standardize_ecosystem_column(df_recov)
recov_col = 'mean_recovery' if 'mean_recovery' in df_recov.columns else 'recovery'

# Statistical p‑values
df_stats = pd.read_csv(stats_file)
metric_to_response = {
    "Plasticity slope (SPEI-3)": "plasticity_slope_SPEI3",
    "Plasticity range": "plasticity_p95_p05",
    "Drought resistance": "resistance",
    "Post-drought recovery": "mean_recovery"
}
p_value_dict = {}
for label, resp_var in metric_to_response.items():
    if "response_variable" in df_stats.columns:
        mask = df_stats["response_variable"].astype(str).str.lower() == resp_var.lower()
    elif "metric" in df_stats.columns:
        mask = df_stats["metric"].astype(str).str.lower() == resp_var.lower()
    else:
        mask = pd.Series(False, index=df_stats.index)
    if mask.any() and "p_value" in df_stats.columns:
        p_value_dict[label] = df_stats.loc[mask, "p_value"].values[0]
    else:
        p_value_dict[label] = np.nan

# ============================================================================
# SAMPLE SIZES
# ============================================================================
n_slope = get_n_per_ecosystem(df_slope)
n_range = get_n_per_ecosystem(df_range)
n_resist = get_n_per_ecosystem(df_resist)
n_recov = get_n_per_ecosystem(df_recov)

# ============================================================================
# CREATE FIGURE (2×2)
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(17, 13))
ax = axes.flatten()

# Panel letters only
panel_letters = ['a)', 'b)', 'c)', 'd)']

# Y‑axis labels – concise
ylabels = [
    r"Plasticity slope (WUE$_T$ ~ SPEI-3)",
    r"Plasticity range",
    r"Drought resistance",
    r"Post-drought recovery"
]

# Data: (DataFrame, column, n_dict, reference y, p-position)
panel_data = [
    (df_slope, slope_col, n_slope, 0, 'bottomright'),
    (df_range, range_col, n_range, None, 'topright'),
    (df_resist, resist_col, n_resist, 1, 'topright'),
    (df_recov, recov_col, n_recov, 1, 'topright')
]

# ============================================================================
# PLOTTING LOOP
# ============================================================================
for i, (df, val_col, n_dict, ref_y, p_pos) in enumerate(panel_data):
    ax_i = ax[i]
    
    df_plot = df[['ecosystem_class', val_col]].dropna()
    df_plot = df_plot.sort_values('ecosystem_class')
    
    data_n = df_plot.groupby('ecosystem_class', observed=True).size().reindex(ECOSYSTEM_ORDER).fillna(0).astype(int)
    x_labels = [f"{eco}\n(N={data_n[eco]})" for eco in ECOSYSTEM_ORDER]
    
    # ---- Styles ----
    if i == 0:  # a) violin + box + jitter
        sns.violinplot(data=df_plot, x='ecosystem_class', y=val_col,
                       order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                       inner=None, cut=0, linewidth=0, alpha=0.3, ax=ax_i)
        sns.boxplot(data=df_plot, x='ecosystem_class', y=val_col,
                    order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                    width=0.2, boxprops={'alpha':0.7, 'edgecolor':'black', 'linewidth':1.5},
                    whiskerprops={'color':'black', 'linewidth':1.5},
                    capprops={'color':'black', 'linewidth':1.5},
                    medianprops={'color':'black', 'linewidth':1.5},
                    fliersize=0, ax=ax_i)
        sns.stripplot(data=df_plot, x='ecosystem_class', y=val_col,
                      order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                      size=8, jitter=0.2, alpha=0.8, ax=ax_i)
    elif i == 1:  # b) box + jitter
        sns.boxplot(data=df_plot, x='ecosystem_class', y=val_col,
                    order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                    width=0.5, boxprops={'alpha':0.7, 'edgecolor':'black', 'linewidth':1.5},
                    whiskerprops={'color':'black', 'linewidth':1.5},
                    capprops={'color':'black', 'linewidth':1.5},
                    medianprops={'color':'black', 'linewidth':1.5},
                    fliersize=0, ax=ax_i)
        sns.stripplot(data=df_plot, x='ecosystem_class', y=val_col,
                      order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                      size=8, jitter=0.25, alpha=0.8, ax=ax_i)
    elif i == 2:  # c) violin + box + jitter
        sns.violinplot(data=df_plot, x='ecosystem_class', y=val_col,
                       order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                       inner=None, cut=0, linewidth=0, alpha=0.3, ax=ax_i)
        sns.boxplot(data=df_plot, x='ecosystem_class', y=val_col,
                    order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                    width=0.2, boxprops={'alpha':0.7, 'edgecolor':'black', 'linewidth':1.5},
                    whiskerprops={'color':'black', 'linewidth':1.5},
                    capprops={'color':'black', 'linewidth':1.5},
                    medianprops={'color':'black', 'linewidth':1.5},
                    fliersize=0, ax=ax_i)
        sns.stripplot(data=df_plot, x='ecosystem_class', y=val_col,
                      order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                      size=8, jitter=0.2, alpha=0.8, ax=ax_i)
    else:  # d) strip + box (point‑heavy)
        sns.stripplot(data=df_plot, x='ecosystem_class', y=val_col,
                      order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                      size=9, jitter=0.15, alpha=0.8, ax=ax_i)
        sns.boxplot(data=df_plot, x='ecosystem_class', y=val_col,
                    order=ECOSYSTEM_ORDER, palette=ECOSYSTEM_COLORS,
                    width=0.3, boxprops={'alpha':0.5, 'edgecolor':'black', 'linewidth':1.5, 'zorder':1},
                    whiskerprops={'color':'black', 'linewidth':1.5, 'zorder':1},
                    capprops={'color':'black', 'linewidth':1.5, 'zorder':1},
                    medianprops={'color':'black', 'linewidth':1.5, 'zorder':1},
                    fliersize=0, ax=ax_i)
    
    # ---- Reference line ----
    if ref_y is not None:
        ax_i.axhline(ref_y, color='gray', linestyle='--', linewidth=2, alpha=0.7)
    
    # ---- X‑axis labels ----
    ax_i.set_xticklabels(x_labels, rotation=0, fontsize=TICK_LABEL_SIZE)
    ax_i.set_xlabel('')
    
    # ---- Y‑axis label ----
    ax_i.set_ylabel(ylabels[i], fontsize=AXIS_LABEL_SIZE, labelpad=12)
    
    # ---- Panel letter only ----
    ax_i.set_title(panel_letters[i], fontsize=TITLE_SIZE, fontweight='bold', loc='left', pad=12)
    
    # ---- Tick parameters ----
    ax_i.tick_params(axis='both', labelsize=TICK_LABEL_SIZE, width=1.8, length=6)
    
    # ---- Light grid ----
    ax_i.grid(axis='y', linestyle=':', alpha=0.4)
    ax_i.set_axisbelow(True)
    
    # ---- Spines ----
    sns.despine(ax=ax_i, top=True, right=True)
    for spine in ax_i.spines.values():
        spine.set_linewidth(1.8)
        spine.set_color('black')
    
    # ---- p‑value annotation ----
    key_map = ["Plasticity slope (SPEI-3)", "Plasticity range",
               "Drought resistance", "Post-drought recovery"]
    p_val = p_value_dict.get(key_map[i], np.nan)
    if not np.isnan(p_val):
        p_text = f"p = {p_val:.2f}" if p_val >= 0.01 else "p < 0.01"
        if p_pos == 'bottomright':
            if i == 0:
                x_pos, y_pos = 0.98, 0.08   # panel a moved up
            else:
                x_pos, y_pos = 0.98, 0.02
            va, ha = 'bottom', 'right'
        else:
            x_pos, y_pos = 0.98, 0.95
            va, ha = 'top', 'right'
        ax_i.text(x_pos, y_pos, p_text, transform=ax_i.transAxes,
                  fontsize=P_VALUE_SIZE, fontweight='bold', va=va, ha=ha,
                  bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

# ============================================================================
# FINAL LAYOUT AND SAVE
# ============================================================================
plt.subplots_adjust(hspace=0.4, wspace=0.25, top=0.95, bottom=0.08)
fig.canvas.draw()
plt.show(block=True)

# ---- Save to the specified output directory ----
output_file = os.path.join(output_dir, "Q2_supplementary_figure.png")
fig.savefig(output_file, dpi=300, facecolor='white', bbox_inches='tight')
print(f"Figure saved to: {output_file}")