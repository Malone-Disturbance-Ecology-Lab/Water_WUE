# -*- coding: utf-8 -*-
"""
Q3_ecosystem_GAM_curves_supplementary.py (FINAL)

Supplementary figure showing ecosystem-specific GAM smooth curves only.
Layout: 3 rows (Upland, Freshwater, Saline) × 4 columns (SPEI-1, -3, -24, -48).

- X-axis ticks: all integers from -5 to 5 (filtered by data range).
- Figure size increased to 54 x 36 inches.
- All fonts enlarged significantly.
- Section labels (a, b, c) removed.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# -----------------------------------------------------------------------------
# PATHS AND CONSTANTS
# -----------------------------------------------------------------------------
input_csv = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_august_update\Q3_ecosystem_GAM_predictions_all_ecosystems_all_months.csv"
output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_figures\talib_manuscript"
os.makedirs(output_dir, exist_ok=True)

# Ecosystem colours
ecosystem_colors = {
    "Upland":    "#800080",   # purple
    "Freshwater": "#0000FF",  # blue
    "Saline":    "#FFA500"    # orange
}

ecosystem_order = ["Upland", "Freshwater", "Saline"]
# No section labels needed

spei_order = ["SPEI_1", "SPEI_3", "SPEI_24", "SPEI_48"]
spei_labels = ["SPEI-1", "SPEI-3", "SPEI-24", "SPEI-48"]

# -----------------------------------------------------------------------------
# LOAD AND AGGREGATE DATA (average over months)
# -----------------------------------------------------------------------------
df = pd.read_csv(input_csv)
df["month_f"] = df["month_f"].astype(int)
df = df[df["SPEI_timescale"].isin(spei_order)]

agg_df = df.groupby(
    ["water_class", "SPEI_timescale", "SPEI_value"],
    observed=True
).agg(
    mean_pct=("predicted_pct_change", "mean"),
    lower_pct=("predicted_lower_pct", "mean"),
    upper_pct=("predicted_upper_pct", "mean"),
    n_months=("month_f", "nunique")
).reset_index()

# -----------------------------------------------------------------------------
# COMMON Y‑AXIS LIMIT
# -----------------------------------------------------------------------------
max_abs = 0
for eco in ecosystem_order:
    eco_data = agg_df[agg_df["water_class"] == eco]
    for spei in spei_order:
        sub = eco_data[eco_data["SPEI_timescale"] == spei]
        if not sub.empty:
            vals = sub[["mean_pct", "lower_pct", "upper_pct"]].to_numpy()
            max_abs = max(max_abs, np.nanmax(np.abs(vals)))

def nice_y_limit(x, base=10, min_lim=20):
    if not np.isfinite(x) or x <= 0:
        return min_lim
    return max(min_lim, int(np.ceil(x / base) * base))

ylim = nice_y_limit(max_abs, base=10, min_lim=20)

# -----------------------------------------------------------------------------
# X‑AXIS: full integer ticks from -5 to 5 (filtered by data range)
# -----------------------------------------------------------------------------
x_min = np.floor(agg_df["SPEI_value"].min())
x_max = np.ceil(agg_df["SPEI_value"].max())
all_integer_ticks = list(range(-5, 6))
x_ticks = [t for t in all_integer_ticks if x_min <= t <= x_max]

# -----------------------------------------------------------------------------
# CREATE FIGURE – FURTHER ENLARGED
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(
    nrows=3, ncols=4,
    figsize=(54, 36),          # increased from 48x32
    sharex=False,
    sharey=True,
    gridspec_kw={'wspace': 0.30, 'hspace': 0.40}
)

# ---- Greatly increased font sizes ----
tick_fontsize      = 50        # was 42
title_fontsize     = 56        # was 48
axis_label_fontsize = 60       # was 50
line_width         = 7.0       # thicker curves
ribbon_alpha       = 0.15

# ---- Plot each panel ----
for row_idx, eco in enumerate(ecosystem_order):
    for col_idx, spei in enumerate(spei_order):
        ax = axes[row_idx, col_idx]
        sub = agg_df[
            (agg_df["water_class"] == eco) &
            (agg_df["SPEI_timescale"] == spei)
        ]
        if not sub.empty:
            sub = sub.sort_values("SPEI_value")
            x = sub["SPEI_value"].values
            y = sub["mean_pct"].values
            lower = sub["lower_pct"].values
            upper = sub["upper_pct"].values
            color = ecosystem_colors[eco]

            ax.fill_between(
                x, lower, upper,
                color=color,
                alpha=ribbon_alpha,
                linewidth=0,
                zorder=1
            )
            ax.plot(
                x, y,
                color=color,
                linewidth=line_width,
                zorder=3
            )

        # Reference lines
        ax.axhline(0, color="black", linestyle="-", linewidth=3.0, alpha=0.6)
        ax.axvline(-1, color="gray", linestyle="--", linewidth=3.0, alpha=0.6)
        ax.axvline(1, color="gray", linestyle="--", linewidth=3.0, alpha=0.6)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-ylim, ylim)
        ax.set_xticks(x_ticks)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

        ax.grid(True, linestyle=":", color="lightgray", alpha=0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for spine in ax.spines.values():
            spine.set_linewidth(5.0)               # even thicker spines

        ax.tick_params(
            axis='both',
            which='major',
            direction='out',
            length=24,           # longer ticks
            width=4.0,
            colors='black',
            bottom=True,
            left=True,
            top=False,
            right=False,
            pad=16
        )

        # Column titles (top row)
        if row_idx == 0:
            ax.set_title(spei_labels[col_idx], fontsize=title_fontsize, pad=20)

        # Row labels (left column) – ecosystem names
        if col_idx == 0:
            ax.set_ylabel(
                eco,
                fontsize=title_fontsize,
                labelpad=20,
                rotation=90,
                color=ecosystem_colors[eco]
            )
            ax.tick_params(axis="y", labelsize=tick_fontsize, pad=14)
        else:
            ax.tick_params(axis="y", labelleft=False)

        # X‑tick labels only on bottom row
        if row_idx == 2:
            ax.tick_params(axis="x", labelsize=tick_fontsize, pad=14, rotation=0)
        else:
            ax.tick_params(axis="x", labelbottom=False)

# -----------------------------------------------------------------------------
# GLOBAL AXIS LABELS
# -----------------------------------------------------------------------------
fig.text(
    0.50, 0.035,
    "SPEI",
    ha="center",
    fontsize=axis_label_fontsize
)

fig.text(
    0.015, 0.50,
    "Predicted WUE$_T$ change from near‑normal (%)",
    va="center",
    rotation=90,
    fontsize=axis_label_fontsize
)

# -----------------------------------------------------------------------------
# SECTION LABELS (a, b, c) – REMOVED as requested
# -----------------------------------------------------------------------------
# (No code for section labels)

# -----------------------------------------------------------------------------
# SAVE FIGURE
# -----------------------------------------------------------------------------
output_file = os.path.join(output_dir, "Q3_ecosystem_GAM_curves_supplementary.png")
fig.savefig(output_file, dpi=500, bbox_inches="tight")
print(f"Ecosystem GAM curves supplementary figure saved to: {output_file}")

plt.show()
print("\nDone.")