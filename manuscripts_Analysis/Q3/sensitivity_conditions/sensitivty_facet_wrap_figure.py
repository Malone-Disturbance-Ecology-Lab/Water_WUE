# -*- coding: utf-8 -*-
"""
Q3_coast_threshold_GAM_all_SPEI_facet_curves.py

Creates a 4-row × 7-column figure of coast-threshold GAM response curves
for all coasts and SPEI timescales. Each panel shows predicted WUE_T
percent change vs. SPEI for Upland, Freshwater, and Saline water classes,
aggregated over growing-season months.

Input:  prediction_curves_all_reference_conditions.csv
Output: Q3_coast_threshold_GAM_all_SPEI_facet_curves.png

No new CSV files are saved.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# -----------------------------------------------------------------------------
# PATHS AND CONSTANTS
# -----------------------------------------------------------------------------
input_csv = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs\sensitivity_conditions_results\prediction_curves_all_reference_conditions.csv"
output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_figures\talib_manuscript"
output_file = os.path.join(output_dir, "Q3_coast_threshold_GAM_all_SPEI_facet_curves.png")
os.makedirs(output_dir, exist_ok=True)

# Water class colours (same as Q3 workflow)
ecosystem_palette = {
    "Upland": "#800080",      # purple
    "Freshwater": "#0000FF",  # blue
    "Saline": "#FFA500"       # orange
}
water_class_order = ["Upland", "Freshwater", "Saline"]

# Coast order (rows)
coast_order = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
coast_labels = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "Alaska Coast"]  # display names

# SPEI timescales (columns)
spei_order = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
spei_labels = ["SPEI-1", "SPEI-3", "SPEI-6", "SPEI-12", "SPEI-24", "SPEI-36", "SPEI-48"]

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
print("Loading data...")
df = pd.read_csv(input_csv)
print(f"Input file: {input_csv}")
print(f"Number of rows: {len(df)}")

# Ensure month_f is integer
df["month_f"] = df["month_f"].astype(int)

# Check unique values
print(f"Unique coast regions: {df['coast_region'].unique().tolist()}")
print(f"Unique water classes: {df['water_class'].unique().tolist()}")
print(f"Unique SPEI timescales: {df['SPEI_timescale'].unique().tolist()}")
months = sorted(df["month_f"].unique())
print(f"Months included: {months} (growing-season months)")
print(f"SPEI value range: {df['SPEI_value'].min():.2f} to {df['SPEI_value'].max():.2f}")

# -----------------------------------------------------------------------------
# AGGREGATE OVER MONTHS (mean & sd per coast, water_class, SPEI_timescale, SPEI_value)
# -----------------------------------------------------------------------------
agg_df = df.groupby(
    ["coast_region", "water_class", "SPEI_timescale", "SPEI_value"],
    observed=True
).agg(
    mean_pct=("predicted_pct_change", "mean"),
    sd_pct=("predicted_pct_change", "std"),
    n_months=("month_f", "nunique")   # count unique months, not rows
).reset_index()

# -----------------------------------------------------------------------------
# COMPUTE DYNAMIC Y-AXIS LIMIT
# -----------------------------------------------------------------------------
# Compute absolute extremes across all panels (mean_pct ± sd_pct)
max_abs_val = 0.0
for _, row in agg_df.iterrows():
    lower = row["mean_pct"] - row["sd_pct"]
    upper = row["mean_pct"] + row["sd_pct"]
    max_abs_val = max(max_abs_val, abs(lower), abs(upper))

# Round up to next multiple of 5, with a minimum of 20
y_lim = int(np.ceil(max_abs_val / 5) * 5)
if y_lim < 20:
    y_lim = 20
print(f"Dynamic y-axis limit: ±{y_lim}")

# -----------------------------------------------------------------------------
# CREATE FIGURE AND SUBPLOTS – INCREASED SPACING & LARGER FONTS
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(
    nrows=4, ncols=7,
    figsize=(26, 16),          # bigger figure for more breathing room
    sharex=True,
    sharey=True,
    constrained_layout=False
)

# Increase space between rows and columns, and adjust margins
fig.subplots_adjust(
    left=0.08,       # space for y-axis label and row labels
    right=0.98,
    top=0.94,
    bottom=0.16,     # space for x-axis label, legend, and tick labels
    wspace=0.25,     # wider gap between columns
    hspace=0.25      # wider gap between rows
)

# Set common x range
x_min, x_max = -3.0, 3.0
x_ticks = np.arange(-3, 4, 1)

# FONT SIZES – much larger for readability
tick_fontsize = 26
axis_fontsize = 28
title_fontsize = 28
legend_fontsize = 28
line_width = 2.5
ribbon_alpha = 0.20

# Loop over rows (coasts) and columns (SPEI timescales)
for i, coast in enumerate(coast_order):
    for j, spei in enumerate(spei_order):
        ax = axes[i, j]

        # Filter data for this coast and SPEI timescale
        sub = agg_df[(agg_df["coast_region"] == coast) &
                     (agg_df["SPEI_timescale"] == spei)]

        # Plot each water class
        for wc in water_class_order:
            wc_data = sub[sub["water_class"] == wc]
            if not wc_data.empty:
                # Sort by SPEI_value for smooth lines
                wc_data = wc_data.sort_values("SPEI_value")
                x = wc_data["SPEI_value"].values
                y = wc_data["mean_pct"].values
                sd = wc_data["sd_pct"].values

                # Line – thicker now
                ax.plot(x, y,
                        color=ecosystem_palette[wc],
                        linewidth=line_width,
                        label=wc if (i == 0 and j == 0) else "")

                # Ribbon (mean ± sd)
                ax.fill_between(x,
                                y - sd,
                                y + sd,
                                color=ecosystem_palette[wc],
                                alpha=ribbon_alpha,
                                linewidth=0)

        # Reference lines
        ax.axhline(0, color="black", linestyle="-", linewidth=1.5, alpha=0.7)
        ax.axvline(-1, color="gray", linestyle="--", linewidth=1.5, alpha=0.7)
        ax.axvline(1, color="gray", linestyle="--", linewidth=1.5, alpha=0.7)

        # Axes limits and ticks
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-y_lim, y_lim)
        ax.set_xticks(x_ticks)
        ax.set_yticks(np.arange(-y_lim, y_lim + 0.1, 10))

        # Grid (very light)
        ax.grid(True, linestyle=":", color="lightgray", alpha=0.5)

        # Remove top and right spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Column titles (only top row) – larger font
        if i == 0:
            ax.set_title(spei_labels[j], fontsize=title_fontsize, pad=12)

        # Row labels (only leftmost column) – with extra padding
        if j == 0:
            ax.set_ylabel(coast_labels[i], fontsize=axis_fontsize, labelpad=25, rotation=90)
            # Show y tick labels only in first column
            ax.tick_params(axis="y", labelsize=tick_fontsize)
        else:
            ax.tick_params(axis="y", labelleft=False)

        # Show x tick labels only in bottom row
        if i == 3:
            ax.tick_params(axis="x", labelsize=tick_fontsize)
        else:
            ax.tick_params(axis="x", labelbottom=False)

# -----------------------------------------------------------------------------
# SHARED AXIS LABELS – placed with fig.text, larger font
# -----------------------------------------------------------------------------
fig.text(0.50, 0.07, "SPEI", ha="center", fontsize=axis_fontsize)
fig.text(0.001, 0.5, r"WUE$_{T}$ change from near-normal (%)",
         va="center", rotation=90, fontsize=axis_fontsize)

# -----------------------------------------------------------------------------
# LEGEND – placed well below panels with larger font
# -----------------------------------------------------------------------------
handles = [
    plt.Line2D([0], [0], color=ecosystem_palette[wc], lw=4, label=wc)
    for wc in water_class_order
]
fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.50, 0.01),
    ncol=3,
    fontsize=legend_fontsize,
    frameon=True
)

# -----------------------------------------------------------------------------
# SAVE AND SHOW
# -----------------------------------------------------------------------------
plt.savefig(output_file, dpi=900, bbox_inches="tight")   # <--- HIGHER RESOLUTION
print(f"Figure saved to: {output_file}")

plt.show(block=True)   # Display in Spyder / console

# -----------------------------------------------------------------------------
# CONSOLE SUMMARIES
# -----------------------------------------------------------------------------
# Top 10 combos by mean maximum absolute predicted WUE_T change
top_combos = (
    agg_df.groupby(["coast_region", "water_class", "SPEI_timescale"], observed=True)
    .apply(lambda g: g["mean_pct"].abs().max())
    .reset_index(name="max_abs_mean_pct")
    .sort_values("max_abs_mean_pct", ascending=False)
    .head(10)
)
print("\nTop 10 coast × water_class × SPEI_timescale combinations by mean max absolute change:")
print(top_combos.to_string(index=False))

# -----------------------------------------------------------------------------
# TAKE-HOME MESSAGE (refined)
# -----------------------------------------------------------------------------
print("\nTake-home message (coast-threshold GAM sensitivity):")
print("  - Pacific Coast shows the strongest short-to-intermediate nonlinear WUE_T response, especially around SPEI-3.")
print("  - Alaska shows dry-side WUE_T decreases at longer SPEI timescales, consistent with the long-timescale impairment-style signal.")
print("  - Gulf Coast shows an attenuated WUE_T response across SPEI timescales.")
print("  - Atlantic Coast shows weaker to moderate responses, with stronger responses at some intermediate and longer SPEI timescales.")
print("\nThese curves come from the coast-threshold GAM evaluated under Upland, Freshwater, and Saline water-class settings; they are not independent ecosystem-specific smooths.")
print("Ribbons show variability across growing-season months, not 95% prediction intervals.")

print("\nAll done.")