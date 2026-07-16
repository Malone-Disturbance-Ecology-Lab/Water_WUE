# -*- coding: utf-8 -*-
"""
Q3_combined_figure_2x4_GridSpec.py

Single figure with:
- Top-left: Upland (4×4)
- Top-right: Freshwater (4×4)
- Bottom (centered): Saline (4×4)

All blocks have the same width and height.
GridSpec layout: Upland = row0, cols0-2; Freshwater = row0, cols2-4; Saline = row1, cols1-3
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
input_csv = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs\sensitivity_conditions_results\prediction_curves_all_reference_conditions.csv"
output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_figures\talib_manuscript"
os.makedirs(output_dir, exist_ok=True)

coast_palette = {
    "Atlantic Coast": "#2E8B57",
    "Pacific Coast":  "#DC143C",
    "Gulf Coast":     "#00CED1",
    "AK Coast":       "#8B4513"
}
coast_order = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
coast_labels_plot = ["Atlantic", "Pacific", "Gulf", "Alaska"]

ecosystem_order = ["Upland", "Freshwater", "Saline"]

spei_order = ["SPEI_3", "SPEI_24", "SPEI_36", "SPEI_48"]
spei_labels = ["SPEI-3", "SPEI-24", "SPEI-36", "SPEI-48"]

# -----------------------------------------------------------------------------
# LOAD AND AGGREGATE DATA
# -----------------------------------------------------------------------------
df = pd.read_csv(input_csv)
df["month_f"] = df["month_f"].astype(int)
df = df[df["SPEI_timescale"].isin(spei_order)]

agg_df = df.groupby(
    ["coast_region", "water_class", "SPEI_timescale", "SPEI_value"],
    observed=True
).agg(
    mean_pct=("predicted_pct_change", "mean"),
    lower_pct=("predicted_lower_pct", "mean"),
    upper_pct=("predicted_upper_pct", "mean"),
    n_months=("month_f", "nunique")
).reset_index()

# -----------------------------------------------------------------------------
# Y-LIMITS PER (ECOSYSTEM, COAST)
# -----------------------------------------------------------------------------
def nice_y_limit(x, base=10, min_lim=20, pad=1.05):
    if not np.isfinite(x) or x <= 0:
        return min_lim
    return max(min_lim, int(np.ceil((x * pad) / base) * base))

y_lims = {}
for eco in ecosystem_order:
    eco_data = agg_df[agg_df["water_class"] == eco]
    for coast in coast_order:
        coast_data = eco_data[eco_data["coast_region"] == coast]
        if coast_data.empty:
            y_lims[(eco, coast)] = 20
            continue
        max_abs = np.nanmax(np.abs(coast_data[["mean_pct", "lower_pct", "upper_pct"]].to_numpy()))
        y_lims[(eco, coast)] = nice_y_limit(max_abs, base=10, min_lim=20)

# -----------------------------------------------------------------------------
# FUNCTION TO CREATE ONE 4×4 ECOSYSTEM BLOCK
# -----------------------------------------------------------------------------
def fill_ecosystem_block(fig, outer_spec, ecosystem, panel_title):
    """
    Create one 4×4 ecosystem block inside a fixed GridSpec slot.
    """
    inner = outer_spec.subgridspec(
        4, 4,
        wspace=0.12,
        hspace=0.16
    )

    axes = np.empty((4, 4), dtype=object)

    # FONT SIZES
    spei_title_fontsize = 26
    tick_fontsize = 33
    row_label_fontsize = 30
    section_fontsize = 30
    line_width = 4.0
    ribbon_alpha = 0.10

    # Section title
    bbox = outer_spec.get_position(fig)
    fig.text(
        bbox.x0,
        bbox.y1 + 0.018,
        panel_title,
        fontsize=section_fontsize,
        weight="bold",
        ha="left",
        va="bottom"
    )

    for coast_idx, coast in enumerate(coast_order):
        for col, spei in enumerate(spei_order):
            ax = fig.add_subplot(inner[coast_idx, col])
            axes[coast_idx, col] = ax

            sub = agg_df[
                (agg_df["coast_region"] == coast) &
                (agg_df["water_class"] == ecosystem) &
                (agg_df["SPEI_timescale"] == spei)
            ]

            if not sub.empty:
                sub = sub.sort_values("SPEI_value")
                x = sub["SPEI_value"].values
                y = sub["mean_pct"].values
                lower = sub["lower_pct"].values
                upper = sub["upper_pct"].values
                color = coast_palette[coast]

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
            ax.axhline(0, color="black", linestyle="-", linewidth=1.4, alpha=0.65)
            ax.axvline(-1, color="gray", linestyle="--", linewidth=1.4, alpha=0.65)
            ax.axvline(1, color="gray", linestyle="--", linewidth=1.4, alpha=0.65)

            # Row-specific y-limit
            ylim = y_lims[(ecosystem, coast)]
            ax.set_ylim(-ylim, ylim)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

            ax.set_xlim(-3.0, 3.0)
            ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])

            ax.grid(True, linestyle=":", color="lightgray", alpha=0.4)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # SPEI titles only on top row
            if coast_idx == 0:
                ax.set_title(spei_labels[col], fontsize=spei_title_fontsize, pad=8)

            # Coast labels only on leftmost column
            if col == 0:
                ax.set_ylabel(
                    coast_labels_plot[coast_idx],
                    fontsize=row_label_fontsize,
                    labelpad=6,          # reduced from 10 to reduce overlap
                    rotation=90
                )
                ax.tick_params(axis="y", labelsize=tick_fontsize, pad=8)
            else:
                ax.tick_params(axis="y", labelleft=False)

            # X tick labels only on bottom row of each block
            if coast_idx == 3:
                ax.tick_params(axis="x", labelsize=tick_fontsize, pad=8)
            else:
                ax.tick_params(axis="x", labelbottom=False)

    return axes

# -----------------------------------------------------------------------------
# CREATE MAIN FIGURE – INCREASED LEFT MARGIN FOR Y-LABEL
# -----------------------------------------------------------------------------
fig = plt.figure(figsize=(36, 28))

outer = fig.add_gridspec(
    2, 4,
    left=0.1,             # increased from 0.05 to give space for y-label
    right=0.98,
    top=0.96,
    bottom=0.07,
    wspace=0.35,
    hspace=0.40
)

upland_axes = fill_ecosystem_block(
    fig,
    outer[0, 0:2],
    "Upland",
    "a) Upland"
)

freshwater_axes = fill_ecosystem_block(
    fig,
    outer[0, 2:4],
    "Freshwater",
    "b) Freshwater"
)

saline_axes = fill_ecosystem_block(
    fig,
    outer[1, 1:3],
    "Saline",
    "c) Saline"
)

# -----------------------------------------------------------------------------
# GLOBAL AXIS LABELS – repositioned to avoid coast labels
# -----------------------------------------------------------------------------
fig.text(
    0.50, 0.035,
    "SPEI",
    ha="center",
    fontsize=32
)

fig.text(
    0.01, 0.73,                # moved slightly right, but still outside axes
    r"WUE$_{T}$ change from near-normal (%)",
    va="center",
    rotation=90,
    fontsize=32
)

# -----------------------------------------------------------------------------
# SAVE WITH DPI = 600
# -----------------------------------------------------------------------------
output_file = os.path.join(output_dir, "Q3_combined_2row_layout_equal_blocks.png")
fig.savefig(output_file, dpi=600, bbox_inches="tight")
print(f"High-res combined figure saved to: {output_file}")

plt.show()
print("\nDone.")