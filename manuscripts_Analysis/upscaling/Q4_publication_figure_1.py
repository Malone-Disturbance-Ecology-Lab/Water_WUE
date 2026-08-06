# -*- coding: utf-8 -*-
"""
Recreate Supplementary_plot_02 with:
- New coast colours
- Shared y‑axis label (bold, closer to figure)
- No main title / footnote
- No panel labels
- At least 5 y‑ticks per panel
- Thicker, longer x‑axis tick marks
- Bold tick labels and axis labels
- Display in Spyder
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
original_output = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
new_output = original_output / "talib_publication"
new_output.mkdir(parents=True, exist_ok=True)

monthly_cum_file = original_output / "Supplementary_table_02_monthly_cumulative_SPEI3_by_coast.csv"
annual_cum_file  = original_output / "Supplementary_table_03_annual_cumulative_SPEI3_by_coast.csv"

# ------------------------------------------------------------------
# Settings – coast order, labels, and colours
# ------------------------------------------------------------------
# ORDER: Alaska, Pacific, Gulf, Atlantic
coast_levels = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]

coast_labels = {
    "AK Coast": "Alaska Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast",
}

# Colours as per your latest request
coast_colors = {
    'Pacific Coast': '#0072B2',   # blue
    'Gulf Coast':   '#B3B300',    # olive/khaki
    'Atlantic Coast': '#CC79A7',  # pink
    'AK Coast':     '#009E73'     # teal
}

# ------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------
monthly_cum = pd.read_csv(monthly_cum_file)
annual_cum  = pd.read_csv(annual_cum_file)

monthly_cum["date"] = pd.to_datetime(
    monthly_cum["year"].astype(str) + "-" +
    monthly_cum["month"].astype(str).str.zfill(2) + "-15"
)
annual_cum["date"] = pd.to_datetime(
    annual_cum["year"].astype(str) + "-07-01"
)

monthly_cum = monthly_cum[monthly_cum["coast_region"].isin(coast_levels)].copy()
annual_cum  = annual_cum[annual_cum["coast_region"].isin(coast_levels)].copy()

monthly_cum["coast_region"] = pd.Categorical(
    monthly_cum["coast_region"], categories=coast_levels, ordered=True
)
annual_cum["coast_region"] = pd.Categorical(
    annual_cum["coast_region"], categories=coast_levels, ordered=True
)

monthly_cum = monthly_cum.sort_values(["coast_region", "date"])
annual_cum  = annual_cum.sort_values(["coast_region", "year"])

date_min = monthly_cum["date"].min()
date_max = monthly_cum["date"].max()

# ------------------------------------------------------------------
# Helper: format axis with ≥5 y‑ticks, clear x‑tick marks, bold ticks
# ------------------------------------------------------------------
def format_axis(ax):
    # X‑axis: 5‑year spacing, longer/thicker ticks
    ax.xaxis.set_major_locator(mdates.YearLocator(base=5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.tick_params(axis='x', length=8, width=2)
    ax.tick_params(axis='y', length=6, width=1.5)
    # Enforce at least 5 y‑ticks (use nbins=6)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune=None))
    ax.grid(axis="y", color="0.90", linewidth=0.7)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # Make tick labels bold and larger
    for label in ax.get_xticklabels():
        label.set_fontsize(14)
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontsize(14)
        label.set_fontweight('bold')

# ------------------------------------------------------------------
# Build figure
# ------------------------------------------------------------------
fig, axes = plt.subplots(
    nrows=4, ncols=1, figsize=(7.2, 8.2), sharex=True
)

for ax, coast in zip(axes, coast_levels):
    sub = monthly_cum[monthly_cum["coast_region"] == coast]
    ann = annual_cum[annual_cum["coast_region"] == coast]

    # Zero line
    ax.axhline(0, color="0.35", linewidth=0.55)

    # Monthly cumulative (thicker)
    ax.plot(
        sub["date"],
        sub["cumulative_SPEI3"],
        color=coast_colors[coast],
        linewidth=1.0,
        alpha=0.88
    )

    # Annual cumulative (thicker)
    ax.plot(
        ann["date"],
        ann["cumulative_annual_SPEI3"],
        color=coast_colors[coast],
        linewidth=2.0,
    )

    ax.set_title(coast_labels[coast], fontweight="bold", fontsize=14)
    ax.set_ylabel("")                         # remove individual y-label
    ax.set_xlim(date_min, date_max)
    format_axis(ax)

# X‑axis label only for bottom panel – bold and large
axes[-1].set_xlabel("Year", fontsize=14, fontweight='bold')

# ------------------------------------------------------------------
# Shared y‑axis label – bold, closer to figure
# ------------------------------------------------------------------
fig.text(
    0.08, 0.5,
    "Cumulative mean SPEI-3",
    va="center",
    ha="center",
    rotation="vertical",
    fontsize=14,
    fontweight='bold'
)

# Adjust margins
plt.subplots_adjust(
    left=0.18,
    right=0.98,
    top=0.94,
    bottom=0.08,
    hspace=0.45
)

# ------------------------------------------------------------------
# Save and show
# ------------------------------------------------------------------
out_file = new_output / "Supplementary_plot_02_cumulative_SPEI3_by_coast.png"
fig.savefig(out_file, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Figure saved to: {out_file}")

# Display in Spyder
plt.show(block=True)
print("Figure displayed. If you do not see it, check Spyder Graphics backend settings.")