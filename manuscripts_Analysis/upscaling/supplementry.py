# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 16:30:00 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
Python replicate of Supplementary.R for SPEI-3 coast supplementary figures.

Uses Python/manuscript updated workflow outputs only.

Input:
  M:/Research/WUE_CUE/WUE_manuscript_version6/upscaling/
    EDI_response_v2_upland_monthly_coast_summary.csv

Fallback:
    EDI_logistic_monthly_coast_summary.csv

Outputs:
  Supplementary_table_01_annual_mean_SPEI3_by_coast.csv
  Supplementary_table_02_monthly_cumulative_SPEI3_by_coast.csv
  Supplementary_table_03_annual_cumulative_SPEI3_by_coast.csv
  Supplementary_plot_01_mean_SPEI3_timeseries_by_coast.png
  Supplementary_plot_02_cumulative_SPEI3_by_coast.png
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ------------------------------------------------------------------
# Python/manuscript updated paths
# ------------------------------------------------------------------
output_dir = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")

input_file = output_dir / "EDI_response_v2_upland_monthly_coast_summary.csv"
fallback_file = output_dir / "EDI_logistic_monthly_coast_summary.csv"

if not input_file.exists():
    if not fallback_file.exists():
        raise FileNotFoundError(
            "Could not find monthly coast summary. Expected one of:\n"
            f"{input_file}\n{fallback_file}\n"
            "Run updated 24-EDI_spatial.py first."
        )
    input_file = fallback_file

print(f"Using input file: {input_file}")

# ------------------------------------------------------------------
# Settings matching Supplementary.R
# ------------------------------------------------------------------
coast_levels = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]

coast_labels = {
    "AK Coast": "Alaska Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast",
}

coast_colors = {
    "AK Coast": "#4D4D4D",
    "Pacific Coast": "#C44E52",
    "Gulf Coast": "#8C564B",
    "Atlantic Coast": "#008B8B",
}

# ------------------------------------------------------------------
# Read and validate input
# ------------------------------------------------------------------
monthly = pd.read_csv(input_file)

required_cols = ["year", "month", "coast_region", "mean_SPEI3"]
missing_cols = [c for c in required_cols if c not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

monthly = monthly.copy()
monthly["date"] = pd.to_datetime(
    monthly["year"].astype(int).astype(str)
    + "-"
    + monthly["month"].astype(int).astype(str).str.zfill(2)
    + "-15"
)

monthly = monthly[
    monthly["coast_region"].isin(coast_levels)
    & np.isfinite(monthly["mean_SPEI3"])
].copy()

monthly["coast_region"] = pd.Categorical(
    monthly["coast_region"],
    categories=coast_levels,
    ordered=True
)

monthly = monthly.sort_values(["coast_region", "date"])

print("Months present:", sorted(monthly["month"].unique()))
print("Rows:", len(monthly))

# ------------------------------------------------------------------
# Annual mean SPEI-3 table
# ------------------------------------------------------------------
annual = (
    monthly
    .groupby(["coast_region", "year"], observed=True, as_index=False)
    .agg(annual_mean_SPEI3=("mean_SPEI3", "mean"))
)

# ------------------------------------------------------------------
# Monthly cumulative table
# ------------------------------------------------------------------
cumulative = monthly.copy()
cumulative["cumulative_SPEI3"] = (
    cumulative
    .groupby("coast_region", observed=True)["mean_SPEI3"]
    .cumsum()
)

# ------------------------------------------------------------------
# Annual cumulative table
# Same logic as Supplementary.R:
# cumulative sum of annual_mean_SPEI3
# ------------------------------------------------------------------
annual_cumulative = annual.copy()
annual_cumulative = annual_cumulative.sort_values(["coast_region", "year"])
annual_cumulative["cumulative_annual_SPEI3"] = (
    annual_cumulative
    .groupby("coast_region", observed=True)["annual_mean_SPEI3"]
    .cumsum()
)

# ------------------------------------------------------------------
# Write tables
# ------------------------------------------------------------------
annual.to_csv(
    output_dir / "Supplementary_table_01_annual_mean_SPEI3_by_coast.csv",
    index=False
)

cumulative.to_csv(
    output_dir / "Supplementary_table_02_monthly_cumulative_SPEI3_by_coast.csv",
    index=False
)

annual_cumulative.to_csv(
    output_dir / "Supplementary_table_03_annual_cumulative_SPEI3_by_coast.csv",
    index=False
)

# ------------------------------------------------------------------
# Helper for stacked coast panels
# ------------------------------------------------------------------
def format_time_axis(ax):
    ax.xaxis.set_major_locator(mdates.YearLocator(base=5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(axis="y", color="0.90", linewidth=0.7)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ------------------------------------------------------------------
# Plot 01: Monthly mean SPEI-3 by coast
# ------------------------------------------------------------------
fig, axes = plt.subplots(
    nrows=4,
    ncols=1,
    figsize=(7.2, 8.2),
    sharex=True
)

date_min = monthly["date"].min()
date_max = monthly["date"].max()

for ax, coast in zip(axes, coast_levels):
    sub = monthly[monthly["coast_region"] == coast]
    ann = annual[annual["coast_region"] == coast].copy()
    ann["date"] = pd.to_datetime(ann["year"].astype(int).astype(str) + "-07-01")

    ax.axhspan(-1, 1, color="gray", alpha=0.18)
    ax.axhline(0, color="0.35", linewidth=0.45)
    ax.axhline(-1, color="0.50", linewidth=0.40, linestyle="--")
    ax.axhline(1, color="0.50", linewidth=0.40, linestyle="--")

    ax.plot(
        sub["date"],
        sub["mean_SPEI3"],
        color=coast_colors[coast],
        linewidth=0.55,
        alpha=0.85
    )

    ax.plot(
        ann["date"],
        ann["annual_mean_SPEI3"],
        color=coast_colors[coast],
        linewidth=1.25
    )

    ax.set_title(coast_labels[coast], fontweight="bold", fontsize=12)
    ax.set_ylabel("Mean SPEI-3")
    ax.set_xlim(date_min, date_max)
    format_time_axis(ax)

axes[-1].set_xlabel("")
fig.suptitle("Monthly Mean SPEI-3 by Coastal Region", fontweight="bold", fontsize=14, y=0.995)
fig.text(
    0.5,
    0.955,
    "Thin lines show monthly raster means; thick lines show annual means. Gray band marks near-normal SPEI-3 (-1 to 1).",
    ha="center",
    fontsize=9
)

fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(
    output_dir / "Supplementary_plot_01_mean_SPEI3_timeseries_by_coast.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)
plt.close(fig)

# ------------------------------------------------------------------
# Plot 02: Cumulative mean SPEI-3 by coast
# ------------------------------------------------------------------
fig, axes = plt.subplots(
    nrows=4,
    ncols=1,
    figsize=(7.2, 8.2),
    sharex=True
)

for ax, coast in zip(axes, coast_levels):
    sub = cumulative[cumulative["coast_region"] == coast]
    ann = annual_cumulative[annual_cumulative["coast_region"] == coast].copy()
    ann["date"] = pd.to_datetime(ann["year"].astype(int).astype(str) + "-07-01")

    ax.axhline(0, color="0.35", linewidth=0.55)

    ax.plot(
        sub["date"],
        sub["cumulative_SPEI3"],
        color=coast_colors[coast],
        linewidth=0.70,
        alpha=0.88
    )

    ax.plot(
        ann["date"],
        ann["cumulative_annual_SPEI3"],
        color=coast_colors[coast],
        linewidth=1.45
    )

    ax.set_title(coast_labels[coast], fontweight="bold", fontsize=12)
    ax.set_ylabel("Cumulative mean SPEI-3")
    ax.set_xlim(date_min, date_max)
    format_time_axis(ax)

axes[-1].set_xlabel("")
fig.suptitle("Cumulative Mean SPEI-3 by Coastal Region", fontweight="bold", fontsize=14, y=0.995)
fig.text(
    0.5,
    0.955,
    "Thin lines show cumulative monthly raster means; thick lines show cumulative annual means.",
    ha="center",
    fontsize=9
)

fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(
    output_dir / "Supplementary_plot_02_cumulative_SPEI3_by_coast.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)
plt.close(fig)

print("=" * 80)
print("Supplementary outputs written to:")
print(output_dir)
print("=" * 80)
print(output_dir / "Supplementary_plot_01_mean_SPEI3_timeseries_by_coast.png")
print(output_dir / "Supplementary_plot_02_cumulative_SPEI3_by_coast.png")
print(output_dir / "Supplementary_table_01_annual_mean_SPEI3_by_coast.csv")
print(output_dir / "Supplementary_table_02_monthly_cumulative_SPEI3_by_coast.csv")
print(output_dir / "Supplementary_table_03_annual_cumulative_SPEI3_by_coast.csv")