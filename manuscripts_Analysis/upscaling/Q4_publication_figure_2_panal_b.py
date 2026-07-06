# -*- coding: utf-8 -*-
"""
Panel B – Annual impacted-area time series (final)
- Four stacked subplots (one per coast)
- Shows dry area, negative events, positive events, any event
- Saves as 'EDI_panel_B_timeseries.png' in talib_publication
- Displays in Spyder
"""

import csv
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
import numpy as np

# ---------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------
SPATIAL_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
OUT_DIR = SPATIAL_DIR / "talib_publication"
OUT_DIR.mkdir(parents=True, exist_ok=True)

annual_path = SPATIAL_DIR / "EDI_response_communication_annual_impacted_area_by_coast_2000_2025.csv"
OUT_FIG = OUT_DIR / "EDI_panel_B_timeseries.png"

# ---------------------------------------------------------------------
# COAST REGIONS
# ---------------------------------------------------------------------
REGIONS = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]
REGION_DISPLAY = {
    "AK Coast": "Alaska Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast",
}
def display_region(region):
    return REGION_DISPLAY.get(region, region)

# Colours for time-series elements
DRY_AREA_COLOR = "#C9B037"
NEGATIVE_EVENT_COLOR = "#2C7FB8"
POSITIVE_EVENT_COLOR = "#D95F0E"
ANY_EVENT_COLOR = "#2B2B2B"

# ---------------------------------------------------------------------
# READ CSV
# ---------------------------------------------------------------------
def read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

annual_rows = read_csv_rows(annual_path)

# ---------------------------------------------------------------------
# HELPER FOR Y‑AXIS LIMIT
# ---------------------------------------------------------------------
def percent_axis_limit(*series, minimum=8, maximum=100):
    vals = []
    for seq in series:
        vals.extend([v for v in seq if np.isfinite(v)])
    if not vals:
        return minimum
    return min(maximum, max(minimum, np.nanmax(vals) * 1.18))

# ---------------------------------------------------------------------
# CREATE PANEL B
# ---------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 5.0), facecolor="white")
gs = gridspec.GridSpec(4, 1, hspace=0.34)
axes = [fig.add_subplot(gs[i, 0]) for i in range(4)]

# Font sizes
FONT_AXIS_LABEL = 12
FONT_TICK = 10
FONT_REGION = 10
FONT_LEGEND = 10
FONT_PANEL = 12

for ax, region in zip(axes, REGIONS):
    sub = [r for r in annual_rows if r["coast_region"] == region]
    years = [int(r["year"]) for r in sub]
    neg_area = [float(r["annual_mean_negative_impacted_area_pct"]) for r in sub]
    pos_area = [float(r["annual_mean_positive_impacted_area_pct"]) for r in sub]
    any_area = [float(r["annual_mean_any_impacted_area_pct"]) for r in sub]
    dry_area = [float(r["annual_mean_dry_area_pct"]) for r in sub]

    ax.fill_between(years, 0, dry_area, color=DRY_AREA_COLOR, alpha=0.32, label="Dry area")
    ax.plot(years, neg_area, color=NEGATIVE_EVENT_COLOR, marker="o", markersize=2.2, linewidth=1.2,
            label="Negative event area")
    ax.plot(years, pos_area, color=POSITIVE_EVENT_COLOR, marker="o", markersize=2.2, linewidth=1.2,
            label="Positive event area")
    ax.plot(years, any_area, color=ANY_EVENT_COLOR, linestyle="--", linewidth=0.9,
            label="Any signed event area")

    ylim = percent_axis_limit(neg_area, pos_area, any_area, dry_area)
    ax.set_ylim(0, ylim)
    ax.set_xlim(1999.5, 2025.5)

    # Y‑axis: numbers only (no % sign), at least 4 ticks
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    ax.set_ylabel("Area (%)", fontsize=FONT_AXIS_LABEL)
    ax.tick_params(labelsize=FONT_TICK)

    # Region label – moved higher, smaller font
    ax.text(
        0.01, 0.94, display_region(region),
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=FONT_REGION, fontweight="bold",
    )

    # Aesthetics
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="0.90")

    if ax is not axes[-1]:
        ax.tick_params(labelbottom=False)

axes[-1].set_xlabel("Year", fontsize=FONT_AXIS_LABEL)

# Shared legend (above the panels)
handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(
    handles, labels,
    loc="upper center",
    bbox_to_anchor=(0.60, 1.42),
    ncol=4,
    fontsize=FONT_LEGEND,
    frameon=False,
    columnspacing=0.8,
    handlelength=1.8,
)

# Panel label "b)" – outside the plot area
axes[0].text(-0.08, 1.22, "b)", transform=axes[0].transAxes,
             fontsize=FONT_PANEL, fontweight="bold")

# No suptitle

# ---------------------------------------------------------------------
# SAVE AND DISPLAY
# ---------------------------------------------------------------------
fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Panel B saved to: {OUT_FIG}")

plt.show(block=True)
print("Panel B displayed.")