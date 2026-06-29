import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import gridspec
from matplotlib.ticker import PercentFormatter
import netCDF4 as nc
import numpy as np

# ---------------------------------------------------------------------
# PATHS — Python/manuscript-only version
# ---------------------------------------------------------------------
SPATIAL_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
ECO_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q4\Q4_Ecological_Impacts_outputs")

OUT = SPATIAL_DIR / "EDI_response_communication_plot_03_composite_maps_timeseries_Q3bars.png"
OUT_MAIN = SPATIAL_DIR / "EDI_response_communication_plot_04_main_maps_timeseries.png"
OUT_SUPP = SPATIAL_DIR / "EDI_response_communication_plot_05_supp_Q3bars.png"

REGIONS = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]
REGION_DISPLAY = {
    "AK Coast": "Alaska Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast",
}
REGION_LABEL_POS = {
    "AK Coast": (-171.0, 73.0, "Alaska"),
    "Pacific Coast": (-132.0, 50.0, "Pacific"),
    "Gulf Coast": (-101.0, 22.0, "Gulf"),
    "Atlantic Coast": (-82.0, 50.0, "Atlantic"),
}
REGION_COLORS = {
    "AK Coast": "#4D4D4D",
    "Pacific Coast": "#C44E52",
    "Gulf Coast": "#8C564B",
    "Atlantic Coast": "#008B8B",
}
DRY_AREA_COLOR = "#C9B037"
NEGATIVE_EVENT_COLOR = "#2C7FB8"
POSITIVE_EVENT_COLOR = "#D95F0E"
ANY_EVENT_COLOR = "#2B2B2B"
CLASS_COLORS = {
    "No meaningful change (<5%)": "#EFEFEF",
    "Watch (5-10%)": "#BFD3E6",
    "Stress (10-20%)": "#F4A582",
    "Impact (>=20%)": "#B2182B",
}
PI_COLORS = {
    "PI supports decrease": "#2166AC",
    "PI overlaps no-change": "#D9D9D9",
    "PI supports increase": "#B2182B",
}


def read_raster(path):
    ds = nc.Dataset(path)
    lat = np.array(ds.variables["lat"][:])
    lon = np.array(ds.variables["lon"][:])
    var_name = [v for v in ds.variables if v not in ("lat", "lon")][0]
    arr = np.array(ds.variables[var_name][:], dtype=float)
    fill = getattr(ds.variables[var_name], "_FillValue", None)
    ds.close()
    if fill is not None:
        arr[arr == fill] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return lon, lat, arr


def read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def proportion_table(rows, group_field, value_field, value_levels):
    counts = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)
    for row in rows:
        group = row[group_field]
        value = row[value_field]
        counts[group][value] += 1
        totals[group] += 1
    props = {
        group: {level: (counts[group][level] / totals[group] * 100 if totals[group] else 0)
                for level in value_levels}
        for group in REGIONS
    }
    return props


def style_map_axis(ax, label_fontsize=11):
    ax.set_xlim(-180, -60)
    ax.set_ylim(17, 77)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_facecolor("white")
    ax.grid(False)
    for region, (x, y, label) in REGION_LABEL_POS.items():
        ax.text(
            x, y, label,
            ha="center", va="center",
            fontsize=label_fontsize, fontweight="bold", color="black",
        )


def display_region(region):
    return REGION_DISPLAY.get(region, region)


def percent_axis_limit(*series, minimum=8, maximum=100):
    vals = []
    for seq in series:
        vals.extend([v for v in seq if np.isfinite(v)])
    if not vals:
        return minimum
    return min(maximum, max(minimum, np.nanmax(vals) * 1.18))


def main():
    neg_path = SPATIAL_DIR / "EDI_response_v2_freq_strong_decrease_upland.nc"
    pos_path = SPATIAL_DIR / "EDI_response_v2_freq_strong_increase_upland.nc"
    annual_path = SPATIAL_DIR / "EDI_response_communication_annual_impacted_area_by_coast_2000_2025.csv"
    monthly_path = ECO_DIR / "EDI_Q3_aligned_monthly_scores.csv"

    lon, lat, neg = read_raster(neg_path)
    _, _, pos = read_raster(pos_path)
    annual_rows = read_csv_rows(annual_path)
    monthly_rows = read_csv_rows(monthly_path)

    class_levels = list(CLASS_COLORS)
    pi_levels = list(PI_COLORS)
    class_props = proportion_table(monthly_rows, "coast_region", "EDI_Q3_class", class_levels)
    pi_props = proportion_table(monthly_rows, "coast_region", "prediction_interval_support", pi_levels)

    neg_cmap = mcolors.LinearSegmentedColormap.from_list(
        "negative_events", ["#FFFFFF", "#D7ECF7", "#73B3D8", "#2879B9", "#08306B"]
    )
    pos_cmap = mcolors.LinearSegmentedColormap.from_list(
        "positive_events", ["#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33", "#7F0000"]
    )
    neg_cmap.set_bad((1, 1, 1, 0))
    pos_cmap.set_bad((1, 1, 1, 0))

    fig = plt.figure(figsize=(13, 18), facecolor="white")
    gs = gridspec.GridSpec(
        4, 2,
        height_ratios=[1.35, 2.6, 1.05, 1.05],
        hspace=0.50,
        wspace=0.25,
    )

    # A. Event-frequency maps.
    ax_neg = fig.add_subplot(gs[0, 0])
    ax_pos = fig.add_subplot(gs[0, 1], sharex=ax_neg, sharey=ax_neg)
    neg_vmin, neg_vmax = np.nanpercentile(neg, [2, 98])
    pos_vmin, pos_vmax = np.nanpercentile(pos, [2, 98])
    im_neg = ax_neg.pcolormesh(lon, lat, neg, cmap=neg_cmap, vmin=neg_vmin, vmax=neg_vmax, shading="auto")
    im_pos = ax_pos.pcolormesh(lon, lat, pos, cmap=pos_cmap, vmin=pos_vmin, vmax=pos_vmax, shading="auto")
    style_map_axis(ax_neg)
    style_map_axis(ax_pos)
    ax_pos.set_ylabel("")
    ax_neg.set_title("Negative events: WUE$_T$ <= -5%", fontweight="bold")
    ax_pos.set_title("Positive events: WUE$_T$ >= +5%", fontweight="bold")
    cb1 = fig.colorbar(im_neg, ax=ax_neg, fraction=0.046, pad=0.02)
    cb2 = fig.colorbar(im_pos, ax=ax_pos, fraction=0.046, pad=0.02)
    cb1.set_label("% months, 2000-2025")
    cb2.set_label("% months, 2000-2025")
    ax_neg.text(-0.10, 1.10, "A", transform=ax_neg.transAxes, fontsize=18, fontweight="bold")

    # B. Annual impacted area.
    ts_gs = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[1, :], hspace=0.38)
    ts_axes = [fig.add_subplot(ts_gs[i, 0]) for i in range(4)]
    for ax, region in zip(ts_axes, REGIONS):
        sub = [r for r in annual_rows if r["coast_region"] == region]
        years = [int(r["year"]) for r in sub]
        neg_area = [float(r["annual_mean_negative_impacted_area_pct"]) for r in sub]
        pos_area = [float(r["annual_mean_positive_impacted_area_pct"]) for r in sub]
        any_area = [float(r["annual_mean_any_impacted_area_pct"]) for r in sub]
        dry_area = [float(r["annual_mean_dry_area_pct"]) for r in sub]
        ax.fill_between(years, 0, dry_area, color=DRY_AREA_COLOR, alpha=0.32, label="Dry area")
        ax.plot(years, neg_area, color=NEGATIVE_EVENT_COLOR, marker="o", markersize=3.4, linewidth=1.8,
                label="Negative event area")
        ax.plot(years, pos_area, color=POSITIVE_EVENT_COLOR, marker="o", markersize=3.4, linewidth=1.8,
                label="Positive event area")
        ax.plot(years, any_area, color=ANY_EVENT_COLOR, linestyle="--", linewidth=1.3,
                label="Any signed event area")
        ax.set_ylim(0, percent_axis_limit(neg_area, pos_area, any_area, dry_area))
        ax.set_xlim(1999.5, 2025.5)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        ax.set_ylabel("Area (%)")
        ax.text(
            0.01, 0.88, display_region(region),
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=11, fontweight="bold",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="0.90")
        if ax is not ts_axes[-1]:
            ax.tick_params(labelbottom=False)
    ts_axes[-1].set_xlabel("Year")
    handles, labels = ts_axes[0].get_legend_handles_labels()
    ts_axes[0].legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.58, 1.38),
        ncol=4,
        fontsize=8.5,
        frameon=False,
    )
    ts_axes[0].text(-0.08, 1.22, "B", transform=ts_axes[0].transAxes, fontsize=18, fontweight="bold")

    # C. Monthly Q3 EDI class proportions.
    ax_c = fig.add_subplot(gs[2:, 0])
    y = np.arange(len(REGIONS))
    left = np.zeros(len(REGIONS))
    for level in class_levels:
        vals = np.array([class_props[region][level] for region in REGIONS])
        ax_c.barh(y, vals, left=left, color=CLASS_COLORS[level], edgecolor="white", linewidth=0.5, label=level)
        left += vals
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([display_region(region) for region in REGIONS])
    ax_c.invert_yaxis()
    ax_c.set_xlim(0, 100)
    ax_c.xaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax_c.set_xlabel("% of upland site-months")
    ax_c.set_title("Monthly SPEI-3 Upland EDI Classes", fontweight="bold")
    ax_c.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=1, fontsize=8, frameon=False)
    ax_c.grid(axis="x", color="0.90")
    ax_c.set_axisbelow(True)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.text(-0.14, 1.06, "C", transform=ax_c.transAxes, fontsize=18, fontweight="bold")

    # D. Prediction-interval support.
    ax_d = fig.add_subplot(gs[2:, 1])
    left = np.zeros(len(REGIONS))
    for level in pi_levels:
        vals = np.array([pi_props[region][level] for region in REGIONS])
        ax_d.barh(y, vals, left=left, color=PI_COLORS[level], edgecolor="white", linewidth=0.5, label=level)
        left += vals
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([display_region(region) for region in REGIONS])
    ax_d.invert_yaxis()
    ax_d.set_xlim(0, 100)
    ax_d.xaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax_d.set_xlabel("% of upland site-months")
    ax_d.set_title("Prediction-Interval Support", fontweight="bold")
    ax_d.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=1, fontsize=8, frameon=False)
    ax_d.grid(axis="x", color="0.90")
    ax_d.set_axisbelow(True)
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)
    ax_d.text(-0.14, 1.06, "D", transform=ax_d.transAxes, fontsize=18, fontweight="bold")

    fig.suptitle(
        "Spatial and Site-Based Upland WUE$_T$ Drought Response",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Main-text version: maps + annual impacted-area time series.
    fig_main = plt.figure(figsize=(7.2, 7.5), facecolor="white")
    main_gs = gridspec.GridSpec(
        2, 2,
        height_ratios=[1.0, 1.55],
        hspace=0.42,
        wspace=0.38,
    )
    ax_main_neg = fig_main.add_subplot(main_gs[0, 0])
    ax_main_pos = fig_main.add_subplot(main_gs[0, 1], sharex=ax_main_neg, sharey=ax_main_neg)
    im_main_neg = ax_main_neg.pcolormesh(lon, lat, neg, cmap=neg_cmap, vmin=neg_vmin, vmax=neg_vmax, shading="auto")
    im_main_pos = ax_main_pos.pcolormesh(lon, lat, pos, cmap=pos_cmap, vmin=pos_vmin, vmax=pos_vmax, shading="auto")
    for ax in (ax_main_neg, ax_main_pos):
        style_map_axis(ax, label_fontsize=7)
        ax.tick_params(labelsize=7)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
    ax_main_pos.set_ylabel("")
    ax_main_neg.set_title("Negative events\nWUE$_T$ <= -5%", fontweight="bold", fontsize=9)
    ax_main_pos.set_title("Positive events\nWUE$_T$ >= +5%", fontweight="bold", fontsize=9)
    cb_main_neg = fig_main.colorbar(im_main_neg, ax=ax_main_neg, fraction=0.045, pad=0.02)
    cb_main_pos = fig_main.colorbar(im_main_pos, ax=ax_main_pos, fraction=0.045, pad=0.02)
    cb_main_neg.set_label("% months", fontsize=7, labelpad=2)
    cb_main_pos.set_label("% months", fontsize=7, labelpad=2)
    cb_main_neg.ax.tick_params(labelsize=7)
    cb_main_pos.ax.tick_params(labelsize=7)
    ax_main_neg.text(-0.16, 1.12, "A", transform=ax_main_neg.transAxes, fontsize=13, fontweight="bold")

    main_ts_gs = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=main_gs[1, :], hspace=0.34)
    main_ts_axes = [fig_main.add_subplot(main_ts_gs[i, 0]) for i in range(4)]
    for ax, region in zip(main_ts_axes, REGIONS):
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
        ax.set_ylim(0, percent_axis_limit(neg_area, pos_area, any_area, dry_area))
        ax.set_xlim(1999.5, 2025.5)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        ax.set_ylabel("Area (%)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.text(
            0.01, 0.88, display_region(region),
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=8.5, fontweight="bold",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="0.90")
        if ax is not main_ts_axes[-1]:
            ax.tick_params(labelbottom=False)
    main_ts_axes[-1].set_xlabel("Year", fontsize=8)
    handles, labels = main_ts_axes[0].get_legend_handles_labels()
    main_ts_axes[0].legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.60, 1.42),
        ncol=4,
        fontsize=6.5,
        frameon=False,
        columnspacing=0.8,
        handlelength=1.8,
    )
    main_ts_axes[0].text(-0.08, 1.22, "B", transform=main_ts_axes[0].transAxes, fontsize=13, fontweight="bold")
    fig_main.suptitle(
        "Upland WUE$_T$ Event Frequency and Impacted Area",
        fontsize=11.5,
        fontweight="bold",
        y=0.995,
    )
    fig_main.savefig(OUT_MAIN, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_main)

    # Supplemental version: Q3 class proportions + uncertainty support.
    fig_supp, (ax_supp_c, ax_supp_d) = plt.subplots(
        1, 2, figsize=(7.2, 5.8), facecolor="white", sharey=False
    )
    y = np.arange(len(REGIONS))
    left = np.zeros(len(REGIONS))
    for level in class_levels:
        vals = np.array([class_props[region][level] for region in REGIONS])
        ax_supp_c.barh(y, vals, left=left, color=CLASS_COLORS[level], edgecolor="white", linewidth=0.5, label=level)
        left += vals
    ax_supp_c.set_yticks(y)
    ax_supp_c.set_yticklabels([display_region(region) for region in REGIONS], fontsize=8)
    ax_supp_c.invert_yaxis()
    ax_supp_c.set_xlim(0, 100)
    ax_supp_c.xaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax_supp_c.set_xlabel("% of upland site-months", fontsize=8)
    ax_supp_c.set_title("Monthly SPEI-3 Upland EDI Classes", fontweight="bold", fontsize=9)
    ax_supp_c.tick_params(axis="x", labelsize=8)
    ax_supp_c.grid(axis="x", color="0.90")
    ax_supp_c.set_axisbelow(True)
    ax_supp_c.spines["top"].set_visible(False)
    ax_supp_c.spines["right"].set_visible(False)
    ax_supp_c.text(-0.18, 1.06, "A", transform=ax_supp_c.transAxes, fontsize=13, fontweight="bold")

    left = np.zeros(len(REGIONS))
    for level in pi_levels:
        vals = np.array([pi_props[region][level] for region in REGIONS])
        ax_supp_d.barh(y, vals, left=left, color=PI_COLORS[level], edgecolor="white", linewidth=0.5, label=level)
        left += vals
    ax_supp_d.set_yticks(y)
    ax_supp_d.set_yticklabels([display_region(region) for region in REGIONS], fontsize=8)
    ax_supp_d.invert_yaxis()
    ax_supp_d.set_xlim(0, 100)
    ax_supp_d.xaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax_supp_d.set_xlabel("% of upland site-months", fontsize=8)
    ax_supp_d.set_title("Prediction-Interval Support", fontweight="bold", fontsize=9)
    ax_supp_d.tick_params(axis="x", labelsize=8)
    ax_supp_d.grid(axis="x", color="0.90")
    ax_supp_d.set_axisbelow(True)
    ax_supp_d.spines["top"].set_visible(False)
    ax_supp_d.spines["right"].set_visible(False)
    ax_supp_d.text(-0.18, 1.06, "B", transform=ax_supp_d.transAxes, fontsize=13, fontweight="bold")
    class_handles, class_labels = ax_supp_c.get_legend_handles_labels()
    pi_handles, pi_labels = ax_supp_d.get_legend_handles_labels()
    fig_supp.legend(
        class_handles, class_labels,
        loc="lower left",
        bbox_to_anchor=(0.12, 0.02),
        ncol=1,
        fontsize=7,
        frameon=False,
    )
    fig_supp.legend(
        pi_handles, pi_labels,
        loc="lower left",
        bbox_to_anchor=(0.60, 0.045),
        ncol=1,
        fontsize=7,
        frameon=False,
    )
    fig_supp.suptitle("Q3-Aligned Upland EDI Classes and Uncertainty", fontsize=11.5, fontweight="bold", y=0.99)
    fig_supp.tight_layout(rect=[0, 0.30, 1, 0.94])
    fig_supp.savefig(OUT_SUPP, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_supp)

    print(OUT)
    print(OUT_MAIN)
    print(OUT_SUPP)


if __name__ == "__main__":
    main()