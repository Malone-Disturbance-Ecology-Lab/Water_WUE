"""
24-EDI_spatial.py
==================
Apply the Ecological Drought Index (EDI) response model to spatially-masked
SPEI-3 rasters and summarise results by US coast region.

Background
----------
The primary response-based EDI uses the Q3 coast-threshold GAM output exported
by Ecological_Impacts.R:

  Malone_Workflow/results/Ecological_Impacts/
  EDI_Q3_aligned_prediction_scores_upland.csv

That table is copied from the final Q3 workflow and gives coast-specific upland
WUE_T percent change from near-normal SPEI. Freshwater and Saline response
curves are retained only as legacy scenario diagnostics because the final Q3
manuscript logic is upland-only and coast-specific.

Primary Q3 response classes:
  No meaningful change < 5% absolute predicted WUE_T change
  Watch                5-10%
  Stress               10-20%
  Impact               >= 20%

Spatial response_v2 event maps use the 5% threshold to separate decrease and
increase events, including dry-only maps and prediction-interval support maps.

The older pooled logistic EDI is retained as a secondary exposure diagnostic.

Coast regions (consistent with Q3 / Q4 workflow):
  AK Coast       lat > 50
  Pacific Coast  lat <= 50 AND lon < -120
  Gulf Coast     lat <= 50 AND -120 <= lon <= -100
  Atlantic Coast lat <= 50 AND lon > -100

Inputs
------
  Masked SPEI-3 files:
    M:/Research/WUE_CUE/spatial_SPEI/SPEI-3/
    SPEI3_*_US_OCEAN_COAST_80km.nc

Outputs
-------
  M:/Research/WUE_CUE/WUE_manuscript_version6/upscaling/
    EDI_response_monthly_coast_summary.csv — monthly response classes by region
    EDI_response_annual_coast_summary.csv  — annual response summaries by region
    EDI_response_pixel_mean_abs_change_*.nc — mean absolute WUE_T change raster
    EDI_response_pixel_freq_strong_*.nc     — frequency (%) of strong change
    EDI_logistic_*                          — secondary pooled logistic outputs
    EDI_response_communication_annual_impacted_area_by_coast_2000_2025.csv
    EDI_response_communication_plot_01_event_frequency_maps_2000_2025.png
    EDI_response_communication_plot_02_annual_impacted_area_by_coast_2000_2025.png
    EDI_plot_01_mean_map.png            — map of time-mean EDI
    EDI_plot_02_impact_freq_map.png     — map of % months in Impact class
    EDI_plot_03_timeseries_by_coast.png — monthly mean EDI time series by region
    EDI_plot_04_severity_stack.png      — stacked area: severity proportions over time
    EDI_plot_05_annual_summary.png      — annual summaries by coast region
    EDI_plot_06_composite.png           — composite figure
"""

import os, sys, glob, re
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.ticker import PercentFormatter
from pathlib import Path
from datetime import datetime
import csv

# ---------------------------------------------------------------------------
# PATHS  — using manuscript version6 directories (Python/manuscript-only)
# ---------------------------------------------------------------------------
SPEI3_DIR   = r"M:\Research\WUE_CUE\spatial_SPEI\SPEI-3"
OUTPUT_DIR  = r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling"
Q3_DIR      = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs"
Q4_DIR      = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q4\Q4_Ecological_Impacts_outputs"

SPATIAL_START_YEAR = int(os.environ.get("EDI_START_YEAR", "2000"))
SPATIAL_END_YEAR = int(os.environ.get("EDI_END_YEAR", "2025"))
EVENT_THRESHOLD_PCT = float(os.environ.get("EDI_EVENT_THRESHOLD_PCT", "5.0"))

# Required response curves from Q4 directory
RESPONSE_CURVES = os.path.join(Q4_DIR, "EDI_response_prediction_curves.csv")
Q3_ALIGNED_CURVES = os.path.join(Q4_DIR, "EDI_Q3_aligned_prediction_scores_upland.csv")

# Required slope summary from Q3 directory
Q3_SLOPE_SUMMARY = os.path.join(Q3_DIR, "Q3_WUE_T_SPEI_coast_region_slope_summary.csv")

# Optional support files (will be skipped if missing)
Q3_THRESHOLD_SUMMARY = os.path.join(Q3_DIR, "Q3_WUE_T_SPEI_coast_predicted_threshold_summary_5_10_20pct.csv")
Q4_DECLINE_THRESHOLDS = os.path.join(Q4_DIR, "Q4_WUE_decline_5pct_probability_thresholds_by_group.csv")

# ---------------------------------------------------------------------------
# EDI MODEL COEFFICIENTS  (from Ecological_Impacts.R — pooled model)
# ---------------------------------------------------------------------------
B0_POOLED = 0.1457163    # intercept
B1_POOLED = 0.0984531    # slope on SPEI_3

FILL       = 9.969209968386869e+36

# EDI thresholds
EDI_WATCH  = 0.30
EDI_STRESS = 0.50
EDI_IMPACT = 0.70

# ---------------------------------------------------------------------------
# COAST REGION DEFINITIONS  (matches Q3 classify_coast_region)
# ---------------------------------------------------------------------------
COAST_REGIONS = {
    "AK Coast":       lambda lat, lon: lat > 50,
    "Pacific Coast":  lambda lat, lon: (lat <= 50) & (lon < -120),
    "Gulf Coast":     lambda lat, lon: (lat <= 50) & (lon >= -120) & (lon <= -100),
    "Atlantic Coast": lambda lat, lon: (lat <= 50) & (lon > -100),
}
COAST_COLORS = {
    "AK Coast":       "#2B6CB0",
    "Pacific Coast":  "#F6AD55",
    "Gulf Coast":     "#6A3D9A",
    "Atlantic Coast": "#2B9348",
}
COAST_DISPLAY_NAMES = {
    "AK Coast":       "Alaska Coast",
    "Pacific Coast":  "Pacific Coast",
    "Gulf Coast":     "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast",
}
SEVERITY_COLORS = {
    "None":   "#2166AC",
    "Watch":  "#92C5DE",
    "Stress": "#F4A582",
    "Impact": "#D6604D",
}
RESPONSE_CLASS_COLORS = {
    "Stable":   "#2166AC",
    "Mild":     "#67A9CF",
    "Moderate": "#FDDC8A",
    "Large":    "#F4A582",
    "Strong":   "#B2182B",
}
RESPONSE_CLASSES = ["Stable", "Mild", "Moderate", "Large", "Strong"]
RESPONSE_CLASS_IDS = {name: i for i, name in enumerate(RESPONSE_CLASSES)}
WATER_COLORS = {
    "Upland": "#7570B3",
    "Freshwater": "#1B9E77",
    "Saline": "#D95F02",
}

# ---------------------------------------------------------------------------
def logistic(x, b0, b1):
    return 1.0 / (1.0 + np.exp(-(b0 + b1 * x)))


def classify_edi(edi):
    """Return integer class: 0=None, 1=Watch, 2=Stress, 3=Impact."""
    c = np.full(edi.shape, -1, dtype=np.int8)
    c[edi <  EDI_WATCH]  = 0
    c[(edi >= EDI_WATCH)  & (edi < EDI_STRESS)] = 1
    c[(edi >= EDI_STRESS) & (edi < EDI_IMPACT)] = 2
    c[edi >= EDI_IMPACT] = 3
    return c


def classify_response(abs_change):
    """Return integer class for absolute predicted WUE_T percent change."""
    c = np.full(abs_change.shape, -1, dtype=np.int8)
    c[abs_change < 0.5] = RESPONSE_CLASS_IDS["Stable"]
    c[(abs_change >= 0.5) & (abs_change < 1.0)] = RESPONSE_CLASS_IDS["Mild"]
    c[(abs_change >= 1.0) & (abs_change < 2.5)] = RESPONSE_CLASS_IDS["Moderate"]
    c[(abs_change >= 2.5) & (abs_change < 5.0)] = RESPONSE_CLASS_IDS["Large"]
    c[abs_change >= 5.0] = RESPONSE_CLASS_IDS["Strong"]
    return c


def load_response_curves(path):
    curves = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            water_class = row["water_class"]
            curves.setdefault(water_class, {"spei": [], "pct": [], "lower": [], "upper": []})
            curves[water_class]["spei"].append(float(row["SPEI_3"]))
            curves[water_class]["pct"].append(float(row["predicted_pct_change"]))
            curves[water_class]["lower"].append(float(row["predicted_pct_change_lower"]))
            curves[water_class]["upper"].append(float(row["predicted_pct_change_upper"]))

    out = {}
    for water_class, vals in curves.items():
        order = np.argsort(vals["spei"])
        spei = np.array(vals["spei"], dtype=np.float64)[order]
        pct = np.array(vals["pct"], dtype=np.float64)[order]
        lower = np.array(vals["lower"], dtype=np.float64)[order]
        upper = np.array(vals["upper"], dtype=np.float64)[order]
        uniq, idx = np.unique(spei, return_index=True)
        out[water_class] = {
            "spei": uniq,
            "pct": pct[idx],
            "lower": lower[idx],
            "upper": upper[idx],
        }
    return out


def load_q3_aligned_upland_curves(path, spei_timescale="SPEI_3"):
    curves = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("SPEI_timescale") != spei_timescale:
                continue
            if row.get("water_class") != "Upland":
                continue
            coast_region = row["coast_region"]
            curves.setdefault(coast_region, {"spei": [], "pct": [], "lower": [], "upper": []})
            curves[coast_region]["spei"].append(float(row["SPEI_value"]))
            curves[coast_region]["pct"].append(float(row["predicted_pct_change"]))
            curves[coast_region]["lower"].append(float(row["predicted_lower_pct"]))
            curves[coast_region]["upper"].append(float(row["predicted_upper_pct"]))

    out = {}
    for coast_region, vals in curves.items():
        order = np.argsort(vals["spei"])
        spei = np.array(vals["spei"], dtype=np.float64)[order]
        pct = np.array(vals["pct"], dtype=np.float64)[order]
        lower = np.array(vals["lower"], dtype=np.float64)[order]
        upper = np.array(vals["upper"], dtype=np.float64)[order]
        uniq, idx = np.unique(spei, return_index=True)
        out[coast_region] = {
            "spei": uniq,
            "pct": pct[idx],
            "lower": lower[idx],
            "upper": upper[idx],
        }
    return out


def predict_response_pct(spei, curve, field="pct"):
    spei_x = curve["spei"]
    pct_y = curve[field]
    flat = spei.ravel()
    pred = np.interp(flat, spei_x, pct_y, left=pct_y[0], right=pct_y[-1])
    pred = pred.reshape(spei.shape)
    pred[~np.isfinite(spei)] = np.nan
    return pred


def predict_q3_coast_response_pct(spei, coast_id, region_names, curves, field="pct"):
    pred = np.full(spei.shape, np.nan, dtype=np.float64)
    for rid, name in enumerate(region_names):
        if name not in curves:
            continue
        pix = (coast_id == rid) & np.isfinite(spei)
        if not np.any(pix):
            continue
        pred[pix] = np.interp(
            spei[pix],
            curves[name]["spei"],
            curves[name][field],
            left=curves[name][field][0],
            right=curves[name][field][-1],
        )
    return pred


def safe_nanmean(vals):
    vals = np.asarray(vals, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float(np.mean(vals))


def pct_true(vals):
    vals = np.asarray(vals)
    if vals.size == 0:
        return np.nan
    return float(np.mean(vals) * 100)


def read_csv_lookup(path, key_fields, filters=None):
    filters = filters or {}
    out = {}
    if not Path(path).exists():
        print(f"  Warning: Optional file not found, skipping: {path}")
        return out
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            keep = True
            for k, v in filters.items():
                if row.get(k) != v:
                    keep = False
                    break
            if not keep:
                continue
            key = tuple(row.get(k, "") for k in key_fields)
            out[key] = row
    return out


def num_or_nan(row, field):
    if not row:
        return np.nan
    val = row.get(field, "")
    if val in ("", "NA", "NaN", "nan", None):
        return np.nan
    try:
        return float(val)
    except ValueError:
        return np.nan


def extract_yyyymm(path):
    m = re.search(r"(\d{6})", os.path.basename(path))
    return int(m.group(1)) if m else 0


def write_csv(rows, header, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def display_coast(name):
    return COAST_DISPLAY_NAMES.get(name, name)


# ---------------------------------------------------------------------------
def main():
    print("=" * 68)
    print("24-EDI_spatial.py — Spatial EDI from SPEI-3")
    print(f"Started: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 68)

    # Create output directory
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    
    # Check if input directories exist
    if not Path(SPEI3_DIR).exists():
        sys.exit(f"SPEI-3 directory not found: {SPEI3_DIR}")
    if not Path(Q3_DIR).exists():
        print(f"Warning: Q3 directory not found: {Q3_DIR}")
        print("Continuing without Q3 validation files...")
    if not Path(Q4_DIR).exists():
        sys.exit(f"Q4 directory not found: {Q4_DIR}")
    
    # Load required response curves
    try:
        response_curves = load_response_curves(RESPONSE_CURVES)
    except FileNotFoundError:
        sys.exit(f"Required file not found: {RESPONSE_CURVES}\n"
                 f"Please ensure EDI_response_prediction_curves.csv exists in {Q4_DIR}")
    
    try:
        q3_upland_curves = load_q3_aligned_upland_curves(Q3_ALIGNED_CURVES, spei_timescale="SPEI_3")
    except FileNotFoundError:
        sys.exit(f"Required file not found: {Q3_ALIGNED_CURVES}\n"
                 f"Please ensure EDI_Q3_aligned_prediction_scores_upland.csv exists in {Q4_DIR}")
    
    water_classes = list(response_curves.keys())
    print("Loaded legacy ecosystem response curves:")
    for water_class in water_classes:
        spei_x = response_curves[water_class]["spei"]
        pct_y = response_curves[water_class]["pct"]
        print(
            f"  {water_class}: SPEI {spei_x.min():.2f} to {spei_x.max():.2f}; "
            f"WUE_T change {pct_y.min():.2f}% to {pct_y.max():.2f}%"
        )
    print("Loaded Q3-aligned coast-specific upland SPEI-3 curves:")
    for coast_region, curve in q3_upland_curves.items():
        print(
            f"  {coast_region}: SPEI {curve['spei'].min():.2f} to {curve['spei'].max():.2f}; "
            f"WUE_T change {curve['pct'].min():.2f}% to {curve['pct'].max():.2f}%"
        )

    # ---- Find input files ---------------------------------------------------
    files = sorted(glob.glob(os.path.join(SPEI3_DIR, "SPEI3_*_US_OCEAN_COAST_80km.nc")),
                   key=extract_yyyymm)
    if not files:
        sys.exit(f"No masked SPEI-3 files found in {SPEI3_DIR}\n"
                 "Run 23-mask_spei3_coastal.py first.")
    # Filter to keep only April-October (months 4-10) for the specified year range
    files = [
        fp for fp in files
        if SPATIAL_START_YEAR <= extract_yyyymm(fp) // 100 <= SPATIAL_END_YEAR
        and 4 <= extract_yyyymm(fp) % 100 <= 10
    ]
    if not files:
        sys.exit(
            f"No masked SPEI-3 files found for {SPATIAL_START_YEAR}-{SPATIAL_END_YEAR}, "
            f"April-October months in {SPEI3_DIR}"
        )
    preview_n = os.environ.get("EDI_PREVIEW_N")
    if preview_n:
        files = files[:int(preview_n)]
        print(f"Preview mode: using first {len(files)} masked SPEI-3 files")
    print(
        f"Found {len(files)} masked SPEI-3 files "
        f"for {SPATIAL_START_YEAR}-{SPATIAL_END_YEAR} (April-October only)"
    )

    # ---- Read grid from first file ------------------------------------------
    ref   = nc.Dataset(files[0])
    lat1d = np.array(ref.variables["lat"][:])
    lon1d = np.array(ref.variables["lon"][:])
    ref.close()
    LAT, LON = np.meshgrid(lat1d, lon1d, indexing="ij")  # (nlat, nlon)
    nlat, nlon = LAT.shape

    # ---- Build coast-region mask array (nlat, nlon) -------------------------
    coast_id = np.full((nlat, nlon), -1, dtype=np.int8)
    region_names = list(COAST_REGIONS.keys())
    for rid, (name, fn) in enumerate(COAST_REGIONS.items()):
        coast_id[fn(LAT, LON)] = rid
    print(f"Coast region pixel counts:")
    for rid, name in enumerate(region_names):
        print(f"  {name}: {(coast_id == rid).sum():,}")

    # ---- Process all files --------------------------------------------------
    n_files      = len(files)
    edi_stack    = np.full((n_files, nlat, nlon), np.nan, dtype=np.float32)
    response_abs_stacks = {
        water_class: np.full((n_files, nlat, nlon), np.nan, dtype=np.float32)
        for water_class in water_classes
    }
    response_pct_stacks = {
        water_class: np.full((n_files, nlat, nlon), np.nan, dtype=np.float32)
        for water_class in water_classes
    }
    valid_mask   = np.zeros((nlat, nlon), dtype=bool)
    monthly_rows = []
    response_monthly_rows = []
    response_v2_monthly_rows = []

    # Upland-focused response_v2 accumulators. These avoid keeping several
    # extra full time stacks while still allowing pixel-level output rasters.
    v2_count = np.zeros((nlat, nlon), dtype=np.float32)
    v2_sum_pct = np.zeros((nlat, nlon), dtype=np.float64)
    v2_sum_abs = np.zeros((nlat, nlon), dtype=np.float64)
    v2_sum_lower = np.zeros((nlat, nlon), dtype=np.float64)
    v2_sum_upper = np.zeros((nlat, nlon), dtype=np.float64)
    v2_sum_ci_width = np.zeros((nlat, nlon), dtype=np.float64)
    v2_count_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_strong_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_strong_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_possible_strong = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_robust_strong = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_possible_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_robust_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_possible_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_count_robust_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_sum_pct = np.zeros((nlat, nlon), dtype=np.float64)
    v2_dry_sum_abs = np.zeros((nlat, nlon), dtype=np.float64)
    v2_dry_count_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_strong_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_strong_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_possible_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_robust_decrease = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_possible_increase = np.zeros((nlat, nlon), dtype=np.float32)
    v2_dry_count_robust_increase = np.zeros((nlat, nlon), dtype=np.float32)

    print(f"\nProcessing {n_files} files ...")
    for i, fp in enumerate(files):
        yyyymm = extract_yyyymm(fp)
        year   = yyyymm // 100
        month  = yyyymm  % 100

        f     = nc.Dataset(fp)
        spei  = np.array(f.variables["SPEI3"][0, :, :], dtype=np.float64)
        f.close()

        # Mask fill / invalid
        bad = (spei > 1e30) | ~np.isfinite(spei)
        spei[bad] = np.nan

        # Compute EDI
        edi = logistic(spei, B0_POOLED, B1_POOLED)
        edi[bad] = np.nan

        edi_stack[i] = edi.astype(np.float32)
        valid_mask |= ~bad

        # Classify
        cls = classify_edi(edi)

        # Per-region summaries
        for rid, name in enumerate(region_names):
            pix = (coast_id == rid) & ~bad
            n   = pix.sum()
            if n == 0:
                continue
            mean_spei = float(np.nanmean(spei[pix]))
            mean_edi  = float(np.nanmean(edi[pix]))
            p_none    = float((cls[pix] == 0).sum() / n * 100)
            p_watch   = float((cls[pix] == 1).sum() / n * 100)
            p_stress  = float((cls[pix] == 2).sum() / n * 100)
            p_impact  = float((cls[pix] == 3).sum() / n * 100)
            monthly_rows.append([year, month, name, n,
                                  round(mean_spei, 4), round(mean_edi, 4),
                                  round(p_none, 2), round(p_watch, 2),
                                  round(p_stress, 2), round(p_impact, 2)])

        # Response-based EDI scenarios by ecosystem class.
        for water_class in water_classes:
            if water_class == "Upland":
                pct_change = predict_q3_coast_response_pct(
                    spei, coast_id, region_names, q3_upland_curves, field="pct"
                )
            else:
                pct_change = predict_response_pct(spei, response_curves[water_class])
            abs_change = np.abs(pct_change)
            response_pct_stacks[water_class][i] = pct_change.astype(np.float32)
            response_abs_stacks[water_class][i] = abs_change.astype(np.float32)
            response_cls = classify_response(abs_change)

            if water_class == "Upland":
                lower = predict_q3_coast_response_pct(
                    spei, coast_id, region_names, q3_upland_curves, field="lower"
                )
                upper = predict_q3_coast_response_pct(
                    spei, coast_id, region_names, q3_upland_curves, field="upper"
                )
                good = ~bad
                v2_count[good] += 1
                v2_sum_pct[good] += pct_change[good]
                v2_sum_abs[good] += abs_change[good]
                v2_sum_lower[good] += lower[good]
                v2_sum_upper[good] += upper[good]
                v2_sum_ci_width[good] += (upper[good] - lower[good])
                v2_count_decrease[good] += pct_change[good] < 0
                v2_count_increase[good] += pct_change[good] > 0
                v2_count_strong_decrease[good] += pct_change[good] <= -5.0
                v2_count_strong_increase[good] += pct_change[good] >= 5.0
                v2_count_possible_strong[good] += (
                    (np.abs(lower[good]) >= 5.0) | (np.abs(upper[good]) >= 5.0)
                )
                v2_count_robust_strong[good] += (
                    (np.minimum(np.abs(lower[good]), np.abs(upper[good])) >= 5.0) &
                    (np.sign(lower[good]) == np.sign(upper[good]))
                )
                v2_count_possible_decrease[good] += lower[good] <= -5.0
                v2_count_robust_decrease[good] += upper[good] <= -5.0
                v2_count_possible_increase[good] += upper[good] >= 5.0
                v2_count_robust_increase[good] += lower[good] >= 5.0

                dry = good & (spei <= -1.0)
                v2_dry_count[dry] += 1
                v2_dry_sum_pct[dry] += pct_change[dry]
                v2_dry_sum_abs[dry] += abs_change[dry]
                v2_dry_count_decrease[dry] += pct_change[dry] < 0
                v2_dry_count_increase[dry] += pct_change[dry] > 0
                v2_dry_count_strong_decrease[dry] += pct_change[dry] <= -5.0
                v2_dry_count_strong_increase[dry] += pct_change[dry] >= 5.0
                v2_dry_count_possible_decrease[dry] += lower[dry] <= -5.0
                v2_dry_count_robust_decrease[dry] += upper[dry] <= -5.0
                v2_dry_count_possible_increase[dry] += upper[dry] >= 5.0
                v2_dry_count_robust_increase[dry] += lower[dry] >= 5.0
            else:
                lower = upper = None

            for rid, name in enumerate(region_names):
                pix = (coast_id == rid) & ~bad
                n = pix.sum()
                if n == 0:
                    continue
                vals = pct_change[pix]
                abs_vals = abs_change[pix]
                response_monthly_rows.append([
                    year, month, name, water_class, n,
                    round(float(np.nanmean(spei[pix])), 4),
                    round(float(np.nanmean(vals)), 4),
                    round(float(np.nanmean(abs_vals)), 4),
                    round(float(np.nanmean(vals > 0) * 100), 2),
                    round(float(np.nanmean(vals < 0) * 100), 2),
                    round(float((response_cls[pix] == RESPONSE_CLASS_IDS["Stable"]).sum() / n * 100), 2),
                    round(float((response_cls[pix] == RESPONSE_CLASS_IDS["Mild"]).sum() / n * 100), 2),
                    round(float((response_cls[pix] == RESPONSE_CLASS_IDS["Moderate"]).sum() / n * 100), 2),
                    round(float((response_cls[pix] == RESPONSE_CLASS_IDS["Large"]).sum() / n * 100), 2),
                    round(float((response_cls[pix] == RESPONSE_CLASS_IDS["Strong"]).sum() / n * 100), 2),
                ])

                if water_class == "Upland":
                    dry_pix = pix & (spei <= -1.0)
                    dry_n = int(dry_pix.sum())
                    vals = pct_change[pix]
                    lower_vals = lower[pix]
                    upper_vals = upper[pix]
                    ci_width = upper_vals - lower_vals
                    dry_vals = pct_change[dry_pix] if dry_n else np.array([])

                    response_v2_monthly_rows.append([
                        year, month, name, n,
                        round(float(np.nanmean(spei[pix])), 4),
                        round(float(np.nanmean(vals)), 4),
                        round(float(np.nanmean(abs_vals)), 4),
                        round(safe_nanmean(lower_vals), 4),
                        round(safe_nanmean(upper_vals), 4),
                        round(safe_nanmean(ci_width), 4),
                        round(pct_true(vals < 0), 2),
                        round(pct_true(vals > 0), 2),
                        round(pct_true(vals <= -5.0), 2),
                        round(pct_true(vals >= 5.0), 2),
                        round(pct_true((np.minimum(np.abs(lower_vals), np.abs(upper_vals)) >= 5.0) &
                                       (np.sign(lower_vals) == np.sign(upper_vals))), 2),
                        round(pct_true((np.abs(lower_vals) >= 5.0) | (np.abs(upper_vals) >= 5.0)), 2),
                        dry_n,
                        round(float(dry_n / n * 100), 2),
                        round(safe_nanmean(dry_vals), 4),
                        round(safe_nanmean(np.abs(dry_vals)), 4),
                        round(pct_true(dry_vals < 0), 2) if dry_n else np.nan,
                        round(pct_true(dry_vals > 0), 2) if dry_n else np.nan,
                        round(pct_true(dry_vals <= -5.0), 2) if dry_n else np.nan,
                        round(pct_true(dry_vals >= 5.0), 2) if dry_n else np.nan,
                    ])

        if (i + 1) % 50 == 0 or (i + 1) == n_files:
            print(f"  {i+1}/{n_files} processed")

    # ---- Time-mean EDI raster -----------------------------------------------
    mean_edi    = np.nanmean(edi_stack, axis=0)
    freq_impact = np.nanmean(edi_stack >= EDI_IMPACT, axis=0) * 100  # % months
    mean_abs_response = {
        water_class: np.nanmean(stack, axis=0)
        for water_class, stack in response_abs_stacks.items()
    }
    mean_pct_response = {
        water_class: np.nanmean(stack, axis=0)
        for water_class, stack in response_pct_stacks.items()
    }
    freq_strong_response = {
        water_class: np.nanmean(stack >= 5.0, axis=0) * 100
        for water_class, stack in response_abs_stacks.items()
    }
    with np.errstate(invalid="ignore", divide="ignore"):
        response_v2_rasters = {
            "mean_pct_change": v2_sum_pct / v2_count,
            "mean_abs_change": v2_sum_abs / v2_count,
            "mean_lower_pct_change": v2_sum_lower / v2_count,
            "mean_upper_pct_change": v2_sum_upper / v2_count,
            "mean_ci_width": v2_sum_ci_width / v2_count,
            "freq_decrease": v2_count_decrease / v2_count * 100,
            "freq_increase": v2_count_increase / v2_count * 100,
            "freq_strong_decrease": v2_count_strong_decrease / v2_count * 100,
            "freq_strong_increase": v2_count_strong_increase / v2_count * 100,
            "freq_possible_strong": v2_count_possible_strong / v2_count * 100,
            "freq_robust_strong": v2_count_robust_strong / v2_count * 100,
            "freq_possible_decrease": v2_count_possible_decrease / v2_count * 100,
            "freq_robust_decrease": v2_count_robust_decrease / v2_count * 100,
            "freq_possible_increase": v2_count_possible_increase / v2_count * 100,
            "freq_robust_increase": v2_count_robust_increase / v2_count * 100,
            "dry_month_frequency": v2_dry_count / v2_count * 100,
            "dry_mean_pct_change": v2_dry_sum_pct / v2_dry_count,
            "dry_mean_abs_change": v2_dry_sum_abs / v2_dry_count,
            "dry_freq_decrease": v2_dry_count_decrease / v2_dry_count * 100,
            "dry_freq_increase": v2_dry_count_increase / v2_dry_count * 100,
            "dry_freq_strong_decrease": v2_dry_count_strong_decrease / v2_dry_count * 100,
            "dry_freq_strong_increase": v2_dry_count_strong_increase / v2_dry_count * 100,
            "dry_freq_possible_decrease": v2_dry_count_possible_decrease / v2_dry_count * 100,
            "dry_freq_robust_decrease": v2_dry_count_robust_decrease / v2_dry_count * 100,
            "dry_freq_possible_increase": v2_dry_count_possible_increase / v2_dry_count * 100,
            "dry_freq_robust_increase": v2_dry_count_robust_increase / v2_dry_count * 100,
        }
        response_v2_rasters["drought_adjustment_index"] = response_v2_rasters["dry_freq_strong_increase"]
        response_v2_rasters["drought_impairment_index"] = response_v2_rasters["dry_freq_strong_decrease"]
        response_v2_rasters["drought_net_adjustment_minus_impairment"] = (
            response_v2_rasters["drought_adjustment_index"] -
            response_v2_rasters["drought_impairment_index"]
        )
        response_v2_rasters["drought_possible_adjustment_index"] = response_v2_rasters["dry_freq_possible_increase"]
        response_v2_rasters["drought_robust_adjustment_index"] = response_v2_rasters["dry_freq_robust_increase"]
        response_v2_rasters["drought_possible_impairment_index"] = response_v2_rasters["dry_freq_possible_decrease"]
        response_v2_rasters["drought_robust_impairment_index"] = response_v2_rasters["dry_freq_robust_decrease"]

    # ---- Annual summary ------------------------------------------------------
    annual_rows = []
    years = sorted(set(r[0] for r in monthly_rows))
    for yr in years:
        yr_rows = [r for r in monthly_rows if r[0] == yr]
        for rid, name in enumerate(region_names):
            sub = [r for r in yr_rows if r[2] == name]
            if not sub:
                continue
            annual_rows.append([
                yr, name,
                round(np.mean([r[5] for r in sub]), 4),   # mean EDI
                round(np.mean([r[4] for r in sub]), 4),   # mean SPEI
                round(np.mean([r[6] for r in sub]), 2),   # % None
                round(np.mean([r[7] for r in sub]), 2),   # % Watch
                round(np.mean([r[8] for r in sub]), 2),   # % Stress
                round(np.mean([r[9] for r in sub]), 2),   # % Impact
            ])

    response_annual_rows = []
    for yr in years:
        yr_rows = [r for r in response_monthly_rows if r[0] == yr]
        for name in region_names:
            for water_class in water_classes:
                sub = [r for r in yr_rows if r[2] == name and r[3] == water_class]
                if not sub:
                    continue
                response_annual_rows.append([
                    yr, name, water_class,
                    round(np.mean([r[6] for r in sub]), 4),   # mean pct change
                    round(np.mean([r[7] for r in sub]), 4),   # mean abs change
                    round(np.mean([r[8] for r in sub]), 2),   # % increase
                    round(np.mean([r[9] for r in sub]), 2),   # % decrease
                    round(np.mean([r[10] for r in sub]), 2),  # % Stable
                    round(np.mean([r[11] for r in sub]), 2),  # % Mild
                    round(np.mean([r[12] for r in sub]), 2),  # % Moderate
                    round(np.mean([r[13] for r in sub]), 2),  # % Large
                    round(np.mean([r[14] for r in sub]), 2),  # % Strong
                ])

    response_v2_annual_rows = []
    for yr in years:
        yr_rows = [r for r in response_v2_monthly_rows if r[0] == yr]
        for name in region_names:
            sub = [r for r in yr_rows if r[2] == name]
            if not sub:
                continue
            response_v2_annual_rows.append([
                yr, name,
                round(np.mean([r[4] for r in sub]), 4),
                round(np.mean([r[5] for r in sub]), 4),
                round(np.mean([r[6] for r in sub]), 4),
                round(np.mean([r[7] for r in sub]), 4),
                round(np.mean([r[8] for r in sub]), 4),
                round(np.mean([r[9] for r in sub]), 4),
                round(np.mean([r[10] for r in sub]), 2),
                round(np.mean([r[11] for r in sub]), 2),
                round(np.mean([r[12] for r in sub]), 2),
                round(np.mean([r[13] for r in sub]), 2),
                round(np.mean([r[14] for r in sub]), 2),
                round(np.mean([r[15] for r in sub]), 2),
                round(np.mean([r[17] for r in sub]), 2),
                round(safe_nanmean([r[18] for r in sub]), 4),
                round(safe_nanmean([r[19] for r in sub]), 4),
                round(safe_nanmean([r[20] for r in sub]), 2),
                round(safe_nanmean([r[21] for r in sub]), 2),
                round(safe_nanmean([r[22] for r in sub]), 2),
                round(safe_nanmean([r[23] for r in sub]), 2),
            ])

    # Load Q3 slopes (required for validation summary)
    q3_slopes = read_csv_lookup(
        Q3_SLOPE_SUMMARY, ["coast_region"], filters={"SPEI_timescale": "SPEI_3"}
    )
    
    # Load optional Q3 thresholds and Q4 decline thresholds (will skip if missing)
    q3_thresholds = read_csv_lookup(
        Q3_THRESHOLD_SUMMARY, ["coast_region"], filters={"SPEI_timescale": "SPEI_3"}
    )
    q4_thresholds = read_csv_lookup(
        Q4_DECLINE_THRESHOLDS, ["coast_region"],
        filters={"metric": "WUE_T", "SPEI_timescale": "SPEI_3", "water_class": "Upland"}
    )
    
    response_v2_validation_rows = []
    response_v3_rows = []
    for name in region_names:
        sub = [r for r in response_v2_annual_rows if r[1] == name]
        if not sub:
            continue
        q3s = q3_slopes.get((name,), {})
        q3t = q3_thresholds.get((name,), {})
        q4t = q4_thresholds.get((name,), {})
        response_v2_validation_rows.append([
            name,
            round(np.mean([r[3] for r in sub]), 4),
            round(np.mean([r[4] for r in sub]), 4),
            round(np.mean([r[8] for r in sub]), 2),
            round(np.mean([r[9] for r in sub]), 2),
            round(np.mean([r[10] for r in sub]), 2),
            round(np.mean([r[11] for r in sub]), 2),
            round(np.mean([r[14] for r in sub]), 2),
            round(safe_nanmean([r[15] for r in sub]), 4),
            round(safe_nanmean([r[16] for r in sub]), 4),
            round(safe_nanmean([r[19] for r in sub]), 2),
            round(safe_nanmean([r[20] for r in sub]), 2),
            num_or_nan(q3s, "n_sites"),
            num_or_nan(q3s, "mean_slope"),
            num_or_nan(q3s, "median_slope"),
            num_or_nan(q3s, "pct_negative"),
            num_or_nan(q3t, "SPEI_threshold_dry_increase_5%"),
            num_or_nan(q3t, "SPEI_threshold_dry_decrease_5%"),
            num_or_nan(q3t, "SPEI_threshold_wet_increase_5%"),
            num_or_nan(q3t, "SPEI_threshold_wet_decrease_5%"),
            num_or_nan(q4t, "dry_SPEI_prob_25"),
            num_or_nan(q4t, "dry_SPEI_prob_50"),
            num_or_nan(q4t, "dry_SPEI_prob_75"),
            num_or_nan(q4t, "min_decline_probability_5"),
            num_or_nan(q4t, "max_decline_probability_5"),
        ])
        dry_adjustment = round(safe_nanmean([r[20] for r in sub]), 2)
        dry_impairment = round(safe_nanmean([r[19] for r in sub]), 2)
        net_adjustment = round(dry_adjustment - dry_impairment, 2)
        q4_max_decline = num_or_nan(q4t, "max_decline_probability_5")
        if dry_impairment >= 25:
            interpretation = "dry_impairment_signal"
        elif dry_adjustment >= 25 and dry_impairment < 5:
            interpretation = "dry_adjustment_signal_not_impairment"
        elif dry_adjustment >= 25 and dry_impairment >= 5:
            interpretation = "mixed_dry_adjustment_and_impairment"
        else:
            interpretation = "weak_dry_response"
        response_v3_rows.append([
            name,
            round(np.mean([r[14] for r in sub]), 2),
            round(safe_nanmean([r[15] for r in sub]), 4),
            round(safe_nanmean([r[16] for r in sub]), 4),
            dry_adjustment,
            dry_impairment,
            net_adjustment,
            q4_max_decline,
            num_or_nan(q3s, "mean_slope"),
            num_or_nan(q3s, "pct_negative"),
            interpretation,
        ])

    communication_annual_rows = []
    for r in response_v2_annual_rows:
        communication_annual_rows.append([
            r[0],
            r[1],
            EVENT_THRESHOLD_PCT,
            r[10],
            r[11],
            r[10] + r[11],
            r[11] - r[10],
            r[14],
            r[19],
            r[20],
        ])

    # ---- Write CSVs ---------------------------------------------------------
    write_csv(monthly_rows,
              ["year","month","coast_region","n_pixels",
               "mean_SPEI3","mean_EDI",
               "pct_None","pct_Watch","pct_Stress","pct_Impact"],
              os.path.join(OUTPUT_DIR, "EDI_logistic_monthly_coast_summary.csv"))

    write_csv(annual_rows,
              ["year","coast_region","mean_EDI","mean_SPEI3",
               "pct_None","pct_Watch","pct_Stress","pct_Impact"],
              os.path.join(OUTPUT_DIR, "EDI_logistic_annual_coast_summary.csv"))

    write_csv(response_monthly_rows,
              ["year","month","coast_region","water_class","n_pixels",
               "mean_SPEI3","mean_WUE_T_pct_change","mean_abs_WUE_T_pct_change",
               "pct_pixels_increase","pct_pixels_decrease",
               "pct_Stable","pct_Mild","pct_Moderate","pct_Large","pct_Strong"],
              os.path.join(OUTPUT_DIR, "EDI_response_monthly_coast_summary.csv"))

    write_csv(response_annual_rows,
              ["year","coast_region","water_class",
               "mean_WUE_T_pct_change","mean_abs_WUE_T_pct_change",
               "pct_pixels_increase","pct_pixels_decrease",
               "pct_Stable","pct_Mild","pct_Moderate","pct_Large","pct_Strong"],
              os.path.join(OUTPUT_DIR, "EDI_response_annual_coast_summary.csv"))

    write_csv(response_v2_monthly_rows,
              ["year","month","coast_region","n_pixels",
               "mean_SPEI3","upland_mean_WUE_T_pct_change","upland_mean_abs_WUE_T_pct_change",
               "upland_mean_lower_pct_change","upland_mean_upper_pct_change",
               "upland_mean_prediction_interval_width",
               "pct_pixels_decrease","pct_pixels_increase",
               "pct_pixels_strong_decrease","pct_pixels_strong_increase",
               "pct_pixels_robust_strong_response","pct_pixels_possible_strong_response",
               "dry_n_pixels","dry_pct_pixels",
               "dry_mean_WUE_T_pct_change","dry_mean_abs_WUE_T_pct_change",
               "dry_pct_pixels_decrease","dry_pct_pixels_increase",
               "dry_pct_pixels_strong_decrease","dry_pct_pixels_strong_increase"],
              os.path.join(OUTPUT_DIR, "EDI_response_v2_upland_monthly_coast_summary.csv"))

    write_csv(response_v2_annual_rows,
              ["year","coast_region",
               "mean_SPEI3","upland_mean_WUE_T_pct_change","upland_mean_abs_WUE_T_pct_change",
               "upland_mean_lower_pct_change","upland_mean_upper_pct_change",
               "upland_mean_prediction_interval_width",
               "pct_pixels_decrease","pct_pixels_increase",
               "pct_pixels_strong_decrease","pct_pixels_strong_increase",
               "pct_pixels_robust_strong_response","pct_pixels_possible_strong_response",
               "dry_pct_pixels",
               "dry_mean_WUE_T_pct_change","dry_mean_abs_WUE_T_pct_change",
               "dry_pct_pixels_decrease","dry_pct_pixels_increase",
               "dry_pct_pixels_strong_decrease","dry_pct_pixels_strong_increase"],
              os.path.join(OUTPUT_DIR, "EDI_response_v2_upland_annual_coast_summary.csv"))

    write_csv(response_v2_validation_rows,
              ["coast_region",
               "spatial_period_mean_WUE_T_pct_change",
               "spatial_period_mean_abs_WUE_T_pct_change",
               "spatial_pct_pixels_decrease",
               "spatial_pct_pixels_increase",
               "spatial_pct_pixels_strong_decrease",
               "spatial_pct_pixels_strong_increase",
               "spatial_dry_pct_pixels",
               "spatial_dry_mean_WUE_T_pct_change",
               "spatial_dry_mean_abs_WUE_T_pct_change",
               "spatial_dry_pct_pixels_strong_decrease",
               "spatial_dry_pct_pixels_strong_increase",
               "Q3_n_sites",
               "Q3_SPEI3_mean_slope",
               "Q3_SPEI3_median_slope",
               "Q3_pct_negative_site_slopes",
               "Q3_SPEI3_threshold_dry_increase_5pct",
               "Q3_SPEI3_threshold_dry_decrease_5pct",
               "Q3_SPEI3_threshold_wet_increase_5pct",
               "Q3_SPEI3_threshold_wet_decrease_5pct",
               "Q4_WUE_T_SPEI3_upland_dry_SPEI_prob25",
               "Q4_WUE_T_SPEI3_upland_dry_SPEI_prob50",
               "Q4_WUE_T_SPEI3_upland_dry_SPEI_prob75",
               "Q4_WUE_T_SPEI3_upland_min_decline_probability_5pct",
               "Q4_WUE_T_SPEI3_upland_max_decline_probability_5pct"],
              os.path.join(OUTPUT_DIR, "EDI_response_v2_regional_validation_Q3_Q4.csv"))

    write_csv(response_v3_rows,
              ["coast_region",
               "dry_month_frequency_pct",
               "dry_mean_WUE_T_pct_change",
               "dry_mean_abs_WUE_T_pct_change",
               "drought_adjustment_index_pct",
               "drought_impairment_index_pct",
               "net_adjustment_minus_impairment_pct",
               "Q4_WUE_T_SPEI3_upland_max_decline_probability_5pct",
               "Q3_SPEI3_mean_slope",
               "Q3_pct_negative_site_slopes",
               "interpretation"],
              os.path.join(OUTPUT_DIR, "EDI_response_v3_drought_adjustment_impairment_summary.csv"))

    write_csv(communication_annual_rows,
              ["year","coast_region","event_threshold_pct",
               "annual_mean_negative_impacted_area_pct",
               "annual_mean_positive_impacted_area_pct",
               "annual_mean_any_impacted_area_pct",
               "annual_mean_net_positive_minus_negative_pct",
               "annual_mean_dry_area_pct",
               "dry_negative_impacted_area_pct",
               "dry_positive_impacted_area_pct"],
              os.path.join(
                  OUTPUT_DIR,
                  f"EDI_response_communication_annual_impacted_area_by_coast_"
                  f"{SPATIAL_START_YEAR}_{SPATIAL_END_YEAR}.csv"
              ))

    # ---- Write rasters -------------------------------------------------------
    def write_raster(data, varname, long_name, outpath):
        ds = nc.Dataset(outpath, "w", format="NETCDF4")
        ds.createDimension("lat", nlat); ds.createDimension("lon", nlon)
        la = ds.createVariable("lat","f4",("lat",)); la.units="degrees_north"; la[:]=lat1d
        lo = ds.createVariable("lon","f4",("lon",)); lo.units="degrees_east";  lo[:]=lon1d
        v  = ds.createVariable(varname,"f4",("lat","lon"),fill_value=np.float32(FILL),
                                zlib=True,complevel=4)
        v.long_name = long_name
        d = data.copy().astype(np.float32)
        d[(~valid_mask) | (~np.isfinite(d))] = np.float32(FILL)
        v[:] = d
        ds.created = datetime.now().isoformat(); ds.close()

    write_raster(mean_edi,    "EDI_mean",        "Time-mean EDI",
                 os.path.join(OUTPUT_DIR, "EDI_logistic_pixel_mean.nc"))
    write_raster(freq_impact, "EDI_freq_impact", "% months in Impact class (EDI > 0.70)",
                 os.path.join(OUTPUT_DIR, "EDI_logistic_pixel_freq_impact.nc"))

    for water_class in water_classes:
        suffix = water_class.lower()
        write_raster(
            mean_pct_response[water_class],
            "WUE_T_pct_change_mean",
            f"Time-mean predicted WUE_T percent change, {water_class} response scenario",
            os.path.join(OUTPUT_DIR, f"EDI_response_pixel_mean_pct_change_{suffix}.nc"),
        )
        write_raster(
            mean_abs_response[water_class],
            "WUE_T_abs_pct_change_mean",
            f"Time-mean absolute predicted WUE_T percent change, {water_class} response scenario",
            os.path.join(OUTPUT_DIR, f"EDI_response_pixel_mean_abs_change_{suffix}.nc"),
        )
        write_raster(
            freq_strong_response[water_class],
            "WUE_T_freq_strong_change",
            f"% months with >=5% absolute predicted WUE_T change, {water_class} response scenario",
            os.path.join(OUTPUT_DIR, f"EDI_response_pixel_freq_strong_{suffix}.nc"),
        )

    response_v2_raster_specs = [
        ("mean_pct_change", "upland_mean_pct_change",
         "Upland response_v2: time-mean predicted WUE_T percent change"),
        ("mean_abs_change", "upland_mean_abs_change",
         "Upland response_v2: time-mean absolute predicted WUE_T percent change"),
        ("freq_decrease", "upland_freq_decrease",
         "Upland response_v2: % months with predicted WUE_T decrease"),
        ("freq_increase", "upland_freq_increase",
         "Upland response_v2: % months with predicted WUE_T increase"),
        ("freq_strong_decrease", "upland_freq_strong_decrease",
         "Upland response_v2: % months with <= -5% predicted WUE_T change"),
        ("freq_strong_increase", "upland_freq_strong_increase",
         "Upland response_v2: % months with >= 5% predicted WUE_T change"),
        ("dry_month_frequency", "upland_dry_month_frequency",
         "Upland response_v2: % months with SPEI-3 <= -1"),
        ("dry_mean_pct_change", "upland_dry_mean_pct_change",
         "Upland response_v2: dry-month mean predicted WUE_T percent change"),
        ("dry_mean_abs_change", "upland_dry_mean_abs_change",
         "Upland response_v2: dry-month mean absolute predicted WUE_T percent change"),
        ("dry_freq_decrease", "upland_dry_freq_decrease",
         "Upland response_v2: dry-month % months with predicted WUE_T decrease"),
        ("dry_freq_increase", "upland_dry_freq_increase",
         "Upland response_v2: dry-month % months with predicted WUE_T increase"),
        ("dry_freq_strong_decrease", "upland_dry_freq_strong_decrease",
         "Upland response_v2: dry-month % months with <= -5% predicted WUE_T change"),
        ("dry_freq_strong_increase", "upland_dry_freq_strong_increase",
         "Upland response_v2: dry-month % months with >= 5% predicted WUE_T change"),
        ("dry_freq_possible_decrease", "upland_dry_freq_possible_decrease",
         "Upland response_v2: dry-month % months where interval allows <= -5% response"),
        ("dry_freq_robust_decrease", "upland_dry_freq_robust_decrease",
         "Upland response_v2: dry-month % months where interval upper bound is <= -5%"),
        ("dry_freq_possible_increase", "upland_dry_freq_possible_increase",
         "Upland response_v2: dry-month % months where interval allows >= 5% response"),
        ("dry_freq_robust_increase", "upland_dry_freq_robust_increase",
         "Upland response_v2: dry-month % months where interval lower bound is >= 5%"),
        ("mean_lower_pct_change", "upland_mean_lower_pct_change",
         "Upland response_v2: mean lower prediction interval WUE_T percent change"),
        ("mean_upper_pct_change", "upland_mean_upper_pct_change",
         "Upland response_v2: mean upper prediction interval WUE_T percent change"),
        ("mean_ci_width", "upland_mean_prediction_interval_width",
         "Upland response_v2: mean prediction interval width"),
        ("freq_possible_strong", "upland_freq_possible_strong",
         "Upland response_v2: % months where interval allows >=5% absolute response"),
        ("freq_robust_strong", "upland_freq_robust_strong",
         "Upland response_v2: % months where interval supports >=5% absolute response with same sign"),
        ("freq_possible_decrease", "upland_freq_possible_decrease",
         "Upland response_v2: % months where interval allows <= -5% response"),
        ("freq_robust_decrease", "upland_freq_robust_decrease",
         "Upland response_v2: % months where interval upper bound is <= -5%"),
        ("freq_possible_increase", "upland_freq_possible_increase",
         "Upland response_v2: % months where interval allows >= 5% response"),
        ("freq_robust_increase", "upland_freq_robust_increase",
         "Upland response_v2: % months where interval lower bound is >= 5%"),
        ("drought_adjustment_index", "upland_drought_adjustment_index",
         "Upland drought adjustment index: dry-month % months with >= 5% WUE_T increase"),
        ("drought_impairment_index", "upland_drought_impairment_index",
         "Upland drought impairment index: dry-month % months with <= -5% WUE_T decrease"),
        ("drought_net_adjustment_minus_impairment", "upland_drought_net_adjustment_minus_impairment",
         "Upland drought adjustment minus impairment index"),
        ("drought_possible_adjustment_index", "upland_drought_possible_adjustment_index",
         "Upland drought adjustment index allowing prediction interval uncertainty"),
        ("drought_robust_adjustment_index", "upland_drought_robust_adjustment_index",
         "Upland drought adjustment index supported by lower prediction interval >= 5%"),
        ("drought_possible_impairment_index", "upland_drought_possible_impairment_index",
         "Upland drought impairment index allowing prediction interval uncertainty"),
        ("drought_robust_impairment_index", "upland_drought_robust_impairment_index",
         "Upland drought impairment index supported by upper prediction interval <= -5%"),
    ]
    for key, varname, long_name in response_v2_raster_specs:
        write_raster(
            response_v2_rasters[key],
            varname,
            long_name,
            os.path.join(OUTPUT_DIR, f"EDI_response_v2_{key}_upland.nc"),
        )

    print(f"\nCSVs and rasters written to {OUTPUT_DIR}")

    # =========================================================================
    # FIGURES
    # =========================================================================

    # helper: date float from year+month
    def ym2date(y, m):
        return y + (m - 0.5) / 12

    dates = [ym2date(r[0], r[1]) for r in monthly_rows
             if r[2] == region_names[0]]
    all_years = sorted(set(r[0] for r in monthly_rows))

    # ---- Response Plot 1: mean absolute WUE_T change maps --------------------
    fig_r1, axes_r1 = plt.subplots(
        len(water_classes), 2,
        figsize=(14, 4.2 * len(water_classes)),
        squeeze=False,
    )
    for row, water_class in enumerate(water_classes):
        abs_plot = np.where(valid_mask, mean_abs_response[water_class], np.nan)
        strong_plot = np.where(valid_mask, freq_strong_response[water_class], np.nan)

        ax_abs = axes_r1[row, 0]
        ax_abs.set_facecolor("#d0e8f5")
        im_abs = ax_abs.pcolormesh(
            lon1d, lat1d, abs_plot, cmap="YlOrRd", vmin=0, vmax=10, shading="auto"
        )
        plt.colorbar(im_abs, ax=ax_abs, label="Mean |WUE_T change| (%)", fraction=0.035, pad=0.02)
        ax_abs.set_xlim(-180, -60); ax_abs.set_ylim(17, 77)
        ax_abs.set_xlabel("Longitude"); ax_abs.set_ylabel("Latitude")
        ax_abs.set_title(f"{water_class}: Mean Absolute WUE_T Change", fontweight="bold")

        ax_freq = axes_r1[row, 1]
        ax_freq.set_facecolor("#d0e8f5")
        im_freq = ax_freq.pcolormesh(
            lon1d, lat1d, strong_plot, cmap="PuRd", vmin=0, vmax=100, shading="auto"
        )
        plt.colorbar(im_freq, ax=ax_freq, label="% months >=5% change", fraction=0.035, pad=0.02)
        ax_freq.set_xlim(-180, -60); ax_freq.set_ylim(17, 77)
        ax_freq.set_xlabel("Longitude"); ax_freq.set_ylabel("Latitude")
        ax_freq.set_title(f"{water_class}: Strong Response Frequency", fontweight="bold")

    fig_r1.suptitle(
        "Response-Based EDI Spatial Upscaling: Ecosystem Response Scenarios",
        fontweight="bold",
        y=1.01,
    )
    fig_r1.tight_layout()
    fig_r1.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_plot_01_spatial_response_scenarios.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_r1)

    # ---- Response Plot 2: monthly response by coast and ecosystem ------------
    fig_r2, axes_r2 = plt.subplots(
        len(region_names), 1,
        figsize=(14, 3.2 * len(region_names)),
        sharex=True,
    )
    if len(region_names) == 1:
        axes_r2 = [axes_r2]
    for ax, name in zip(axes_r2, region_names):
        for water_class in water_classes:
            sub = [r for r in response_monthly_rows if r[2] == name and r[3] == water_class]
            if not sub:
                continue
            ax.plot(
                [ym2date(r[0], r[1]) for r in sub],
                [r[7] for r in sub],
                color=WATER_COLORS.get(water_class, "black"),
                linewidth=1.0,
                alpha=0.9,
                label=water_class,
            )
        ax.axhline(0.5, color="#999999", linewidth=0.8, linestyle=":")
        ax.axhline(2.5, color="#666666", linewidth=0.8, linestyle="--")
        ax.axhline(5.0, color="#333333", linewidth=0.8, linestyle="--")
        ax.set_ylabel("Mean |WUE_T change| (%)")
        ax.set_title(display_coast(name), fontweight="bold")
    axes_r2[-1].set_xlabel("Year")
    axes_r2[0].legend(loc="upper left", ncol=len(water_classes), fontsize=9)
    fig_r2.suptitle(
        "Response-Based EDI: Mean Predicted WUE_T Change by Coast",
        fontweight="bold",
        y=1.01,
    )
    fig_r2.tight_layout()
    fig_r2.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_plot_02_timeseries_by_coast.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_r2)

    # ---- Response v3 Plot: model-output time series by coast ----------------
    fig_ts, axes_ts = plt.subplots(
        len(region_names), 1,
        figsize=(15, 3.6 * len(region_names)),
        sharex=True,
    )
    if len(region_names) == 1:
        axes_ts = [axes_ts]
    for ax, name in zip(axes_ts, region_names):
        sub = [r for r in response_v2_monthly_rows if r[2] == name]
        if not sub:
            continue
        xs = [ym2date(r[0], r[1]) for r in sub]
        mean_change = [r[5] for r in sub]
        strong_decrease = [r[12] for r in sub]
        strong_increase = [r[13] for r in sub]
        dry_pct = [r[17] for r in sub]
        dry_impairment = [0 if not np.isfinite(r[22]) else r[22] for r in sub]
        dry_adjustment = [0 if not np.isfinite(r[23]) else r[23] for r in sub]

        ax.axhline(0, color="#555555", linewidth=0.8)
        ax.fill_between(xs, 0, dry_pct, color="#D9A441", alpha=0.18, label="Dry pixel frequency")
        ax.plot(xs, mean_change, color="#222222", linewidth=1.1, label="Mean WUE_T change")
        ax.plot(xs, strong_increase, color="#B83280", linewidth=0.9, alpha=0.8,
                label="Strong increase")
        ax.plot(xs, strong_decrease, color="#2B6CB0", linewidth=0.9, alpha=0.8,
                label="Strong decrease")
        ax.plot(xs, dry_adjustment, color="#C2185B", linewidth=1.2, linestyle="--",
                label="Dry adjustment")
        ax.plot(xs, dry_impairment, color="#08519C", linewidth=1.2, linestyle="--",
                label="Dry impairment")
        ax.set_ylim(-20, 105)
        ax.set_ylabel("% or WUE_T %")
        ax.set_title(display_coast(name), fontweight="bold")
    axes_ts[-1].set_xlabel("Year")
    handles, labels = axes_ts[0].get_legend_handles_labels()
    fig_ts.legend(handles, labels, loc="lower center", ncol=3, fontsize=9,
                  bbox_to_anchor=(0.5, -0.01))
    fig_ts.suptitle(
        "Response v3 Model Outputs by Coast Region: Upland Adjustment vs Impairment",
        fontweight="bold",
        y=1.01,
    )
    fig_ts.tight_layout()
    fig_ts.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_v3_plot_02_model_output_timeseries_by_coast.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_ts)

    fig_ann, axes_ann = plt.subplots(
        2, 2,
        figsize=(14, 9),
        sharex=True,
        sharey=True,
    )
    for ax, name in zip(axes_ann.flat, region_names):
        sub = [r for r in response_v2_annual_rows if r[1] == name]
        if not sub:
            continue
        xs = [r[0] for r in sub]
        dry_pct = [r[14] for r in sub]
        dry_mean = [r[15] for r in sub]
        dry_impairment = [0 if not np.isfinite(r[19]) else r[19] for r in sub]
        dry_adjustment = [0 if not np.isfinite(r[20]) else r[20] for r in sub]
        ax.axhline(0, color="#555555", linewidth=0.8)
        ax.fill_between(xs, 0, dry_pct, color="#D9A441", alpha=0.22, label="Dry pixel frequency")
        ax.plot(xs, dry_mean, color="#222222", marker="o", markersize=3,
                linewidth=1.0, label="Dry mean WUE_T change")
        ax.plot(xs, dry_adjustment, color="#C2185B", marker="o", markersize=3,
                linewidth=1.0, label="Dry adjustment")
        ax.plot(xs, dry_impairment, color="#08519C", marker="o", markersize=3,
                linewidth=1.0, label="Dry impairment")
        ax.set_title(display_coast(name), fontweight="bold")
        ax.set_ylim(-20, 105)
        ax.set_ylabel("% or WUE_T %")
    for ax in axes_ann[-1, :]:
        ax.set_xlabel("Year")
    handles, labels = axes_ann.flat[0].get_legend_handles_labels()
    fig_ann.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
                   bbox_to_anchor=(0.5, -0.02))
    fig_ann.suptitle(
        "Annual Upland Drought Model Outputs by Coast Region",
        fontweight="bold",
        y=1.01,
    )
    fig_ann.tight_layout()
    fig_ann.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_v3_plot_03_annual_model_outputs_by_coast.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_ann)

    # ---- Response v2 Plot: Upland direction, dry-only, and uncertainty -------
    fig_v2, axes_v2 = plt.subplots(3, 2, figsize=(14, 13))
    v2_panels = [
        ("mean_abs_change", "Mean Absolute WUE_T Change", "YlOrRd", 0, 15, "%"),
        ("freq_strong_decrease", "Strong Decrease Frequency", "Blues", 0, 100, "% months"),
        ("freq_strong_increase", "Strong Increase Frequency", "PuRd", 0, 100, "% months"),
        ("dry_freq_strong_decrease", "Dry-Only Strong Decrease", "Blues", 0, 100, "% dry months"),
        ("dry_freq_strong_increase", "Dry-Only Strong Increase", "PuRd", 0, 100, "% dry months"),
        ("mean_ci_width", "Prediction Interval Width", "Greys", 0, 35, "%"),
    ]
    for ax, (key, title, cmap, vmin, vmax, label) in zip(axes_v2.flat, v2_panels):
        ax.set_facecolor("#d0e8f5")
        data = np.where(valid_mask, response_v2_rasters[key], np.nan)
        im = ax.pcolormesh(lon1d, lat1d, data, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        plt.colorbar(im, ax=ax, label=label, fraction=0.035, pad=0.02)
        ax.set_xlim(-180, -60); ax.set_ylim(17, 77)
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
        ax.set_title(title, fontweight="bold")
    fig_v2.suptitle(
        "Response v2: Upland WUE_T Spatial Test Against SPEI-3",
        fontweight="bold",
        y=1.01,
    )
    fig_v2.tight_layout()
    fig_v2.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_v2_plot_01_upland_direction_dry_uncertainty.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_v2)

    # ---- Communication Plot: 2000-2025 positive/negative event frequencies --
    fig_comm_map, axes_comm_map = plt.subplots(1, 2, figsize=(15, 6.4), sharex=True, sharey=True)
    comm_map_specs = [
        (
            "freq_strong_decrease",
            f"Negative events: WUE_T <= -{EVENT_THRESHOLD_PCT:g}%",
            "Blues",
            f"% months, {SPATIAL_START_YEAR}-{SPATIAL_END_YEAR}",
        ),
        (
            "freq_strong_increase",
            f"Positive events: WUE_T >= +{EVENT_THRESHOLD_PCT:g}%",
            "Reds",
            f"% months, {SPATIAL_START_YEAR}-{SPATIAL_END_YEAR}",
        ),
    ]
    for ax, (key, title, cmap, cbar_label) in zip(axes_comm_map, comm_map_specs):
        ax.set_facecolor("#d0e8f5")
        data = np.where(valid_mask, response_v2_rasters[key], np.nan)
        im = ax.pcolormesh(lon1d, lat1d, data, cmap=cmap, vmin=0, vmax=100, shading="auto")
        plt.colorbar(im, ax=ax, label=cbar_label, fraction=0.035, pad=0.02)
        ax.set_xlim(-180, -60)
        ax.set_ylim(17, 77)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title, fontweight="bold")
    fig_comm_map.suptitle(
        f"Upland WUE_T Event Frequency from SPEI-3 Spatial Upscaling ({SPATIAL_START_YEAR}-{SPATIAL_END_YEAR})",
        fontweight="bold",
        y=1.02,
    )
    fig_comm_map.tight_layout()
    fig_comm_map.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"EDI_response_communication_plot_01_event_frequency_maps_"
            f"{SPATIAL_START_YEAR}_{SPATIAL_END_YEAR}.png",
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_comm_map)

    fig_comm_ts, axes_comm_ts = plt.subplots(
        len(region_names),
        1,
        figsize=(13, 3.1 * len(region_names)),
        sharex=True,
        sharey=True,
    )
    if len(region_names) == 1:
        axes_comm_ts = [axes_comm_ts]
    for ax, name in zip(axes_comm_ts, region_names):
        sub = [r for r in communication_annual_rows if r[1] == name]
        if not sub:
            continue
        xs = [r[0] for r in sub]
        neg = [r[3] for r in sub]
        pos = [r[4] for r in sub]
        any_impacted = [r[5] for r in sub]
        dry_area = [r[7] for r in sub]
        ax.fill_between(xs, 0, dry_area, color="#D9A441", alpha=0.16, label="Dry area")
        ax.plot(xs, neg, color="#08519C", linewidth=1.8, marker="o", markersize=3.2,
                label=f"Negative event area (<= -{EVENT_THRESHOLD_PCT:g}%)")
        ax.plot(xs, pos, color="#B2182B", linewidth=1.8, marker="o", markersize=3.2,
                label=f"Positive event area (>= +{EVENT_THRESHOLD_PCT:g}%)")
        ax.plot(xs, any_impacted, color="#222222", linewidth=1.1, linestyle="--",
                label="Any signed event area")
        ax.set_title(display_coast(name), fontweight="bold")
        ax.set_ylabel("Mean annual\narea (%)")
        ax.set_ylim(0, 100)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100))
    axes_comm_ts[-1].set_xlabel("Year")
    handles, labels = axes_comm_ts[0].get_legend_handles_labels()
    fig_comm_ts.legend(handles, labels, loc="lower center", ncol=2, fontsize=9,
                       bbox_to_anchor=(0.5, -0.015))
    fig_comm_ts.suptitle(
        f"Mean Annual Impacted Area by Coastline ({SPATIAL_START_YEAR}-{SPATIAL_END_YEAR})",
        fontweight="bold",
        y=1.01,
    )
    fig_comm_ts.tight_layout()
    fig_comm_ts.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"EDI_response_communication_plot_02_annual_impacted_area_by_coast_"
            f"{SPATIAL_START_YEAR}_{SPATIAL_END_YEAR}.png",
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_comm_ts)

    fig_val, ax_val = plt.subplots(figsize=(12, 6))
    x = np.arange(len(region_names))
    validation_by_region = {r[0]: r for r in response_v2_validation_rows}
    strong_dec = [validation_by_region[name][10] for name in region_names]
    strong_inc = [validation_by_region[name][11] for name in region_names]
    q3_slope = [validation_by_region[name][13] for name in region_names]
    width = 0.36
    ax_val.bar(x - width/2, strong_dec, width, color="#2B6CB0", label="Dry strong decrease")
    ax_val.bar(x + width/2, strong_inc, width, color="#B83280", label="Dry strong increase")
    ax_val.set_xticks(x)
    ax_val.set_xticklabels([display_coast(name) for name in region_names], rotation=20, ha="right")
    ax_val.set_ylabel("% dry pixels/months")
    ax_val.set_title("Response v2 Validation Summary: Spatial Dry Signal vs Q3 Site Slopes",
                     fontweight="bold")
    ax_val.legend(loc="upper left")
    ax_slope = ax_val.twinx()
    ax_slope.plot(x, q3_slope, color="#222222", marker="o", linewidth=1.5, label="Q3 SPEI-3 mean slope")
    ax_slope.axhline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax_slope.set_ylabel("Q3 site-level mean slope")
    ax_slope.legend(loc="upper right")
    fig_val.tight_layout()
    fig_val.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_v2_plot_02_validation_Q3_Q4.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_val)

    # ---- Response v3 Plot: adjustment versus impairment framing --------------
    fig_v3, axes_v3 = plt.subplots(2, 3, figsize=(16, 9))
    v3_panels = [
        ("drought_adjustment_index", "Drought Adjustment Index\nDry SPEI + WUE_T Increase", "PuRd", 0, 100, "% dry months"),
        ("drought_impairment_index", "Drought Impairment Index\nDry SPEI + WUE_T Decrease", "Blues", 0, 100, "% dry months"),
        ("drought_net_adjustment_minus_impairment", "Net Adjustment - Impairment", "RdBu_r", -100, 100, "percentage points"),
        ("drought_robust_adjustment_index", "Robust Adjustment\nLower PI >= +5%", "PuRd", 0, 100, "% dry months"),
        ("drought_robust_impairment_index", "Robust Impairment\nUpper PI <= -5%", "Blues", 0, 100, "% dry months"),
        ("dry_month_frequency", "Dry Month Frequency\nSPEI-3 <= -1", "YlOrBr", 0, 50, "% months"),
    ]
    for ax, (key, title, cmap, vmin, vmax, label) in zip(axes_v3.flat, v3_panels):
        ax.set_facecolor("#d0e8f5")
        data = np.where(valid_mask, response_v2_rasters[key], np.nan)
        im = ax.pcolormesh(lon1d, lat1d, data, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        plt.colorbar(im, ax=ax, label=label, fraction=0.04, pad=0.02)
        ax.set_xlim(-180, -60); ax.set_ylim(17, 77)
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
        ax.set_title(title, fontweight="bold")
    fig_v3.suptitle(
        "Response v3: Separate Drought Adjustment from Drought Impairment",
        fontweight="bold",
        y=1.02,
    )
    fig_v3.tight_layout()
    fig_v3.savefig(
        os.path.join(OUTPUT_DIR, "EDI_response_v3_plot_01_adjustment_vs_impairment.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig_v3)

    # ---- Plot 1: time-mean EDI map ------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    ax1.set_facecolor("#d0e8f5")
    edi_plot = np.where(valid_mask, mean_edi, np.nan)
    im = ax1.pcolormesh(lon1d, lat1d, edi_plot,
                        cmap="RdBu_r", vmin=0.3, vmax=0.8, shading="auto")
    plt.colorbar(im, ax=ax1, label="Time-mean EDI", fraction=0.03, pad=0.02)

    # Coast region boundaries (approximate box outlines)
    region_boxes = {
        "AK Coast":       ([-180,-60],[50,76]),
        "Pacific Coast":  ([-130,-116],[31,50]),
        "Gulf Coast":     ([-100,-78],[24,32]),
        "Atlantic Coast": ([-82,-65],[24,50]),
    }
    for name, color in COAST_COLORS.items():
        lons, lats = region_boxes[name]
        rect = mpatches.Rectangle((lons[0], lats[0]),
                                   lons[1]-lons[0], lats[1]-lats[0],
                                   linewidth=1.2, edgecolor=color,
                                   facecolor="none", linestyle="--", label=display_coast(name))
        ax1.add_patch(rect)

    ax1.set_xlim(-180, -60); ax1.set_ylim(17, 77)
    ax1.set_xlabel("Longitude"); ax1.set_ylabel("Latitude")
    ax1.set_title("Time-Mean Ecological Drought Index (EDI) — SPEI-3 Calibrated",
                  fontweight="bold")
    ax1.legend(loc="lower right", fontsize=9, title="Coast region")
    ax1.axhline(50, color="gray", linewidth=0.5, linestyle=":")
    ax1.axvline(-120, color="gray", linewidth=0.5, linestyle=":")
    ax1.axvline(-100, color="gray", linewidth=0.5, linestyle=":")
    for th, label in [(EDI_WATCH,"Watch"), (EDI_STRESS,"Stress"), (EDI_IMPACT,"Impact")]:
        pass  # thresholds on colorbar already
    fig1.tight_layout()
    fig1.savefig(os.path.join(OUTPUT_DIR,"EDI_plot_01_mean_map.png"),
                 dpi=300, bbox_inches="tight")
    plt.close(fig1)

    # ---- Plot 2: % months in Impact class -----------------------------------
    fig2, ax2 = plt.subplots(figsize=(12, 6))
    ax2.set_facecolor("#d0e8f5")
    fi_plot = np.where(valid_mask, freq_impact, np.nan)
    im2 = ax2.pcolormesh(lon1d, lat1d, fi_plot,
                         cmap="OrRd", vmin=0, vmax=60, shading="auto")
    plt.colorbar(im2, ax=ax2, label="% months EDI > 0.70 (Impact)", fraction=0.03, pad=0.02)
    ax2.set_xlim(-180,-60); ax2.set_ylim(17,77)
    ax2.set_xlabel("Longitude"); ax2.set_ylabel("Latitude")
    ax2.set_title("Frequency of EDI 'Impact' Class (EDI > 0.70) — SPEI-3 Calibrated",
                  fontweight="bold")
    fig2.tight_layout()
    fig2.savefig(os.path.join(OUTPUT_DIR,"EDI_plot_02_impact_freq_map.png"),
                 dpi=300, bbox_inches="tight")
    plt.close(fig2)

    # ---- Plot 3: monthly mean EDI time series by coast region ----------------
    fig3, ax3 = plt.subplots(figsize=(14, 6))
    for name, color in COAST_COLORS.items():
        sub   = [r for r in monthly_rows if r[2] == name]
        if not sub: continue
        xvals = [ym2date(r[0], r[1]) for r in sub]
        yvals = [r[5] for r in sub]          # mean_EDI
        ax3.plot(xvals, yvals, color=color, linewidth=1.0, label=display_coast(name), alpha=0.85)

    for th, col, lbl in [(EDI_WATCH,"#92C5DE","Watch 0.30"),
                          (EDI_STRESS,"#F4A582","Stress 0.50"),
                          (EDI_IMPACT,"#D6604D","Impact 0.70")]:
        ax3.axhline(th, color=col, linewidth=0.8, linestyle="--")
        ax3.text(all_years[-1] + 0.15, th, lbl, va="center", fontsize=8, color=col)

    ax3.set_xlabel("Year")
    ax3.set_ylabel("Mean EDI")
    ax3.set_title("Monthly Mean EDI by Coast Region", fontweight="bold")
    ax3.legend(loc="upper left", fontsize=9)
    ax3.set_ylim(0, 1)
    fig3.tight_layout()
    fig3.savefig(os.path.join(OUTPUT_DIR,"EDI_plot_03_timeseries_by_coast.png"),
                 dpi=300, bbox_inches="tight")
    plt.close(fig3)

    # ---- Plot 4: stacked severity proportions over time (all regions) -------
    fig4, axes4 = plt.subplots(len(region_names), 1,
                               figsize=(14, 3.5*len(region_names)), sharex=True)
    if len(region_names) == 1:
        axes4 = [axes4]
    sev_cols  = ["pct_None","pct_Watch","pct_Stress","pct_Impact"]
    sev_names = ["None","Watch","Stress","Impact"]
    sev_cidx  = [5+1, 5+2, 5+3, 5+4]   # column indices in monthly_rows: 6,7,8,9

    for ax, name in zip(axes4, region_names):
        sub = [r for r in monthly_rows if r[2] == name]
        if not sub:
            ax.set_title(display_coast(name)); continue
        xs   = [ym2date(r[0], r[1]) for r in sub]
        base = np.zeros(len(sub))
        for sname, cidx in zip(sev_names, [6,7,8,9]):
            vals = np.array([r[cidx] for r in sub])
            ax.fill_between(xs, base, base + vals,
                            color=SEVERITY_COLORS[sname],
                            alpha=0.85, label=sname)
            base += vals
        ax.set_ylim(0,100); ax.set_ylabel("% pixels")
        ax.set_title(display_coast(name), fontweight="bold")
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100))

    axes4[-1].set_xlabel("Year")
    handles = [mpatches.Patch(color=SEVERITY_COLORS[s], label=s) for s in sev_names]
    axes4[0].legend(handles=handles, loc="upper right", ncol=4, fontsize=9)
    fig4.suptitle("EDI Severity Class Proportions by Coast Region",
                  fontweight="bold", y=1.01)
    fig4.tight_layout()
    fig4.savefig(os.path.join(OUTPUT_DIR,"EDI_plot_04_severity_stack.png"),
                 dpi=300, bbox_inches="tight")
    plt.close(fig4)

    # ---- Plot 5: annual mean EDI by coast region (dot + line) ---------------
    fig5, axes5 = plt.subplots(1, 2, figsize=(14, 6))

    # Left: annual mean EDI
    ax5a = axes5[0]
    for name, color in COAST_COLORS.items():
        sub = [r for r in annual_rows if r[1] == name]
        if not sub: continue
        ax5a.plot([r[0] for r in sub], [r[2] for r in sub],
                  marker="o", markersize=4, color=color, label=display_coast(name), linewidth=1.2)
    for th, col in [(EDI_WATCH,"#92C5DE"),(EDI_STRESS,"#F4A582"),(EDI_IMPACT,"#D6604D")]:
        ax5a.axhline(th, color=col, linewidth=0.8, linestyle="--")
    ax5a.set_xlabel("Year"); ax5a.set_ylabel("Annual mean EDI")
    ax5a.set_title("Annual Mean EDI by Coast Region", fontweight="bold")
    ax5a.legend(fontsize=9); ax5a.set_ylim(0,1)

    # Right: annual % Impact
    ax5b = axes5[1]
    for name, color in COAST_COLORS.items():
        sub = [r for r in annual_rows if r[1] == name]
        if not sub: continue
        ax5b.plot([r[0] for r in sub], [r[7] for r in sub],
                  marker="o", markersize=4, color=color, label=display_coast(name), linewidth=1.2)
    ax5b.set_xlabel("Year"); ax5b.set_ylabel("% pixels in Impact class")
    ax5b.set_title("Annual % EDI 'Impact' Pixels by Coast Region", fontweight="bold")
    ax5b.legend(fontsize=9); ax5b.set_ylim(0,100)
    ax5b.yaxis.set_major_formatter(PercentFormatter(xmax=100))

    fig5.tight_layout()
    fig5.savefig(os.path.join(OUTPUT_DIR,"EDI_plot_05_annual_summary.png"),
                 dpi=300, bbox_inches="tight")
    plt.close(fig5)

    # ---- Plot 6: composite (maps + time series) -----------------------------
    fig6 = plt.figure(figsize=(16, 20))
    gs   = fig6.add_gridspec(4, 2, hspace=0.45, wspace=0.3)

    # Row 0: mean EDI map + impact freq map
    ax_m1 = fig6.add_subplot(gs[0, 0])
    ax_m2 = fig6.add_subplot(gs[0, 1])
    for ax, data, cmap, vmax, title in [
        (ax_m1, edi_plot,  "RdBu_r", 0.8,  "Time-Mean EDI"),
        (ax_m2, fi_plot,   "OrRd",   60,   "% Months in Impact class"),
    ]:
        ax.set_facecolor("#d0e8f5")
        im = ax.pcolormesh(lon1d, lat1d, data, cmap=cmap,
                           vmin=0.3 if cmap=="RdBu_r" else 0,
                           vmax=vmax, shading="auto")
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        ax.set_xlim(-180,-60); ax.set_ylim(17,77)
        ax.set_title(title, fontweight="bold", fontsize=11)
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")

    # Row 1: monthly time series
    ax_ts = fig6.add_subplot(gs[1, :])
    for name, color in COAST_COLORS.items():
        sub = [r for r in monthly_rows if r[2] == name]
        if not sub: continue
        ax_ts.plot([ym2date(r[0],r[1]) for r in sub],
                   [r[5] for r in sub],
                   color=color, linewidth=1.0, label=display_coast(name), alpha=0.85)
    for th, col in [(EDI_WATCH,"#92C5DE"),(EDI_STRESS,"#F4A582"),(EDI_IMPACT,"#D6604D")]:
        ax_ts.axhline(th, color=col, linewidth=0.7, linestyle="--")
    ax_ts.set_ylim(0,1); ax_ts.set_xlabel("Year"); ax_ts.set_ylabel("Mean EDI")
    ax_ts.set_title("Monthly Mean EDI by Coast Region", fontweight="bold", fontsize=11)
    ax_ts.legend(fontsize=9, ncol=4)

    # Rows 2–3: severity stacks (2 regions per row)
    for i, name in enumerate(region_names):
        row = 2 + i // 2; col = i % 2
        ax_sv = fig6.add_subplot(gs[row, col])
        sub   = [r for r in monthly_rows if r[2] == name]
        if not sub: ax_sv.set_title(display_coast(name)); continue
        xs = [ym2date(r[0],r[1]) for r in sub]
        base = np.zeros(len(sub))
        for sname, cidx in zip(sev_names, [6,7,8,9]):
            vals = np.array([r[cidx] for r in sub])
            ax_sv.fill_between(xs, base, base+vals,
                               color=SEVERITY_COLORS[sname], alpha=0.85, label=sname)
            base += vals
        ax_sv.set_ylim(0,100); ax_sv.set_ylabel("% pixels"); ax_sv.set_title(display_coast(name), fontweight="bold")
        ax_sv.yaxis.set_major_formatter(PercentFormatter(xmax=100))
        if row == 3: ax_sv.set_xlabel("Year")

    handles = [mpatches.Patch(color=SEVERITY_COLORS[s], label=s) for s in sev_names]
    fig6.legend(handles=handles, loc="lower center", ncol=4, fontsize=10,
                bbox_to_anchor=(0.5, -0.01))
    fig6.suptitle(
        "Ecological Drought Index (EDI) — SPEI-3 Spatial Application\n"
        "US Coastline + 80 km Buffer  |  Pooled Logistic Model  |  EDI = P(WUE_T impairment | SPEI-3)",
        fontweight="bold", fontsize=13, y=1.01
    )
    fig6.savefig(os.path.join(OUTPUT_DIR,"EDI_plot_06_composite.png"),
                 dpi=300, bbox_inches="tight")
    plt.close(fig6)

    print(f"\nAll figures saved to {OUTPUT_DIR}")
    print(f"Done: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 68)


if __name__ == "__main__":
    main()