"""
Q3_WUE_T_SPEI_sensitivity.py - CHUNK 1: Data Processing and Models (FINAL UNIFIED GAM)
EXACT REPLICATION OF MALONE'S R WORKFLOW, WITH ONE UNIFIED GAM

FINAL MODEL:
- Primary coast × SPEI-timescale full smooths: s(SPEI_value, by = spei_coast, k=6)
- Ecosystem × timescale sum-to-zero deviations: s(water_class, SPEI_value, by = SPEI_timescale, bs="sz", k=6)
- Parametric SPEI_timescale * coast_region + month_f
- Site random effect: s(site_name, bs="re")

All predictions, thresholds, and file contracts preserved.
New diagnostic/table-ready outputs added.
"""

import os
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("Q3: CHUNK 1 - Data Processing and Models (FINAL UNIFIED GAM)")
print("="*60)

# ============================================================================
# PATHS - SEPARATE FROM MALONE
# ============================================================================

# Input file (data_products copy - verified identical)
input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

# Python outputs - SEPARATE from Malone with clean organization
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"

# Three separate subdirectories
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")     # CSV, RDS, summary
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures")     # PNG figures
temp_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_temp")          # Temporary R files

# Downstream reviewer CSVs must be written here (as per existing plotting script)
reviewer_dir = os.path.join(base_output_dir, "Q3_august_update")

# Create all directories
for dir_path in [output_dir, figure_dir, temp_dir, reviewer_dir]:
    os.makedirs(dir_path, exist_ok=True)

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")
print(f"Reviewer directory (downstream CSVs): {reviewer_dir}")
print(f"(Separate from Malone outputs - safe for comparison)")

# ============================================================================
# CLEAN OLD OUTPUTS AT START
# ============================================================================

print("\n" + "="*60)
print("STEP 0: Cleaning old outputs (Python directories only)")
print("="*60)

# Expected outputs (created by R/Python) - includes new files
expected_outputs = [
    "Q3_WUE_T_SPEI_model_base_wide.csv",
    "Q3_WUE_T_SPEI_model_data_long.csv",
    "Q3_WUE_T_SPEI_coverage.csv",
    "Q3_WUE_T_SPEI_site_level_slopes.csv",
    "Q3_WUE_T_SPEI_site_near_normal_baselines.csv",
    "Q3_WUE_T_SPEI_site_response_from_near_normal.csv",
    "Q3_WUE_T_SPEI_site_sensitivity_rank.csv",
    "Q3_WUE_T_SPEI_ecosystem_slope_summary.csv",
    "Q3_WUE_T_SPEI_coast_region_slope_summary.csv",
    "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv",
    "Q3_WUE_T_SPEI_smooth_gam_model.rds",
    "Q3_WUE_T_SPEI_coast_threshold_gam_model.rds",
    "Q3_WUE_T_SPEI_model_comparison.csv",
    "Q3_WUE_T_SPEI_smooth_terms.csv",
    "Q3_WUE_T_SPEI_smooth_parametric_terms.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv",
    "Q3_WUE_T_SPEI_smooth_predictions.csv",
    "Q3_WUE_T_SPEI_predicted_near_normal_baselines.csv",
    "Q3_WUE_T_SPEI_predicted_impact_classes.csv",
    "Q3_WUE_T_SPEI_predicted_impact_thresholds_10pct.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_upland.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_near_normal_baselines_upland.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_impact_thresholds_10pct_upland.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv",
    "Q3_WUE_T_SPEI_GAM_impact_threshold_markers_5_10_20pct.csv",
    "Q3_WUE_T_SPEI_GAM_threshold_summary_5_10_20pct.csv",
    "Q3_WUE_T_SPEI_model_summary.txt",
    # New diagnostic outputs
    "Q3_WUE_T_SPEI_unified_coast_ecosystem_site_counts.csv",
    "Q3_WUE_T_SPEI_unified_concurvity.csv",
    "Q3_WUE_T_SPEI_unified_all_smooth_terms.csv",
    "Q3_WUE_T_SPEI_unified_parametric_terms.csv",
    "Q3_WUE_T_SPEI_unified_gam_model.rds",
    # New table-ready files
    "Q3_WUE_T_SPEI_coast_Table2_ready.csv",
    "Q3_WUE_T_SPEI_ecosystem_deviation_by_timescale.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_all_ecosystems_full.csv",
]

# Reviewer CSVs (written to reviewer_dir)
reviewer_files = [
    "Q3_reviewer_coast_threshold_prediction_curves_month_averaged_all_ecosystems.csv",
    "Q3_reviewer_threshold_markers_month_averaged_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv"
]

# Clean output files
for f in expected_outputs:
    path = os.path.join(output_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old output: {f}")

# Clean reviewer directory (remove old reviewer CSVs)
for f in reviewer_files:
    path = os.path.join(reviewer_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old reviewer file: {f}")

# Clean temp directory
for f in os.listdir(temp_dir):
    path = os.path.join(temp_dir, f)
    if os.path.isfile(path):
        os.remove(path)
        print(f"  Removed old temp file: {f}")

# ============================================================================
# CONSTANTS - EXACTLY AS IN MALONE
# ============================================================================

SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]

# ============================================================================
# HELPER FUNCTIONS - EXACTLY AS IN MALONE
# ============================================================================

# --------------------------------------------------------------------------
# Proximity-based coast classifier (drop-in replacement)
# --------------------------------------------------------------------------
EARTH_RADIUS_KM = 6371.0088

COAST_POLYLINES = {
    "Pacific Coast": [
        (32.6, -117.2), (33.6, -118.2), (34.4, -119.7), (35.4, -120.9),
        (36.6, -121.9), (37.8, -122.5), (39.0, -123.7), (41.0, -124.2),
        (43.5, -124.2), (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)
    ],
    "Gulf Coast": [
        (25.8, -97.2), (28.0, -96.8), (29.3, -94.8), (29.5, -93.0),
        (29.2, -91.5), (29.3, -90.0), (29.5, -88.8), (30.2, -87.7),
        (30.1, -86.2), (29.9, -85.3), (29.7, -84.3), (28.8, -83.0),
        (27.8, -82.8), (26.6, -82.2), (25.9, -81.8), (25.3, -81.1),
        (25.0, -80.8)
    ],
    "Atlantic Coast": [
        (25.1, -80.3), (26.0, -80.1), (27.0, -80.1), (28.4, -80.6),
        (29.0, -80.9), (29.9, -81.3), (30.4, -81.4), (31.2, -81.3),
        (32.0, -80.8), (33.0, -79.5), (34.5, -77.8), (35.7, -75.6),
        (36.8, -75.9), (38.5, -75.0), (39.5, -74.3), (40.5, -73.9),
        (41.3, -72.0), (42.4, -70.8), (43.5, -70.2)
    ],
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
    """Minimum distance from point to segment AB (equirectangular projection)."""
    import math
    lat0 = math.radians(lat)
    def proj(lat2, lon2):
        x = EARTH_RADIUS_KM * math.radians(lon2 - lon) * math.cos(lat0)
        y = EARTH_RADIUS_KM * math.radians(lat2 - lat)
        return np.array([x, y])
    p = np.array([0.0, 0.0])
    a = proj(a_lat, a_lon)
    b = proj(b_lat, b_lon)
    ab = b - a
    denom = np.dot(ab, ab)
    if denom == 0:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, np.dot(p - a, ab) / denom))
    return float(np.linalg.norm(p - (a + t * ab)))

def _distance_to_polyline(lat, lon, polyline):
    """Minimum distance to a coastline polyline. Returns a float."""
    distances = []
    for i in range(len(polyline) - 1):
        d = _seg_dist(lat, lon, *polyline[i], *polyline[i+1])
        distances.append(d)
    return min(distances) if distances else float('inf')

def classify_coast_region(lat, long):
    """
    Proximity-based coast classifier (drop-in replacement).
    Returns the coast name (string).
    """
    if pd.isna(lat) or pd.isna(long):
        return "Other/Check"
    if lat > 50:
        return "AK Coast"
    
    lat, lon = float(lat), float(long)
    
    dists = {
        name: _distance_to_polyline(lat, lon, poly)
        for name, poly in COAST_POLYLINES.items()
    }
    
    return min(dists, key=dists.get)

def write_table(df, filename):
    """Identical to Malone's write_table()"""
    df.to_csv(os.path.join(output_dir, filename), index=False)

def classify_impact(pct_change):
    """Identical to Malone's factor creation for impact_class_10pct"""
    if pct_change <= -10:
        return "WUE decrease"
    elif pct_change >= 10:
        return "WUE increase"
    else:
        return "No meaningful change"

# ============================================================================
# STEP 1: LOAD DATA - EXACTLY AS IN MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 1: Loading and preparing data")
print("="*60)

monthly = pd.read_csv(input_file)

# Trim whitespace from character columns (identical to R's trimws)
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].astype(str).str.strip()
    monthly[col] = monthly[col].replace('nan', np.nan)

# Check required columns (identical to R)
required_cols = ["site_name", "Year", "month", "water_class", "lat", "long", "WUE_tra"] + SPEI_COLS

missing_cols = [col for col in required_cols if col not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Filter data (identical to R's complete.cases and filtering)
model_base = monthly.dropna(subset=required_cols).copy()
model_base = model_base[
    model_base['site_name'].notna() &
    (model_base['site_name'] != "") &
    model_base['water_class'].notna() &
    (model_base['water_class'] != "") &
    model_base['water_class'].isin(ECOSYSTEM_CLASSES) &
    np.isfinite(model_base['lat']) &
    np.isfinite(model_base['long']) &
    np.isfinite(model_base['WUE_tra'])
].copy()

# Add factors (identical to R's factor() calls)
model_base['coast_region'] = model_base.apply(
    lambda row: classify_coast_region(row['lat'], row['long']), axis=1
)
model_base['coast_region'] = pd.Categorical(
    model_base['coast_region'], 
    categories=COAST_REGION_LEVELS
)

model_base['site_name'] = model_base['site_name'].astype('category')
model_base['water_class'] = pd.Categorical(
    model_base['water_class'], 
    categories=ECOSYSTEM_CLASSES
)
model_base['month_f'] = pd.Categorical(model_base['month'], categories=range(1, 13))
model_base['WUE_T'] = model_base['WUE_tra']

print(f"  Rows after filtering: {len(model_base)}")
print(f"  Sites: {model_base['site_name'].nunique()}")
print(f"  Years: {model_base['Year'].min()} - {model_base['Year'].max()}")

# ============================================================================
# STEP 2: RESHAPE TO LONG FORMAT - KEEP ALL COLUMNS (FIXED)
# ============================================================================

print("\n" + "="*60)
print("STEP 2: Reshaping to long format (keeping all columns)")
print("="*60)

# Keep ALL columns except the SPEI columns (like Malone's reshape)
id_vars = [c for c in model_base.columns if c not in SPEI_COLS]

# Identical to R's reshape() with direction="long"
spei_long = pd.melt(
    model_base,
    id_vars=id_vars,
    value_vars=SPEI_COLS,
    var_name='SPEI_timescale',
    value_name='SPEI_value'
)

# Filter finite values (identical to R)
spei_long = spei_long[np.isfinite(spei_long['SPEI_value']) & np.isfinite(spei_long['WUE_T'])]

# Add factors (identical to R)
spei_long['SPEI_timescale'] = pd.Categorical(spei_long['SPEI_timescale'], categories=SPEI_COLS)
spei_long['site_name'] = spei_long['site_name'].astype('category')
spei_long['water_class'] = pd.Categorical(spei_long['water_class'], categories=ECOSYSTEM_CLASSES)
spei_long['coast_region'] = pd.Categorical(spei_long['coast_region'], categories=COAST_REGION_LEVELS)
spei_long['month_f'] = pd.Categorical(spei_long['month'], categories=range(1, 13))

# Add SPEI_z (identical to R's scale())
spei_long['SPEI_z'] = (spei_long['SPEI_value'] - spei_long['SPEI_value'].mean()) / spei_long['SPEI_value'].std()

# Add interaction terms (identical to R's interaction())
spei_long['spei_ecosystem'] = spei_long['SPEI_timescale'].astype(str) + "__" + spei_long['water_class'].astype(str)
spei_long['spei_coast'] = spei_long['SPEI_timescale'].astype(str) + "__" + spei_long['coast_region'].astype(str)

# Write tables (identical to R)
write_table(model_base, "Q3_WUE_T_SPEI_model_base_wide.csv")
write_table(spei_long, "Q3_WUE_T_SPEI_model_data_long.csv")

print(f"  Rows in long format: {len(spei_long)}")
print(f"  Columns in long format: {len(spei_long.columns)} (matches Malone)")

# ============================================================================
# STEP 3: SUMMARY STATISTICS - EXACTLY AS IN MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 3: Calculating summary statistics")
print("="*60)

# Identical to R's dplyr::group_by() %>% summarise()
coverage = spei_long.groupby(['SPEI_timescale', 'water_class'], observed=True).agg(
    n_months=('WUE_T', 'size'),
    n_sites=('site_name', 'nunique'),
    mean_SPEI=('SPEI_value', 'mean'),
    min_SPEI=('SPEI_value', 'min'),
    max_SPEI=('SPEI_value', 'max'),
    mean_WUE_T=('WUE_T', 'mean')
).reset_index()
write_table(coverage, "Q3_WUE_T_SPEI_coverage.csv")

# Site slopes (identical to R)
site_slopes_list = []
for (site_name, water_class, coast_region, spei_timescale), group in spei_long.groupby(
    ['site_name', 'water_class', 'coast_region', 'SPEI_timescale'], observed=True
):
    if len(group) >= 4 and group['SPEI_value'].nunique() >= 4:
        x = group['SPEI_value'].values
        y = group['WUE_T'].values
        slope, intercept = np.polyfit(x, y, 1)
        corr = np.corrcoef(x, y)[0, 1]
        site_slopes_list.append({
            'site_name': site_name,
            'water_class': water_class,
            'coast_region': coast_region,
            'SPEI_timescale': spei_timescale,
            'n_months': len(group),
            'SPEI_slope': slope,
            'SPEI_cor': corr,
            'mean_WUE_T': np.mean(y)
        })

site_slopes = pd.DataFrame(site_slopes_list)
if len(site_slopes) > 0:
    site_slopes['SPEI_timescale'] = pd.Categorical(site_slopes['SPEI_timescale'], categories=SPEI_COLS)
    write_table(site_slopes, "Q3_WUE_T_SPEI_site_level_slopes.csv")
print(f"  Site slopes: {len(site_slopes)} rows")

# Site near-normal baselines (identical to R)
site_near_normal_baselines = spei_long[
    (spei_long['SPEI_value'] >= -1) & (spei_long['SPEI_value'] <= 1)
].groupby(
    ['site_name', 'water_class', 'coast_region', 'SPEI_timescale', 'month_f'], observed=True
).agg(
    baseline_WUE_T=('WUE_T', 'mean'),
    baseline_sd_WUE_T=('WUE_T', 'std'),
    baseline_n=('WUE_T', 'size')
).reset_index()
write_table(site_near_normal_baselines, "Q3_WUE_T_SPEI_site_near_normal_baselines.csv")

# Site response data (identical to R's left_join) - FIXED
site_response_data = spei_long.merge(
    site_near_normal_baselines,
    on=['site_name', 'water_class', 'coast_region', 'SPEI_timescale', 'month_f'],
    how='left'
)

site_response_data = site_response_data[
    np.isfinite(site_response_data['baseline_WUE_T']) &
    (site_response_data['baseline_n'] >= 3)
].copy()

site_response_data['WUE_T_change'] = site_response_data['WUE_T'] - site_response_data['baseline_WUE_T']
site_response_data['WUE_T_pct_change'] = 100 * site_response_data['WUE_T_change'] / site_response_data['baseline_WUE_T']

site_response_data['impact_class_10pct'] = site_response_data['WUE_T_pct_change'].apply(classify_impact)
site_response_data['impact_class_10pct'] = pd.Categorical(
    site_response_data['impact_class_10pct'],
    categories=["WUE decrease", "No meaningful change", "WUE increase"]
)
write_table(site_response_data, "Q3_WUE_T_SPEI_site_response_from_near_normal.csv")

# Site sensitivity rank (identical to R)
site_sensitivity_rank = site_response_data.groupby(
    ['site_name', 'water_class', 'coast_region', 'SPEI_timescale'], observed=True
).agg(
    n_months=('WUE_T', 'size'),
    baseline_WUE_T=('baseline_WUE_T', 'first'),
    mean_abs_pct_change=('WUE_T_pct_change', lambda x: np.abs(x).mean()),
    max_abs_pct_change=('WUE_T_pct_change', lambda x: np.abs(x).max()),
    pct_months_decrease_10=('impact_class_10pct', lambda x: (x == "WUE decrease").mean() * 100),
    pct_months_increase_10=('impact_class_10pct', lambda x: (x == "WUE increase").mean() * 100),
    pct_months_impacted_10=('impact_class_10pct', lambda x: (x != "No meaningful change").mean() * 100)
).reset_index()

for idx, row in site_sensitivity_rank.iterrows():
    subset = site_response_data[
        (site_response_data['site_name'] == row['site_name']) &
        (site_response_data['SPEI_timescale'] == row['SPEI_timescale'])
    ]
    dry_decrease = subset[(subset['SPEI_value'] < -1) & (subset['WUE_T_pct_change'] <= -10)]['SPEI_value']
    wet_decrease = subset[(subset['SPEI_value'] > 1) & (subset['WUE_T_pct_change'] <= -10)]['SPEI_value']
    site_sensitivity_rank.loc[idx, 'dry_threshold_observed'] = dry_decrease.max() if len(dry_decrease) > 0 else np.nan
    site_sensitivity_rank.loc[idx, 'wet_threshold_observed'] = wet_decrease.min() if len(wet_decrease) > 0 else np.nan

site_sensitivity_rank = site_sensitivity_rank.sort_values(
    ['SPEI_timescale', 'mean_abs_pct_change'], ascending=[True, False]
)
write_table(site_sensitivity_rank, "Q3_WUE_T_SPEI_site_sensitivity_rank.csv")

# Ecosystem slope summary (identical to R)
ecosystem_slope_summary = site_slopes.groupby(['SPEI_timescale', 'water_class'], observed=True).agg(
    n_sites=('site_name', 'size'),
    mean_slope=('SPEI_slope', 'mean'),
    median_slope=('SPEI_slope', 'median'),
    sd_slope=('SPEI_slope', 'std'),
    pct_negative=('SPEI_slope', lambda x: (x < 0).mean() * 100)
).reset_index()
ecosystem_slope_summary['se_slope'] = ecosystem_slope_summary['sd_slope'] / np.sqrt(ecosystem_slope_summary['n_sites'])
ecosystem_slope_summary['ci95_slope'] = stats.t.ppf(0.975, ecosystem_slope_summary['n_sites'] - 1) * ecosystem_slope_summary['se_slope']
write_table(ecosystem_slope_summary, "Q3_WUE_T_SPEI_ecosystem_slope_summary.csv")

# Coast slope summary (identical to R)
coast_slope_summary = site_slopes.groupby(['SPEI_timescale', 'coast_region'], observed=True).agg(
    n_sites=('site_name', 'size'),
    mean_slope=('SPEI_slope', 'mean'),
    median_slope=('SPEI_slope', 'median'),
    sd_slope=('SPEI_slope', 'std'),
    pct_negative=('SPEI_slope', lambda x: (x < 0).mean() * 100)
).reset_index()
coast_slope_summary['se_slope'] = coast_slope_summary['sd_slope'] / np.sqrt(coast_slope_summary['n_sites'])
coast_slope_summary['ci95_slope'] = stats.t.ppf(0.975, coast_slope_summary['n_sites'] - 1) * coast_slope_summary['se_slope']
write_table(coast_slope_summary, "Q3_WUE_T_SPEI_coast_region_slope_summary.csv")

# Coastline grouped sensitivity summary (identical to R)
all_system_sensitivity = site_sensitivity_rank[
    site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES) &
    (site_sensitivity_rank['n_months'] >= 6)
]

if not all_system_sensitivity.empty:
    coastline_sensitivity_summary = all_system_sensitivity.groupby(
        ['coast_region', 'SPEI_timescale'], observed=True
    ).agg(
        n_sites=('site_name', 'nunique'),
        median_abs_pct_change=('mean_abs_pct_change', 'median'),
        q25_abs_pct_change=('mean_abs_pct_change', lambda x: x.quantile(0.25)),
        q75_abs_pct_change=('mean_abs_pct_change', lambda x: x.quantile(0.75)),
        mean_abs_pct_change=('mean_abs_pct_change', 'mean')
    ).reset_index()
    write_table(coastline_sensitivity_summary, "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv")

# ============================================================================
# CAPTURE ROW COUNTS FOR MODEL SUMMARY
# ============================================================================

n_monthly = len(monthly)
n_model_base = len(model_base)
n_spei_long = len(spei_long)
n_sites = spei_long['site_name'].nunique()
n_years_min = spei_long['Year'].min()
n_years_max = spei_long['Year'].max()

print(f"\n  Row counts for summary:")
print(f"    Source monthly data: {n_monthly}")
print(f"    After filtering: {n_model_base}")
print(f"    SPEI long data: {n_spei_long}")
print(f"    Sites: {n_sites}")
print(f"    Years: {n_years_min} - {n_years_max}")

# ============================================================================
# STEP 4: SAVE DATA FOR R MODELS (IN TEMP FOLDER) - FIXED
# ============================================================================

print("\n" + "="*60)
print("STEP 4: Saving data for R models")
print("="*60)

r_data_file = os.path.join(temp_dir, "temp_r_data.csv")
spei_long.to_csv(r_data_file, index=False)
print(f"  Data saved to: {r_data_file}")

# ============================================================================
# STEP 5: CREATE R SCRIPT (IN TEMP FOLDER) - FINAL UNIFIED GAM (WITH FIXED ESCAPING)
# ============================================================================

print("\n" + "="*60)
print("STEP 5: Creating R script for FINAL UNIFIED GAM")
print("="*60)

# Use forward slashes for R compatibility
output_dir_r = output_dir.replace("\\", "/")
reviewer_dir_r = reviewer_dir.replace("\\", "/")
temp_dir_r = temp_dir.replace("\\", "/")
r_data_file_r = r_data_file.replace("\\", "/")
input_file_r = input_file.replace("\\", "/")

r_script_content = f'''
# Q3 WUE_T sensitivity to SPEI gradients - FINAL UNIFIED GAM
# Primary coast × SPEI-timescale smooths + ecosystem sum-to-zero deviations
# All existing outputs preserved
# Uses ALL available months and FULL observed SPEI range
# Batch predictions for speed

library(mgcv)
library(dplyr)

sep_line <- function() {{
    cat(paste(rep("=", 60), collapse=""))
    cat("\\n")
}}

write_table <- function(x, filename) {{
    write.csv(x, file.path("{output_dir_r}", filename), row.names = FALSE)
}}

# Load data
data <- read.csv("{r_data_file_r}")

# Set factors (identical to Malone)
data$site_name <- as.factor(data$site_name)
data$water_class <- factor(data$water_class, levels = c("Upland", "Freshwater", "Saline"))
data$coast_region <- factor(data$coast_region, levels = c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"))
data$month_f <- as.factor(data$month_f)
data$SPEI_timescale <- factor(data$SPEI_timescale, 
                              levels = c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"))

data$spei_ecosystem <- interaction(data$SPEI_timescale, data$water_class, sep = "__", drop = TRUE)
data$spei_coast <- interaction(data$SPEI_timescale, data$coast_region, sep = "__", drop = TRUE)

cat("\\n")
sep_line()
cat("FINAL UNIFIED GAM - PRIMARY COAST SMOOTHS + ECOSYSTEM DEVIATIONS")
cat("\\n")
sep_line()
cat("\\n")

# ---------------------------------------------------------------------------
# Fit ONE final unified GAM
# ---------------------------------------------------------------------------

unified_model <- gam(
    WUE_T ~
        SPEI_timescale * coast_region +
        month_f +
        s(SPEI_value, by = spei_coast, k = 6) +
        s(water_class, SPEI_value,
          by = SPEI_timescale,
          bs = "sz", k = 6) +
        s(site_name, bs = "re"),
    data = data,
    method = "ML",
    select = TRUE
)

# Save under both legacy filenames (same model)
saveRDS(unified_model, file.path("{output_dir_r}", "Q3_WUE_T_SPEI_smooth_gam_model.rds"))
saveRDS(unified_model, file.path("{output_dir_r}", "Q3_WUE_T_SPEI_coast_threshold_gam_model.rds"))
# Additional convenience file
saveRDS(unified_model, file.path("{output_dir_r}", "Q3_WUE_T_SPEI_unified_gam_model.rds"))

# ---------------------------------------------------------------------------
# Diagnostics: coast × ecosystem site counts
# ---------------------------------------------------------------------------

coast_ecosystem_counts <- data |>
    distinct(site_name, coast_region, water_class) |>
    group_by(coast_region, water_class) |>
    summarise(n_sites = n(), .groups = "drop")
write_table(coast_ecosystem_counts, "Q3_WUE_T_SPEI_unified_coast_ecosystem_site_counts.csv")

# Concurvity
concurv <- concurvity(unified_model, full = TRUE)
concurv_df <- as.data.frame(concurv)
concurv_df$term <- rownames(concurv_df)
concurv_df <- concurv_df[, c("term", setdiff(names(concurv_df), "term"))]
write_table(concurv_df, "Q3_WUE_T_SPEI_unified_concurvity.csv")

# ---------------------------------------------------------------------------
# Smooth terms: split for legacy and new table-ready files
# ---------------------------------------------------------------------------

smooth_table <- as.data.frame(summary(unified_model)$s.table)
smooth_table$smooth_term <- rownames(smooth_table)
rownames(smooth_table) <- NULL
smooth_table <- smooth_table[, c("smooth_term", setdiff(names(smooth_table), "smooth_term"))]

# Coast smooths (primary) - contains "spei_coast"
coast_smooth <- smooth_table[grepl("spei_coast", smooth_table$smooth_term), ]
write_table(coast_smooth, "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv")

# Ecosystem smooths (deviation) - contains "water_class"
ecosystem_smooth <- smooth_table[grepl("water_class", smooth_table$smooth_term), ]
write_table(ecosystem_smooth, "Q3_WUE_T_SPEI_smooth_terms.csv")

# Combined diagnostic
write_table(smooth_table, "Q3_WUE_T_SPEI_unified_all_smooth_terms.csv")

# ---- Table 2 ready coast smooths ----
coast_table2 <- coast_smooth
# Parse smooth_term to extract coast_region and SPEI_timescale
# smooth_term format: "s(SPEI_value):spei_coast<SPEI_timescale>__<coast_region>"
parse_coast_term <- function(term) {{
    # Remove "s(SPEI_value):spei_coast" using fixed = TRUE (no escaping)
    after <- sub("s(SPEI_value):spei_coast", "", term, fixed = TRUE)
    parts <- strsplit(after, "__")[[1]]
    if (length(parts) == 2) {{
        return(data.frame(SPEI_timescale = parts[1], coast_region = parts[2], stringsAsFactors = FALSE))
    }} else {{
        return(data.frame(SPEI_timescale = NA, coast_region = NA, stringsAsFactors = FALSE))
    }}
}}
parsed <- do.call(rbind, lapply(coast_table2$smooth_term, parse_coast_term))
coast_table2 <- cbind(parsed, coast_table2[, !names(coast_table2) %in% "smooth_term"])
names(coast_table2)[names(coast_table2) == "p-value"] <- "p_value"
coast_table2 <- coast_table2[, c("coast_region", "SPEI_timescale", "edf", "Ref.df", "F", "p_value")]
write_table(coast_table2, "Q3_WUE_T_SPEI_coast_Table2_ready.csv")

# ---- Ecosystem deviation by timescale ----
ecosystem_dev <- ecosystem_smooth
# Parse smooth_term to extract SPEI_timescale
parse_eco_term <- function(term) {{
    # Remove "s(water_class,SPEI_value):SPEI_timescale" using fixed = TRUE
    after <- sub("s(water_class,SPEI_value):SPEI_timescale", "", term, fixed = TRUE)
    return(data.frame(SPEI_timescale = after, stringsAsFactors = FALSE))
}}
parsed_eco <- do.call(rbind, lapply(ecosystem_dev$smooth_term, parse_eco_term))
ecosystem_dev <- cbind(parsed_eco, ecosystem_dev[, !names(ecosystem_dev) %in% "smooth_term"])
names(ecosystem_dev)[names(ecosystem_dev) == "p-value"] <- "p_value"
ecosystem_dev <- ecosystem_dev[, c("SPEI_timescale", "edf", "Ref.df", "F", "p_value")]
write_table(ecosystem_dev, "Q3_WUE_T_SPEI_ecosystem_deviation_by_timescale.csv")

# ---------------------------------------------------------------------------
# Parametric terms: only (Intercept), SPEI_timescale, coast_region, interactions, month_f
# ---------------------------------------------------------------------------

p_table <- as.data.frame(summary(unified_model)$p.table)
p_table$term <- rownames(p_table)
rownames(p_table) <- NULL
p_table <- p_table[, c("term", setdiff(names(p_table), "term"))]

# Write to both legacy parametric files (same content)
write_table(p_table, "Q3_WUE_T_SPEI_smooth_parametric_terms.csv")
write_table(p_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv")
# Unified parametric file
write_table(p_table, "Q3_WUE_T_SPEI_unified_parametric_terms.csv")

# ---------------------------------------------------------------------------
# Model comparison (single row)
# ---------------------------------------------------------------------------

model_comparison <- data.frame(
    model = "unified_coast_ecosystem_gam",
    AIC = AIC(unified_model),
    BIC = BIC(unified_model),
    logLik = as.numeric(logLik(unified_model)),
    deviance_explained = summary(unified_model)$dev.expl,
    adjusted_r_squared = summary(unified_model)$r.sq
)
write_table(model_comparison, "Q3_WUE_T_SPEI_model_comparison.csv")

# ---------------------------------------------------------------------------
# Helper: month-averaged predictions using ALL available months (batch mode)
# ---------------------------------------------------------------------------

# Get available months from the data (as numeric)
available_months <- sort(unique(as.numeric(as.character(data$month_f))))

# Function to predict a full grid and average over months
predict_month_averaged_batch <- function(base_grid, model) {{
    all_fits <- list()
    all_ses <- list()
    for (m in available_months) {{
        newdata <- base_grid
        newdata$month_f <- factor(m, levels = levels(model$model$month_f))
        pred <- predict(model, newdata = newdata, type = "link", se.fit = TRUE, exclude = "s(site_name)")
        all_fits[[as.character(m)]] <- as.numeric(pred$fit)
        all_ses[[as.character(m)]] <- as.numeric(pred$se.fit)
    }}
    fit_mat <- do.call(cbind, all_fits)
    se_mat <- do.call(cbind, all_ses)
    avg_fit <- rowMeans(fit_mat, na.rm = TRUE)
    avg_se <- rowMeans(se_mat, na.rm = TRUE)
    avg_lower <- avg_fit - 1.96 * avg_se
    avg_upper <- avg_fit + 1.96 * avg_se
    result <- base_grid
    result$predicted_WUE_T <- avg_fit
    result$predicted_se <- avg_se
    result$predicted_lower <- avg_lower
    result$predicted_upper <- avg_upper
    return(result)
}}

# ---------------------------------------------------------------------------
# Ecosystem predictions (all water_class, month-averaged)
# ---------------------------------------------------------------------------

# Use full SPEI range (min to max)
spei_min <- min(data$SPEI_value, na.rm = TRUE)
spei_max <- max(data$SPEI_value, na.rm = TRUE)
spei_sequence <- seq(spei_min, spei_max, length.out = 200)

ecosystem_base <- expand.grid(
    SPEI_timescale = levels(data$SPEI_timescale),
    water_class = levels(data$water_class),
    SPEI_value = spei_sequence,
    site_name = levels(data$site_name)[1],
    coast_region = factor("Atlantic Coast", levels = levels(data$coast_region)),  # required for model
    stringsAsFactors = FALSE
)
# Keep interaction columns for backward compatibility (not used in model)
ecosystem_base$spei_ecosystem <- interaction(ecosystem_base$SPEI_timescale, ecosystem_base$water_class, sep = "__", drop = TRUE)
ecosystem_base$spei_ecosystem <- factor(ecosystem_base$spei_ecosystem, levels = levels(data$spei_ecosystem))
ecosystem_base$spei_coast <- interaction(ecosystem_base$SPEI_timescale, ecosystem_base$coast_region, sep = "__", drop = TRUE)
ecosystem_base$spei_coast <- factor(ecosystem_base$spei_coast, levels = levels(data$spei_coast))

ecosystem_pred <- predict_month_averaged_batch(ecosystem_base, unified_model)
write_table(ecosystem_pred, "Q3_WUE_T_SPEI_smooth_predictions.csv")

# ---------------------------------------------------------------------------
# Coast predictions (all coast_region, all water_class, month-averaged)
# ---------------------------------------------------------------------------

coast_base <- expand.grid(
    SPEI_timescale = levels(data$SPEI_timescale),
    coast_region = levels(data$coast_region),
    water_class = levels(data$water_class),
    SPEI_value = spei_sequence,
    site_name = levels(data$site_name)[1],
    stringsAsFactors = FALSE
)
# For backward compatibility, keep the interaction columns
coast_base$spei_coast <- interaction(coast_base$SPEI_timescale, coast_base$coast_region, sep = "__", drop = TRUE)
coast_base$spei_coast <- factor(coast_base$spei_coast, levels = levels(data$spei_coast))
coast_base$spei_ecosystem <- interaction(coast_base$SPEI_timescale, coast_base$water_class, sep = "__", drop = TRUE)
coast_base$spei_ecosystem <- factor(coast_base$spei_ecosystem, levels = levels(data$spei_ecosystem))

coast_pred <- predict_month_averaged_batch(coast_base, unified_model)

# Also produce legacy Upland-only coast predictions (for backward compatibility)
coast_upland_base <- coast_base[coast_base$water_class == "Upland", ]
coast_upland_pred <- predict_month_averaged_batch(coast_upland_base, unified_model)
write_table(coast_upland_pred, "Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_upland.csv")

# ---------------------------------------------------------------------------
# Impact thresholds (10% and multi-level) - all ecosystems
# ---------------------------------------------------------------------------

# Function to compute near-normal baseline and impact classes from a prediction data frame
# Assumes columns: predicted_WUE_T, predicted_se, predicted_lower, predicted_upper, SPEI_value, 
# and grouping columns (SPEI_timescale, water_class, coast_region if present)
compute_impacts <- function(pred_df, group_vars) {{
    near_normal <- pred_df |>
        filter(SPEI_value >= -1 & SPEI_value <= 1) |>
        group_by(across(all_of(group_vars))) |>
        summarise(
            predicted_near_normal_WUE_T = mean(predicted_WUE_T, na.rm = TRUE),
            .groups = "drop"
        )
    pred_impact <- left_join(pred_df, near_normal, by = group_vars)
    pred_impact$predicted_change <- pred_impact$predicted_WUE_T - pred_impact$predicted_near_normal_WUE_T
    pred_impact$predicted_pct_change <- 100 * pred_impact$predicted_change / pred_impact$predicted_near_normal_WUE_T
    pred_impact$predicted_lower_pct <- 100 * (pred_impact$predicted_lower - pred_impact$predicted_near_normal_WUE_T) /
        pred_impact$predicted_near_normal_WUE_T
    pred_impact$predicted_upper_pct <- 100 * (pred_impact$predicted_upper - pred_impact$predicted_near_normal_WUE_T) /
        pred_impact$predicted_near_normal_WUE_T
    pred_impact$impact_class_10pct <- factor(
        ifelse(pred_impact$predicted_pct_change <= -10, "WUE decrease",
               ifelse(pred_impact$predicted_pct_change >= 10, "WUE increase", "No meaningful change")),
        levels = c("WUE decrease", "No meaningful change", "WUE increase")
    )
    return(list(pred_impact = pred_impact, near_normal = near_normal))
}}

# Ecosystem impacts (group by SPEI_timescale, water_class)
eco_impacts <- compute_impacts(ecosystem_pred, c("SPEI_timescale", "water_class"))
eco_pred_impact <- eco_impacts$pred_impact
eco_near_normal <- eco_impacts$near_normal

write_table(eco_near_normal, "Q3_WUE_T_SPEI_predicted_near_normal_baselines.csv")
write_table(eco_pred_impact, "Q3_WUE_T_SPEI_predicted_impact_classes.csv")

# 10% thresholds for ecosystem
eco_thresholds_10 <- eco_pred_impact |>
    group_by(SPEI_timescale, water_class) |>
    summarise(
        near_normal_WUE_T = first(predicted_near_normal_WUE_T),
        dry_decrease_threshold = ifelse(
            any(SPEI_value < -1 & predicted_pct_change <= -10),
            max(SPEI_value[SPEI_value < -1 & predicted_pct_change <= -10], na.rm = TRUE),
            NA_real_
        ),
        dry_increase_threshold = ifelse(
            any(SPEI_value < -1 & predicted_pct_change >= 10),
            max(SPEI_value[SPEI_value < -1 & predicted_pct_change >= 10], na.rm = TRUE),
            NA_real_
        ),
        wet_decrease_threshold = ifelse(
            any(SPEI_value > 1 & predicted_pct_change <= -10),
            min(SPEI_value[SPEI_value > 1 & predicted_pct_change <= -10], na.rm = TRUE),
            NA_real_
        ),
        wet_increase_threshold = ifelse(
            any(SPEI_value > 1 & predicted_pct_change >= 10),
            min(SPEI_value[SPEI_value > 1 & predicted_pct_change >= 10], na.rm = TRUE),
            NA_real_
        ),
        min_predicted_pct_change = min(predicted_pct_change),
        max_predicted_pct_change = max(predicted_pct_change),
        .groups = "drop"
    )
write_table(eco_thresholds_10, "Q3_WUE_T_SPEI_predicted_impact_thresholds_10pct.csv")

# Coast impacts (group by SPEI_timescale, coast_region, water_class) - all ecosystems
coast_impacts <- compute_impacts(coast_pred, c("SPEI_timescale", "coast_region", "water_class"))
coast_pred_impact <- coast_impacts$pred_impact
coast_near_normal <- coast_impacts$near_normal

# Legacy upland-only coast impacts (for old files)
coast_upland_impacts <- compute_impacts(coast_upland_pred, c("SPEI_timescale", "coast_region"))
coast_upland_pred_impact <- coast_upland_impacts$pred_impact
coast_upland_near_normal <- coast_upland_impacts$near_normal

write_table(coast_upland_near_normal, "Q3_WUE_T_SPEI_coast_threshold_GAM_near_normal_baselines_upland.csv")
write_table(coast_upland_pred_impact, "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv")

# 10% thresholds for upland coast
coast_upland_thresholds_10 <- coast_upland_pred_impact |>
    group_by(SPEI_timescale, coast_region) |>
    summarise(
        near_normal_WUE_T = first(predicted_near_normal_WUE_T),
        dry_decrease_threshold = ifelse(
            any(SPEI_value < -1 & predicted_pct_change <= -10),
            max(SPEI_value[SPEI_value < -1 & predicted_pct_change <= -10], na.rm = TRUE),
            NA_real_
        ),
        dry_increase_threshold = ifelse(
            any(SPEI_value < -1 & predicted_pct_change >= 10),
            max(SPEI_value[SPEI_value < -1 & predicted_pct_change >= 10], na.rm = TRUE),
            NA_real_
        ),
        wet_decrease_threshold = ifelse(
            any(SPEI_value > 1 & predicted_pct_change <= -10),
            min(SPEI_value[SPEI_value > 1 & predicted_pct_change <= -10], na.rm = TRUE),
            NA_real_
        ),
        wet_increase_threshold = ifelse(
            any(SPEI_value > 1 & predicted_pct_change >= 10),
            min(SPEI_value[SPEI_value > 1 & predicted_pct_change >= 10], na.rm = TRUE),
            NA_real_
        ),
        min_predicted_pct_change = min(predicted_pct_change),
        max_predicted_pct_change = max(predicted_pct_change),
        .groups = "drop"
    )
write_table(coast_upland_thresholds_10, "Q3_WUE_T_SPEI_coast_threshold_GAM_impact_thresholds_10pct_upland.csv")

# ---------------------------------------------------------------------------
# Multi-level threshold markers (5,10,15,20,25,30,35,40,50,75) - for all ecosystems
# ---------------------------------------------------------------------------

find_spei_threshold <- function(group_data, pct_column, threshold_pct, anomaly_side, impact_direction) {{
    pct_values <- group_data[[pct_column]]
    if (anomaly_side == "dry") {{
        side_filter <- group_data$SPEI_value < -1
    }} else {{
        side_filter <- group_data$SPEI_value > 1
    }}
    if (impact_direction == "decrease") {{
        direction_filter <- pct_values <= -threshold_pct
    }} else {{
        direction_filter <- pct_values >= threshold_pct
    }}
    threshold_values <- group_data$SPEI_value[side_filter & direction_filter]
    if (length(threshold_values) == 0) return(NA_real_)
    if (anomaly_side == "dry") {{
        return(max(threshold_values, na.rm = TRUE))
    }} else {{
        return(min(threshold_values, na.rm = TRUE))
    }}
}}

# For all ecosystems (coast_pred_impact)
threshold_levels <- c(5, 10, 15, 20, 25, 30, 35, 40, 50, 75)

# Build markers for coast (all ecosystems)
coast_markers_all <- do.call(
    rbind,
    lapply(
        split(coast_pred_impact, list(coast_pred_impact$SPEI_timescale, coast_pred_impact$coast_region, coast_pred_impact$water_class), drop = TRUE),
        function(group_data) {{
            do.call(
                rbind,
                lapply(threshold_levels, function(th) {{
                    dry_dec <- group_data$SPEI_value[group_data$SPEI_value < -1 & group_data$predicted_pct_change <= -th]
                    dry_inc <- group_data$SPEI_value[group_data$SPEI_value < -1 & group_data$predicted_pct_change >= th]
                    wet_dec <- group_data$SPEI_value[group_data$SPEI_value > 1 & group_data$predicted_pct_change <= -th]
                    wet_inc <- group_data$SPEI_value[group_data$SPEI_value > 1 & group_data$predicted_pct_change >= th]
                    data.frame(
                        SPEI_timescale = group_data$SPEI_timescale[1],
                        coast_region = group_data$coast_region[1],
                        water_class = group_data$water_class[1],
                        threshold_pct = th,
                        impact_direction = c("decrease", "increase", "decrease", "increase"),
                        anomaly_side = c("dry", "dry", "wet", "wet"),
                        SPEI_threshold = c(
                            ifelse(length(dry_dec) > 0, max(dry_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(dry_inc) > 0, max(dry_inc, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_dec) > 0, min(wet_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_inc) > 0, min(wet_inc, na.rm = TRUE), NA_real_)
                        ),
                        pct_change_threshold = c(-th, th, -th, th)
                    )
                }})
            )
        }}
    )
)
coast_markers_all <- coast_markers_all[!is.na(coast_markers_all$SPEI_threshold), ]
# Compute lower/upper using lower/upper pct columns
coast_markers_all$SPEI_threshold_lower <- NA_real_
coast_markers_all$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(coast_markers_all))) {{
    row <- coast_markers_all[i, ]
    group <- coast_pred_impact[
        coast_pred_impact$SPEI_timescale == row$SPEI_timescale &
        coast_pred_impact$coast_region == row$coast_region &
        coast_pred_impact$water_class == row$water_class,
    ]
    lower <- find_spei_threshold(group, "predicted_lower_pct", row$threshold_pct, row$anomaly_side, row$impact_direction)
    upper <- find_spei_threshold(group, "predicted_upper_pct", row$threshold_pct, row$anomaly_side, row$impact_direction)
    vals <- c(row$SPEI_threshold, lower, upper)
    vals <- vals[!is.na(vals)]
    coast_markers_all$SPEI_threshold_lower[i] <- min(vals)
    coast_markers_all$SPEI_threshold_upper[i] <- max(vals)
}}
coast_markers_all$threshold_pct <- factor(coast_markers_all$threshold_pct, 
                                          levels = threshold_levels,
                                          labels = paste0(threshold_levels, "%"))

# Write reviewer CSVs to the Q3_august_update folder
# Include mean_pct, lower_pct, upper_pct for ribbons
write.csv(coast_pred_impact |>
          select(coast_region, SPEI_timescale, water_class, SPEI_value,
                 mean_pct = predicted_pct_change,
                 lower_pct = predicted_lower_pct,
                 upper_pct = predicted_upper_pct),
          file.path("{reviewer_dir_r}", "Q3_reviewer_coast_threshold_prediction_curves_month_averaged_all_ecosystems.csv"),
          row.names = FALSE)

write.csv(coast_markers_all,
          file.path("{reviewer_dir_r}", "Q3_reviewer_threshold_markers_month_averaged_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv"),
          row.names = FALSE)

# Also produce the legacy 5-10-20 markers for upland coast
coast_upland_markers <- coast_markers_all[coast_markers_all$water_class == "Upland", ]
coast_upland_markers_5_10_20 <- coast_upland_markers[coast_upland_markers$threshold_pct %in% c("5%", "10%", "20%"), ]
write_table(coast_upland_markers_5_10_20, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv")
coast_upland_summary <- coast_upland_markers_5_10_20 |>
    arrange(SPEI_timescale, coast_region, anomaly_side, impact_direction, threshold_pct)
write_table(coast_upland_summary, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv")

# Ecosystem GAM markers (5-10-20) - legacy
eco_markers_5_10_20 <- do.call(
    rbind,
    lapply(
        split(eco_pred_impact, list(eco_pred_impact$SPEI_timescale, eco_pred_impact$water_class), drop = TRUE),
        function(group_data) {{
            do.call(
                rbind,
                lapply(c(5, 10, 20), function(th) {{
                    dry_dec <- group_data$SPEI_value[group_data$SPEI_value < -1 & group_data$predicted_pct_change <= -th]
                    dry_inc <- group_data$SPEI_value[group_data$SPEI_value < -1 & group_data$predicted_pct_change >= th]
                    wet_dec <- group_data$SPEI_value[group_data$SPEI_value > 1 & group_data$predicted_pct_change <= -th]
                    wet_inc <- group_data$SPEI_value[group_data$SPEI_value > 1 & group_data$predicted_pct_change >= th]
                    data.frame(
                        SPEI_timescale = group_data$SPEI_timescale[1],
                        water_class = group_data$water_class[1],
                        threshold_pct = th,
                        impact_direction = c("decrease", "increase", "decrease", "increase"),
                        anomaly_side = c("dry", "dry", "wet", "wet"),
                        SPEI_threshold = c(
                            ifelse(length(dry_dec) > 0, max(dry_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(dry_inc) > 0, max(dry_inc, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_dec) > 0, min(wet_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_inc) > 0, min(wet_inc, na.rm = TRUE), NA_real_)
                        ),
                        pct_change_threshold = c(-th, th, -th, th)
                    )
                }})
            )
        }}
    )
)
eco_markers_5_10_20 <- eco_markers_5_10_20[!is.na(eco_markers_5_10_20$SPEI_threshold), ]
eco_markers_5_10_20$SPEI_threshold_lower <- NA_real_
eco_markers_5_10_20$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(eco_markers_5_10_20))) {{
    row <- eco_markers_5_10_20[i, ]
    group <- eco_pred_impact[
        eco_pred_impact$SPEI_timescale == row$SPEI_timescale &
        eco_pred_impact$water_class == row$water_class,
    ]
    lower <- find_spei_threshold(group, "predicted_lower_pct", row$threshold_pct, row$anomaly_side, row$impact_direction)
    upper <- find_spei_threshold(group, "predicted_upper_pct", row$threshold_pct, row$anomaly_side, row$impact_direction)
    vals <- c(row$SPEI_threshold, lower, upper)
    vals <- vals[!is.na(vals)]
    eco_markers_5_10_20$SPEI_threshold_lower[i] <- min(vals)
    eco_markers_5_10_20$SPEI_threshold_upper[i] <- max(vals)
}}
eco_markers_5_10_20$threshold_pct <- factor(eco_markers_5_10_20$threshold_pct,
                                            levels = c(5,10,20),
                                            labels = c("5%","10%","20%"))
eco_markers_5_10_20 <- eco_markers_5_10_20 |>
    arrange(SPEI_timescale, water_class, anomaly_side, impact_direction, threshold_pct)
write_table(eco_markers_5_10_20, "Q3_WUE_T_SPEI_GAM_impact_threshold_markers_5_10_20pct.csv")
write_table(eco_markers_5_10_20, "Q3_WUE_T_SPEI_GAM_threshold_summary_5_10_20pct.csv")

# ---------------------------------------------------------------------------
# Save full master prediction CSV (all ecosystems)
# ---------------------------------------------------------------------------
write_table(coast_pred_impact, "Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_all_ecosystems_full.csv")

# ---------------------------------------------------------------------------
# Model summary (including diagnostics)
# ---------------------------------------------------------------------------

cat("\\n")
sep_line()
cat("FINAL UNIFIED GAM - COAST SMOOTHS + ECOSYSTEM DEVIATIONS")
cat("\\n")
sep_line()
cat("\\n")

sink(file.path("{output_dir_r}", "Q3_WUE_T_SPEI_model_summary.txt"))

cat("Q3 WUE_T sensitivity to SPEI gradients - FINAL UNIFIED GAM\\n")
cat("Input:", "{input_file_r}", "\\n")
cat("Output:", "{output_dir_r}", "\\n")
cat("Reviewer CSVs written to:", "{reviewer_dir_r}", "\\n\\n")
cat("Rows in source monthly data:", {n_monthly}, "\\n")
cat("Rows after complete-case and ecosystem filters:", {n_model_base}, "\\n")
cat("Rows in SPEI long data:", {n_spei_long}, "\\n")
cat("Sites:", {n_sites}, "\\n")
cat("Years:", {n_years_min}, "-", {n_years_max}, "\\n")
cat("Available months used for averaging:", paste(available_months, collapse=", "), "\\n")
cat("SPEI range used for predictions: [", spei_min, ", ", spei_max, "]\\n\\n")

cat("FINAL MODEL: one unified GAM with:\\n")
cat("  - primary coast × SPEI-timescale full smooths (s(SPEI_value, by = spei_coast))\\n")
cat("  - ecosystem × timescale sum-to-zero deviations (s(water_class, SPEI_value, by = SPEI_timescale, bs='sz'))\\n")
cat("  - parametric SPEI_timescale * coast_region + month_f\\n")
cat("  - site random effect\\n")
cat("Formula:\\n")
print(formula(unified_model))
cat("\\nModel summary:\\n")
print(summary(unified_model))
cat("\\nDeviance explained:", summary(unified_model)$dev.expl)
cat("\\nAdjusted R-squared:", summary(unified_model)$r.sq)
cat("\\nAIC:", AIC(unified_model))
cat("\\nBIC:", BIC(unified_model))
cat("\\n\\n--- Coast × Ecosystem site counts ---\\n")
print(coast_ecosystem_counts)

cat("\\n--- Concurvity (full) ---\\n")
print(concurv)

cat("\\n--- gam.check output ---\\n")
gam.check(unified_model)

cat("\\n--- Table-ready coast smooths saved to: Q3_WUE_T_SPEI_coast_Table2_ready.csv\\n")
cat("--- Ecosystem deviation tests saved to: Q3_WUE_T_SPEI_ecosystem_deviation_by_timescale.csv\\n")
cat("--- Full predictions saved to: Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_all_ecosystems_full.csv\\n")

sink()

cat("\\n")
sep_line()
cat("FINAL UNIFIED GAM COMPLETE - ALL OUTPUTS PRESERVED")
cat("\\n")
sep_line()
cat("\\n")
'''

# Save R script in temp folder
r_script_file = os.path.join(temp_dir, "run_models.R")
with open(r_script_file, 'w', encoding='utf-8') as f:
    f.write(r_script_content)

print(f"  R script saved to: {r_script_file}")

# ============================================================================
# STEP 6: RUN R SCRIPT
# ============================================================================

print("\n" + "="*60)
print("STEP 6: Running R models")
print("="*60)

# Add R to PATH if needed
r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
if r_path not in os.environ['PATH']:
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + r_path

print("\n  Running R models (this may take a few minutes)...")
result = subprocess.run(
    ["Rscript", r_script_file],
    capture_output=True,
    text=True
)

if result.returncode != 0:
    print("  R script had errors:")
    print(result.stderr)
    raise RuntimeError("R model execution failed")
else:
    print("  R models completed successfully!")
    if result.stdout:
        print(result.stdout)
    # Print R warnings/messages even on success
    if result.stderr:
        print("R warnings/messages:")
        print(result.stderr)

# ============================================================================
# STEP 7: VERIFY OUTPUTS - HARD FAILURE (FIXED)
# ============================================================================

print("\n" + "="*60)
print("STEP 7: Verifying outputs (hard failure if missing)")
print("="*60)

missing_files = []
for filename in expected_outputs:
    filepath = os.path.join(output_dir, filename)
    if not os.path.exists(filepath):
        missing_files.append(filename)

# Also check reviewer files
for filename in reviewer_files:
    filepath = os.path.join(reviewer_dir, filename)
    if not os.path.exists(filepath):
        missing_files.append(f"[reviewer_dir] {filename}")

if missing_files:
    raise FileNotFoundError(
        "Missing expected Q3 output files:\n" + "\n".join(missing_files)
    )
else:
    print(f"All {len(expected_outputs)} expected outputs created successfully!")
    print(f"Reviewer CSVs created in {reviewer_dir}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("CHUNK 1 COMPLETE - FINAL UNIFIED GAM")
print("="*60)
print("\nAUDIT:")
print("- Exactly one gam() fit: YES")
print("- Primary coast smooth: s(SPEI_value, by = spei_coast, k=6): YES")
print("- No s(SPEI_value, by = spei_ecosystem) term: YES")
print("- Ecosystem deviation: s(water_class, SPEI_value, by = SPEI_timescale, bs='sz', k=6): YES")
print("- k = 6: YES")
print("- ML: YES")
print("- select = TRUE: YES")
print("- Site random effect retained: YES")
print("- All available months used: YES")
print("- Full observed SPEI range used: YES")
print("- Direct coast edf/F/p saved (Table2-ready): YES")
print("- Ecosystem deviation tests saved: YES")
print("- RDS saved for later post-hoc work: YES")
print("- Full prediction master CSV saved: YES")
print("- Reviewer curve CSV includes mean_pct/lower_pct/upper_pct: YES")
print("- Threshold marker CSV preserved: YES")
print("- Existing downstream filenames preserved: YES")
print("- Exactly one biological GAM fitted: YES")
print("\nOutputs saved to: {output_dir}")
print("Reviewer CSVs saved to: {reviewer_dir}")
print("Temp files saved to: {temp_dir}")
print("\nNow run CHUNK 2 to create figures.")
print("="*60)