"""
Q3_WUE_T_SPEI_sensitivity.py - CHUNK 1: Data Processing and Models (FIXED)
EXACT REPLICATION OF MALONE'S R WORKFLOW

FIXES:
1. model_data_long keeps ALL original columns (like Malone)
2. site_response_data uses LEFT JOIN then filter (like Malone)
3. Temp files go to temp_dir (not output_dir)
4. Hard failure for missing outputs
5. Clean directory structure (outputs/figures/temp)
"""

import os
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("Q3: CHUNK 1 - Data Processing and Models (FIXED)")
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

# Create all directories
for dir_path in [output_dir, figure_dir, temp_dir]:
    os.makedirs(dir_path, exist_ok=True)

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")
print(f"(Separate from Malone outputs - safe for comparison)")

# ============================================================================
# CLEAN OLD OUTPUTS AT START
# ============================================================================

print("\n" + "="*60)
print("STEP 0: Cleaning old outputs (Python directories only)")
print("="*60)

# Expected outputs (created by R/Python)
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
    "Q3_WUE_T_SPEI_model_summary.txt"
]

# Clean output files
for f in expected_outputs:
    path = os.path.join(output_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old output: {f}")

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
    # This explicitly returns a float, never None
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
    # Identical to R's filter(n_distinct(SPEI_value) >= 4)
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
    how='left'  # FIXED: left join, then filter (matches Malone)
)

# Filter finite values (identical to R)
site_response_data = site_response_data[
    np.isfinite(site_response_data['baseline_WUE_T']) &
    (site_response_data['baseline_n'] >= 3)
].copy()

# Calculate changes (identical to R)
site_response_data['WUE_T_change'] = site_response_data['WUE_T'] - site_response_data['baseline_WUE_T']
site_response_data['WUE_T_pct_change'] = 100 * site_response_data['WUE_T_change'] / site_response_data['baseline_WUE_T']

# Add impact class (identical to R's factor)
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

# Add thresholds (identical to R)
for idx, row in site_sensitivity_rank.iterrows():
    subset = site_response_data[
        (site_response_data['site_name'] == row['site_name']) &
        (site_response_data['SPEI_timescale'] == row['SPEI_timescale'])
    ]
    dry_decrease = subset[(subset['SPEI_value'] < -1) & (subset['WUE_T_pct_change'] <= -10)]['SPEI_value']
    wet_decrease = subset[(subset['SPEI_value'] > 1) & (subset['WUE_T_pct_change'] <= -10)]['SPEI_value']
    site_sensitivity_rank.loc[idx, 'dry_threshold_observed'] = dry_decrease.max() if len(dry_decrease) > 0 else np.nan
    site_sensitivity_rank.loc[idx, 'wet_threshold_observed'] = wet_decrease.min() if len(wet_decrease) > 0 else np.nan

# Arrange (identical to R's arrange)
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

# Use temp_dir for R data (FIXED)
r_data_file = os.path.join(temp_dir, "temp_r_data.csv")
spei_long.to_csv(r_data_file, index=False)
print(f"  Data saved to: {r_data_file}")

# ============================================================================
# STEP 5: CREATE R SCRIPT (IN TEMP FOLDER) - FIXED
# ============================================================================

print("\n" + "="*60)
print("STEP 5: Creating R script for GAM models")
print("="*60)

# Use forward slashes for R compatibility
output_dir_r = output_dir.replace("\\", "/")
temp_dir_r = temp_dir.replace("\\", "/")
r_data_file_r = r_data_file.replace("\\", "/")
input_file_r = input_file.replace("\\", "/")

r_script_content = f'''
# Q3 WUE_T sensitivity to SPEI gradients - EXACT MALONE REPLICATION
# This R script replicates Malone's GAM models exactly

library(mgcv)
library(dplyr)

# Helper functions (identical to Malone)
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
cat("GAM MODELS - EXACT MALONE REPLICATION")
cat("\\n")
sep_line()
cat("\\n")

# ---------------------------------------------------------------------------
# Models (identical to Malone)
# ---------------------------------------------------------------------------

smooth_model <- gam(
    WUE_T ~
        SPEI_timescale * water_class +
        month_f +
        s(SPEI_value, by = spei_ecosystem, k = 6) +
        s(site_name, bs = "re"),
    data = data,
    method = "ML",
    select = TRUE
)

saveRDS(smooth_model, file.path("{output_dir_r}", "Q3_WUE_T_SPEI_smooth_gam_model.rds"))

smooth_coast_model <- gam(
    WUE_T ~
        SPEI_timescale * coast_region +
        water_class +
        month_f +
        s(SPEI_value, by = spei_coast, k = 6) +
        s(site_name, bs = "re"),
    data = data,
    method = "ML",
    select = TRUE
)

saveRDS(smooth_coast_model, file.path("{output_dir_r}", "Q3_WUE_T_SPEI_coast_threshold_gam_model.rds"))

# ---------------------------------------------------------------------------
# Model summaries (identical to Malone)
# ---------------------------------------------------------------------------

model_comparison <- data.frame(
    model = c("ecosystem_smooth_gam", "coast_threshold_gam"),
    AIC = c(AIC(smooth_model), AIC(smooth_coast_model)),
    BIC = c(BIC(smooth_model), BIC(smooth_coast_model)),
    logLik = c(as.numeric(logLik(smooth_model)), as.numeric(logLik(smooth_coast_model))),
    deviance_explained = c(summary(smooth_model)$dev.expl, summary(smooth_coast_model)$dev.expl),
    adjusted_r_squared = c(summary(smooth_model)$r.sq, summary(smooth_coast_model)$r.sq)
)
write_table(model_comparison, "Q3_WUE_T_SPEI_model_comparison.csv")

smooth_table <- as.data.frame(summary(smooth_model)$s.table)
smooth_table$smooth_term <- rownames(smooth_table)
rownames(smooth_table) <- NULL
smooth_table <- smooth_table[, c("smooth_term", setdiff(names(smooth_table), "smooth_term"))]
write_table(smooth_table, "Q3_WUE_T_SPEI_smooth_terms.csv")

parametric_table <- as.data.frame(summary(smooth_model)$p.table)
parametric_table$term <- rownames(parametric_table)
rownames(parametric_table) <- NULL
parametric_table <- parametric_table[, c("term", setdiff(names(parametric_table), "term"))]
write_table(parametric_table, "Q3_WUE_T_SPEI_smooth_parametric_terms.csv")

coast_smooth_table <- as.data.frame(summary(smooth_coast_model)$s.table)
coast_smooth_table$smooth_term <- rownames(coast_smooth_table)
rownames(coast_smooth_table) <- NULL
coast_smooth_table <- coast_smooth_table[, c("smooth_term", setdiff(names(coast_smooth_table), "smooth_term"))]
write_table(coast_smooth_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv")

coast_parametric_table <- as.data.frame(summary(smooth_coast_model)$p.table)
coast_parametric_table$term <- rownames(coast_parametric_table)
rownames(coast_parametric_table) <- NULL
coast_parametric_table <- coast_parametric_table[, c("term", setdiff(names(coast_parametric_table), "term"))]
write_table(coast_parametric_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv")

# ---------------------------------------------------------------------------
# Predictions (identical to Malone)
# ---------------------------------------------------------------------------

cat("\\n--- Generating predictions ---\\n")

spei_sequence <- seq(quantile(data$SPEI_value, 0.02), quantile(data$SPEI_value, 0.98), length.out = 120)

prediction_grid <- expand.grid(
    SPEI_timescale = levels(data$SPEI_timescale),
    water_class = levels(data$water_class),
    SPEI_value = spei_sequence,
    month_f = factor("7", levels = levels(data$month_f)),
    site_name = levels(data$site_name)[1]
)
prediction_grid$spei_ecosystem <- interaction(prediction_grid$SPEI_timescale, prediction_grid$water_class, sep = "__", drop = TRUE)
prediction_grid$spei_ecosystem <- factor(prediction_grid$spei_ecosystem, levels = levels(data$spei_ecosystem))

smooth_pred <- predict(smooth_model, newdata = prediction_grid, type = "link", se.fit = TRUE, exclude = "s(site_name)")
prediction_grid$predicted_WUE_T <- as.numeric(smooth_pred$fit)
prediction_grid$predicted_se <- as.numeric(smooth_pred$se.fit)
prediction_grid$predicted_lower <- prediction_grid$predicted_WUE_T - 1.96 * prediction_grid$predicted_se
prediction_grid$predicted_upper <- prediction_grid$predicted_WUE_T + 1.96 * prediction_grid$predicted_se
write_table(prediction_grid, "Q3_WUE_T_SPEI_smooth_predictions.csv")

coast_prediction_grid <- expand.grid(
    SPEI_timescale = levels(data$SPEI_timescale),
    coast_region = levels(data$coast_region),
    SPEI_value = spei_sequence,
    water_class = factor("Upland", levels = levels(data$water_class)),
    month_f = factor("7", levels = levels(data$month_f)),
    site_name = levels(data$site_name)[1]
)
coast_prediction_grid$spei_coast <- interaction(coast_prediction_grid$SPEI_timescale, coast_prediction_grid$coast_region, sep = "__", drop = TRUE)
coast_prediction_grid$spei_coast <- factor(coast_prediction_grid$spei_coast, levels = levels(data$spei_coast))

coast_smooth_pred <- predict(smooth_coast_model, newdata = coast_prediction_grid, type = "link", se.fit = TRUE, exclude = "s(site_name)")
coast_prediction_grid$predicted_WUE_T <- as.numeric(coast_smooth_pred$fit)
coast_prediction_grid$predicted_se <- as.numeric(coast_smooth_pred$se.fit)
coast_prediction_grid$predicted_lower <- coast_prediction_grid$predicted_WUE_T - 1.96 * coast_prediction_grid$predicted_se
coast_prediction_grid$predicted_upper <- coast_prediction_grid$predicted_WUE_T + 1.96 * coast_prediction_grid$predicted_se
write_table(coast_prediction_grid, "Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_upland.csv")

# ---------------------------------------------------------------------------
# Impact classes and thresholds (identical to Malone)
# ---------------------------------------------------------------------------

near_normal_prediction_baseline <- prediction_grid[
    prediction_grid$SPEI_value >= -1 & prediction_grid$SPEI_value <= 1,
] |>
    group_by(SPEI_timescale, water_class) |>
    summarise(
        predicted_near_normal_WUE_T = mean(predicted_WUE_T),
        .groups = "drop"
    )

prediction_impact <- left_join(
    prediction_grid,
    near_normal_prediction_baseline,
    by = c("SPEI_timescale", "water_class")
)
prediction_impact$predicted_change <- prediction_impact$predicted_WUE_T -
    prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_pct_change <- 100 * prediction_impact$predicted_change /
    prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_lower_pct <- 100 *
    (prediction_impact$predicted_lower - prediction_impact$predicted_near_normal_WUE_T) /
    prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_upper_pct <- 100 *
    (prediction_impact$predicted_upper - prediction_impact$predicted_near_normal_WUE_T) /
    prediction_impact$predicted_near_normal_WUE_T
prediction_impact$impact_class_10pct <- factor(
    ifelse(
        prediction_impact$predicted_pct_change <= -10, "WUE decrease",
        ifelse(prediction_impact$predicted_pct_change >= 10, "WUE increase", "No meaningful change")
    ),
    levels = c("WUE decrease", "No meaningful change", "WUE increase")
)

impact_thresholds <- prediction_impact |>
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

write_table(near_normal_prediction_baseline, "Q3_WUE_T_SPEI_predicted_near_normal_baselines.csv")
write_table(prediction_impact, "Q3_WUE_T_SPEI_predicted_impact_classes.csv")
write_table(impact_thresholds, "Q3_WUE_T_SPEI_predicted_impact_thresholds_10pct.csv")

coast_near_normal_prediction_baseline <- coast_prediction_grid[
    coast_prediction_grid$SPEI_value >= -1 & coast_prediction_grid$SPEI_value <= 1,
] |>
    group_by(SPEI_timescale, coast_region) |>
    summarise(
        predicted_near_normal_WUE_T = mean(predicted_WUE_T),
        .groups = "drop"
    )

coast_prediction_impact <- left_join(
    coast_prediction_grid,
    coast_near_normal_prediction_baseline,
    by = c("SPEI_timescale", "coast_region")
)
coast_prediction_impact$predicted_change <- coast_prediction_impact$predicted_WUE_T -
    coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$predicted_pct_change <- 100 * coast_prediction_impact$predicted_change /
    coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$predicted_lower_pct <- 100 *
    (coast_prediction_impact$predicted_lower - coast_prediction_impact$predicted_near_normal_WUE_T) /
    coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$predicted_upper_pct <- 100 *
    (coast_prediction_impact$predicted_upper - coast_prediction_impact$predicted_near_normal_WUE_T) /
    coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$impact_class_10pct <- factor(
    ifelse(
        coast_prediction_impact$predicted_pct_change <= -10, "WUE decrease",
        ifelse(coast_prediction_impact$predicted_pct_change >= 10, "WUE increase", "No meaningful change")
    ),
    levels = c("WUE decrease", "No meaningful change", "WUE increase")
)

coast_impact_thresholds <- coast_prediction_impact |>
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

write_table(coast_near_normal_prediction_baseline, "Q3_WUE_T_SPEI_coast_threshold_GAM_near_normal_baselines_upland.csv")
write_table(coast_prediction_impact, "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv")
write_table(coast_impact_thresholds, "Q3_WUE_T_SPEI_coast_threshold_GAM_impact_thresholds_10pct_upland.csv")

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

# Coast threshold markers (identical to Malone)
coast_impact_threshold_markers <- do.call(
    rbind,
    lapply(
        split(coast_prediction_impact, list(coast_prediction_impact$SPEI_timescale, coast_prediction_impact$coast_region), drop = TRUE),
        function(group_data) {{
            do.call(
                rbind,
                lapply(c(5, 10, 20), function(threshold_pct) {{
                    dry_decrease <- group_data$SPEI_value[
                        group_data$SPEI_value < -1 & group_data$predicted_pct_change <= -threshold_pct
                    ]
                    dry_increase <- group_data$SPEI_value[
                        group_data$SPEI_value < -1 & group_data$predicted_pct_change >= threshold_pct
                    ]
                    wet_decrease <- group_data$SPEI_value[
                        group_data$SPEI_value > 1 & group_data$predicted_pct_change <= -threshold_pct
                    ]
                    wet_increase <- group_data$SPEI_value[
                        group_data$SPEI_value > 1 & group_data$predicted_pct_change >= threshold_pct
                    ]
                    data.frame(
                        SPEI_timescale = group_data$SPEI_timescale[1],
                        coast_region = group_data$coast_region[1],
                        water_class = group_data$water_class[1],
                        threshold_pct = threshold_pct,
                        impact_direction = c("decrease", "increase", "decrease", "increase"),
                        anomaly_side = c("dry", "dry", "wet", "wet"),
                        SPEI_threshold = c(
                            ifelse(length(dry_decrease) > 0, max(dry_decrease, na.rm = TRUE), NA_real_),
                            ifelse(length(dry_increase) > 0, max(dry_increase, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_decrease) > 0, min(wet_decrease, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_increase) > 0, min(wet_increase, na.rm = TRUE), NA_real_)
                        ),
                        pct_change_threshold = c(
                            -threshold_pct,
                            threshold_pct,
                            -threshold_pct,
                            threshold_pct
                        )
                    )
                }})
            )
        }}
    )
)
coast_impact_threshold_markers <- coast_impact_threshold_markers[
    !is.na(coast_impact_threshold_markers$SPEI_threshold),
]
coast_impact_threshold_markers$SPEI_threshold_lower <- NA_real_
coast_impact_threshold_markers$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(coast_impact_threshold_markers))) {{
    marker_row <- coast_impact_threshold_markers[i, ]
    marker_group <- coast_prediction_impact[
        coast_prediction_impact$SPEI_timescale == marker_row$SPEI_timescale &
            coast_prediction_impact$coast_region == marker_row$coast_region,
    ]
    lower_threshold <- find_spei_threshold(
        marker_group,
        "predicted_lower_pct",
        marker_row$threshold_pct,
        marker_row$anomaly_side,
        marker_row$impact_direction
    )
    upper_threshold <- find_spei_threshold(
        marker_group,
        "predicted_upper_pct",
        marker_row$threshold_pct,
        marker_row$anomaly_side,
        marker_row$impact_direction
    )
    threshold_range <- c(marker_row$SPEI_threshold, lower_threshold, upper_threshold)
    threshold_range <- threshold_range[!is.na(threshold_range)]
    coast_impact_threshold_markers$SPEI_threshold_lower[i] <- min(threshold_range)
    coast_impact_threshold_markers$SPEI_threshold_upper[i] <- max(threshold_range)
}}
coast_impact_threshold_markers$threshold_pct <- factor(
    coast_impact_threshold_markers$threshold_pct,
    levels = c(5, 10, 20),
    labels = c("5%", "10%", "20%")
)

coast_threshold_summary_table <- coast_impact_threshold_markers |>
    arrange(SPEI_timescale, coast_region, anomaly_side, impact_direction, threshold_pct)

write_table(coast_impact_threshold_markers, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv")
write_table(coast_threshold_summary_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv")

# GAM threshold markers (identical to Malone)
gam_impact_threshold_markers <- do.call(
    rbind,
    lapply(
        split(prediction_impact, list(prediction_impact$SPEI_timescale, prediction_impact$water_class), drop = TRUE),
        function(group_data) {{
            do.call(
                rbind,
                lapply(c(5, 10, 20), function(threshold_pct) {{
                    dry_decrease <- group_data$SPEI_value[
                        group_data$SPEI_value < -1 & group_data$predicted_pct_change <= -threshold_pct
                    ]
                    dry_increase <- group_data$SPEI_value[
                        group_data$SPEI_value < -1 & group_data$predicted_pct_change >= threshold_pct
                    ]
                    wet_decrease <- group_data$SPEI_value[
                        group_data$SPEI_value > 1 & group_data$predicted_pct_change <= -threshold_pct
                    ]
                    wet_increase <- group_data$SPEI_value[
                        group_data$SPEI_value > 1 & group_data$predicted_pct_change >= threshold_pct
                    ]
                    data.frame(
                        SPEI_timescale = group_data$SPEI_timescale[1],
                        water_class = group_data$water_class[1],
                        threshold_pct = threshold_pct,
                        impact_direction = c("decrease", "increase", "decrease", "increase"),
                        anomaly_side = c("dry", "dry", "wet", "wet"),
                        SPEI_threshold = c(
                            ifelse(length(dry_decrease) > 0, max(dry_decrease, na.rm = TRUE), NA_real_),
                            ifelse(length(dry_increase) > 0, max(dry_increase, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_decrease) > 0, min(wet_decrease, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_increase) > 0, min(wet_increase, na.rm = TRUE), NA_real_)
                        ),
                        pct_change_threshold = c(
                            -threshold_pct,
                            threshold_pct,
                            -threshold_pct,
                            threshold_pct
                        )
                    )
                }})
            )
        }}
    )
)
gam_impact_threshold_markers <- gam_impact_threshold_markers[
    !is.na(gam_impact_threshold_markers$SPEI_threshold),
]
gam_impact_threshold_markers$SPEI_threshold_lower <- NA_real_
gam_impact_threshold_markers$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(gam_impact_threshold_markers))) {{
    marker_row <- gam_impact_threshold_markers[i, ]
    marker_group <- prediction_impact[
        prediction_impact$SPEI_timescale == marker_row$SPEI_timescale &
            prediction_impact$water_class == marker_row$water_class,
    ]
    lower_threshold <- find_spei_threshold(
        marker_group,
        "predicted_lower_pct",
        marker_row$threshold_pct,
        marker_row$anomaly_side,
        marker_row$impact_direction
    )
    upper_threshold <- find_spei_threshold(
        marker_group,
        "predicted_upper_pct",
        marker_row$threshold_pct,
        marker_row$anomaly_side,
        marker_row$impact_direction
    )
    threshold_range <- c(marker_row$SPEI_threshold, lower_threshold, upper_threshold)
    threshold_range <- threshold_range[!is.na(threshold_range)]
    gam_impact_threshold_markers$SPEI_threshold_lower[i] <- min(threshold_range)
    gam_impact_threshold_markers$SPEI_threshold_upper[i] <- max(threshold_range)
}}
gam_impact_threshold_markers$threshold_pct <- factor(
    gam_impact_threshold_markers$threshold_pct,
    levels = c(5, 10, 20),
    labels = c("5%", "10%", "20%")
)

gam_threshold_summary_table <- gam_impact_threshold_markers |>
    arrange(SPEI_timescale, water_class, anomaly_side, impact_direction, threshold_pct)

write_table(gam_impact_threshold_markers, "Q3_WUE_T_SPEI_GAM_impact_threshold_markers_5_10_20pct.csv")
write_table(gam_threshold_summary_table, "Q3_WUE_T_SPEI_GAM_threshold_summary_5_10_20pct.csv")

# ---------------------------------------------------------------------------
# Model summary (identical to Malone's capture_to_file)
# ---------------------------------------------------------------------------

cat("\\n")
sep_line()
cat("MODEL SUMMARY - EXACT MALONE REPLICATION")
cat("\\n")
sep_line()
cat("\\n")

sink(file.path("{output_dir_r}", "Q3_WUE_T_SPEI_model_summary.txt"))

cat("Q3 WUE_T sensitivity to SPEI gradients\\n")
cat("Input:", "{input_file_r}", "\\n")
cat("Output:", "{output_dir_r}", "\\n\\n")
cat("Rows in source monthly data:", {n_monthly}, "\\n")
cat("Rows after complete-case and ecosystem filters:", {n_model_base}, "\\n")
cat("Rows in SPEI long data:", {n_spei_long}, "\\n")
cat("Sites:", {n_sites}, "\\n")
cat("Years:", {n_years_min}, "-", {n_years_max}, "\\n\\n")

cat("Primary model formula:\\n")
print(formula(smooth_model))

cat("\\nGAM model information:\\n")
print(model_comparison)

cat("\\nSmooth GAM summary:\\n")
print(summary(smooth_model))

cat("\\nCoast-threshold GAM summary:\\n")
print(summary(smooth_coast_model))

cat("\\nGAM-predicted ecological-impact threshold markers:\\n")
print(gam_threshold_summary_table)

cat("\\nCoast-specific GAM-predicted upland ecological-impact threshold markers:\\n")
print(coast_threshold_summary_table)

sink()

cat("\\n")
sep_line()
cat("R MODELS COMPLETE - EXACT MALONE REPLICATION")
cat("\\n")
sep_line()
cat("\\n")
'''

# Save R script in temp folder (FIXED)
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

if missing_files:
    raise FileNotFoundError(
        "Missing expected Q3 output files:\n" + "\n".join(missing_files)
    )
else:
    print(f"All {len(expected_outputs)} expected outputs created successfully!")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("CHUNK 1 COMPLETE!")
print("="*60)
print(f"\nOutputs saved to: {output_dir}")
print(f"Temp files saved to: {temp_dir}")
print("\nRow counts used in summary:")
print(f"  Source monthly data: {n_monthly}")
print(f"  After filtering: {n_model_base}")
print(f"  SPEI long data: {n_spei_long}")
print(f"  Sites: {n_sites}")
print(f"  Years: {n_years_min} - {n_years_max}")
print("\nFIXES APPLIED:")
print("  - model_data_long keeps ALL columns (like Malone)")
print("  - site_response_data uses LEFT JOIN then filter (like Malone)")
print("  - Temp files in temp_dir (not output_dir)")
print("  - Hard failure for missing outputs")
print("\nNow run CHUNK 2 to create figures.")
print("="*60)