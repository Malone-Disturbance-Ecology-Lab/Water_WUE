"""
Q3_COMPLETE_UPDATED_WORKFLOW.py
Complete Q3 workflow with confirmed values and refined interpretation
Based on helper script outputs and corrected data extraction
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime

print("="*80)
print("Q3 COMPLETE UPDATED WORKFLOW - METHODS & RESULTS EXTRACTION")
print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# ============================================================================
# PATHS
# ============================================================================

base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures")

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")

# ============================================================================
# CONSTANTS - CONFIRMED VALUES FROM ANALYSIS
# ============================================================================

# Confirmed dataset stats
N_OBSERVATIONS = 1872
N_SITES = 64
YEAR_START = 1994
YEAR_END = 2025
"""
CHUNK 3A: Flux-Tower Disturbance Exposure Screening
Purpose: For each flux tower, identify IDS disturbance polygons that contain the tower point 
or overlap simple tower-centered windows.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import box
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

DISTURBANCE_GPKG = r"M:\Research\NE_temperate_forest_resilience\Data\Disturbance_data\results\NE_AOI_all_disturbances_master.gpkg"
TOWER_CSV = r"M:\Research\NE_temperate_forest_resilience\Data\tower_sites.csv"
OUTPUT_FOLDER = r"M:\Research\NE_temperate_forest_resilience\Data\Disturbance_data\results\tower_screening"

# Create output folder if it doesn't exist
Path(OUTPUT_FOLDER).mkdir(parents=True, exist_ok=True)

# Window sizes in meters (radial distance from tower)
WINDOW_SIZES = [500, 1500, 2500, 4000]

# Use EPSG:5070 (NAD83 / Conus Albers) - good equal-area projected CRS for CONUS in meters
TARGET_CRS = "EPSG:5070"


# ============================================================================
# Step 1: Load data
# ============================================================================

print("=" * 80)
print("Step 1: Loading data")
print("=" * 80)

# Load disturbance polygons
print("Loading disturbance polygons...")
disturbance_gdf = gpd.read_file(DISTURBANCE_GPKG)
print(f"  Loaded {len(disturbance_gdf)} disturbance polygons")
print(f"  Original CRS: {disturbance_gdf.crs}")

# Force projected CRS in meters
print(f"  Reprojecting to {TARGET_CRS}...")
disturbance_gdf = disturbance_gdf.to_crs(TARGET_CRS)
print(f"  New CRS: {disturbance_gdf.crs}")

# Load tower locations
print("Loading tower locations...")
tower_df = pd.read_csv(TOWER_CSV)
print(f"  Loaded {len(tower_df)} tower records")
print(f"  Columns: {list(tower_df.columns)}")

# Detect tower ID column
id_col = None
for col in ['Site ID', 'Site Name', 'site_code', 'tower_id', 'SITE_ID', 'SITE_NAME']:
    if col in tower_df.columns:
        id_col = col
        break
if id_col is None:
    id_col = tower_df.columns[0]  # fallback to first column
print(f"  Using ID column: {id_col}")

# Detect latitude column
lat_col = None
for col in ['Latitude', 'latitude', 'lat', 'LATITUDE', 'LAT']:
    if col in tower_df.columns:
        lat_col = col
        break
if lat_col is None:
    lat_col = tower_df.columns[1]  # fallback
print(f"  Using latitude column: {lat_col}")

# Detect longitude column
lon_col = None
for col in ['Longitude', 'longitude', 'lon', 'LONGITUDE', 'LON']:
    if col in tower_df.columns:
        lon_col = col
        break
if lon_col is None:
    lon_col = tower_df.columns[0]  # fallback
print(f"  Using longitude column: {lon_col}")

# Create tower GeoDataFrame
print("Creating tower GeoDataFrame...")
tower_gdf = gpd.GeoDataFrame(
    tower_df,
    geometry=gpd.points_from_xy(tower_df[lon_col], tower_df[lat_col]),
    crs="EPSG:4326"
)

# Reproject to match disturbance CRS
tower_gdf = tower_gdf.to_crs(TARGET_CRS)
print(f"  Reprojected towers to {TARGET_CRS}")


# ============================================================================
# Step 2: Create tower windows
# ============================================================================

print("\n" + "=" * 80)
print("Step 2: Creating tower windows (circular buffers)")
print("=" * 80)

def create_window(point, radius_m):
    """Create circular buffer around tower point with given radius in meters."""
    return point.buffer(radius_m)

# Create windows for each tower and window size
tower_windows = []
for idx, tower in tower_gdf.iterrows():
    tower_id = tower[id_col]
    tower_name = tower.get('tower_name', tower_id)
    
    for size in WINDOW_SIZES:
        window_geom = create_window(tower.geometry, size)
        window_area_m2 = window_geom.area
        window_area_km2 = window_area_m2 / 1_000_000
        window_area_acres = window_area_m2 * 0.000247105
        
        tower_windows.append({
            'tower_id': tower_id,
            'tower_name': tower_name,
            'latitude': tower[lat_col],
            'longitude': tower[lon_col],
            'window_size_m': size,
            'window_geometry': window_geom,
            'tower_point_geometry': tower.geometry,  # Store real tower point for later use
            'window_area_m2': window_area_m2,
            'window_area_km2': window_area_km2,
            'window_area_acres': window_area_acres
        })

# Create GeoDataFrame with proper geometry column
tower_windows_df = gpd.GeoDataFrame(
    tower_windows,
    geometry="window_geometry",
    crs=TARGET_CRS
)
print(f"  Created {len(tower_windows_df)} tower-window combinations")


# ============================================================================
# Step 3: Find disturbance polygons near each tower
# ============================================================================

print("\n" + "=" * 80)
print("Step 3: Finding disturbance polygon overlaps")
print("=" * 80)

# Build spatial index for disturbance polygons
print("Building spatial index...")
disturbance_sindex = disturbance_gdf.sindex

all_overlaps = []

for tower_idx, tower_window in tower_windows_df.iterrows():
    tower_id = tower_window['tower_id']
    tower_name = tower_window['tower_name']
    window_size = tower_window['window_size_m']
    window_geom = tower_window['window_geometry']
    window_area_m2 = tower_window['window_area_m2']
    window_area_km2 = tower_window['window_area_km2']
    window_area_acres = tower_window['window_area_acres']
    tower_point = tower_window['tower_point_geometry']
    
    # Find candidate polygons using spatial index
    possible_matches_idx = list(disturbance_sindex.intersection(window_geom.bounds))
    
    if not possible_matches_idx:
        continue
    
    # Get candidate polygons
    candidates = disturbance_gdf.iloc[possible_matches_idx]
    
    # Check for actual intersections
    for dist_idx, disturbance in candidates.iterrows():
        if not window_geom.intersects(disturbance.geometry):
            continue
        
        # Calculate overlap
        overlap = window_geom.intersection(disturbance.geometry)
        if overlap.is_empty:
            continue
        
        overlap_area_m2 = overlap.area
        overlap_area_km2 = overlap_area_m2 / 1_000_000
        overlap_area_acres = overlap_area_m2 * 0.000247105
        percent_of_window_overlapped = (overlap_area_m2 / window_area_m2) * 100
        
        # Calculate percent of disturbance polygon overlapped
        disturbance_area_m2 = disturbance.geometry.area
        percent_of_polygon_overlapped = (
            overlap_area_m2 / disturbance_area_m2 * 100
            if disturbance_area_m2 > 0 else np.nan
        )
        
        # Check if tower point is inside disturbance polygon using .covers()
        tower_inside = disturbance.geometry.covers(tower_point)
        
        # Get disturbance attributes
        all_overlaps.append({
            'tower_id': tower_id,
            'tower_name': tower_name,
            'latitude': tower_window['latitude'],
            'longitude': tower_window['longitude'],
            'window_size_m': window_size,
            'window_area_km2': window_area_km2,
            'SURVEY_YEAR': disturbance.get('SURVEY_YEAR', None),
            'original_agent': disturbance.get('original_agent', None),
            'original_damage_type': disturbance.get('original_damage_type', None),
            'original_host': disturbance.get('original_host', None),
            'original_host_group': disturbance.get('original_host_group', None),
            'DAMAGE_AREA_ID': disturbance.get('DAMAGE_AREA_ID', None),
            'OBSERVATION_ID': disturbance.get('OBSERVATION_ID', None),
            'AOI_state': disturbance.get('AOI_state', None),
            'AOI_state_abbr': disturbance.get('AOI_state_abbr', None),
            'AREA_TYPE': disturbance.get('AREA_TYPE', None),
            'DATA_SOURCE_NAME': disturbance.get('DATA_SOURCE_NAME', None),
            'PERCENT_AFFECTED': disturbance.get('PERCENT_AFFECTED', None),
            'PERCENT_MID': disturbance.get('PERCENT_MID', None),
            'LEGACY_SEVERITY': disturbance.get('LEGACY_SEVERITY', None),
            'tower_inside_disturbance_polygon': tower_inside,
            'overlap_area_m2': overlap_area_m2,
            'overlap_area_km2': overlap_area_km2,
            'overlap_area_acres': overlap_area_acres,
            'percent_of_window_overlapped': percent_of_window_overlapped,
            'disturbance_polygon_area_km2': disturbance_area_m2 / 1_000_000,
            'percent_of_polygon_overlapped_by_window': percent_of_polygon_overlapped
        })

overlaps_df = pd.DataFrame(all_overlaps)
print(f"  Found {len(overlaps_df)} tower-window-disturbance overlaps")


# ============================================================================
# Step 4: Create simple damage-effect flags
# ============================================================================

print("\n" + "=" * 80)
print("Step 4: Creating damage-effect flags")
print("=" * 80)

def create_damage_flags(damage_type):
    """Create damage effect flags based on original_damage_type."""
    if pd.isna(damage_type):
        damage_type = ''
    else:
        damage_type = str(damage_type).lower()
    
    return {
        'mortality_flag': any(x in damage_type for x in ['mortality', 'dead', 'dieback', 'topkill']),
        'defoliation_flag': any(x in damage_type for x in ['defoliation', 'defol']),
        'discoloration_flag': 'discoloration' in damage_type,
        'crown_dieback_flag': 'crown dieback' in damage_type,
        'breakage_or_uprooted_flag': any(x in damage_type for x in ['branch breakage', 'main stem broken', 'uprooted'])
    }

if len(overlaps_df) > 0:
    flags = overlaps_df['original_damage_type'].apply(create_damage_flags)
    flags_df = pd.DataFrame(flags.tolist())
    overlaps_df = pd.concat([overlaps_df, flags_df], axis=1)
else:
    # Add empty flag columns if no overlaps
    for col in ['mortality_flag', 'defoliation_flag', 'discoloration_flag', 
                'crown_dieback_flag', 'breakage_or_uprooted_flag']:
        overlaps_df[col] = False

print(f"  Added damage flags to {len(overlaps_df)} records")


# ============================================================================
# Step 5: Create reported damage intensity field (FIXED - case insensitive)
# ============================================================================

print("\n" + "=" * 80)
print("Step 5: Creating damage intensity field")
print("=" * 80)

def parse_percent_affected(value):
    """
    Parse PERCENT_AFFECTED text to numeric value.
    Case insensitive - handles both uppercase and lowercase text.
    """
    if pd.isna(value):
        return None
    value = str(value).strip().lower()
    
    if 'very severe' in value or '>50' in value:
        return 75
    elif 'severe' in value or '30-50' in value:
        return 40
    elif 'moderate' in value or '11-29' in value:
        return 20
    elif 'very light' in value or '1-3' in value:
        return 2
    elif 'light' in value or '4-10' in value:
        return 7
    else:
        return None

def calculate_damage_intensity(row):
    """Calculate damage intensity using PERCENT_MID, then PERCENT_AFFECTED, then LEGACY_SEVERITY."""
    # Try PERCENT_MID first
    percent_mid = row.get('PERCENT_MID')
    if pd.notna(percent_mid):
        try:
            return float(percent_mid)
        except (ValueError, TypeError):
            pass
    
    # Try PERCENT_AFFECTED
    percent_affected = row.get('PERCENT_AFFECTED')
    if pd.notna(percent_affected):
        parsed = parse_percent_affected(percent_affected)
        if parsed is not None:
            return parsed
    
    # Try LEGACY_SEVERITY
    legacy = row.get('LEGACY_SEVERITY')
    if pd.notna(legacy):
        parsed = parse_percent_affected(legacy)
        if parsed is not None:
            return parsed
    
    return None

if len(overlaps_df) > 0:
    overlaps_df['damage_intensity_mid'] = overlaps_df.apply(calculate_damage_intensity, axis=1)
    overlaps_df['high_damage_50_flag'] = overlaps_df['damage_intensity_mid'] >= 50
else:
    overlaps_df['damage_intensity_mid'] = None
    overlaps_df['high_damage_50_flag'] = False

print(f"  Calculated damage intensity for {len(overlaps_df)} records")


# ============================================================================
# Step 6: Save polygon-level tower overlap table
# ============================================================================

print("\n" + "=" * 80)
print("Step 6: Saving polygon-level overlap table")
print("=" * 80)

# Ensure all required columns exist (add any missing ones with NaN)
required_columns = [
    'tower_id', 'tower_name', 'latitude', 'longitude', 'window_size_m', 'window_area_km2',
    'SURVEY_YEAR', 'original_agent', 'original_damage_type', 'original_host', 
    'original_host_group', 'DAMAGE_AREA_ID', 'OBSERVATION_ID', 'AOI_state', 
    'AOI_state_abbr', 'AREA_TYPE', 'DATA_SOURCE_NAME', 'PERCENT_AFFECTED', 
    'PERCENT_MID', 'LEGACY_SEVERITY', 'damage_intensity_mid', 'high_damage_50_flag',
    'mortality_flag', 'defoliation_flag', 'discoloration_flag', 'crown_dieback_flag',
    'breakage_or_uprooted_flag', 'tower_inside_disturbance_polygon', 'overlap_area_m2',
    'overlap_area_km2', 'overlap_area_acres', 'percent_of_window_overlapped',
    'disturbance_polygon_area_km2', 'percent_of_polygon_overlapped_by_window'
]

# Add any missing columns
for col in required_columns:
    if col not in overlaps_df.columns:
        overlaps_df[col] = None

# Save
output_file = Path(OUTPUT_FOLDER) / "NE_AOI_tower_disturbance_polygon_overlaps.csv"
overlaps_df[required_columns].to_csv(output_file, index=False)
print(f"  Saved to: {output_file}")
print(f"  Records: {len(overlaps_df)}")


# ============================================================================
# Step 7: Save tower-year-agent summary table (FIXED)
# ============================================================================

print("\n" + "=" * 80)
print("Step 7: Saving tower-year-agent summary table")
print("=" * 80)

if len(overlaps_df) > 0:
    # Create damage-specific area columns BEFORE grouping
    overlaps_df["mortality_overlap_area_km2"] = np.where(
        overlaps_df["mortality_flag"], overlaps_df["overlap_area_km2"], 0
    )
    overlaps_df["defoliation_overlap_area_km2"] = np.where(
        overlaps_df["defoliation_flag"], overlaps_df["overlap_area_km2"], 0
    )
    overlaps_df["discoloration_overlap_area_km2"] = np.where(
        overlaps_df["discoloration_flag"], overlaps_df["overlap_area_km2"], 0
    )
    overlaps_df["crown_dieback_overlap_area_km2"] = np.where(
        overlaps_df["crown_dieback_flag"], overlaps_df["overlap_area_km2"], 0
    )

    # Define grouping columns
    summary_cols = [
        "tower_id", "tower_name", "window_size_m",
        "SURVEY_YEAR", "original_agent", "original_damage_type"
    ]

    # Aggregate using named aggregation (fixes the duplicate column bug)
    summary_df = (
        overlaps_df
        .groupby(summary_cols, dropna=False)
        .agg(
            overlap_record_count=("DAMAGE_AREA_ID", "count"),
            unique_damage_area_count=("DAMAGE_AREA_ID", "nunique"),
            unique_observation_count=("OBSERVATION_ID", "nunique"),
            total_overlap_area_km2=("overlap_area_km2", "sum"),
            total_overlap_area_acres=("overlap_area_acres", "sum"),
            percent_of_window_overlapped_total=("percent_of_window_overlapped", "sum"),
            mean_damage_intensity_mid=("damage_intensity_mid", "mean"),
            max_damage_intensity_mid=("damage_intensity_mid", "max"),
            high_damage_50_record_count=("high_damage_50_flag", "sum"),
            tower_inside_any_polygon=("tower_inside_disturbance_polygon", "any"),
            mortality_overlap_area_km2=("mortality_overlap_area_km2", "sum"),
            defoliation_overlap_area_km2=("defoliation_overlap_area_km2", "sum"),
            discoloration_overlap_area_km2=("discoloration_overlap_area_km2", "sum"),
            crown_dieback_overlap_area_km2=("crown_dieback_overlap_area_km2", "sum"),
            dominant_host=("original_host", lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else None),
            states_present=("AOI_state_abbr", lambda x: ", ".join(sorted(set(x.dropna()))))
        )
        .reset_index()
    )

    # Save
    output_file = Path(OUTPUT_FOLDER) / "NE_AOI_tower_disturbance_year_agent_summary.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"  Saved to: {output_file}")
    print(f"  Records: {len(summary_df)}")
else:
    print("  No overlaps found, skipping summary table")


# ============================================================================
# Step 8: Save top 10 tower-specific disturbances
# ============================================================================

print("\n" + "=" * 80)
print("Step 8: Saving top 10 disturbances by window")
print("=" * 80)

if len(overlaps_df) > 0 and len(summary_df) > 0:
    # Group by tower, window
    top10_list = []
    
    for (tower_id, tower_name, window_size), group in summary_df.groupby(['tower_id', 'tower_name', 'window_size_m']):
        # Sort by total overlap area and take top 10
        top10 = group.nlargest(10, 'total_overlap_area_km2')
        top10 = top10.copy()
        top10['rank'] = range(1, len(top10) + 1)
        top10_list.append(top10)
    
    if top10_list:
        top10_df = pd.concat(top10_list, ignore_index=True)
        
        # Select columns
        top10_cols = [
            'tower_id', 'tower_name', 'window_size_m', 'rank',
            'SURVEY_YEAR', 'original_agent', 'original_damage_type',
            'total_overlap_area_km2', 'percent_of_window_overlapped_total',
            'mean_damage_intensity_mid', 'high_damage_50_record_count',
            'tower_inside_any_polygon',
            'mortality_overlap_area_km2', 'defoliation_overlap_area_km2',
            'discoloration_overlap_area_km2', 'crown_dieback_overlap_area_km2',
            'dominant_host', 'states_present'
        ]
        top10_df = top10_df[top10_cols]
        
        # Save
        output_file = Path(OUTPUT_FOLDER) / "NE_AOI_tower_top10_disturbances_by_window.csv"
        top10_df.to_csv(output_file, index=False)
        print(f"  Saved to: {output_file}")
        print(f"  Records: {len(top10_df)}")
    else:
        print("  No top 10 records to save")
else:
    print("  No overlaps found, skipping top 10 table")


# ============================================================================
# Step 9: Save simple tower screening recommendation table (FIXED - safer)
# ============================================================================

print("\n" + "=" * 80)
print("Step 9: Saving tower screening recommendations")
print("=" * 80)

recommendations = []

for tower_id in tower_gdf[id_col].unique():
    tower_name = tower_gdf[tower_gdf[id_col] == tower_id]['tower_name'].iloc[0] if 'tower_name' in tower_gdf.columns else tower_id
    
    rec = {
        'tower_id': tower_id,
        'tower_name': tower_name,
        'has_disturbance_500m': False,
        'has_disturbance_1500m': False,
        'has_disturbance_2500m': False,
        'has_disturbance_4000m': False,
        'years_with_disturbance_500m': '',
        'years_with_disturbance_1500m': '',
        'years_with_disturbance_2500m': '',
        'years_with_disturbance_4000m': '',
        'top_disturbance_500m': '',
        'top_disturbance_1500m': '',
        'top_disturbance_2500m': '',
        'top_disturbance_4000m': ''
    }
    
    # Check each window size
    for size in WINDOW_SIZES:
        size_col = f'has_disturbance_{size}m'
        years_col = f'years_with_disturbance_{size}m'
        top_col = f'top_disturbance_{size}m'
        
        # Get overlaps for this tower and window
        tower_overlaps = overlaps_df[(overlaps_df['tower_id'] == tower_id) & 
                                     (overlaps_df['window_size_m'] == size)]
        
        if len(tower_overlaps) > 0:
            rec[size_col] = True
            
            # Get unique years
            years = tower_overlaps['SURVEY_YEAR'].dropna().unique()
            rec[years_col] = ', '.join(sorted([str(int(y)) for y in years if pd.notna(y)]))
            
            # Get top disturbance (by overlap area) - SAFER VERSION with missing value handling
            top_group = (
                tower_overlaps
                .groupby(['original_agent', 'original_damage_type', 'SURVEY_YEAR'], dropna=False)['overlap_area_km2']
                .sum()
            )
            
            if len(top_group) > 0:
                top_idx = top_group.idxmax()
                if isinstance(top_idx, tuple):
                    # Handle potential None/NaN values in the tuple
                    agent = top_idx[0] if top_idx[0] is not None else 'Unknown Agent'
                    damage_type = top_idx[1] if top_idx[1] is not None else 'Unknown Damage'
                    year = top_idx[2] if top_idx[2] is not None else 'Unknown Year'
                    rec[top_col] = f"{agent} - {damage_type} ({year})"
    
    # Determine recommended first window
    if rec['has_disturbance_500m']:
        rec['recommended_first_window_to_test'] = '500m'
    elif rec['has_disturbance_1500m']:
        rec['recommended_first_window_to_test'] = '1500m'
    elif rec['has_disturbance_2500m']:
        rec['recommended_first_window_to_test'] = '2500m'
    elif rec['has_disturbance_4000m']:
        rec['recommended_first_window_to_test'] = '4000m'
    else:
        rec['recommended_first_window_to_test'] = 'no nearby mapped disturbance'
    
    recommendations.append(rec)

recommendations_df = pd.DataFrame(recommendations)

# Save
output_file = Path(OUTPUT_FOLDER) / "NE_AOI_tower_screening_recommendations.csv"
recommendations_df.to_csv(output_file, index=False)
print(f"  Saved to: {output_file}")
print(f"  Records: {len(recommendations_df)}")


# ============================================================================
# Step 10: Console summary
# ============================================================================

print("\n" + "=" * 80)
print("Step 10: Console Summary")
print("=" * 80)

print(f"\nNumber of towers loaded: {len(tower_gdf)}")
print(f"Number of disturbance polygons loaded: {len(disturbance_gdf)}")
print(f"Number of tower-window polygon overlaps: {len(overlaps_df)}")

print("\n" + "-" * 80)
print("Tower-specific summaries:")
print("-" * 80)

for tower_id in tower_gdf[id_col].unique():
    tower_name = tower_gdf[tower_gdf[id_col] == tower_id]['tower_name'].iloc[0] if 'tower_name' in tower_gdf.columns else tower_id
    
    print(f"\nTower: {tower_id} ({tower_name})")
    
    # Check disturbances in each window
    for size in WINDOW_SIZES:
        count = len(overlaps_df[(overlaps_df['tower_id'] == tower_id) & 
                               (overlaps_df['window_size_m'] == size)])
        print(f"  {size}m window: {'Yes' if count > 0 else 'No'} ({count} overlaps)")
    
    # Get top 5 disturbances for 1500m and 4000m
    for size in [1500, 4000]:
        tower_overlaps = overlaps_df[(overlaps_df['tower_id'] == tower_id) & 
                                     (overlaps_df['window_size_m'] == size)]
        
        if len(tower_overlaps) > 0:
            top5 = tower_overlaps.groupby(['original_agent', 'original_damage_type', 'SURVEY_YEAR'])['overlap_area_km2'].sum().nlargest(5)
            print(f"\n  Top 5 for {size}m window:")
            for idx, area in top5.items():
                print(f"    - {idx[0]} / {idx[1]} ({idx[2]}) - {area:.2f} km²")
    
    # Get recommendation
    rec = recommendations_df[recommendations_df['tower_id'] == tower_id]
    if len(rec) > 0:
        print(f"\n  Recommended first window: {rec['recommended_first_window_to_test'].iloc[0]}")

print("\n" + "=" * 80)
print("Processing complete!")
print("=" * 80)N_YEARS = 32
N_LONG_OBS = 13104
N_SPEI_TIMESCALES = 7

# Confirmed ecosystem distribution
ECOSYSTEM_DIST = {
    "Upland": {"obs": 793, "sites": 29},
    "Freshwater": {"obs": 629, "sites": 20},
    "Saline": {"obs": 450, "sites": 15}
}

# Confirmed coast distribution
COAST_DIST = {
    "Atlantic Coast": {"obs": 1017, "sites": 29},
    "Pacific Coast": {"obs": 560, "sites": 18},
    "Gulf Coast": {"obs": 110, "sites": 5},
    "AK Coast": {"obs": 185, "sites": 12}
}

# Confirmed model comparison values
MODEL_COMP = {
    "ecosystem": {
        "AIC": 55467.76,
        "BIC": 56263.82,
        "logLik": -27627.46,
        "deviance_explained": 0.579,
        "r_squared": 0.575
    },
    "coast": {
        "AIC": 55347.65,
        "BIC": 56361.12,
        "logLik": -27538.35,
        "deviance_explained": 0.584,
        "r_squared": 0.580
    }
}

# Confirmed strongest signal
STRONGEST_SIGNAL = {
    "coast": "Gulf Coast",
    "timescale": "SPEI_3",
    "edf": 4.34,
    "f_stat": 29.44,
    "p_value": "<0.001"
}

# Confirmed coastline sensitivity medians
COAST_SENSITIVITY = {
    "Alaska": {"SPEI_1": 37.5, "SPEI_3": 43.6, "SPEI_48": 44.7},
    "Atlantic": {"SPEI_1": 18.3, "SPEI_3": 18.0, "SPEI_48": 20.2},
    "Gulf": {"SPEI_1": 39.3, "SPEI_3": 42.5, "SPEI_48": 37.8},
    "Pacific": {"SPEI_1": 17.2, "SPEI_3": 17.7, "SPEI_48": 16.3}
}

# Confirmed threshold ranges
THRESHOLDS = {
    "Gulf_SPEI3_dry_decrease": {"min": -2.13, "max": -2.02},
    "Gulf_SPEI3_wet_increase": {"min": 1.37, "max": 1.81},
    "AK_SPEI48_dry_decrease": {"min": -1.44, "max": -1.44},
    "AK_SPEI48_wet_increase": {"min": 1.44, "max": 1.44},
    "Atlantic_SPEI48_wet_decrease": {"min": 1.20, "max": 1.20},
    "Pacific_SPEI1_dry_increase": {"min": -1.78, "max": -1.03}
}

# ============================================================================
# SECTION 1: METHODS - DATASET DESCRIPTION
# ============================================================================

print("\n" + "="*80)
print("SECTION 1: METHODS - DATASET DESCRIPTION (CONFIRMED VALUES)")
print("="*80)

print(f"\nFinal dataset after Q3 complete-case filtering:")
print(f"  • {N_OBSERVATIONS:,} site-month observations")
print(f"  • {N_SITES} sites")
print(f"  • {YEAR_START} to {YEAR_END} ({N_YEARS} years)")

print(f"\nLong-format dataset (reshaped across SPEI timescales):")
print(f"  • {N_LONG_OBS:,} site-month-timescale observations")
print(f"  • {N_SPEI_TIMESCALES} SPEI timescales (1-48 months)")

print(f"\nObservations by ecosystem class:")
for eco, stats in ECOSYSTEM_DIST.items():
    print(f"  • {eco}: {stats['obs']:,} observations ({stats['sites']} sites)")

print(f"\nObservations by coast region:")
for coast, stats in COAST_DIST.items():
    print(f"  • {coast}: {stats['obs']:,} observations ({stats['sites']} sites)")

# ============================================================================
# SECTION 2: METHODS - MODEL COMPARISON
# ============================================================================

print("\n" + "="*80)
print("SECTION 2: METHODS - MODEL COMPARISON (TABLE Q3.1)")
print("="*80)

print("\nTable Q3.1: Comparison of ecosystem and coast-threshold GAMs")
print("-" * 80)
print(f"{'Model':<20} {'AIC':>10} {'BIC':>10} {'logLik':>10} {'Deviance Expl.':>14} {'Adj R²':>10}")
print("-" * 80)

eco = MODEL_COMP["ecosystem"]
coast = MODEL_COMP["coast"]

print(f"{'Ecosystem GAM':<20} {eco['AIC']:>10.2f} {eco['BIC']:>10.2f} {eco['logLik']:>10.2f} {eco['deviance_explained']*100:>13.1f}% {eco['r_squared']:>10.3f}")
print(f"{'Coast-threshold GAM':<20} {coast['AIC']:>10.2f} {coast['BIC']:>10.2f} {coast['logLik']:>10.2f} {coast['deviance_explained']*100:>13.1f}% {coast['r_squared']:>10.3f}")

print("\nModel Selection Rationale:")
print(f"  • Lower AIC: {coast['AIC']:.2f} vs {eco['AIC']:.2f} (ΔAIC = {eco['AIC'] - coast['AIC']:.2f})")
print(f"  • Higher deviance explained: {coast['deviance_explained']*100:.1f}% vs {eco['deviance_explained']*100:.1f}%")
print(f"  • Higher adjusted R²: {coast['r_squared']:.3f} vs {eco['r_squared']:.3f}")

print("\nMETHODS SENTENCE TEMPLATE:")
print(f"  \"The coast-threshold GAM performed better than the ecosystem GAM (AIC: {coast['AIC']:.2f} vs {eco['AIC']:.2f}, deviance explained: {coast['deviance_explained']*100:.1f}% vs {eco['deviance_explained']*100:.1f}%, adjusted R²: {coast['r_squared']:.3f} vs {eco['r_squared']:.3f}) and was used as the primary model.\"")

# ============================================================================
# SECTION 3: RESULTS - STRONGEST SIGNAL
# ============================================================================

print("\n" + "="*80)
print("SECTION 3: RESULTS - STRONGEST SIGNAL")
print("="*80)

strong = STRONGEST_SIGNAL
print(f"\nStrongest nonlinear SPEI signal:")
print(f"  • Coast: {strong['coast']}")
print(f"  • Timescale: {strong['timescale']}")
print(f"  • EDF (Effective Degrees of Freedom): {strong['edf']:.2f}")
print(f"  • F-statistic: {strong['f_stat']:.2f}")
print(f"  • p-value: {strong['p_value']}")

print("\nRESULTS SENTENCE TEMPLATE:")
print(f"  \"The strongest nonlinear SPEI signal occurred in the {strong['coast']} at {strong['timescale']} (EDF={strong['edf']:.2f}, F={strong['f_stat']:.2f}, p<0.001).\"")

# ============================================================================
# SECTION 4: RESULTS - COASTLINE SENSITIVITY MEDIANS (TABLE Q3.3)
# ============================================================================

print("\n" + "="*80)
print("SECTION 4: RESULTS - COASTLINE SENSITIVITY MEDIANS (TABLE Q3.3)")
print("="*80)

print("\nTable Q3.3: Coastline sensitivity medians for selected SPEI timescales")
print("-" * 80)
print(f"{'Coast':<10} {'SPEI-1 (%)':>12} {'SPEI-3 (%)':>12} {'SPEI-48 (%)':>13}")
print("-" * 80)

for coast, values in COAST_SENSITIVITY.items():
    print(f"{coast:<10} {values['SPEI_1']:>12.1f} {values['SPEI_3']:>12.1f} {values['SPEI_48']:>13.1f}")

# Calculate averages
print("\nCoast patterns (average median across selected timescales):")
for coast, values in COAST_SENSITIVITY.items():
    avg = sum(values.values()) / 3
    print(f"  • {coast}: {avg:.1f}%")

# Highest values
print("\nHighest observed sensitivity:")
max_overall = max(max(v.values()) for v in COAST_SENSITIVITY.values())
for coast, values in COAST_SENSITIVITY.items():
    if max(values.values()) == max_overall:
        ts = max(values, key=values.get)
        print(f"  • Overall: {coast} at {ts} ({max_overall:.1f}%)")

print("\nRESULTS SENTENCE TEMPLATE:")
print(f"  \"Observed site sensitivity was highest in the Alaska Coast (SPEI-48: {COAST_SENSITIVITY['Alaska']['SPEI_48']:.1f}%, SPEI-3: {COAST_SENSITIVITY['Alaska']['SPEI_3']:.1f}%) and Gulf Coast (SPEI-3: {COAST_SENSITIVITY['Gulf']['SPEI_3']:.1f}%, SPEI-1: {COAST_SENSITIVITY['Gulf']['SPEI_1']:.1f}%), while Atlantic and Pacific coasts showed lower sensitivity (approximately 16-20%).\"")

# ============================================================================
# SECTION 5: RESULTS - THRESHOLD SUMMARIES
# ============================================================================

print("\n" + "="*80)
print("SECTION 5: RESULTS - IMPACT THRESHOLDS")
print("="*80)

print("\nGulf Coast SPEI-3 thresholds (strongest signal):")
print(f"  • Dry-side decreases: SPEI {THRESHOLDS['Gulf_SPEI3_dry_decrease']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_dry_decrease']['max']:.2f}")
print(f"  • Wet-side increases: SPEI {THRESHOLDS['Gulf_SPEI3_wet_increase']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_wet_increase']['max']:.2f}")

print("\nAK Coast SPEI-48 thresholds:")
print(f"  • Dry-side decreases: SPEI {THRESHOLDS['AK_SPEI48_dry_decrease']['min']:.2f}")
print(f"  • Wet-side increases: SPEI {THRESHOLDS['AK_SPEI48_wet_increase']['min']:.2f}")

print("\nAtlantic Coast SPEI-48 thresholds:")
print(f"  • Wet-side decreases: SPEI {THRESHOLDS['Atlantic_SPEI48_wet_decrease']['min']:.2f}")

print("\nPacific Coast SPEI-1 thresholds:")
print(f"  • Dry-side increases: SPEI {THRESHOLDS['Pacific_SPEI1_dry_increase']['min']:.2f} to {THRESHOLDS['Pacific_SPEI1_dry_increase']['max']:.2f}")

# ============================================================================
# SECTION 6: RESULTS - SITE-LEVEL SUMMARY
# ============================================================================

print("\n" + "="*80)
print("SECTION 6: RESULTS - SITE-LEVEL SUMMARY")
print("="*80)

# Load site slopes for additional stats
try:
    site_slopes = pd.read_csv(os.path.join(output_dir, "Q3_WUE_T_SPEI_site_level_slopes.csv"))
    
    print(f"\nSite-level slopes across all sites and timescales:")
    print(f"  • Mean slope: {site_slopes['SPEI_slope'].mean():.4f}")
    print(f"  • Median slope: {site_slopes['SPEI_slope'].median():.4f}")
    print(f"  • % negative slopes: {(site_slopes['SPEI_slope'] < 0).mean()*100:.1f}%")
    
    # Top sensitive sites
    sensitivity = pd.read_csv(os.path.join(output_dir, "Q3_WUE_T_SPEI_site_sensitivity_rank.csv"))
    top5 = sensitivity.nlargest(5, 'mean_abs_pct_change')[['site_name', 'water_class', 'coast_region', 'SPEI_timescale', 'mean_abs_pct_change']]
    
    print(f"\nTop 5 most sensitive sites overall:")
    for _, row in top5.iterrows():
        print(f"  • {row['site_name']} ({row['water_class']}, {row['coast_region']}) - {row['mean_abs_pct_change']:.1f}% at {row['SPEI_timescale']}")
    
except Exception as e:
    print(f"  Note: Could not load site-level data: {e}")

# ============================================================================
# SECTION 7: FIGURE 4 INTERPRETATION (UPDATED)
# ============================================================================

print("\n" + "="*80)
print("SECTION 7: FIGURE 4 INTERPRETATION (UPDATED)")
print("="*80)

print("\nPanel A: Modeled coast-specific upland WUE_T response")
print("  • Shows GAM-predicted WUE_T percent change across SPEI gradients")
print("  • Displays 5%, 10%, and 20% ecological impact thresholds")
print("  • Gulf Coast SPEI-3 shows strongest threshold behavior")
print("  • AK Coast shows thresholds under wet conditions at SPEI-3 and SPEI-48")
print("  • Atlantic/Pacific coasts show flatter responses")
print("  • Direction is mixed - WUE_T can increase OR decrease under dry/wet conditions")

print("\nPanel B: Observed site-level sensitivity by coastline")
print("  • Shows mean absolute WUE_T change from near-normal conditions")
print("  • Groups sites by coastline and SPEI timescale")
print("  • Gulf and AK coasts show highest observed sensitivity")
print("  • Atlantic/Pacific coasts show lower observed sensitivity")

print("\nPanel B provides observed site-level support for the coast-specific threshold patterns shown in Panel A.")
print("The GAMs were fit using all ecosystem classes, but Figure 4A predictions were standardized to")
print("upland conditions, July, and population-level site effects to isolate coast-specific SPEI threshold behavior.")

print("\nKEY INSIGHT (UPDATED):")
print("  The strongest modeled nonlinear response was Gulf Coast SPEI-3, while observed site-level")
print("  sensitivity was highest in the Gulf and Alaska coasts. This means SPEI thresholds are most")
print("  informative when interpreted by coastline and drought-memory timescale, not as one universal")
print("  SPEI cutoff. Coast-specific thresholds should be used explicitly in drought-impact workflows,")
print("  with threshold direction and uncertainty treated explicitly.")

# ============================================================================
# SECTION 8: MANUSCRIPT-READY VALUE SUMMARY
# ============================================================================

print("\n" + "="*80)
print("SECTION 8: MANUSCRIPT-READY VALUE SUMMARY")
print("="*80)

print("\nMETHODS:")
print(f"  • Dataset: {N_OBSERVATIONS:,} observations, {N_SITES} sites, {YEAR_START}-{YEAR_END} ({N_YEARS} years)")
print(f"  • Long format: {N_LONG_OBS:,} site-month-timescale observations")
print(f"  • {N_SPEI_TIMESCALES} SPEI timescales (1-48 months)")
print(f"  • 3 ecosystem classes: Upland ({ECOSYSTEM_DIST['Upland']['obs']} obs, {ECOSYSTEM_DIST['Upland']['sites']} sites), Freshwater ({ECOSYSTEM_DIST['Freshwater']['obs']} obs, {ECOSYSTEM_DIST['Freshwater']['sites']} sites), Saline ({ECOSYSTEM_DIST['Saline']['obs']} obs, {ECOSYSTEM_DIST['Saline']['sites']} sites)")
print(f"  • 4 coast regions: Atlantic ({COAST_DIST['Atlantic Coast']['obs']} obs, {COAST_DIST['Atlantic Coast']['sites']} sites), Pacific ({COAST_DIST['Pacific Coast']['obs']} obs, {COAST_DIST['Pacific Coast']['sites']} sites), Gulf ({COAST_DIST['Gulf Coast']['obs']} obs, {COAST_DIST['Gulf Coast']['sites']} sites), AK ({COAST_DIST['AK Coast']['obs']} obs, {COAST_DIST['AK Coast']['sites']} sites)")

print("\nMODEL PERFORMANCE:")
print(f"  • Coast-threshold GAM: {MODEL_COMP['coast']['deviance_explained']*100:.1f}% deviance explained, R²={MODEL_COMP['coast']['r_squared']:.3f}")
print(f"  • Ecosystem GAM: {MODEL_COMP['ecosystem']['deviance_explained']*100:.1f}% deviance explained, R²={MODEL_COMP['ecosystem']['r_squared']:.3f}")
print(f"  • ΔAIC: {MODEL_COMP['ecosystem']['AIC'] - MODEL_COMP['coast']['AIC']:.2f} (coast-threshold model better)")

print("\nSTRONGEST SIGNAL:")
print(f"  • {STRONGEST_SIGNAL['coast']} {STRONGEST_SIGNAL['timescale']}: EDF={STRONGEST_SIGNAL['edf']:.2f}, F={STRONGEST_SIGNAL['f_stat']:.2f}, p<0.001")

print("\nTHRESHOLDS:")
print(f"  • Gulf Coast SPEI-3 dry-side decreases: {THRESHOLDS['Gulf_SPEI3_dry_decrease']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_dry_decrease']['max']:.2f}")
print(f"  • Gulf Coast SPEI-3 wet-side increases: {THRESHOLDS['Gulf_SPEI3_wet_increase']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_wet_increase']['max']:.2f}")

print("\nSENSITIVITY:")
for coast, values in COAST_SENSITIVITY.items():
    print(f"  • {coast}: {values['SPEI_1']:.1f}% (SPEI-1), {values['SPEI_3']:.1f}% (SPEI-3), {values['SPEI_48']:.1f}% (SPEI-48)")

print("\nOVERALL INTERPRETATION:")
print("  • SPEI thresholds are most informative when interpreted by coastline and drought-memory timescale")
print("  • The strongest modeled response was Gulf Coast SPEI-3 (EDF=4.34, F=29.44, p<0.001)")
print("  • Observed sensitivity was highest in Gulf and Alaska coasts (∼40-45%)")
print("  • Atlantic and Pacific coasts showed lower sensitivity (∼16-20%)")
print("  • Coast-specific thresholds should be used in drought-impact workflows")
print("  • Threshold direction and uncertainty should be treated explicitly")

# ============================================================================
# SECTION 9: SAVE COMPLETE REPORT
# ============================================================================

print("\n" + "="*80)
print("SECTION 9: SAVING COMPLETE REPORT")
print("="*80)

report_file = os.path.join(output_dir, "Q3_COMPLETE_UPDATED_WORKFLOW_report.txt")
print(f"\nSaving report to: {report_file}")

try:
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("Q3 COMPLETE UPDATED WORKFLOW - METHODS & RESULTS EXTRACTION\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n\n")
        
        # Methods section
        f.write("METHODS:\n")
        f.write("-"*80 + "\n")
        f.write(f"Dataset: {N_OBSERVATIONS:,} observations, {N_SITES} sites, {YEAR_START}-{YEAR_END} ({N_YEARS} years)\n")
        f.write(f"Long format: {N_LONG_OBS:,} site-month-timescale observations\n")
        f.write(f"SPEI timescales: {N_SPEI_TIMESCALES} (1-48 months)\n\n")
        
        f.write("Ecosystem distribution:\n")
        for eco, stats in ECOSYSTEM_DIST.items():
            f.write(f"  {eco}: {stats['obs']} obs, {stats['sites']} sites\n")
        
        f.write("\nCoast distribution:\n")
        for coast, stats in COAST_DIST.items():
            f.write(f"  {coast}: {stats['obs']} obs, {stats['sites']} sites\n")
        
        f.write("\nModel Comparison (Table Q3.1):\n")
        f.write(f"  Ecosystem GAM: {MODEL_COMP['ecosystem']['deviance_explained']*100:.1f}% deviance, R²={MODEL_COMP['ecosystem']['r_squared']:.3f}, AIC={MODEL_COMP['ecosystem']['AIC']:.2f}\n")
        f.write(f"  Coast-threshold GAM: {MODEL_COMP['coast']['deviance_explained']*100:.1f}% deviance, R²={MODEL_COMP['coast']['r_squared']:.3f}, AIC={MODEL_COMP['coast']['AIC']:.2f}\n")
        f.write(f"  ΔAIC: {MODEL_COMP['ecosystem']['AIC'] - MODEL_COMP['coast']['AIC']:.2f}\n\n")
        
        # Results section
        f.write("RESULTS:\n")
        f.write("-"*80 + "\n")
        f.write(f"Strongest signal: {STRONGEST_SIGNAL['coast']} {STRONGEST_SIGNAL['timescale']} (EDF={STRONGEST_SIGNAL['edf']:.2f}, F={STRONGEST_SIGNAL['f_stat']:.2f}, p<0.001)\n\n")
        
        f.write("Coastline sensitivity medians (Table Q3.3):\n")
        for coast, values in COAST_SENSITIVITY.items():
            f.write(f"  {coast}: {values['SPEI_1']:.1f}% (SPEI-1), {values['SPEI_3']:.1f}% (SPEI-3), {values['SPEI_48']:.1f}% (SPEI-48)\n")
        
        f.write("\nThreshold summaries:\n")
        f.write(f"  Gulf Coast SPEI-3 dry-side decreases: {THRESHOLDS['Gulf_SPEI3_dry_decrease']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_dry_decrease']['max']:.2f}\n")
        f.write(f"  Gulf Coast SPEI-3 wet-side increases: {THRESHOLDS['Gulf_SPEI3_wet_increase']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_wet_increase']['max']:.2f}\n")
        f.write(f"  AK Coast SPEI-48 dry-side decreases: {THRESHOLDS['AK_SPEI48_dry_decrease']['min']:.2f}\n")
        f.write(f"  AK Coast SPEI-48 wet-side increases: {THRESHOLDS['AK_SPEI48_wet_increase']['min']:.2f}\n\n")
        
        # Interpretation
        f.write("INTERPRETATION:\n")
        f.write("-"*80 + "\n")
        f.write("Panel B provides observed site-level support for the coast-specific threshold patterns shown in Panel A.\n")
        f.write("The GAMs were fit using all ecosystem classes, but Figure 4A predictions were standardized to\n")
        f.write("upland conditions, July, and population-level site effects to isolate coast-specific SPEI threshold behavior.\n\n")
        
        f.write("KEY INSIGHT:\n")
        f.write("  The strongest modeled nonlinear response was Gulf Coast SPEI-3, while observed site-level\n")
        f.write("  sensitivity was highest in the Gulf and Alaska coasts. This means SPEI thresholds are most\n")
        f.write("  informative when interpreted by coastline and drought-memory timescale, not as one universal\n")
        f.write("  SPEI cutoff. Coast-specific thresholds should be used explicitly in drought-impact workflows,\n")
        f.write("  with threshold direction and uncertainty treated explicitly.\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"  ✓ Complete report saved successfully")
except Exception as e:
    print(f"  ERROR saving report: {e}")

print("\n" + "="*80)
print("Q3 COMPLETE UPDATED WORKFLOW COMPLETE")
print("="*80)
print(f"\nReport saved to: {report_file}")
print("\nKey findings confirmed:")
print(f"  • Dataset: {N_OBSERVATIONS:,} observations, {N_SITES} sites, {YEAR_START}-{YEAR_END}")
print(f"  • Best model: Coast-threshold GAM ({MODEL_COMP['coast']['deviance_explained']*100:.1f}% deviance, R²={MODEL_COMP['coast']['r_squared']:.3f})")
print(f"  • Strongest signal: {STRONGEST_SIGNAL['coast']} {STRONGEST_SIGNAL['timescale']} (EDF={STRONGEST_SIGNAL['edf']:.2f}, F={STRONGEST_SIGNAL['f_stat']:.2f})")
print(f"  • Highest sensitivity: Alaska Coast ({COAST_SENSITIVITY['Alaska']['SPEI_48']:.1f}% at SPEI-48)")
print(f"  • Key insight: Coast-specific SPEI thresholds are most informative for drought-impact workflows")
print("="*80)