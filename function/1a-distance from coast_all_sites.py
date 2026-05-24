# -*- coding: utf-8 -*-
"""
Corrected coastline distance calculation for AmeriFlux sites
Shows ONLY coastal sites (≤80km) on maps - with clear visual distinction
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, LineString, MultiLineString, box
from shapely.ops import nearest_points, unary_union
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os
from shapely.validation import make_valid
import chardet
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# STEP 1: Geographic filters
# ============================================================================

def in_conus_or_alaska(lat, lon):
    """Filter to CONUS and Alaska only using coordinate bounds"""
    in_conus = (24 <= lat <= 50) and (-125 <= lon <= -66)
    in_alaska = (51 <= lat <= 72) and (-180 <= lon <= -130)
    return in_conus or in_alaska

def get_region(lat, lon):
    """Return region for projection selection"""
    if 24 <= lat <= 50 and -125 <= lon <= -66:
        return 'CONUS'
    elif 51 <= lat <= 72 and -180 <= lon <= -130:
        return 'Alaska'
    else:
        return 'Other'

# ============================================================================
# STEP 2: Read CSV with encoding detection
# ============================================================================

def read_csv_with_encoding(file_path):
    """
    Read CSV file with automatic encoding detection
    """
    # Try common encodings first
    encodings_to_try = ['utf-8', 'cp1252', 'latin-1', 'iso-8859-1', 'cp1250', 'utf-16']
    
    for encoding in encodings_to_try:
        try:
            df = pd.read_csv(file_path, encoding=encoding)
            print(f"   Successfully read file with {encoding} encoding")
            return df
        except (UnicodeDecodeError, UnicodeError):
            continue
        except Exception as e:
            continue
    
    # If all fail, detect encoding with chardet
    try:
        print("   Attempting encoding detection with chardet...")
        with open(file_path, 'rb') as f:
            raw_data = f.read()
            result = chardet.detect(raw_data)
            detected_encoding = result['encoding']
            confidence = result['confidence']
            print(f"   Detected encoding: {detected_encoding} (confidence: {confidence:.2f})")
            
            df = pd.read_csv(file_path, encoding=detected_encoding)
            return df
    except Exception as e:
        raise Exception(f"Could not read CSV file. Tried multiple encodings. Last error: {str(e)}")

# ============================================================================
# STEP 3: Load and clean coastline data (CRITICAL: Remove Great Lakes)
# ============================================================================

def remove_great_lakes_from_coastline(coastline_gdf):
    """
    Remove Great Lakes shoreline from coastline data.
    """
    great_lakes_bboxes = [
        box(-92, 46, -84, 49),      # Lake Superior
        box(-88, 41, -86, 46),      # Lake Michigan
        box(-84, 43, -80, 46),      # Lake Huron
        box(-83, 41, -79, 43),      # Lake Erie
        box(-80, 43, -76, 44),      # Lake Ontario
        box(-88, 44, -82, 48)       # Combined northern lakes region
    ]
    
    keep_mask = pd.Series(True, index=coastline_gdf.index)
    
    for idx, geom in coastline_gdf.geometry.items():
        intersects_gl = False
        for gl_bbox in great_lakes_bboxes:
            try:
                if geom.intersects(gl_bbox):
                    intersects_gl = True
                    break
            except:
                continue
        
        if intersects_gl:
            keep_mask[idx] = False
    
    coastline_filtered = coastline_gdf[keep_mask].copy()
    
    gl_region_mask = ~((coastline_filtered.geometry.bounds.minx > -95) &
                       (coastline_filtered.geometry.bounds.maxx < -75) &
                       (coastline_filtered.geometry.bounds.miny > 40) &
                       (coastline_filtered.geometry.bounds.maxy < 50))
    
    coastline_filtered = coastline_filtered[gl_region_mask].copy()
    
    print(f"   Removed {len(coastline_gdf) - len(coastline_filtered)} Great Lakes shoreline features")
    
    return coastline_filtered


def load_and_clean_coastline(coastline_folder_path, region='CONUS'):
    """Load coastline shapefile and remove Great Lakes contamination."""
    shapefile_path = None
    for file in os.listdir(coastline_folder_path):
        if file.endswith('.shp'):
            shapefile_path = os.path.join(coastline_folder_path, file)
            print(f"   Found shapefile: {file}")
            break
    
    if shapefile_path is None:
        raise FileNotFoundError(f"No shapefile in {coastline_folder_path}")
    
    coastline = gpd.read_file(shapefile_path)
    print(f"   Loaded {len(coastline)} features")
    
    coastline = remove_great_lakes_from_coastline(coastline)
    coastline = coastline[~coastline.geometry.is_empty]
    coastline['geometry'] = coastline['geometry'].apply(
        lambda x: make_valid(x) if not x.is_valid else x
    )
    
    if region == 'CONUS':
        coastline = coastline.cx[-125:-66, 24:50]
        print(f"   Filtered to CONUS: {len(coastline)} features")
    elif region == 'Alaska':
        coastline = coastline.cx[-180:-130, 51:72]
        print(f"   Filtered to Alaska: {len(coastline)} features")
    
    return coastline


# ============================================================================
# STEP 4: Project to appropriate CRS
# ============================================================================

def get_projection_for_region(region):
    """Return appropriate projected CRS for accurate distance calculations"""
    if region == 'CONUS':
        return 'EPSG:5070'  # USA Contiguous Albers Equal Area Conic
    elif region == 'Alaska':
        return 'EPSG:3338'  # Alaska Albers Equal Area
    else:
        return 'EPSG:4326'


def calculate_distances_by_region(sites_gdf, coastline_gdf, region):
    """Calculate distances for ALL sites using union of coastline geometries."""
    print(f"\n   Processing {region} sites...")
    
    proj_crs = get_projection_for_region(region)
    
    sites_proj = sites_gdf.to_crs(proj_crs)
    coastline_proj = coastline_gdf.to_crs(proj_crs)
    
    print(f"      Creating coastline union...")
    coast_union = unary_union(coastline_proj.geometry.values)
    print(f"      Coastline union created")
    
    distances_km = []
    nearest_points_list = []
    
    for i, (idx, site) in enumerate(sites_proj.iterrows(), 1):
        site_point = site.geometry
        
        try:
            nearest = nearest_points(site_point, coast_union)[1]
            dist_m = site_point.distance(nearest)
            dist_km = dist_m / 1000.0
            
            distances_km.append(dist_km)
            nearest_points_list.append(nearest)
            
        except Exception as e:
            distances_km.append(np.nan)
            nearest_points_list.append(None)
        
        if i % 50 == 0 or i == len(sites_proj):
            print(f"      Processed {i}/{len(sites_proj)} sites")
    
    return distances_km, nearest_points_list


# ============================================================================
# STEP 5: Plot maps showing ONLY coastal sites with clear visual distinction
# ============================================================================

def plot_coastal_maps_only(coastal_sites, conus_coastline, alaska_coastline, max_distance_km=80):
    """
    Create separate maps for CONUS and Alaska showing ONLY coastal sites (≤80km).
    Uses a single color (red) to avoid confusion between distance colors and coastline.
    """
    
    print("\n" + "=" * 80)
    print("CREATING MAPS - SHOWING ONLY COASTAL SITES (≤80km)")
    print("=" * 80)
    
    # Split coastal sites by region (these are already filtered to ≤80km)
    conus_coastal = coastal_sites[coastal_sites['region'] == 'CONUS'].copy()
    alaska_coastal = coastal_sites[coastal_sites['region'] == 'Alaska'].copy()
    
    # DEBUG: Verify only coastal sites are being plotted
    print("\n   DEBUG CHECK BEFORE PLOTTING:")
    print(f"   CONUS coastal sites to plot: {len(conus_coastal)}")
    print(f"   Alaska coastal sites to plot: {len(alaska_coastal)}")
    if len(conus_coastal) > 0:
        print(f"   Any CONUS sites > {max_distance_km} km?: {(conus_coastal['distance_km'] > max_distance_km).any()}")
        if len(conus_coastal[conus_coastal['distance_km'] >= 70]) > 0:
            print(f"   CONUS sites near threshold (70-80km):")
            near_threshold = conus_coastal[conus_coastal['distance_km'] >= 70].sort_values('distance_km', ascending=False)
            for _, site in near_threshold.head(10).iterrows():
                print(f"      {site['site_id']}: {site['distance_km']:.2f} km")
    
    if len(alaska_coastal) > 0:
        print(f"   Any Alaska sites > {max_distance_km} km?: {(alaska_coastal['distance_km'] > max_distance_km).any()}")
        if len(alaska_coastal[alaska_coastal['distance_km'] >= 70]) > 0:
            print(f"   Alaska sites near threshold (70-80km):")
            near_threshold = alaska_coastal[alaska_coastal['distance_km'] >= 70].sort_values('distance_km', ascending=False)
            for _, site in near_threshold.head(10).iterrows():
                print(f"      {site['site_id']}: {site['distance_km']:.2f} km")
    
    # Create figure with two subplots side by side
    fig, (ax_conus, ax_alaska) = plt.subplots(1, 2, figsize=(20, 8))
    
    # ========================================================================
    # CONUS MAP - ONLY COASTAL SITES (all in red)
    # ========================================================================
    print("\n   Plotting CONUS map (coastal sites only in RED)...")
    
    # Plot coastline in light gray (less visually dominant)
    if len(conus_coastline) > 0:
        conus_coastline.boundary.plot(ax=ax_conus, color='lightgray', linewidth=0.6, alpha=0.7, label='Coastline')
    
    # Plot ONLY coastal sites (within 80km) - ALL IN RED to avoid confusion
    if len(conus_coastal) > 0:
        ax_conus.scatter(conus_coastal['lon'], conus_coastal['lat'],
                        color='red', s=120, alpha=0.9, 
                        edgecolors='darkred', linewidth=1.5, zorder=5)
        
        # Add labels only for sites to show they are coastal
        # But limit labels to avoid clutter - only label if less than 50 sites
        if len(conus_coastal) <= 50:
            for _, site in conus_coastal.iterrows():
                ax_conus.annotate(f"{site['site_id']}\n{site['distance_km']:.1f}km", 
                                 (site['lon'], site['lat']), 
                                 xytext=(8, 8), textcoords='offset points',
                                 fontsize=7, alpha=0.8,
                                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor='gray'))
        else:
            # Only label closest and farthest to avoid clutter
            closest = conus_coastal.nsmallest(5, 'distance_km')
            farthest = conus_coastal.nlargest(5, 'distance_km')
            for _, site in pd.concat([closest, farthest]).iterrows():
                ax_conus.annotate(f"{site['site_id']}\n{site['distance_km']:.1f}km", 
                                 (site['lon'], site['lat']), 
                                 xytext=(8, 8), textcoords='offset points',
                                 fontsize=7, alpha=0.8,
                                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor='gray'))
    
    # Set CONUS map properties
    ax_conus.set_xlim(-125, -66)
    ax_conus.set_ylim(24, 50)
    ax_conus.set_title(f'CONTIGUOUS USA\nCOASTAL SITES ONLY (≤ {max_distance_km} km)\nTotal: {len(conus_coastal)} sites', 
                      fontsize=14, fontweight='bold', pad=20, color='darkred')
    ax_conus.set_xlabel('Longitude', fontsize=11)
    ax_conus.set_ylabel('Latitude', fontsize=11)
    ax_conus.grid(True, alpha=0.3, linestyle='--')
    
    # Add legend
    legend_elements = [Line2D([0], [0], color='lightgray', linewidth=1.5, label='Coastline'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                            markersize=10, label=f'Coastal Sites (≤{max_distance_km}km)', 
                            markeredgecolor='darkred')]
    ax_conus.legend(handles=legend_elements, loc='lower left', fontsize=10, framealpha=0.9)
    
    # ========================================================================
    # ALASKA MAP - ONLY COASTAL SITES (all in red)
    # ========================================================================
    print("   Plotting Alaska map (coastal sites only in RED)...")
    
    # Plot coastline in light gray
    if len(alaska_coastline) > 0:
        alaska_coastline.boundary.plot(ax=ax_alaska, color='lightgray', linewidth=0.6, alpha=0.7, label='Coastline')
    
    # Plot ONLY coastal sites (within 80km) - ALL IN RED
    if len(alaska_coastal) > 0:
        ax_alaska.scatter(alaska_coastal['lon'], alaska_coastal['lat'],
                         color='red', s=120, alpha=0.9, 
                         edgecolors='darkred', linewidth=1.5, zorder=5)
        
        # Add labels for all Alaska sites (usually fewer)
        for _, site in alaska_coastal.iterrows():
            ax_alaska.annotate(f"{site['site_id']}\n{site['distance_km']:.1f}km", 
                              (site['lon'], site['lat']), 
                              xytext=(8, 8), textcoords='offset points',
                              fontsize=8, alpha=0.9,
                              bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor='gray'))
    
    # Set Alaska map properties
    ax_alaska.set_xlim(-180, -130)
    ax_alaska.set_ylim(51, 72)
    ax_alaska.set_title(f'ALASKA\nCOASTAL SITES ONLY (≤ {max_distance_km} km)\nTotal: {len(alaska_coastal)} sites', 
                       fontsize=14, fontweight='bold', pad=20, color='darkred')
    ax_alaska.set_xlabel('Longitude', fontsize=11)
    ax_alaska.set_ylabel('Latitude', fontsize=11)
    ax_alaska.grid(True, alpha=0.3, linestyle='--')
    
    # Add legend for Alaska
    ax_alaska.legend(handles=legend_elements, loc='lower left', fontsize=10, framealpha=0.9)
    
    # Add overall title
    plt.suptitle(f'COASTAL SITE ANALYSIS - FINAL RETAINED SITES\nOcean Coastline Only (Great Lakes Excluded) | Distance ≤ {max_distance_km} km\nAll sites shown in RED', 
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.show()
    
    print(f"\n   ✓ Maps created showing ONLY coastal sites (≤{max_distance_km} km):")
    print(f"      - CONUS: {len(conus_coastal)} coastal sites displayed (ALL IN RED)")
    print(f"      - Alaska: {len(alaska_coastal)} coastal sites displayed (ALL IN RED)")
    print(f"      - Coastline shown in light gray for clarity")
    print(f"      - NO non-coastal sites (> {max_distance_km} km) are shown")
    print(f"      - Blue colors in previous maps were coastal sites with distances near {max_distance_km} km")
    
    return fig


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def analyze_coastal_sites(site_csv_path, coastline_folder_path, max_distance_km=80):
    """
    Identify sites within max_distance_km of ocean coastline.
    Returns only coastal sites for mapping.
    """
    
    print("=" * 80)
    print("COASTAL SITE ANALYSIS - FINAL RETAINED SITES ONLY")
    print("=" * 80)
    print(f"Input CSV: {site_csv_path}")
    print(f"Distance threshold: {max_distance_km} km")
    print("Geographic domain: CONUS + Alaska only (excludes Great Lakes)")
    print("=" * 80)
    
    # ------------------------------------------------------------------------
    # Step 1: Load site data with encoding detection
    # ------------------------------------------------------------------------
    print("\n1. LOADING SITE DATA")
    print("-" * 40)
    
    # Read CSV with encoding detection
    sites_df = read_csv_with_encoding(site_csv_path)
    print(f"   Total sites in file: {len(sites_df)}")
    print(f"   Columns: {list(sites_df.columns)}")
    
    # Identify columns (handle your exact naming)
    lat_col = None
    lon_col = None
    site_col = None
    
    # Look for latitude column
    for col in sites_df.columns:
        col_lower = col.lower()
        if 'lat' in col_lower and not 'long' in col_lower and not 'lon' in col_lower:
            lat_col = col
            break
    
    # Look for longitude column
    for col in sites_df.columns:
        col_lower = col.lower()
        if 'lon' in col_lower or 'long' in col_lower:
            lon_col = col
            break
    
    # Look for site ID column
    for col in sites_df.columns:
        col_lower = col.lower()
        if 'site' in col_lower and ('id' in col_lower or 'name' in col_lower):
            site_col = col
            break
    
    if lat_col is None or lon_col is None:
        raise ValueError(f"Cannot find lat/lon columns. Available columns: {list(sites_df.columns)}")
    
    # Rename for consistency
    sites_df = sites_df.rename(columns={
        lat_col: 'lat',
        lon_col: 'lon'
    })
    
    if site_col:
        sites_df = sites_df.rename(columns={site_col: 'site_id'})
    else:
        sites_df['site_id'] = sites_df.index
    
    print(f"   Using: lat='{lat_col}', lon='{lon_col}'")
    print(f"   Site ID column: '{site_col if site_col else 'auto-generated'}'")
    
    # Show first few rows to verify
    print(f"\n   First 5 rows of data:")
    print(sites_df[['site_id', 'lat', 'lon']].head())
    
    # ------------------------------------------------------------------------
    # Step 2: Filter to CONUS + Alaska by coordinates
    # ------------------------------------------------------------------------
    print("\n2. FILTERING TO CONUS + ALASKA")
    print("-" * 40)
    
    mask = sites_df.apply(lambda row: in_conus_or_alaska(row['lat'], row['lon']), axis=1)
    original_count = len(sites_df)
    sites_df = sites_df[mask].copy()
    
    print(f"   Sites in CONUS/Alaska: {len(sites_df)} (removed {original_count - len(sites_df)})")
    
    if len(sites_df) == 0:
        print("   No sites in target region. Exiting.")
        return None, None, None, None, None
    
    # Add region column
    sites_df['region'] = sites_df.apply(lambda row: get_region(row['lat'], row['lon']), axis=1)
    print(f"   Regions: {sites_df['region'].value_counts().to_dict()}")
    
    # Create geometry
    geometry = [Point(xy) for xy in zip(sites_df['lon'], sites_df['lat'])]
    sites_gdf = gpd.GeoDataFrame(sites_df, geometry=geometry, crs="EPSG:4326")
    
    # ------------------------------------------------------------------------
    # Step 3: Load and clean coastline data (remove Great Lakes)
    # ------------------------------------------------------------------------
    print("\n3. LOADING AND CLEANING COASTLINE DATA")
    print("-" * 40)
    
    conus_coastline = load_and_clean_coastline(coastline_folder_path, region='CONUS')
    alaska_coastline = load_and_clean_coastline(coastline_folder_path, region='Alaska')
    
    print(f"   CONUS coastline features: {len(conus_coastline)}")
    print(f"   Alaska coastline features: {len(alaska_coastline)}")
    
    # ------------------------------------------------------------------------
    # Step 4: Calculate distances for ALL sites
    # ------------------------------------------------------------------------
    print("\n4. CALCULATING DISTANCES FOR ALL SITES")
    print("-" * 40)
    
    conus_sites = sites_gdf[sites_gdf['region'] == 'CONUS'].copy()
    alaska_sites = sites_gdf[sites_gdf['region'] == 'Alaska'].copy()
    
    print(f"   CONUS sites to analyze: {len(conus_sites)}")
    print(f"   Alaska sites to analyze: {len(alaska_sites)}")
    
    all_distances = {}
    all_nearest = {}
    
    if len(conus_sites) > 0:
        dist, nearest = calculate_distances_by_region(conus_sites, conus_coastline, 'CONUS')
        all_distances['CONUS'] = dist
        all_nearest['CONUS'] = nearest
    
    if len(alaska_sites) > 0:
        dist, nearest = calculate_distances_by_region(alaska_sites, alaska_coastline, 'Alaska')
        all_distances['Alaska'] = dist
        all_nearest['Alaska'] = nearest
    
    # Combine results
    sites_gdf['distance_km'] = np.nan
    sites_gdf['nearest_point'] = None
    
    if len(conus_sites) > 0:
        sites_gdf.loc[conus_sites.index, 'distance_km'] = all_distances['CONUS']
        sites_gdf.loc[conus_sites.index, 'nearest_point'] = all_nearest['CONUS']
    
    if len(alaska_sites) > 0:
        sites_gdf.loc[alaska_sites.index, 'distance_km'] = all_distances['Alaska']
        sites_gdf.loc[alaska_sites.index, 'nearest_point'] = all_nearest['Alaska']
    
    nan_count = sites_gdf['distance_km'].isna().sum()
    if nan_count > 0:
        print(f"\n   WARNING: {nan_count} sites had calculation issues")
        sites_gdf = sites_gdf.dropna(subset=['distance_km'])
    
    print(f"\n   ✓ Successfully calculated distances for {len(sites_gdf)} sites")
    
    # ------------------------------------------------------------------------
    # Step 5: Filter to coastal sites within threshold (RETAINED SITES)
    # ------------------------------------------------------------------------
    print("\n5. FILTERING TO RETAINED COASTAL SITES")
    print("-" * 40)
    
    coastal_sites = sites_gdf[sites_gdf['distance_km'] <= max_distance_km].copy()
    coastal_sites = coastal_sites.sort_values('distance_km')
    
    non_coastal_sites = sites_gdf[sites_gdf['distance_km'] > max_distance_km].copy()
    
    print(f"   Total sites analyzed: {len(sites_gdf)}")
    print(f"   ✓ RETAINED coastal sites (≤ {max_distance_km} km): {len(coastal_sites)}")
    print(f"   ✗ EXCLUDED sites (> {max_distance_km} km): {len(non_coastal_sites)}")
    if len(sites_gdf) > 0:
        print(f"   Retention rate: {len(coastal_sites)/len(sites_gdf)*100:.1f}%")
    
    # Show distribution of distances among retained sites
    if len(coastal_sites) > 0:
        print(f"\n   Distance distribution among RETAINED coastal sites:")
        print(f"      Min: {coastal_sites['distance_km'].min():.2f} km")
        print(f"      Max: {coastal_sites['distance_km'].max():.2f} km")
        print(f"      Mean: {coastal_sites['distance_km'].mean():.2f} km")
        print(f"      Median: {coastal_sites['distance_km'].median():.2f} km")
        print(f"      Sites near threshold (70-80 km): {len(coastal_sites[coastal_sites['distance_km'] >= 70])}")
    
    # ------------------------------------------------------------------------
    # Step 6: Display results summary
    # ------------------------------------------------------------------------
    if len(coastal_sites) > 0:
        print("\n" + "=" * 80)
        print("RETAINED COASTAL SITES SUMMARY")
        print("=" * 80)
        
        print("\nBy Region:")
        print("-" * 40)
        for region in ['CONUS', 'Alaska']:
            region_sites = coastal_sites[coastal_sites['region'] == region]
            if len(region_sites) > 0:
                print(f"\n{region}:")
                print(f"   Count: {len(region_sites)} sites (RETAINED)")
                print(f"   Min distance: {region_sites['distance_km'].min():.2f} km")
                print(f"   Max distance: {region_sites['distance_km'].max():.2f} km")
                print(f"   Mean distance: {region_sites['distance_km'].mean():.2f} km")
                print(f"   Median distance: {region_sites['distance_km'].median():.2f} km")
        
        print("\n" + "-" * 40)
        print("Distance Bands (RETAINED Coastal Sites):")
        print("-" * 40)
        for band in [10, 20, 40, 60, 80]:
            count = len(coastal_sites[coastal_sites['distance_km'] <= band])
            print(f"  ≤ {band:2d} km: {count:3d} sites ({count/len(coastal_sites)*100:5.1f}% of retained sites)")
        
        print("\n" + "-" * 40)
        print("RETAINED SITES NEAR THRESHOLD (70-80 km):")
        print("-" * 40)
        near_threshold = coastal_sites[coastal_sites['distance_km'] >= 70].sort_values('distance_km', ascending=False)
        if len(near_threshold) > 0:
            near_display = near_threshold[['site_id', 'region', 'lat', 'lon', 'distance_km']].copy()
            near_display['distance_km'] = near_display['distance_km'].round(2)
            print(near_display.to_string(index=False))
        else:
            print("   No sites near the threshold")
    
    # ------------------------------------------------------------------------
    # Step 7: Create maps showing ONLY coastal sites (retained sites) - ALL IN RED
    # ------------------------------------------------------------------------
    if len(coastal_sites) > 0:
        plot_coastal_maps_only(coastal_sites, conus_coastline, alaska_coastline, max_distance_km)
    else:
        print("\n   No coastal sites to display in maps (all sites > 80km)")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nFINAL SUMMARY:")
    print(f"   Total sites in CONUS + Alaska: {len(sites_gdf)}")
    print(f"   ✓ RETAINED coastal sites (≤ {max_distance_km} km): {len(coastal_sites)}")
    print(f"   ✗ EXCLUDED sites (> {max_distance_km} km): {len(non_coastal_sites)}")
    if len(sites_gdf) > 0:
        print(f"   Retention rate: {len(coastal_sites)/len(sites_gdf)*100:.1f}%")
    print("\n   Great Lakes shoreline REMOVED from coastline data")
    print("   ✓ Maps show ONLY retained coastal sites (≤80km) - ALL IN RED")
    print("   ✗ Non-coastal sites (>80km) are NOT displayed on maps")
    print("   ✓ Coastline shown in light gray for clarity")
    print("\n   NOTE: Previous 'blue' sites were coastal sites with distances near 80km")
    
    return sites_gdf, coastal_sites, non_coastal_sites, conus_coastline, alaska_coastline


# ============================================================================
# EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Paths
    #site_csv_path = r"M:\Research\WUE_CUE\data_products\info\all_sites.csv"
    site_csv_path = r"M:\Research\WUE_CUE\data_products\info\site_lat_long.csv"
    coastline_folder_path = r"M:\Research\WUE_CUE\spatial_SPEI\coastline"
    
    try:
        all_sites, coastal_sites, excluded_sites, conus_coast, ak_coast = analyze_coastal_sites(
            site_csv_path, 
            coastline_folder_path,
            max_distance_km=80
        )
        
        # Final validation
        if coastal_sites is not None and len(coastal_sites) > 0:
            print("\n" + "=" * 80)
            print("VALIDATION CHECK")
            print("=" * 80)
            print("✓ Maps show ONLY coastal sites (≤80km) - ALL IN RED")
            print(f"✓ Non-coastal sites ({len(excluded_sites)} sites) are NOT displayed")
            print("✓ Coastline shown in light gray")
            
            # Check for Great Lakes in retained sites
            gl_sites = coastal_sites[
                (coastal_sites['lat'] > 40) & (coastal_sites['lat'] < 50) &
                (coastal_sites['lon'] > -95) & (coastal_sites['lon'] < -75)
            ]
            if len(gl_sites) > 0:
                print(f"\n⚠ WARNING: {len(gl_sites)} Great Lakes sites found in retained coastal results:")
                print(gl_sites[['site_id', 'lat', 'lon', 'distance_km']])
            else:
                print("\n✓ No Great Lakes sites found in retained coastal results")
            
            print("\n" + "=" * 80)
            print("INTERPRETATION NOTE:")
            print("=" * 80)
            print("If you previously saw blue site markers, those were coastal sites")
            print("with distances close to 80km (blue in the RdYlBu_r colormap).")
            print("Now ALL coastal sites are shown in RED to avoid confusion.")
            print("The maps now clearly show ONLY the retained coastal sites.")
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()