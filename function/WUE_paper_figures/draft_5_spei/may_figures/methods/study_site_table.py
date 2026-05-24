# -*- coding: utf-8 -*-
"""
TABLE 1: Characteristics of EC flux tower sites
Uses SAME filtering as upstream workflow (triple intersection)
Extracts actual year ranges and T:ET ranges from monthly data
UPDATED: Shortened climate names for cleaner display
"""

import pandas as pd
import numpy as np
from tabulate import tabulate
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# FILE PATHS
# ============================================================================

MONTHLY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
NN_MEDIANS_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"
SUMMARY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\site_summary_with_SPEI48.csv"
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# CLIMATE NAME MAPPING (SHORTENED VERSIONS)
# ============================================================================

def shorten_climate_names(df):
    """Convert long climate names to shorter, cleaner versions"""
    
    if 'climate' not in df.columns:
        return df
    
    # Define mapping from long names to short names
    climate_map = {
        'Tundra': 'Tundra',
        'Humid Subtropical: mild, no dry season, hot summer': 'Humid Subtropical (mild)',
        'Humid Subtropical: dry winter, hot summer': 'Humid Subtropical (dry winter)',
        'Mediterranean: mild, dry, hot summer': 'Mediterranean (hot summer)',
        'Mediterranean: mild, dry, warm summer': 'Mediterranean (warm summer)',
        'Marine West Coast: mild, no dry season, cool summer': 'Marine West Coast',
        'Dry Continental: cool summer': 'Dry Continental',
        'Humid Continental: humid, severe winter, no dry season, hot summer': 'Humid Continental'
    }
    
    # Apply mapping
    df['climate_short'] = df['climate'].map(climate_map).fillna(df['climate'])
    
    # Replace the original climate column with shortened version
    df['climate'] = df['climate_short']
    df = df.drop('climate_short', axis=1)
    
    return df

# ============================================================================
# FUNCTION: GET TRIPLE INTERSECTION SITES (MATCHES UPSTREAM)
# ============================================================================

def get_triple_intersection_sites():
    """Get the strict triple intersection sites matching upstream workflow"""
    
    print("\n📊 Getting strict triple intersection sites...")
    
    nn_data = pd.read_csv(NN_MEDIANS_PATH)
    
    wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
    eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
    tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
    
    shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
    
    print(f"  WUE sites: {len(wue_sites)}")
    print(f"  WUE_eva sites: {len(eva_sites)}")
    print(f"  WUE_tra sites: {len(tra_sites)}")
    print(f"  Triple intersection: {len(shared_sites)} sites")
    
    return shared_sites

# ============================================================================
# FUNCTION: EXTRACT YEAR RANGES FROM MONTHLY DATA (USING YEAR COLUMN)
# ============================================================================

def extract_year_ranges_from_monthly(df_monthly, shared_sites):
    """
    Extract actual year ranges from monthly data using 'year' column
    Handles discontinuous ranges: e.g., "2012-2014, 2021-2025"
    """
    
    print("\n📅 Extracting actual year ranges from monthly data...")
    
    # Filter to shared sites only
    df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
    print(f"  Filtered to {df_monthly['site_name'].nunique()} shared sites")
    
    # Check if 'year' column exists
    if 'year' not in df_monthly.columns:
        print("  ⚠️ No 'year' column found!")
        return {}
    
    year_ranges = {}
    
    for site_name in shared_sites:
        site_data = df_monthly[df_monthly['site_name'] == site_name]
        
        # Get unique years, sorted
        unique_years = sorted(site_data['year'].dropna().unique())
        
        if len(unique_years) == 0:
            year_ranges[site_name] = 'N/A'
            continue
        
        # Find consecutive year ranges
        ranges = []
        start_year = unique_years[0]
        end_year = unique_years[0]
        
        for i in range(1, len(unique_years)):
            if unique_years[i] == unique_years[i-1] + 1:
                # Consecutive, extend the range
                end_year = unique_years[i]
            else:
                # Gap found, save current range
                if start_year == end_year:
                    ranges.append(f"{start_year}")
                else:
                    ranges.append(f"{start_year}-{end_year}")
                start_year = unique_years[i]
                end_year = unique_years[i]
        
        # Add the last range
        if start_year == end_year:
            ranges.append(f"{start_year}")
        else:
            ranges.append(f"{start_year}-{end_year}")
        
        # Join ranges with commas
        year_ranges[site_name] = ', '.join(ranges)
    
    print(f"  ✓ Extracted year ranges for {len(year_ranges)} sites")
    
    # Show examples
    print("\n  Example year ranges:")
    for i, (site, yr_range) in enumerate(list(year_ranges.items())[:5]):
        print(f"    {site}: {yr_range}")
    
    return year_ranges

# ============================================================================
# FUNCTION: EXTRACT T:ET RANGE (MIN AND MAX) FROM MONTHLY DATA
# ============================================================================

def extract_tet_range_from_monthly(df_monthly, shared_sites):
    """
    Extract T:ET ratio min and max from monthly data using 'Trans_ratio' column
    """
    
    print("\n📈 Extracting T:ET range (min and max) from monthly data...")
    
    # Filter to shared sites only
    df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
    
    tet_min = {}
    tet_max = {}
    
    # Check for Trans_ratio column
    if 'Trans_ratio' not in df_monthly.columns:
        print("  ⚠️ No 'Trans_ratio' column found! Available columns:")
        for col in df_monthly.columns:
            if 'trans' in col.lower() or 'ratio' in col.lower():
                print(f"    - {col}")
        return tet_min, tet_max
    
    for site_name in shared_sites:
        site_data = df_monthly[df_monthly['site_name'] == site_name]
        tet_values = site_data['Trans_ratio'].dropna()
        
        if len(tet_values) > 0:
            tet_min[site_name] = tet_values.min()
            tet_max[site_name] = tet_values.max()
        else:
            tet_min[site_name] = np.nan
            tet_max[site_name] = np.nan
    
    print(f"  ✓ Extracted T:ET ranges for {len(tet_min)} sites")
    
    # Show examples
    print("\n  Example T:ET ranges (min - max):")
    for i, (site, tmin) in enumerate(list(tet_min.items())[:5]):
        tmax = tet_max.get(site, np.nan)
        if not np.isnan(tmin) and not np.isnan(tmax):
            print(f"    {site}: {tmin:.3f} - {tmax:.3f}")
        else:
            print(f"    {site}: N/A")
    
    return tet_min, tet_max

# ============================================================================
# FUNCTION: FORMAT T:ET RANGE WITH 2 SIGNIFICANT DIGITS
# ============================================================================

def format_tet_range(tmin, tmax):
    """Format T:ET range with 2 significant digits"""
    if pd.isna(tmin) or pd.isna(tmax):
        return 'N/A'
    
    # Format with 2 significant digits
    def sig_fig(x, sig=2):
        if x == 0:
            return 0
        return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)
    
    tmin_fmt = sig_fig(tmin, 2)
    tmax_fmt = sig_fig(tmax, 2)
    
    # Handle formatting
    if tmin_fmt == int(tmin_fmt):
        tmin_fmt = int(tmin_fmt)
    if tmax_fmt == int(tmax_fmt):
        tmax_fmt = int(tmax_fmt)
    
    return f"{tmin_fmt} - {tmax_fmt}"

# ============================================================================
# FUNCTION: STANDARDIZE CLIMATE NAMES
# ============================================================================

def standardize_climate_names(df):
    """Standardize AND SHORTEN climate names"""
    
    if 'climate' not in df.columns:
        return df
    
    # First remove trailing parentheses and clean up
    df['climate'] = df['climate'].astype(str).str.rstrip(')')
    
    # Replace ' with ' with ', '
    df['climate'] = df['climate'].str.replace(' with ', ', ')
    
    # Clean up
    df['climate'] = df['climate'].str.replace('  ', ' ').str.strip()
    
    # THEN shorten the climate names
    df = shorten_climate_names(df)
    
    return df

# ============================================================================
# BIOME MAPPING
# ============================================================================

IGBP_MAP = {
    'ENF': 'Evergreen Needleleaf Forest',
    'EBF': 'Evergreen Broadleaf Forest',
    'DNF': 'Deciduous Needleleaf Forest',
    'DBF': 'Deciduous Broadleaf Forest',
    'MF': 'Mixed Forest',
    'CSH': 'Closed Shrublands',
    'OSH': 'Open Shrublands',
    'WSA': 'Woody Savannas',
    'SAV': 'Savannas',
    'GRA': 'Grasslands',
    'WET': 'Wetlands',
    'CRO': 'Croplands',
    'URB': 'Urban',
    'SNO': 'Snow and Ice',
    'BSV': 'Barren or Sparse Vegetation'
}

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Create Table 1 with actual year ranges and T:ET ranges"""
    
    print("="*80)
    print("📊 TABLE 1: Site Characteristics with Actual Year Ranges")
    print("="*80)
    
    # 1. Get triple intersection sites (matching upstream)
    shared_sites = get_triple_intersection_sites()
    
    # 2. Load monthly data
    print(f"\n📂 Loading monthly data: {MONTHLY_DATA_PATH}")
    df_monthly = pd.read_csv(MONTHLY_DATA_PATH)
    print(f"  ✓ Loaded {len(df_monthly):,} rows, {df_monthly['site_name'].nunique()} sites")
    
    # 3. Load summary data
    print(f"\n📂 Loading summary data: {SUMMARY_DATA_PATH}")
    df_summary = pd.read_csv(SUMMARY_DATA_PATH)
    print(f"  ✓ Loaded {len(df_summary)} sites")
    
    # 4. Filter summary to shared sites that are in both datasets
    # Get sites that are in monthly data AND summary data
    monthly_sites = set(df_monthly[df_monthly['site_name'].isin(shared_sites)]['site_name'].unique())
    summary_sites = set(df_summary['site_name'].unique())
    final_sites = monthly_sites.intersection(summary_sites)
    
    print(f"\n  Monthly shared sites: {len(monthly_sites)}")
    print(f"  Summary sites: {len(summary_sites)}")
    print(f"  Final sites for table: {len(final_sites)}")
    
    # Filter both dataframes to final sites
    df_monthly = df_monthly[df_monthly['site_name'].isin(final_sites)]
    df_summary = df_summary[df_summary['site_name'].isin(final_sites)]
    
    # 5. Extract year ranges from monthly data
    year_ranges = extract_year_ranges_from_monthly(df_monthly, final_sites)
    
    # 6. Extract T:ET ranges from monthly data
    tet_min_dict, tet_max_dict = extract_tet_range_from_monthly(df_monthly, final_sites)
    
    # 7. Standardize and shorten climate names
    df_summary = standardize_climate_names(df_summary)
    
    # 8. Create Table 1 with correct column order
    print("\n📊 Creating Table 1...")
    
    table_df = pd.DataFrame()
    
    # Column 1: Ameriflux Site
    table_df['Ameriflux Site'] = df_summary['site_name']
    
    # Column 2: Latitude
    table_df['Latitude'] = df_summary['lat'].round(2)
    
    # Column 3: Longitude
    table_df['Longitude'] = df_summary['long'].round(2)
    
    # Column 4: Climate (now with shortened names)
    table_df['Climate'] = df_summary['climate']
    
    # Column 5: Biome (from IGBP)
    if 'IGBP' in df_summary.columns:
        table_df['Biome'] = df_summary['IGBP'].map(IGBP_MAP).fillna(df_summary['IGBP'])
    else:
        table_df['Biome'] = 'N/A'
    
    # Column 6: Ecosystem Type
    if 'Salinity_Category' in df_summary.columns:
        table_df['Ecosystem Type'] = df_summary['Salinity_Category']
    else:
        table_df['Ecosystem Type'] = 'N/A'
    
    # Column 7: Data Range (actual years from monthly data)
    table_df['Data Range'] = table_df['Ameriflux Site'].map(year_ranges).fillna('N/A')
    
    # Column 8: T:ET Range (formatted with 2 significant digits)
    tet_range_list = []
    for site in table_df['Ameriflux Site']:
        tmin = tet_min_dict.get(site, np.nan)
        tmax = tet_max_dict.get(site, np.nan)
        tet_range_list.append(format_tet_range(tmin, tmax))
    table_df['T:ET Range'] = tet_range_list
    
    # 9. Sort by site name
    table_df = table_df.sort_values('Ameriflux Site').reset_index(drop=True)
    
    # 10. Print table to console
    print("\n" + "="*160)
    print("TABLE 1. Characteristics of EC flux tower sites used in this study")
    print("="*160)
    print(f"\n📍 Total sites included: {len(table_df)}")
    
    # Print unique climate types after shortening
    print(f"\n🌍 Climate types (after shortening):")
    for climate in sorted(table_df['Climate'].unique()):
        count = len(table_df[table_df['Climate'] == climate])
        print(f"   {climate}: {count} sites")
    
    print("\n" + tabulate(table_df, 
                          headers='keys', 
                          tablefmt='grid', 
                          maxcolwidths=[18, 10, 12, 25, 30, 16, 20, 18],
                          stralign='left'))
    
    # 11. Save to files
    csv_path = os.path.join(OUTPUT_DIR, "Table1_Site_Characteristics.csv")
    table_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    print(f"\n✅ CSV saved: {csv_path}")
    
    excel_path = os.path.join(OUTPUT_DIR, "Table1_Site_Characteristics.xlsx")
    try:
        table_df.to_excel(excel_path, index=False, engine='openpyxl')
        print(f"✅ Excel saved: {excel_path}")
    except Exception as e:
        print(f"⚠️ Could not save Excel: {e}")
    
    # 12. Summary statistics
    print("\n" + "="*80)
    print("📊 SUMMARY STATISTICS")
    print("="*80)
    
    print(f"\n🏞️ Ecosystem Types:")
    for eco, count in table_df['Ecosystem Type'].value_counts().items():
        print(f"   {eco}: {count} sites")
    
    print(f"\n🌿 Biomes:")
    for biome, count in table_df['Biome'].value_counts().head(10).items():
        print(f"   {biome}: {count} sites")
    
    print(f"\n🌍 Climate Zones (shortened):")
    for climate, count in table_df['Climate'].value_counts().items():
        print(f"   {climate}: {count} sites")
    
    # Year range coverage
    valid_ranges = table_df[table_df['Data Range'] != 'N/A']
    print(f"\n📅 Year Range Coverage:")
    print(f"   Sites with valid year ranges: {len(valid_ranges)} / {len(table_df)}")
    
    # Show sites with discontinuous ranges
    discontinuous = table_df[table_df['Data Range'].str.contains(',', na=False)]
    if len(discontinuous) > 0:
        print(f"\n⚠️ Sites with discontinuous data (multiple ranges):")
        for _, row in discontinuous.iterrows():
            print(f"   {row['Ameriflux Site']}: {row['Data Range']}")
    
    print("\n" + "="*80)
    print("✅ TABLE 1 CREATION COMPLETE!")
    print("="*80)
    
    return table_df

# ============================================================================
# RUN
# ============================================================================

if __name__ == "__main__":
    table = main()