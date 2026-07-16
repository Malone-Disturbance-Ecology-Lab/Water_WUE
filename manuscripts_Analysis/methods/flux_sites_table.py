# -*- coding: utf-8 -*-
"""
SITE TABLE – 64 SITES, SHORT CLIMATE (SAMPLE STYLE), FULL BIOME NAMES
- Exact Study_area filtering (no SPEI filter)
- Climate shortened to readable format
- Biome: IGBP codes mapped to full names
- Output columns match the provided sample image
"""

import os
import pandas as pd
import numpy as np

# =============================================================================
# CONFIGURATION
# =============================================================================

MONTHLY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\tables"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "Site_Table_64sites.csv")

ECOSYSTEM_TYPES = ['Freshwater', 'Saline', 'Upland']

# =============================================================================
# CLIMATE SHORTENING (SAMPLE STYLE)
# =============================================================================

def shorten_climate(climate_str):
    """Convert verbose climate descriptions to a compact, readable format."""
    if pd.isna(climate_str):
        return ""
    s = str(climate_str).strip()
    
    patterns = [
        ("Humid Subtropical: mild with no dry season, hot summer", "Humid Subtropical (mild)"),
        ("Humid Subtropical: dry winter, hot summer", "Humid Subtropical (dry winter)"),
        ("Mediterranean: mild with dry, hot summer", "Mediterranean (hot summer)"),
        ("Mediterranean: mild with dry, warm summer", "Mediterranean (warm summer)"),
        ("Dry Continental: cool summer", "Dry Continental (cool summer)"),
        ("Humid Continental: humid with severe winter, no dry season, hot summer", "Humid Continental (hot summer)"),
        ("Subarctic: severe winter, no dry season, cool summer", "Subarctic"),
        ("Marine West Coast: mild with no dry season, cool summer", "Marine West Coast"),
        ("Tundra", "Tundra"),
    ]
    
    for pattern, short in patterns:
        if pattern in s:
            return short
    
    return s if len(s) <= 30 else s[:27] + "..."

# =============================================================================
# BIOME CODE MAPPING (IGBP → FULL NAME)
# =============================================================================

BIOME_MAP = {
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
    'URB': 'Urban and Built-up',
    'SNO': 'Snow and Ice',
    'BSV': 'Barren or Sparse Vegetation',
    # Add any other codes that appear
}

def full_biome_name(code):
    """Convert IGBP code to full biome name; return code if not found."""
    if pd.isna(code):
        return ""
    code = str(code).strip().upper()
    return BIOME_MAP.get(code, code)  # return code itself if not mapped

# =============================================================================
# DATA LOADING & FILTERING (EXACT STUDY_AREA)
# =============================================================================

def load_and_filter_data(filepath):
    print("\n" + "=" * 80)
    print("LOADING MONTHLY DATA – EXACT STUDY AREA FILTERS")
    print("=" * 80)
    
    df = pd.read_csv(filepath)
    print(f"✓ Loaded {len(df)} rows")
    
    # Clean strings
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].where(df[col].isna(), df[col].astype(str).str.strip())
        df[col] = df[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})
    
    required = ['site_name', 'Year', 'month', 'water_class', 'lat', 'long',
                'Trans_ratio', 'WUE_tra']
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found.")
    
    df = df.dropna(subset=required).copy()
    df = df[df['water_class'].isin(ECOSYSTEM_TYPES)].copy()
    df = df[np.isfinite(df['lat']) & np.isfinite(df['long']) &
            np.isfinite(df['Trans_ratio']) & (df['Trans_ratio'] >= 0) & (df['Trans_ratio'] <= 1) &
            np.isfinite(df['WUE_tra'])].copy()
    
    print(f"✓ After all filters: {len(df)} rows, {df['site_name'].nunique()} unique sites")
    return df

# =============================================================================
# SITE SUMMARY
# =============================================================================

def compute_site_summary(df):
    sites = df['site_name'].unique()
    rows = []
    
    # Identify climate and biome columns
    climate_col = 'climate' if 'climate' in df.columns else None
    biome_col = 'IGBP' if 'IGBP' in df.columns else None
    
    if climate_col:
        print(f"Using climate column: {climate_col}")
    else:
        print("No climate column – leaving blank.")
    if biome_col:
        print(f"Using biome column: {biome_col} (codes will be mapped to full names)")
    else:
        print("No biome column – leaving blank.")
    
    for site in sites:
        site_df = df[df['site_name'] == site]
        lat = site_df['lat'].iloc[0]
        lon = site_df['long'].iloc[0]
        eco = site_df['water_class'].iloc[0]
        
        # Climate (shortened)
        raw_climate = site_df[climate_col].dropna().iloc[0] if climate_col and climate_col in site_df.columns else np.nan
        climate_short = shorten_climate(raw_climate) if not pd.isna(raw_climate) else ""
        
        # Biome (full name from IGBP code)
        raw_biome = site_df[biome_col].dropna().iloc[0] if biome_col and biome_col in site_df.columns else np.nan
        biome_full = full_biome_name(raw_biome) if not pd.isna(raw_biome) else ""
        
        min_year = site_df['Year'].min()
        max_year = site_df['Year'].max()
        data_range = f"{min_year}-{max_year}"
        
        te_min = site_df['Trans_ratio'].min()
        te_max = site_df['Trans_ratio'].max()
        te_range = f"{te_min:.2f} - {te_max:.2f}"
        
        rows.append({
            'Ameriflux Site': site,
            'Latitude': lat,
            'Longitude': lon,
            'Climate': climate_short,
            'Biome': biome_full,
            'Ecosystem Type': eco,
            'Data Range': data_range,
            'T:ET Range': te_range
        })
    
    return pd.DataFrame(rows)

# =============================================================================
# OUTPUT
# =============================================================================

def print_table(df):
    print("\n" + "=" * 140)
    print(f"SITE TABLE – {len(df)} sites (matching study area map)")
    print("=" * 140)
    print(df.to_string(index=False))

def save_table(df, output_path):
    df.to_csv(output_path, index=False, encoding='utf-8')
    print(f"\n✅ Site table saved to: {output_path}")

# =============================================================================
# MAIN
# =============================================================================

def main():
    df_monthly = load_and_filter_data(MONTHLY_DATA_PATH)
    if df_monthly.empty:
        print("❌ No data after filtering.")
        return
    
    site_df = compute_site_summary(df_monthly)
    print(f"\nNumber of sites: {len(site_df)}")
    print_table(site_df)
    save_table(site_df, OUTPUT_CSV)
    
    print("\n" + "=" * 80)
    print("SITE TABLE GENERATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()