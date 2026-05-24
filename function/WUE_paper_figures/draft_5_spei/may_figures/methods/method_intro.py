# -*- coding: utf-8 -*-
"""
COMPLETE METHODS TABLE with T:ET Ratio
Includes T:ET by:
- Coastal Region (with IQR)
- Ecosystem Type / Salinity Class (Upland, Freshwater, Saline, Brackish) with IQR
- IGBP Biome (WET, GRA, ENF, CSH, etc.) with IQR
"""

import pandas as pd
import numpy as np
import os

results_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"

# Load the final methods dataset
input_file = os.path.join(results_dir, "methods_site_summary_COMPLETE.csv")
df = pd.read_csv(input_file)

print("="*80)
print("COMPLETE METHODS TABLE with T:ET Ratio BY ECOSYSTEM TYPE AND REGION")
print("="*80)

# =============================================================================
# 1. VERIFY REGION COUNTS
# =============================================================================

print("\n" + "="*80)
print("1. REGION COUNTS VERIFICATION")
print("="*80)

region_counts = df['coast_region'].value_counts()
print(f"\n  Pacific Coast:          n = {region_counts.get('Pacific', 0)}")
print(f"  Atlantic North Coast:   n = {region_counts.get('Atlantic North', 0)}")
print(f"  Alaska Coast:           n = {region_counts.get('Alaska', 0)}")
print(f"  Gulf of Mexico:         n = {region_counts.get('Gulf of Mexico', 0)}")
print(f"  Southeast Atlantic:     n = {region_counts.get('Southeast Atlantic', 0)}")
print(f"\n  TOTAL:                  n = {len(df)}")

# =============================================================================
# 2. T:ET BY COASTAL REGION (with IQR)
# =============================================================================

print("\n" + "="*80)
print("2. T:ET RATIO BY COASTAL REGION (with IQR)")
print("="*80)

def calc_stats_with_iqr(data):
    """Calculate statistics including IQR"""
    return pd.Series({
        'n_sites': len(data),
        'median': data.median(),
        'mean': data.mean(),
        'std': data.std(),
        'min': data.min(),
        'max': data.max(),
        'Q1 (25th)': data.quantile(0.25),
        'Q3 (75th)': data.quantile(0.75),
        'IQR': data.quantile(0.75) - data.quantile(0.25)
    })

# Region T:ET statistics
region_tet_stats = df.groupby('coast_region')['Trans_ratio_median'].apply(calc_stats_with_iqr).unstack()

# Reorder regions
region_order = ['Pacific', 'Atlantic North', 'Alaska', 'Gulf of Mexico', 'Southeast Atlantic']
region_tet_stats = region_tet_stats.reindex([r for r in region_order if r in region_tet_stats.index])

print("\n  T:ET Ratio by Coastal Region (Median and IQR):")
print("  " + "-" * 90)
print(f"  {'Region':<20} {'N':<6} {'Median':<10} {'IQR':<10} {'Q1':<10} {'Q3':<10} {'Min':<8} {'Max':<8}")
print("  " + "-" * 90)

for region, row in region_tet_stats.iterrows():
    print(f"  {region:<20} {int(row['n_sites']):<6} {row['median']:<10.3f} {row['IQR']:<10.3f} {row['Q1 (25th)']:<10.3f} {row['Q3 (75th)']:<10.3f} {row['min']:<8.3f} {row['max']:<8.3f}")

print("  " + "-" * 90)

# =============================================================================
# 3. T:ET BY ECOSYSTEM TYPE (Salinity Class) with IQR
# =============================================================================

print("\n" + "="*80)
print("3. T:ET RATIO BY ECOSYSTEM TYPE (Upland, Freshwater, Saline, Brackish)")
print("="*80)

ecosystem_tet_stats = df.groupby('Salinity_Category')['Trans_ratio_median'].apply(calc_stats_with_iqr).unstack()

# Reorder ecosystem types
ecosystem_order = ['Upland', 'Freshwater', 'Saline', 'Brackish']
ecosystem_tet_stats = ecosystem_tet_stats.reindex([e for e in ecosystem_order if e in ecosystem_tet_stats.index])

print("\n  T:ET Ratio by Ecosystem Type (Median and IQR):")
print("  " + "-" * 90)
print(f"  {'Ecosystem Type':<20} {'N':<6} {'Median':<10} {'IQR':<10} {'Q1':<10} {'Q3':<10} {'Min':<8} {'Max':<8}")
print("  " + "-" * 90)

for eco, row in ecosystem_tet_stats.iterrows():
    print(f"  {eco:<20} {int(row['n_sites']):<6} {row['median']:<10.3f} {row['IQR']:<10.3f} {row['Q1 (25th)']:<10.3f} {row['Q3 (75th)']:<10.3f} {row['min']:<8.3f} {row['max']:<8.3f}")

print("  " + "-" * 90)

# =============================================================================
# 4. T:ET BY IGBP BIOME (with IQR)
# =============================================================================

print("\n" + "="*80)
print("4. T:ET RATIO BY IGBP BIOME")
print("="*80)

# IGBP display names
igbp_display = {
    'WET': 'Wetlands',
    'GRA': 'Grasslands',
    'ENF': 'Evergreen Needleleaf Forest',
    'CSH': 'Closed Shrublands',
    'CRO': 'Croplands',
    'DBF': 'Deciduous Broadleaf Forest',
    'MF': 'Mixed Forest',
    'OSH': 'Open Shrublands',
    'BSV': 'Barren/Sparse Vegetation'
}

# Add display name column
df['IGBP_display'] = df['IGBP'].map(igbp_display)

igbp_tet_stats = df.groupby('IGBP_display')['Trans_ratio_median'].apply(calc_stats_with_iqr).unstack()

# Sort by median T:ET
igbp_tet_stats = igbp_tet_stats.sort_values('median', ascending=False)

print("\n  T:ET Ratio by IGBP Biome (Median and IQR):")
print("  " + "-" * 90)
print(f"  {'IGBP Biome':<30} {'N':<6} {'Median':<10} {'IQR':<10} {'Q1':<10} {'Q3':<10} {'Min':<8} {'Max':<8}")
print("  " + "-" * 90)

for biome, row in igbp_tet_stats.iterrows():
    print(f"  {biome:<30} {int(row['n_sites']):<6} {row['median']:<10.3f} {row['IQR']:<10.3f} {row['Q1 (25th)']:<10.3f} {row['Q3 (75th)']:<10.3f} {row['min']:<8.3f} {row['max']:<8.3f}")

print("  " + "-" * 90)

# =============================================================================
# 5. CREATE FORMATTED SITE TABLE
# =============================================================================

print("\n" + "="*80)
print("5. FORMATTED SITE TABLE (Alphabetical by Site Name)")
print("="*80)

# Get data range per site from monthly data
monthly_file = os.path.join(results_dir, "monthly_data_after_outlier_removal.csv")
df_monthly = pd.read_csv(monthly_file)

data_range = {}
for site in df['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site]
    if len(site_data) > 0:
        min_year = int(site_data['year'].min())
        max_year = int(site_data['year'].max())
        data_range[site] = f"{min_year}-{max_year}"
    else:
        data_range[site] = "N/A"

# Simplify climate names
def simplify_climate(climate):
    if pd.isna(climate):
        return "Unknown"
    climate_lower = climate.lower()
    if 'tundra' in climate_lower:
        return "Polar / Tundra"
    elif 'mediterranean' in climate_lower:
        return "Mediterranean"
    elif 'humid subtropical' in climate_lower:
        return "Humid Subtropical"
    elif 'humid continental' in climate_lower:
        return "Humid Continental"
    elif 'marine west coast' in climate_lower:
        return "Marine West Coast"
    elif 'dry continental' in climate_lower:
        return "Dry Continental"
    else:
        return climate[:30]

# Create table
table_rows = []
for _, row in df.sort_values('site_name').iterrows():
    table_rows.append({
        'Ameriflux Site': row['site_name'],
        'Latitude': f"{row['lat']:.2f}",
        'Longitude': f"{row['long']:.2f}",
        'Climate': simplify_climate(row['climate']),
        'Biome': igbp_display.get(row['IGBP'], row['IGBP']),
        'Ecosystem Type': row['Salinity_Category'],
        'Data range': data_range.get(row['site_name'], "N/A"),
        'T:ET Ratio': f"{row['Trans_ratio_median']:.3f}"
    })

df_table = pd.DataFrame(table_rows)

# Print header
print("\n")
header = f"{'Site':<12} {'Lat':<8} {'Lon':<9} {'Climate':<20} {'Biome':<28} {'Eco Type':<12} {'Data Range':<12} {'T:ET':<6}"
print(header)
print("-" * 120)

for _, row in df_table.iterrows():
    print(f"{row['Ameriflux Site']:<12} {row['Latitude']:<8} {row['Longitude']:<9} {row['Climate']:<20} {row['Biome']:<28} {row['Ecosystem Type']:<12} {row['Data range']:<12} {row['T:ET Ratio']:<6}")

print("\n" + "="*80)

# =============================================================================
# 6. SAVE TABLES TO CSV
# =============================================================================

# Save site table
output_table = os.path.join(results_dir, "methods_site_table_with_TET.csv")
df_table.to_csv(output_table, index=False)
print(f"\n✅ Saved site table to: {output_table}")

# Save region T:ET stats
region_stats_output = os.path.join(results_dir, "TET_by_region.csv")
region_tet_stats.to_csv(region_stats_output)
print(f"✅ Saved T:ET by region to: {region_stats_output}")

# Save ecosystem T:ET stats
ecosystem_stats_output = os.path.join(results_dir, "TET_by_ecosystem_type.csv")
ecosystem_tet_stats.to_csv(ecosystem_stats_output)
print(f"✅ Saved T:ET by ecosystem type to: {ecosystem_stats_output}")

# Save IGBP T:ET stats
igbp_stats_output = os.path.join(results_dir, "TET_by_IGBP.csv")
igbp_tet_stats.to_csv(igbp_stats_output)
print(f"✅ Saved T:ET by IGBP to: {igbp_stats_output}")

# =============================================================================
# 7. PRINT SUMMARY STATISTICS FOR METHODS SECTION
# =============================================================================

print("\n" + "="*80)
print("6. SUMMARY STATISTICS FOR METHODS SECTION")
print("="*80)

tet_min = df['Trans_ratio_median'].min()
tet_max = df['Trans_ratio_median'].max()
tet_median = df['Trans_ratio_median'].median()
tet_iqr = df['Trans_ratio_median'].quantile(0.75) - df['Trans_ratio_median'].quantile(0.25)

print(f"""
  T:ET RATIO (Transpiration / Evapotranspiration):
  -------------------------------------------------
  Network-wide median:     {tet_median:.3f}
  Network-wide IQR:        {tet_iqr:.3f}
  Range across sites:      {tet_min:.3f} to {tet_max:.3f}
  
  T:ET BY COASTAL REGION (median [IQR]):
  -------------------------------------------------""")

for region, row in region_tet_stats.iterrows():
    print(f"    {region}:           {row['median']:.3f} [{row['IQR']:.3f}] (n={int(row['n_sites'])})")

print(f"""
  T:ET BY ECOSYSTEM TYPE (median [IQR]):
  -------------------------------------------------""")

for eco, row in ecosystem_tet_stats.iterrows():
    print(f"    {eco}:        {row['median']:.3f} [{row['IQR']:.3f}] (n={int(row['n_sites'])})")

print(f"""
  T:ET BY IGBP BIOME (median [IQR]):
  -------------------------------------------------""")

for biome, row in igbp_tet_stats.head(5).iterrows():
    print(f"    {biome}:           {row['median']:.3f} [{row['IQR']:.3f}] (n={int(row['n_sites'])})")

# =============================================================================
# 8. FINAL METHODS SECTION TEXT
# =============================================================================

print("\n" + "="*80)
print("7. METHODS SECTION TEXT — COPY THIS DIRECTLY")
print("="*80)

print(f"""
We obtained eddy covariance (EC) data from {len(df)} towers within the 
U.S. coastal zone (Fig. 1) across a broad latitudinal and climatic gradient, from 
high-latitude tundra to warm-temperate and subtropical coastal zones in the United 
States. The site network encompasses sites along the Alaska (n = {region_counts.get('Alaska', 0)}), 
Pacific (n = {region_counts.get('Pacific', 0)}), Gulf of Mexico (n = {region_counts.get('Gulf of Mexico', 0)}), 
Southeast Atlantic (n = {region_counts.get('Southeast Atlantic', 0)}), and Atlantic North 
(n = {region_counts.get('Atlantic North', 0)}) coasts. 

The network includes a diverse set of vegetation types based on the International 
Geosphere-Biosphere Programme (IGBP) classification, including wetlands 
(WET, n = 33), grasslands (GRA, n = 6), evergreen needleleaf forests (ENF, n = 5), 
closed shrublands (CSH, n = 5), croplands (CRO, n = 2), deciduous broadleaf forests 
(DBF, n = 2), mixed forests (MF, n = 1), open shrublands (OSH, n = 1), and barren 
or sparsely vegetated systems (BSV, n = 2) (Table 1). 

Sites were distributed across polar (tundra), Mediterranean-type, and humid temperate 
to subtropical climates. Sites were classified as upland (n = 23), freshwater (n = 16), 
saline (n = 13), or brackish (n = 5) based on hydrological and salinity characteristics.

The median T:ET ratio (transpiration ÷ evapotranspiration) across all sites was 
{tet_median:.3f} (IQR = {tet_iqr:.3f}), ranging from {tet_min:.3f} to {tet_max:.3f}. 
Regional median T:ET values ranged from {region_tet_stats['median'].min():.3f} to {region_tet_stats['median'].max():.3f}, 
with the highest values in {region_tet_stats['median'].idxmax()} ({region_tet_stats.loc[region_tet_stats['median'].idxmax(), 'median']:.3f}) 
and the lowest in {region_tet_stats['median'].idxmin()} ({region_tet_stats.loc[region_tet_stats['median'].idxmin(), 'median']:.3f}). 
By ecosystem type, median T:ET was highest in {ecosystem_tet_stats['median'].idxmax()} ecosystems ({ecosystem_tet_stats.loc[ecosystem_tet_stats['median'].idxmax(), 'median']:.3f}) 
and lowest in {ecosystem_tet_stats['median'].idxmin()} ({ecosystem_tet_stats.loc[ecosystem_tet_stats['median'].idxmin(), 'median']:.3f}).
""")

print("\n" + "="*80)
print("✅ COMPLETE — All tables and summaries created")
print("="*80)

print(f"\n📁 Files saved:")
print(f"   1. {output_table}")
print(f"   2. {region_stats_output}")
print(f"   3. {ecosystem_stats_output}")
print(f"   4. {igbp_stats_output}")