# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 17:35:46 2026

@author: ammar
"""

"""
CHUNK 1: Northeast AOI Disturbance Context

Purpose:
Create a clean Northeast disturbance master file and supporting summary tables/figures
from the original IDS geodatabase.

INPUTS:
- IDS Geodatabase: CONUS_Region9_AllYears.gdb (layer: DAMAGE_AREAS_FLAT_AllYears_CONUS_Rgn9)
- State boundaries: cb_2018_us_state_500k.shp
- AOI States: CT, ME, MA, NH, NJ, NY, PA, RI, VT

OUTPUTS:
- Master Files: NE_AOI_all_disturbances_master.gpkg, .csv
- Summary Tables: agent, damage type, agent-damage, year summaries
- Figures: 4 context figures (overview, agents, categories, states)
"""

import os
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import MultiPolygon
from shapely.validation import make_valid
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR = r'M:\Research\NE_temperate_forest_resilience\Data'
DISTURBANCE_GDB = os.path.join(BASE_DIR, r'Disturbance_data\CONUS_Region9_AllYears.gdb')
STATE_BOUNDARY = os.path.join(BASE_DIR, r'Boundaries\cb_2018_us_state_500k.shp')
OUTPUT_DIR = os.path.join(BASE_DIR, r'Disturbance_data\results')
FIGURE_DIR = os.path.join(OUTPUT_DIR, 'figures')

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# AOI States
AOI_STATES = ['CT', 'ME', 'MA', 'NH', 'NJ', 'NY', 'PA', 'RI', 'VT']
STATE_NAMES = {'CT': 'Connecticut', 'ME': 'Maine', 'MA': 'Massachusetts',
               'NH': 'New Hampshire', 'NJ': 'New Jersey', 'NY': 'New York',
               'PA': 'Pennsylvania', 'RI': 'Rhode Island', 'VT': 'Vermont'}

print("=" * 80)
print("CHUNK 1: Northeast AOI Disturbance Context")
print("=" * 80)

# ============================================================================
# STEP 1: Load and Prepare Data
# ============================================================================

print("\nStep 1: Loading data...")

# Read disturbance layer - SAFER VERSION
layers = gpd.list_layers(DISTURBANCE_GDB)

if isinstance(layers, pd.DataFrame):
    layer_names = layers["name"].tolist()
else:
    layer_names = [l.name if hasattr(l, "name") else str(l) for l in layers]

layer = next((l for l in layer_names if "DAMAGE_AREAS_FLAT" in l), None)
if layer is None:
    raise ValueError("Could not find DAMAGE_AREAS_FLAT layer")

gdf = gpd.read_file(DISTURBANCE_GDB, layer=layer)
print(f"  Loaded {len(gdf):,} records")

# Filter to MODIS era (>= 2000)
if 'SURVEY_YEAR' in gdf.columns:
    gdf = gdf[gdf['SURVEY_YEAR'] >= 2000].copy()
    print(f"  Filtered to >=2000: {len(gdf):,} records")

# Fix geometries
def fix_geometry(geom):
    if geom is None or geom.is_empty:
        return None
    if geom.is_valid:
        return geom
    try:
        valid = make_valid(geom)
        if valid and not valid.is_empty:
            return valid
    except:
        pass
    try:
        valid = geom.buffer(0)
        if valid and not valid.is_empty:
            return valid
    except:
        pass
    return None

gdf['geometry'] = gdf['geometry'].apply(fix_geometry)
gdf = gdf[gdf['geometry'].notna()].copy()
gdf = gdf[gdf['geometry'].apply(lambda x: x.geom_type in ['Polygon', 'MultiPolygon'])].copy()
print(f"  Valid polygons: {len(gdf):,}")

# Load and prepare state boundaries
states = gpd.read_file(STATE_BOUNDARY)
states = states[states['STUSPS'].isin(AOI_STATES)].copy()
states['AOI_state_abbr'] = states['STUSPS']
states['AOI_state'] = states['STUSPS'].map(STATE_NAMES)
states['geometry'] = states['geometry'].apply(fix_geometry)
states = states[states['geometry'].notna()].copy()

# Align CRS
if gdf.crs != states.crs:
    states = states.to_crs(gdf.crs)

# ============================================================================
# STEP 2: Clip to AOI
# ============================================================================

print("\nStep 2: Clipping to AOI states...")

# Spatial join to assign states
gdf = gpd.sjoin(gdf, states[['AOI_state', 'AOI_state_abbr', 'geometry']], 
                how='inner', predicate='intersects')
print(f"  After state assignment: {len(gdf):,}")

# Clip each state
clipped = []
for abbr in AOI_STATES:
    state_geom = states[states['AOI_state_abbr'] == abbr].geometry.iloc[0]
    state_gdf = gdf[gdf['AOI_state_abbr'] == abbr].copy()
    
    def clip_geom(geom):
        try:
            clipped = geom.intersection(state_geom)
            if clipped and not clipped.is_empty:
                if clipped.geom_type in ['Polygon', 'MultiPolygon']:
                    return clipped
                elif clipped.geom_type == 'GeometryCollection':
                    polys = [g for g in clipped.geoms if g.geom_type in ['Polygon', 'MultiPolygon']]
                    if polys:
                        return MultiPolygon(polys) if len(polys) > 1 else polys[0]
            return None
        except:
            return None
    
    state_gdf['geometry'] = state_gdf['geometry'].apply(clip_geom)
    state_gdf = state_gdf[state_gdf['geometry'].notna()].copy()
    if len(state_gdf) > 0:
        clipped.append(state_gdf)

gdf = gpd.GeoDataFrame(pd.concat(clipped, ignore_index=True), geometry='geometry', crs=gdf.crs)
print(f"  Clipped records: {len(gdf):,}")

# Remove sjoin helper column
if 'index_right' in gdf.columns:
    gdf = gdf.drop(columns=['index_right'])

# ============================================================================
# STEP 3: Calculate Areas and Add Fields
# ============================================================================

print("\nStep 3: Calculating areas and adding fields...")

# Project to EPSG:5070 for accurate area
gdf_proj = gdf.to_crs('EPSG:5070')

# Area fields
gdf_proj['area_m2'] = gdf_proj.geometry.area
gdf_proj['area_km2'] = gdf_proj['area_m2'] / 1_000_000
gdf_proj['area_ha'] = gdf_proj['area_m2'] / 10_000
gdf_proj['area_acres_clipped'] = gdf_proj['area_m2'] / 4046.8564224
gdf_proj['modis500_pixel_equiv'] = gdf_proj['area_km2'] / 0.25
gdf_proj['modis250_pixel_equiv'] = gdf_proj['area_km2'] / 0.0625

# Clean fields
gdf_proj['original_agent'] = gdf_proj.get('DCA_COMMON_NAME', gdf_proj.get('AGENT_NAME', ''))
gdf_proj['original_damage_type'] = gdf_proj.get('DAMAGE_TYPE', '')
gdf_proj['original_host'] = gdf_proj.get('HOST', '')
gdf_proj['original_host_group'] = gdf_proj.get('HOST_GROUP', '')

# Damage intensity
def parse_intensity(text):
    if pd.isna(text) or text == '':
        return np.nan
    t = str(text).lower().strip()
    if 'very severe' in t or '>50%' in t:
        return 75
    elif 'severe' in t and '30-50' in t:
        return 40
    elif 'severe' in t:
        return 40
    elif 'moderate' in t and '11-29' in t:
        return 20
    elif 'moderate' in t:
        return 20
    elif 'very light' in t or '1-3' in t:
        return 2
    elif 'light' in t and '4-10' in t:
        return 7
    elif 'light' in t:
        return 7
    elif 'trace' in t:
        return 1
    return np.nan

gdf_proj['damage_intensity_mid'] = gdf_proj.apply(
    lambda r: float(r['PERCENT_MID']) if pd.notna(r.get('PERCENT_MID')) else 
              parse_intensity(r.get('PERCENT_AFFECTED')) if pd.notna(r.get('PERCENT_AFFECTED')) else
              parse_intensity(r.get('LEGACY_SEVERITY')) if pd.notna(r.get('LEGACY_SEVERITY')) else np.nan,
    axis=1
)

gdf_proj['high_damage_50_flag'] = gdf_proj['damage_intensity_mid'] >= 50
gdf_proj['moderate_damage_30_flag'] = gdf_proj['damage_intensity_mid'] >= 30

# Suggested category
def suggest_category(agent):
    if pd.isna(agent) or agent == '':
        return 'unknown/no-data'
    t = str(agent).lower()
    if any(k in t for k in ['moth', 'caterpillar', 'budworm', 'leafminer', 'beetle', 'borer', 'adelgid', 'scale']):
        return 'insect'
    if any(k in t for k in ['disease', 'wilt', 'blight', 'canker', 'rust', 'fungus', 'needlecast']):
        return 'disease/pathogen'
    if any(k in t for k in ['decline', 'dieback', 'scorch']):
        return 'decline/complex'
    if any(k in t for k in ['drought', 'frost', 'flood', 'wind', 'hurricane', 'fire']):
        return 'abiotic/weather'
    if any(k in t for k in ['harvest', 'timber', 'logging', 'cut', 'development']):
        return 'anthropogenic/land-use'
    return 'other/unclassified'

gdf_proj['suggested_category'] = gdf_proj['original_agent'].apply(suggest_category)

print(f"  Total area: {gdf_proj['area_km2'].sum():.2f} km²")
print(f"  MODIS 500m pixels: {gdf_proj['modis500_pixel_equiv'].sum():.0f}")

# ============================================================================
# STEP 4: Save Master Files
# ============================================================================

print("\nStep 4: Saving master files...")

# Save projected
gdf_proj.to_file(os.path.join(OUTPUT_DIR, 'NE_AOI_all_disturbances_master.gpkg'), driver='GPKG')

# Save WGS84
gdf_wgs84 = gdf_proj.to_crs('EPSG:4326')
gdf_wgs84.to_file(os.path.join(OUTPUT_DIR, 'NE_AOI_all_disturbances_master_WGS84.gpkg'), driver='GPKG')

# Save CSV
gdf_wgs84.drop(columns=['geometry']).to_csv(
    os.path.join(OUTPUT_DIR, 'NE_AOI_all_disturbances_master.csv'), index=False
)

# Save AOI boundary
states.to_file(os.path.join(OUTPUT_DIR, 'AOI_Northeast_Census_states.gpkg'), driver='GPKG')
print("  Master files saved")

# ============================================================================
# STEP 5: Create Summary Tables
# ============================================================================

print("\nStep 5: Creating summary tables...")

df = gdf_proj

# Agent summary
agent_summary = df.groupby('original_agent').agg(
    record_count=('DAMAGE_AREA_ID', 'count'),
    total_area_km2=('area_km2', 'sum'),
    total_modis500_pixel_equiv=('modis500_pixel_equiv', 'sum'),
    high_damage_record_count=('high_damage_50_flag', 'sum'),
    states_present=('AOI_state_abbr', lambda x: '; '.join(sorted(x.dropna().unique())))
).reset_index()
agent_summary = agent_summary.sort_values('total_area_km2', ascending=False)
agent_summary.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_original_agent_summary.csv'), index=False)

# Damage type summary
damage_summary = df.groupby('original_damage_type').agg(
    record_count=('DAMAGE_AREA_ID', 'count'),
    total_area_km2=('area_km2', 'sum'),
    total_modis500_pixel_equiv=('modis500_pixel_equiv', 'sum'),
    high_damage_record_count=('high_damage_50_flag', 'sum')
).reset_index()
damage_summary = damage_summary.sort_values('total_area_km2', ascending=False)
damage_summary.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_damage_type_summary.csv'), index=False)

# Agent-damage summary
agent_damage = df.groupby(['original_agent', 'original_damage_type']).agg(
    record_count=('DAMAGE_AREA_ID', 'count'),
    total_area_km2=('area_km2', 'sum'),
    total_modis500_pixel_equiv=('modis500_pixel_equiv', 'sum'),
    high_damage_record_count=('high_damage_50_flag', 'sum'),
    states_present=('AOI_state_abbr', lambda x: '; '.join(sorted(x.dropna().unique())))
).reset_index()
agent_damage = agent_damage.sort_values('total_area_km2', ascending=False)
agent_damage.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_agent_damage_type_summary.csv'), index=False)

# Agent-year regional
agent_year = df.groupby(['SURVEY_YEAR', 'original_agent', 'original_damage_type']).agg(
    annual_area_km2=('area_km2', 'sum'),
    annual_modis500_pixel_equiv=('modis500_pixel_equiv', 'sum'),
    states_present=('AOI_state_abbr', lambda x: '; '.join(sorted(x.dropna().unique())))
).reset_index()
agent_year = agent_year.sort_values(['SURVEY_YEAR', 'annual_area_km2'], ascending=[True, False])
agent_year.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_agent_year_summary_regional.csv'), index=False)

# Agent-year by state
state_year = df.groupby(['AOI_state', 'AOI_state_abbr', 'SURVEY_YEAR', 'original_agent', 'original_damage_type']).agg(
    annual_area_km2=('area_km2', 'sum'),
    annual_modis500_pixel_equiv=('modis500_pixel_equiv', 'sum')
).reset_index()
state_year = state_year.sort_values(['AOI_state_abbr', 'SURVEY_YEAR', 'annual_area_km2'], ascending=[True, True, False])
state_year.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_agent_year_summary_by_state.csv'), index=False)

# Major disturbance years (>= 25 MODIS pixels)
major_years = agent_year[agent_year['annual_modis500_pixel_equiv'] >= 25].copy()
major_years = major_years.sort_values('annual_modis500_pixel_equiv', ascending=False)
major_years.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_major_disturbance_years.csv'), index=False)

# High damage summary
high_damage = df[df['high_damage_50_flag']].copy()
if len(high_damage) > 0:
    high_summary = high_damage.groupby(['original_agent', 'original_damage_type']).agg(
        high_damage_area_km2=('area_km2', 'sum'),
        high_damage_modis500_pixel_equiv=('modis500_pixel_equiv', 'sum'),
        mean_intensity=('damage_intensity_mid', 'mean')
    ).reset_index()
    high_summary = high_summary.sort_values('high_damage_area_km2', ascending=False)
    high_summary.to_csv(os.path.join(OUTPUT_DIR, 'NE_AOI_high_damage_summary.csv'), index=False)

print("  Summary tables saved")

# ============================================================================
# STEP 6: Create Figures
# ============================================================================

print("\nStep 6: Creating figures...")

# Set plot style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 9

# Color palette for categories
cat_colors = {
    'insect': '#2ecc71',
    'disease/pathogen': '#3498db',
    'decline/complex': '#9b59b6',
    'abiotic/weather': '#f39c12',
    'anthropogenic/land-use': '#e74c3c',
    'unknown/no-data': '#95a5a6',
    'other/unclassified': '#bdc3c7'
}
cat_order = ['insect', 'disease/pathogen', 'decline/complex', 'abiotic/weather', 
             'anthropogenic/land-use', 'unknown/no-data', 'other/unclassified']

# Add categories to agent_year for figures
agent_year_cat = agent_year.merge(
    df[['original_agent', 'suggested_category']].drop_duplicates('original_agent'),
    on='original_agent', how='left'
)
agent_year_cat['suggested_category'] = agent_year_cat['suggested_category'].fillna('other/unclassified')

# Figure 1: Overview
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))
fig1.suptitle('Northeast AOI Disturbance Overview', fontsize=16, fontweight='bold')

# A: Map
states.boundary.plot(ax=ax1, edgecolor='black', linewidth=0.5)
states.plot(ax=ax1, facecolor='lightgray', alpha=0.3)
for _, row in states.iterrows():
    ax1.text(row.geometry.centroid.x, row.geometry.centroid.y, 
             row['AOI_state_abbr'], ha='center', va='center', fontweight='bold')
ax1.set_title('A) Study Area')
ax1.axis('off')

# B: State totals
state_totals = state_year.groupby('AOI_state_abbr')['annual_modis500_pixel_equiv'].sum().sort_values()
state_totals.plot(kind='barh', ax=ax2, color='steelblue')
ax2.set_title('B) Mapped Disturbance by State')
ax2.set_xlabel('MODIS 500m pixels')
ax2.grid(axis='x', alpha=0.3)

# C: Top agents
top_agents = agent_summary.nlargest(15, 'total_modis500_pixel_equiv').sort_values('total_modis500_pixel_equiv')
ax3.barh(top_agents['original_agent'], top_agents['total_modis500_pixel_equiv']/1000, color='forestgreen')
ax3.set_title('C) Top 15 Disturbance Agents')
ax3.set_xlabel('MODIS 500m pixels (thousands)')
ax3.grid(axis='x', alpha=0.3)

# D: Annual time series
annual = agent_year.groupby('SURVEY_YEAR')['annual_modis500_pixel_equiv'].sum()
ax4.bar(annual.index, annual/1000, color='darkorange')
ax4.set_title('D) Annual Mapped Disturbance')
ax4.set_xlabel('Year')
ax4.set_ylabel('MODIS 500m pixels (thousands)')
ax4.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure1_disturbance_overview.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure1_disturbance_overview.pdf'), bbox_inches='tight')
plt.close()
print("  Figure 1 saved")

# Figure 2: Agents, damage types, events
fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))
fig2.suptitle('Northeast AOI Disturbance Identity', fontsize=16, fontweight='bold')

# A: Agent-damage combinations
top_combos = agent_damage.nlargest(15, 'total_modis500_pixel_equiv').sort_values('total_modis500_pixel_equiv')
labels = [f"{row['original_agent'][:20]}…{row['original_damage_type'][:15]}" 
          for _, row in top_combos.iterrows()]
ax1.barh(labels, top_combos['total_modis500_pixel_equiv']/1000, color='darkcyan')
ax1.set_title('A) Top Agent × Damage Type Combinations')
ax1.set_xlabel('MODIS 500m pixels (thousands)')
ax1.grid(axis='x', alpha=0.3)

# B: Damage types
top_damage = damage_summary.nlargest(15, 'total_modis500_pixel_equiv').sort_values('total_modis500_pixel_equiv')
ax2.barh(top_damage['original_damage_type'], top_damage['total_modis500_pixel_equiv']/1000, color='coral')
ax2.set_title('B) Top Damage Types')
ax2.set_xlabel('MODIS 500m pixels (thousands)')
ax2.grid(axis='x', alpha=0.3)

# C: Major events
top_events = major_years.nlargest(15, 'annual_modis500_pixel_equiv').sort_values('annual_modis500_pixel_equiv')
event_labels = [f"{row['SURVEY_YEAR']} - {row['original_agent'][:20]}" for _, row in top_events.iterrows()]
ax3.barh(event_labels, top_events['annual_modis500_pixel_equiv']/1000, color='mediumpurple')
ax3.set_title('C) Major Disturbance Events (>= 25 pixels)')
ax3.set_xlabel('MODIS 500m pixels (thousands)')
ax3.grid(axis='x', alpha=0.3)

# D: High damage
if len(high_damage) > 0:
    top_high = high_summary.nlargest(15, 'high_damage_modis500_pixel_equiv').sort_values('high_damage_modis500_pixel_equiv')
    high_labels = [f"{row['original_agent'][:20]}…{row['original_damage_type'][:15]}" for _, row in top_high.iterrows()]
    ax4.barh(high_labels, top_high['high_damage_modis500_pixel_equiv']/1000, color='crimson')
    ax4.set_title('D) High-Damage Events (>= 50% intensity)')
    ax4.set_xlabel('MODIS 500m pixels (thousands)')
    ax4.grid(axis='x', alpha=0.3)
else:
    ax4.text(0.5, 0.5, 'No high-damage records', transform=ax4.transAxes, ha='center', va='center')
    ax4.set_title('D) High-Damage Events')

plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure2_agents_damage_years.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure2_agents_damage_years.pdf'), bbox_inches='tight')
plt.close()
print("  Figure 2 saved")

# Figure 3: Category time series
fig3, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))
fig3.suptitle('Northeast AOI Disturbance Categories Through Time', fontsize=16, fontweight='bold')

# Pivot categories by year
cat_pivot = agent_year_cat.pivot_table(
    index='SURVEY_YEAR', columns='suggested_category', 
    values='annual_modis500_pixel_equiv', aggfunc='sum', fill_value=0
)
for cat in cat_order:
    if cat not in cat_pivot.columns:
        cat_pivot[cat] = 0
cat_pivot = cat_pivot[cat_order]

# A: Stacked bars - FIXED
(cat_pivot / 1000).plot(kind='bar', stacked=True, ax=ax1, color=[cat_colors[c] for c in cat_order])
ax1.set_title('A) Annual Disturbance by Category')
ax1.set_xlabel('Year')
ax1.set_ylabel('MODIS 500m pixels (thousands)')
ax1.legend(title='Category', bbox_to_anchor=(1.02, 1))
ax1.grid(axis='y', alpha=0.3)

# B: Total by category
cat_totals = cat_pivot.sum().sort_values()
ax2.barh(cat_totals.index, cat_totals/1000, color=[cat_colors[c] for c in cat_totals.index])
ax2.set_title('B) Total Disturbance by Category')
ax2.set_xlabel('MODIS 500m pixels (thousands)')
ax2.grid(axis='x', alpha=0.3)

# C: High damage by category
if len(high_damage) > 0:
    high_cat = high_damage.merge(
        df[['original_agent', 'suggested_category']].drop_duplicates('original_agent'),
        on='original_agent', how='left'
    )
    high_cat['suggested_category'] = high_cat['suggested_category'].fillna('other/unclassified')
    high_cat_totals = high_cat.groupby('suggested_category')['area_km2'].sum().reindex(cat_order, fill_value=0)
    ax3.barh(high_cat_totals.index, high_cat_totals, color=[cat_colors[c] for c in high_cat_totals.index])
    ax3.set_title('C) High-Damage Area by Category')
    ax3.set_xlabel('Area (km²)')
    ax3.grid(axis='x', alpha=0.3)
else:
    ax3.text(0.5, 0.5, 'No high-damage data', transform=ax3.transAxes, ha='center', va='center')
    ax3.set_title('C) High-Damage Area by Category')

# D: Category share over time
cat_pct = cat_pivot.div(cat_pivot.sum(axis=1), axis=0) * 100
cat_pct.plot(kind='bar', stacked=True, ax=ax4, color=[cat_colors[c] for c in cat_order])
ax4.set_title('D) Annual Category Share')
ax4.set_xlabel('Year')
ax4.set_ylabel('Percentage (%)')
ax4.legend(title='Category', bbox_to_anchor=(1.02, 1))
ax4.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure3_category_time_series.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure3_category_time_series.pdf'), bbox_inches='tight')
plt.close()
print("  Figure 3 saved")

# Figure 4: State patterns
fig4, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))
fig4.suptitle('Northeast AOI State-Level Disturbance Patterns', fontsize=16, fontweight='bold')

# A: State totals (same as Figure 1B but with names)
state_totals = state_year.groupby('AOI_state_abbr')['annual_modis500_pixel_equiv'].sum().sort_values()
state_labels = [f"{STATE_NAMES.get(s, s)} ({s})" for s in state_totals.index]
ax1.barh(state_labels, state_totals/1000, color='teal')
ax1.set_title('A) Total Mapped Disturbance by State')
ax1.set_xlabel('MODIS 500m pixels (thousands)')
ax1.grid(axis='x', alpha=0.3)

# B: Top agent by state
state_agent = state_year.groupby(['AOI_state_abbr', 'original_agent'])['annual_modis500_pixel_equiv'].sum().reset_index()
top_state_agent = state_agent.sort_values(['AOI_state_abbr', 'annual_modis500_pixel_equiv'], ascending=[True, False])
top_state_agent = top_state_agent.drop_duplicates('AOI_state_abbr').sort_values('annual_modis500_pixel_equiv')
state_labels2 = [f"{STATE_NAMES.get(s, s)} ({s})" for s in top_state_agent['AOI_state_abbr']]
ax2.barh(state_labels2, top_state_agent['annual_modis500_pixel_equiv']/1000, color='darkgreen')
ax2.set_title('B) Dominant Agent by State')
ax2.set_xlabel('MODIS 500m pixels (thousands)')
ax2.grid(axis='x', alpha=0.3)
for i, (_, row) in enumerate(top_state_agent.iterrows()):
    ax2.text(row['annual_modis500_pixel_equiv']/1000 + 5, i, row['original_agent'][:20], va='center', fontsize=7)

# C: State time series - FIXED
state_annual = state_year.pivot_table(
    index='SURVEY_YEAR', columns='AOI_state_abbr', 
    values='annual_modis500_pixel_equiv', aggfunc='sum', fill_value=0
)
state_annual = state_annual[state_annual.sum().sort_values(ascending=False).index[:8]]
(state_annual / 1000).plot(kind='bar', stacked=True, ax=ax3, colormap='viridis')
ax3.set_title('C) Annual Disturbance by State')
ax3.set_xlabel('Year')
ax3.set_ylabel('MODIS 500m pixels (thousands)')
ax3.legend(title='State', bbox_to_anchor=(1.02, 1))
ax3.grid(axis='y', alpha=0.3)

# D: State category composition - FIXED
state_cat = state_year.merge(
    df[['original_agent', 'suggested_category']].drop_duplicates('original_agent'),
    on='original_agent', how='left'
)
state_cat['suggested_category'] = state_cat['suggested_category'].fillna('other/unclassified')
state_cat_pivot = state_cat.pivot_table(
    index='AOI_state_abbr', columns='suggested_category',
    values='annual_modis500_pixel_equiv', aggfunc='sum', fill_value=0
)
for cat in cat_order:
    if cat not in state_cat_pivot.columns:
        state_cat_pivot[cat] = 0
state_cat_pivot = state_cat_pivot[cat_order]
state_cat_pivot = state_cat_pivot.loc[state_cat_pivot.sum(axis=1).sort_values(ascending=False).index]
(state_cat_pivot / 1000).plot(kind='barh', stacked=True, ax=ax4, color=[cat_colors[c] for c in cat_order])
ax4.set_title('D) State Composition by Category')
ax4.set_xlabel('MODIS 500m pixels (thousands)')
ax4.set_ylabel('State')
ax4.legend(title='Category', bbox_to_anchor=(1.02, 1))
ax4.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure4_state_patterns.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURE_DIR, 'NE_AOI_Figure4_state_patterns.pdf'), bbox_inches='tight')
plt.close()
print("  Figure 4 saved")

# ============================================================================
# STEP 7: Create Excel Workbook
# ============================================================================

print("\nStep 7: Creating Excel workbook...")

with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'NE_AOI_chunk1_disturbance_context_workbook.xlsx'), 
                    engine='openpyxl') as writer:
    agent_summary.to_excel(writer, sheet_name='Agent Summary', index=False)
    damage_summary.to_excel(writer, sheet_name='Damage Summary', index=False)
    agent_damage.to_excel(writer, sheet_name='Agent-Damage', index=False)
    agent_year.to_excel(writer, sheet_name='Agent-Year Regional', index=False)
    state_year.to_excel(writer, sheet_name='Agent-Year by State', index=False)
    major_years.to_excel(writer, sheet_name='Major Events', index=False)
    if len(high_damage) > 0:
        high_summary.to_excel(writer, sheet_name='High Damage', index=False)
    
    # Auto-adjust widths
    for sheet_name in writer.sheets:
        worksheet = writer.sheets[sheet_name]
        for column in worksheet.columns:
            max_len = max([len(str(cell.value)) for cell in column] + [len(column[0].column_letter)])
            worksheet.column_dimensions[column[0].column_letter].width = min(max_len + 2, 40)

print("  Workbook saved")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("CHUNK 1 COMPLETE")
print("=" * 80)
print(f"\nTotal records: {len(df):,}")
print(f"Total area: {df['area_km2'].sum():.2f} km²")
print(f"Year range: {df['SURVEY_YEAR'].min()} - {df['SURVEY_YEAR'].max()}")
print(f"\nOutput directory: {OUTPUT_DIR}")
print("=" * 80)