# -*- coding: utf-8 -*-
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 1. PATHS (Use exact paths)
# ============================================================================
MONTHLY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Removed PDF, added JPG
OUTPUT_JPG = os.path.join(OUTPUT_DIR, "coverage_figure.jpg")

# ============================================================================
# 2. EXACT COAST CLASSIFIER (From August figure)
# ============================================================================
EARTH_RADIUS_KM = 6371.0088
COAST_POLYLINES = {
    "Pacific Coast": [(32.6, -117.2), (33.6, -118.2), (34.4, -119.7), (35.4, -120.9), (36.6, -121.9), (37.8, -122.5), (39.0, -123.7), (41.0, -124.2), (43.5, -124.2), (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)],
    "Gulf Coast": [(25.8, -97.2), (28.0, -96.8), (29.3, -94.8), (29.5, -93.0), (29.2, -91.5), (29.3, -90.0), (29.5, -88.8), (30.2, -87.7), (30.1, -86.2), (29.9, -85.3), (29.7, -84.3), (28.8, -83.0), (27.8, -82.8), (26.6, -82.2), (25.9, -81.8), (25.3, -81.1), (25.0, -80.8)],
    "Atlantic Coast": [(25.1, -80.3), (26.0, -80.1), (27.0, -80.1), (28.4, -80.6), (29.0, -80.9), (29.9, -81.3), (30.4, -81.4), (31.2, -81.3), (32.0, -80.8), (33.0, -79.5), (34.5, -77.8), (35.7, -75.6), (36.8, -75.9), (38.5, -75.0), (39.5, -74.3), (40.5, -73.9), (41.3, -72.0), (42.4, -70.8), (43.5, -70.2)]
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
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
    return min(_seg_dist(lat, lon, *polyline[i], *polyline[i + 1]) for i in range(len(polyline) - 1))

def assign_coast_region(lat, lon):
    if pd.isna(lat) or pd.isna(lon):
        return np.nan
    lat = float(lat); lon = float(lon)
    if lat > 50:
        return "AK Coast"
    dists = {coast: _distance_to_polyline(lat, lon, polyline) for coast, polyline in COAST_POLYLINES.items()}
    return min(dists, key=dists.get)

# ============================================================================
# 3. LOAD DATA (Using exact rules from August study)
# ============================================================================
print(f"📂 Loading data from {MONTHLY_DATA_PATH} ...")
df = pd.read_csv(MONTHLY_DATA_PATH)

for col in df.select_dtypes(include="object").columns:
    df[col] = df[col].where(df[col].isna(), df[col].astype(str).str.strip())
    df[col] = df[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

required = ['site_name', 'Year', 'month', 'water_class', 'lat', 'long', 'Trans_ratio', 'WUE_tra']
for c in required:
    if c not in df.columns:
        raise ValueError(f"Column '{c}' not found")

df = df.dropna(subset=required).copy()
df = df[df['water_class'].isin(['Upland', 'Freshwater', 'Saline'])].copy()
df = df[np.isfinite(df['lat']) & np.isfinite(df['long']) & np.isfinite(df['Trans_ratio']) & (df['Trans_ratio'] >= 0) & (df['Trans_ratio'] <= 1) & np.isfinite(df['WUE_tra'])].copy()

print(f"  Filtered to {len(df)} monthly rows, {df['site_name'].nunique()} unique sites")

df['coast_region'] = df.apply(lambda row: assign_coast_region(row['lat'], row['long']), axis=1)
df = df.dropna(subset=['coast_region']).copy()

mapping = {'AK Coast': 'Alaska', 'Pacific Coast': 'Pacific', 'Gulf Coast': 'Gulf', 'Atlantic Coast': 'Atlantic'}
df['coast_final'] = df['coast_region'].map(mapping)

# ============================================================================
# 4. BUILD MATRIX AUTOMATICALLY
# ============================================================================
COAST_ORDER = ['Alaska', 'Pacific', 'Gulf', 'Atlantic']
COAST_COLORS = {
    'Alaska':  '#009E73',
    'Pacific': '#0072B2',
    'Gulf':    '#4D4D4D',
    'Atlantic':'#A50F15'
}
coast_val_map = {coast: i + 1 for i, coast in enumerate(COAST_ORDER)}

df_sites = df.groupby('site_name')['coast_final'].first().reset_index()
df_sites['coast_order'] = df_sites['coast_final'].map({c: i for i, c in enumerate(COAST_ORDER)})
df_sites = df_sites.sort_values(by=['coast_order', 'site_name'])

site_labels = df_sites['site_name'].tolist()
coast_labels = df_sites['coast_final'].tolist()
num_sites = df_sites['coast_final'].value_counts().to_dict()
all_years = sorted(df['Year'].unique())

data_matrix = []
for site in site_labels:
    years_present = set(df[df['site_name'] == site]['Year'].tolist())
    row = []
    for y in all_years:
        if y in years_present:
            row.append(coast_val_map[df_sites[df_sites['site_name'] == site]['coast_final'].iloc[0]])
        else:
            row.append(0)
    data_matrix.append(row)

matrix = np.array(data_matrix)

# ============================================================================
# 5. PLOT THE FIGURE (Updated formatting, larger fonts, JPG only)
# ============================================================================
print("Creating coverage figure ...")
fig, ax = plt.subplots(figsize=(28, 20)) 

colors_list = ['#E5E5E5'] + [COAST_COLORS[c] for c in COAST_ORDER]
cmap = ListedColormap(colors_list)

sns.heatmap(matrix, cmap=cmap, linewidths=1.5, linecolor='white', 
            xticklabels=all_years, yticklabels=site_labels,
            cbar=False, vmin=0, vmax=len(COAST_ORDER), ax=ax)

# --- Fix Separator Lines (Exactly at the boundary between rows) ---
cum_sums = np.cumsum([num_sites.get(coast, 0) for coast in COAST_ORDER])
for c in cum_sums[:-1]:
    ax.axhline(c, color='black', linewidth=3)

# --- Fix Tick Centers and Increase Font Sizes ---
# X-axis (Years): Centered + Font 26
ax.set_xticks(np.arange(len(all_years)) + 0.5)
ax.set_xticklabels(all_years, rotation=90, fontsize=26)
ax.set_xlabel('Year', fontsize=32, fontweight='bold', labelpad=15)

# Y-axis (Sites): Centered + Font 26 + Remove dot/line
ax.set_yticks(np.arange(len(site_labels)) + 0.5)
ax.set_yticklabels(site_labels, rotation=0, fontsize=24)
ax.tick_params(axis='y', length=0)  # Removes the dot/line next to site names

# Color site labels based on coast
for label, coast in zip(ax.get_yticklabels(), coast_labels):
    label.set_color(COAST_COLORS[coast])
    label.set_fontweight('bold')

# --- Reduce Space between Coast Names and Years ---
ax.set_xlim(-3, len(all_years))
ax.set_ylim(len(site_labels) - 0.5, -1.5)

# --- Top labels (Increased font size) ---
ax.text(-2.8, -1.0, "b)", fontsize=36, fontweight='bold', ha='left', va='center')
ax.text(-0.8, -1.0, f"N = {len(site_labels)}", fontsize=30, fontweight='bold', ha='left', va='center')

# --- Coast names (Increased font size) ---
centers = []
start = 0
for c in [num_sites.get(coast, 0) for coast in COAST_ORDER]:
    centers.append(start + (c / 2) - 0.5)
    start += c

for i, coast in enumerate(COAST_ORDER):
    ax.text(-2.8, centers[i], f"{coast}\n(N = {num_sites.get(coast, 0)})", 
            ha='left', va='center', fontsize=30, fontweight='bold', 
            color=COAST_COLORS[coast])

# Clean Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()

# Save JPG only (Higher DPI for better resolution)
fig.savefig(OUTPUT_JPG, dpi=900, bbox_inches='tight', facecolor='white')
print(f"✅ Saved JPG:  {OUTPUT_JPG}")

plt.show()