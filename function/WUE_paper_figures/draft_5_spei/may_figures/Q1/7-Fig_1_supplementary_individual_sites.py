# -*- coding: utf-8 -*-
"""
SUPPLEMENTARY FIGURE: Site-level Boxplots by Ecosystem Type
USING SAME STRICT TRIPLE INTERSECTION AS FIGURE 1
- LOADS pre-computed strict triple intersection from Figure 1
- Uses same wue_site_level_NN_medians_SPEI_1.csv for site selection
- Applies SAME strict metric filter as CHUNK 3 (WUE, WUE_eva, WUE_tra ALL have data)
- Panel A: T:ET ratio
- Panel B: WUE_ET
- Panel C: WUE_T
- X-axis site names ONLY on Panel C (bottom)
- Sites ordered by T:ET ratio within each ecosystem group
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input files
MONTHLY_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
NN_MEDIANS_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"
METADATA_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\site_metadata_with_salinity_SPEIinfo.csv"

# Output directory
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\figures\supplementary"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_PNG = os.path.join(OUTPUT_DIR, "Figure_S1_StrictMonths_ThreePanel_Boxplots.png")
OUTPUT_PDF = os.path.join(OUTPUT_DIR, "Figure_S1_StrictMonths_ThreePanel_Boxplots.pdf")

# =============================================================================
# FIGURE DIMENSIONS
# =============================================================================
FIGURE_WIDTH = 55
FIGURE_HEIGHT = 32
DPI = 300

# =============================================================================
# FONT SIZES
# =============================================================================
PANEL_LABEL_FONT = 80
AXIS_LABEL_FONT = 65
Y_TICK_FONT = 75
X_TICK_FONT = 56
ECOSYSTEM_LABEL_FONT = 65
LEGEND_FONT = 50

# =============================================================================
# Y-AXIS SETTINGS
# =============================================================================
Y_LIMITS = {
    'Trans_ratio': (0, 1.0),
    'WUE': (0, 15),
    'WUE_tra': (0, 15)
}

Y_TICK_INTERVALS = {
    'Trans_ratio': 0.2,
    'WUE': 3,
    'WUE_tra': 3
}

# Y-LABELS WITH TWO ROWS TO PREVENT OVERLAP BETWEEN PANELS
Y_LABELS = {
    'Trans_ratio': 'T:ET ratio',
    'WUE': 'WUE$_{ET}$\n(g C kg$^{-1}$ H$_2$O$^{-1}$)',
    'WUE_tra': 'WUE$_T$\n(g C kg$^{-1}$ H$_2$O$^{-1}$)'
}

# =============================================================================
# BOXPLOT STYLING
# =============================================================================
BOX_WIDTH = 0.7
BOX_ALPHA = 0.7
BOX_LINEWIDTH = 3
MEDIAN_LINEWIDTH = 4
WHISKER_WIDTH = 2.5

# =============================================================================
# COLORS
# =============================================================================
ECOSYSTEM_COLORS = {
    'Upland': '#800080',
    'Freshwater': '#0000FF',
    'Brackish': '#008080',
    'Saline': '#FFA500'
}

ECOSYSTEM_ORDER = ['Upland', 'Freshwater', 'Brackish', 'Saline']

# =============================================================================
# PANEL ORDER (top to bottom)
# =============================================================================
PANELS = [
    {'col': 'Trans_ratio', 'label': 'T:ET ratio'},
    {'col': 'WUE', 'label': 'WUE_ET'},
    {'col': 'WUE_tra', 'label': 'WUE_T'}
]

# =============================================================================
# STEP 1: GET STRICT TRIPLE INTERSECTION SITES FROM FIGURE 1
# =============================================================================

print("="*80)
print("SUPPLEMENTARY FIGURE: USING SAME STRICT TRIPLE INTERSECTION AS FIGURE 1")
print("="*80)

# Load the NN medians file from Figure 1
nn_medians = pd.read_csv(NN_MEDIANS_FILE)
print(f"Loaded NN medians file: {len(nn_medians)} rows")

# Get sites with WUE available
wue_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())

# Get sites with WUE_eva available
eva_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())

# Get sites with WUE_tra available
tra_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())

# STRICT TRIPLE INTERSECTION: sites with ALL THREE metrics (SAME AS FIGURE 1)
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

print(f"\nStrict triple intersection sites (from Figure 1): {len(shared_sites)}")

# =============================================================================
# STEP 2: LOAD FILTERED MONTHLY DATA AND APPLY STRICT METRIC FILTER
# =============================================================================

print("\n" + "="*60)
print("LOADING FILTERED MONTHLY DATA (upstream CHUNK 3 already applied)")
print("="*60)

df_monthly = pd.read_csv(MONTHLY_FILE)
print(f"Loaded monthly data: {len(df_monthly)} rows")

# Filter to NN conditions using SPEI_1_Cat (already in file)
if 'SPEI_1_Cat' in df_monthly.columns:
    df_nn = df_monthly[df_monthly['SPEI_1_Cat'] == 'NN'].copy()
else:
    # Fallback if column doesn't exist
    def classify_spei(spei_value):
        if pd.isna(spei_value):
            return np.nan
        if spei_value >= 2.0:
            return "EW"
        elif 1.5 <= spei_value < 2.0:
            return "SW"
        elif 1.0 <= spei_value < 1.5:
            return "MoW"
        elif 0.5 < spei_value < 1.0:
            return "MW"
        elif -0.5 <= spei_value <= 0.5:
            return "NN"
        elif -1.0 < spei_value < -0.5:
            return "MD"
        elif -1.5 < spei_value <= -1.0:
            return "MoD"
        elif -2.0 < spei_value <= -1.5:
            return "SD"
        elif spei_value <= -2.0:
            return "ED"
        else:
            return np.nan
    df_nn = df_monthly.copy()
    df_nn['SPEI_1_Cat'] = df_nn['SPEI_1'].apply(classify_spei)
    df_nn = df_nn[df_nn['SPEI_1_Cat'] == 'NN'].copy()

print(f"NN rows: {len(df_nn)}")

# Filter to strict triple intersection sites ONLY
df_nn = df_nn[df_nn['site_name'].isin(shared_sites)].copy()
print(f"NN rows after filtering to strict triple intersection: {len(df_nn)}")
print(f"Unique sites: {df_nn['site_name'].nunique()} (expected: {len(shared_sites)})")

# =============================================================================
# CRITICAL FIX: Apply SAME strict metric filter as CHUNK 3
# Only keep rows where WUE, WUE_eva, AND WUE_tra ALL have data
# =============================================================================
print("\n" + "="*60)
print("APPLYING STRICT METRIC FILTER (SAME AS CHUNK 3)")
print("="*60)

strict_mask = (
    df_nn['WUE'].notna() &
    df_nn['WUE_eva'].notna() &
    df_nn['WUE_tra'].notna()
)

df_nn = df_nn[strict_mask].copy()
print(f"NN rows after strict metric filter: {len(df_nn):,}")
print(f"Rows removed: {len(df_nn) - len(df_nn) if 'df_nn_before' not in dir() else 'see above'}")
print(f"Unique sites after strict metric filter: {df_nn['site_name'].nunique()}")

# =============================================================================
# STEP 3: LOAD METADATA AND ADD ECOSYSTEM TYPE
# =============================================================================

print("\n" + "="*60)
print("LOADING METADATA")
print("="*60)

df_meta = pd.read_csv(METADATA_FILE)
site_ecosystem = df_meta[['site_name', 'Salinity_Category']].drop_duplicates()

# Add ecosystem to nn data
df_nn = df_nn.merge(site_ecosystem, on='site_name', how='left')

# =============================================================================
# STEP 4: CALCULATE T:ET MEDIANS FROM STRICT MONTHS
# =============================================================================

print("\n" + "="*60)
print("CALCULATING T:ET MEDIANS")
print("="*60)

tet_medians = df_nn.groupby('site_name')['Trans_ratio'].median().to_dict()

print(f"\nT:ET MEDIAN RANGE FOR {len(tet_medians)} SITES:")
tet_values_list = list(tet_medians.values())
print(f"  Min T:ET: {min(tet_values_list):.4f}")
print(f"  Max T:ET: {max(tet_values_list):.4f}")
print(f"  Mean T:ET: {np.mean(tet_values_list):.4f}")
print(f"  Median T:ET: {np.median(tet_values_list):.4f}")

# Low T:ET sites (<0.1)
low_tet_sites = [(site, med) for site, med in tet_medians.items() if med < 0.1]
low_tet_sites.sort(key=lambda x: x[1])
print(f"\nSites with T:ET < 0.1: {len(low_tet_sites)}")
for site, med in low_tet_sites:
    eco = site_ecosystem[site_ecosystem['site_name'] == site]['Salinity_Category'].iloc[0]
    print(f"  {site}: {eco}, T:ET = {med:.4f}")

# High T:ET sites (>0.7)
high_tet_sites = [(site, med) for site, med in tet_medians.items() if med > 0.7]
high_tet_sites.sort(key=lambda x: x[1], reverse=True)
print(f"\nSites with T:ET > 0.7: {len(high_tet_sites)}")
for site, med in high_tet_sites:
    eco = site_ecosystem[site_ecosystem['site_name'] == site]['Salinity_Category'].iloc[0]
    print(f"  {site}: {eco}, T:ET = {med:.4f}")

# =============================================================================
# STEP 5: ORDER SITES BY T:ET RATIO WITHIN ECOSYSTEM GROUP
# =============================================================================

print("\n" + "="*60)
print("ORDERING SITES BY T:ET RATIO WITHIN ECOSYSTEM")
print("="*60)

ordered_sites = []
ordered_ecosystems = []

for eco in ECOSYSTEM_ORDER:
    eco_sites = site_ecosystem[site_ecosystem['Salinity_Category'] == eco]['site_name'].tolist()
    eco_sites_in_shared = [s for s in eco_sites if s in shared_sites]
    # Sort by T:ET median
    eco_sites_sorted = sorted(eco_sites_in_shared, key=lambda x: tet_medians.get(x, 0))
    ordered_sites.extend(eco_sites_sorted)
    ordered_ecosystems.extend([eco] * len(eco_sites_sorted))
    print(f"  {eco}: {len(eco_sites_sorted)} sites")

print(f"\nTotal ordered sites: {len(ordered_sites)} (should match {len(shared_sites)})")

# =============================================================================
# STEP 6: PREPARE BOXPLOT DATA (use ALL monthly values from filtered data)
# =============================================================================

print("\n" + "="*60)
print("PREPARING BOXPLOT DATA")
print("="*60)

# Collect data for each site
trans_data = []
wue_data = []
wue_t_data = []

for site in ordered_sites:
    site_data = df_nn[df_nn['site_name'] == site]
    
    # Trans_ratio values
    trans_vals = site_data['Trans_ratio'].dropna().values
    trans_data.append(trans_vals if len(trans_vals) > 0 else [np.nan])
    
    # WUE_ET values (cap at 15 for visualization)
    wue_vals = site_data['WUE'].dropna().values
    wue_vals = np.where(wue_vals > 15, 15, wue_vals)
    wue_data.append(wue_vals if len(wue_vals) > 0 else [np.nan])
    
    # WUE_T values (cap at 15 for visualization)
    wue_t_vals = site_data['WUE_tra'].dropna().values
    wue_t_vals = np.where(wue_t_vals > 15, 15, wue_t_vals)
    wue_t_data.append(wue_t_vals if len(wue_t_vals) > 0 else [np.nan])

panel_data = {
    'Trans_ratio': trans_data,
    'WUE': wue_data,
    'WUE_tra': wue_t_data
}

print(f"\nTrans_ratio: {len(trans_data)} boxes")
print(f"WUE_ET: {len(wue_data)} boxes")
print(f"WUE_T: {len(wue_t_data)} boxes")

# =============================================================================
# STEP 7: CREATE THREE-PANEL FIGURE
# =============================================================================

print("\n" + "="*80)
print("CREATING THREE-PANEL BOXPLOT FIGURE")
print("="*80)

fig, axes = plt.subplots(3, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), 
                          sharex=False, gridspec_kw={'hspace': 0.15})

# Calculate ecosystem boundaries
ecosystem_boundaries = []
current_eco = None
start_idx = 0

for i, eco in enumerate(ordered_ecosystems):
    if current_eco is None:
        current_eco = eco
        start_idx = i
    elif eco != current_eco:
        ecosystem_boundaries.append((current_eco, start_idx, i - 1))
        current_eco = eco
        start_idx = i
if current_eco:
    ecosystem_boundaries.append((current_eco, start_idx, len(ordered_sites) - 1))

positions = np.arange(len(ordered_sites))

for idx, panel in enumerate(PANELS):
    ax = axes[idx]
    col = panel['col']
    boxplot_data = panel_data[col]
    box_colors = [ECOSYSTEM_COLORS[eco] for eco in ordered_ecosystems]
    
    # Create boxplot
    bp = ax.boxplot(boxplot_data, positions=positions, widths=BOX_WIDTH,
                    patch_artist=True, showfliers=False, vert=True)
    
    # Color boxes
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(box_colors[i])
        patch.set_alpha(BOX_ALPHA)
        patch.set_edgecolor('black')
        patch.set_linewidth(BOX_LINEWIDTH)
    
    # Style medians
    for median in bp['medians']:
        median.set_color('red')
        median.set_linewidth(MEDIAN_LINEWIDTH)
    
    # Style whiskers and caps
    for whisker in bp['whiskers']:
        whisker.set_color('black')
        whisker.set_linewidth(WHISKER_WIDTH)
    for cap in bp['caps']:
        cap.set_color('black')
        cap.set_linewidth(WHISKER_WIDTH)
    
    # X-AXIS: Site names ONLY on bottom panel (Panel C, index 2)
    if idx == 2:
        ax.set_xticks(positions)
        ax.set_xticklabels(ordered_sites, rotation=90, fontsize=X_TICK_FONT, ha='center')
    else:
        ax.set_xticks([])
        ax.set_xlabel('')
    
    # ECOSYSTEM LABELS ONLY IN TOP PANEL
    if idx == 0:
        y_top = ax.get_ylim()[1]
        y_text_pos = y_top * 1.10
        for eco, start_idx, end_idx in ecosystem_boundaries:
            center_x = (start_idx + end_idx) / 2
            ax.text(center_x, y_text_pos, eco, ha='center', va='bottom',
                    fontsize=ECOSYSTEM_LABEL_FONT, fontweight='bold', 
                    color=ECOSYSTEM_COLORS[eco])
    
    # VERTICAL SEPARATOR LINES
    for eco, start_idx, end_idx in ecosystem_boundaries:
        if start_idx > 0:
            ax.axvline(x=start_idx - 0.5, color='black', linestyle='--', alpha=0.5, linewidth=2.5)
    
    # Y-AXIS - SIMPLE FORMATTING FOR ALL PANELS
    y_min, y_max = Y_LIMITS.get(col, (0, None))
    if y_max:
        ax.set_ylim(bottom=y_min, top=y_max)
        y_interval = Y_TICK_INTERVALS.get(col, 0.2)
        if col == 'Trans_ratio':
            y_ticks = np.arange(y_min, y_max + y_interval, y_interval)
            ax.set_yticks(y_ticks)
            ax.set_yticklabels([f'{tick:.1f}' for tick in y_ticks], fontsize=Y_TICK_FONT)
        else:
            y_ticks = np.arange(y_min, y_max + y_interval, y_interval)
            ax.set_yticks(y_ticks)
            ax.set_yticklabels([f'{tick:.0f}' for tick in y_ticks], fontsize=Y_TICK_FONT)
    
    # Y-AXIS LABEL
    ax.set_ylabel(Y_LABELS[col], fontsize=AXIS_LABEL_FONT, fontweight='bold', labelpad=30)
    
    # SPINES
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    
    # GRID
    ax.grid(True, linestyle=':', alpha=0.3, axis='y')
    
    # PANEL LABEL
    panel_label = chr(97 + idx)
    ax.text(-0.10, 1.12, f'({panel_label})', transform=ax.transAxes,
            fontsize=PANEL_LABEL_FONT, fontweight='bold', va='top', ha='left')

# LEGEND
legend_elements = []
for eco in ECOSYSTEM_ORDER:
    legend_elements.append(
        plt.Rectangle((0, 0), 1, 1, facecolor=ECOSYSTEM_COLORS[eco], 
                     edgecolor='black', alpha=BOX_ALPHA, label=eco)
    )

fig.legend(handles=legend_elements, loc='center right', 
           bbox_to_anchor=(1.02, 0.5), ncol=1, fontsize=LEGEND_FONT + 10,
           frameon=True, fancybox=True, shadow=True, title='Ecosystem Type',
           title_fontsize=LEGEND_FONT + 10, handlelength=3.0, handleheight=2.0)

# ADJUST LAYOUT
plt.subplots_adjust(left=0.12, right=0.86, top=0.93, bottom=0.10, hspace=0.18)

# Save
plt.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_PDF, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {OUTPUT_PNG}")
print(f"✅ Saved: {OUTPUT_PDF}")

plt.show()

# =============================================================================
# PRINT FINAL SUMMARY
# =============================================================================

print("\n" + "="*80)
print("SUPPLEMENTARY FIGURE COMPLETE")
print("="*80)
print(f"\nSTRICT TRIPLE INTERSECTION SUMMARY (SAME AS FIGURE 1):")
print(f"  Total sites: {len(shared_sites)}")
print(f"  T:ET range: {min(tet_values_list):.4f} - {max(tet_values_list):.4f}")
print(f"  Sites with T:ET < 0.1: {len(low_tet_sites)}")
print(f"  Sites with T:ET > 0.7: {len(high_tet_sites)}")
print(f"\nEcosystem distribution:")
for eco in ECOSYSTEM_ORDER:
    count = ordered_ecosystems.count(eco)
    print(f"  {eco}: {count} sites")

print(f"\nFILTERING SUMMARY:")
print(f"  Strict triple intersection sites: {len(shared_sites)}")
print(f"  NN rows after strict metric filter: {len(df_nn):,}")
print(f"  Strict metric filter applied: WUE, WUE_eva, WUE_tra ALL have data")