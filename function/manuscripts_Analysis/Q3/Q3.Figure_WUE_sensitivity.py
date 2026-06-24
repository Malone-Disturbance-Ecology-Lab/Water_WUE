"""
Q3_Figure_WUE_sensitivity.py - CHUNK 2: Create Figures (FULLY FIXED)
EXACT REPLICATION OF MALONE'S FIGURE CREATION

FIXES:
1. Hard failure for Figure 4 composite creation
2. Hard failure if any expected figure is missing
3. Clean figure directory at start
4. 8 expected figures (Q3_plot_02, 03, 04, 06, 07, 09, 10, 11)
5. PATHS MATCH CHUNK 1 - outputs and figures separated
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("Q3: CHUNK 2 - Create Figures (FULLY FIXED)")
print("="*60)

# ============================================================================
# PATHS - MATCH CHUNK 1 (FIXED)
# ============================================================================

base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"

# Directories matching Chunk 1
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures")
temp_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_temp")

os.makedirs(figure_dir, exist_ok=True)

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")

# ============================================================================
# CONSTANTS - EXACTLY AS IN MALONE
# ============================================================================

SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
ECOSYSTEM_COLORS = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500"
}
COAST_COLORS = {
    "Atlantic Coast": "#008B8B",
    "Pacific Coast": "#C44E52",
    "Gulf Coast": "#8C564B",
    "AK Coast": "#4D4D4D"
}
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]
GRAY_COLOR = "#595959"

# Marker map for threshold_pct labels (Malone-style "5%", "10%", "20%")
THRESHOLD_MARKER_MAP = {"5%": "o", "10%": "s", "20%": "^"}

# Expected figures (8 total - Malone)
EXPECTED_FIGURES = [
    "Q3_plot_02_site_level_sensitivity_slopes.png",
    "Q3_plot_03_smooth_WUE_T_response_by_SPEI.png",
    "Q3_plot_04_WUE_T_SPEI_sensitivity_multipanel.png",
    "Q3_plot_06_WUE_T_selected_SPEI_timescales.png",
    "Q3_plot_07_most_sensitive_sites.png",
    "Q3_plot_09_all_system_sensitivity_heatmap.png",
    "Q3_plot_10_coast_specific_GAM_SPEI_thresholds_upland.png",
    "Q3_plot_11_coastline_grouped_site_sensitivity.png",
]

# ============================================================================
# CLEAN OLD FIGURES AT START
# ============================================================================

print("\n" + "="*60)
print("STEP 0: Cleaning old figures")
print("="*60)

for f in EXPECTED_FIGURES:
    path = os.path.join(figure_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old figure: {f}")

# Also clean temp panel files if they exist
for f in ["temp_panel_A.png", "temp_panel_B.png"]:
    path = os.path.join(figure_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed temp file: {f}")

# ============================================================================
# HELPER FUNCTIONS - MATCH MALONE'S SAVE_PLOT AND THEME_WUE
# ============================================================================

def save_plot(fig, filename, width=10, height=7, dpi=300):
    """Identical to Malone's save_plot()"""
    fig.set_size_inches(width, height)
    fig.savefig(os.path.join(figure_dir, filename), dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def theme_wue(base_size=13):
    """Identical to Malone's theme_wue()"""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': base_size,
        'axes.labelsize': base_size,
        'axes.titlesize': base_size + 2,
        'xtick.labelsize': base_size - 1,
        'ytick.labelsize': base_size - 1,
        'legend.fontsize': base_size - 1,
        'figure.titlesize': base_size + 4,
        'axes.titleweight': 'bold',
        'legend.title_fontsize': base_size - 1,
        'grid.alpha': 0.3,
    })

def get_threshold_marker(threshold_pct):
    """Get marker symbol for threshold percentage label"""
    threshold_label = str(threshold_pct)
    if threshold_label in THRESHOLD_MARKER_MAP:
        return THRESHOLD_MARKER_MAP[threshold_label]
    try:
        pct_val = int(float(threshold_label))
        label_map = {5: "5%", 10: "10%", 20: "20%"}
        return THRESHOLD_MARKER_MAP.get(label_map.get(pct_val, "10%"), "s")
    except:
        return "s"

# ============================================================================
# LOAD DATA
# ============================================================================

print("\n" + "="*60)
print("STEP 1: Loading data files")
print("="*60)

data_vars = {}
files_to_load = {
    "site_slopes": "Q3_WUE_T_SPEI_site_level_slopes.csv",
    "site_sensitivity_rank": "Q3_WUE_T_SPEI_site_sensitivity_rank.csv",
    "prediction_impact": "Q3_WUE_T_SPEI_predicted_impact_classes.csv",
    "gam_threshold_markers": "Q3_WUE_T_SPEI_GAM_impact_threshold_markers_5_10_20pct.csv",
    "coast_prediction_impact": "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv",
    "coast_threshold_markers": "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv",
    "smooth_predictions": "Q3_WUE_T_SPEI_smooth_predictions.csv",
}

for var_name, filename in files_to_load.items():
    filepath = os.path.join(output_dir, filename)
    if os.path.exists(filepath):
        data_vars[var_name] = pd.read_csv(filepath)
        print(f"  Loaded: {filename}")
    else:
        data_vars[var_name] = pd.DataFrame()
        print(f"  WARNING: {filename} not found")

# Extract variables for convenience
site_slopes = data_vars.get("site_slopes", pd.DataFrame())
site_sensitivity_rank = data_vars.get("site_sensitivity_rank", pd.DataFrame())
prediction_impact = data_vars.get("prediction_impact", pd.DataFrame())
gam_threshold_markers = data_vars.get("gam_threshold_markers", pd.DataFrame())
coast_prediction_impact = data_vars.get("coast_prediction_impact", pd.DataFrame())
coast_threshold_markers = data_vars.get("coast_threshold_markers", pd.DataFrame())
smooth_predictions = data_vars.get("smooth_predictions", pd.DataFrame())

# Convert threshold_pct to string for consistent handling
for df in [gam_threshold_markers, coast_threshold_markers]:
    if not df.empty and 'threshold_pct' in df.columns:
        if df['threshold_pct'].dtype in ['int64', 'float64']:
            df['threshold_pct'] = df['threshold_pct'].astype(int).astype(str) + '%'

print("\n  Data loaded successfully!")

# ============================================================================
# CREATE ALL INDIVIDUAL FIGURES
# ============================================================================

print("\n" + "="*60)
print("STEP 2: Creating individual figures")
print("="*60)
theme_wue()

# ============================================================================
# FIGURE 2: Site-level sensitivity slopes - DODGED BOXPLOTS
# ============================================================================

if not site_slopes.empty:
    print("\n  Creating Figure 2: Site-level sensitivity slopes...")
    
    fig, ax = plt.subplots(figsize=(11, 7))
    
    x_positions = np.arange(len(SPEI_COLS))
    width = 0.25
    
    for j, ecosystem in enumerate(ECOSYSTEM_CLASSES):
        subset = site_slopes[site_slopes['water_class'] == ecosystem]
        data_by_timescale = []
        positions = x_positions + (j - 1) * width
        
        for i, timescale in enumerate(SPEI_COLS):
            ts_data = subset[subset['SPEI_timescale'] == timescale]['SPEI_slope'].dropna().values
            if len(ts_data) > 0:
                data_by_timescale.append(ts_data)
            else:
                data_by_timescale.append([])
        
        if any(len(d) > 0 for d in data_by_timescale):
            bp = ax.boxplot(data_by_timescale, positions=positions, widths=width, patch_artist=True,
                          boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.5),
                          whiskerprops=dict(color=GRAY_COLOR),
                          capprops=dict(color=GRAY_COLOR),
                          medianprops=dict(color='black', linewidth=1.5))

    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(SPEI_COLS, rotation=30, ha='right')
    ax.set_xlabel('SPEI timescale')
    ax.set_ylabel('Site-level slope of WUE_T vs SPEI')
    ax.set_title('Site-Level WUE_T Sensitivity Slopes', fontweight='bold')
    
    legend_elements = [Patch(facecolor=ECOSYSTEM_COLORS[eco], alpha=0.5, label=eco) for eco in ECOSYSTEM_CLASSES]
    ax.legend(handles=legend_elements, loc='best', title='Ecosystem class')
    ax.grid(True, alpha=0.3)
    save_plot(fig, "Q3_plot_02_site_level_sensitivity_slopes.png", width=11, height=7)
    print("    Saved: Q3_plot_02_site_level_sensitivity_slopes.png")
else:
    print("  Skipped: Figure 2 - site_slopes data not available")

# ============================================================================
# FIGURE 3: Smooth WUE_T response by SPEI - EMPTY PANEL REMOVED
# ============================================================================

if not smooth_predictions.empty:
    print("\n  Creating Figure 3: Smooth WUE_T response by SPEI...")
    prediction_grid = smooth_predictions
    
    fig, axes = plt.subplots(4, 2, figsize=(12, 12))
    axes = axes.flatten()
    
    for idx, timescale in enumerate(SPEI_COLS):
        ax = axes[idx]
        subset = prediction_grid[prediction_grid['SPEI_timescale'] == timescale]
        
        for ecosystem in ECOSYSTEM_CLASSES:
            eco_subset = subset[subset['water_class'] == ecosystem]
            if len(eco_subset) > 0:
                ax.fill_between(
                    eco_subset['SPEI_value'],
                    eco_subset['predicted_lower'],
                    eco_subset['predicted_upper'],
                    color=ECOSYSTEM_COLORS[ecosystem],
                    alpha=0.14
                )
                ax.plot(
                    eco_subset['SPEI_value'],
                    eco_subset['predicted_WUE_T'],
                    color=ECOSYSTEM_COLORS[ecosystem],
                    linewidth=0.85,
                    label=ecosystem
                )
        
        ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.5)
        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.5)
        ax.set_title(timescale, fontweight='bold')
        ax.set_xlabel('SPEI')
        ax.set_ylabel('Predicted WUE_T')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Remove empty 8th panel (Malone-style)
    fig.delaxes(axes[-1])
    
    plt.suptitle('GAM WUE_T Response Across SPEI Gradients', fontsize=16, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "Q3_plot_03_smooth_WUE_T_response_by_SPEI.png", width=12, height=12)
    print("    Saved: Q3_plot_03_smooth_WUE_T_response_by_SPEI.png")
else:
    print("  Skipped: Figure 3 - smooth_predictions data not available")

# ============================================================================
# FIGURE 6: Selected SPEI timescales with thresholds
# ============================================================================

if not prediction_impact.empty:
    print("\n  Creating Figure 6: Selected SPEI timescales...")
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 6))
    
    for idx, timescale in enumerate(SELECTED_TIMESCALES):
        ax = axes[idx]
        subset = prediction_impact[prediction_impact['SPEI_timescale'] == timescale]
        
        for ecosystem in ECOSYSTEM_CLASSES:
            eco_subset = subset[subset['water_class'] == ecosystem]
            if len(eco_subset) > 0:
                ax.fill_between(
                    eco_subset['SPEI_value'],
                    eco_subset['predicted_lower_pct'],
                    eco_subset['predicted_upper_pct'],
                    color=ECOSYSTEM_COLORS[ecosystem],
                    alpha=0.14
                )
                ax.plot(
                    eco_subset['SPEI_value'],
                    eco_subset['predicted_pct_change'],
                    color=ECOSYSTEM_COLORS[ecosystem],
                    linewidth=0.9,
                    label=ecosystem
                )
        
        ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.5)
        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.5)
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
        ax.axhline(y=-10, color='gray', linestyle=':', linewidth=0.5)
        ax.axhline(y=10, color='gray', linestyle=':', linewidth=0.5)
        
        if not gam_threshold_markers.empty:
            markers_subset = gam_threshold_markers[gam_threshold_markers['SPEI_timescale'] == timescale]
            for _, row in markers_subset.iterrows():
                if row['water_class'] in ECOSYSTEM_COLORS:
                    color = ECOSYSTEM_COLORS[row['water_class']]
                    marker = get_threshold_marker(row['threshold_pct'])
                    
                    ax.hlines(
                        y=row['pct_change_threshold'],
                        xmin=row['SPEI_threshold_lower'],
                        xmax=row['SPEI_threshold_upper'],
                        color=color,
                        linewidth=1.5,
                        alpha=0.7
                    )
                    ax.scatter(
                        x=row['SPEI_threshold'],
                        y=row['pct_change_threshold'],
                        color=color,
                        s=80,
                        marker=marker,
                        edgecolor='black',
                        linewidth=1,
                        zorder=5
                    )
        
        ax.set_title(timescale, fontweight='bold')
        ax.set_xlabel('SPEI')
        ax.set_ylabel('Predicted WUE_T change from near-normal (%)')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('GAM-Predicted WUE_T Change from Near-Normal Conditions', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "Q3_plot_06_WUE_T_selected_SPEI_timescales.png", width=14, height=6)
    print("    Saved: Q3_plot_06_WUE_T_selected_SPEI_timescales.png")
else:
    print("  Skipped: Figure 6 - prediction_impact data not available")

# ============================================================================
# FIGURE 7: Most sensitive sites
# ============================================================================

if not site_sensitivity_rank.empty:
    print("\n  Creating Figure 7: Most sensitive sites...")
    top_sensitive_sites = site_sensitivity_rank[
        site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES)
    ].groupby('SPEI_timescale', observed=True).apply(
        lambda x: x.nlargest(12, 'mean_abs_pct_change')
    ).reset_index(drop=True)
    
    if len(top_sensitive_sites) > 0:
        top_sensitive_sites['site_name'] = pd.Categorical(
            top_sensitive_sites['site_name'],
            categories=top_sensitive_sites.sort_values('mean_abs_pct_change')['site_name'].unique()
        )
        
        fig, axes = plt.subplots(1, 3, figsize=(13, 8))
        
        for idx, timescale in enumerate(SELECTED_TIMESCALES):
            ax = axes[idx]
            subset = top_sensitive_sites[top_sensitive_sites['SPEI_timescale'] == timescale]
            
            if len(subset) > 0:
                subset = subset.sort_values('mean_abs_pct_change')
                ax.barh(
                    subset['site_name'],
                    subset['mean_abs_pct_change'],
                    color=[COAST_COLORS.get(c, '#808080') for c in subset['coast_region']],
                    height=0.72
                )
                ax.set_title(timescale, fontweight='bold')
                ax.set_xlabel('Mean absolute WUE_T change (%)')
                ax.grid(True, alpha=0.3, axis='x')
        
        plt.suptitle('Most Sensitive Sites - Ranked by Mean Absolute WUE_T Change', fontsize=14, fontweight='bold')
        plt.tight_layout()
        save_plot(fig, "Q3_plot_07_most_sensitive_sites.png", width=13, height=8)
        print("    Saved: Q3_plot_07_most_sensitive_sites.png")
    else:
        print("    Skipped: Figure 7 - no top sensitive sites found")
else:
    print("  Skipped: Figure 7 - site_sensitivity_rank data not available")

# ============================================================================
# FIGURE 9: All system sensitivity
# ============================================================================

if not site_sensitivity_rank.empty:
    print("\n  Creating Figure 9: All system sensitivity...")
    all_system_sensitivity = site_sensitivity_rank[
        site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES) &
        (site_sensitivity_rank['n_months'] >= 6)
    ]
    
    if len(all_system_sensitivity) > 0:
        site_order = all_system_sensitivity.groupby('site_name')['mean_abs_pct_change'].max().sort_values().index
        all_system_sensitivity['site_name'] = pd.Categorical(
            all_system_sensitivity['site_name'],
            categories=site_order
        )
        
        fig, ax = plt.subplots(figsize=(13, 13))
        
        shapes = {'SPEI_1': 'o', 'SPEI_3': 's', 'SPEI_48': '^'}
        
        for timescale in SELECTED_TIMESCALES:
            subset = all_system_sensitivity[all_system_sensitivity['SPEI_timescale'] == timescale]
            ax.scatter(
                subset['mean_abs_pct_change'],
                subset['site_name'],
                color=[ECOSYSTEM_COLORS.get(c, '#808080') for c in subset['water_class']],
                s=80,
                marker=shapes.get(timescale, 'o'),
                alpha=0.88,
                label=timescale,
                edgecolor='black',
                linewidth=0.5
            )
        
        ax.axvline(x=5, color='gray', linestyle=':', linewidth=1, alpha=0.7)
        ax.axvline(x=10, color='gray', linestyle=':', linewidth=1, alpha=0.7)
        ax.axvline(x=20, color='gray', linestyle=':', linewidth=1, alpha=0.7)
        
        ax.set_xlabel('Mean absolute WUE_T change from near-normal (%)')
        ax.set_ylabel('Site')
        ax.set_title('Sensitivity of All Systems', fontweight='bold')
        
        eco_handles = [Patch(facecolor=ECOSYSTEM_COLORS.get(c, '#808080'), label=c) for c in ECOSYSTEM_CLASSES]
        shape_handles = []
        for ts, marker in shapes.items():
            shape_handles.append(plt.Line2D([0], [0], marker=marker, color='black', linestyle='None', 
                                            markersize=10, label=ts))
        
        legend1 = ax.legend(handles=eco_handles, loc='upper right', title='Ecosystem class')
        ax.add_artist(legend1)
        ax.legend(handles=shape_handles, loc='lower right', title='SPEI timescale')
        
        ax.grid(True, alpha=0.3, axis='x')
        plt.tight_layout()
        save_plot(fig, "Q3_plot_09_all_system_sensitivity_heatmap.png", width=13, height=13)
        print("    Saved: Q3_plot_09_all_system_sensitivity_heatmap.png")
    else:
        print("    Skipped: Figure 9 - no data for selected timescales")
else:
    print("  Skipped: Figure 9 - site_sensitivity_rank data not available")

# ============================================================================
# FIGURE 10: Coast-specific GAM thresholds
# ============================================================================

print("\n  Creating Figure 10: Coast-specific GAM thresholds...")
if not coast_prediction_impact.empty and not coast_threshold_markers.empty:
    
    fig, axes = plt.subplots(4, 3, figsize=(12, 10))
    
    for idx_c, coast in enumerate(COAST_REGION_LEVELS):
        for idx_t, timescale in enumerate(SELECTED_TIMESCALES):
            ax = axes[idx_c, idx_t]
            
            subset = coast_prediction_impact[
                (coast_prediction_impact['coast_region'] == coast) &
                (coast_prediction_impact['SPEI_timescale'] == timescale)
            ]
            
            if len(subset) > 0:
                ax.fill_between(
                    subset['SPEI_value'],
                    subset['predicted_lower_pct'],
                    subset['predicted_upper_pct'],
                    color=COAST_COLORS.get(coast, '#808080'),
                    alpha=0.1
                )
                ax.plot(
                    subset['SPEI_value'],
                    subset['predicted_pct_change'],
                    color=COAST_COLORS.get(coast, '#808080'),
                    linewidth=0.9,
                    label=coast
                )
                
                markers_subset = coast_threshold_markers[
                    (coast_threshold_markers['coast_region'] == coast) &
                    (coast_threshold_markers['SPEI_timescale'] == timescale)
                ]
                for _, row in markers_subset.iterrows():
                    marker = get_threshold_marker(row['threshold_pct'])
                    
                    ax.hlines(
                        y=row['pct_change_threshold'],
                        xmin=row['SPEI_threshold_lower'],
                        xmax=row['SPEI_threshold_upper'],
                        color=COAST_COLORS.get(coast, '#808080'),
                        linewidth=1.5,
                        alpha=0.7
                    )
                    ax.scatter(
                        x=row['SPEI_threshold'],
                        y=row['pct_change_threshold'],
                        color=COAST_COLORS.get(coast, '#808080'),
                        s=60,
                        marker=marker,
                        edgecolor='black',
                        linewidth=1,
                        zorder=5
                    )
            
            ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.5)
            ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.5)
            ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            ax.axhline(y=-20, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=-10, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=-5, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=5, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=10, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=20, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            
            if idx_c == 0:
                ax.set_title(timescale, fontweight='bold')
            if idx_t == 0:
                ax.set_ylabel(coast, fontweight='bold')
            if idx_c == 3:
                ax.set_xlabel('SPEI')
    
    plt.suptitle('Coast-Specific GAM Thresholds for Upland WUE_T Response', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "Q3_plot_10_coast_specific_GAM_SPEI_thresholds_upland.png", width=12, height=10)
    print("    Saved: Q3_plot_10_coast_specific_GAM_SPEI_thresholds_upland.png")
else:
    print("  Skipped: Figure 10 - coast data not available")

# ============================================================================
# FIGURE 11: Coastline grouped sensitivity - FLIPPED TO MATCH MALONE
# ============================================================================

print("\n  Creating Figure 11: Coastline grouped sensitivity (Malone-style flipped)...")
if not site_sensitivity_rank.empty:
    all_system_sensitivity = site_sensitivity_rank[
        site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES) &
        (site_sensitivity_rank['n_months'] >= 6)
    ]
    
    if len(all_system_sensitivity) > 0:
        fig, axes = plt.subplots(1, 3, figsize=(14, 6))
        
        for idx, timescale in enumerate(SELECTED_TIMESCALES):
            ax = axes[idx]
            subset = all_system_sensitivity[all_system_sensitivity['SPEI_timescale'] == timescale]
            
            if len(subset) > 0:
                data_for_box = []
                coast_order = COAST_REGION_LEVELS
                for coast in coast_order:
                    coast_subset = subset[subset['coast_region'] == coast]
                    if len(coast_subset) > 0:
                        data_for_box.append(coast_subset['mean_abs_pct_change'].values)
                    else:
                        data_for_box.append([])
                
                positions = range(len(coast_order))
                bp = ax.boxplot(data_for_box, positions=positions, widths=0.6, patch_artist=True,
                              vert=False,
                              boxprops=dict(facecolor='lightgray', alpha=0.5),
                              whiskerprops=dict(color=GRAY_COLOR),
                              capprops=dict(color=GRAY_COLOR),
                              medianprops=dict(color='black', linewidth=1.5))
                
                for ecosystem in ECOSYSTEM_CLASSES:
                    eco_subset = subset[subset['water_class'] == ecosystem]
                    if len(eco_subset) > 0:
                        for i, coast in enumerate(coast_order):
                            coast_eco_subset = eco_subset[eco_subset['coast_region'] == coast]
                            if len(coast_eco_subset) > 0:
                                y_jitter = np.random.normal(i, 0.08, len(coast_eco_subset))
                                ax.scatter(coast_eco_subset['mean_abs_pct_change'], y_jitter,
                                         color=ECOSYSTEM_COLORS[ecosystem], s=60, alpha=0.7,
                                         edgecolor='black', linewidth=0.5, label=ecosystem if i == 0 else "")
                
                ax.axvline(x=5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
                ax.axvline(x=10, color='gray', linestyle=':', linewidth=1, alpha=0.5)
                ax.axvline(x=20, color='gray', linestyle=':', linewidth=1, alpha=0.5)
                
                ax.set_yticks(range(len(coast_order)))
                ax.set_yticklabels(coast_order)
                ax.set_title(timescale, fontweight='bold')
                ax.set_xlabel('Mean absolute WUE_T change from near-normal (%)')
                if idx == 0:
                    ax.set_ylabel('Coast region')
                ax.grid(True, alpha=0.3, axis='x')
        
        handles = [Patch(facecolor=ECOSYSTEM_COLORS[eco], alpha=0.7, label=eco) for eco in ECOSYSTEM_CLASSES]
        fig.legend(handles=handles, loc='lower center', ncol=3, title='Ecosystem class', fontsize=11)
        
        plt.suptitle('Site Sensitivity Grouped by Coastline', fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0.05, 1, 0.95])
        save_plot(fig, "Q3_plot_11_coastline_grouped_site_sensitivity.png", width=14, height=6)
        print("    Saved: Q3_plot_11_coastline_grouped_site_sensitivity.png")
    else:
        print("    Skipped: Figure 11 - no data for selected timescales")
else:
    print("  Skipped: Figure 11 - site_sensitivity_rank data not available")

# ============================================================================
# FIGURE 4: COMPOSITE MULTIPANEL - WITH HARD FAILURE
# ============================================================================

print("\n  Creating Figure 04: Multipanel composite (Malone-style)...")

# Check if we have the needed data for the composite
if not coast_prediction_impact.empty and not coast_threshold_markers.empty and not site_sensitivity_rank.empty:
    
    # Panel A: Coast-specific GAM thresholds (top panel)
    fig_panel_a, axes_a = plt.subplots(4, 3, figsize=(12, 10))
    
    for idx_c, coast in enumerate(COAST_REGION_LEVELS):
        for idx_t, timescale in enumerate(SELECTED_TIMESCALES):
            ax = axes_a[idx_c, idx_t]
            
            subset = coast_prediction_impact[
                (coast_prediction_impact['coast_region'] == coast) &
                (coast_prediction_impact['SPEI_timescale'] == timescale)
            ]
            
            if len(subset) > 0:
                ax.fill_between(
                    subset['SPEI_value'],
                    subset['predicted_lower_pct'],
                    subset['predicted_upper_pct'],
                    color=COAST_COLORS.get(coast, '#808080'),
                    alpha=0.1
                )
                ax.plot(
                    subset['SPEI_value'],
                    subset['predicted_pct_change'],
                    color=COAST_COLORS.get(coast, '#808080'),
                    linewidth=0.9,
                    label=coast
                )
                
                markers_subset = coast_threshold_markers[
                    (coast_threshold_markers['coast_region'] == coast) &
                    (coast_threshold_markers['SPEI_timescale'] == timescale)
                ]
                for _, row in markers_subset.iterrows():
                    marker = get_threshold_marker(row['threshold_pct'])
                    
                    ax.hlines(
                        y=row['pct_change_threshold'],
                        xmin=row['SPEI_threshold_lower'],
                        xmax=row['SPEI_threshold_upper'],
                        color=COAST_COLORS.get(coast, '#808080'),
                        linewidth=1.5,
                        alpha=0.7
                    )
                    ax.scatter(
                        x=row['SPEI_threshold'],
                        y=row['pct_change_threshold'],
                        color=COAST_COLORS.get(coast, '#808080'),
                        s=60,
                        marker=marker,
                        edgecolor='black',
                        linewidth=1,
                        zorder=5
                    )
            
            ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.5)
            ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.5)
            ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            ax.axhline(y=-20, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=-10, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=-5, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=5, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=10, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            ax.axhline(y=20, color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
            
            if idx_c == 0:
                ax.set_title(timescale, fontweight='bold')
            if idx_t == 0:
                ax.set_ylabel(coast, fontweight='bold')
            if idx_c == 3:
                ax.set_xlabel('SPEI')
    
    plt.suptitle('Coast-Specific GAM Thresholds for Upland WUE_T Response', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Panel B: Coastline grouped sensitivity (bottom panel) - FLIPPED
    fig_panel_b, axes_b = plt.subplots(1, 3, figsize=(14, 6))
    
    all_system_sensitivity = site_sensitivity_rank[
        site_sensitivity_rank['SPEI_timescale'].isin(SELECTED_TIMESCALES) &
        (site_sensitivity_rank['n_months'] >= 6)
    ]
    
    if len(all_system_sensitivity) > 0:
        for idx, timescale in enumerate(SELECTED_TIMESCALES):
            ax = axes_b[idx]
            subset = all_system_sensitivity[all_system_sensitivity['SPEI_timescale'] == timescale]
            
            if len(subset) > 0:
                data_for_box = []
                coast_order = COAST_REGION_LEVELS
                for coast in coast_order:
                    coast_subset = subset[subset['coast_region'] == coast]
                    if len(coast_subset) > 0:
                        data_for_box.append(coast_subset['mean_abs_pct_change'].values)
                    else:
                        data_for_box.append([])
                
                positions = range(len(coast_order))
                bp = ax.boxplot(data_for_box, positions=positions, widths=0.6, patch_artist=True,
                              vert=False,
                              boxprops=dict(facecolor='lightgray', alpha=0.5),
                              whiskerprops=dict(color=GRAY_COLOR),
                              capprops=dict(color=GRAY_COLOR),
                              medianprops=dict(color='black', linewidth=1.5))
                
                for ecosystem in ECOSYSTEM_CLASSES:
                    eco_subset = subset[subset['water_class'] == ecosystem]
                    if len(eco_subset) > 0:
                        for i, coast in enumerate(coast_order):
                            coast_eco_subset = eco_subset[eco_subset['coast_region'] == coast]
                            if len(coast_eco_subset) > 0:
                                y_jitter = np.random.normal(i, 0.08, len(coast_eco_subset))
                                ax.scatter(coast_eco_subset['mean_abs_pct_change'], y_jitter,
                                         color=ECOSYSTEM_COLORS[ecosystem], s=60, alpha=0.7,
                                         edgecolor='black', linewidth=0.5)
                
                ax.axvline(x=5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
                ax.axvline(x=10, color='gray', linestyle=':', linewidth=1, alpha=0.5)
                ax.axvline(x=20, color='gray', linestyle=':', linewidth=1, alpha=0.5)
                
                ax.set_yticks(range(len(coast_order)))
                ax.set_yticklabels(coast_order)
                ax.set_title(timescale, fontweight='bold')
                ax.set_xlabel('Mean absolute WUE_T change from near-normal (%)')
                if idx == 0:
                    ax.set_ylabel('Coast region')
                ax.grid(True, alpha=0.3, axis='x')
    
    plt.suptitle('Site Sensitivity Grouped by Coastline', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save panel A and B as temp files
    temp_a_path = os.path.join(figure_dir, "temp_panel_A.png")
    temp_b_path = os.path.join(figure_dir, "temp_panel_B.png")
    fig_panel_a.savefig(temp_a_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig_panel_b.savefig(temp_b_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig_panel_a)
    plt.close(fig_panel_b)
    
    # Create composite using PIL
    try:
        from PIL import Image, ImageDraw, ImageFont
        
        if os.path.exists(temp_a_path) and os.path.exists(temp_b_path):
            img_a = Image.open(temp_a_path)
            img_b = Image.open(temp_b_path)
            
            composite_width = max(img_a.width, img_b.width)
            img_a_resized = img_a.resize((composite_width, img_a.height))
            img_b_resized = img_b.resize((composite_width, img_b.height))
            
            composite_height = img_a_resized.height + img_b_resized.height + 50
            composite = Image.new('RGB', (composite_width, composite_height), 'white')
            composite.paste(img_a_resized, (0, 0))
            composite.paste(img_b_resized, (0, img_a_resized.height + 20))
            
            draw = ImageDraw.Draw(composite)
            try:
                font = ImageFont.truetype("arial.ttf", 40)
            except:
                font = ImageFont.load_default()
            
            draw.text((20, 20), "A", fill='black', font=font)
            draw.text((20, img_a_resized.height + 30), "B", fill='black', font=font)
            
            composite_path = os.path.join(figure_dir, "Q3_plot_04_WUE_T_SPEI_sensitivity_multipanel.png")
            composite.save(composite_path, dpi=(300, 300))
            print("    Saved: Q3_plot_04_WUE_T_SPEI_sensitivity_multipanel.png")
            
            # Clean up temp files
            os.remove(temp_a_path)
            os.remove(temp_b_path)
        else:
            raise RuntimeError("Temporary panel files not found")
            
    except ImportError:
        print("    WARNING: PIL not installed - cannot create composite")
        print("    Saved individual panels instead:")
        print("      - temp_panel_A.png")
        print("      - temp_panel_B.png")
        raise RuntimeError("PIL not installed - cannot create Figure 4 composite")
else:
    print("  Skipped: Figure 04 - required data not available")

# ============================================================================
# STEP 3: VERIFY ALL FIGURES - HARD FAILURE
# ============================================================================

print("\n" + "="*60)
print("STEP 3: Verifying all figures (hard failure if missing)")
print("="*60)

missing_figures = []
for filename in EXPECTED_FIGURES:
    filepath = os.path.join(figure_dir, filename)
    if not os.path.exists(filepath):
        missing_figures.append(filename)

if missing_figures:
    raise FileNotFoundError(
        "Missing expected Q3 figure files:\n" + "\n".join(missing_figures)
    )
else:
    print(f"All {len(EXPECTED_FIGURES)} expected figures created successfully!")

# ============================================================================
# CLEAN UP TEMP FILES
# ============================================================================

print("\n" + "="*60)
print("STEP 4: Cleaning up temporary files")
print("="*60)

for f in ["temp_panel_A.png", "temp_panel_B.png"]:
    path = os.path.join(figure_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed temp file: {f}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("CHUNK 2 COMPLETE! All figures saved.")
print("="*60)
print(f"\nFigures saved to: {figure_dir}")
print("\n8 figures created (Malone-style):")
print("  - Q3_plot_02_site_level_sensitivity_slopes.png")
print("  - Q3_plot_03_smooth_WUE_T_response_by_SPEI.png")
print("  - Q3_plot_04_WUE_T_SPEI_sensitivity_multipanel.png")
print("  - Q3_plot_06_WUE_T_selected_SPEI_timescales.png")
print("  - Q3_plot_07_most_sensitive_sites.png")
print("  - Q3_plot_09_all_system_sensitivity_heatmap.png")
print("  - Q3_plot_10_coast_specific_GAM_SPEI_thresholds_upland.png")
print("  - Q3_plot_11_coastline_grouped_site_sensitivity.png")
print("\nFIXES APPLIED:")
print("  - Paths match Chunk 1 (outputs/figures separated)")
print("  - Hard failure for Figure 4 composite")
print("  - Hard failure if any expected figure missing")
print("  - Clean figure directory at start")
print("="*60)