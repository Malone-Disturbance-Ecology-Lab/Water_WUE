# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:58:37 2026
@author: ammar

UPDATED: Now uses SAME filtering criteria as Figure 1
- Reads monthly_data_after_outlier_removal.csv (already filtered by CHUNK 3)
- Strict month filter: only months where WUE, WUE_eva, WUE_tra ALL have data
- Sites must have ALL THREE metrics (WUE, WUE_eva, WUE_tra) - strict triple intersection
- Global 3-month NN filter applied (same as Figure 1)
- Uses wue_site_level_NN_medians_SPEI_1.csv for site selection (same as Figure 1)
- No downstream output file names changed
- FIXED: Removed duplicate function with old drop_duplicates() logic
- FIXED: All NN observation counts now use len(site_data) not unique month combos
- UPDATED: New color gradient (vivid purple → blue-purple → muted brown-gold → golden yellow)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, theilslopes
import os
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

# ------------------------------
# Configuration
# ------------------------------
# UPDATED: Use FILTERED monthly data from upstream CHUNK 3 (not raw data)
MONTHLY_FILTERED_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"

# NN medians file (same as Figure 1 - for site selection)
NN_MEDIANS_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"

# Output directory
OUTPUT_DIR = r"M:\Research\WUE_CUE\data_products\results\T_ET_ratio"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files (SAME NAMES - no downstream impact)
OUTPUT_SITE_CSV = os.path.join(OUTPUT_DIR, "TET_gradient_site_NN_medians.csv")
OUTPUT_GROUP_SUMMARY = os.path.join(OUTPUT_DIR, "TET_gradient_group_summary.csv")
OUTPUT_FIGURE_PNG = os.path.join(OUTPUT_DIR, "Figure_TET_gradient_WUE_divergence.png")
OUTPUT_FIGURE_PDF = os.path.join(OUTPUT_DIR, "Figure_TET_gradient_WUE_divergence.pdf")

# Near-normal baseline (SAME SPEI RANGE as NN definition)
SPEI_MIN = -0.5
SPEI_MAX = 0.5

# Minimum NN months required per site (SAME as Figure 1)
MIN_NN_MONTHS = 3

# Set font sizes
plt.rcParams.update({
    'font.size': 21,
    'axes.titlesize': 0,
    'axes.labelsize': 22.5,
    'xtick.labelsize': 21,
    'ytick.labelsize': 21,
    'legend.fontsize': 19.5,
})

# =============================================================================
# STEP 1: GET STRICT TRIPLE INTERSECTION SITES (SAME AS FIGURE 1)
# =============================================================================

def get_strict_triple_intersection_sites():
    """
    Get the SAME strict triple intersection sites used in Figure 1
    This ensures consistency across all analyses
    
    Returns:
        shared_sites: set of site names that have ALL THREE metrics
        obs_counts: dict with NN observation counts per site
    """
    print("\n" + "="*80)
    print("STEP 1: LOADING STRICT TRIPLE INTERSECTION SITES (SAME AS FIGURE 1)")
    print("="*80)
    
    if not os.path.exists(NN_MEDIANS_FILE):
        print(f"❌ ERROR: {NN_MEDIANS_FILE} not found!")
        print("   Please run Figure 1 workflow first to create this file.")
        return None, None
    
    nn_data = pd.read_csv(NN_MEDIANS_FILE)
    
    # Get sites with WUE available
    wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
    
    # Get sites with WUE_eva available
    eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
    
    # Get sites with WUE_tra available
    tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
    
    # STRICT TRIPLE INTERSECTION: sites with ALL THREE metrics
    shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
    
    print(f"\n  Sites with WUE (bulk): {len(wue_sites)}")
    print(f"  Sites with WUE_E (evaporation): {len(eva_sites)}")
    print(f"  Sites with WUE_T (transpiration): {len(tra_sites)}")
    print(f"  Strict triple intersection (ALL metrics): {len(shared_sites)} sites")
    
    # =========================================================================
    # STEP 2: APPLY GLOBAL 3-MONTH NN FILTER (SAME AS FIGURE 1)
    # =========================================================================
    print(f"\n" + "="*80)
    print(f"STEP 2: APPLYING GLOBAL {MIN_NN_MONTHS}-MONTH NN FILTER (SAME AS FIGURE 1)")
    print("="*80)
    
    # Load monthly filtered data (already has strict month filter from CHUNK 3)
    if not os.path.exists(MONTHLY_FILTERED_FILE):
        print(f"⚠️ Warning: {MONTHLY_FILTERED_FILE} not found!")
        print("   Cannot apply global 3-month NN filter.")
        return shared_sites, {}
    
    df_monthly = pd.read_csv(MONTHLY_FILTERED_FILE)
    
    # Filter to NN conditions at SPEI-1 (same as Figure 1)
    if "SPEI_1_Cat" in df_monthly.columns:
        df_nn = df_monthly[df_monthly["SPEI_1_Cat"] == "NN"].copy()
    else:
        print(f"⚠️ Warning: SPEI_1_Cat column not found")
        df_nn = df_monthly[(df_monthly["SPEI_1"] >= SPEI_MIN) & (df_monthly["SPEI_1"] <= SPEI_MAX)].copy()
    
    # Count retained NN observations per site - FIXED: count ALL retained observations
    obs_counts = {}
    for site in shared_sites:
        site_data = df_nn[df_nn["site_name"] == site]
        n_obs = len(site_data)  # FIXED: Count ALL retained NN observations (not unique month combos)
        obs_counts[site] = n_obs
    
    # Filter sites to those with >= MIN_NN_MONTHS observations
    sites_to_keep = {site for site in shared_sites if obs_counts.get(site, 0) >= MIN_NN_MONTHS}
    sites_to_remove = shared_sites - sites_to_keep
    
    print(f"\n  Original shared sites: {len(shared_sites)}")
    print(f"  Sites removed (<{MIN_NN_MONTHS} observations): {len(sites_to_remove)}")
    print(f"  Sites kept (≥{MIN_NN_MONTHS} observations): {len(sites_to_keep)}")
    
    if sites_to_remove:
        print(f"\n  Sites REMOVED (insufficient NN observations):")
        for site in sorted(sites_to_remove):
            obs = obs_counts.get(site, 0)
            print(f"      - {site}: {obs} observations")
    
    return sites_to_keep, obs_counts

# =============================================================================
# STEP 3: LOAD AND FILTER MONTHLY DATA (ALREADY FILTERED BY UPSTREAM CHUNK 3)
# =============================================================================

def load_filtered_monthly_data(kept_sites):
    """
    Load monthly data that has already been filtered by upstream CHUNK 3
    (strict month filter: only months where WUE, WUE_eva, WUE_tra ALL have data)
    """
    print("\n" + "="*80)
    print("STEP 3: LOADING FILTERED MONTHLY DATA (FROM UPSTREAM CHUNK 3)")
    print("="*80)
    print("  NOTE: Data already has strict month filter applied")
    print("  (only months where WUE, WUE_eva, WUE_tra ALL have data)")
    
    if not os.path.exists(MONTHLY_FILTERED_FILE):
        raise FileNotFoundError(f"Filtered monthly data not found: {MONTHLY_FILTERED_FILE}")
    
    df = pd.read_csv(MONTHLY_FILTERED_FILE)
    print(f"  Loaded {len(df):,} rows from filtered monthly data")
    
    # Filter to NN conditions (SPEI between -0.5 and 0.5)
    nn_mask = (df["SPEI_1"] >= SPEI_MIN) & (df["SPEI_1"] <= SPEI_MAX)
    df_nn = df[nn_mask].copy()
    print(f"  NN rows (SPEI between {SPEI_MIN} and {SPEI_MAX}): {len(df_nn)}")
    
    # Filter to strict triple intersection sites ONLY
    df_nn = df_nn[df_nn['site_name'].isin(kept_sites)].copy()
    print(f"  NN rows after filtering to kept sites: {len(df_nn)}")
    print(f"  Unique sites retained: {df_nn['site_name'].nunique()}")
    print(f"  Expected sites: {len(kept_sites)}")
    
    # Verify we have the correct number of sites
    if df_nn['site_name'].nunique() != len(kept_sites):
        print(f"⚠️ WARNING: Expected {len(kept_sites)} sites but found {df_nn['site_name'].nunique()}")
    
    return df_nn

# =============================================================================
# STEP 4: CALCULATE SITE-LEVEL MEDIANS
# =============================================================================

def calculate_site_medians(df_nn, kept_sites):
    """Calculate site-level medians for all metrics"""
    print("\n" + "="*80)
    print("STEP 4: CALCULATING SITE-LEVEL MEDIANS")
    print("="*80)
    
    site_medians = df_nn.groupby("site_name").agg(
        median_Trans_ratio=("Trans_ratio", "median"),
        median_WUE=("WUE", "median"),
        median_WUE_eva=("WUE_eva", "median"),
        median_WUE_tra=("WUE_tra", "median"),
        n_nn_observations=("Trans_ratio", "count"),  # Renamed from n_nn_months
        water_class=("water_class", "first"),
        IGBP=("IGBP", "first"),
        climate=("climate", "first"),
        lat=("lat", "first"),
        long=("long", "first"),
        Salinity_Category=("Salinity_Category", "first")
    ).reset_index()
    
    print(f"  Site-level medians calculated for {len(site_medians)} sites")
    
    # Check for any missing values
    required_vars = ["median_Trans_ratio", "median_WUE", "median_WUE_eva", "median_WUE_tra"]
    before_drop = len(site_medians)
    for var in required_vars:
        missing = site_medians[var].isna().sum()
        if missing > 0:
            print(f"  Warning: {var} has {missing} missing values")
        site_medians = site_medians.dropna(subset=[var])
    
    n_sites = len(site_medians)
    print(f"\n  Number of retained sites (all 4 metrics valid): {n_sites}")
    
    if n_sites != len(kept_sites):
        print(f"⚠️ WARNING: Expected {len(kept_sites)} sites but found {n_sites} after dropping NaN")
    
    return site_medians

# =============================================================================
# STEP 5: MAIN EXECUTION
# =============================================================================

print("="*80)
print("T:ET GRADIENT ANALYSIS WITH SAME FILTERING AS FIGURE 1")
print("="*80)
print(f"SPEI range for NN: {SPEI_MIN} to {SPEI_MAX}")
print(f"Minimum NN observations per site: {MIN_NN_MONTHS}")
print("Filtering rules (SAME as Figure 1):")
print("  1. Strict month filter: months where WUE, WUE_eva, WUE_tra ALL have data")
print("  2. Strict triple intersection: sites with ALL THREE metrics")
print(f"  3. Global {MIN_NN_MONTHS}-observation NN filter")
print("="*80)

# Get strict triple intersection sites (same as Figure 1)
kept_sites, obs_counts = get_strict_triple_intersection_sites()

if kept_sites is None or len(kept_sites) == 0:
    raise ValueError("Could not load strict triple intersection sites. Please run Figure 1 first.")

print(f"\n✅ Using {len(kept_sites)} sites from strict triple intersection (same as Figure 1)")

# Load filtered monthly data (already has strict month filter from CHUNK 3)
df_nn = load_filtered_monthly_data(kept_sites)

# Calculate site-level medians
site_medians = calculate_site_medians(df_nn, kept_sites)

# =============================================================================
# STEP 6: ASSIGN T:ET GROUPS BASED ON TERTILES
# =============================================================================
print("\n" + "="*80)
print("STEP 5: ASSIGNING T:ET GROUPS (TERTILES)")
print("="*80)
print("\n  EXPLANATION: How Low/Intermediate/High are defined:")
print("  ----------------------------------------------------")
print("  • Groups are based on TERTILES (33rd and 67th percentiles)")
print("  • Low T:ET:     sites below the 33rd percentile (evaporation-dominated)")
print("  • Intermediate: sites between 33rd and 67th percentiles")
print("  • High T:ET:    sites above the 67th percentile (transpiration-dominated)")
print("  • This ensures equal group sizes (approximately 1/3 of sites each)")
print("  • Threshold values are data-driven (percentiles, not fixed numbers)")
print("="*80)

site_medians = site_medians.sort_values("median_Trans_ratio").reset_index(drop=True)
q1 = site_medians["median_Trans_ratio"].quantile(1/3)
q2 = site_medians["median_Trans_ratio"].quantile(2/3)

def assign_group(x):
    if x <= q1:
        return "Low T:ET (evap-dominated)"
    elif x <= q2:
        return "Intermediate T:ET"
    else:
        return "High T:ET (transp-dominated)"

site_medians["TET_group"] = site_medians["median_Trans_ratio"].apply(assign_group)

print(f"\n  T:ET tertile thresholds:")
print(f"    33rd percentile (Low/Intermediate boundary): {q1:.3f}")
print(f"    67th percentile (Intermediate/High boundary): {q2:.3f}")

print(f"\n  T:ET group sizes:")
for group in ["Low T:ET (evap-dominated)", "Intermediate T:ET", "High T:ET (transp-dominated)"]:
    count = len(site_medians[site_medians["TET_group"] == group])
    print(f"    {group}: {count} sites")

# =============================================================================
# STEP 7: SAVE OUTPUTS (SAME NAMES - no downstream impact)
# =============================================================================
print("\n" + "="*80)
print("STEP 6: SAVING OUTPUTS")
print("="*80)

# Save site-level CSV
site_medians.to_csv(OUTPUT_SITE_CSV, index=False)
print(f"  Saved: {OUTPUT_SITE_CSV}")

# Summary table by T:ET group
def summarize_group(df_group):
    return pd.Series({
        "n_sites": len(df_group),
        "WUE_ET_median": df_group["median_WUE"].median(),
        "WUE_ET_IQR": df_group["median_WUE"].quantile(0.75) - df_group["median_WUE"].quantile(0.25),
        "WUE_E_median": df_group["median_WUE_eva"].median(),
        "WUE_E_IQR": df_group["median_WUE_eva"].quantile(0.75) - df_group["median_WUE_eva"].quantile(0.25),
        "WUE_T_median": df_group["median_WUE_tra"].median(),
        "WUE_T_IQR": df_group["median_WUE_tra"].quantile(0.75) - df_group["median_WUE_tra"].quantile(0.25),
    })

summary = site_medians.groupby("TET_group", group_keys=False).apply(summarize_group).reset_index()
summary.to_csv(OUTPUT_GROUP_SUMMARY, index=False)
print(f"  Saved: {OUTPUT_GROUP_SUMMARY}")

# =============================================================================
# STEP 8: SPEARMAN CORRELATIONS
# =============================================================================
print("\n" + "="*80)
print("STEP 7: SPEARMAN CORRELATIONS")
print("="*80)

corr_WUE, p_WUE = spearmanr(site_medians["median_Trans_ratio"], site_medians["median_WUE"])
corr_WUE_eva, p_eva = spearmanr(site_medians["median_Trans_ratio"], site_medians["median_WUE_eva"])
corr_WUE_tra, p_tra = spearmanr(site_medians["median_Trans_ratio"], site_medians["median_WUE_tra"])

print(f"  WUE_ET : rho = {corr_WUE:.3f}, p = {p_WUE:.4f}")
print(f"  WUE_E  : rho = {corr_WUE_eva:.3f}, p = {p_eva:.4f}")
print(f"  WUE_T  : rho = {corr_WUE_tra:.3f}, p = {p_tra:.4f}")

# =============================================================================
# STEP 9: CONVERGENCE METRICS
# =============================================================================
print("\n" + "="*80)
print("STEP 8: CONVERGENCE METRICS")
print("="*80)

# Calculate absolute difference for each site
site_medians["abs_diff_WUE_ET_T"] = np.abs(site_medians["median_WUE"] - site_medians["median_WUE_tra"])

# Calculate group-wise medians for annotation
low_group = site_medians[site_medians["TET_group"] == "Low T:ET (evap-dominated)"]
high_group = site_medians[site_medians["TET_group"] == "High T:ET (transp-dominated)"]

# Get median absolute differences (rounded to 1 significant digit)
diff_low_median = low_group["abs_diff_WUE_ET_T"].median() if len(low_group) > 0 else np.nan
diff_high_median = high_group["abs_diff_WUE_ET_T"].median() if len(high_group) > 0 else np.nan

# Round to 1 significant digit
def round_to_1_sigfig(x):
    if np.isnan(x):
        return np.nan
    return round(x, -int(np.floor(np.log10(abs(x)))) + 1)

diff_low_rounded = round_to_1_sigfig(diff_low_median)
diff_high_rounded = round_to_1_sigfig(diff_high_median)

print(f"  Median |WUE_ET - WUE_T| in low T:ET group: {diff_low_median:.3f} (rounded to {diff_low_rounded})" if not np.isnan(diff_low_median) else "  Low group: insufficient data")
print(f"  Median |WUE_ET - WUE_T| in high T:ET group: {diff_high_median:.3f} (rounded to {diff_high_rounded})" if not np.isnan(diff_high_median) else "  High group: insufficient data")

# IQR comparison
wue_et_iqr = site_medians["median_WUE"].quantile(0.75) - site_medians["median_WUE"].quantile(0.25)
wue_t_iqr = site_medians["median_WUE_tra"].quantile(0.75) - site_medians["median_WUE_tra"].quantile(0.25)
print(f"  IQR of WUE_ET: {wue_et_iqr:.3f}")
print(f"  IQR of WUE_T: {wue_t_iqr:.3f}")

# =============================================================================
# STEP 10: ROBUST THRESHOLD DETECTION (Theil-Sen)
# =============================================================================
print("\n" + "="*80)
print("STEP 9: ROBUST THRESHOLD DETECTION (Theil-Sen segmented regression)")
print("="*80)

x = site_medians["median_Trans_ratio"].values
y = site_medians["abs_diff_WUE_ET_T"].values

valid = np.isfinite(x) & np.isfinite(y)
x = x[valid]
y = y[valid]

min_sites_per_side = 8

x_low = np.quantile(x, 0.15)
x_high = np.quantile(x, 0.85)

candidate_breaks = np.linspace(x_low, x_high, 100)

best_score = np.inf
best_break = np.nan
best_models = {}

for bp in candidate_breaks:
    left_mask = x < bp
    right_mask = x >= bp

    n_left = np.sum(left_mask)
    n_right = np.sum(right_mask)

    if n_left < min_sites_per_side or n_right < min_sites_per_side:
        continue

    x_left, y_left = x[left_mask], y[left_mask]
    x_right, y_right = x[right_mask], y[right_mask]

    slope_left, intercept_left, _, _ = theilslopes(y_left, x_left)
    slope_right, intercept_right, _, _ = theilslopes(y_right, x_right)

    pred_left = slope_left * x_left + intercept_left
    pred_right = slope_right * x_right + intercept_right

    # Robust loss: sum of absolute residuals
    score = np.sum(np.abs(y_left - pred_left)) + np.sum(np.abs(y_right - pred_right))

    if score < best_score:
        best_score = score
        best_break = bp
        best_models = {
            "slope_low": slope_left,
            "intercept_low": intercept_left,
            "slope_high": slope_right,
            "intercept_high": intercept_right,
            "n_left": n_left,
            "n_right": n_right,
            "score": score
        }

if np.isfinite(best_break):
    print(f"\n  Robust estimated T:ET threshold = {best_break:.3f}")
    print(f"  Sites below threshold: {best_models['n_left']}")
    print(f"  Sites at/above threshold: {best_models['n_right']}")
    print(f"  Robust slope (T:ET < {best_break:.3f}): {best_models['slope_low']:.4f}")
    print(f"  Robust slope (T:ET ≥ {best_break:.3f}): {best_models['slope_high']:.4f}")

    if best_models["n_left"] < 10 or best_break <= np.quantile(x, 0.20):
        threshold_interpretation = "boundary-like behavior; interpret cautiously"
    elif np.sign(best_models["slope_low"]) != np.sign(best_models["slope_high"]):
        threshold_interpretation = "threshold-like behavior with a change in slope direction"
    else:
        threshold_interpretation = "threshold-like behavior with weaker evidence"

    print(f"  Interpretation: {threshold_interpretation}")
else:
    print("\n  No robust threshold could be estimated with the current minimum site requirement.")
    threshold_interpretation = "no robust threshold detected"

# =============================================================================
# STEP 11: CREATE FIGURE
# =============================================================================
print("\n" + "="*80)
print("STEP 10: CREATING FIGURE")
print("="*80)

fig, ax = plt.subplots(1, 1, figsize=(15, 7))
fig.subplots_adjust(top=0.78) 
# Sort for plotting
plot_df = site_medians.sort_values("median_Trans_ratio").reset_index(drop=True)

# =============================================================================
# UPDATED COLOR GRADIENT (vivid purple → blue-purple → muted brown-gold → golden yellow)
# =============================================================================
# Top (high T:ET):    #5B1CFF (vivid purple)
# Middle:             #4A36D8 (blue-purple)
# Lower:              #9B8A71 (muted brown-gold)
# Bottom (low T:ET):  #D7A900 (golden yellow)
tet_cmap = LinearSegmentedColormap.from_list(
    "arrow_gradient_grouped",
    [
        (0.0, "#D7A900"),   # Bottom (low T:ET) - golden yellow
        (0.33, "#9B8A71"),  # Lower - muted brown-gold
        (0.66, "#0066FF"),  # Middle - blue-purple
        (1.0, "#5B1CFF")    # Top (high T:ET) - vivid purple
    ]
)
norm = Normalize(vmin=plot_df["median_Trans_ratio"].min(), vmax=plot_df["median_Trans_ratio"].max())

# Draw vertical connecting lines for each site
for idx, row in plot_df.iterrows():
    tet = row["median_Trans_ratio"]
    wue_et = row["median_WUE"]
    wue_t = row["median_WUE_tra"]
    
    color = tet_cmap(norm(tet))
    ax.plot([tet, tet], [wue_et, wue_t], color=color, linewidth=1.8, alpha=0.75)

# Add WUE_ET points (open circles)
ax.scatter(plot_df["median_Trans_ratio"], plot_df["median_WUE"], 
           facecolor='none', edgecolor='black', s=42, linewidth=1.0, 
           label='WUE$_{ET}$', zorder=3)

# Add WUE_T points (filled circles)
ax.scatter(plot_df["median_Trans_ratio"], plot_df["median_WUE_tra"], 
           facecolor='black', edgecolor='black', s=38, linewidth=0.8, 
           label='WUE$_T$', alpha=0.75, zorder=3)

# =============================================================================
# ADD ANNOTATION LABELS AT TOP LEFT AND TOP RIGHT
# =============================================================================

# Row 1 titles
fig.text(0.15, 1, "Low T:ET", 
         ha='center', va='bottom', fontsize=16, fontweight='bold', color='#4a4a4a')
fig.text(0.55, 1, "High T:ET", 
         ha='center', va='bottom', fontsize=16, fontweight='bold', color='#4a4a4a')

# Row 2 values with 1 significant digit
fig.text(0.15, 0.96, f"|WUE$_{{ET}}$ − WUE$_T$| = {diff_low_rounded}", 
         ha='center', va='bottom', fontsize=14, color='#4a4a4a')
fig.text(0.55, 0.96, f"|WUE$_{{ET}}$ − WUE$_T$| = {diff_high_rounded}", 
         ha='center', va='bottom', fontsize=14, color='#4a4a4a')

ax.set_xlabel("Median T:ET ratio", fontsize=22.5)
ax.set_ylabel("Median WUE (g C kg$^{-1}$ H$_2$O$^{-1}$)", fontsize=22.5)

# Add colorbar (automatically uses tet_cmap with updated colors)
sm = ScalarMappable(norm=norm, cmap=tet_cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=30, pad=0.02)
cbar.set_label('T:ET ratio', fontsize=19.5)
cbar.ax.tick_params(labelsize=18)

# Legend
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='none', markerfacecolor='none', 
               markeredgecolor='black', markersize=10, linewidth=1.0, label='WUE$_{ET}$'),
    plt.Line2D([0], [0], marker='o', color='none', markerfacecolor='black', 
               markeredgecolor='black', markersize=10, label='WUE$_T$'),
    plt.Line2D([0], [0], color='gray', linewidth=1.8, alpha=0.75, label='|WUE$_{ET}$ − WUE$_T$|'),
]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
          frameon=True, fontsize=17)

ax.grid(True, linestyle=':', alpha=0.15, axis='both')

plt.tight_layout(rect=[0, 0, 0.82, 1])

# Save figure (SAME NAMES)
plt.savefig(OUTPUT_FIGURE_PNG, dpi=600, bbox_inches='tight')
plt.savefig(OUTPUT_FIGURE_PDF, dpi=600, bbox_inches='tight')
print(f"\n  Saved figure: {OUTPUT_FIGURE_PNG}")
print(f"  Saved figure: {OUTPUT_FIGURE_PDF}")

plt.show()

# =============================================================================
# STEP 12: PRINT FINAL SUMMARY
# =============================================================================
print("\n" + "="*80)
print("FINAL SUMMARY - T:ET GRADIENT ANALYSIS")
print(f"Using SAME filtering criteria as Figure 1")
print("="*80)

print(f"\n1. FILTERING SUMMARY:")
print(f"   Strict month filter (upstream CHUNK 3): months where WUE, WUE_E, WUE_T ALL have data")
print(f"   Strict triple intersection sites: {len(kept_sites)} (SAME as Figure 1)")
print(f"   Global {MIN_NN_MONTHS}-observation NN filter applied")
print(f"   NN rows used: {len(df_nn):,}")
print(f"   Retained sites (all metrics valid): {len(site_medians)}")

print(f"\n2. SPEARMAN CORRELATIONS (Trans_ratio vs WUE metrics):")
print(f"   WUE_ET : rho = {corr_WUE:.3f}, p = {p_WUE:.4f}")
print(f"   WUE_E  : rho = {corr_WUE_eva:.3f}, p = {p_eva:.4f}")
print(f"   WUE_T  : rho = {corr_WUE_tra:.3f}, p = {p_tra:.4f}")

print(f"\n3. CONVERGENCE CHECK:")
if not np.isnan(diff_low_median) and not np.isnan(diff_high_median):
    print(f"   |WUE_ET - WUE_T| at low T:ET: {diff_low_median:.3f} → {diff_low_rounded} (1 sig fig)")
    print(f"   |WUE_ET - WUE_T| at high T:ET: {diff_high_median:.3f} → {diff_high_rounded} (1 sig fig)")
    if diff_high_median < diff_low_median:
        print(f"   ✓ WUE_ET and WUE_T CONVERGE at high T:ET")
    else:
        print(f"   ✗ WUE_ET and WUE_T do NOT converge at high T:ET")

print(f"\n4. ROBUST THRESHOLD:")
if np.isfinite(best_break):
    print(f"   Estimated T:ET threshold = {best_break:.3f}")
    print(f"   Interpretation: {threshold_interpretation}")

print("\n" + "="*80)
print("✅ ANALYSIS COMPLETE - Using SAME filtering as Figure 1")
print("="*80)