"""
Q2_WUE_performance.py - EXACT MALONE REPLICATION
Malone-style Q2 WUE_T performance metrics with Python

OUTPUT STRUCTURE (separate from Malone):
  Q2_WUE_performance_outputs/   - CSV files, summary text
  Q2_WUE_performance_figures/   - PNG figure files
  Q2_WUE_performance_temp/      - Temporary files (optional, for future use)

EXACT MALONE OUTPUT NAMES:
  CSV Files:
    Q2_site_stability.csv
    Q2_site_plasticity_slope.csv
    Q2_site_plasticity_range.csv
    Q2_site_resistance.csv
    Q2_site_recovery_events.csv
    Q2_site_recovery.csv
    Q2_site_performance_summary.csv
    Q2_ecosystem_class_statistical_tests.csv
    Q2_WUE_performance_summary.txt

  Figures:
    Q2_plot_01_stability.png
    Q2_plot_02_plasticity_slope.png
    Q2_plot_03_plasticity_range.png
    Q2_plot_04_resistance.png
    Q2_plot_05_recovery.png
    Q2_plot_06_performance_heatmap.png
    Q2_plot_07_composite.png
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from scipy.stats import false_discovery_control
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("Q2: WUE_T Performance Metrics - EXACT MALONE REPLICATION")
print("="*60)

# ============================================================================
# PATHS - SEPARATE FROM MALONE
# ============================================================================

# IMPORTANT: Use the SAME input file as Malone
input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"


# Python outputs - SEPARATE from Malone with clean organization
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2"

# Three separate subdirectories for organization
output_dir = os.path.join(base_output_dir, "Q2_WUE_performance_outputs")     # CSV, summary
figure_dir = os.path.join(base_output_dir, "Q2_WUE_performance_figures")     # PNG figures
temp_dir = os.path.join(base_output_dir, "Q2_WUE_performance_temp")          # Temporary files (optional)

# Create all directories
for dir_path in [output_dir, figure_dir, temp_dir]:
    os.makedirs(dir_path, exist_ok=True)

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")
print(f"(Separate from Malone outputs - safe for comparison)")

# ============================================================================
# CLEAN OLD OUTPUTS AT START
# ============================================================================

print("\n" + "="*60)
print("STEP 0: Cleaning old outputs (Python directories only)")
print("="*60)

expected_outputs = [
    "Q2_site_stability.csv",
    "Q2_site_plasticity_slope.csv",
    "Q2_site_plasticity_range.csv",
    "Q2_site_resistance.csv",
    "Q2_site_recovery_events.csv",
    "Q2_site_recovery.csv",
    "Q2_site_performance_summary.csv",
    "Q2_ecosystem_class_statistical_tests.csv",
    "Q2_WUE_performance_summary.txt",
]

expected_figures = [
    "Q2_plot_01_stability.png",
    "Q2_plot_02_plasticity_slope.png",
    "Q2_plot_03_plasticity_range.png",
    "Q2_plot_04_resistance.png",
    "Q2_plot_05_recovery.png",
    "Q2_plot_06_performance_heatmap.png",
    "Q2_plot_07_composite.png",
]

# Clean output files
for f in expected_outputs:
    path = os.path.join(output_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old output: {f}")

# Clean figure files
for f in expected_figures:
    path = os.path.join(figure_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old figure: {f}")

# Clean temp directory
for f in os.listdir(temp_dir):
    path = os.path.join(temp_dir, f)
    if os.path.isfile(path):
        os.remove(path)
        print(f"  Removed old temp file: {f}")

# ============================================================================
# CONSTANTS - EXACT MALONE
# ============================================================================

SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
SINGLE_SPEI = "SPEI_3"
SINGLE_SPEI_LABEL = "SPEI-3"
ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]

# EXACT Malone colors - CONSISTENT across ALL plots
ECOSYSTEM_COLORS = {
    "Upland": "#800080",     # Purple
    "Freshwater": "#0000FF", # Blue
    "Saline": "#FFA500"      # Orange - CONSISTENT with Malone
}

def get_ecosystem_color(ecosystem):
    """Return the correct color for each ecosystem"""
    return ECOSYSTEM_COLORS.get(ecosystem, "#808080")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def identify_drought_events(df, spei_col="SPEI_3", threshold=-1, min_duration=2):
    """Identify consecutive drought events per site - EXACT MALONE METHOD"""
    df = df.sort_values(['site_name', 'Year', 'month']).copy()
    events = []
    
    for site_name, site_df in df.groupby('site_name'):
        site_df = site_df.sort_values(['Year', 'month']).reset_index(drop=True)
        spei_values = site_df[spei_col].values
        in_drought = (spei_values < threshold) & np.isfinite(spei_values)
        
        i = 0
        event_id = 1
        while i < len(in_drought):
            if in_drought[i]:
                start = i
                while i < len(in_drought) and in_drought[i]:
                    i += 1
                end = i - 1
                duration = end - start + 1
                if duration >= min_duration:
                    events.append({
                        'site_name': site_name,
                        'event_id': event_id,
                        'start_idx': start,
                        'end_idx': end,
                        'start_year': site_df.loc[start, 'Year'],
                        'start_month': site_df.loc[start, 'month'],
                        'end_year': site_df.loc[end, 'Year'],
                        'end_month': site_df.loc[end, 'month'],
                        'duration_months': duration
                    })
                    event_id += 1
            else:
                i += 1
    
    if not events:
        return pd.DataFrame()
    return pd.DataFrame(events)

def sig(p):
    """Significance stars - EXACT MALONE"""
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.1:
        return "."
    return ""

def ecosystem_metric_test(df, metric_col, metric_label):
    """Perform Kruskal-Wallis test with pairwise Wilcoxon - EXACT MALONE"""
    test_df = df[df[metric_col].notna() & df['water_class'].notna()].copy()
    test_df = test_df[['water_class', metric_col]].dropna()
    test_df['water_class'] = pd.Categorical(test_df['water_class'], categories=ECOSYSTEM_CLASSES)
    test_df = test_df.dropna(subset=['water_class'])
    
    group_counts = test_df.groupby('water_class').size().to_dict()
    
    result = {
        'metric': metric_label,
        'response_variable': metric_col,
        'test': 'Kruskal-Wallis',
        'n_sites': len(test_df),
        'n_upland': group_counts.get('Upland', 0),
        'n_freshwater': group_counts.get('Freshwater', 0),
        'n_saline': group_counts.get('Saline', 0),
        'statistic': np.nan,
        'df': np.nan,
        'p_value': np.nan,
        'p_adjust_method': 'BH/FDR for pairwise Wilcoxon',
        'pairwise_upland_vs_freshwater_p_adj': np.nan,
        'pairwise_upland_vs_saline_p_adj': np.nan,
        'pairwise_freshwater_vs_saline_p_adj': np.nan,
    }
    
    if len(test_df) >= 3 and len(test_df['water_class'].unique()) >= 2:
        groups = [group[metric_col].values for name, group in test_df.groupby('water_class') if len(group) >= 2]
        if len(groups) >= 2:
            h_stat, p_val = kruskal(*groups)
            result['statistic'] = h_stat
            result['df'] = len(groups) - 1
            result['p_value'] = p_val
            
            # Pairwise Wilcoxon with BH correction
            classes = test_df['water_class'].unique()
            pair_p_values = []
            pair_names = []
            for i, c1 in enumerate(classes):
                for c2 in classes[i+1:]:
                    data1 = test_df[test_df['water_class'] == c1][metric_col].values
                    data2 = test_df[test_df['water_class'] == c2][metric_col].values
                    if len(data1) >= 2 and len(data2) >= 2:
                        _, p_pair = mannwhitneyu(data1, data2, alternative='two-sided')
                        pair_p_values.append(p_pair)
                        pair_names.append(f'{c1}_vs_{c2}')
            
            # Apply BH correction
            if len(pair_p_values) > 0:
                p_adj = false_discovery_control(pair_p_values, method='bh')
                
                for name, p_adj_val in zip(pair_names, p_adj):
                    if name == 'Upland_vs_Freshwater' or name == 'Freshwater_vs_Upland':
                        result['pairwise_upland_vs_freshwater_p_adj'] = p_adj_val
                    elif name == 'Upland_vs_Saline' or name == 'Saline_vs_Upland':
                        result['pairwise_upland_vs_saline_p_adj'] = p_adj_val
                    elif name == 'Freshwater_vs_Saline' or name == 'Saline_vs_Freshwater':
                        result['pairwise_freshwater_vs_saline_p_adj'] = p_adj_val
    
    result['signif'] = sig(result['p_value'])
    return pd.DataFrame([result])

def theme_wue(base_size=13):
    """Malone-style theme"""
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
    })

def save_plot(fig, filename, width=11, height=7, dpi=300):
    """Save figure to figure_dir - EXACT MALONE"""
    fig.set_size_inches(width, height)
    fig.savefig(os.path.join(figure_dir, filename), dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def write_table(df, filename):
    """Write CSV to output_dir - EXACT MALONE"""
    df.to_csv(os.path.join(output_dir, filename), index=False)

def capture_to_file(content, filename):
    """Write summary text - EXACT MALONE"""
    with open(os.path.join(output_dir, filename), 'w') as f:
        f.write(content)

# ============================================================================
# DATA LOADING - EXACT MALONE FILTERING
# ============================================================================

print("\n" + "="*60)
print("STEP 1: Data Loading and Filtering")
print("="*60)

print(f"\nLoading data from: {input_file}")
monthly = pd.read_csv(input_file)

# Trim whitespace - exact Malone
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].astype(str).str.strip()
    monthly[col] = monthly[col].replace('nan', np.nan)

# Filter - exact Malone
base_df = monthly[
    monthly['site_name'].notna() &
    (monthly['site_name'] != "") &
    monthly['water_class'].notna() &
    (monthly['water_class'] != "") &
    monthly['water_class'].isin(ECOSYSTEM_CLASSES) &
    monthly['lat'].notna() &
    np.isfinite(monthly['lat']) &
    monthly['long'].notna() &
    np.isfinite(monthly['long']) &
    monthly['WUE_tra'].notna() &
    np.isfinite(monthly['WUE_tra']) &
    monthly['Trans_ratio'].notna() &
    np.isfinite(monthly['Trans_ratio']) &
    (monthly['Trans_ratio'] >= 0) &
    (monthly['Trans_ratio'] <= 1)
].copy()

base_df['WUE_T'] = base_df['WUE_tra']
base_df['water_class'] = pd.Categorical(base_df['water_class'], categories=ECOSYSTEM_CLASSES)
base_df = base_df.sort_values(['site_name', 'Year', 'month']).reset_index(drop=True)

print(f"\nData loaded:")
print(f"  Rows: {len(base_df):,}")
print(f"  Sites: {base_df['site_name'].nunique():,}")
print(f"  Years: {base_df['Year'].min()} - {base_df['Year'].max()}")

# ============================================================================
# 1. STABILITY - EXACT MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 2: Calculating Stability")
print("="*60)

stability = base_df.groupby(['site_name', 'water_class']).agg(
    n_months=('WUE_T', 'size'),
    mean_TET=('Trans_ratio', 'mean'),
    mean_WUE_T=('WUE_T', 'mean'),
    var_WUE_T=('WUE_T', 'var'),
    sd_WUE_T=('WUE_T', 'std')
).reset_index()

stability = stability[stability['n_months'] >= 6]
stability['cv_WUE_T'] = stability['sd_WUE_T'] / stability['mean_WUE_T']
stability['stability'] = stability['mean_WUE_T'] / stability['var_WUE_T']

write_table(stability, "Q2_site_stability.csv")
print(f"  Saved: {len(stability)} sites")

# ============================================================================
# 2. PLASTICITY SLOPE - EXACT MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 3: Calculating Plasticity Slopes")
print("="*60)

plasticity_slope_list = []

for spei_col in SPEI_COLS:
    temp_df = base_df[base_df[spei_col].notna() & np.isfinite(base_df[spei_col])].copy()
    temp_df = temp_df[['site_name', 'water_class', 'WUE_T', spei_col]].dropna()
    temp_df.rename(columns={spei_col: 'SPEI_value'}, inplace=True)
    
    for (site_name, water_class), group in temp_df.groupby(['site_name', 'water_class']):
        if len(group) >= 6 and group['SPEI_value'].nunique() >= 4:
            x = group['SPEI_value'].values
            y = group['WUE_T'].values
            slope, intercept = np.polyfit(x, y, 1)
            y_pred = slope * x + intercept
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else np.nan
            correlation = np.corrcoef(x, y)[0, 1] if len(x) > 1 else np.nan
            
            plasticity_slope_list.append({
                'site_name': site_name,
                'water_class': water_class,
                'SPEI_timescale': spei_col,
                'n_months': len(group),
                'slope': slope,
                'r_squared': r_squared,
                'correlation': correlation
            })

plasticity_slope = pd.DataFrame(plasticity_slope_list)
if len(plasticity_slope) > 0:
    plasticity_slope['SPEI_timescale'] = pd.Categorical(plasticity_slope['SPEI_timescale'], categories=SPEI_COLS)
    write_table(plasticity_slope, "Q2_site_plasticity_slope.csv")
print(f"  Saved: {len(plasticity_slope)} site × timescale rows")

# ============================================================================
# 3. PLASTICITY RANGE - EXACT MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 4: Calculating Plasticity Range")
print("="*60)

plasticity_range = base_df.groupby(['site_name', 'water_class']).agg(
    n_months=('WUE_T', 'size'),
    WUE_T_max=('WUE_T', 'max'),
    WUE_T_min=('WUE_T', 'min'),
    WUE_T_p95=('WUE_T', lambda x: np.percentile(x, 95) if len(x) >= 2 else np.nan),
    WUE_T_p05=('WUE_T', lambda x: np.percentile(x, 5) if len(x) >= 2 else np.nan)
).reset_index()

plasticity_range = plasticity_range[
    (plasticity_range['n_months'] >= 6) &
    (plasticity_range['WUE_T_min'] > 0)
]

plasticity_range['plasticity_range'] = plasticity_range['WUE_T_max'] / plasticity_range['WUE_T_min']
plasticity_range['plasticity_p95_p05'] = plasticity_range['WUE_T_p95'] / plasticity_range['WUE_T_p05']

write_table(plasticity_range, "Q2_site_plasticity_range.csv")
print(f"  Saved: {len(plasticity_range)} sites")

# ============================================================================
# 4. RESISTANCE - EXACT MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 5: Calculating Resistance")
print("="*60)

resistance_df = base_df[base_df[SINGLE_SPEI].notna() & np.isfinite(base_df[SINGLE_SPEI])].copy()
resistance_df['drought_class'] = np.where(
    resistance_df[SINGLE_SPEI] < -1, 'drought',
    np.where(resistance_df[SINGLE_SPEI] <= 1, 'near_normal', 'wet')
)

resistance_wide = resistance_df[
    resistance_df['drought_class'].isin(['drought', 'near_normal'])
].groupby(['site_name', 'water_class', 'drought_class']).agg(
    n_months=('WUE_T', 'size'),
    mean_WUE_T=('WUE_T', 'mean')
).reset_index()

resistance_pivot = resistance_wide.pivot_table(
    index=['site_name', 'water_class'],
    columns='drought_class',
    values=['n_months', 'mean_WUE_T']
).reset_index()

resistance_pivot.columns = ['_'.join(col).strip() if col[1] else col[0] for col in resistance_pivot.columns.values]
resistance_pivot = resistance_pivot.rename(columns={'site_name_': 'site_name', 'water_class_': 'water_class'})

if 'n_months_drought' in resistance_pivot.columns and 'n_months_near_normal' in resistance_pivot.columns:
    resistance = resistance_pivot[
        (resistance_pivot['n_months_drought'] >= 3) &
        (resistance_pivot['n_months_near_normal'] >= 3) &
        (resistance_pivot['mean_WUE_T_near_normal'] > 0)
    ].copy()
    
    resistance['resistance'] = resistance['mean_WUE_T_drought'] / resistance['mean_WUE_T_near_normal']
    write_table(resistance, "Q2_site_resistance.csv")
    print(f"  Saved: {len(resistance)} sites")
else:
    resistance = pd.DataFrame()
    print("  Warning: Could not calculate resistance")

# ============================================================================
# 5. RECOVERY - EXACT MALONE (FIXED INDEXING)
# ============================================================================

print("\n" + "="*60)
print("STEP 6: Calculating Recovery")
print("="*60)

# Identify drought events
events = identify_drought_events(
    base_df[base_df[SINGLE_SPEI].notna() & np.isfinite(base_df[SINGLE_SPEI])],
    spei_col=SINGLE_SPEI,
    threshold=-1,
    min_duration=2
)
print(f"  Drought events identified: {len(events)}")

WINDOW_MONTHS = 12
MIN_WINDOW_N = 3
recovery_list = []

if len(events) > 0:
    for idx, event in events.iterrows():
        site_df = base_df[base_df['site_name'] == event['site_name']].sort_values(['Year', 'month']).reset_index(drop=True)
        n = len(site_df)
        
        start_idx = int(event['start_idx'])
        end_idx = int(event['end_idx'])
        
        # Pre-drought: 12 months BEFORE event start (EXACT MALONE)
        pre_end = start_idx - 1
        pre_start = max(0, pre_end - WINDOW_MONTHS + 1)
        
        # Post-drought: 12 months AFTER event end (EXACT MALONE)
        post_start = end_idx + 1
        post_end = min(n, post_start + WINDOW_MONTHS)
        
        # Check windows exist
        if pre_end < pre_start or post_start >= n:
            continue
        
        # Extract WUE_T values
        pre_wue = site_df.iloc[pre_start:pre_end + 1]['WUE_T'].values
        post_wue = site_df.iloc[post_start:post_end]['WUE_T'].values
        
        pre_wue = pre_wue[np.isfinite(pre_wue)]
        post_wue = post_wue[np.isfinite(post_wue)]
        
        if len(pre_wue) < MIN_WINDOW_N or len(post_wue) < MIN_WINDOW_N:
            continue
        
        mean_pre = np.mean(pre_wue)
        mean_post = np.mean(post_wue)
        
        if mean_pre <= 0:
            continue
        
        recovery_list.append({
            'site_name': event['site_name'],
            'water_class': site_df.loc[0, 'water_class'],
            'event_id': event['event_id'],
            'start_year': event['start_year'],
            'start_month': event['start_month'],
            'end_year': event['end_year'],
            'end_month': event['end_month'],
            'duration_months': event['duration_months'],
            'n_pre': len(pre_wue),
            'n_post': len(post_wue),
            'mean_WUE_T_pre': mean_pre,
            'mean_WUE_T_post': mean_post,
            'recovery': mean_post / mean_pre
        })

recovery_events = pd.DataFrame(recovery_list)

if len(recovery_events) > 0:
    write_table(recovery_events, "Q2_site_recovery_events.csv")
    
    recovery_site = recovery_events.groupby(['site_name', 'water_class']).agg(
        n_events=('event_id', 'size'),
        mean_recovery=('recovery', 'mean'),
        median_recovery=('recovery', 'median'),
        sd_recovery=('recovery', 'std'),
        min_recovery=('recovery', 'min'),
        max_recovery=('recovery', 'max'),
        pct_full_recovery=('recovery', lambda x: np.mean(x >= 1.0) * 100 if len(x) > 0 else 0)
    ).reset_index()
    
    write_table(recovery_site, "Q2_site_recovery.csv")
    print(f"  Saved: {len(recovery_events)} events across {len(recovery_site)} sites")
else:
    recovery_site = pd.DataFrame()
    print("  No recovery events found")

# ============================================================================
# 6. SUMMARY TABLE - EXACT MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 7: Creating Summary Table")
print("="*60)

plasticity_slope_3 = plasticity_slope[plasticity_slope['SPEI_timescale'] == SINGLE_SPEI][
    ['site_name', 'slope', 'r_squared']
].rename(columns={'slope': 'plasticity_slope_SPEI3', 'r_squared': 'plasticity_r2_SPEI3'})

summary_table = stability[['site_name', 'water_class', 'n_months', 'mean_TET', 
                           'mean_WUE_T', 'var_WUE_T', 'cv_WUE_T', 'stability']].copy()

if len(plasticity_slope_3) > 0:
    summary_table = summary_table.merge(plasticity_slope_3, on='site_name', how='left')

if len(plasticity_range) > 0:
    summary_table = summary_table.merge(
        plasticity_range[['site_name', 'plasticity_range', 'plasticity_p95_p05', 
                          'WUE_T_max', 'WUE_T_min']],
        on='site_name', how='left'
    )

if len(resistance) > 0:
    summary_table = summary_table.merge(
        resistance[['site_name', 'resistance', 'mean_WUE_T_drought', 'mean_WUE_T_near_normal']],
        on='site_name', how='left'
    )

if len(recovery_site) > 0:
    summary_table = summary_table.merge(
        recovery_site[['site_name', 'n_events', 'mean_recovery', 'pct_full_recovery']],
        on='site_name', how='left'
    )

write_table(summary_table, "Q2_site_performance_summary.csv")
print(f"  Saved: {len(summary_table)} sites")

# ============================================================================
# 7. ECOSYSTEM-CLASS STATISTICAL TESTS - EXACT MALONE
# ============================================================================

print("\n" + "="*60)
print("STEP 8: Calculating Ecosystem-Class Statistical Tests")
print("="*60)

test_metrics = [
    ('stability', 'Stability index'),
    ('plasticity_slope_SPEI3', 'Plasticity slope (SPEI-3)'),
    ('plasticity_p95_p05', 'Plasticity range (95th/5th)'),
    ('resistance', 'Drought resistance'),
    ('mean_recovery', 'Drought recovery')
]

test_results = []
for metric_col, metric_label in test_metrics:
    if metric_col in summary_table.columns:
        result = ecosystem_metric_test(summary_table, metric_col, metric_label)
        test_results.append(result)

if test_results:
    ecosystem_tests = pd.concat(test_results, ignore_index=True)
    write_table(ecosystem_tests, "Q2_ecosystem_class_statistical_tests.csv")
    print(f"  Saved: {len(ecosystem_tests)} metrics")

# ============================================================================
# 8. CREATE FIGURES - EXACT MALONE NAMES AND STYLES
# ============================================================================

print("\n" + "="*60)
print("STEP 9: Creating Figures (Malone-style)")
print("="*60)

theme_wue()

# -------------------------------------------------------------------------
# FIGURE 1: Stability - Malone: Q2_plot_01_stability.png
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))
for ecosystem in ECOSYSTEM_CLASSES:
    subset = stability[stability['water_class'] == ecosystem]
    ax.scatter(subset['mean_TET'], subset['stability'], 
               color=ECOSYSTEM_COLORS[ecosystem], alpha=0.75, s=50, label=ecosystem)
ax.set_xlabel('Mean T:ET ratio')
ax.set_ylabel('Stability index: mean / variance (log scale)')
ax.set_xlim(0, 1)
ax.set_yscale('log')
ax.set_title('WUE_T Stability by Site', fontweight='bold')
ax.legend(loc='best', title='Ecosystem class')
ax.grid(True, alpha=0.3)
save_plot(fig, "Q2_plot_01_stability.png", width=10, height=6)
print("  Saved: Q2_plot_01_stability.png")

# -------------------------------------------------------------------------
# FIGURE 2: Plasticity slope - Malone: Q2_plot_02_plasticity_slope.png
# -------------------------------------------------------------------------
if len(plasticity_slope) > 0:
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x_positions = np.arange(len(SPEI_COLS))
    width = 0.25
    
    for i, ecosystem in enumerate(ECOSYSTEM_CLASSES):
        positions = x_positions + (i - 1) * width
        subset = plasticity_slope[plasticity_slope['water_class'] == ecosystem]
        data = [subset[subset['SPEI_timescale'] == s]['slope'].dropna().values for s in SPEI_COLS]
        
        bp = ax.boxplot(data, positions=positions, widths=width, patch_artist=True,
                        boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.5),
                        whiskerprops=dict(color='gray'),
                        capprops=dict(color='gray'),
                        medianprops=dict(color='black', linewidth=1.5))
    
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(SPEI_COLS, rotation=30, ha='right')
    ax.set_xlabel('SPEI timescale')
    ax.set_ylabel('Slope (g C mm⁻¹ per SPEI unit)')
    ax.set_title('WUE_T Plasticity: Slope of WUE_T ~ SPEI', fontweight='bold')
    
    legend_elements = [Patch(facecolor=ECOSYSTEM_COLORS[eco], alpha=0.5, label=eco) for eco in ECOSYSTEM_CLASSES]
    ax.legend(handles=legend_elements, loc='best', title='Ecosystem class')
    ax.grid(True, alpha=0.3)
    save_plot(fig, "Q2_plot_02_plasticity_slope.png", width=12, height=6)
    print("  Saved: Q2_plot_02_plasticity_slope.png")

# -------------------------------------------------------------------------
# FIGURE 3-5: Violin + Boxplot - Malone: Q2_plot_03/04/05
# -------------------------------------------------------------------------
def create_violin_boxplot(data, y_col, y_label, title, filename):
    fig, ax = plt.subplots(figsize=(8, 7))
    for i, ecosystem in enumerate(ECOSYSTEM_CLASSES):
        subset = data[data['water_class'] == ecosystem][y_col].dropna()
        if len(subset) > 0:
            parts = ax.violinplot(subset, positions=[i], widths=0.6, showmeans=False, showmedians=False)
            for pc in parts['bodies']:
                pc.set_facecolor(ECOSYSTEM_COLORS[ecosystem])
                pc.set_alpha(0.5)
            bp = ax.boxplot(subset, positions=[i], widths=0.18, patch_artist=True,
                           boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.7),
                           whiskerprops=dict(color='gray'),
                           capprops=dict(color='gray'),
                           medianprops=dict(color='black', linewidth=1.5))
    ax.set_xticks(range(len(ECOSYSTEM_CLASSES)))
    ax.set_xticklabels(ECOSYSTEM_CLASSES)
    ax.set_ylabel(y_label)
    ax.set_title(title, fontweight='bold')
    ax.grid(True, alpha=0.3)
    save_plot(fig, filename, width=8, height=7)
    print(f"  Saved: {filename}")

if len(plasticity_range) > 0:
    create_violin_boxplot(plasticity_range, 'plasticity_p95_p05', 
                         'WUE_T 95th / 5th percentile', 
                         'WUE_T Plasticity: Dynamic Range',
                         'Q2_plot_03_plasticity_range.png')

if len(resistance) > 0:
    create_violin_boxplot(resistance, 'resistance',
                         'Resistance (drought WUE_T / near-normal WUE_T)',
                         'WUE_T Resistance to Drought',
                         'Q2_plot_04_resistance.png')

if len(recovery_site) > 0:
    create_violin_boxplot(recovery_site, 'mean_recovery',
                         'Recovery (post-drought / pre-drought)',
                         'WUE_T Recovery Post-Drought',
                         'Q2_plot_05_recovery.png')

# -------------------------------------------------------------------------
# FIGURE 6: Heatmap - Malone: Q2_plot_06_performance_heatmap.png
# -------------------------------------------------------------------------
metrics_for_heatmap = ['stability', 'plasticity_slope_SPEI3', 'plasticity_range', 
                       'resistance', 'mean_recovery']
metric_labels = ['Stability index\n(mean/variance)', 'Plasticity\n(slope~SPEI-3)', 
                 'Plasticity\n(max/min)', 'Resistance\n(drought/normal)', 
                 'Recovery\n(post/pre)']

if all(col in summary_table.columns for col in metrics_for_heatmap):
    heatmap_df = summary_table[['site_name', 'water_class'] + metrics_for_heatmap].dropna()
    
    if len(heatmap_df) > 0:
        # Calculate z-scores
        for metric in metrics_for_heatmap:
            heatmap_df[f'{metric}_z'] = (heatmap_df[metric] - heatmap_df[metric].mean()) / heatmap_df[metric].std()
        
        # Create faceted heatmap by ecosystem
        fig, axes = plt.subplots(len(ECOSYSTEM_CLASSES), 1, figsize=(10, 14), squeeze=False)
        
        for idx, ecosystem in enumerate(ECOSYSTEM_CLASSES):
            eco_df = heatmap_df[heatmap_df['water_class'] == ecosystem].copy()
            if len(eco_df) > 0:
                # Sort by stability z-score
                eco_df = eco_df.sort_values('stability_z')
                
                # Prepare data for heatmap
                z_cols = [f'{m}_z' for m in metrics_for_heatmap]
                z_data = eco_df[z_cols].values
                
                # Plot on appropriate axis
                ax = axes[idx, 0]
                im = ax.imshow(z_data.T, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
                
                # Labels
                ax.set_xticks(range(len(eco_df)))
                ax.set_xticklabels(eco_df['site_name'], rotation=90, ha='center', fontsize=6)
                ax.set_yticks(range(len(metric_labels)))
                ax.set_yticklabels(metric_labels, fontsize=8)
                ax.set_title(f'{ecosystem} (n={len(eco_df)})', fontsize=12, fontweight='bold')
                ax.grid(False)
        
        # Add colorbar
        plt.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax, label='Z-score')
        
        fig.suptitle('WUE_T Performance Profile: All Sites', fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 0.85, 0.95])
        save_plot(fig, "Q2_plot_06_performance_heatmap.png", width=10, height=14)
        print("  Saved: Q2_plot_06_performance_heatmap.png")

# -------------------------------------------------------------------------
# FIGURE 7: Composite Panel - Malone: Q2_plot_07_composite.png
# -------------------------------------------------------------------------
fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1], hspace=0.3, wspace=0.3)

# Panel A: Stability
ax_a = fig.add_subplot(gs[0, 0])
for ecosystem in ECOSYSTEM_CLASSES:
    subset = stability[stability['water_class'] == ecosystem]
    ax_a.scatter(subset['mean_TET'], subset['stability'], 
                 color=ECOSYSTEM_COLORS[ecosystem], alpha=0.75, s=30)
ax_a.set_xlabel('T:ET ratio')
ax_a.set_ylabel('Stability index')
ax_a.set_yscale('log')
ax_a.set_title('A: Stability', fontweight='bold')
ax_a.grid(True, alpha=0.3)

# Panel B: Plasticity slope (ALL SPEI timescales - boxplot)
ax_b = fig.add_subplot(gs[0, 1])
x_positions = np.arange(len(SPEI_COLS))
width = 0.25
for i, ecosystem in enumerate(ECOSYSTEM_CLASSES):
    positions = x_positions + (i - 1) * width
    subset = plasticity_slope[plasticity_slope['water_class'] == ecosystem]
    data = [subset[subset['SPEI_timescale'] == s]['slope'].dropna().values for s in SPEI_COLS]
    bp = ax_b.boxplot(data, positions=positions, widths=width, patch_artist=True,
                     boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.5),
                     whiskerprops=dict(color='gray'),
                     capprops=dict(color='gray'),
                     medianprops=dict(color='black', linewidth=1.5))
ax_b.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
ax_b.set_xticks(x_positions)
ax_b.set_xticklabels(SPEI_COLS, rotation=30, ha='right', fontsize=8)
ax_b.set_ylabel('Slope')
ax_b.set_title('B: Plasticity slope', fontweight='bold')
ax_b.grid(True, alpha=0.3)

# Panel C: Plasticity range
ax_c = fig.add_subplot(gs[1, 0])
for i, ecosystem in enumerate(ECOSYSTEM_CLASSES):
    subset = plasticity_range[plasticity_range['water_class'] == ecosystem]['plasticity_p95_p05'].dropna().values
    if len(subset) > 0:
        parts = ax_c.violinplot(subset, positions=[i], widths=0.6, showmeans=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(ECOSYSTEM_COLORS[ecosystem])
            pc.set_alpha(0.5)
        bp = ax_c.boxplot(subset, positions=[i], widths=0.18, patch_artist=True,
                         boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.7),
                         whiskerprops=dict(color='gray'),
                         capprops=dict(color='gray'),
                         medianprops=dict(color='black', linewidth=1.5))
ax_c.set_xticks(range(len(ECOSYSTEM_CLASSES)))
ax_c.set_xticklabels(ECOSYSTEM_CLASSES, rotation=30, ha='right', fontsize=8)
ax_c.set_ylabel('95th / 5th')
ax_c.set_title('C: Plasticity range', fontweight='bold')
ax_c.grid(True, alpha=0.3)

# Panel D: Resistance
ax_d = fig.add_subplot(gs[1, 1])
for i, ecosystem in enumerate(ECOSYSTEM_CLASSES):
    subset = resistance[resistance['water_class'] == ecosystem]['resistance'].dropna().values
    if len(subset) > 0:
        parts = ax_d.violinplot(subset, positions=[i], widths=0.6, showmeans=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(ECOSYSTEM_COLORS[ecosystem])
            pc.set_alpha(0.5)
        bp = ax_d.boxplot(subset, positions=[i], widths=0.18, patch_artist=True,
                         boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.7),
                         whiskerprops=dict(color='gray'),
                         capprops=dict(color='gray'),
                         medianprops=dict(color='black', linewidth=1.5))
ax_d.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
ax_d.set_xticks(range(len(ECOSYSTEM_CLASSES)))
ax_d.set_xticklabels(ECOSYSTEM_CLASSES, rotation=30, ha='right', fontsize=8)
ax_d.set_ylabel('Resistance')
ax_d.set_title('D: Resistance', fontweight='bold')
ax_d.grid(True, alpha=0.3)

# Panel E: Recovery
ax_e = fig.add_subplot(gs[1, 2])
for i, ecosystem in enumerate(ECOSYSTEM_CLASSES):
    subset = recovery_site[recovery_site['water_class'] == ecosystem]['mean_recovery'].dropna().values
    if len(subset) > 0:
        parts = ax_e.violinplot(subset, positions=[i], widths=0.6, showmeans=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(ECOSYSTEM_COLORS[ecosystem])
            pc.set_alpha(0.5)
        bp = ax_e.boxplot(subset, positions=[i], widths=0.18, patch_artist=True,
                         boxprops=dict(facecolor=ECOSYSTEM_COLORS[ecosystem], alpha=0.7),
                         whiskerprops=dict(color='gray'),
                         capprops=dict(color='gray'),
                         medianprops=dict(color='black', linewidth=1.5))
ax_e.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
ax_e.set_xticks(range(len(ECOSYSTEM_CLASSES)))
ax_e.set_xticklabels(ECOSYSTEM_CLASSES, rotation=30, ha='right', fontsize=8)
ax_e.set_ylabel('Recovery')
ax_e.set_title('E: Recovery', fontweight='bold')
ax_e.grid(True, alpha=0.3)

fig.suptitle('Q2 WUE_T Performance Metrics', fontsize=16, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_07_composite.png", width=12, height=10)
print("  Saved: Q2_plot_07_composite.png")

# ============================================================================
# 9. SUMMARY TEXT - EXACT MALONE (Q2_WUE_performance_summary.txt)
# ============================================================================

print("\n" + "="*60)
print("STEP 10: Creating Summary Text")
print("="*60)

summary_text = f"""Q2_WUE_performance.py — run summary
Input:  {input_file}
Output: {output_dir}

Single-SPEI drought metric: {SINGLE_SPEI_LABEL}

--- Data ---
Rows in cleaned monthly data: {len(base_df):,}
Sites: {base_df['site_name'].nunique():,}
Years: {base_df['Year'].min()} - {base_df['Year'].max()}

--- 1. Stability ---
Sites with stability estimates: {len(stability)}
Mean stability: {stability['stability'].mean():.4f}
Median stability: {stability['stability'].median():.4f}
By ecosystem class:
{stability.groupby('water_class').agg(n=('site_name', 'size'), mean_stability=('stability', 'mean')).to_string()}

--- 2. Plasticity: slope ---
Site × timescale rows: {len(plasticity_slope)}
By SPEI timescale and ecosystem:
{plasticity_slope.groupby(['SPEI_timescale', 'water_class']).agg(n=('site_name', 'size'), mean_slope=('slope', 'mean')).to_string()}

--- 3. Plasticity: range ---
Sites with range estimates: {len(plasticity_range)}
Mean plasticity range (95th/5th): {plasticity_range['plasticity_p95_p05'].mean():.4f}
Median plasticity range (95th/5th): {plasticity_range['plasticity_p95_p05'].median():.4f}

--- 4. Resistance ---
Sites with resistance estimates: {len(resistance)}
Mean resistance: {resistance['resistance'].mean():.4f}
Median resistance: {resistance['resistance'].median():.4f}
By ecosystem class:
{resistance.groupby('water_class').agg(n=('site_name', 'size'), mean_resistance=('resistance', 'mean')).to_string()}

--- 5. Recovery ---
Drought events: {len(recovery_events)}
Sites with recovery estimates: {len(recovery_site)}
Mean recovery: {recovery_site['mean_recovery'].mean():.4f}
Median recovery: {recovery_site['mean_recovery'].median():.4f}
By ecosystem class:
{recovery_site.groupby('water_class').agg(n=('site_name', 'size'), mean_recovery=('mean_recovery', 'mean')).to_string()}

--- Summary table coverage ---
Sites in summary table: {len(summary_table)}
Sites with all five metrics: {len(summary_table.dropna(subset=['stability', 'plasticity_slope_SPEI3', 'plasticity_range', 'resistance', 'mean_recovery']))}

--- Ecosystem-class statistical tests ---
{ecosystem_tests.to_string() if len(test_results) > 0 else 'No tests performed'}
"""

capture_to_file(summary_text, "Q2_WUE_performance_summary.txt")
print("  Saved: Q2_WUE_performance_summary.txt")

# ============================================================================
# 10. VERIFY OUTPUTS
# ============================================================================

print("\n" + "="*60)
print("STEP 11: Verifying Outputs")
print("="*60)

missing_outputs = [
    f for f in expected_outputs
    if not os.path.exists(os.path.join(output_dir, f))
]

missing_figures = [
    f for f in expected_figures
    if not os.path.exists(os.path.join(figure_dir, f))
]

if missing_outputs:
    print(f"\nWARNING: Missing expected output files:")
    for f in missing_outputs:
        print(f"  - {f}")
else:
    print(f"\nAll {len(expected_outputs)} expected outputs created.")

if missing_figures:
    print(f"\nWARNING: Missing expected figure files:")
    for f in missing_figures:
        print(f"  - {f}")
else:
    print(f"All {len(expected_figures)} expected figures created.")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("Q2 COMPLETE!")
print("="*60)
print(f"\nAll outputs saved to:")
print(f"  Outputs (CSV/summary): {output_dir}")
print(f"  Figures (PNG):         {figure_dir}")
print(f"  Temp files:            {temp_dir}")

print("\nFiles generated:")
print("  CSV Files:")
print("  - Q2_site_stability.csv")
print("  - Q2_site_plasticity_slope.csv")
print("  - Q2_site_plasticity_range.csv")
print("  - Q2_site_resistance.csv")
print("  - Q2_site_recovery_events.csv")
print("  - Q2_site_recovery.csv")
print("  - Q2_site_performance_summary.csv")
print("  - Q2_ecosystem_class_statistical_tests.csv")
print("  - Q2_WUE_performance_summary.txt")

print("\n  Figures (Malone-style):")
print("  - Q2_plot_01_stability.png")
print("  - Q2_plot_02_plasticity_slope.png")
print("  - Q2_plot_03_plasticity_range.png")
print("  - Q2_plot_04_resistance.png")
print("  - Q2_plot_05_recovery.png")
print("  - Q2_plot_06_performance_heatmap.png")
print("  - Q2_plot_07_composite.png")

print("\n" + "="*60)
print("✓ Q2 workflow complete - EXACT MALONE REPLICATION")
print("  - Same output file names as Malone")
print("  - Same figure file names as Malone")
print("  - Separate directories (safe from Malone)")
print("  - Cleaned old outputs before run")
print("  - Verified all expected outputs created")
print("="*60)