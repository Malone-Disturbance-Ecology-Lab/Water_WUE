"""
FIGURE 2: Q3 Sensitivity Results with Significance Stars
- * added to sites with p < 0.05 (statistically significant slope)
- BRIGHT RED stars (size reduced for better fit)
- Negative slopes: star BELOW the dot
- Positive slopes: star ABOVE the dot
- Pacific region changed to LIGHTER color for better distinction
- Stars placed CLOSER to dots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from scipy.stats import kruskal
from pathlib import Path
from sklearn.linear_model import TheilSenRegressor
import warnings
warnings.filterwarnings('ignore')

# ============================================
# SETUP OUTPUT DIRECTORY
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("FIGURE 2: Q3 Sensitivity Results (with Significance Stars)")
print("="*70)

# ============================================
# STEP 1 — LOAD DATA AND COMPUTE P-VALUES
# ============================================

print("\n📂 STEP 1: Loading data and computing p-values...")

# Load the Q3 final dataset (has slopes)
input_file = output_dir / 'Q3_final_dataset.csv'
df = pd.read_csv(input_file)

# Load monthly data to compute p-values for each site
monthly_filtered_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
nn_medians_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"

df_monthly = pd.read_csv(monthly_filtered_path)
nn_data = pd.read_csv(nn_medians_path)

# Get strict triple intersection sites
wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()

# Compute p-values for each site (using linear regression on SPEI-6 vs WUE_tra)
site_pvalues = {}
for site in df['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site].dropna(subset=['SPEI_6', 'WUE_tra'])
    if len(site_data) >= 5:
        X = site_data['SPEI_6'].values
        y = site_data['WUE_tra'].values
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(X, y)
            site_pvalues[site] = p_value
        except:
            site_pvalues[site] = np.nan
    else:
        site_pvalues[site] = np.nan

# Add p-value and significance flag to dataframe
df['p_value'] = df['site_name'].map(site_pvalues)
df['is_significant'] = (df['p_value'] < 0.05) & (df['p_value'].notna())

# Drop missing slopes
df = df.dropna(subset=['slope_theilsen'])
print(f"  Loaded {len(df)} sites with valid slopes")
print(f"  Significant sites (p < 0.05): {df['is_significant'].sum()} / {len(df)} ({df['is_significant'].sum()/len(df)*100:.1f}%)")

# ============================================
# STEP 2 — DEFINE COLORS
# ============================================
coast_colors = {
    'AK_coast': '#F06292',           # Deep Pink/Maroon
    'West_Coast': '#66C2A5',         # Lighter teal/blue-green (was #00ACC1)
    'Gulf_of_America': '#000000',    # Black (Gulf)
    'Southeast_Atlantic': '#795548', # Brown (Southeast)
    'Atlantic_Coast_North': '#1E88E5' # Bright Blue (Atlantic)
}

coast_full_names = {
    'AK_coast': 'Alaska',
    'West_Coast': 'Pacific',
    'Gulf_of_America': 'Gulf',
    'Southeast_Atlantic': 'Southeast Atlantic',
    'Atlantic_Coast_North': 'Atlantic North'
}

# ============================================
# STEP 3 — GLOBAL CAPPING
# ============================================

global_pos_95 = df['slope_theilsen'].quantile(0.95)
global_neg_95 = df['slope_theilsen'].quantile(0.05)

print(f"\n  GLOBAL capping thresholds:")
print(f"    Lower bound (5th percentile): {global_neg_95:.4f}")
print(f"    Upper bound (95th percentile): {global_pos_95:.4f}")

df['slope_capped'] = df['slope_theilsen'].clip(lower=global_neg_95, upper=global_pos_95)
df['is_extreme_low'] = df['slope_theilsen'] < global_neg_95
df['is_extreme_high'] = df['slope_theilsen'] > global_pos_95

print(f"\n  Extreme values:")
print(f"    Below lower bound: {df['is_extreme_low'].sum()} sites")
print(f"    Above upper bound: {df['is_extreme_high'].sum()} sites")

# Sort by capped slope
df_sorted = df.sort_values('slope_capped').reset_index(drop=True)

# ============================================
# STEP 4 — REORDER COASTS FOR PANEL B
# ============================================

coast_medians = df.groupby('coast_region_analysis')['slope_theilsen'].median().sort_values()
coast_order_ascending = coast_medians.index.tolist()

print(f"\n  Coast order based on median slope (ascending):")
for coast in coast_order_ascending:
    median_val = coast_medians[coast]
    print(f"    {coast}: {median_val:.4f}")

coast_short_labels_ascending = {
    'AK_coast': 'AK',
    'West_Coast': 'Pacific',
    'Gulf_of_America': 'Gulf',
    'Southeast_Atlantic': 'SE Atl',
    'Atlantic_Coast_North': 'Atl N'
}

# Calculate y-axis range for star placement
y_range = global_pos_95 - global_neg_95

# ============================================
# STEP 5 — CREATE FIGURE
# ============================================

print("\n🎨 STEP 5: Creating 2-panel figure...")

fig, axes = plt.subplots(1, 2, figsize=(30, 18), dpi=300)

# ============================================
# PANEL A: SITE-LEVEL DISTRIBUTION WITH STARS
# ============================================

ax1 = axes[0]
slopes = df_sorted['slope_capped'].values
n_sites = len(slopes)
n_negative = (slopes < 0).sum()
n_positive = (slopes > 0).sum()

for i, (slope, row) in enumerate(zip(slopes, df_sorted.iterrows())):
    row_data = row[1]
    color = coast_colors.get(row_data['coast_region_analysis'], '#9E9E9E')
    is_sig = row_data['is_significant']
    
    # Vertical line
    ax1.plot([i, i], [0, slope], color='black', linewidth=2.5, alpha=0.4)
    
    # Plot point
    if row_data['is_extreme_low']:
        y_plot = global_neg_95
        marker = 'v'
    elif row_data['is_extreme_high']:
        y_plot = global_pos_95
        marker = '^'
    else:
        y_plot = slope
        marker = 'o'
    
    ax1.scatter(i, y_plot, color=color, s=200, marker=marker,
               zorder=5, edgecolors='none', alpha=0.85)
    
    # Add bright red star (half size, closer to dot)
    if is_sig:
        if slope < 0:
            # Negative slope: star BELOW the dot
            star_y = y_plot - y_range * 0.06  # Closer (0.06 instead of 0.15)
            va = 'top'
        else:
            # Positive slope: star ABOVE the dot
            star_y = y_plot + y_range * 0.06  # Closer (0.06 instead of 0.15)
            va = 'bottom'
        
        ax1.text(i, star_y, '*', fontsize=60, fontweight='bold',  # Half size (60 instead of 120)
                color='red', ha='center', va=va, zorder=10)

# Reference lines
ax1.axhline(y=0, color='#666666', linestyle=':', linewidth=2.5, alpha=0.7)

# Vertical separator
if n_negative > 0 and n_positive > 0:
    separator_pos = n_negative - 0.5
    ax1.axvline(x=separator_pos, color='black', linestyle='--', linewidth=2, alpha=0.6)
    ax1.text(n_negative/2.5, ax1.get_ylim()[1] * 0.96, 'Negative\nsensitivity', 
            ha='center', va='top', fontsize=38, fontweight='bold', color='#666666')
    ax1.text(n_negative + n_positive/2.5, ax1.get_ylim()[1] * 0.96, 'Positive\nsensitivity', 
            ha='center', va='top', fontsize=38, fontweight='bold', color='#666666')

# Y-axis limits with extra space for stars
y_margin = y_range * 0.05
ax1.set_ylim(global_neg_95 - y_margin - y_range * 0.08, 
             global_pos_95 + y_margin + y_range * 0.08)

# Axes formatting
ax1.set_xlabel('Sites ordered by WUE$_T$ sensitivity', fontsize=48, fontweight='bold', labelpad=15)
ax1.set_ylabel('Sensitivity slope (ΔWUE$_T$ / ΔSPEI-6)', fontsize=48, fontweight='bold', labelpad=15)
ax1.set_xticks([])
ax1.tick_params(axis='both', labelsize=34, width=2, length=10)
ax1.grid(True, alpha=0.3, axis='y', linestyle='--')
ax1.set_axisbelow(True)

for spine in ax1.spines.values():
    spine.set_linewidth(3)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=54, fontweight='bold', va='bottom')

# ============================================
# PANEL A STATISTICS
# ============================================

n_sites_total = len(df)
n_pos = (df['slope_theilsen'] > 0).sum()
n_neg = (df['slope_theilsen'] < 0).sum()
median_slope = df['slope_theilsen'].median()
iqr_slope = df['slope_theilsen'].quantile(0.75) - df['slope_theilsen'].quantile(0.25)
wilcoxon_stat, wilcoxon_p = stats.wilcoxon(df['slope_theilsen'])

def format_p(p_val, sig=2):
    if p_val < 0.001:
        return '<0.001'
    else:
        return f'{p_val:.{sig}g}'

def format_num(x, sig=2):
    if abs(x) < 0.001:
        return '≈0'
    else:
        return f'{x:.{sig}g}'

n_sig_total = df['is_significant'].sum()

stats_text_a = (f"N = {n_sites_total}\n"
                f"Positive: {n_pos} ({n_pos/n_sites_total*100:.0f}%)\n"
                f"Negative: {n_neg} ({n_neg/n_sites_total*100:.0f}%)\n"
                f"Significant: {n_sig_total} ({n_sig_total/n_sites_total*100:.0f}%)\n"
                f"Median = {format_num(median_slope)}\n"
                f"IQR = {format_num(iqr_slope)}\n"
                f"p = {format_p(wilcoxon_p)}")

ax1.text(0.98, 0.05, stats_text_a, transform=ax1.transAxes, fontsize=42,
         verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# ============================================
# PANEL B: COAST EFFECT WITH STARS
# ============================================

ax2 = axes[1]

coast_positions = range(len(coast_order_ascending))

for i, coast in enumerate(coast_order_ascending):
    coast_df = df[df['coast_region_analysis'] == coast]
    coast_slopes_true = coast_df['slope_theilsen'].values
    coast_is_sig = coast_df['is_significant'].values
    coast_extreme_low = coast_df['is_extreme_low'].values
    coast_extreme_high = coast_df['is_extreme_high'].values
    
    if len(coast_slopes_true) > 0:
        median_val = np.median(coast_slopes_true)
        jitter = np.random.normal(0, 0.12, size=len(coast_slopes_true))
        x_positions = np.ones(len(coast_slopes_true)) * i + jitter
        
        for j, (x_pos, slope_true, is_sig, is_low, is_high) in enumerate(zip(
                x_positions, coast_slopes_true, coast_is_sig, 
                coast_extreme_low, coast_extreme_high)):
            
            # Determine y position and marker
            if is_low:
                y_plot = global_neg_95
                marker = 'v'
            elif is_high:
                y_plot = global_pos_95
                marker = '^'
            else:
                y_plot = slope_true
                marker = 'o'
            
            ax2.scatter(x_pos, y_plot, color=coast_colors[coast], s=220, 
                       marker=marker, alpha=0.7, edgecolors='none', zorder=3)
            
            # Add bright red star (half size, closer to dot)
            if is_sig:
                if slope_true < 0:
                    # Negative slope: star BELOW the dot
                    star_y = y_plot - y_range * 0.05  # Closer
                    va = 'top'
                else:
                    # Positive slope: star ABOVE the dot
                    star_y = y_plot + y_range * 0.05  # Closer
                    va = 'bottom'
                
                ax2.text(x_pos, star_y, '*', fontsize=50, fontweight='bold',  # Half size (50 instead of 100)
                        color='red', ha='center', va=va, zorder=10)
        
        ax2.hlines(y=median_val, xmin=i - 0.4, xmax=i + 0.4, 
                  color=coast_colors[coast], linewidth=7, zorder=4, alpha=0.9)

ax2.axhline(y=0, color='#666666', linestyle=':', linewidth=2.5, alpha=0.7)
ax2.set_ylim(global_neg_95 - y_margin - y_range * 0.08, 
             global_pos_95 + y_margin + y_range * 0.08)

ax2.set_xlabel('Coastal region (ordered by median sensitivity)', fontsize=48, fontweight='bold', labelpad=40, x=0.4)
ax2.set_ylabel('Sensitivity slope (ΔWUE$_T$ / ΔSPEI-6)', fontsize=48, fontweight='bold', labelpad=15)

ax2.set_xticks(coast_positions)
ax2.set_xticklabels([coast_short_labels_ascending[coast] for coast in coast_order_ascending], 
                    fontsize=44, fontweight='bold')

for i, coast in enumerate(coast_order_ascending):
    ax2.get_xticklabels()[i].set_color(coast_colors[coast])

ax2.tick_params(axis='both', labelsize=34, width=2, length=10)
ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
ax2.set_axisbelow(True)

for spine in ax2.spines.values():
    spine.set_linewidth(3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, fontsize=54, fontweight='bold', va='bottom')

# ============================================
# PANEL B STATISTICS
# ============================================

coast_groups = [df[df['coast_region_analysis'] == coast]['slope_theilsen'].values 
                for coast in coast_order_ascending if len(df[df['coast_region_analysis'] == coast]) > 0]
kw_stat, kw_p = kruskal(*coast_groups)

if kw_p < 0.0001:
    p_text = 'p < 0.0001'
else:
    p_text = f'p = {kw_p:.4f}'

stats_text_b = f"Kruskal-Wallis {p_text}"

ax2.text(0.5, 0.98, stats_text_b, transform=ax2.transAxes, fontsize=44,
         verticalalignment='top', horizontalalignment='center',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# ============================================
# LEGEND
# ============================================

legend_elements = []
for coast in coast_order_ascending:
    if coast in df['coast_region_analysis'].values:
        count = len(df[df['coast_region_analysis'] == coast])
        n_sig_coast = df[(df['coast_region_analysis'] == coast) & (df['is_significant'] == True)].shape[0]
        label = f"{coast_full_names[coast]} (N={count})"
        if n_sig_coast > 0:
            label += f"  {n_sig_coast}*"
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', 
                   label=label,
                   markerfacecolor=coast_colors[coast], markersize=24,
                   markeredgecolor='none')
        )

legend = fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.06),
                    ncol=3, fontsize=38, title='Coastal Regions', title_fontsize=42, 
                    framealpha=0.95, edgecolor='gray', handlelength=1.2, handleheight=1, 
                    handletextpad=0.05, columnspacing=0.01)

for i, coast in enumerate(coast_order_ascending):
    if coast in df['coast_region_analysis'].values:
        legend_text = legend.get_texts()[i]
        legend_text.set_color(coast_colors[coast])
        legend_text.set_fontsize(38)
        legend_text.set_fontweight('bold')

# ============================================
# ADJUST LAYOUT AND SAVE
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, wspace=0.3)

# Save
png_file = output_dir / 'Figure2_Q3_sensitivity.png'
pdf_file = output_dir / 'Figure2_Q3_sensitivity.pdf'

plt.savefig(png_file, dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n  ✓ Saved: {png_file} (600 DPI)")
print(f"  ✓ Saved: {pdf_file}")

plt.show()

# ============================================
# CONSOLE SUMMARY WITH SIGNIFICANCE DETAILS
# ============================================

print("\n" + "="*70)
print("FIGURE 2: CONSOLE SUMMARY")
print("="*70)

print("\n📊 PANEL A: Site-level sensitivity distribution")
print(f"  • Total sites: {n_sites_total}")
print(f"  • Positive sensitivity: {n_pos} ({n_pos/n_sites_total*100:.1f}%)")
print(f"  • Negative sensitivity: {n_neg} ({n_neg/n_sites_total*100:.1f}%)")
print(f"  • Significant (p<0.05): {n_sig_total} ({n_sig_total/n_sites_total*100:.1f}%)")
print(f"  • Median slope: {format_num(median_slope)}")
print(f"  • IQR: {format_num(iqr_slope)}")
print(f"  • Wilcoxon p-value: {format_p(wilcoxon_p)}")

print("\n📊 PANEL B: Coast region effect (ordered by median slope)")
print(f"  • Kruskal-Wallis p-value: {p_text}")

print("\n  Significant sites by coast region:")
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    n_total_coast = len(coast_df)
    n_sig_coast = coast_df['is_significant'].sum()
    pct_sig = (n_sig_coast / n_total_coast * 100) if n_total_coast > 0 else 0
    median_val = coast_medians[coast]
    print(f"    {coast_full_names[coast]}: n={n_total_coast}, significant={n_sig_coast} ({pct_sig:.0f}%), median={format_num(median_val)}")

print("\n✅ UPDATES APPLIED:")
print("  • * = significant slope (p < 0.05) - BRIGHT RED")
print("  • Star size reduced by half (60pt Panel A, 50pt Panel B)")
print("  • Stars placed CLOSER to dots (offset 0.06 instead of 0.15)")
print("  • Pacific region changed to LIGHTER color (#66C2A5)")
print("  • Atlantic Coast North remains BRIGHT BLUE (#1E88E5)")
print("  • Negative slopes: star BELOW the dot")
print("  • Positive slopes: star ABOVE the dot")

print("\n" + "="*70)
print("✅ FIGURE 2 COMPLETE")
print(f"📁 Files saved to: {output_dir}")
print("="*70)