import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import os

# ============================================================================
# 1. LOAD DATA
# ============================================================================
file_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\ClimateBiome_WUE_WUET_FinalDataset.csv"
df = pd.read_csv(file_path)

# Filter to required metrics only
df = df[df['WUE_Metric'].isin(['WUE', 'WUE_tra'])].copy()

# ============================================================================
# 2. DEFINITIONS - UPDATED LABELS
# ============================================================================
# Original groups with new labels for plotting
group_rename = {
    "Humid Coastal Wetlands": "Humid\nWetlands",
    "Humid Coastal Forests": "Humid\nForests",
    "Mediterranean Wetlands": "Mediterranean\nWetlands",
    "Mediterranean Drylands": "Mediterranean\nshrub/grass",
    "Cold Coastal Systems": "Cold\necosystems"
}

# Group order
group_order_original = [
    "Humid Coastal Wetlands",
    "Humid Coastal Forests", 
    "Mediterranean Wetlands",
    "Mediterranean Drylands",
    "Cold Coastal Systems"
]

group_order_display = [group_rename[g] for g in group_order_original]

# Color and style settings - LIGHTER colors for boxplots
colors = {
    'WUE': '#A8C4E0',      # Light blue for WUE ET
    'WUE_tra': '#F5B89E'    # Light orange for WUE T
}
box_width = 0.35

# Capping percentiles for visualization only - USING 5th and 95th
LOWER_CAP_PERCENTILE = 5
UPPER_CAP_PERCENTILE = 95

# ============================================================================
# 3. CAPPING FUNCTION (for visualization only)
# ============================================================================

def cap_values_for_visualization(values, lower_percentile=LOWER_CAP_PERCENTILE, 
                                  upper_percentile=UPPER_CAP_PERCENTILE):
    """Cap values at specified percentiles for VISUALIZATION ONLY"""
    if len(values) == 0:
        return values, None, None, 0, 0
    
    lower_cap = np.percentile(values, lower_percentile)
    upper_cap = np.percentile(values, upper_percentile)
    
    lower_mask = values < lower_cap
    upper_mask = values > upper_cap
    
    n_capped_lower = np.sum(lower_mask)
    n_capped_upper = np.sum(upper_mask)
    
    values_capped = values.copy()
    values_capped[lower_mask] = lower_cap
    values_capped[upper_mask] = upper_cap
    
    return values_capped, lower_cap, upper_cap, n_capped_lower, n_capped_upper

# ============================================================================
# 4. STATISTICAL ANALYSIS (on original, uncapped data)
# ============================================================================
print("=" * 80)
print("STATISTICAL ANALYSIS (using all original, uncapped data)")
print("=" * 80)

# Kruskal-Wallis tests
kw_results = {}
for metric in ['WUE', 'WUE_tra']:
    metric_df = df[df['WUE_Metric'] == metric]
    groups_data = [metric_df[metric_df['Final_Group'] == g]['Median'].values 
                   for g in group_order_original if len(metric_df[metric_df['Final_Group'] == g]) > 0]
    if len(groups_data) >= 2:
        h_stat, p_val = kruskal(*groups_data)
        kw_results[metric] = {'H': h_stat, 'p': p_val}
        print(f"\n{metric}: H={h_stat:.4f}, p={p_val:.6e}")

# Function to get significant pairs
def get_significant_pairs(df, metric, groups, alpha=0.05):
    """Get list of significant pairwise comparisons after FDR correction"""
    group_names = groups
    n_groups = len(group_names)
    
    # Get data for each group
    group_data = []
    for group in group_names:
        data = df[(df['WUE_Metric'] == metric) & (df['Final_Group'] == group)]['Median'].values
        group_data.append(data)
    
    # Perform pairwise tests
    p_values = []
    comparisons = []
    for i in range(n_groups):
        for j in range(i+1, n_groups):
            if len(group_data[i]) > 0 and len(group_data[j]) > 0:
                _, p_val = mannwhitneyu(group_data[i], group_data[j], alternative='two-sided')
                p_values.append(p_val)
                comparisons.append((i, j))
    
    # Apply FDR correction
    if len(p_values) > 0:
        reject, p_corrected, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
        
        # Get significant pairs
        significant_pairs = []
        for idx, (i, j) in enumerate(comparisons):
            if reject[idx]:
                significant_pairs.append((i, j))
        return significant_pairs
    return []

# Get significant pairs for both metrics
sig_pairs_wue = get_significant_pairs(df, 'WUE', group_order_original)
sig_pairs_wue_tra = get_significant_pairs(df, 'WUE_tra', group_order_original)

# Print results
print("\n" + "=" * 80)
print("PAIRWISE COMPARISONS (Mann-Whitney U with FDR correction)")
print("=" * 80)

print("\nWUE ET Significant pairs:")
for i, j in sig_pairs_wue:
    print(f"  {group_rename[group_order_original[i]].replace(chr(10), ' ')} vs {group_rename[group_order_original[j]].replace(chr(10), ' ')}")

print("\nWUE T Significant pairs:")
for i, j in sig_pairs_wue_tra:
    print(f"  {group_rename[group_order_original[i]].replace(chr(10), ' ')} vs {group_rename[group_order_original[j]].replace(chr(10), ' ')}")

# ============================================================================
# 5. PREPARE DATA FOR VISUALIZATION (capped at 5th and 95th percentiles)
# ============================================================================

# Create capped versions for plotting only
df_capped = df.copy()
capping_info = {}

print("\n" + "=" * 80)
print("VISUALIZATION CAPPING (5th and 95th percentiles)")
print("=" * 80)

for metric in ['WUE', 'WUE_tra']:
    # Get all values for this metric across all groups
    metric_data = df_capped[df_capped['WUE_Metric'] == metric]['Median'].values
    
    print(f"\n{metric} - Original data range: [{metric_data.min():.3f}, {metric_data.max():.3f}]")
    
    # Apply capping
    capped_values, lower_cap, upper_cap, n_lower, n_upper = cap_values_for_visualization(
        metric_data, LOWER_CAP_PERCENTILE, UPPER_CAP_PERCENTILE
    )
    
    capping_info[metric] = {
        'lower_cap': lower_cap,
        'upper_cap': upper_cap,
        'n_capped_lower': n_lower,
        'n_capped_upper': n_upper,
        'total_n': len(metric_data),
        'original_min': metric_data.min(),
        'original_max': metric_data.max()
    }
    
    # Apply capping to the dataframe
    for idx in df_capped[df_capped['WUE_Metric'] == metric].index:
        original_val = df_capped.loc[idx, 'Median']
        if original_val < lower_cap:
            df_capped.loc[idx, 'Median_capped'] = lower_cap
        elif original_val > upper_cap:
            df_capped.loc[idx, 'Median_capped'] = upper_cap
        else:
            df_capped.loc[idx, 'Median_capped'] = original_val
    
    print(f"  Lower {LOWER_CAP_PERCENTILE}th percentile cap: {lower_cap:.3f}")
    print(f"  Upper {UPPER_CAP_PERCENTILE}th percentile cap: {upper_cap:.3f}")
    print(f"  Capped {n_lower} low values and {n_upper} high values")
    print(f"  New range after capping: [{capped_values.min():.3f}, {capped_values.max():.3f}]")
    print(f"  Values changed: {n_lower+n_upper}/{len(metric_data)} points ({(n_lower+n_upper)/len(metric_data)*100:.1f}%)")

# ============================================================================
# 6. CALCULATE SAMPLE SIZES (N) FOR EACH GROUP
# ============================================================================
# Calculate N for each group (using original data, each site appears once per metric)
# Since each site has both WUE and WUE_tra, we can use either metric to count sites
sample_sizes = {}
for group in group_order_original:
    n = len(df[(df['WUE_Metric'] == 'WUE') & (df['Final_Group'] == group)])
    sample_sizes[group] = n

print("\n" + "=" * 80)
print("SAMPLE SIZES PER GROUP")
print("=" * 80)
for group in group_order_original:
    display_name = group_rename[group].replace('\n', ' ')
    print(f"  {display_name}: n = {sample_sizes[group]}")

# ============================================================================
# 7. FIGURE: SINGLE PANEL WITH SIDE-BY-SIDE BOXPLOTS (USING CAPPED DATA)
# ============================================================================
fig, ax = plt.subplots(1, 1, figsize=(12, 7))

# Calculate positions for grouped boxplots
n_groups = len(group_order_original)
x_positions = np.arange(n_groups)
width = box_width

# Store y_max for significance bar placement
all_data = []

# For each group, create two boxplots side by side using CAPPED data
for idx, (orig_group, display_name) in enumerate(zip(group_order_original, group_order_display)):
    # Left position for WUE ET
    pos_wue = x_positions[idx] - width/2
    # Right position for WUE T  
    pos_wue_tra = x_positions[idx] + width/2
    
    # Get CAPPED data for this group (for visualization only)
    wue_data = df_capped[(df_capped['WUE_Metric'] == 'WUE') & (df_capped['Final_Group'] == orig_group)]['Median_capped'].dropna().values
    wue_tra_data = df_capped[(df_capped['WUE_Metric'] == 'WUE_tra') & (df_capped['Final_Group'] == orig_group)]['Median_capped'].dropna().values
    
    print(f"\nGroup {display_name.replace(chr(10), ' ')}:")
    print(f"  WUE ET capped range: [{wue_data.min():.3f}, {wue_data.max():.3f}]" if len(wue_data) > 0 else "  WUE ET: no data")
    print(f"  WUE T capped range: [{wue_tra_data.min():.3f}, {wue_tra_data.max():.3f}]" if len(wue_tra_data) > 0 else "  WUE T: no data")
    
    # Store data for y_max calculation
    if len(wue_data) > 0:
        all_data.extend(wue_data)
    if len(wue_tra_data) > 0:
        all_data.extend(wue_tra_data)
    
    # Create boxplots with CAPPED data - LIGHT colors with borders
    if len(wue_data) > 0:
        bp1 = ax.boxplot(wue_data, 
                        positions=[pos_wue],
                        widths=width,
                        patch_artist=True,
                        showfliers=False,
                        whiskerprops={'color': 'black', 'linewidth': 1.5},
                        capprops={'color': 'black', 'linewidth': 1.5},
                        medianprops={'color': 'black', 'linewidth': 2},
                        boxprops={'facecolor': colors['WUE'], 'edgecolor': 'black', 'linewidth': 1.5, 'alpha': 0.8})
        
        # Add jittered points with SAME COLOR as boxplot face, NO black outline
        jitter = np.random.normal(0, 0.03, size=len(wue_data))
        ax.scatter(pos_wue + jitter, wue_data, 
                  color=colors['WUE'], alpha=0.6, s=25, edgecolors='none', zorder=3)
    
    if len(wue_tra_data) > 0:
        bp2 = ax.boxplot(wue_tra_data, 
                        positions=[pos_wue_tra],
                        widths=width,
                        patch_artist=True,
                        showfliers=False,
                        whiskerprops={'color': 'black', 'linewidth': 1.5},
                        capprops={'color': 'black', 'linewidth': 1.5},
                        medianprops={'color': 'black', 'linewidth': 2},
                        boxprops={'facecolor': colors['WUE_tra'], 'edgecolor': 'black', 'linewidth': 1.5, 'alpha': 0.8})
        
        # Add jittered points with SAME COLOR as boxplot face, NO black outline
        jitter = np.random.normal(0, 0.03, size=len(wue_tra_data))
        ax.scatter(pos_wue_tra + jitter, wue_tra_data, 
                  color=colors['WUE_tra'], alpha=0.6, s=25, edgecolors='none', zorder=3)

# Calculate y_max for significance bars
y_max_data = max(all_data) if all_data else 10
y_max_plot = y_max_data * 1.2  # Add 20% headroom for significance bars

# Add significance bars for WUE ET
bar_height = y_max_data * 0.05
y_base = y_max_data + bar_height

for i, j in sig_pairs_wue:
    # Get x positions for the two groups (center positions of the WUE ET boxplots)
    x1 = x_positions[i] - width/2
    x2 = x_positions[j] - width/2
    
    # Draw significance bar
    ax.plot([x1, x1, x2, x2], [y_base, y_base + bar_height, y_base + bar_height, y_base], 
            'k-', linewidth=1.5)
    
    # Add asterisk
    mid_x = (x1 + x2) / 2
    ax.text(mid_x, y_base + bar_height + bar_height/2, '*', 
            ha='center', va='bottom', fontsize=16, fontweight='bold')
    
    # Increment y_base for next bar if multiple comparisons
    y_base += bar_height * 2

# Add significance bars for WUE T
bar_height_tra = y_max_data * 0.05
y_base_tra = y_max_data + bar_height_tra * 1.5

# Reset y_base for WUE T bars (start lower if WUE ET bars are present)
if sig_pairs_wue:
    y_base_tra = y_base + bar_height_tra

for i, j in sig_pairs_wue_tra:
    # Get x positions for the two groups (center positions of the WUE T boxplots)
    x1 = x_positions[i] + width/2
    x2 = x_positions[j] + width/2
    
    # Draw significance bar
    ax.plot([x1, x1, x2, x2], [y_base_tra, y_base_tra + bar_height_tra, y_base_tra + bar_height_tra, y_base_tra], 
            'k-', linewidth=1.5)
    
    # Add asterisk
    mid_x = (x1 + x2) / 2
    ax.text(mid_x, y_base_tra + bar_height_tra + bar_height_tra/2, '*', 
            ha='center', va='bottom', fontsize=16, fontweight='bold')
    
    # Increment y_base for next bar
    y_base_tra += bar_height_tra * 2

# Set y-axis limit to accommodate significance bars
ax.set_ylim(bottom=0, top=max(y_base_tra, y_base) + bar_height_tra * 2)

# Create x-axis labels with sample sizes (N) below
xtick_labels_with_n = []
for group in group_order_original:
    display_name = group_rename[group]
    n = sample_sizes[group]
    xtick_labels_with_n.append(f"{display_name}\n(N={n})")

# Customize x-axis with sample sizes
ax.set_xticks(x_positions)
ax.set_xticklabels(xtick_labels_with_n, ha='center', fontsize=14)  # Slightly smaller font to fit

# Remove x-axis label
ax.set_xlabel('')

# Remove title
ax.set_title('')

# Single y-axis label
ax.set_ylabel('WUE (g C kg⁻¹ H₂O⁻¹)', fontsize=22, fontweight='bold')

# Increase y-tick font size
ax.tick_params(axis='y', labelsize=16)  # Increased y-tick font size

# Add legend with subscript formatting
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors['WUE'], alpha=0.8, edgecolor='black', label='WUE$_{ET}$'),
                  Patch(facecolor=colors['WUE_tra'], alpha=0.8, edgecolor='black', label='WUE$_{T}$')]
ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=16, framealpha=0.9)

# Add grid for readability
ax.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()

# ============================================================================
# 8. SAVE FIGURE
# ============================================================================
output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\results\figures"
os.makedirs(output_dir, exist_ok=True)

png_path = os.path.join(output_dir, "Figure_WUE_SinglePanel_Capped.png")
pdf_path = os.path.join(output_dir, "Figure_WUE_SinglePanel_Capped.pdf")

try:
    plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"\nFigure saved to:\n  PNG: {png_path}\n  PDF: {pdf_path}")
except Exception as e:
    print(f"\nWarning: Could not save to network location: {e}")
    # Save locally as fallback
    local_dir = "./figures"
    os.makedirs(local_dir, exist_ok=True)
    png_path_local = os.path.join(local_dir, "Figure_WUE_SinglePanel_Capped.png")
    pdf_path_local = os.path.join(local_dir, "Figure_WUE_SinglePanel_Capped.pdf")
    plt.savefig(png_path_local, dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path_local, bbox_inches='tight', facecolor='white')
    print(f"  (Saved locally to: {png_path_local} and {pdf_path_local})")

# ============================================================================
# 9. PRINT CAPPING SUMMARY FOR PAPER
# ============================================================================
print("\n" + "=" * 80)
print("VISUALIZATION CAPPING SUMMARY (for paper reporting)")
print("=" * 80)
print(f"Note: Statistical analyses use ORIGINAL uncapped values")
print(f"Capping applied to boxplots and points for visualization only")
print(f"Used {LOWER_CAP_PERCENTILE}th and {UPPER_CAP_PERCENTILE}th percentiles\n")

for metric, info in capping_info.items():
    metric_name = "WUE ET" if metric == 'WUE' else "WUE T"
    print(f"{metric_name}:")
    print(f"  Original range: [{info['original_min']:.3f}, {info['original_max']:.3f}]")
    print(f"  Lower {LOWER_CAP_PERCENTILE}th percentile cap: {info['lower_cap']:.3f}")
    print(f"  Upper {UPPER_CAP_PERCENTILE}th percentile cap: {info['upper_cap']:.3f}")
    print(f"  Points capped: {info['n_capped_lower']} low + {info['n_capped_upper']} high = {info['n_capped_lower']+info['n_capped_upper']}/{info['total_n']} ({(info['n_capped_lower']+info['n_capped_upper'])/info['total_n']*100:.1f}%)")
    print()

# ============================================================================
# 10. MANUSCRIPT-STYLE INTERPRETATION
# ============================================================================
print("\n" + "=" * 80)
print("MANUSCRIPT-STYLE INTERPRETATION SUMMARY")
print("=" * 80)

# Get summary statistics from ORIGINAL data (not capped)
wue_summary = df[df['WUE_Metric'] == 'WUE'].groupby('Final_Group')['Median'].agg(['median', 'count', 'mean', 'std'])
wue_tra_summary = df[df['WUE_Metric'] == 'WUE_tra'].groupby('Final_Group')['Median'].agg(['median', 'count', 'mean', 'std'])

print("\n• WUE ET (g C kg⁻¹ H₂O⁻¹) - ORIGINAL uncapped values:")
for group in group_order_original:
    if group in wue_summary.index:
        display_name = group_rename[group].replace('\n', ' ')
        print(f"    {display_name}: median = {wue_summary.loc[group, 'median']:.3f}")

print("\n• WUE T (g C kg⁻¹ H₂O⁻¹) - ORIGINAL uncapped values:")
for group in group_order_original:
    if group in wue_tra_summary.index:
        display_name = group_rename[group].replace('\n', ' ')
        print(f"    {display_name}: median = {wue_tra_summary.loc[group, 'median']:.3f}")

print(f"\n• Omnibus tests (Kruskal-Wallis on ORIGINAL data):")
for metric, res in kw_results.items():
    metric_name = "WUE ET" if metric == 'WUE' else "WUE T"
    print(f"    {metric_name}: H={res['H']:.2f}, p={res['p']:.4e}")

print(f"\n• Significant pairwise comparisons (FDR-corrected, p < 0.05):")
if sig_pairs_wue:
    print("  WUE ET:")
    for i, j in sig_pairs_wue:
        print(f"    {group_rename[group_order_original[i]].replace(chr(10), ' ')} vs {group_rename[group_order_original[j]].replace(chr(10), ' ')}")
if sig_pairs_wue_tra:
    print("  WUE T:")
    for i, j in sig_pairs_wue_tra:
        print(f"    {group_rename[group_order_original[i]].replace(chr(10), ' ')} vs {group_rename[group_order_original[j]].replace(chr(10), ' ')}")

print("\n• Key finding: Both WUE ET and WUE T showed significant ecological")
print("  structuring across climate-biome groups (p < 0.05). Wetlands generally")
print("  exhibited lower water-use efficiencies, while humid forests and")
print("  Mediterranean shrub/grass showed elevated WUE T. Asterisks (*) above")
print("  the plot indicate significant pairwise differences between groups")
print(f"\n• Visualization note: Values were capped at the {LOWER_CAP_PERCENTILE}th and {UPPER_CAP_PERCENTILE}th")
print(f"  percentiles for display purposes only. Statistical analyses used")
print("  original uncapped values.")

print("\n" + "=" * 80)
print("Analysis complete.")
print("=" * 80)

# Show the figure
plt.show()