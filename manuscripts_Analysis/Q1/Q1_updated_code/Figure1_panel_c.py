import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde, kruskal

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
fig_output = os.path.join(figure_dir, "Q1_panels_g_h_TET_coast.png")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
# Monthly data for Panel g (dots) and Panel h (densities)
monthly_file = os.path.join(base_dir, "outputs", "T_ET_ratio_coast", "Q1_NN_TET_monthly_by_coast_waterclass.csv")
df_monthly = pd.read_csv(monthly_file)

# Site‑level data for coast means and p‑value
site_file = os.path.join(base_dir, "outputs", "T_ET_ratio_coast", "Q1_NN_TET_site_level_by_coast_waterclass.csv")
df_site = pd.read_csv(site_file)

# ------------------------------------------------------------------------------
# Determine coast order from site‑level means (low → high)
# ------------------------------------------------------------------------------
coast_means = df_site.groupby('coast_region', observed=True)['mean_TET'].mean().sort_values(ascending=True)
coast_order = coast_means.index.tolist()

print("\nCoast order (low → high mean T:ET):")
for coast, mean_val in coast_means.items():
    print(f"  {coast}: {mean_val:.4f}")

# Convert to categorical for consistent ordering
df_site['coast_region'] = pd.Categorical(df_site['coast_region'], categories=coast_order, ordered=True)
df_monthly['coast_region'] = pd.Categorical(df_monthly['coast_region'], categories=coast_order, ordered=True)

# ------------------------------------------------------------------------------
# Custom colours (Atlantic=green, Gulf=cyan, Pacific=red, AK=brown)
# ------------------------------------------------------------------------------
custom_colors = {
    'Atlantic Coast': '#2E8B57',   # SeaGreen
    'Gulf Coast':     '#00CED1',   # DarkTurquoise
    'Pacific Coast':  '#DC143C',   # Crimson
    'AK Coast':       '#8B4513'    # SaddleBrown
}

# Short names for legend (remove "Coast")
short_names = {
    'Atlantic Coast': 'Atlantic',
    'Pacific Coast': 'Pacific',
    'Gulf Coast': 'Gulf',
    'AK Coast': 'AK'
}

# ------------------------------------------------------------------------------
# Compute Kruskal‑Wallis p‑value on site‑level means (for the p‑value box)
# ------------------------------------------------------------------------------
groups = [df_site[df_site['coast_region'] == coast]['mean_TET'].dropna().values for coast in coast_order]
h_stat, p_val = kruskal(*groups)

# Format p‑value text (only the p‑value itself)
if p_val < 0.001:
    p_text = "p < 0.001"
else:
    p_text = f"p = {p_val:.4f}"

# ------------------------------------------------------------------------------
# Plotting – one row, two panels
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(24, 6))

plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 28,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
})

# Format axes
for ax in axes:
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# ------------------------------------------------------------
# Panel g – Jittered points (using monthly data), median = black, mean label above
# ------------------------------------------------------------
ax_g = axes[0]
positions = np.arange(1, len(coast_order) + 1)

for i, coast in enumerate(coast_order):
    # Use monthly data for the dots (many points)
    vals = df_monthly[df_monthly['coast_region'] == coast]['TET_ratio'].dropna().values
    x_jitter = np.random.normal(i + 1, 0.06, size=len(vals))
    ax_g.scatter(x_jitter, vals, alpha=0.4, s=25,
                 color=custom_colors[coast], label='_nolegend_',
                 edgecolors='black', linewidth=0.5)
    # Median line – black
    median_val = np.median(vals)
    ax_g.plot([i + 1 - 0.15, i + 1 + 0.15], [median_val, median_val],
              color='black', linewidth=2.5, solid_capstyle='butt')
    # Mean label above group – site‑level mean (from coast_means)
    mean_val = coast_means[coast]
    ax_g.text(i + 1, 0.92, f"{mean_val:.2f}",
              ha='center', va='bottom', fontsize=24,
              color='black', fontweight='bold')

ax_g.axhline(0.5, color='gray', linestyle='--', alpha=0.4, linewidth=2)
ax_g.set_xticks(positions)
ax_g.set_xticklabels(coast_order)
for tick, coast in zip(ax_g.get_xticklabels(), coast_order):
    tick.set_color(custom_colors[coast])

# Y‑axis label: simply "T:ET ratio"
ax_g.set_ylabel('T:ET ratio', fontsize=28, labelpad=10)
ax_g.set_xlabel('')
ax_g.set_ylim(0, 1)
ax_g.set_title('')
ax_g.text(-0.15, 1.02, 'g)', transform=ax_g.transAxes,
          fontsize=34, fontweight='bold', va='bottom', ha='left')

# Add p‑value text box (only the p‑value)
ax_g.text(0.97, 0.05, p_text, transform=ax_g.transAxes,
          fontsize=24, ha='right', va='bottom',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

# ------------------------------------------------------------
# Panel h – Density curves (monthly data, same coast order)
# ------------------------------------------------------------
ax_h = axes[1]
x_density = np.linspace(0, 1, 300)

for coast in coast_order:
    data = df_monthly[df_monthly['coast_region'] == coast]['TET_ratio'].dropna()
    if len(data) > 2:
        kde = gaussian_kde(data)
        y = kde(x_density)
        ax_h.plot(x_density, y, label=short_names[coast],
                  color=custom_colors[coast], linewidth=2.5)
        ax_h.fill_between(x_density, 0, y, color=custom_colors[coast], alpha=0.15)

ax_h.set_xlabel('T:ET ratio (monthly)', fontsize=28)
ax_h.set_ylabel('Probability density', fontsize=28, labelpad=10)
ax_h.set_ylim(0, None)
ax_h.set_xlim(0, 1)
ax_h.set_title('')

# Legend – top‑left, two columns, coloured text, reduced column spacing
legend = ax_h.legend(loc='upper left', ncol=2, fontsize=24,
                     framealpha=0.9, columnspacing=0.5)
for text, coast in zip(legend.get_texts(), coast_order):
    text.set_color(custom_colors[coast])

ax_h.text(-0.15, 1.02, 'h)', transform=ax_h.transAxes,
          fontsize=34, fontweight='bold', va='bottom', ha='left')

# ------------------------------------------------------------------------------
# Adjust spacing and save
# ------------------------------------------------------------------------------
plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.15, wspace=0.3)
fig.savefig(fig_output, dpi=600, bbox_inches='tight', pad_inches=0.5, facecolor='white')
print(f"Panels g and h saved to: {fig_output}")