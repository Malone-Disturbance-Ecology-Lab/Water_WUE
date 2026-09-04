# -*- coding: utf-8 -*-
"""
Combined Q1 panels: A (a–c), B (d–f), and G/H (g–h)
Updated with:
- mathematical reference line in d–f
- panel labels outside axes (top-left)
- brackets with significance stars for post-hoc comparisons
- increased row spacing and bracket separation
- bbox_inches='tight' to prevent cut-off when saving
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal
import re

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
fig_output = os.path.join(figure_dir, "Q1_combined_panels_A_B_G_H_with_brackets_final.png")

monthly_data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
summary_file = os.path.join(output_dir, "Q1_final_near_normal_summary_by_difference_type.csv")
df_summary = pd.read_csv(summary_file)

smooth_file = os.path.join(output_dir, "Q1_final_smooth_gam_predictions_TET.csv")
gam_comp_file = os.path.join(output_dir, "Q1_final_gam_model_comparison.csv")
edf_file = os.path.join(output_dir, "Q1_final_smooth_gam_smooth_terms.csv")
df_smooth = pd.read_csv(smooth_file)
gam_comp = pd.read_csv(gam_comp_file)
edf_df = pd.read_csv(edf_file)

wue_monthly = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                       "Q1_NN_WUE_T_monthly_by_coast_waterclass.csv"))
wue_site = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                    "Q1_NN_WUE_T_site_level_by_coast_waterclass.csv"))
tet_monthly = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                       "Q1_NN_TET_monthly_by_coast_waterclass.csv"))
tet_site = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                    "Q1_NN_TET_site_level_by_coast_waterclass.csv"))

# ------------------------------------------------------------------------------
# Common settings
# ------------------------------------------------------------------------------
ecosystem_order = ["Upland", "Freshwater", "Saline"]
ecosystem_colours = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500",
}
coast_colors = {
    'Atlantic Coast': '#A50F15',
    'Pacific Coast':  '#0072B2',
    'Gulf Coast':     '#4D4D4D',
    'AK Coast':       '#009E73'
}
coast_short = {
    'Atlantic Coast': 'Atlantic',
    'Pacific Coast': 'Pacific',
    'Gulf Coast': 'Gulf',
    'AK Coast': 'Alaska'
}

diff_label_mapping = {
    "WUE_ET - WUE_T": r"WUE$_{ET}$ – WUE$_{T}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_ET - WUE_E": r"WUE$_{ET}$ – WUE$_{E}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_E - WUE_T": r"WUE$_{E}$ – WUE$_{T}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
}
diff_labels = list(diff_label_mapping.keys())

# ------------------------------------------------------------------------------
# Data preparation for Panel A
# ------------------------------------------------------------------------------
df_summary = df_summary[df_summary["difference_label"].isin(diff_labels)]
df_summary["difference_label"] = pd.Categorical(df_summary["difference_label"],
                                                categories=diff_labels, ordered=True)
df_summary["water_class"] = pd.Categorical(df_summary["water_class"],
                                           categories=ecosystem_order, ordered=True)

# ------------------------------------------------------------------------------
# Data preparation for Panel B (smooth GAM)
# ------------------------------------------------------------------------------
df_smooth = df_smooth[df_smooth["difference_label"].isin(diff_labels)]
df_smooth["difference_label"] = pd.Categorical(df_smooth["difference_label"],
                                               categories=diff_labels, ordered=True)
df_smooth["ecosystem_label"] = pd.Categorical(df_smooth["ecosystem_label"],
                                              categories=ecosystem_order, ordered=True)

# ------------------------------------------------------------------------------
# Data preparation for Panels g/h (capping, ordering)
# ------------------------------------------------------------------------------
def cap_by_group(df, group_col, value_col, percentile=99):
    capped = df.copy()
    capped['WUE_T_capped'] = np.nan
    for group in df[group_col].unique():
        mask = df[group_col] == group
        vals = df.loc[mask, value_col].dropna()
        if len(vals) > 0:
            p = np.percentile(vals, percentile)
            capped.loc[mask, 'WUE_T_capped'] = np.minimum(df.loc[mask, value_col], p)
    return capped

wue_monthly = cap_by_group(wue_monthly, 'water_class', 'WUE_T', 99)

eco_order_gh = ["Upland", "Freshwater", "Saline"]

coast_means_order = tet_site.groupby('coast_region')['mean_TET'].mean().sort_values()
coast_order = coast_means_order.index.tolist()

wue_monthly['water_class'] = pd.Categorical(wue_monthly['water_class'], categories=eco_order_gh, ordered=True)
wue_site['water_class'] = pd.Categorical(wue_site['water_class'], categories=eco_order_gh, ordered=True)
tet_monthly['coast_region'] = pd.Categorical(tet_monthly['coast_region'], categories=coast_order, ordered=True)
tet_site['coast_region'] = pd.Categorical(tet_site['coast_region'], categories=coast_order, ordered=True)

# Omnibus tests
eco_groups = [wue_site[wue_site['water_class'] == e]['mean_WUE_T'].dropna().values for e in eco_order_gh]
h_eco, p_eco = kruskal(*eco_groups)
p_text_eco = f"p = {p_eco:.4f}" if p_eco >= 0.001 else "p < 0.001"

coast_groups = [tet_site[tet_site['coast_region'] == c]['mean_TET'].dropna().values for c in coast_order]
h_coast, p_coast = kruskal(*coast_groups)
p_text_coast = f"p = {p_coast:.3f}" if p_coast >= 0.001 else "p < 0.001"

# Group statistics for plotting
eco_medians = wue_monthly.groupby('water_class', observed=False)['WUE_T'].median().reindex(eco_order_gh)
eco_means = wue_site.groupby('water_class', observed=False)['mean_WUE_T'].mean().reindex(eco_order_gh)

coast_stats = tet_site.groupby('coast_region', observed=False)['mean_TET'].agg(['mean', 'std', 'count']).reindex(coast_order)
coast_means = coast_stats['mean']
coast_sem = coast_stats['std'] / np.sqrt(coast_stats['count'])
coast_ci = 1.96 * coast_sem

# ------------------------------------------------------------------------------
# Compute mathematical reference curves
# ------------------------------------------------------------------------------
print("\nLoading monthly data for reference calculation...")
monthly = pd.read_csv(monthly_data_path)
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].where(
        monthly[col].isna(),
        monthly[col].astype(str).str.strip()
    )
    monthly[col] = monthly[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

ecosystem_classes = ["Upland", "Freshwater", "Saline"]
required_cols = ["site_name", "Year", "month", "water_class", "Trans_ratio",
                 "WUE", "WUE_tra", "WUE_eva", "SPEI_1"]
model_base = monthly.dropna(subset=required_cols).copy()
model_base = model_base[
    (model_base['site_name'] != "") &
    (model_base['water_class'].isin(ecosystem_classes)) &
    np.isfinite(model_base['Trans_ratio']) &
    (model_base['Trans_ratio'] >= 0) & (model_base['Trans_ratio'] <= 1) &
    np.isfinite(model_base['WUE']) &
    np.isfinite(model_base['WUE_tra']) &
    np.isfinite(model_base['WUE_eva']) &
    np.isfinite(model_base['SPEI_1'])
].copy()

nn_data = model_base[(model_base['SPEI_1'] >= -1) & (model_base['SPEI_1'] <= 1)]
n_nn = len(nn_data)
mean_wue_et = nn_data['WUE'].mean()
print(f"N used for mathematical reference = {n_nn}")
print(f"Mean WUE_ET used for mathematical reference = {mean_wue_et:.4f}")

if n_nn != 1202:
    raise ValueError(f"Expected N=1202, got {n_nn}. Check filtering.")

tet_range = (df_smooth['Trans_ratio'].min(), df_smooth['Trans_ratio'].max())
tet_seq = np.linspace(tet_range[0], tet_range[1], 200)
ref_curves = {}
for diff_label in diff_labels:
    r = tet_seq
    if diff_label == "WUE_ET - WUE_T":
        ref = mean_wue_et * (1 - 1/r)
    elif diff_label == "WUE_ET - WUE_E":
        ref = mean_wue_et * (1 - 1/(1-r))
    elif diff_label == "WUE_E - WUE_T":
        ref = mean_wue_et * (1/(1-r) - 1/r)
    else:
        raise ValueError(f"Unknown diff label: {diff_label}")
    ref_curves[diff_label] = (tet_seq, ref)

# ------------------------------------------------------------------------------
# Print removed statistics to console
# ------------------------------------------------------------------------------
print("\n" + "="*60)
print("GAM STATISTICS (removed from figure, now printed to console)")
print("="*60)

smooth_row = gam_comp[gam_comp["model"] == "q1_smooth_gam"].iloc[0]
print(f"\nSmooth GAM adjusted R²: {smooth_row['adjusted_r_squared']:.4f}")
print(f"Smooth GAM deviance explained: {smooth_row['deviance_explained']*100:.2f}%")
print(f"Smooth GAM AIC: {smooth_row['AIC']:.2f}")

print("\nSmooth terms for s(Trans_ratio) by diff_ecosystem:")
f_col_candidates = [c for c in edf_df.columns if c.strip().lower() in ['f', 'f statistic', 'f_statistic']]
if not f_col_candidates:
    f_col_candidates = [c for c in edf_df.columns if c not in ['term', 'edf', 'ref.df'] and 'p' not in c.lower()]
f_col = f_col_candidates[0] if f_col_candidates else None
p_col_candidates = [c for c in edf_df.columns if 'p' in c.lower()]
p_col = p_col_candidates[0] if p_col_candidates else None

if f_col is None or p_col is None:
    print("Could not locate F or p columns in edf_df; check column names:")
    print(edf_df.columns.tolist())
else:
    for _, row in edf_df.iterrows():
        term = row['term']
        match = re.search(r"diff_ecosystem(.+)__(.+)", term)
        if match:
            diff_type_raw = match.group(1)
            ecosystem = match.group(2)
            diff_map = {
                "WUE_ET_minus_WUE_T": "WUE_ET - WUE_T",
                "WUE_ET_minus_WUE_E": "WUE_ET - WUE_E",
                "WUE_E_minus_WUE_T": "WUE_E - WUE_T",
            }
            diff_label = diff_map.get(diff_type_raw, diff_type_raw)
        else:
            parts = term.split("__")
            if len(parts) == 2:
                diff_type_raw = parts[0].replace("s(Trans_ratio):diff_ecosystem", "")
                ecosystem = parts[1]
                diff_label = diff_map.get(diff_type_raw, diff_type_raw)
            else:
                diff_label = term
                ecosystem = "unknown"
        edf_val = row['edf']
        f_val = row[f_col]
        p_val = row[p_col]
        print(f"  {diff_label} – {ecosystem}: edf = {edf_val:.3f}, F = {f_val:.3f}, p = {p_val:.3e}")

summary_txt = os.path.join(output_dir, "Q1_final_model_summary.txt")
if os.path.exists(summary_txt):
    with open(summary_txt, 'r') as f:
        txt = f.read()
    site_match = re.search(r"s\(site_name\).*?edf\s*=\s*([\d.]+).*?p-value\s*=\s*([\d.e+-]+)", txt, re.DOTALL)
    month_match = re.search(r"s\(month_f\).*?edf\s*=\s*([\d.]+).*?p-value\s*=\s*([\d.e+-]+)", txt, re.DOTALL)
    if site_match:
        print(f"\nSite random effect: edf = {site_match.group(1)}, p = {site_match.group(2)}")
    else:
        print("\nSite random effect: not found in summary.")
    if month_match:
        print(f"Month random effect: edf = {month_match.group(1)}, p = {month_match.group(2)}")
    else:
        print("Month random effect: not found in summary.")
else:
    print("\nModel summary file not found; cannot extract random effect statistics.")
print("="*60 + "\n")

# ------------------------------------------------------------------------------
# Font settings
# ------------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 28,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
})

# ------------------------------------------------------------------------------
# Function to draw a bracket with significance star
# ------------------------------------------------------------------------------
def add_significance_bracket(ax, x1, x2, y, star, color='black', linewidth=1.5, star_size=20):
    """
    Draw a horizontal bracket with a star above it.
    x1, x2: positions of the two groups (in data coordinates)
    y: vertical position of the bracket (in data coordinates)
    star: string, e.g. '*', '**'
    """
    ax.plot([x1, x2], [y, y], color=color, linewidth=linewidth, clip_on=False)
    tick_len = 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([x1, x1], [y - tick_len/2, y + tick_len/2], color=color, linewidth=linewidth, clip_on=False)
    ax.plot([x2, x2], [y - tick_len/2, y + tick_len/2], color=color, linewidth=linewidth, clip_on=False)
    mid_x = (x1 + x2) / 2
    ax.text(mid_x, y + tick_len, star, ha='center', va='bottom', fontsize=star_size,
            fontweight='bold', color=color)

# ------------------------------------------------------------------------------
# Create figure with increased row spacing and larger top margin
# ------------------------------------------------------------------------------
fig_height = 6 + 7 + 6 + 0.5 * 2
fig = plt.figure(figsize=(24, fig_height))

outer = fig.add_gridspec(
    3, 1,
    height_ratios=[6, 7, 6],
    hspace=0.50,
    left=0.12,
    right=0.95,
    top=0.97,        # increased from 0.96 to give more room at top
    bottom=0.06
)

gs_top = outer[0].subgridspec(1, 3, wspace=0.30)
gs_mid = outer[1].subgridspec(1, 3, wspace=0.35)
gs_bot = outer[2].subgridspec(1, 2, wspace=0.30)

axes_a = [fig.add_subplot(gs_top[0, i]) for i in range(3)]
axes_b = [fig.add_subplot(gs_mid[0, i]) for i in range(3)]
ax_g = fig.add_subplot(gs_bot[0, 0])
ax_h = fig.add_subplot(gs_bot[0, 1])

def add_panel_label(ax, label):
    pos = ax.get_position()
    x0, y0, x1, y1 = pos.x0, pos.y0, pos.x1, pos.y1
    fig.text(x0 - 0.015, y1 + 0.01, label, fontsize=32, fontweight='bold',
             va='bottom', ha='left')

# ------------------------------------------------------------------------------
# Panel A (a, b, c) – no changes
# ------------------------------------------------------------------------------
for i, diff_label in enumerate(diff_labels):
    ax = axes_a[i]
    subset = df_summary[df_summary["difference_label"] == diff_label]
    x_positions = np.arange(len(ecosystem_order))
    means, cis = [], []
    for eco in ecosystem_order:
        row = subset[subset["water_class"] == eco]
        if len(row) > 0:
            means.append(row["mean_difference"].values[0])
            cis.append(row["ci95_difference"].values[0])
        else:
            means.append(np.nan)
            cis.append(np.nan)
    for j, eco in enumerate(ecosystem_order):
        if not np.isnan(means[j]):
            ax.errorbar(x_positions[j], means[j], yerr=cis[j], fmt='o',
                        color=ecosystem_colours[eco], capsize=8, elinewidth=3,
                        markersize=12, markeredgecolor='none', capthick=3)
    ax.axhline(0, color='red', linestyle='--', linewidth=3, alpha=0.8)
    ax.set_xticks(x_positions)
    tick_labels = ax.set_xticklabels(ecosystem_order, fontsize=26)
    for tick, eco in zip(tick_labels, ecosystem_order):
        tick.set_color(ecosystem_colours[eco])
    ax.set_ylabel(diff_label_mapping[diff_label], fontsize=24, labelpad=10)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)
    add_panel_label(ax, f"{chr(97+i)})")

# ------------------------------------------------------------------------------
# Panel B (d, e, f) – with reference line
# ------------------------------------------------------------------------------
for i, diff_label in enumerate(diff_labels):
    ax = axes_b[i]
    subset = df_smooth[df_smooth["difference_label"] == diff_label]
    for eco in ecosystem_order:
        eco_sub = subset[subset["ecosystem_label"] == eco]
        if eco_sub.empty:
            continue
        eco_sub = eco_sub.sort_values("Trans_ratio")
        x = eco_sub["Trans_ratio"].values
        y = eco_sub["prediction"].values
        lower = eco_sub["lower"].values
        upper = eco_sub["upper"].values
        color = ecosystem_colours[eco]
        ax.fill_between(x, lower, upper, color=color, alpha=0.15)
        ax.plot(x, y, color=color, linewidth=3, label=eco)
    x_ref, y_ref = ref_curves[diff_label]
    ax.plot(x_ref, y_ref, color='black', linestyle='--', linewidth=3, label='Reference')
    ax.axhline(0, color='red', linestyle='--', linewidth=3, alpha=0.8)
    ax.set_xlabel("T:ET ratio", fontsize=28)
    ax.set_ylabel(diff_label_mapping[diff_label], fontsize=28, labelpad=10)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)
    add_panel_label(ax, f"{chr(100 + i)})")

# Shared legend for d-f
leg = axes_b[2].legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=24,
                       frameon=True, edgecolor='black', title='')
leg.get_frame().set_linewidth(1.5)

# ------------------------------------------------------------------------------
# Panel g: Ecosystem – violin plot with bracket for Upland-Freshwater
# ------------------------------------------------------------------------------
positions = np.arange(len(eco_order_gh))
data_violin = [wue_monthly[wue_monthly['water_class'] == e]['WUE_T_capped'].dropna().values for e in eco_order_gh]

violin_parts = ax_g.violinplot(data_violin, positions=positions, showmeans=False, showmedians=False, widths=0.8)
for i, pc in enumerate(violin_parts['bodies']):
    pc.set_facecolor(ecosystem_colours[eco_order_gh[i]])
    pc.set_edgecolor(ecosystem_colours[eco_order_gh[i]])
    pc.set_alpha(0.7)

for i, eco in enumerate(eco_order_gh):
    vals = wue_monthly[wue_monthly['water_class'] == eco]['WUE_T_capped'].dropna().values
    x_jitter = np.random.normal(i, 0.06, size=len(vals))
    ax_g.scatter(x_jitter, vals, alpha=0.5, s=20, color=ecosystem_colours[eco],
                 edgecolors='none')

ax_g.set_xlim(-0.5, len(eco_order_gh) - 0.5)
max_data = max([np.max(vals) for vals in data_violin]) if data_violin else 1
max_label = eco_means.max() if not eco_means.empty else 1
ylim_top = max(max_data, max_label) * 1.15
ax_g.set_ylim(0, ylim_top)
ax_g.locator_params(axis='y', nbins=3)

# Medians and means
for i, eco in enumerate(eco_order_gh):
    med = eco_medians[eco]
    ax_g.plot([i - 0.15, i + 0.15], [med, med],
              color='black', linewidth=2.5, solid_capstyle='butt')

ymin_g, ymax_g = ax_g.get_ylim()
for i, eco in enumerate(eco_order_gh):
    mean_val = eco_means[eco]
    ax_g.text(i, ymax_g - 0.02 * (ymax_g - ymin_g), f"{mean_val:.2f}",
              ha='center', va='top', fontsize=24, color='black', fontweight='bold')

# Omnibus p-value box
ax_g.text(0.45, 0.80, p_text_eco, transform=ax_g.transAxes,
          fontsize=24, ha='right', va='top',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

# Bracket for Upland-Freshwater with '*'
bracket_y_g = ymax_g + 0.03 * (ymax_g - ymin_g)
add_significance_bracket(ax_g, x1=0, x2=1, y=bracket_y_g, star='*', color='black', linewidth=1.5, star_size=22)

ax_g.set_xticks(positions)
ax_g.set_xticklabels(eco_order_gh)
for tick, eco in zip(ax_g.get_xticklabels(), eco_order_gh):
    tick.set_color(ecosystem_colours[eco])

ax_g.set_ylabel('WUE$_T$ (g C kg$^{-1}$ H$_2$O$^{-1}$)', fontsize=28, labelpad=10)
for spine in ax_g.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)
ax_g.tick_params(axis='both', colors='black', width=2, length=8)
ax_g.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
add_panel_label(ax_g, "g)")

# ------------------------------------------------------------------------------
# Panel h: Coast – mean ± CI with brackets (separated vertically)
# ------------------------------------------------------------------------------
x_pos_h = np.arange(len(coast_order))

for i, coast in enumerate(coast_order):
    vals = tet_site[tet_site['coast_region'] == coast]['mean_TET'].dropna().values
    x_jitter = np.random.normal(i, 0.06, size=len(vals))
    ax_h.scatter(x_jitter, vals, alpha=0.5, s=30, color=coast_colors[coast],
                 edgecolors='none')

for i, coast in enumerate(coast_order):
    ax_h.errorbar(i, coast_means[coast], yerr=coast_ci[coast],
                  fmt='o', color=coast_colors[coast],
                  markersize=16, capsize=10, elinewidth=3,
                  markeredgecolor='black', markeredgewidth=1.2)

ymin_h, ymax_h = 0.2, 0.8
ax_h.set_ylim(ymin_h, ymax_h)
for i, coast in enumerate(coast_order):
    ax_h.text(i, ymax_h - 0.1*(ymax_h-ymin_h), f"{coast_means[coast]:.2f}",
              ha='center', va='bottom', fontsize=24, color='black', fontweight='bold')

ax_h.set_xticks(x_pos_h)
ax_h.set_xticklabels([coast_short[c] for c in coast_order])
for tick, coast in zip(ax_h.get_xticklabels(), coast_order):
    tick.set_color(coast_colors[coast])

ax_h.axhline(0.5, color='gray', linestyle='--', linewidth=2, alpha=0.6)
ax_h.set_ylabel('T:ET ratio', fontsize=28, labelpad=10)

# Omnibus p-value box
ax_h.text(0.97, 0.05, p_text_coast, transform=ax_h.transAxes,
          fontsize=24, ha='right', va='bottom',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

# Brackets for coast comparisons – separated vertically
# Alaska–Atlantic (**) – upper bracket
bracket_y_h = ymax_h + 0.10 * (ymax_h - ymin_h)
add_significance_bracket(ax_h, x1=0, x2=3, y=bracket_y_h, star='**', color='black', linewidth=1.5, star_size=22)

# Pacific–Atlantic (*) – lower bracket
bracket_y_h2 = ymax_h + 0.05 * (ymax_h - ymin_h)
add_significance_bracket(ax_h, x1=1, x2=3, y=bracket_y_h2, star='*', color='black', linewidth=1.5, star_size=22)

for spine in ax_h.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)
ax_h.tick_params(axis='both', colors='black', width=2, length=8)
ax_h.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
add_panel_label(ax_h, "h)")

# ------------------------------------------------------------------------------
# Save figure with bbox_inches='tight' to prevent cut-off
# ------------------------------------------------------------------------------
fig.savefig(fig_output, dpi=600, bbox_inches='tight', facecolor='white')
print(f"\nCombined figure saved to: {fig_output}")

# ------------------------------------------------------------------------------
# Suggested caption wording
# ------------------------------------------------------------------------------
print("\nSuggested Figure 3 caption wording for panels d–f:")
print("Panels d–f show GAM-predicted relationships between T:ET and the observed")
print("WUE-metric differences, with shaded bands indicating 95% confidence intervals.")
print("Black dashed lines show the mathematically expected WUE difference across T:ET")
print("when WUE_ET is held at its overall mean under near-normal conditions.")
print("\nPost-hoc comparisons (Holm-adjusted):")
print("Panel g: * Upland vs Freshwater (p = 0.0173)")
print("Panel h: ** Alaska vs Atlantic (p = 0.0090), * Pacific vs Atlantic (p = 0.0107)")