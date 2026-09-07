#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q3_ERA5_reviewer_chunk2_figure.py

Loads the cached intermediate data from Chunk 1 and regenerates:
- Figure_Sx_ERA5_bias_correction.jpg (with improved formatting, percent bias)
- Table_Sx_ERA5_correction_performance_by_coast.csv
- Table_Sy_ERA5_fill_fraction_by_coast.csv

Run this after Chunk 1 to adjust figure formatting without reprocessing.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================
OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods\reviewer_ERA5_correction"
CACHE_DIR = os.path.join(OUTPUT_DIR, "cache")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# LOAD CACHED DATA
# ============================================================================
print("Loading cached data...")
df_tair = pd.read_csv(os.path.join(CACHE_DIR, 'df_tair.csv'))
df_rg = pd.read_csv(os.path.join(CACHE_DIR, 'df_rg.csv'))
df_vpd = pd.read_csv(os.path.join(CACHE_DIR, 'df_vpd.csv'))
df_fill = pd.read_csv(os.path.join(CACHE_DIR, 'df_fill.csv'))
final_domain = pd.read_csv(os.path.join(CACHE_DIR, 'final_domain.csv'))

# Load pooled arrays
tower_tair = np.load(os.path.join(CACHE_DIR, 'tower_tair.npy'))
raw_tair = np.load(os.path.join(CACHE_DIR, 'raw_tair.npy'))
corr_tair = np.load(os.path.join(CACHE_DIR, 'corr_tair.npy'))
tower_rg = np.load(os.path.join(CACHE_DIR, 'tower_rg.npy'))
raw_rg = np.load(os.path.join(CACHE_DIR, 'raw_rg.npy'))
corr_rg = np.load(os.path.join(CACHE_DIR, 'corr_rg.npy'))
tower_vpd = np.load(os.path.join(CACHE_DIR, 'tower_vpd.npy'))
raw_vpd = np.load(os.path.join(CACHE_DIR, 'raw_vpd.npy'))
corr_vpd = np.load(os.path.join(CACHE_DIR, 'corr_vpd.npy'))

# Add variable columns to dataframes
for df in [df_tair, df_rg, df_vpd]:
    if 'variable' not in df.columns:
        if df is df_tair:
            df['variable'] = 'Tair'
        elif df is df_rg:
            df['variable'] = 'Rg'
        else:
            df['variable'] = 'VPD'

# ============================================================================
# HELPER: compute metrics (for pooled text only)
# ============================================================================
from scipy.stats import linregress
def compute_metrics(obs, pred):
    mask = np.isfinite(obs) & np.isfinite(pred)
    if mask.sum() < 2:
        return np.nan, np.nan, np.nan, 0, np.nan
    obs = obs[mask]; pred = pred[mask]
    slope, intercept, r_value, _, _ = linregress(obs, pred)
    r2 = r_value**2
    bias = np.mean(pred - obs)
    rmse = np.sqrt(np.mean((pred - obs)**2))
    mean_obs = np.mean(obs)
    percent_bias = 100 * bias / mean_obs if mean_obs != 0 else np.nan
    return r2, bias, rmse, len(obs), percent_bias

r2_tair_raw, bias_tair_raw, _, _, pct_bias_tair_raw = compute_metrics(tower_tair, raw_tair)
r2_tair_corr, bias_tair_corr, _, _, pct_bias_tair_corr = compute_metrics(tower_tair, corr_tair)
r2_rg_raw, bias_rg_raw, _, _, pct_bias_rg_raw = compute_metrics(tower_rg, raw_rg)
r2_rg_corr, bias_rg_corr, _, _, pct_bias_rg_corr = compute_metrics(tower_rg, corr_rg)
r2_vpd_raw, bias_vpd_raw, _, _, pct_bias_vpd_raw = compute_metrics(tower_vpd, raw_vpd)
r2_vpd_corr, bias_vpd_corr, _, _, pct_bias_vpd_corr = compute_metrics(tower_vpd, corr_vpd)

# Format numbers to 2 significant digits (for figure)
def fmt_sig(x):
    if np.isnan(x):
        return "NA"
    return f"{x:.2g}"

# Format numbers to 3 significant digits (for tables)
def fmt_sig3(x):
    if np.isnan(x):
        return "NA"
    return f"{x:.3g}"

# ============================================================================
# GENERATE FIGURE WITH IMPROVED FORMATTING (PERCENT BIAS)
# ============================================================================
print("Generating figure...")

# Use colorblind-friendly colors
COLOR_RAW = '#0072B2'    # blue
COLOR_CORR = '#D55E00'   # orange

# Larger font sizes
FONT_TICK = 22
FONT_LABEL = 26
FONT_PANEL = 26
FONT_LEGEND = 18
FONT_STATS = 16

fig = plt.figure(figsize=(20, 18))
gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.3, hspace=0.3)

def scatter_panel(ax, obs, raw, corr, xlabel, ylabel, panel_label,
                  r2_raw, pct_bias_raw, r2_corr, pct_bias_corr):
    # Scatter points with transparency
    ax.scatter(raw, obs, s=5, alpha=0.15, color=COLOR_RAW, label='Raw ERA5')
    ax.scatter(corr, obs, s=5, alpha=0.15, color=COLOR_CORR, label='Corrected ERA5')
    # 1:1 line
    lims = [min(obs.min(), raw.min(), corr.min()), max(obs.max(), raw.max(), corr.max())]
    ax.plot(lims, lims, 'k--', linewidth=1.5, alpha=0.6)
    ax.set_xlabel(xlabel, fontsize=FONT_LABEL)
    ax.set_ylabel(ylabel, fontsize=FONT_LABEL)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    # Increase tick label size
    ax.tick_params(axis='both', labelsize=FONT_TICK)
    # Panel label
    ax.text(0.03, 0.95, panel_label, transform=ax.transAxes, fontsize=FONT_PANEL,
            fontweight='bold', va='top', ha='left')
    # Stats text with percent bias – moved higher (y=0.97) and larger font
    stats_text = (f"Before: R²={fmt_sig(r2_raw)}, Bias={fmt_sig(pct_bias_raw)}%\n"
                  f"After:  R²={fmt_sig(r2_corr)}, Bias={fmt_sig(pct_bias_corr)}%")
    ax.text(0.97, 0.97, stats_text, transform=ax.transAxes,
            fontsize=FONT_STATS, va='top', ha='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.85, edgecolor='black', linewidth=1.5))

# Panel (a): Tair
ax_a = fig.add_subplot(gs[0, 0])
scatter_panel(ax_a, tower_tair, raw_tair, corr_tair,
              'ERA5 Tair (°C)', 'AmeriFlux Tair (°C)', '(a)',
              r2_tair_raw, pct_bias_tair_raw, r2_tair_corr, pct_bias_tair_corr)

# Panel (b): Rg
ax_b = fig.add_subplot(gs[0, 1])
scatter_panel(ax_b, tower_rg, raw_rg, corr_rg,
              'ERA5 Rg (W m⁻²)', 'AmeriFlux Rg (W m⁻²)', '(b)',
              r2_rg_raw, pct_bias_rg_raw, r2_rg_corr, pct_bias_rg_corr)

# Panel (c): VPD
ax_c = fig.add_subplot(gs[0, 2])
scatter_panel(ax_c, tower_vpd, raw_vpd, corr_vpd,
              'ERA5 VPD (hPa)', 'AmeriFlux VPD (hPa)', '(c)',
              r2_vpd_raw, pct_bias_vpd_raw, r2_vpd_corr, pct_bias_vpd_corr)

# Panel (d): Percent RMSE reduction
ax_d = fig.add_subplot(gs[1, 0])
reductions_tair = []
for _, row in df_tair.iterrows():
    rmse_b = row['rmse_before']
    rmse_a = row['rmse_after']
    if not np.isnan(rmse_b) and not np.isnan(rmse_a) and rmse_b > 0:
        reductions_tair.append(100 * (rmse_b - rmse_a) / rmse_b)

reductions_rg = []
for _, row in df_rg.iterrows():
    rmse_b = row['rmse_before']
    rmse_a = row['rmse_after']
    if not np.isnan(rmse_b) and not np.isnan(rmse_a) and rmse_b > 0:
        reductions_rg.append(100 * (rmse_b - rmse_a) / rmse_b)

reductions_vpd = []
for _, row in df_vpd.iterrows():
    rmse_b = row['rmse_before']
    rmse_a = row['rmse_after']
    if not np.isnan(rmse_b) and not np.isnan(rmse_a) and rmse_b > 0:
        reductions_vpd.append(100 * (rmse_b - rmse_a) / rmse_b)

data_d = [reductions_tair, reductions_rg, reductions_vpd]
labels_d = ['Tair', 'Rg', 'VPD']
bp_d = ax_d.boxplot(data_d, labels=labels_d, patch_artist=True,
                    notch=True, widths=0.6,
                    boxprops=dict(facecolor='lightgray', edgecolor='black', linewidth=1.5),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5))
ax_d.axhline(y=0, color='k', linestyle='--', linewidth=1.5, alpha=0.6)
ax_d.set_ylabel('RMSE reduction (%)', fontsize=FONT_LABEL)
ax_d.text(0.03, 0.95, '(d)', transform=ax_d.transAxes, fontsize=FONT_PANEL,
          fontweight='bold', va='top', ha='left')
ax_d.tick_params(axis='both', labelsize=FONT_TICK)
ax_d.grid(True, axis='y', alpha=0.3)

# Panel (e): Correction slopes
ax_e = fig.add_subplot(gs[1, 1])
slopes_tair = df_tair['slope'].dropna().values
slopes_rg = df_rg['slope'].dropna().values
slopes_vpd = df_vpd['slope'].dropna().values
data_e = [slopes_tair, slopes_rg, slopes_vpd]
bp_e = ax_e.boxplot(data_e, labels=labels_d, patch_artist=True,
                    notch=True, widths=0.6,
                    boxprops=dict(facecolor='lightgray', edgecolor='black', linewidth=1.5),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5))
ax_e.axhline(y=1, color='k', linestyle='--', linewidth=1.5, alpha=0.6)
ax_e.set_ylabel('Correction slope', fontsize=FONT_LABEL)
ax_e.text(0.03, 0.95, '(e)', transform=ax_e.transAxes, fontsize=FONT_PANEL,
          fontweight='bold', va='top', ha='left')
ax_e.tick_params(axis='both', labelsize=FONT_TICK)
ax_e.grid(True, axis='y', alpha=0.3)

# Panel (f): Fraction filled by ERA5 (site-level) – cap visually at 98th percentile
ax_f = fig.add_subplot(gs[1, 2])
fill_tair = df_fill['Tair_fill_frac'].dropna().values
fill_rg = df_fill['Rg_fill_frac'].dropna().values
fill_vpd = df_fill['VPD_fill_frac'].dropna().values
data_f = [fill_tair, fill_rg, fill_vpd]
# Compute 98th percentile for each to set ylim
p98_tair = np.percentile(fill_tair, 98) if len(fill_tair)>0 else 1
p98_rg = np.percentile(fill_rg, 98) if len(fill_rg)>0 else 1
p98_vpd = np.percentile(fill_vpd, 98) if len(fill_vpd)>0 else 1
ylim_upper = max(p98_tair, p98_rg, p98_vpd) * 1.05

bp_f = ax_f.boxplot(data_f, labels=labels_d, patch_artist=True,
                    notch=True, widths=0.6,
                    boxprops=dict(facecolor='lightgray', edgecolor='black', linewidth=1.5),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5))
ax_f.set_ylabel('Fraction filled by ERA5', fontsize=FONT_LABEL)
ax_f.set_ylim(0, ylim_upper)
ax_f.text(0.03, 0.95, '(f)', transform=ax_f.transAxes, fontsize=FONT_PANEL,
          fontweight='bold', va='top', ha='left')
ax_f.tick_params(axis='both', labelsize=FONT_TICK)
ax_f.grid(True, axis='y', alpha=0.3)

# Legend for scatter points
handles_scatter = [Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_RAW, markersize=8, label='Raw ERA5'),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_CORR, markersize=8, label='Corrected ERA5')]
ax_a.legend(handles=handles_scatter, loc='lower right', fontsize=FONT_LEGEND, framealpha=0.9, edgecolor='black')

plt.tight_layout()
fig_path = os.path.join(OUTPUT_DIR, 'Figure_Sx_ERA5_bias_correction.jpg')
plt.savefig(fig_path, dpi=600, bbox_inches='tight', format='jpg')
plt.close()
print(f"Figure saved: {fig_path}")

# ============================================================================
# GENERATE COAST-BASED TABLES WITH 3 SIGNIFICANT DIGITS
# ============================================================================
print("Generating coast-based tables...")

def format_med_iqr(x):
    """Format a string like 'median [lower, upper]' with 3 significant digits."""
    if pd.isna(x) or x == "NA":
        return "NA"
    # Expect format: "1.234 [0.123, 0.456]"
    import re
    parts = re.split(r' \[|\]|, ', x)
    if len(parts) >= 3:
        median = float(parts[0])
        lower = float(parts[1])
        upper = float(parts[2])
        return f"{fmt_sig3(median)} [{fmt_sig3(lower)}, {fmt_sig3(upper)}]"
    else:
        return x

def format_range(x):
    """Format a range like 'min–max' with 3 significant digits."""
    if pd.isna(x) or x == "NA":
        return "NA"
    if '–' in x:
        parts = x.split('–')
        if len(parts) == 2:
            try:
                a = float(parts[0])
                b = float(parts[1])
                return f"{fmt_sig3(a)}–{fmt_sig3(b)}"
            except:
                return x
    return x

def agg_coast_performance(df, var_name):
    sub = df[df['variable'] == var_name].copy()
    if sub.empty:
        return None
    grouped = sub.groupby('coast')
    rows = []
    for coast, grp in grouped:
        if len(grp) == 0:
            continue
        def med_iqr_str(x):
            q1 = x.quantile(0.25)
            q3 = x.quantile(0.75)
            return f"{x.median():.3g} [{q1:.3g}, {q3:.3g}]"
        row = {
            'Coast': coast,
            'Variable': var_name,
            'N sites': grp['site'].nunique() if 'site' in grp.columns else len(grp),
            'N units': len(grp),
            'Correction slope, median [IQR]': med_iqr_str(grp['slope']),
            'Correction intercept, median [IQR]': med_iqr_str(grp['intercept']),
            'R² before, median [IQR]': med_iqr_str(grp['r2_before']),
            'R² after, median [IQR]': med_iqr_str(grp['r2_after']),
            'Bias before, median [IQR]': med_iqr_str(grp['bias_before']),
            'Bias after, median [IQR]': med_iqr_str(grp['bias_after']),
            'RMSE before, median [IQR]': med_iqr_str(grp['rmse_before']),
            'RMSE after, median [IQR]': med_iqr_str(grp['rmse_after']),
        }
        rows.append(row)
    return pd.DataFrame(rows)

table_sx_parts = []
for var in ['Tair', 'Rg', 'VPD']:
    if var == 'Tair':
        df_use = df_tair
    elif var == 'Rg':
        df_use = df_rg
    else:
        df_use = df_vpd
    part = agg_coast_performance(df_use, var)
    if part is not None:
        table_sx_parts.append(part)
table_sx = pd.concat(table_sx_parts, ignore_index=True)

def fill_table_by_coast():
    rows = []
    for coast in ['AK Coast', 'Pacific Coast', 'Gulf Coast', 'Atlantic Coast']:
        sub = df_fill[df_fill['coast'] == coast]
        if sub.empty:
            continue
        n_sites = sub['site'].nunique()
        sy = final_domain[final_domain['site_name'].isin(sub['site'])][['site_name','Year']].drop_duplicates()
        n_site_years = len(sy)
        total_records = sub['total_records'].sum()
        for var in ['Tair', 'Rg', 'VPD']:
            fill_col = f'{var}_fill_count'
            frac_col = f'{var}_fill_frac'
            fill_total = sub[fill_col].sum()
            fracs = sub[frac_col].dropna()
            if len(fracs) > 0:
                med = fracs.median()
                q1 = fracs.quantile(0.25)
                q3 = fracs.quantile(0.75)
                med_iqr = f"{fmt_sig3(med)} [{fmt_sig3(q1)}, {fmt_sig3(q3)}]"
                min_frac = fracs.min()
                max_frac = fracs.max()
                range_str = f"{fmt_sig3(min_frac)}-{fmt_sig3(max_frac)}"
            else:
                med_iqr = "NA"
                range_str = "NA"
            zero_sites = sub[sub[frac_col] == 0]['site'].nunique() if len(fracs) > 0 else 0
            complete_sites = sub[sub[frac_col] == 1]['site'].nunique() if len(fracs) > 0 else 0
            rows.append({
                'Coast': coast,
                'Variable': var,
                'N sites': n_sites,
                'N final site-years': n_site_years,
                'N final half-hourly records': total_records,
                'N records filled directly by ERA5': fill_total,
                'ERA5 fill fraction, median [IQR]': med_iqr,
                'ERA5 fill fraction, range': range_str,
                'N sites with zero ERA5 filling': zero_sites,
                'N sites where ERA5 supplied complete variable': complete_sites
            })
    return pd.DataFrame(rows)

table_sy = fill_table_by_coast()

# Save tables with UTF-8 encoding to preserve en dash
table_sx_path = os.path.join(OUTPUT_DIR, 'Table_Sx_ERA5_correction_performance_by_coast.csv')
table_sy_path = os.path.join(OUTPUT_DIR, 'Table_Sy_ERA5_fill_fraction_by_coast.csv')
table_sx.to_csv(table_sx_path, index=False, encoding='utf-8')
table_sy.to_csv(table_sy_path, index=False, encoding='utf-8')
print(f"Table Sx saved: {table_sx_path}")
print(f"Table Sy saved: {table_sy_path}")

print("\nChunk 2 complete.")