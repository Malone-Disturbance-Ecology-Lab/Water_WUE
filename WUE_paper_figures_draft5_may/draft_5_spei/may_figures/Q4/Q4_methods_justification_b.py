# -*- coding: utf-8 -*-
"""
CHUNK 2: PERCENTILE THRESHOLD SENSITIVITY - FIGURE CREATION (FINAL COMPLETE)
=============================================================================
Requirements:
- Panel (a): x-axis ticks ~5-6 values, y-axis limit 0.35 to 0.65, legend lower left
- Panel (b): Order thresholds by lowest to highest p-value (5/95, 10/90, 15/85, 20/80, 25/75)
- Panel (a): Same order for curves (lowest to highest p-value)
- Increased height of both panels
=============================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# =============================================================================
# LOAD RESULTS FROM CHUNK 1
# =============================================================================

OUTPUT_BASE_DIR = r"M:\Research\WUE_CUE\data_products\results\sensitivity_analysis_percentiles"

# Load master summary (now includes Intercept)
master_df = pd.read_csv(os.path.join(OUTPUT_BASE_DIR, "Sensitivity_Summary_AllThresholds.csv"))

# Check if Intercept column exists
if 'Intercept' not in master_df.columns:
    raise ValueError("Intercept column not found! Please re-run Chunk 1 with the updated version.")

# Define thresholds in original order
thresholds_original = ['5/95', '10/90', '15/85', '20/80', '25/75']

# Extract data into dataframe for sorting
data_list = []
for thresh in thresholds_original:
    row = master_df[master_df['Threshold'] == thresh].iloc[0]
    data_list.append({
        'Threshold': thresh,
        'Coefficient': row['Coefficient_Standardized'],
        'Intercept': row['Intercept'],
        'CI_Lower': row['CI_Lower'],
        'CI_Upper': row['CI_Upper'],
        'P_Value': row['P_Value'],
        'ROC_AUC': row['ROC_AUC'],
        'Pct_Decrease': row['Pct_Decrease'],
        'SPEI_Mean': row['SPEI_Mean'],
        'SPEI_Std': row['SPEI_Std']
    })

# Create dataframe and sort by p-value (lowest to highest)
df_results = pd.DataFrame(data_list)
df_results = df_results.sort_values('P_Value').reset_index(drop=True)

# Get sorted thresholds
thresholds_order = df_results['Threshold'].tolist()

print("\n" + "="*60)
print("THRESHOLDS ORDERED BY P-VALUE (lowest to highest):")
print("="*60)
for i, thresh in enumerate(thresholds_order):
    p_val = df_results[df_results['Threshold'] == thresh]['P_Value'].values[0]
    print(f"  {i+1}. {thresh}: p = {p_val:.2e}")

# Extract sorted data
threshold_labels = thresholds_order
coefs = df_results['Coefficient'].tolist()
intercepts = df_results['Intercept'].tolist()
ci_lowers = df_results['CI_Lower'].tolist()
ci_uppers = df_results['CI_Upper'].tolist()
p_values = df_results['P_Value'].tolist()
roc_aucs = df_results['ROC_AUC'].tolist()
pct_decrease = df_results['Pct_Decrease'].tolist()
spei_means = df_results['SPEI_Mean'].tolist()
spei_stds = df_results['SPEI_Std'].tolist()

# Print intercepts for verification
print("\n" + "="*60)
print("VERIFICATION: Actual fitted intercepts from GLM")
print("="*60)
for i, thresh in enumerate(threshold_labels):
    print(f"  {thresh}: Intercept = {intercepts[i]:.4f}")

# =============================================================================
# DEFINE COLORS AND LINE STYLES (consistent across both panels)
# =============================================================================
colors = {
    '5/95': '#1f77b4',   # blue
    '10/90': '#ff7f0e',  # orange
    '15/85': '#2ca02c',  # green
    '20/80': '#d62728',  # red
    '25/75': '#9467bd'   # purple
}

line_styles = {
    '5/95': '-',
    '10/90': '--',
    '15/85': '-.',
    '20/80': ':',
    '25/75': (0, (3, 1, 1, 1))
}

# =============================================================================
# FUNCTION: Get probability curve from ACTUAL fitted model parameters
# =============================================================================

def get_probability_correct(coef, intercept, spei_mean, spei_std, n_points=200):
    """
    Generate probability curve using ACTUAL fitted GLM parameters.
    P = 1 / (1 + exp(-(intercept + coef * z)))
    where z = (SPEI - spei_mean) / spei_std
    """
    spei_range = np.linspace(-3.5, 3.5, n_points)
    z_spei = (spei_range - spei_mean) / spei_std
    log_odds = intercept + coef * z_spei
    probs = 1 / (1 + np.exp(-log_odds))
    return spei_range, probs

# =============================================================================
# CREATE FIGURE - INCREASED HEIGHT
# =============================================================================

print("\n" + "="*60)
print("CREATING SUPPLEMENTARY FIGURE")
print("="*60)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))  # Increased height from 7 to 8

# =============================================================================
# PANEL (a) - Overlay probability curves (sorted by p-value)
# =============================================================================

print("\nGenerating Panel (a): Probability curves (ordered by p-value)...")

for i, thresh in enumerate(thresholds_order):
    # Get index for this threshold in the sorted lists
    idx = threshold_labels.index(thresh)
    spei_range, probs = get_probability_correct(
        coefs[idx], intercepts[idx], spei_means[idx], spei_stds[idx]
    )
    
    # Make 10/90 line thicker for emphasis
    line_width = 3.0 if thresh == '10/90' else 2.5
    
    ax1.plot(spei_range, probs, 
             color=colors[thresh], 
             linestyle=line_styles[thresh],
             linewidth=line_width, 
             label=f'{thresh}',
             alpha=0.9)
    
    # Print verification at SPEI=0
    z0 = (0 - spei_means[idx]) / spei_stds[idx]
    p_at_zero = 1 / (1 + np.exp(-(intercepts[idx] + coefs[idx] * z0)))
    print(f"  {thresh}: P(Decrease|SPEI=0) = {p_at_zero:.3f} (intercept={intercepts[idx]:.4f})")

ax1.axvspan(-0.5, 0.5, alpha=0.08, color='gray', label='Near-normal region')
ax1.axhline(0.5, color='black', linestyle='--', alpha=0.5, linewidth=1.2)
ax1.set_xlabel('SPEI-48', fontsize=16, fontweight='bold')
ax1.set_ylabel('P(WUE$_T$ decrease)', fontsize=16, fontweight='bold')
ax1.set_title('(a) Sensitivity of WUE$_T$ decline probability\nto threshold', 
              fontsize=14, fontweight='bold', loc='left')
ax1.set_ylim(0.35, 0.65)  # Changed from 0.30-0.70 to 0.35-0.65
ax1.set_xlim(-3.2, 3.2)

# Set x-axis ticks to 7 values
ax1.set_xticks([-3, -2, -1, 0, 1, 2, 3])
ax1.set_xticklabels(['-3', '-2', '-1', '0', '1', '2', '3'], fontsize=24)

ax1.grid(True, alpha=0.15, linestyle='--')
ax1.legend(loc='lower left', fontsize=14, framealpha=0.9, ncol=2)
ax1.tick_params(axis='both', labelsize=24, width=1.5, length=8)

# =============================================================================
# PANEL (b) - Coefficient robustness (sorted by p-value)
# =============================================================================

print("\nGenerating Panel (b): Coefficient robustness (ordered by p-value)...")

x_pos = np.arange(len(threshold_labels))

# Plot error bars with middle dots
for i, thresh in enumerate(threshold_labels):
    idx = threshold_labels.index(thresh)
    ax2.errorbar(x_pos[i], coefs[idx], 
                 yerr=[[coefs[idx] - ci_lowers[idx]], [ci_uppers[idx] - coefs[idx]]],
                 fmt='o', 
                 color=colors[thresh],
                 capsize=8, 
                 capthick=2.5,
                 markersize=12, 
                 markeredgecolor='black', 
                 markeredgewidth=1.5,
                 elinewidth=2.5, 
                 ecolor=colors[thresh],
                 alpha=0.9,
                 zorder=3)

ax2.axhline(0, color='red', linestyle='--', alpha=0.6, linewidth=1.5, zorder=1)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(threshold_labels, fontsize=24)
ax2.set_ylabel('Standardized β coefficient\n(per 1 SD change in SPEI-48)', fontsize=16, fontweight='bold')
ax2.set_title('(b) Coefficient robustness across thresholds', fontsize=14, fontweight='bold', loc='left')
ax2.grid(True, alpha=0.15, linestyle='--', axis='y', zorder=0)
ax2.tick_params(axis='both', labelsize=24, width=1.5, length=8)

# Set y-axis limits with extra space
y_min = min(ci_lowers) - 0.22
y_max = max(ci_uppers) + 0.35
ax2.set_ylim(y_min, y_max)

# =============================================================================
# ADD VALUE LABELS WITH CUSTOM POSITIONING
# =============================================================================

for i, thresh in enumerate(threshold_labels):
    idx = threshold_labels.index(thresh)
    x = x_pos[i]
    coef = coefs[idx]
    p_val = p_values[idx]
    roc_auc = roc_aucs[idx]
    pct = pct_decrease[idx]
    
    # Format β with 3 significant digits for 5/95, 2 digits for others
    if thresh == '5/95':
        beta_text = f'{coef:.3f}'  # 3 significant digits for 5/95
    else:
        beta_text = f'{coef:.2f}'  # 2 digits for others
    
    # Add asterisk for significance (p < 0.05)
    if p_val < 0.05:
        beta_text += '*'
    
    # Format AUC text
    if thresh == '10/90':
        auc_text = f'AUC = {roc_auc:.1f}'
    else:
        auc_text = f'AUC = {roc_auc:.2f}'
    
    # Format p-value text with specific requirements
    if thresh == '5/95':
        if p_val < 0.001:
            p_text = 'p < 0.001'
        elif p_val < 0.01:
            p_text = f'p = {p_val:.3f}'
        else:
            p_text = f'p = {p_val:.2f}'
    elif thresh == '10/90':
        if p_val < 0.001:
            p_text = 'p < 0.001'
        else:
            p_text = f'p = {p_val:.3f}'
    else:
        if p_val < 0.001:
            p_text = 'p < 0.001'
        elif p_val < 0.01:
            p_text = f'p = {p_val:.3f}'
        else:
            p_text = f'p = {p_val:.2f}'
    
    # Custom positioning based on threshold
    if thresh == '5/95':
        x_offset = +0.15
        y_beta = coef + 0.1
        y_auc = coef + 0.16
        y_p = coef - 0.14
        ha_align = 'right'
    elif thresh == '25/75':
        x_offset = +0.18
        y_beta = coef + 0.12
        y_auc = coef + 0.16
        y_p = coef - 0.1
        ha_align = 'right'
    elif thresh == '10/90':
        x_offset = -0.15
        y_beta = coef + 0.14
        y_auc = coef + 0.18
        y_p = coef - 0.12
        ha_align = 'left'
    else:
        x_offset = 0
        y_beta = coef + 0.12
        y_auc = coef + 0.18
        y_p = coef - 0.1
        ha_align = 'center'
    
    # Beta label
    ax2.annotate(beta_text, (x + x_offset, y_beta), 
                 ha=ha_align, va='bottom', 
                 fontsize=13, fontweight='bold',
                 color=colors[thresh],
                 zorder=4)
    
    # AUC label (black)
    ax2.annotate(auc_text, (x + x_offset, y_auc), 
                 ha=ha_align, va='bottom',
                 fontsize=12, color='black', fontweight='bold',
                 zorder=4)
    
    # p-value label (black)
    ax2.annotate(p_text, (x + x_offset, y_p), 
                 ha=ha_align, va='top',
                 fontsize=12, color='black', fontweight='bold',
                 zorder=4)

plt.tight_layout()

# =============================================================================
# SAVE FIGURE
# =============================================================================

fig_path_png = os.path.join(OUTPUT_BASE_DIR, "Supplementary_Figure_PercentileSensitivity.png")
fig_path_pdf = os.path.join(OUTPUT_BASE_DIR, "Supplementary_Figure_PercentileSensitivity.pdf")
fig.savefig(fig_path_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(fig_path_pdf, bbox_inches='tight', facecolor='white')
print(f"\n  ✅ Saved: {fig_path_png}")
print(f"  ✅ Saved: {fig_path_pdf}")

# Display in console
plt.show()

# =============================================================================
# PRINT SUMMARY
# =============================================================================

print("\n" + "="*60)
print("FIGURE SUMMARY")
print("="*60)
print("\nThreshold | β (std) | Intercept | 95% CI | p-value | ROC-AUC | % Decrease")
print("-" * 85)
for i, thresh in enumerate(threshold_labels):
    idx = threshold_labels.index(thresh)
    sig_star = '*' if p_values[idx] < 0.05 else ''
    
    # Format beta for display
    if thresh == '5/95':
        beta_display = f"{coefs[idx]:.3f}"
    else:
        beta_display = f"{coefs[idx]:.2f}"
    
    # Format p-value for display
    if thresh == '5/95':
        if p_values[idx] < 0.001:
            p_display = "<0.001"
        else:
            p_display = f"{p_values[idx]:.2f}"
    elif thresh == '10/90':
        if p_values[idx] < 0.001:
            p_display = "<0.001"
        else:
            p_display = f"{p_values[idx]:.3f}"
    else:
        if p_values[idx] < 0.001:
            p_display = "<0.001"
        elif p_values[idx] < 0.01:
            p_display = f"{p_values[idx]:.3f}"
        else:
            p_display = f"{p_values[idx]:.2f}"
    
    print(f"{thresh:8} | {beta_display}{sig_star} | {intercepts[idx]:8.4f} | [{ci_lowers[idx]:5.3f}, {ci_uppers[idx]:5.3f}] | {p_display:10} | {roc_aucs[idx]:5.3f} | {pct_decrease[idx]:5.1f}%")

print("\n" + "="*60)
print("✅ CHUNK 2 COMPLETE - All formatting requirements applied")
print("="*60)
print("\nFORMATTING REQUIREMENTS MET:")
print("   PANEL (a):")
print("   1. ✓ X-axis ticks: 7 values (-3, -2, -1, 0, 1, 2, 3)")
print("   2. ✓ Y-axis limit: 0.35 to 0.65")
print("   3. ✓ Legend position: lower left")
print("   4. ✓ Curves ordered by p-value (lowest to highest)")
print("\n   PANEL (b):")
print("   5. ✓ Thresholds ordered by p-value (lowest to highest)")
print("   6. ✓ 5/95: Beta with 3 significant digits, p-value with 2 digits")
print("   7. ✓ 10/90: 3 significant digits for p-value")
print("   8. ✓ Error bars with middle dots colored by threshold")
print("   9. ✓ AUC and p-value text in BLACK")
print("   10. ✓ X and Y tick font size DOUBLED (24pt)")
print("   11. ✓ Figure height increased (8 instead of 7)")
print("="*60)