"""
Q3_COMPLETE_UPDATED_WORKFLOW.py
Complete Q3 workflow with confirmed values and refined interpretation
Based on helper script outputs and corrected data extraction
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime

print("="*80)
print("Q3 COMPLETE UPDATED WORKFLOW - METHODS & RESULTS EXTRACTION")
print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# ============================================================================
# PATHS
# ============================================================================

base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")
figure_dir = os.path.join(base_output_dir, "Q3_WUE_T_SPEI_sensitivity_figures")

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")

# ============================================================================
# CONSTANTS - CONFIRMED VALUES FROM ANALYSIS
# ============================================================================

# Confirmed dataset stats
N_OBSERVATIONS = 1872
N_SITES = 64
YEAR_START = 1994
YEAR_END = 2025
N_YEARS = 32
N_LONG_OBS = 13104
N_SPEI_TIMESCALES = 7

# Confirmed ecosystem distribution
ECOSYSTEM_DIST = {
    "Upland": {"obs": 793, "sites": 29},
    "Freshwater": {"obs": 629, "sites": 20},
    "Saline": {"obs": 450, "sites": 15}
}

# Confirmed coast distribution
COAST_DIST = {
    "Atlantic Coast": {"obs": 1017, "sites": 29},
    "Pacific Coast": {"obs": 560, "sites": 18},
    "Gulf Coast": {"obs": 110, "sites": 5},
    "AK Coast": {"obs": 185, "sites": 12}
}

# Confirmed model comparison values
MODEL_COMP = {
    "ecosystem": {
        "AIC": 55467.76,
        "BIC": 56263.82,
        "logLik": -27627.46,
        "deviance_explained": 0.579,
        "r_squared": 0.575
    },
    "coast": {
        "AIC": 55347.65,
        "BIC": 56361.12,
        "logLik": -27538.35,
        "deviance_explained": 0.584,
        "r_squared": 0.580
    }
}

# Confirmed strongest signal
STRONGEST_SIGNAL = {
    "coast": "Gulf Coast",
    "timescale": "SPEI_3",
    "edf": 4.34,
    "f_stat": 29.44,
    "p_value": "<0.001"
}

# Confirmed coastline sensitivity medians
COAST_SENSITIVITY = {
    "Alaska": {"SPEI_1": 37.5, "SPEI_3": 43.6, "SPEI_48": 44.7},
    "Atlantic": {"SPEI_1": 18.3, "SPEI_3": 18.0, "SPEI_48": 20.2},
    "Gulf": {"SPEI_1": 39.3, "SPEI_3": 42.5, "SPEI_48": 37.8},
    "Pacific": {"SPEI_1": 17.2, "SPEI_3": 17.7, "SPEI_48": 16.3}
}

# Confirmed threshold ranges
THRESHOLDS = {
    "Gulf_SPEI3_dry_decrease": {"min": -2.13, "max": -2.02},
    "Gulf_SPEI3_wet_increase": {"min": 1.37, "max": 1.81},
    "AK_SPEI48_dry_decrease": {"min": -1.44, "max": -1.44},
    "AK_SPEI48_wet_increase": {"min": 1.44, "max": 1.44},
    "Atlantic_SPEI48_wet_decrease": {"min": 1.20, "max": 1.20},
    "Pacific_SPEI1_dry_increase": {"min": -1.78, "max": -1.03}
}

# ============================================================================
# SECTION 1: METHODS - DATASET DESCRIPTION
# ============================================================================

print("\n" + "="*80)
print("SECTION 1: METHODS - DATASET DESCRIPTION (CONFIRMED VALUES)")
print("="*80)

print(f"\nFinal dataset after Q3 complete-case filtering:")
print(f"  • {N_OBSERVATIONS:,} site-month observations")
print(f"  • {N_SITES} sites")
print(f"  • {YEAR_START} to {YEAR_END} ({N_YEARS} years)")

print(f"\nLong-format dataset (reshaped across SPEI timescales):")
print(f"  • {N_LONG_OBS:,} site-month-timescale observations")
print(f"  • {N_SPEI_TIMESCALES} SPEI timescales (1-48 months)")

print(f"\nObservations by ecosystem class:")
for eco, stats in ECOSYSTEM_DIST.items():
    print(f"  • {eco}: {stats['obs']:,} observations ({stats['sites']} sites)")

print(f"\nObservations by coast region:")
for coast, stats in COAST_DIST.items():
    print(f"  • {coast}: {stats['obs']:,} observations ({stats['sites']} sites)")

# ============================================================================
# SECTION 2: METHODS - MODEL COMPARISON
# ============================================================================

print("\n" + "="*80)
print("SECTION 2: METHODS - MODEL COMPARISON (TABLE Q3.1)")
print("="*80)

print("\nTable Q3.1: Comparison of ecosystem and coast-threshold GAMs")
print("-" * 80)
print(f"{'Model':<20} {'AIC':>10} {'BIC':>10} {'logLik':>10} {'Deviance Expl.':>14} {'Adj R²':>10}")
print("-" * 80)

eco = MODEL_COMP["ecosystem"]
coast = MODEL_COMP["coast"]

print(f"{'Ecosystem GAM':<20} {eco['AIC']:>10.2f} {eco['BIC']:>10.2f} {eco['logLik']:>10.2f} {eco['deviance_explained']*100:>13.1f}% {eco['r_squared']:>10.3f}")
print(f"{'Coast-threshold GAM':<20} {coast['AIC']:>10.2f} {coast['BIC']:>10.2f} {coast['logLik']:>10.2f} {coast['deviance_explained']*100:>13.1f}% {coast['r_squared']:>10.3f}")

print("\nModel Selection Rationale:")
print(f"  • Lower AIC: {coast['AIC']:.2f} vs {eco['AIC']:.2f} (ΔAIC = {eco['AIC'] - coast['AIC']:.2f})")
print(f"  • Higher deviance explained: {coast['deviance_explained']*100:.1f}% vs {eco['deviance_explained']*100:.1f}%")
print(f"  • Higher adjusted R²: {coast['r_squared']:.3f} vs {eco['r_squared']:.3f}")

print("\nMETHODS SENTENCE TEMPLATE:")
print(f"  \"The coast-threshold GAM performed better than the ecosystem GAM (AIC: {coast['AIC']:.2f} vs {eco['AIC']:.2f}, deviance explained: {coast['deviance_explained']*100:.1f}% vs {eco['deviance_explained']*100:.1f}%, adjusted R²: {coast['r_squared']:.3f} vs {eco['r_squared']:.3f}) and was used as the primary model.\"")

# ============================================================================
# SECTION 3: RESULTS - STRONGEST SIGNAL
# ============================================================================

print("\n" + "="*80)
print("SECTION 3: RESULTS - STRONGEST SIGNAL")
print("="*80)

strong = STRONGEST_SIGNAL
print(f"\nStrongest nonlinear SPEI signal:")
print(f"  • Coast: {strong['coast']}")
print(f"  • Timescale: {strong['timescale']}")
print(f"  • EDF (Effective Degrees of Freedom): {strong['edf']:.2f}")
print(f"  • F-statistic: {strong['f_stat']:.2f}")
print(f"  • p-value: {strong['p_value']}")

print("\nRESULTS SENTENCE TEMPLATE:")
print(f"  \"The strongest nonlinear SPEI signal occurred in the {strong['coast']} at {strong['timescale']} (EDF={strong['edf']:.2f}, F={strong['f_stat']:.2f}, p<0.001).\"")

# ============================================================================
# SECTION 4: RESULTS - COASTLINE SENSITIVITY MEDIANS (TABLE Q3.3)
# ============================================================================

print("\n" + "="*80)
print("SECTION 4: RESULTS - COASTLINE SENSITIVITY MEDIANS (TABLE Q3.3)")
print("="*80)

print("\nTable Q3.3: Coastline sensitivity medians for selected SPEI timescales")
print("-" * 80)
print(f"{'Coast':<10} {'SPEI-1 (%)':>12} {'SPEI-3 (%)':>12} {'SPEI-48 (%)':>13}")
print("-" * 80)

for coast, values in COAST_SENSITIVITY.items():
    print(f"{coast:<10} {values['SPEI_1']:>12.1f} {values['SPEI_3']:>12.1f} {values['SPEI_48']:>13.1f}")

# Calculate averages
print("\nCoast patterns (average median across selected timescales):")
for coast, values in COAST_SENSITIVITY.items():
    avg = sum(values.values()) / 3
    print(f"  • {coast}: {avg:.1f}%")

# Highest values
print("\nHighest observed sensitivity:")
max_overall = max(max(v.values()) for v in COAST_SENSITIVITY.values())
for coast, values in COAST_SENSITIVITY.items():
    if max(values.values()) == max_overall:
        ts = max(values, key=values.get)
        print(f"  • Overall: {coast} at {ts} ({max_overall:.1f}%)")

print("\nRESULTS SENTENCE TEMPLATE:")
print(f"  \"Observed site sensitivity was highest in the Alaska Coast (SPEI-48: {COAST_SENSITIVITY['Alaska']['SPEI_48']:.1f}%, SPEI-3: {COAST_SENSITIVITY['Alaska']['SPEI_3']:.1f}%) and Gulf Coast (SPEI-3: {COAST_SENSITIVITY['Gulf']['SPEI_3']:.1f}%, SPEI-1: {COAST_SENSITIVITY['Gulf']['SPEI_1']:.1f}%), while Atlantic and Pacific coasts showed lower sensitivity (approximately 16-20%).\"")

# ============================================================================
# SECTION 5: RESULTS - THRESHOLD SUMMARIES
# ============================================================================

print("\n" + "="*80)
print("SECTION 5: RESULTS - IMPACT THRESHOLDS")
print("="*80)

print("\nGulf Coast SPEI-3 thresholds (strongest signal):")
print(f"  • Dry-side decreases: SPEI {THRESHOLDS['Gulf_SPEI3_dry_decrease']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_dry_decrease']['max']:.2f}")
print(f"  • Wet-side increases: SPEI {THRESHOLDS['Gulf_SPEI3_wet_increase']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_wet_increase']['max']:.2f}")

print("\nAK Coast SPEI-48 thresholds:")
print(f"  • Dry-side decreases: SPEI {THRESHOLDS['AK_SPEI48_dry_decrease']['min']:.2f}")
print(f"  • Wet-side increases: SPEI {THRESHOLDS['AK_SPEI48_wet_increase']['min']:.2f}")

print("\nAtlantic Coast SPEI-48 thresholds:")
print(f"  • Wet-side decreases: SPEI {THRESHOLDS['Atlantic_SPEI48_wet_decrease']['min']:.2f}")

print("\nPacific Coast SPEI-1 thresholds:")
print(f"  • Dry-side increases: SPEI {THRESHOLDS['Pacific_SPEI1_dry_increase']['min']:.2f} to {THRESHOLDS['Pacific_SPEI1_dry_increase']['max']:.2f}")

# ============================================================================
# SECTION 6: RESULTS - SITE-LEVEL SUMMARY
# ============================================================================

print("\n" + "="*80)
print("SECTION 6: RESULTS - SITE-LEVEL SUMMARY")
print("="*80)

# Load site slopes for additional stats
try:
    site_slopes = pd.read_csv(os.path.join(output_dir, "Q3_WUE_T_SPEI_site_level_slopes.csv"))
    
    print(f"\nSite-level slopes across all sites and timescales:")
    print(f"  • Mean slope: {site_slopes['SPEI_slope'].mean():.4f}")
    print(f"  • Median slope: {site_slopes['SPEI_slope'].median():.4f}")
    print(f"  • % negative slopes: {(site_slopes['SPEI_slope'] < 0).mean()*100:.1f}%")
    
    # Top sensitive sites
    sensitivity = pd.read_csv(os.path.join(output_dir, "Q3_WUE_T_SPEI_site_sensitivity_rank.csv"))
    top5 = sensitivity.nlargest(5, 'mean_abs_pct_change')[['site_name', 'water_class', 'coast_region', 'SPEI_timescale', 'mean_abs_pct_change']]
    
    print(f"\nTop 5 most sensitive sites overall:")
    for _, row in top5.iterrows():
        print(f"  • {row['site_name']} ({row['water_class']}, {row['coast_region']}) - {row['mean_abs_pct_change']:.1f}% at {row['SPEI_timescale']}")
    
except Exception as e:
    print(f"  Note: Could not load site-level data: {e}")

# ============================================================================
# SECTION 7: FIGURE 4 INTERPRETATION (UPDATED)
# ============================================================================

print("\n" + "="*80)
print("SECTION 7: FIGURE 4 INTERPRETATION (UPDATED)")
print("="*80)

print("\nPanel A: Modeled coast-specific upland WUE_T response")
print("  • Shows GAM-predicted WUE_T percent change across SPEI gradients")
print("  • Displays 5%, 10%, and 20% ecological impact thresholds")
print("  • Gulf Coast SPEI-3 shows strongest threshold behavior")
print("  • AK Coast shows thresholds under wet conditions at SPEI-3 and SPEI-48")
print("  • Atlantic/Pacific coasts show flatter responses")
print("  • Direction is mixed - WUE_T can increase OR decrease under dry/wet conditions")

print("\nPanel B: Observed site-level sensitivity by coastline")
print("  • Shows mean absolute WUE_T change from near-normal conditions")
print("  • Groups sites by coastline and SPEI timescale")
print("  • Gulf and AK coasts show highest observed sensitivity")
print("  • Atlantic/Pacific coasts show lower observed sensitivity")

print("\nPanel B provides observed site-level support for the coast-specific threshold patterns shown in Panel A.")
print("The GAMs were fit using all ecosystem classes, but Figure 4A predictions were standardized to")
print("upland conditions, July, and population-level site effects to isolate coast-specific SPEI threshold behavior.")

print("\nKEY INSIGHT (UPDATED):")
print("  The strongest modeled nonlinear response was Gulf Coast SPEI-3, while observed site-level")
print("  sensitivity was highest in the Gulf and Alaska coasts. This means SPEI thresholds are most")
print("  informative when interpreted by coastline and drought-memory timescale, not as one universal")
print("  SPEI cutoff. Coast-specific thresholds should be used explicitly in drought-impact workflows,")
print("  with threshold direction and uncertainty treated explicitly.")

# ============================================================================
# SECTION 8: MANUSCRIPT-READY VALUE SUMMARY
# ============================================================================

print("\n" + "="*80)
print("SECTION 8: MANUSCRIPT-READY VALUE SUMMARY")
print("="*80)

print("\nMETHODS:")
print(f"  • Dataset: {N_OBSERVATIONS:,} observations, {N_SITES} sites, {YEAR_START}-{YEAR_END} ({N_YEARS} years)")
print(f"  • Long format: {N_LONG_OBS:,} site-month-timescale observations")
print(f"  • {N_SPEI_TIMESCALES} SPEI timescales (1-48 months)")
print(f"  • 3 ecosystem classes: Upland ({ECOSYSTEM_DIST['Upland']['obs']} obs, {ECOSYSTEM_DIST['Upland']['sites']} sites), Freshwater ({ECOSYSTEM_DIST['Freshwater']['obs']} obs, {ECOSYSTEM_DIST['Freshwater']['sites']} sites), Saline ({ECOSYSTEM_DIST['Saline']['obs']} obs, {ECOSYSTEM_DIST['Saline']['sites']} sites)")
print(f"  • 4 coast regions: Atlantic ({COAST_DIST['Atlantic Coast']['obs']} obs, {COAST_DIST['Atlantic Coast']['sites']} sites), Pacific ({COAST_DIST['Pacific Coast']['obs']} obs, {COAST_DIST['Pacific Coast']['sites']} sites), Gulf ({COAST_DIST['Gulf Coast']['obs']} obs, {COAST_DIST['Gulf Coast']['sites']} sites), AK ({COAST_DIST['AK Coast']['obs']} obs, {COAST_DIST['AK Coast']['sites']} sites)")

print("\nMODEL PERFORMANCE:")
print(f"  • Coast-threshold GAM: {MODEL_COMP['coast']['deviance_explained']*100:.1f}% deviance explained, R²={MODEL_COMP['coast']['r_squared']:.3f}")
print(f"  • Ecosystem GAM: {MODEL_COMP['ecosystem']['deviance_explained']*100:.1f}% deviance explained, R²={MODEL_COMP['ecosystem']['r_squared']:.3f}")
print(f"  • ΔAIC: {MODEL_COMP['ecosystem']['AIC'] - MODEL_COMP['coast']['AIC']:.2f} (coast-threshold model better)")

print("\nSTRONGEST SIGNAL:")
print(f"  • {STRONGEST_SIGNAL['coast']} {STRONGEST_SIGNAL['timescale']}: EDF={STRONGEST_SIGNAL['edf']:.2f}, F={STRONGEST_SIGNAL['f_stat']:.2f}, p<0.001")

print("\nTHRESHOLDS:")
print(f"  • Gulf Coast SPEI-3 dry-side decreases: {THRESHOLDS['Gulf_SPEI3_dry_decrease']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_dry_decrease']['max']:.2f}")
print(f"  • Gulf Coast SPEI-3 wet-side increases: {THRESHOLDS['Gulf_SPEI3_wet_increase']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_wet_increase']['max']:.2f}")

print("\nSENSITIVITY:")
for coast, values in COAST_SENSITIVITY.items():
    print(f"  • {coast}: {values['SPEI_1']:.1f}% (SPEI-1), {values['SPEI_3']:.1f}% (SPEI-3), {values['SPEI_48']:.1f}% (SPEI-48)")

print("\nOVERALL INTERPRETATION:")
print("  • SPEI thresholds are most informative when interpreted by coastline and drought-memory timescale")
print("  • The strongest modeled response was Gulf Coast SPEI-3 (EDF=4.34, F=29.44, p<0.001)")
print("  • Observed sensitivity was highest in Gulf and Alaska coasts (∼40-45%)")
print("  • Atlantic and Pacific coasts showed lower sensitivity (∼16-20%)")
print("  • Coast-specific thresholds should be used in drought-impact workflows")
print("  • Threshold direction and uncertainty should be treated explicitly")

# ============================================================================
# SECTION 9: SAVE COMPLETE REPORT
# ============================================================================

print("\n" + "="*80)
print("SECTION 9: SAVING COMPLETE REPORT")
print("="*80)

report_file = os.path.join(output_dir, "Q3_COMPLETE_UPDATED_WORKFLOW_report.txt")
print(f"\nSaving report to: {report_file}")

try:
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("Q3 COMPLETE UPDATED WORKFLOW - METHODS & RESULTS EXTRACTION\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n\n")
        
        # Methods section
        f.write("METHODS:\n")
        f.write("-"*80 + "\n")
        f.write(f"Dataset: {N_OBSERVATIONS:,} observations, {N_SITES} sites, {YEAR_START}-{YEAR_END} ({N_YEARS} years)\n")
        f.write(f"Long format: {N_LONG_OBS:,} site-month-timescale observations\n")
        f.write(f"SPEI timescales: {N_SPEI_TIMESCALES} (1-48 months)\n\n")
        
        f.write("Ecosystem distribution:\n")
        for eco, stats in ECOSYSTEM_DIST.items():
            f.write(f"  {eco}: {stats['obs']} obs, {stats['sites']} sites\n")
        
        f.write("\nCoast distribution:\n")
        for coast, stats in COAST_DIST.items():
            f.write(f"  {coast}: {stats['obs']} obs, {stats['sites']} sites\n")
        
        f.write("\nModel Comparison (Table Q3.1):\n")
        f.write(f"  Ecosystem GAM: {MODEL_COMP['ecosystem']['deviance_explained']*100:.1f}% deviance, R²={MODEL_COMP['ecosystem']['r_squared']:.3f}, AIC={MODEL_COMP['ecosystem']['AIC']:.2f}\n")
        f.write(f"  Coast-threshold GAM: {MODEL_COMP['coast']['deviance_explained']*100:.1f}% deviance, R²={MODEL_COMP['coast']['r_squared']:.3f}, AIC={MODEL_COMP['coast']['AIC']:.2f}\n")
        f.write(f"  ΔAIC: {MODEL_COMP['ecosystem']['AIC'] - MODEL_COMP['coast']['AIC']:.2f}\n\n")
        
        # Results section
        f.write("RESULTS:\n")
        f.write("-"*80 + "\n")
        f.write(f"Strongest signal: {STRONGEST_SIGNAL['coast']} {STRONGEST_SIGNAL['timescale']} (EDF={STRONGEST_SIGNAL['edf']:.2f}, F={STRONGEST_SIGNAL['f_stat']:.2f}, p<0.001)\n\n")
        
        f.write("Coastline sensitivity medians (Table Q3.3):\n")
        for coast, values in COAST_SENSITIVITY.items():
            f.write(f"  {coast}: {values['SPEI_1']:.1f}% (SPEI-1), {values['SPEI_3']:.1f}% (SPEI-3), {values['SPEI_48']:.1f}% (SPEI-48)\n")
        
        f.write("\nThreshold summaries:\n")
        f.write(f"  Gulf Coast SPEI-3 dry-side decreases: {THRESHOLDS['Gulf_SPEI3_dry_decrease']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_dry_decrease']['max']:.2f}\n")
        f.write(f"  Gulf Coast SPEI-3 wet-side increases: {THRESHOLDS['Gulf_SPEI3_wet_increase']['min']:.2f} to {THRESHOLDS['Gulf_SPEI3_wet_increase']['max']:.2f}\n")
        f.write(f"  AK Coast SPEI-48 dry-side decreases: {THRESHOLDS['AK_SPEI48_dry_decrease']['min']:.2f}\n")
        f.write(f"  AK Coast SPEI-48 wet-side increases: {THRESHOLDS['AK_SPEI48_wet_increase']['min']:.2f}\n\n")
        
        # Interpretation
        f.write("INTERPRETATION:\n")
        f.write("-"*80 + "\n")
        f.write("Panel B provides observed site-level support for the coast-specific threshold patterns shown in Panel A.\n")
        f.write("The GAMs were fit using all ecosystem classes, but Figure 4A predictions were standardized to\n")
        f.write("upland conditions, July, and population-level site effects to isolate coast-specific SPEI threshold behavior.\n\n")
        
        f.write("KEY INSIGHT:\n")
        f.write("  The strongest modeled nonlinear response was Gulf Coast SPEI-3, while observed site-level\n")
        f.write("  sensitivity was highest in the Gulf and Alaska coasts. This means SPEI thresholds are most\n")
        f.write("  informative when interpreted by coastline and drought-memory timescale, not as one universal\n")
        f.write("  SPEI cutoff. Coast-specific thresholds should be used explicitly in drought-impact workflows,\n")
        f.write("  with threshold direction and uncertainty treated explicitly.\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"  ✓ Complete report saved successfully")
except Exception as e:
    print(f"  ERROR saving report: {e}")

print("\n" + "="*80)
print("Q3 COMPLETE UPDATED WORKFLOW COMPLETE")
print("="*80)
print(f"\nReport saved to: {report_file}")
print("\nKey findings confirmed:")
print(f"  • Dataset: {N_OBSERVATIONS:,} observations, {N_SITES} sites, {YEAR_START}-{YEAR_END}")
print(f"  • Best model: Coast-threshold GAM ({MODEL_COMP['coast']['deviance_explained']*100:.1f}% deviance, R²={MODEL_COMP['coast']['r_squared']:.3f})")
print(f"  • Strongest signal: {STRONGEST_SIGNAL['coast']} {STRONGEST_SIGNAL['timescale']} (EDF={STRONGEST_SIGNAL['edf']:.2f}, F={STRONGEST_SIGNAL['f_stat']:.2f})")
print(f"  • Highest sensitivity: Alaska Coast ({COAST_SENSITIVITY['Alaska']['SPEI_48']:.1f}% at SPEI-48)")
print(f"  • Key insight: Coast-specific SPEI thresholds are most informative for drought-impact workflows")
print("="*80)