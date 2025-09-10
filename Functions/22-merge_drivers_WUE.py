# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 10:06:34 2025

@author: ammar
"""
#################################################################################################################

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 20:56:21 2025

@author: ammar
"""

#et-surface evaporation (i.e., evaporation from open water + saturated soil + intercepted water on leaves) across the whole flux footprint.

## calculate evaporation based on penmontheith model

import os
import numpy as np
import pandas as pd

# ---------- helper: add evaporation (PM, rs=0) and transpiration columns ----------
def add_evap_and_trans_pm(
    df,
    # column names in your files
    DateTime_col='DateTime',  # only to infer timestep if not given
    Rn_col='NETRAD_f',        # W m-2 (net radiation above surface/canopy)
    Ta_col='Tair_f',          # °C
    P_col='PA_f',             # kPa  (your data)
    VPD_col='VPD_f',          # hPa  (your data)
    RH_col='RH_f',            # %    (fallback if VPD missing)
    WS_col='WS',              # m s-1 (at height z)
    # ET input: if you already have ET (mm per step), give its column; else compute from LE_f
    ET_mm_col=None,           # e.g., 'ET'; if None, ET is computed from LE_f
    LE_col='LE_f',            # W m-2 (used only when ET_mm_col is None)
    # site/physics params
    z=3.0,                    # wind measurement height (m)
    z0m_water=2e-4,           # open-water momentum roughness (m)
    kB_inv=4.0,               # scalar roughness param → z0h = z0m*exp(-kB_inv)
    storage_frac=0.2,         # sub-daily storage: S = storage_frac * Rn_surf  (set 0 for daily)
    Rn_scale_to_surface=1.0,  # scale NETRAD_f to surface if needed (leave 1.0)
    timestep_seconds=None,    # infer from DateTime if None
    clip_nonnegative=True
):
    """
    Returns df with two added columns:
      Evap_pen  : evaporation (mm per timestep) using Penman-Monteith, capped so Evap_pen ≤ ET_mm
      Trans_pen : transpiration by difference = ET_mm - Evap_pen (≥ 0)

    Uses Penman–Monteith with r_s = 0 (open-water / wet-surface).
    Assumes: VPD_f in hPa, PA_f in kPa (converted internally).
    """
    out = df.copy()

    # --- timestep (s) ---
    if timestep_seconds is None:
        try:
            dt = pd.to_datetime(out[DateTime_col]).sort_values().diff().dt.total_seconds().median()
            timestep_seconds = 1800.0 if (dt is None or np.isnan(dt)) else float(dt)
        except Exception:
            timestep_seconds = 1800.0

    # --- inputs & unit conversions ---
    Ta = out[Ta_col].to_numpy(dtype=float)           # °C
    T_k = Ta + 273.15                                # K
    P_pa = out[P_col].to_numpy(dtype=float) * 1000.0 # kPa -> Pa

    # --- VPD calculation ---
    RH = out[RH_col].to_numpy(dtype=float) if RH_col in out.columns else np.full(len(out), np.nan)
    RH = np.clip(RH, 0, 100)
    
    es_pa = 610.78 * np.exp(17.2694 * Ta / (Ta + 237.3))
    ea_pa = es_pa * (RH / 100.0)
    
    if VPD_col in out.columns:
        vpd_hpa = out[VPD_col].to_numpy(dtype=float)
        VPD_pa = np.where(np.isfinite(vpd_hpa), np.clip(vpd_hpa, 0, None) * 100.0,
                          np.maximum(0.0, es_pa - ea_pa))
    else:
        VPD_pa = np.maximum(0.0, es_pa - ea_pa)

    # Δ (Pa/K), λ (J/kg), γ (Pa/K)
    Delta = 4098.0 * es_pa / (Ta + 237.3)**2  # Pa/K

    # CORRECTED: Proper latent heat calculation
    lambda_v = (2.501 - 0.002361 * Ta) * 1e6  # J/kg (MJ/kg to J/kg conversion)
    
    cp = 1004.0  # J/kg/K, specific heat of air at constant pressure
    
    # CORRECTED: Psychrometric constant
    gamma = (cp * P_pa) / (0.622 * lambda_v)  # Pa/K
    
    # --- CORRECTED: Air density with VPD-consistent vapor pressure ---
    e = np.maximum(0.0, es_pa - VPD_pa)  # Pa (use VPD-consistent vapor pressure)
    R_d = 287.058  # J/kg/K, gas constant for dry air
    rho_a = (P_pa - 0.378 * e) / (R_d * T_k)  # kg/m³

    # --- aerodynamic resistance over open water (neutral) ---
    k = 0.4  # von Karman constant
    d = 0.0  # zero displacement for open water
    z0m = float(z0m_water)
    z0h = z0m * np.exp(-float(kB_inv))

    WS = np.clip(out[WS_col].to_numpy(dtype=float), 0.1, None)
    num_m = np.maximum((z - d) / z0m, 1.01)
    num_h = np.maximum((z - d) / z0h, 1.01)
    
    # Aerodynamic resistance formula
    ra = (np.log(num_m) * np.log(num_h)) / (k**2 * WS)  # s/m

    # Apply reasonable bounds for coastal wetlands
    ra = np.clip(ra, 5.0, 500.0)  # More realistic bounds for wetlands

    # --- surface net radiation & storage ---
    Rn_surf = Rn_scale_to_surface * out[Rn_col].to_numpy(dtype=float)  # W m-2
    S = storage_frac * Rn_surf  # W m-2

    # --- CORRECTED: Penman–Monteith evaporation (r_s=0) with proper aerodynamic term ---
    num = Delta * (Rn_surf - S) + (rho_a * cp / gamma) * (VPD_pa / ra)
    den = Delta + gamma
    lambdaE = num / den  # W m-2

    # --- CORRECTED: Unit conversion from m to mm ---
    # Evap per step (mm): (λE * dt) / (λ * ρw) * 1000
    rho_w = 1000.0  # kg/m³, density of water
    Evap_pen_raw = (lambdaE * timestep_seconds) / (lambda_v * rho_w) * 1000.0  # mm

    # --- ET for capping ---
    if ET_mm_col is not None and ET_mm_col in out.columns:
        ET_mm = out[ET_mm_col].to_numpy(dtype=float)
    else:
        if LE_col not in out.columns:
            raise ValueError("Provide ET_mm_col or LE_f to compute ET from energy.")
        LE = out[LE_col].to_numpy(dtype=float)  # W m-2
        # CORRECTED: Unit conversion from m to mm
        ET_mm = (LE * timestep_seconds) / (lambda_v * rho_w) * 1000.0  # mm

    # --- hygiene, cap Evap ≤ ET, compute Trans ---
    if clip_nonnegative:
        lambdaE = np.where(lambdaE < 0, 0.0, lambdaE)
        Evap_pen_raw = np.where(Evap_pen_raw < 0, 0.0, Evap_pen_raw)
        ET_mm = np.where(ET_mm < 0, 0.0, ET_mm)

    Evap_pen = np.where(np.isfinite(ET_mm), np.minimum(Evap_pen_raw, ET_mm), Evap_pen_raw)
    Trans_pen = np.where(np.isfinite(ET_mm), np.maximum(ET_mm - Evap_pen, 0.0), np.nan)

    # Add only the two new columns to the original dataframe
    out = out.copy()  # Ensure we don't modify the original
    out['Evap_pen'] = Evap_pen
    out['Trans_pen'] = Trans_pen
    
    return out


# ---------- batch processor: run on all CSVs and save with same filename ----------
def process_ameri_files(input_folder, output_folder,
                        et_mm_column=None,   # name if ET already exists in mm/step; else None
                        wind_height_m=3.0,
                        storage_frac=0.2,
                        rn_scale_to_surface=1.0):
    """
    Processes all *.csv in input_folder, adds Evap_pen and Trans_pen, and
    saves each to output_folder with the SAME filename as input.

    Parameters
    ----------
    input_folder : str
        Folder containing input CSV files.
    output_folder : str
        Destination folder for processed CSVs.
    et_mm_column : str or None
        If your ET (mm per step) exists, give its column name (e.g., 'ET').
        If None, ET is computed internally from LE_f and Tair_f.
    wind_height_m : float
        Anemometer height (m) used for aerodynamic resistance.
    storage_frac : float
        Sub-daily water heat storage fraction; set 0 for daily data.
    rn_scale_to_surface : float
        Scale factor to map NETRAD_f to surface net radiation.
    """
    os.makedirs(output_folder, exist_ok=True)

    # guard against raw strings ending with a backslash (Windows)
    if input_folder.endswith("\\"):
        input_folder = input_folder[:-1]
    if output_folder.endswith("\\"):
        output_folder = output_folder[:-1]

    files = [f for f in os.listdir(input_folder) if f.lower().endswith(".csv")]
    if not files:
        print(f"[info] No CSV files found in: {input_folder}")
        return

    print(f"[info] Found {len(files)} CSV file(s). Processing...")

    # columns we rely on (some are alternative paths depending on ET availability)
    required = {'DateTime', 'NETRAD_f', 'Tair_f', 'PA_f', 'WS'}
    # At least one of these humidity options must be present
    humidity_ok = lambda cols: ('VPD_f' in cols) or ('RH_f' in cols)

    for fname in files:
        in_path = os.path.join(input_folder, fname)
        try:
            df = pd.read_csv(in_path)
            cols = set(df.columns)

            # quick checks
            missing = required - cols
            if missing or not humidity_ok(cols):
                print(f"[warn] Skipping {fname}: missing {missing or ''} "
                      f"{'(need VPD_f or RH_f)' if not humidity_ok(cols) else ''}")
                continue

            # run evaporation + transpiration add-on
            df_out = add_evap_and_trans_pm(
                df,
                DateTime_col='DateTime',
                Rn_col='NETRAD_f',
                Ta_col='Tair_f',
                P_col='PA_f',
                VPD_col='VPD_f',      # expected in hPa
                RH_col='RH_f',        # fallback if VPD_f missing
                WS_col='WS',
                ET_mm_col=et_mm_column,   # pass your ET column name if available
                LE_col='LE_f',
                z=wind_height_m,
                z0m_water=2e-4,
                kB_inv=4.0,
                storage_frac=storage_frac,
                Rn_scale_to_surface=rn_scale_to_surface,
                timestep_seconds=None,     # infer from DateTime spacing
            )

            out_path = os.path.join(output_folder, fname)
            df_out.to_csv(out_path, index=False)
            print(f"[ok] Saved: {out_path}")

        except Exception as e:
            print(f"[error] {fname}: {e}")

    print("[done] All eligible files processed.")


# ---------- example usage ----------
if __name__ == "__main__":
    # UNC paths: avoid trailing backslash in Python strings.
    input_folder  = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated"
    output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split"

    # If you already have ET (mm per step) in a column, put its name here, e.g., 'ET'.
    # Otherwise leave as None and the code computes ET from LE_f.
    process_ameri_files(input_folder, output_folder, et_mm_column=None,
                        wind_height_m=3.0, storage_frac=0.2, rn_scale_to_surface=1.0)

#################################################################################################################


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def analyze_et_components(input_folder, output_folder_name="3-et_component_analysis"):
    """
    Analyze evaporation and transpiration components on annual basis.
    Focuses only on Penman-Monteith components (Evap_pen/Trans_pen).
    
    Parameters
    ----------
    input_folder : str
        Folder containing CSV files with ET components
    output_folder_name : str
        Name of subfolder to create for analysis results
    """
    
    # Create output folder within the input directory structure
    input_path = Path(input_folder)
    output_path = input_path.parent / output_folder_name
    os.makedirs(output_path, exist_ok=True)
    
    # Get all CSV files in input folder
    csv_files = list(input_path.glob("*.csv"))
    
    if not csv_files:
        print(f"No CSV files found in: {input_folder}")
        return
    
    print(f"Found {len(csv_files)} CSV file(s). Analyzing...")
    print(f"Output will be saved to: {output_path}")
    
    # Initialize results storage
    all_results = []
    site_summaries = []
    
    for csv_file in csv_files:
        try:
            site_name = csv_file.stem
            print(f"Analyzing {site_name}...")
            
            # Read CSV file
            df = pd.read_csv(csv_file)
            
            # Convert DateTime to datetime if needed
            if 'DateTime' in df.columns:
                df['DateTime'] = pd.to_datetime(df['DateTime'])
            
            # Ensure required columns exist
            required_cols = ['Year', 'ET', 'Evap_pen', 'Trans_pen']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                print(f"  Skipping {site_name}: missing columns {missing_cols}")
                continue
            
            # Filter out rows with missing ET data
            df_clean = df.dropna(subset=['ET', 'Evap_pen', 'Trans_pen'])
            
            if df_clean.empty:
                print(f"  Skipping {site_name}: no valid ET data")
                continue
            
            # Group by year and calculate annual totals
            annual_totals = df_clean.groupby('Year').agg({
                'ET': 'sum',
                'Evap_pen': 'sum',
                'Trans_pen': 'sum'
            }).reset_index()
            
            # Calculate fractions
            annual_totals['Evap_pen_frac'] = annual_totals['Evap_pen'] / annual_totals['ET']
            annual_totals['Trans_pen_frac'] = annual_totals['Trans_pen'] / annual_totals['ET']
            
            # Calculate closure error (should be close to 1 if components sum to ET)
            annual_totals['Pen_closure'] = (annual_totals['Evap_pen'] + annual_totals['Trans_pen']) / annual_totals['ET']
            
            # Add site name
            annual_totals['Site'] = site_name
            
            # Store results
            all_results.append(annual_totals)
            
            # Calculate site-wide summary statistics
            site_summary = {
                'Site': site_name,
                'Years': len(annual_totals),
                'Mean_ET_mm': annual_totals['ET'].mean(),
                'Mean_Evap_pen_frac': annual_totals['Evap_pen_frac'].mean(),
                'Mean_Trans_pen_frac': annual_totals['Trans_pen_frac'].mean(),
                'Mean_Pen_closure': annual_totals['Pen_closure'].mean(),
                'Pen_correlation': annual_totals['Evap_pen_frac'].corr(annual_totals['Trans_pen_frac'])
            }
            site_summaries.append(site_summary)
            
            print(f"  Processed {len(annual_totals)} years of data")
                
        except Exception as e:
            print(f"Error processing {csv_file.name}: {str(e)}")
            continue
    
    if not all_results:
        print("No valid data found in any files.")
        return
    
    # Combine all results
    combined_results = pd.concat(all_results, ignore_index=True)
    summary_df = pd.DataFrame(site_summaries)
    
    # Generate summary statistics
    print("\n" + "="*60)
    print("PENMAN-MONTEITH COMPONENT ANALYSIS")
    print("="*60)
    
    # Method comparison
    print("\nMean Fractions of ET:")
    print(f"Evap_pen/ET: {summary_df['Mean_Evap_pen_frac'].mean():.3f} ± {summary_df['Mean_Evap_pen_frac'].std():.3f}")
    print(f"Trans_pen/ET: {summary_df['Mean_Trans_pen_frac'].mean():.3f} ± {summary_df['Mean_Trans_pen_frac'].std():.3f}")
    
    print(f"\nClosure Error (should be close to 1.0):")
    print(f"Evap_pen + Trans_pen: {summary_df['Mean_Pen_closure'].mean():.3f} ± {summary_df['Mean_Pen_closure'].std():.3f}")
    
    print(f"\nCorrelation between Evap_pen and Trans_pen fractions:")
    print(f"Correlation coefficient: {summary_df['Pen_correlation'].mean():.3f}")
    
    # Create visualizations
    create_visualizations(combined_results, summary_df, output_path)
    
    return combined_results, summary_df


def create_visualizations(combined_results, summary_df, output_path):
    """Create visualizations of the ET component analysis focusing only on Penman-Monteith."""
    
    # Set plot style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create a 2x2 panel plot focusing only on Penman-Monteith components
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Penman-Monteith ET Component Analysis', fontsize=16, fontweight='bold')
    
    # 1. Box plot of Penman-Monteith fractions
    pen_data = [combined_results['Evap_pen_frac'], combined_results['Trans_pen_frac']]
    ax1.boxplot(pen_data, labels=['Evaporation\n(Evap_pen)', 'Transpiration\n(Trans_pen)'])
    ax1.set_ylabel('Fraction of ET')
    ax1.set_title('Penman-Monteith ET Components')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)
    
    # 2. Closure error distribution
    ax2.hist(combined_results['Pen_closure'], bins=20, alpha=0.7, edgecolor='black')
    ax2.axvline(x=1.0, color='r', linestyle='--', linewidth=2, label='Perfect closure')
    ax2.set_xlabel('Closure (Evap_pen + Trans_pen) / ET')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Energy Closure Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Scatter plot: Evap_pen vs Trans_pen fractions
    scatter = ax3.scatter(combined_results['Evap_pen_frac'], combined_results['Trans_pen_frac'], 
                         alpha=0.6, c=combined_results['ET'], cmap='viridis', s=50)
    ax3.set_xlabel('Evaporation Fraction (Evap_pen/ET)')
    ax3.set_ylabel('Transpiration Fraction (Trans_pen/ET)')
    ax3.set_title('Evaporation vs Transpiration Fractions\n(color = Total ET)')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    
    # Add colorbar for total ET
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Total ET (mm)')
    
    # 4. Site-wise comparison of mean fractions
    sites = summary_df['Site']
    x_pos = np.arange(len(sites))
    width = 0.35
    
    ax4.bar(x_pos - width/2, summary_df['Mean_Evap_pen_frac'], width, 
            label='Evaporation', alpha=0.8)
    ax4.bar(x_pos + width/2, summary_df['Mean_Trans_pen_frac'], width, 
            label='Transpiration', alpha=0.8)
    
    ax4.set_xlabel('Site')
    ax4.set_ylabel('Mean Fraction of ET')
    ax4.set_title('Site-wise Mean Fractions')
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(sites, rotation=45, ha='right')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_path = output_path / "penman_monteith_analysis.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Penman-Monteith analysis plot: {plot_path}")


# Example usage
if __name__ == "__main__":
    input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split"
    
    # Run the analysis - output will be created in the parent directory
    results, summary = analyze_et_components(input_folder)
    
    if results is not None:
        print("\nAnalysis completed successfully!")
        print(f"Processed {len(summary)} sites with {len(results)} annual records")



#################################################################################################################


import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns

def filter_growing_season_data(input_folder, output_folder, growing_season_folder):
    """
    Filter data to include only growing season periods based on phenofit growing season data
    
    Parameters:
    input_folder: folder containing input CSV files
    output_folder: folder to save filtered CSV files
    growing_season_folder: folder containing phenofit_growing_season.csv
    """
    
    # Read growing season data
    growing_season_file = os.path.join(growing_season_folder, 'phenofit_growing_season.csv')
    growing_season_df = pd.read_csv(growing_season_file)
    
    # Get list of input files
    input_files = glob.glob(os.path.join(input_folder, '*.csv'))
    
    # Set up plotting style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    for input_file in input_files:
        # Extract site name from filename
        site_name = os.path.basename(input_file).split('.')[0]
        
        # Get growing season data for this site
        site_gs_data = growing_season_df[growing_season_df['site_name'] == site_name]
        
        if site_gs_data.empty:
            print(f"No growing season data found for site: {site_name}")
            continue
        
        # Calculate average sos and eos for this site (across all years)
        avg_sos_site = site_gs_data['avg_sos'].mean()
        avg_eos_site = site_gs_data['avg_eos'].mean()
        avg_length_site = site_gs_data['avg_length'].mean()
        
        print(f"\nSite: {site_name}")
        print(f"Using average growing season: DoY {avg_sos_site:.1f} to {avg_eos_site:.1f}")
        print(f"Average length: {avg_length_site:.1f} days")
        
        # Read the input data
        df = pd.read_csv(input_file)
        
        # Convert DateTime to datetime if it's not already
        if 'DateTime' in df.columns:
            df['DateTime'] = pd.to_datetime(df['DateTime'], errors='coerce')
        
        # Add year column if not present
        if 'Year' not in df.columns and 'DateTime' in df.columns:
            df['Year'] = df['DateTime'].dt.year
        
        # Remove rows with invalid years
        df = df.dropna(subset=['Year'])
        df['Year'] = df['Year'].astype(int)
        
        # Filter data for growing season using site-average values
        growing_season_mask = (df['DoY'] >= avg_sos_site) & (df['DoY'] <= avg_eos_site)
        gs_df = df[growing_season_mask].copy()
        
        # Add growing season metadata
        gs_df['avg_sos'] = avg_sos_site
        gs_df['avg_eos'] = avg_eos_site
        gs_df['avg_length'] = avg_length_site
        
        if gs_df.empty:
            print(f"No growing season data found for site: {site_name}")
            continue
            
        # Now apply the 50% data availability filter on growing season data
        valid_years = []
        yearly_gs_data = []
        
        for year, year_gs_data in gs_df.groupby('Year'):
            # Check data availability for NEE_f and ET in growing season
            nee_available = year_gs_data['NEE_f'].notna() if 'NEE_f' in year_gs_data.columns else pd.Series([False] * len(year_gs_data))
            et_available = year_gs_data['ET'].notna() if 'ET' in year_gs_data.columns else pd.Series([False] * len(year_gs_data))
            
            nee_available_pct = nee_available.mean() * 100
            et_available_pct = et_available.mean() * 100
            
            # Keep year only if both columns have >=50% available data in growing season
            if nee_available_pct >= 50 and et_available_pct >= 50:
                valid_years.append(year)
                yearly_gs_data.append(year_gs_data)
                print(f"Keeping year {year}: "
                      f"NEE_f available={nee_available_pct:.1f}%, ET available={et_available_pct:.1f}%")
            else:
                print(f"Excluding year {year}: "
                      f"NEE_f available={nee_available_pct:.1f}%, ET available={et_available_pct:.1f}%")
        
        if not yearly_gs_data:
            print(f"No valid years with sufficient data in growing season for site: {site_name}")
            continue
            
        # Combine valid growing season data
        filtered_df = pd.concat(yearly_gs_data, ignore_index=True)
        
        # Print ranges for specified columns
        print(f"\nRanges for key columns in growing season:")
        for col in ['ET', 'Evap_pen', 'Trans_pen']:
            if col in filtered_df.columns:
                col_data = filtered_df[col].dropna()
                if len(col_data) > 0:
                    print(f"  {col}: {col_data.min():.4f} - {col_data.max():.4f} "
                          f"(n={len(col_data)} observations)")
                else:
                    print(f"  {col}: No valid data")
            else:
                print(f"  {col}: Column not found")
        
        # Calculate yearly sums and ratios - only using rows where all three columns are present
        yearly_data = []
        yearly_ratios = []
        
        for year, year_data in filtered_df.groupby('Year'):
            # Create a mask for rows where all three columns are present
            valid_mask = (year_data['ET'].notna() & 
                         year_data['Evap_pen'].notna() & 
                         year_data['Trans_pen'].notna())
            
            valid_data = year_data[valid_mask]
            
            if len(valid_data) > 0:
                year_sum = valid_data[['ET', 'Evap_pen', 'Trans_pen']].sum()
                
                # Calculate ratios
                if year_sum['ET'] > 0:
                    evap_ratio = year_sum['Evap_pen'] / year_sum['ET']
                    trans_ratio = year_sum['Trans_pen'] / year_sum['ET']
                else:
                    evap_ratio = trans_ratio = np.nan
                
                yearly_data.append({
                    'Year': year,
                    'ET_sum': year_sum['ET'],
                    'Evap_pen_sum': year_sum['Evap_pen'],
                    'Trans_pen_sum': year_sum['Trans_pen'],
                    'Evap_pen/ET': evap_ratio,
                    'Trans_pen/ET': trans_ratio,
                    'n_valid_observations': len(valid_data)
                })
                
                yearly_ratios.append(evap_ratio)
        
        if not yearly_data:
            print(f"No valid yearly data for site: {site_name}")
            continue
            
        yearly_df = pd.DataFrame(yearly_data)
        print("\nYearly sums and ratios (using only rows where all three columns are present):")
        print(yearly_df.to_string(index=False))
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Site: {site_name} - Growing Season Analysis\n(DoY {avg_sos_site:.0f}-{avg_eos_site:.0f})', fontsize=16)
        
        # Plot 1: Time series of key variables
        if 'DateTime' in filtered_df.columns:
            # Plot each year separately for better visualization
            for year, year_data in filtered_df.groupby('Year'):
                axes[0, 0].plot(year_data['DoY'], year_data['ET'], 'b-', alpha=0.7, label=f'ET {year}' if year == filtered_df['Year'].iloc[0] else "")
                axes[0, 0].plot(year_data['DoY'], year_data['Evap_pen'], 'r-', alpha=0.7, label=f'Evap_pen {year}' if year == filtered_df['Year'].iloc[0] else "")
                axes[0, 0].plot(year_data['DoY'], year_data['Trans_pen'], 'g-', alpha=0.7, label=f'Trans_pen {year}' if year == filtered_df['Year'].iloc[0] else "")
            
            axes[0, 0].set_xlabel('Day of Year')
            axes[0, 0].set_ylabel('Flux (mm)')
            axes[0, 0].set_title('Time Series of ET Components (by DoY)')
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Yearly ratios
        years = yearly_df['Year']
        ratios = yearly_df['Evap_pen/ET']
        bars = axes[0, 1].bar(range(len(years)), ratios, alpha=0.7)
        axes[0, 1].set_xlabel('Year')
        axes[0, 1].set_ylabel('Evap_pen/ET Ratio')
        axes[0, 1].set_title('Yearly Evap_pen/ET Ratio')
        axes[0, 1].set_xticks(range(len(years)))
        axes[0, 1].set_xticklabels([str(int(y)) for y in years])
        axes[0, 1].grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, (bar, ratio) in enumerate(zip(bars, ratios)):
            axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                           f'{ratio:.2f}', ha='center', va='bottom')
        
        # Plot 3: Distribution of ET components
        et_components = filtered_df[['ET', 'Evap_pen', 'Trans_pen']].dropna()
        boxplot_data = [et_components['ET'], et_components['Evap_pen'], et_components['Trans_pen']]
        box = axes[1, 0].boxplot(boxplot_data, patch_artist=True)
        
        # Add colors to boxplot
        colors = ['lightblue', 'lightcoral', 'lightgreen']
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
        
        axes[1, 0].set_xticklabels(['ET', 'Evap_pen', 'Trans_pen'])
        axes[1, 0].set_ylabel('Flux (mm)')
        axes[1, 0].set_title('Distribution of ET Components')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Scatter plot of Evap_pen vs Trans_pen with year colors
        unique_years = filtered_df['Year'].unique()
        colors = plt.cm.viridis(np.linspace(0, 1, len(unique_years)))
        
        for i, year in enumerate(unique_years):
            year_data = filtered_df[filtered_df['Year'] == year]
            axes[1, 1].scatter(year_data['Evap_pen'], year_data['Trans_pen'], 
                              alpha=0.6, color=colors[i], label=str(year))
        
        axes[1, 1].set_xlabel('Evap_pen (mm)')
        axes[1, 1].set_ylabel('Trans_pen (mm)')
        axes[1, 1].set_title('Evap_pen vs Trans_pen (colored by year)')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(output_folder, f"{site_name}_analysis_plot.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Saved analysis plot to: {plot_file}")
        
        # Save filtered data
        output_file = os.path.join(output_folder, f"{site_name}.csv")
        filtered_df.to_csv(output_file, index=False)
        print(f"Saved filtered data to: {output_file}")
        
        # Print summary statistics
        print(f"\nSummary for {site_name}:")
        print(f"Total years with valid data: {len(yearly_df)}")
        print(f"Average Evap_pen/ET ratio: {yearly_df['Evap_pen/ET'].mean():.3f} ± {yearly_df['Evap_pen/ET'].std():.3f}")
        print(f"Average Trans_pen/ET ratio: {yearly_df['Trans_pen/ET'].mean():.3f} ± {yearly_df['Trans_pen/ET'].std():.3f}")
        print(f"Growing season: DoY {avg_sos_site:.1f} to {avg_eos_site:.1f} ({avg_length_site:.1f} days)")


# Example usage
if __name__ == "__main__":
    input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split"
    output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split\3-growing_season_transp_Evap"
    growing_season_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
    

    
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    filter_growing_season_data(input_folder, output_folder, growing_season_folder)
    
    

##########################################################################################################

## monthly driver data


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import re
from pathlib import Path
import warnings

def plot_numeric_variables_by_site(df, plot_path):
    """
    Create boxplots for all numeric variables across different sites.
    """
    # Identify numeric columns (exclude text/identifier columns)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    
    # Remove identifier columns that happen to be numeric
    exclude_cols = ['year', 'month', 'lat', 'long']
    numeric_cols = [col for col in numeric_cols if col not in exclude_cols]
    
    if not numeric_cols:
        print("No numeric variables found to plot")
        return
    
    # Create subplots - 4 plots per row
    n_cols = 4
    n_rows = (len(numeric_cols) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
    fig.suptitle('Distribution of Numeric Variables Across Sites', fontsize=16, y=0.98)
    
    # Flatten axes array for easy indexing
    if n_rows > 1:
        axes = axes.flatten()
    else:
        axes = [axes] if n_cols == 1 else axes
    
    for i, col in enumerate(numeric_cols):
        if i < len(axes):
            ax = axes[i]
            
            # Create boxplot
            sns.boxplot(data=df, x='site_name', y=col, ax=ax)
            ax.set_title(f'{col} by Site')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.tick_params(axis='x', labelsize=8)
            
            # Remove empty subplots
        else:
            break
    
    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'numeric_variables_by_site.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created boxplot for {len(numeric_cols)} numeric variables")

def monthly_drivers(min_coverage=0.50):
    """
    Aggregate driver CSVs to monthly resolution, add monthly LAI, join with reference data,
    compute WUE_tra and WUE_eva, and create diagnostic boxplots.
    """
    
    # Define paths
    driver_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split\3-growing_season_transp_Evap"
    reference_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
    lai_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai\growing_season_LAI"
    plot_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\plots_drivers"
    output_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_with_drivers.csv"
    
    os.makedirs(plot_path, exist_ok=True)
    log_messages = []
    
    # 1. Load reference data first
    log_messages.append("Loading reference data...")
    try:
        reference_df = pd.read_csv(reference_path)
        target_months = reference_df[['site_name', 'year', 'month']].drop_duplicates()
        log_messages.append(f"Reference data has {len(target_months)} unique site-month combinations")
        log_messages.append(f"Reference sites: {reference_df['site_name'].unique()}")
    except Exception as e:
        raise ValueError(f"Error loading reference data: {str(e)}")
    
    # 2. Process driver files
    log_messages.append("Processing driver files...")
    driver_files = glob.glob(os.path.join(driver_path, "*.csv"))
    
    if not driver_files:
        raise ValueError(f"No CSV files found in {driver_path}")
    
    all_driver_data = []
    
    for file in driver_files:
        filename = os.path.basename(file)
        
        # Fix site name extraction - remove .csv extension and keep full name
        site_name = filename.replace('.csv', '')
        
        try:
            # Read the file
            df = pd.read_csv(file)
            
            if 'DateTime' not in df.columns:
                log_messages.append(f"Skipping {filename}: No DateTime column")
                continue
            
            # Parse DateTime - try multiple formats
            try:
                # First try the format that matches your data: YYYY-MM-DD HH:MM:SS
                df['datetime'] = pd.to_datetime(df['DateTime'], format='%Y-%m-%d %H:%M:%S')
            except:
                try:
                    # Try without seconds
                    df['datetime'] = pd.to_datetime(df['DateTime'], format='%Y-%m-%d %H:%M')
                except:
                    try:
                        # Try another common format
                        df['datetime'] = pd.to_datetime(df['DateTime'], format='%m/%d/%Y %H:%M')
                    except:
                        # Fall back to pandas automatic parsing
                        df['datetime'] = pd.to_datetime(df['DateTime'])
            
            df['year'] = df['datetime'].dt.year
            df['month'] = df['datetime'].dt.month
            df['site_name'] = site_name
            
            log_messages.append(f"Parsed dates for {filename}: {df['datetime'].min()} to {df['datetime'].max()}")
            
            # Filter to only reference months for this site
            site_targets = target_months[target_months['site_name'] == site_name]
            if not site_targets.empty:
                # Merge to keep only months that exist in reference
                merged = pd.merge(df, site_targets, on=['site_name', 'year', 'month'])
                if not merged.empty:
                    all_driver_data.append(merged)
                    log_messages.append(f"Added {len(merged)} rows from {filename} for site {site_name}")
                else:
                    log_messages.append(f"No matching months between driver {filename} and reference for site {site_name}")
            else:
                log_messages.append(f"No reference months found for site {site_name}")
                
        except Exception as e:
            log_messages.append(f"Error processing {filename}: {str(e)}")
            continue
    
    if not all_driver_data:
        log_messages.append("WARNING: No driver data found for any reference months!")
        driver_monthly = pd.DataFrame(columns=['site_name', 'year', 'month'])
    else:
        # Combine all driver data
        driver_df = pd.concat(all_driver_data, ignore_index=True)
        log_messages.append(f"Total driver rows: {len(driver_df)}")
        log_messages.append(f"Unique site-month combinations in driver data: {len(driver_df[['site_name', 'year', 'month']].drop_duplicates())}")
        
        # Define aggregation rules
        sum_cols = ['precip_mm', 'Evap_pen', 'Trans_pen', 'Co2']
        mean_cols = ['NETRAD_f', 'WS', 'WD', 'H_f', 'PAR_f', 'RH_f', 'PA_f']
        median_cols = ['gs_iPM', 'gs_FG']
        
        # Check which columns exist
        available_sum_cols = [col for col in sum_cols if col in driver_df.columns]
        available_mean_cols = [col for col in mean_cols if col in driver_df.columns]
        available_median_cols = [col for col in median_cols if col in driver_df.columns]
        
        log_messages.append(f"Columns to aggregate - Sum: {available_sum_cols}, Mean: {available_mean_cols}, Median: {available_median_cols}")
        
        # Create aggregation dictionary
        agg_dict = {}
        for col in available_sum_cols:
            agg_dict[col] = 'sum'
        for col in available_mean_cols:
            agg_dict[col] = 'mean'
        for col in available_median_cols:
            agg_dict[col] = 'median'
        
        # Group and aggregate
        grouped = driver_df.groupby(['site_name', 'year', 'month'])
        driver_monthly = grouped.agg(agg_dict).reset_index()
        
        # Apply coverage threshold
        for col in agg_dict.keys():
            for (site, year, month), group in grouped:
                valid_count = group[col].count()
                total_count = len(group)
                
                if total_count > 0 and valid_count / total_count < min_coverage:
                    mask = (driver_monthly['site_name'] == site) & \
                           (driver_monthly['year'] == year) & \
                           (driver_monthly['month'] == month)
                    driver_monthly.loc[mask, col] = np.nan
        
        log_messages.append(f"Aggregated driver data with {len(driver_monthly)} monthly records")
    
    # 3. Process LAI files
    log_messages.append("Processing LAI files...")
    lai_files = glob.glob(os.path.join(lai_path, "*.csv"))
    
    all_lai_data = []
    for file in lai_files:
        # Extract site name without extension
        site_name = os.path.basename(file).split('.')[0]
        
        try:
            df = pd.read_csv(file)
            
            if 'Date' not in df.columns:
                log_messages.append(f"Skipping LAI file {file}: No Date column")
                continue
            
            df['date'] = pd.to_datetime(df['Date'])
            df['year'] = df['date'].dt.year
            df['month'] = df['date'].dt.month
            df['site_name'] = site_name
            
            # Filter to reference months
            site_targets = target_months[target_months['site_name'] == site_name]
            if not site_targets.empty:
                merged = pd.merge(df, site_targets, on=['site_name', 'year', 'month'])
                if not merged.empty:
                    all_lai_data.append(merged)
        except Exception as e:
            log_messages.append(f"Error reading LAI file {file}: {str(e)}")
    
    if all_lai_data:
        lai_df = pd.concat(all_lai_data, ignore_index=True)
        lai_grouped = lai_df.groupby(['site_name', 'year', 'month'])
        lai_monthly = lai_grouped['lai'].mean().reset_index()
        log_messages.append(f"Aggregated LAI data for {len(lai_monthly)} monthly records")
    else:
        log_messages.append("No LAI data found")
        lai_monthly = pd.DataFrame(columns=['site_name', 'year', 'month', 'lai'])
    
    # 4. Merge all data
    log_messages.append("Merging data...")
    
    # Start with reference, left join drivers, then LAI
    merged_df = reference_df.copy()
    
    if not driver_monthly.empty:
        merged_df = pd.merge(merged_df, driver_monthly, on=['site_name', 'year', 'month'], how='left')
    
    merged_df = pd.merge(merged_df, lai_monthly, on=['site_name', 'year', 'month'], how='left')
    
    # 5. Calculate new variables
    if 'Trans_pen' in merged_df.columns and 'GPP' in merged_df.columns:
        merged_df['WUE_tra'] = np.where(
            (merged_df['Trans_pen'] > 0) & merged_df['Trans_pen'].notna() & merged_df['GPP'].notna(),
            merged_df['GPP'] / merged_df['Trans_pen'],
            np.nan
        )
    
    if 'Evap_pen' in merged_df.columns and 'GPP' in merged_df.columns:
        merged_df['WUE_eva'] = np.where(
            (merged_df['Evap_pen'] > 0) & merged_df['Evap_pen'].notna() & merged_df['GPP'].notna(),
            merged_df['GPP'] / merged_df['Evap_pen'],
            np.nan
        )
    
    # 6. Final column order
    ref_cols = ['site_name', 'year', 'month', 'NEE', 'GPP', 'Reco', 'ET', 'NEP', 'State', 'biome', 
                'salinity_ppt', 'salinity_fine', 'salini_coarse', 'WUE', 'CUE', 'Tair_f', 'VPD_f', 
                'Rg_f', 'avg_length', 'avg_sos', 'avg_eos', 'lat', 'long', 'climate', 'climate_2', 'salinity_con']
    
    driver_cols = ['NETRAD_f', 'WS', 'WD', 'Co2', 'H_f', 'PAR_f', 'RH_f', 'PA_f', 'precip_mm', 
                  'gs_iPM', 'gs_FG', 'Evap_pen', 'Trans_pen']
    
    extra_cols = ['lai', 'WUE_tra', 'WUE_eva']
    
    # Reorder columns
    final_cols = []
    for col_list in [ref_cols, driver_cols, extra_cols]:
        for col in col_list:
            if col in merged_df.columns and col not in final_cols:
                final_cols.append(col)
    
    # Add any remaining columns
    for col in merged_df.columns:
        if col not in final_cols and col not in ['NEE_f', 'LE_f']:
            final_cols.append(col)
    
    merged_df = merged_df[final_cols]
    
    # 7. Create diagnostic plots
    log_messages.append("Creating diagnostic plots...")
    try:
        plot_numeric_variables_by_site(merged_df, plot_path)
        log_messages.append("Successfully created numeric variables boxplot")
    except Exception as e:
        log_messages.append(f"Error creating plots: {str(e)}")
    
    # 8. Save and log results
    log_messages.append("Saving final dataset...")
    merged_df.to_csv(output_path, index=False)
    
    # Data availability summary
    driver_vars = ['NETRAD_f', 'WS', 'WD', 'Co2', 'H_f', 'PAR_f', 'RH_f', 'PA_f', 'precip_mm', 
                  'gs_iPM', 'gs_FG', 'Evap_pen', 'Trans_pen']
    extra_vars = ['lai', 'WUE_tra', 'WUE_eva']
    
    log_messages.append("\nData availability summary:")
    for col in driver_vars + extra_vars:
        if col in merged_df.columns:
            non_na = merged_df[col].notna().sum()
            total = len(merged_df)
            log_messages.append(f"{col}: {non_na}/{total} ({non_na/total:.1%})")
        else:
            log_messages.append(f"{col}: Column not found")
    
    # Print all logs
    for msg in log_messages:
        print(msg)
    
    return merged_df

# Execute
if __name__ == "__main__":
    monthly_data = monthly_drivers(min_coverage=0.50)



##################################################################################################################

# recent code that cause error

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import re
from pathlib import Path
import warnings


def remove_outliers_by_site(df):
    """
    Replace values outside the 1% and 99% percentiles with NaNs for numeric columns,
    calculated separately for each site.
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    
    # Exclude identifier columns that shouldn't be filtered
    exclude_cols = ['year', 'month', 'lat', 'long']
    numeric_cols = [col for col in numeric_cols if col not in exclude_cols]
    
    print(f"Applying 1-99% outlier removal to {len(numeric_cols)} numeric columns by site")
    
    total_outliers_removed = 0
    
    for col in numeric_cols:
        if col in df.columns:
            col_outliers = 0
            # Process each site separately
            for site in df['site_name'].unique():
                site_mask = df['site_name'] == site
                site_data = df.loc[site_mask, col]
                
                if len(site_data) > 10:  # Only process if we have enough data
                    # Calculate 1st and 99th percentiles for this site
                    lower_bound = site_data.quantile(0.01)
                    upper_bound = site_data.quantile(0.99)
                    
                    # Create mask for outliers in this site and column
                    outlier_mask = site_mask & (
                        (df[col] < lower_bound) | 
                        (df[col] > upper_bound)
                    ) & df[col].notna()
                    
                    # Count and remove outliers
                    site_outliers = outlier_mask.sum()
                    df.loc[outlier_mask, col] = np.nan
                    
                    col_outliers += site_outliers
            
            total_outliers_removed += col_outliers
            print(f"  {col}: Removed {col_outliers} outliers")
    
    print(f"Total outliers removed: {total_outliers_removed}")
    return df




def plot_numeric_variables_by_site(df, plot_path):
    """
    Create boxplots for all numeric variables across different sites.
    """
    # Identify numeric columns (exclude text/identifier columns)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    
    # Remove identifier columns that happen to be numeric
    exclude_cols = ['year', 'month', 'lat', 'long']
    numeric_cols = [col for col in numeric_cols if col not in exclude_cols]
    
    if not numeric_cols:
        print("No numeric variables found to plot")
        return
    
    # Create subplots - 4 plots per row
    n_cols = 4
    n_rows = (len(numeric_cols) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
    fig.suptitle('Distribution of Numeric Variables Across Sites', fontsize=16, y=0.98)
    
    # Flatten axes array for easy indexing
    if n_rows > 1:
        axes = axes.flatten()
    else:
        axes = [axes] if n_cols == 1 else axes
    
    for i, col in enumerate(numeric_cols):
        if i < len(axes):
            ax = axes[i]
            
            # Create boxplot
            sns.boxplot(data=df, x='site_name', y=col, ax=ax)
            ax.set_title(f'{col} by Site')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.tick_params(axis='x', labelsize=8)
    
    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'numeric_variables_by_site.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created boxplot for {len(numeric_cols)} numeric variables")

def monthly_drivers(min_coverage=0.50, gs_coverage=0.25):
    """
    Aggregate driver CSVs to monthly resolution, add monthly LAI, join with reference data,
    compute WUE_tra and WUE_eva, and create diagnostic boxplots.
    Apply different coverage thresholds: 50% for most variables, 25% for gs_iPM and gs_FG.
    """
    
    # Define paths
    driver_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split\3-growing_season_transp_Evap"
    reference_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
    lai_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai\growing_season_LAI"
    plot_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\plots_drivers"
    output_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_with_drivers.csv"
    
    os.makedirs(plot_path, exist_ok=True)
    log_messages = []
    
    # 1. Load reference data first
    log_messages.append("Loading reference data...")
    try:
        reference_df = pd.read_csv(reference_path)
        target_months = reference_df[['site_name', 'year', 'month']].drop_duplicates()
        log_messages.append(f"Reference data has {len(target_months)} unique site-month combinations")
        log_messages.append(f"Reference sites: {reference_df['site_name'].unique()}")
    except Exception as e:
        raise ValueError(f"Error loading reference data: {str(e)}")
    
    # 2. Process driver files
    log_messages.append("Processing driver files...")
    driver_files = glob.glob(os.path.join(driver_path, "*.csv"))
    
    if not driver_files:
        raise ValueError(f"No CSV files found in {driver_path}")
    
    all_driver_data = []
    
    for file in driver_files:
        filename = os.path.basename(file)
        
        # Fix site name extraction - remove .csv extension and keep full name
        site_name = filename.replace('.csv', '')
        
        try:
            # Read the file
            df = pd.read_csv(file)
            
            if 'DateTime' not in df.columns:
                log_messages.append(f"Skipping {filename}: No DateTime column")
                continue
            
            # Parse DateTime - try multiple formats
            try:
                # First try the format that matches your data: YYYY-MM-DD HH:MM:SS
                df['datetime'] = pd.to_datetime(df['DateTime'], format='%Y-%m-%d %H:%M:%S')
            except:
                try:
                    # Try without seconds
                    df['datetime'] = pd.to_datetime(df['DateTime'], format='%Y-%m-%d %H:%M')
                except:
                    try:
                        # Try another common format
                        df['datetime'] = pd.to_datetime(df['DateTime'], format='%m/%d/%Y %H:%M')
                    except:
                        # Fall back to pandas automatic parsing
                        df['datetime'] = pd.to_datetime(df['DateTime'])
            
            df['year'] = df['datetime'].dt.year
            df['month'] = df['datetime'].dt.month
            df['site_name'] = site_name
            
            log_messages.append(f"Parsed dates for {filename}: {df['datetime'].min()} to {df['datetime'].max()}")
            
            # Filter to only reference months for this site
            site_targets = target_months[target_months['site_name'] == site_name]
            if not site_targets.empty:
                # Merge to keep only months that exist in reference
                merged = pd.merge(df, site_targets, on=['site_name', 'year', 'month'])
                if not merged.empty:
                    all_driver_data.append(merged)
                    log_messages.append(f"Added {len(merged)} rows from {filename} for site {site_name}")
                else:
                    log_messages.append(f"No matching months between driver {filename} and reference for site {site_name}")
            else:
                log_messages.append(f"No reference months found for site {site_name}")
                
        except Exception as e:
            log_messages.append(f"Error processing {filename}: {str(e)}")
            continue
    
    if not all_driver_data:
        log_messages.append("WARNING: No driver data found for any reference months!")
        driver_monthly = pd.DataFrame(columns=['site_name', 'year', 'month'])
    else:
        # Combine all driver data
        driver_df = pd.concat(all_driver_data, ignore_index=True)
        log_messages.append(f"Total driver rows: {len(driver_df)}")
        log_messages.append(f"Unique site-month combinations in driver data: {len(driver_df[['site_name', 'year', 'month']].drop_duplicates())}")
        
        # Define aggregation rules
        sum_cols = ['precip_mm', 'Evap_pen', 'Trans_pen', 'Co2']
        mean_cols = ['NETRAD_f', 'WS', 'WD', 'H_f', 'PAR_f', 'RH_f', 'PA_f']
        median_cols = ['gs_iPM', 'gs_FG']
        
        # Check which columns exist
        available_sum_cols = [col for col in sum_cols if col in driver_df.columns]
        available_mean_cols = [col for col in mean_cols if col in driver_df.columns]
        available_median_cols = [col for col in median_cols if col in driver_df.columns]
        
        log_messages.append(f"Columns to aggregate - Sum: {available_sum_cols}, Mean: {available_mean_cols}, Median: {available_median_cols}")
        
        # Create aggregation dictionary
        agg_dict = {}
        for col in available_sum_cols:
            agg_dict[col] = 'sum'
        for col in available_mean_cols:
            agg_dict[col] = 'mean'
        for col in available_median_cols:
            agg_dict[col] = 'median'
        
        # Group and aggregate
        grouped = driver_df.groupby(['site_name', 'year', 'month'])
        driver_monthly = grouped.agg(agg_dict).reset_index()
        
        # Apply coverage thresholds - different for gs variables vs others
        for col in agg_dict.keys():
            for (site, year, month), group in grouped:
                valid_count = group[col].count()
                total_count = len(group)
                
                if total_count > 0:
                    # Use 25% threshold for gs variables, 50% for others
                    coverage_threshold = gs_coverage if col in ['gs_iPM', 'gs_FG'] else min_coverage
                    
                    if valid_count / total_count < coverage_threshold:
                        mask = (driver_monthly['site_name'] == site) & \
                               (driver_monthly['year'] == year) & \
                               (driver_monthly['month'] == month)
                        driver_monthly.loc[mask, col] = np.nan
        
        log_messages.append(f"Aggregated driver data with {len(driver_monthly)} monthly records")
        log_messages.append(f"Applied coverage thresholds: {min_coverage*100}% for most variables, {gs_coverage*100}% for gs_iPM and gs_FG")
    
    # 3. Process LAI files
    log_messages.append("Processing LAI files...")
    lai_files = glob.glob(os.path.join(lai_path, "*.csv"))
    
    all_lai_data = []
    for file in lai_files:
        # Extract site name without extension
        site_name = os.path.basename(file).split('.')[0]
        
        try:
            df = pd.read_csv(file)
            
            if 'Date' not in df.columns:
                log_messages.append(f"Skipping LAI file {file}: No Date column")
                continue
            
            df['date'] = pd.to_datetime(df['Date'])
            df['year'] = df['date'].dt.year
            df['month'] = df['date'].dt.month
            df['site_name'] = site_name
            
            # Filter to reference months
            site_targets = target_months[target_months['site_name'] == site_name]
            if not site_targets.empty:
                merged = pd.merge(df, site_targets, on=['site_name', 'year', 'month'])
                if not merged.empty:
                    all_lai_data.append(merged)
        except Exception as e:
            log_messages.append(f"Error reading LAI file {file}: {str(e)}")
    
    if all_lai_data:
        lai_df = pd.concat(all_lai_data, ignore_index=True)
        lai_grouped = lai_df.groupby(['site_name', 'year', 'month'])
        lai_monthly = lai_grouped['lai'].mean().reset_index()
        log_messages.append(f"Aggregated LAI data for {len(lai_monthly)} monthly records")
    else:
        log_messages.append("No LAI data found")
        lai_monthly = pd.DataFrame(columns=['site_name', 'year', 'month', 'lai'])
    
    # 4. Merge all data
    log_messages.append("Merging data...")
    
    # Start with reference, left join drivers, then LAI
    merged_df = reference_df.copy()
    
    if not driver_monthly.empty:
        merged_df = pd.merge(merged_df, driver_monthly, on=['site_name', 'year', 'month'], how='left')
    
    merged_df = pd.merge(merged_df, lai_monthly, on=['site_name', 'year', 'month'], how='left')
    
    # 5. Calculate new variables
    if 'Trans_pen' in merged_df.columns and 'GPP' in merged_df.columns:
        merged_df['WUE_tra'] = np.where(
            (merged_df['Trans_pen'] > 0) & merged_df['Trans_pen'].notna() & merged_df['GPP'].notna(),
            merged_df['GPP'] / merged_df['Trans_pen'],
            np.nan
        )
    
    if 'Evap_pen' in merged_df.columns and 'GPP' in merged_df.columns:
        merged_df['WUE_eva'] = np.where(
            (merged_df['Evap_pen'] > 0) & merged_df['Evap_pen'].notna() & merged_df['GPP'].notna(),
            merged_df['GPP'] / merged_df['Evap_pen'],
            np.nan
        )
    
    # 6. Remove outliers using 1-99% percentiles by site
    log_messages.append("Removing outliers using 1-99% percentiles by site...")
    merged_df = remove_outliers_by_site(merged_df)
    
    # 7. Final column order
    ref_cols = ['site_name', 'year', 'month', 'NEE', 'GPP', 'Reco', 'ET', 'NEP', 'State', 'biome', 
                'salinity_ppt', 'salinity_fine', 'salini_coarse', 'WUE', 'CUE', 'Tair_f', 'VPD_f', 
                'Rg_f', 'avg_length', 'avg_sos', 'avg_eos', 'lat', 'long', 'climate', 'climate_2', 'salinity_con']
    
    driver_cols = ['NETRAD_f', 'WS', 'WD', 'Co2', 'H_f', 'PAR_f', 'RH_f', 'PA_f', 'precip_mm', 
                  'gs_iPM', 'gs_FG', 'Evap_pen', 'Trans_pen']
    
    extra_cols = ['lai', 'WUE_tra', 'WUE_eva']
    
    # Reorder columns
    final_cols = []
    for col_list in [ref_cols, driver_cols, extra_cols]:
        for col in col_list:
            if col in merged_df.columns and col not in final_cols:
                final_cols.append(col)
    
    # Add any remaining columns
    for col in merged_df.columns:
        if col not in final_cols and col not in ['NEE_f', 'LE_f']:
            final_cols.append(col)
    
    merged_df = merged_df[final_cols]
    
    # 8. Create diagnostic plots
    log_messages.append("Creating diagnostic plots...")
    try:
        plot_numeric_variables_by_site(merged_df, plot_path)
        log_messages.append("Successfully created numeric variables boxplot")
    except Exception as e:
        log_messages.append(f"Error creating plots: {str(e)}")
    
    # 9. Save and log results
    log_messages.append("Saving final dataset...")
    merged_df.to_csv(output_path, index=False)
    
    # Data availability summary
    driver_vars = ['NETRAD_f', 'WS', 'WD', 'Co2', 'H_f', 'PAR_f', 'RH_f', 'PA_f', 'precip_mm', 
                  'gs_iPM', 'gs_FG', 'Evap_pen', 'Trans_pen']
    extra_vars = ['lai', 'WUE_tra', 'WUE_eva']
    
    log_messages.append("\nData availability summary:")
    for col in driver_vars + extra_vars:
        if col in merged_df.columns:
            non_na = merged_df[col].notna().sum()
            total = len(merged_df)
            log_messages.append(f"{col}: {non_na}/{total} ({non_na/total:.1%})")
        else:
            log_messages.append(f"{col}: Column not found")
    
    # Print all logs
    for msg in log_messages:
        print(msg)
    
    return merged_df

# Execute
if __name__ == "__main__":
    monthly_data = monthly_drivers(min_coverage=0.50, gs_coverage=0.25)

###############################################################################################################


















































write a function to create monthly driver data
 
call function monthly_drivers

path to folder for driver data
\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split\3-growing_season_transp_Evap

this folder have mutiple csv files with unique name. all files have same variable names


variable that will be sum by month
precip_mm, Co2, Evap_pen, 	Trans_pen


variable that will be average by month: dont change name of these variable even when process on monthly scale
NETRAD_f	WS	WD	Co2	H_f	PAR_f	RH_f	PA_f	

for these variable calculate median by month: dont change name of variables even when median is comupted
    gs_iPM	gs_FG


 next step merging data with data in another folder

refernce data in another folder
\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products
name of csv file for reference data WUE_CUE_monthly

reference data is already on monthly scale so no need to change that

include these variables  from reference data at the start of final df same order


site_name	year	month	NEE	GPP	Reco	ET	NEP	State	biome	salinity_ppt	salinity_fine	salini_coarse	WUE	CUE	Tair_f	VPD_f	Rg_f	avg_length	avg_sos	avg_eos	lat	long	climate	climate_2	salinity_con

look at site_name column in reference data to match with name of csv files in other folders in this function such as driver folder and lai data folder


path for lai data\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai\growing_season_LAI

from these folder read csv file . all files have unique site name . take column name lai and merge 
it with refrence data. lai data is weekly so need to calculate mean monthly lai 
lai folder path


new variables 

WUE_tra=monthly['GPP'] / monthly['Trans_pen']
WUE_eva=monthly['GPP'] / monthly['Evap_pen']


site_name	year	month	NEE	GPP	Reco	ET	NEP	State	biome	salinity_ppt	salinity_fine	salini_coarse
	WUE	CUE	Tair_f	VPD_f	Rg_f	avg_length	avg_sos	avg_eos	lat	long	climate	climate_2	salinity_con

NETRAD_f	WS	WD	Co2	H_f	PAR_f	RH_f	PA_f	
precip_mm	NEE_f	LE_f
gs_iPM	gs_FG	Evap_pen	Trans_pen lai
WUE_tra
WUE_eva





































write a function to filter data that belongs to growing season

look at these columns in csv files 
DateTime	Year	DoY	

Hour	Ustar	NETRAD_f	WS	WD	Co2	NEE	H_f	PAR_f	RH_f	PA_f	precip_mm	NEE_f	LE_f	Tair_f	VPD_f	Rg_f	GPP_DT	GPP_nt	Reco_DT	Reco_nt	ET	gs_iPM	gs_FG	LE_tr	Trans_dir	Evap_resid	Evap_pen	Trans_pen

use DoY to calulate growing season
to caluclaute growing season, you need data from another folder.
here is the folder for growing seaosn data

growing season data folder

\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products
here is name of csv file in this folder that has info for growing season 
phenofit_growing_season


site_name	year	avg_length	avg_sos	avg_eos

site_name matches CSV file name

growing season statt with avg_sos and end with avg_eos.
only include data that is part of growing season.

use average of column avg_sos and column avg_eos so that same growing season length for 
all year for a unique site

        # === Add avg_length unchanged per site-year ===
        monthly['avg_length'] = site_df['avg_length'].iloc[0]

save output files in this folder
\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split\1-growing_season_transp_Evap

same name of site as in input folder
print range of these column ET Evap_pen	Trans_pen
yearly sum data  Evap_pen/ET
and Trans_pen/ET





final df should have following columns

DateTime	Year	DoY	Hour	Ustar	NETRAD_f	WS	WD	Co2	NEE	H_f	PAR_f	RH_f	PA_f	precip_mm	NEE_f	LE_f	Tair_f	VPD_f	Rg_f	GPP_DT	GPP_nt	Reco_DT	Reco_nt	ET	gs_iPM	gs_FG	LE_tr	Trans_dir	Evap_resid	Evap_pen	Trans_pen























