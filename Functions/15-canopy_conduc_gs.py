# -*- coding: utf-8 -*-
"""
Created on Sat Sep  6 08:45:34 2025

@author: ammar
"""
## first add precip data into your filled variable set so that canopy conductance can be calculated 

import pandas as pd
import os
from pathlib import Path

def merge_precip_with_gapfilled_data(precip_folder, gap_filled_folder, output_folder):
    """
    Merge precipitation data with gap-filled data and save results.
    
    Parameters:
    precip_folder (str): Path to folder containing precipitation CSV files
    gap_filled_folder (str): Path to folder containing gap-filled CSV files  
    output_folder (str): Path to folder where output files will be saved
    """
    
    # Convert paths to Path objects for easier handling
    precip_path = Path(precip_folder)
    gap_filled_path = Path(gap_filled_folder)
    output_path = Path(output_folder)
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get list of site files from both folders
    precip_files = [f for f in precip_path.glob('*.csv') if f.is_file()]
    
    for precip_file in precip_files:
        site_name = precip_file.stem  # Get site name from filename
        
        # Check if corresponding gap-filled file exists
        gap_filled_file = gap_filled_path / f"{site_name}.csv"
        
        if not gap_filled_file.exists():
            print(f"Warning: No gap-filled file found for site {site_name}")
            continue
        
        try:
            # Read precipitation data
            precip_df = pd.read_csv(precip_file)
            
            # Check if required columns exist
            if not all(col in precip_df.columns for col in ['site', 'date', 'precip_mm']):
                print(f"Warning: Missing required columns in {precip_file.name}")
                continue
            
            # Convert date column to datetime and set as index
            precip_df['date'] = pd.to_datetime(precip_df['date'])
            precip_df = precip_df.set_index('date')
            
            # Read gap-filled data
            gap_filled_df = pd.read_csv(gap_filled_file)
            
            # Convert DateTime column to datetime
            gap_filled_df['DateTime'] = pd.to_datetime(gap_filled_df['DateTime'])
            
            # Create a date column for merging
            gap_filled_df['date'] = gap_filled_df['DateTime'].dt.date
            
            # Convert date column to datetime for merging
            gap_filled_df['date'] = pd.to_datetime(gap_filled_df['date'])
            
            # Calculate number of half-hour intervals per day
            half_hours_per_day = 48  # 24 hours * 2 (half-hour intervals)
            
            # Distribute daily precipitation evenly across half-hour intervals
            precip_daily = precip_df[['precip_mm']].copy()
            precip_daily['precip_half_hour'] = precip_daily['precip_mm'] / half_hours_per_day
            
            # Merge precipitation data with gap-filled data
            merged_df = gap_filled_df.merge(
                precip_daily[['precip_half_hour']], 
                left_on='date', 
                right_index=True, 
                how='left'
            )
            
            # Rename the precipitation column
            merged_df = merged_df.rename(columns={'precip_half_hour': 'precip_mm'})
            
            # Drop the temporary date column
            merged_df = merged_df.drop(columns=['date'])
            
            # Save the merged data
            output_file = output_path / f"{site_name}.csv"
            merged_df.to_csv(output_file, index=False)
            
            print(f"Successfully processed and saved: {site_name}")
            
        except Exception as e:
            print(f"Error processing site {site_name}: {str(e)}")
            continue

# Usage example
if __name__ == "__main__":
    # Define your folder paths
    precip_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\precip"
    gap_filled_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\reddy_proc\gap_filled"
    output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc"
    
    # Run the function
    merge_precip_with_gapfilled_data(precip_folder, gap_filled_folder, output_folder)

#################################################################################################################





# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:13:38 2025

@author: ammar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import os
import glob

# Set plot style
plt.style.use('default')
sns.set_palette("husl")

def lambda_Jkg(Tair):
    """Calculate latent heat of vaporization in J/kg from air temperature."""
    return (2.501 - 0.002361 * Tair) * 1e6  # Convert MJ/kg to J/kg

def growseas_gs(flux_dat, month_min, month_max, LE_min, PAR_min, output_format):
    """
    Calculate canopy conductance (gs) from eddy flux data using 
    inverted Penman-Monteith (iPM) and Flux Gradient (FG) methods.
    Follows paper's approach: use total LE for energy closure, 
    and apply dry-canopy filtering (no precipitation) BEFORE calculations.
    """
    
    # Create a copy to avoid modifying original data
    flux_subset = flux_dat.copy()
    
    # Convert datetime and extract time components
    flux_subset['DateTime'] = pd.to_datetime(flux_subset['DateTime'])
    flux_subset['year'] = flux_subset['DateTime'].dt.year
    flux_subset['month'] = flux_subset['DateTime'].dt.month
    flux_subset['hour'] = flux_subset['DateTime'].dt.hour
    
    # Apply dry-canopy filtering FIRST: only keep hours with no precipitation
    dry_canopy_mask = (flux_subset['precip_mm'] == 0)
    flux_subset = flux_subset[dry_canopy_mask].copy()
    
    # VPD & RH hygiene - SAFETY CHECK
    flux_subset['VPD_f'] = np.clip(flux_subset['VPD_f'], 0, None)       # hPa, prevent negative VPD
    flux_subset['RH_f'] = flux_subset['RH_f'].clip(0, 100)              # %, ensure RH within bounds
    
    # Filter by LE and PAR thresholds (using total LE_f for filtering)
    filter_conditions = [
        (~flux_subset['LE_f'].isna()) & 
        (flux_subset['LE_f'] >= LE_min),
        (~flux_subset['PAR_f'].isna()) & 
        (flux_subset['PAR_f'] >= PAR_min),
        (~flux_subset['WS'].isna()) &  # Use 'WS' instead of 'WS_f'
        (flux_subset['WS'] > 0)        # Avoid division by zero
    ]
    
    for condition in filter_conditions:
        flux_subset = flux_subset[condition]
    
    # Define constants
    sidesOfLeafWithStomata = 1
    gasConstR = 8.314472  # J/mol/K
    c_p = 1012  # J/kg/K, specific heat of air
    airMolarMass_dry = 0.02897  # kg/mol, molar mass of dry air
    
    # Set lambda (latent heat of vaporization) - FIXED: always use J/kg
    if 'lambda' in flux_subset.columns:
        lambda_val = flux_subset['lambda'] * 1e6  # Convert MJ/kg → J/kg
    else:
        lambda_val = 2.45e6  # J/kg
    
    Sc_CO2 = 1.05  # Schmidt number for CO2
    Sc_H2O = Sc_CO2 / 1.57278  # Schmidt number for H2O
    Pr = 0.71  # Prandtl number
    
    # Energy budget calculations using TOTAL LE_f (not partitioned)
    G_guess = flux_subset['NETRAD_f'] * 15 / 700  # W/m², rough estimate of ground heat flux
    TotalTurbHeatFlux = flux_subset['H_f'] + flux_subset['LE_f']  # W/m², use total LE_f
    AvailableEnergy = flux_subset['NETRAD_f'] - G_guess  # W/m²
    
    # Use consistent G estimate for PM inversion
    Rn_minus_G = flux_subset['NETRAD_f'] - G_guess  # W/m², for use in PM equation
    
    # Handle NaN values in the regression input
    valid_mask = (~np.isnan(AvailableEnergy)) & (~np.isnan(TotalTurbHeatFlux))
    
    if np.sum(valid_mask) > 10:  # Need sufficient data points for regression
        # Estimate slope using total LE_f (not partitioned)
        X = AvailableEnergy[valid_mask].values.reshape(-1, 1)
        y = TotalTurbHeatFlux[valid_mask].values
        
        hourlyEbudgratio = LinearRegression().fit(X, y)
        hourlyEnergyBudgetRatio = hourlyEbudgratio.coef_[0]
    else:
        # Fallback: use 1.0 if not enough data for regression
        print(f"Warning: Insufficient data for energy balance regression. Using default ratio of 1.0")
        hourlyEnergyBudgetRatio = 1.0
    
    # Flux correction using total LE_f
    H_hCorr = flux_subset['H_f'] / hourlyEnergyBudgetRatio  # W/m²
    LE_hCorr = flux_subset['LE_f'] / hourlyEnergyBudgetRatio  # W/m², use total LE_f
    
    # Air physics calculations - use consistent Pa units throughout
    P_pa = flux_subset['PA_f'] * 1000  # Pa, convert kPa to Pa
    T_K = flux_subset['Tair_f'] + 273.15  # K
    
    SatVP = 100 * 6.112 * np.exp(17.62 * flux_subset['Tair_f'] / (243.12 + flux_subset['Tair_f']))  # Pa
    
    # FIXED: Correct SatVPslope calculation for Pa units
    SatVPslope = 4098.0 * SatVP / (flux_subset['Tair_f'] + 237.3)**2  # Pa/K (CORRECTED)
    
    # Convert VPD from hPa to Pa (VPD_f is in hPa)
    VPD_pa = flux_subset['VPD_f'] * 100.0  # Convert hPa to Pa
    
    # SAFETY CHECK: Ensure VP is non-negative and doesn't exceed SatVP
    VP = np.maximum(0.0, SatVP - VPD_pa)  # Pa, prevent negative vapor pressure
    VP = np.minimum(VP, SatVP)  # Ensure VP doesn't exceed saturation
    
    VP_n = VP  # Pa
    
    AirConc_dry = (P_pa - VP) / (gasConstR * T_K)  # mol/m³
    AirConc_wet = P_pa / (gasConstR * T_K)  # mol/m³
    
    HeatCapacity_dry = 1003 + (1008 - 1003) * ((flux_subset['Tair_f'] + 23.16) / 100)  # J/kg/K
    HeatCapacity_waterVapor = HeatCapacity_dry * (1 + 0.84 * flux_subset['RH_f'] / 100)  # J/kg/K
    AirDensity_dry = AirConc_dry * airMolarMass_dry  # kg/m³
    WaterVaporDensity = (AirConc_wet - AirConc_dry) * 0.018  # kg/m³
    AirDensity_wet = AirDensity_dry + WaterVaporDensity  # kg/m³
    HeatCapacity_wet = (AirDensity_dry * HeatCapacity_dry + WaterVaporDensity * HeatCapacity_waterVapor) / AirDensity_wet  # J/kg/K
    Gamma = HeatCapacity_wet * P_pa / (lambda_val * 0.62198)  # Pa/K, psychrometric constant
    
    # COASTAL WETLAND AERODYNAMIC PARAMETERIZATION
    # Using physically-based scaling for emergent coastal wetlands
    h_c = 0.6  # m - mean canopy height for coastal salt marsh (adjust based on your site!)
    
    # Scaling relationships for emergent wetlands
    d = 0.67 * h_c                    # displacement height
    z0m = 0.115 * h_c                 # momentum roughness length
    kB_inv = 4.0                      # kB⁻¹ parameter
    z0h = z0m * np.exp(-kB_inv)       # thermal roughness length
    z = h_c + 3.0                     # measurement height
    
    k = 0.4     # von Karman constant
    
    # Calculate terms for logarithmic wind profile with numerical safety
    num_m = (z - d) / z0m  # momentum term
    num_h = (z - d) / z0h  # heat/scalar term
    
    # Ensure valid logarithmic arguments
    valid_log = (num_m > 1.01) & (num_h > 1.01)
    
    # Calculate eddy conductance to heat (G_eh) in m/s using correct formulation
    G_eh = np.where(
        valid_log,
        (k**2 * flux_subset['WS']) / (np.log(num_m) * np.log(num_h)),
        np.nan
    )
    
    # Mild clipping to avoid exploding values - FIXED: more realistic bounds
    G_eh = np.clip(G_eh, 0.001, 1.0)  # Realistic bounds for wetland eddy conductance
    
    # Convert to molar units (mol m⁻² s⁻¹)
    air_molar_density = P_pa / (gasConstR * T_K)  # mol/m³
    G_eh_molar = air_molar_density * G_eh  # mol m⁻² s⁻¹
    
    # Eddy resistance to heat (s/m) - handle division by zero/NaN
    R_eh = np.where(G_eh > 0, 1.0 / G_eh, np.nan)
    
    # Boundary resistance to heat
    R_bH = 10  # s/m, boundary resistance to heat
    R_b = (2 / sidesOfLeafWithStomata) * R_bH * ((Sc_CO2 / Pr)**(2/3))  # s/m, boundary layer resistance for CO2
    R_bV = (2 / sidesOfLeafWithStomata) * R_bH * ((Sc_H2O / Pr)**(2/3))  # s/m, boundary layer resistance for H2O
    
    # Boundary conductance from resistance - FIXED: all use Pa
    G_bH = P_pa / (gasConstR * T_K * R_bH)  # mol/m²/s
    G_b = P_pa / (gasConstR * T_K * R_b)  # mol/m²/s
    G_bV = P_pa / (gasConstR * T_K * R_bV)  # mol/m²/s
    
    # Water flux calculation using LE_hCorr (DO NOT overwrite LE_hCorr)
    E = LE_hCorr / (lambda_val * 0.018)  # mol/m²/s, calculate E from LE_hCorr
    
    # Leaf temperature and VPD calculations - UPDATED to use R_eh
    # Handle NaN values in R_eh
    valid_R_eh = ~np.isnan(R_eh)
    Tair_n = np.where(
        valid_R_eh,
        (H_hCorr * R_eh / (AirDensity_wet * HeatCapacity_wet)) + flux_subset['Tair_f'],
        np.nan
    )  # °C
    
    Tair_leaf = np.where(
        valid_R_eh,
        (H_hCorr * R_bH / (AirDensity_wet * HeatCapacity_wet)) + Tair_n,
        np.nan
    )  # °C
    
    SatVP_leaf = np.where(
        ~np.isnan(Tair_leaf),
        100 * 6.112 * np.exp(17.62 * Tair_leaf / (243.12 + Tair_leaf)),
        np.nan
    )  # Pa
    
    LeafAirVPD = np.where(
        ~np.isnan(SatVP_leaf),
        SatVP_leaf - VP_n,
        np.nan
    )  # Pa
    
    LeafAirConcDiff = np.where(
        ~np.isnan(Tair_n) & ~np.isnan(SatVP_leaf),
        (1 / (gasConstR * (Tair_n + 273.15))) * (SatVP_leaf - VP_n),
        np.nan
    )  # mol/m³
    
    # Canopy resistance to water vapor (Flux Gradient method)
    # Use calculated E to derive resistance, but don't modify original LE_hCorr
    R_sV = np.where(
        ~np.isnan(LeafAirConcDiff) & ~np.isnan(E) & (E > 0),
        LeafAirConcDiff / E - R_bV,
        np.nan
    )  # s/m, canopy resistance to water vapor
    
    R_s = np.where(
        ~np.isnan(R_sV),
        R_sV * 1.57278,
        np.nan
    )  # s/m, canopy resistance to CO2
    
    # Canopy conductance to water vapor (Flux Gradient method) - FIXED: uses Pa
    G_sV_hCorr = np.where(
        ~np.isnan(R_sV) & ~np.isnan(Tair_leaf) & (R_sV > 0),
        P_pa / (gasConstR * (Tair_leaf + 273.15) * R_sV),
        np.nan
    )  # mol/m²/s
    
    G_s = np.where(
        ~np.isnan(R_s) & ~np.isnan(Tair_leaf) & (R_s > 0),
        P_pa / (gasConstR * (Tair_leaf + 273.15) * R_s),
        np.nan
    )  # mol/m²/s
    
    # Penman-Monteith inversion method using total LE_f - FIXED: use R_eh and LE_hCorr
    # Handle potential division by zero or NaN in PM calculation
    valid_pm_mask = (~np.isnan(LE_hCorr)) & (LE_hCorr > 0) & (~np.isnan(Gamma)) & (Gamma > 0) & (~np.isnan(R_eh))
    
    R_sV_PM = np.zeros(len(flux_subset))  # s/m
    if np.any(valid_pm_mask):
        VPD_term = SatVP[valid_pm_mask] - VP[valid_pm_mask]  # Pa
        raH = R_eh[valid_pm_mask]  # aerodynamic resistance to heat
        
        R_sV_PM[valid_pm_mask] = (
            (SatVPslope[valid_pm_mask] * (Rn_minus_G[valid_pm_mask] - LE_hCorr[valid_pm_mask]) 
             + AirDensity_wet[valid_pm_mask] * HeatCapacity_wet[valid_pm_mask] * VPD_term / raH)
            / (Gamma[valid_pm_mask] * LE_hCorr[valid_pm_mask])
            - 1.0
        ) * raH - R_bV
    
    G_sV_PM_hCorr = np.zeros(len(flux_subset))  # mol/m²/s
    valid_g_mask = (~np.isnan(R_sV_PM)) & (R_sV_PM > 0)
    if np.any(valid_g_mask):
        G_sV_PM_hCorr[valid_g_mask] = P_pa[valid_g_mask] / (gasConstR * T_K[valid_g_mask] * R_sV_PM[valid_g_mask])  # FIXED: uses Pa
    
    # Compile results with units in comments
    flux_subset['gs_iPM'] = G_sV_PM_hCorr  # mol/m²/s, canopy conductance (iPM method)
    flux_subset['gs_FG'] = G_sV_hCorr      # mol/m²/s, canopy conductance (FG method)
    flux_subset['E'] = E                   # mol/m²/s, water flux
    flux_subset['LE_tr'] = LE_hCorr        # W/m², transpiration component (≈ LE_f during dry canopy)
    
    # Filter gs values: set to NaN if less than 0 or greater than 10
    flux_subset.loc[(flux_subset['gs_iPM'] < 0) | (flux_subset['gs_iPM'] > 10), 'gs_iPM'] = np.nan
    flux_subset.loc[(flux_subset['gs_FG'] < 0) | (flux_subset['gs_FG'] > 10), 'gs_FG'] = np.nan
    
    # Remove temporary columns
    cols_to_drop = ['date'] if 'date' in flux_subset.columns else []
    flux_subset = flux_subset.drop(columns=cols_to_drop)
    
    # Calculate statistics based on output format
    if output_format == "all":
        return flux_subset
    
    elif output_format == "monthly":
        gs_stats = flux_subset.groupby(['year', 'month']).agg({
            'gs_iPM': ['median', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
            'gs_FG': ['median', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
            'LE_f': 'median',
            'E': 'median',
            'LE_tr': 'median'
        }).reset_index()
        
        gs_stats.columns = ['year', 'month', 'gs_iPM_med', 'gs_iPM_25', 'gs_iPM_75', 
                           'gs_FG_med', 'gs_FG_25', 'gs_FG_75', 'LE_f_med', 'E_med', 
                           'LE_tr_med']
        return gs_stats
    
    elif output_format == "daily":
        flux_subset['date'] = flux_subset['DateTime'].dt.date
        gs_stats = flux_subset.groupby('date').agg({
            'gs_iPM': ['median', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
            'gs_FG': ['median', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
            'LE_f': 'median',
            'E': 'median',
            'LE_tr': 'median'
        }).reset_index()
        
        gs_stats.columns = ['date', 'gs_iPM_med', 'gs_iPM_25', 'gs_iPM_75', 
                           'gs_FG_med', 'gs_FG_25', 'gs_FG_75', 'LE_f_med', 'E_med', 
                           'LE_tr_med']
        return gs_stats
    
    else:
        raise ValueError("No output file format (all, daily, monthly) given.")

def apply_transpiration_mask_advanced(df):
    """
    Apply advanced masking to transpiration values based on precipitation events.
    Sets transpiration to NaN during and after precipitation events.
    """
    # Create a copy to avoid modifying the original
    df = df.copy()
        
    df['transpiration_flag'] = 0  # 0 = valid, 1 = masked due to precipitation
    
    # Find precipitation events
    precip_events = df['precip_mm'] > 0
    
    # Mask transpiration during precipitation events
    df.loc[precip_events, 'transpiration_flag'] = 1
    
    # Also mask for a period after precipitation (e.g., 6 hours)
    post_precip_hours = 6
    for i in range(len(df)):
        if precip_events.iloc[i]:
            # Mask the next 'post_precip_hours' hours
            end_idx = min(i + post_precip_hours, len(df))
            df.iloc[i:end_idx, df.columns.get_loc('transpiration_flag')] = 1
    
    return df

def plot_comprehensive_gs_analysis(gs_data, site_name, output_folder):
    """Create comprehensive plots for canopy conductance analysis."""
    
    # Calculate yearly statistics
    yearly_stats = gs_data.groupby('year').agg({
        'gs_iPM': 'median',
        'gs_FG': 'median',
        'LE_f': 'median',
        'LE_tr': 'median',
        'E': 'median'
    }).reset_index()
    
    # Create comprehensive figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'Canopy Conductance Analysis - {site_name}', fontsize=16, fontweight='bold')
    
    # Plot 1: Yearly median canopy conductance (both methods)
    ax1.plot(yearly_stats['year'], yearly_stats['gs_iPM'], 'o-', linewidth=2, markersize=8, label='iPM Method')
    ax1.plot(yearly_stats['year'], yearly_stats['gs_FG'], 's-', linewidth=2, markersize=8, label='FG Method')
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Canopy Conductance (mol m⁻² s⁻¹)', fontsize=12)
    ax1.set_title('Yearly Median Canopy Conductance\n(Dry Canopy Conditions Only)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Yearly median latent heat flux components
    ax2.plot(yearly_stats['year'], yearly_stats['LE_f'], 'o-', linewidth=2, markersize=8, label='Total LE', color='blue')
    ax2.plot(yearly_stats['year'], yearly_stats['LE_tr'], 's-', linewidth=2, markersize=8, label='LE_tr', color='green')
    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Latent Heat Flux (W m⁻²)', fontsize=12)
    ax2.set_title('Yearly Median Latent Heat Flux\n(Dry Canopy Conditions)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Relationship between LE_tr and canopy conductance
    valid_data = gs_data[(gs_data['gs_iPM'] > 0) & (gs_data['LE_tr'] > 0)]
    if len(valid_data) > 0:
        # Use VPD_f instead of VPDl since we removed VPDl
        scatter = ax3.scatter(valid_data['LE_tr'], valid_data['gs_iPM'], alpha=0.5, s=10, 
                            c=valid_data['VPD_f'], cmap='viridis')
        ax3.set_xlabel('Transpiration (LE_tr, W m⁻²)', fontsize=12)
        ax3.set_ylabel('Canopy Conductance (mol m⁻² s⁻¹)', fontsize=12)
        ax3.set_title('LE_tr vs Canopy Conductance\n(Color: VPD_f)')
        ax3.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax3, label='VPD_f (hPa)')
    
    # Plot 4: Diurnal pattern
    diurnal_stats = gs_data.groupby('hour').agg({
        'gs_iPM': 'median',
        'LE_tr': 'median'
    }).reset_index()
    
    ax4.plot(diurnal_stats['hour'], diurnal_stats['gs_iPM'], 'o-', linewidth=2, color='blue')
    ax4.set_xlabel('Hour of Day', fontsize=12)
    ax4.set_ylabel('Canopy Conductance (mol m⁻² s⁻¹)', fontsize=12, color='blue')
    ax4.tick_params(axis='y', labelcolor='blue')
    ax4.grid(True, alpha=0.3)
    
    ax4_2 = ax4.twinx()
    ax4_2.plot(diurnal_stats['hour'], diurnal_stats['LE_tr'], 's-', color='green', linewidth=2)
    ax4_2.set_ylabel('LE_tr (W m⁻²)', fontsize=12, color='green')
    ax4_2.tick_params(axis='y', labelcolor='green')
    
    ax4.set_title('Diurnal Patterns\n(Dry Canopy Conditions)')
    
    plt.tight_layout()
    
    # Save plot
    plot_path = Path(output_folder) / f"{site_name}_dry_canopy_analysis.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved comprehensive analysis plot: {plot_path}")

def process_canopy_conductance_for_sites(input_folder, output_folder, month_min=6, month_max=8, LE_min=50, PAR_min=500):
    """
    Process canopy conductance for all CSV files in the specified folder.
    Uses dry-canopy filtering approach (no precipitation) as in the paper.
    """
    
    input_path = Path(input_folder)
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    csv_files = list(input_path.glob('*.csv'))
    
    for csv_file in csv_files:
        try:
            print(f"Processing file: {csv_file.name}")
            site_name = csv_file.stem
            
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Process canopy conductance using dry-canopy approach
            result = growseas_gs(df, month_min=month_min, month_max=month_max,
                                LE_min=LE_min, PAR_min=PAR_min, output_format="all")
            
            # Save the results with the same filename in output folder
            output_file = output_path / f"{site_name}.csv"
            result.to_csv(output_file, index=False)
            
            print(f"Successfully processed and saved: {output_file.name}")
            print(f"Number of dry-canopy records: {len(result)}")
            print(f"Canopy conductance range - iPM: {result['gs_iPM'].min():.3f} to {result['gs_iPM'].max():.3f} mol m⁻² s⁻¹")
            print(f"Canopy conductance range - FG: {result['gs_FG'].min():.3f} to {result['gs_FG'].max():.3f} mol m⁻² s⁻¹")
            print(f"LE_tr range: {result['LE_tr'].min():.1f} to {result['LE_tr'].max():.1f} W m⁻²")
            
            # Create comprehensive plots
            plot_comprehensive_gs_analysis(result, site_name, output_path)
            
            print("-" * 50)
            
        except Exception as e:
            print(f"Error processing file {csv_file.name}: {str(e)}")
            continue

def write_log(message):
    """Helper function to write log messages."""
    print(message)

def process_ameri_files():
    # Define paths
    input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc"
    fill_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_fill"
    output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated"
    
    os.makedirs(output_folder, exist_ok=True)
    input_files = glob.glob(os.path.join(input_folder, "*.csv"))
    
    for input_file in input_files:
        try:
            site_name = os.path.basename(input_file).split('.')[0]
            write_log(f"Processing site: {site_name}")
            
            # Read input file
            df = pd.read_csv(input_file)
            
            if 'DateTime' not in df.columns:
                write_log(f"DateTime column not found in {input_file}, skipping")
                continue
            
            df['DateTime'] = pd.to_datetime(df['DateTime'])
            
            # Add Year, DoY, Hour columns
            datetime_cols = ['Year', 'DoY', 'Hour']
            existing_cols = df.columns.tolist()
            
            for col in datetime_cols:
                if col not in existing_cols:
                    if col == 'Year':
                        df.insert(1, 'Year', df['DateTime'].dt.year)
                    elif col == 'DoY':
                        pos = 2 if 'Year' in df.columns else 1
                        df.insert(pos, 'DoY', df['DateTime'].dt.dayofyear)
                    elif col == 'Hour':
                        pos = 3 if all(x in df.columns for x in ['Year', 'DoY']) else \
                              2 if 'Year' in df.columns or 'DoY' in df.columns else 1
                        df.insert(pos, 'Hour', df['DateTime'].dt.hour + df['DateTime'].dt.minute/60)
                    write_log(f"Added column: {col}")
            
            # Find fill file
            fill_file_pattern = os.path.join(fill_folder, f"{site_name}_fill.csv")
            fill_files = glob.glob(fill_file_pattern)
            
            if fill_files:
                fill_file = fill_files[0]
                write_log(f"Found fill file: {fill_file}")
                
                fill_df = pd.read_csv(fill_file)
                
                if 'DateTime' in fill_df.columns:
                    fill_df['DateTime'] = pd.to_datetime(fill_df['DateTime'])
                    
                    columns_to_merge = ['NEE_f', 'LE_f', 'Tair_f', 'VPD_f', 'Rg_f', 'GPP_DT', 'GPP_nt', 'Reco_DT', 'Reco_nt', 'precip_mm']
                    
                    for col in columns_to_merge:
                        if col in fill_df.columns:
                            if col in df.columns:
                                df = df.drop(columns=[col])
                                write_log(f"Replaced existing column: {col}")
                            temp_df = fill_df[['DateTime', col]].copy()
                            df = df.merge(temp_df, on='DateTime', how='left')
                            write_log(f"Merged column: {col}")
            
            # Fill missing columns
            required_cols = ['NEE_f', 'GPP_DT', 'Reco_DT', 'GPP_nt', 'Reco_nt', 'Tair_f', 'LE_f', 'precip_mm']
            for col in required_cols:
                if col not in df.columns:
                    df[col] = np.nan
                    write_log(f"Missing column '{col}' in {site_name}, filled with NaNs")
            
            # Unit conversions - FIXED: use lambda_Jkg instead of lambda_MJkg
            df['lambda'] = lambda_Jkg(df['Tair_f']) / 1e6  # Store as MJ/kg for consistency with ET calculation
            df['NEE'] = df['NEE_f'] * ((12 / 10**6) * 1800)
            df['GPP'] = df['GPP_DT'] * ((12 / 10**6) * 1800)
            df['Reco'] = df['Reco_DT'] * ((12 / 10**6) * 1800)
            df['GPP_nt'] = df['GPP_nt'] * ((12 / 10**6) * 1800)
            df['Reco_nt'] = df['Reco_nt'] * ((12 / 10**6) * 1800)
            df['ET'] = (df['LE_f'] * 1800.0) / (df['lambda'] * 1e6)  # lambda in MJ/kg
            
            # Gap-filling
            df['GPP'] = df['GPP'].fillna(df['GPP_nt'])
            df['Reco'] = df['Reco'].fillna(df['Reco_nt'])
            
            gap_limit = 48 * 14
            for col in ['GPP', 'Reco', 'ET']:
                if col in df.columns:
                    mask = df[col].isna()
                    filled = df[col].interpolate(limit=gap_limit, limit_direction='both')
                    df[col] = np.where(mask & filled.notna(), filled, df[col])
            
            # Calculate NPP
            df['NPP'] = df['GPP'] - df['Reco']
            
            # Calculate canopy conductance using the dry-canopy approach
            gs_result = growseas_gs(df, month_min=6, month_max=8, LE_min=50, PAR_min=500, output_format="all")
            
            # Add gs results to the main dataframe
            gs_columns = ['DateTime', 'gs_iPM', 'gs_FG', 'LE_tr']
            for col in gs_columns:
                if col != 'DateTime':
                    df[col] = np.nan
            
            # Merge gs results with the main dataframe
            for _, row in gs_result.iterrows():
                mask = df['DateTime'] == row['DateTime']
                if mask.any():
                    idx = df[mask].index[0]
                    for col in gs_columns:
                        if col != 'DateTime':
                            df.at[idx, col] = row[col]
            
            # Calculate Trans_dir from LE_tr (transpiration calculated from LE_tr)
            df['Trans_dir'] = (df['LE_tr'] * 1800.0) / (df['lambda'] * 1e6)  # lambda in MJ/kg
            
            # Calculate evaporation as residual
            df['Evap_resid'] = df['ET'] - df['Trans_dir']
            
            # Apply transpiration masking for precipitation events
            df = apply_transpiration_mask_advanced(df)
            
            # Set transpiration to NaN during precipitation events
            precip_mask = df['precip_mm'] > 0
            df.loc[precip_mask, 'LE_tr'] = np.nan
            df.loc[precip_mask, 'Trans_dir'] = np.nan
            df.loc[precip_mask, 'Evap_resid'] = df.loc[precip_mask, 'ET']
            
            # Also mask for a period after precipitation (e.g., 6 hours)
            post_precip_hours = 6
            for i in range(len(df)):
                if precip_mask.iloc[i]:
                    # Mask the next 'post_precip_hours' hours
                    end_idx = min(i + post_precip_hours, len(df))
                    df.iloc[i:end_idx, df.columns.get_loc('LE_tr')] = np.nan
                    df.iloc[i:end_idx, df.columns.get_loc('Trans_dir')] = np.nan
                    df.iloc[i:end_idx, df.columns.get_loc('Evap_resid')] = df.iloc[i:end_idx, df.columns.get_loc('ET')]
            
            # Remove the flag column if it exists
            if 'transpiration_flag' in df.columns:
                df = df.drop(columns=['transpiration_flag'])
                write_log(f"Removed column: transpiration_flag")
            
            # Remove duplicate columns
            df = df.loc[:, ~df.columns.duplicated()]
            
            # Select and reorder the final columns
            final_columns = [
                'DateTime', 'Year', 'DoY', 'Hour', 'Ustar', 'NETRAD_f', 'WS', 'WD', 
                'Co2', 'NEE', 'H_f', 'PAR_f', 'RH_f', 'PA_f', 'precip_mm', 'NEE_f', 
                'LE_f', 'Tair_f', 'VPD_f', 'Rg_f', 'GPP_DT', 'GPP_nt', 'Reco_DT', 
                'Reco_nt', 'ET', 'gs_iPM', 'gs_FG', 'LE_tr', 'Trans_dir', 'Evap_resid'
            ]
            
            # Keep only the columns that exist in the dataframe
            existing_columns = [col for col in final_columns if col in df.columns]
            df = df[existing_columns]
            
            # Save processed file
            output_file = os.path.join(output_folder, f"{site_name}.csv")
            df.to_csv(output_file, index=False)
            write_log(f"Saved file: {output_file}")
            
        except Exception as e:
            write_log(f"Error processing {input_file}: {str(e)}")
            continue

# Example usage
if __name__ == "__main__":
    # Define your folder paths
    input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc"
    output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated"
    
    # Process all CSV files using dry-canopy approach
    process_canopy_conductance_for_sites(input_folder, output_folder)
    
    # Process AMERI files with evaporation/transpiration split
    process_ameri_files()


#################################################################################################################

# next code to calculate canopy conductance
####

### write another function for growing season canopy conductance 