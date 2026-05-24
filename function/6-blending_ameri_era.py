# -*- coding: utf-8 -*-
"""
Created for blending AmeriFlux data with ERA5 data
Uses direct calibration: AmeriFlux = f(ERA5) for gap filling
Includes comprehensive diagnostic plots and NETRAD outlier removal
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings('ignore')

def blended_save(df_b, site_name, save_path):
    """
    Save the blended DataFrame as a CSV file with a dynamic name
    
    Parameters:
    - df_b: DataFrame to be saved (blended data)
    - site_name: Site identifier (e.g., 'US-Skr')
    - save_path: The path where the blended CSV file will be saved
    
    Returns:
    - None
    """
    # Create the new file name based on the site name
    output_file_name = f'gaps_blend_{site_name}.csv'
    
    # Full path to save the file
    output_file_path = os.path.join(save_path, output_file_name)
    
    # Save the blended DataFrame as CSV
    df_b.to_csv(output_file_path, index=False)
    
    print(f"  Blended data saved to: {output_file_path}")

def calculate_metrics(observed, predicted):
    """
    Calculate regression metrics
    
    Parameters:
    - observed: Array of observed values
    - predicted: Array of predicted values
    
    Returns:
    - dict: Dictionary with slope, intercept, r2, rmse
    """
    # Remove NaN values
    mask = ~(np.isnan(observed) | np.isnan(predicted))
    obs_clean = observed[mask]
    pred_clean = predicted[mask]
    
    if len(obs_clean) < 2:
        return None
    
    slope, intercept, r_value, _, _ = linregress(obs_clean, pred_clean)
    r2 = r_value ** 2
    rmse = np.sqrt(mean_squared_error(obs_clean, pred_clean))
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r2': r2,
        'rmse': rmse,
        'n': len(obs_clean)
    }

def clean_netrad(df):
    """
    Remove extreme outliers from NETRAD variable
    Values greater than 1500 are removed (set to NaN)
    No lower bound filtering, no gap filling
    
    Parameters:
    - df: DataFrame with NETRAD column
    
    Returns:
    - df: DataFrame with cleaned NETRAD
    """
    if 'NETRAD' in df.columns:
        netrad_before = df['NETRAD'].count()
        df.loc[df['NETRAD'] > 1500, 'NETRAD'] = np.nan
        netrad_after = df['NETRAD'].count()
        removed = netrad_before - netrad_after
        if removed > 0:
            print(f"    NETRAD cleaned: removed {removed} values > 1500")
    return df

def create_diagnostic_plots(df, site_name, save_path):
    """
    Create comprehensive diagnostic plots for blended data
    
    For each variable (Tair, Rg, VPD), creates:
    1. Time series plot with original, raw ERA5, calibrated ERA5, and final blended
    2. Scatter plot (AmeriFlux vs raw ERA5) with metrics
    3. Scatter plot (AmeriFlux vs calibrated ERA5) with metrics
    
    Parameters:
    - df: Blended DataFrame
    - site_name: Site identifier
    - save_path: Path to save plots
    """
    try:
        # Create plots directory
        plots_path = os.path.join(save_path, 'diagnostic_plots', site_name)
        os.makedirs(plots_path, exist_ok=True)
        
        # Define variables to plot
        variables = [
            {'name': 'Tair', 'ylabel': 'Air Temperature (°C)', 'ylim': [-30, 40]},
            {'name': 'Rg', 'ylabel': 'Shortwave Radiation (W m⁻²)', 'ylim': [-50, 1200]},
            {'name': 'VPD', 'ylabel': 'Vapor Pressure Deficit (hPa)', 'ylim': [0, 70]}
        ]
        
        for var in variables:
            var_name = var['name']
            
            # Check if required columns exist
            if var_name not in df.columns:
                print(f"    Skipping {var_name} plots - column not found")
                continue
            
            # Check if ERA5 columns exist
            era_col = f'{var_name}_era'
            era_c_col = f'{var_name}_era_c'
            
            if era_col not in df.columns:
                print(f"    Skipping {var_name} plots - {era_col} not found")
                continue
            
            # Get data
            df_plot = df.copy()
            df_plot['DateTime'] = pd.to_datetime(df_plot['DateTime'])
            
            # Create figure with 2 rows: time series (top), scatter plots (bottom)
            fig = plt.figure(figsize=(16, 12))
            fig.suptitle(f'{site_name} - {var_name} Diagnostic Plots', fontsize=16, fontweight='bold')
            
            # ========== 1. Time Series Plot (top) ==========
            ax1 = plt.subplot(2, 1, 1)
            
            # Plot original AmeriFlux data
            ax1.plot(df_plot['DateTime'], df_plot[var_name], 
                    color='blue', alpha=0.6, linewidth=0.5, label='AmeriFlux (original)')
            
            # Plot raw ERA5 data
            if df_plot[era_col].notna().any():
                ax1.plot(df_plot['DateTime'], df_plot[era_col], 
                        color='red', alpha=0.4, linewidth=0.5, label='ERA5 (raw)')
            
            # Plot calibrated ERA5 data
            if era_c_col in df_plot.columns and df_plot[era_c_col].notna().any():
                ax1.plot(df_plot['DateTime'], df_plot[era_c_col], 
                        color='orange', alpha=0.4, linewidth=0.5, label='ERA5 (calibrated)')
            
            ax1.set_ylabel(var['ylabel'])
            ax1.set_xlabel('Date')
            ax1.set_title(f'{var_name} Time Series - Comparison of Sources')
            ax1.legend(loc='upper right', fontsize=10)
            ax1.grid(True, alpha=0.3)
            
            # Set y-axis limits
            if var['ylim']:
                ax1.set_ylim(var['ylim'])
            
            # ========== 2. Scatter Plots (bottom) ==========
            # Create two subplots for scatter plots
            ax2 = plt.subplot(2, 2, 3)  # Raw ERA5 vs AmeriFlux
            ax3 = plt.subplot(2, 2, 4)  # Calibrated ERA5 vs AmeriFlux
            
            # Plot Raw ERA5 vs AmeriFlux
            mask_raw = df_plot[var_name].notna() & df_plot[era_col].notna()
            if mask_raw.sum() > 1:
                metrics_raw = calculate_metrics(
                    df_plot.loc[mask_raw, var_name].values,
                    df_plot.loc[mask_raw, era_col].values
                )
                
                ax2.scatter(df_plot.loc[mask_raw, var_name], 
                           df_plot.loc[mask_raw, era_col], 
                           alpha=0.3, s=5, color='red')
                
                # Add 1:1 line
                min_val = min(df_plot[var_name].min(), df_plot[era_col].min())
                max_val = max(df_plot[var_name].max(), df_plot[era_col].max())
                ax2.plot([min_val, max_val], [min_val, max_val], 
                        'k--', alpha=0.5, linewidth=1, label='1:1 line')
                
                # Add regression line
                if metrics_raw:
                    x_range = np.array([min_val, max_val])
                    y_pred = metrics_raw['slope'] * x_range + metrics_raw['intercept']
                    ax2.plot(x_range, y_pred, 'b-', alpha=0.7, linewidth=1.5, 
                            label=f'Fit: y={metrics_raw["slope"]:.2f}x+{metrics_raw["intercept"]:.1f}')
                    
                    # Add metrics text
                    text = f'n={metrics_raw["n"]}\nR²={metrics_raw["r2"]:.3f}\nRMSE={metrics_raw["rmse"]:.2f}'
                    ax2.text(0.05, 0.95, text, transform=ax2.transAxes, 
                            verticalalignment='top', fontsize=10,
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                
                ax2.set_xlabel(f'AmeriFlux {var_name}')
                ax2.set_ylabel(f'ERA5 Raw {var_name}')
                ax2.set_title(f'{var_name}: Raw ERA5 vs AmeriFlux')
                ax2.legend(loc='lower right', fontsize=8)
                ax2.grid(True, alpha=0.3)
            
            # Plot Calibrated ERA5 vs AmeriFlux
            if era_c_col in df_plot.columns:
                mask_cal = df_plot[var_name].notna() & df_plot[era_c_col].notna()
                if mask_cal.sum() > 1:
                    metrics_cal = calculate_metrics(
                        df_plot.loc[mask_cal, var_name].values,
                        df_plot.loc[mask_cal, era_c_col].values
                    )
                    
                    ax3.scatter(df_plot.loc[mask_cal, var_name], 
                               df_plot.loc[mask_cal, era_c_col], 
                               alpha=0.3, s=5, color='green')
                    
                    # Add 1:1 line
                    min_val = min(df_plot[var_name].min(), df_plot[era_c_col].min())
                    max_val = max(df_plot[var_name].max(), df_plot[era_c_col].max())
                    ax3.plot([min_val, max_val], [min_val, max_val], 
                            'k--', alpha=0.5, linewidth=1, label='1:1 line')
                    
                    # Add regression line
                    if metrics_cal:
                        x_range = np.array([min_val, max_val])
                        y_pred = metrics_cal['slope'] * x_range + metrics_cal['intercept']
                        ax3.plot(x_range, y_pred, 'b-', alpha=0.7, linewidth=1.5, 
                                label=f'Fit: y={metrics_cal["slope"]:.2f}x+{metrics_cal["intercept"]:.1f}')
                        
                        # Add metrics text
                        text = f'n={metrics_cal["n"]}\nR²={metrics_cal["r2"]:.3f}\nRMSE={metrics_cal["rmse"]:.2f}'
                        ax3.text(0.05, 0.95, text, transform=ax3.transAxes, 
                                verticalalignment='top', fontsize=10,
                                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                    
                    ax3.set_xlabel(f'AmeriFlux {var_name}')
                    ax3.set_ylabel(f'ERA5 Calibrated {var_name}')
                    ax3.set_title(f'{var_name}: Calibrated ERA5 vs AmeriFlux')
                    ax3.legend(loc='lower right', fontsize=8)
                    ax3.grid(True, alpha=0.3)
            
            # Adjust layout and save
            plt.tight_layout(rect=[0, 0, 1, 0.97])  # Make room for suptitle
            plt.savefig(os.path.join(plots_path, f'{site_name}_{var_name}_diagnostics.png'), dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"    Saved {var_name} diagnostic plots")
            
            # Print metrics summary
            if mask_raw.sum() > 1 and metrics_raw:
                print(f"      Raw ERA5 - R²={metrics_raw['r2']:.3f}, RMSE={metrics_raw['rmse']:.2f}, n={metrics_raw['n']}")
            if mask_cal.sum() > 1 and metrics_cal:
                print(f"      Calibrated - R²={metrics_cal['r2']:.3f}, RMSE={metrics_cal['rmse']:.2f}, n={metrics_cal['n']}")
        
        # Create a summary plot showing all three variables together (time series)
        create_summary_timeseries(df, site_name, plots_path)
        
        # Create NETRAD diagnostic plot if available
        if 'NETRAD' in df.columns:
            create_netrad_diagnostic_plot(df, site_name, plots_path)
        
    except Exception as e:
        print(f"  Warning: Could not create diagnostic plots: {e}")
        import traceback
        traceback.print_exc()

def create_summary_timeseries(df, site_name, plots_path):
    """
    Create a summary time series plot with all three variables
    
    Parameters:
    - df: Blended DataFrame
    - site_name: Site identifier
    - plots_path: Path to save plots
    """
    try:
        df_plot = df.copy()
        df_plot['DateTime'] = pd.to_datetime(df_plot['DateTime'])
        
        # Create figure with 3 subplots
        fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
        fig.suptitle(f'{site_name} - Summary Time Series (Blended Data)', fontsize=16, fontweight='bold')
        
        # Variables to plot
        vars_to_plot = [
            {'name': 'Tair', 'ylabel': 'Air Temperature (°C)', 'ylim': [-30, 40], 'color': 'blue'},
            {'name': 'Rg', 'ylabel': 'Shortwave Radiation (W m⁻²)', 'ylim': [-50, 1200], 'color': 'orange'},
            {'name': 'VPD', 'ylabel': 'Vapor Pressure Deficit (hPa)', 'ylim': [0, 70], 'color': 'green'}
        ]
        
        for idx, var in enumerate(vars_to_plot):
            var_name = var['name']
            
            if var_name not in df_plot.columns:
                axes[idx].text(0.5, 0.5, f'{var_name} data not available', 
                              ha='center', va='center', transform=axes[idx].transAxes)
                continue
            
            # Plot the blended data
            axes[idx].plot(df_plot['DateTime'], df_plot[var_name], 
                          color=var['color'], alpha=0.7, linewidth=0.5)
            
            # Add a rolling average line (optional)
            if len(df_plot) > 100:
                rolling = df_plot[var_name].rolling(window=48*7, center=True).mean()  # 7-day rolling
                axes[idx].plot(df_plot['DateTime'], rolling, 
                              color='red', alpha=0.8, linewidth=1.5, label='7-day rolling mean')
            
            axes[idx].set_ylabel(var['ylabel'])
            axes[idx].set_ylim(var['ylim'])
            axes[idx].grid(True, alpha=0.3)
            axes[idx].legend(loc='upper right', fontsize=8)
            
            # Add data completeness text
            pct = df_plot[var_name].count() / len(df_plot) * 100
            axes[idx].text(0.02, 0.95, f'Completeness: {pct:.1f}%', 
                          transform=axes[idx].transAxes, fontsize=9,
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        axes[-1].set_xlabel('Date')
        
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(os.path.join(plots_path, f'{site_name}_summary_timeseries.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"    Saved summary time series plot")
        
    except Exception as e:
        print(f"  Warning: Could not create summary timeseries: {e}")

def create_netrad_diagnostic_plot(df, site_name, plots_path):
    """
    Create a diagnostic plot for NETRAD showing before/after cleaning
    
    Parameters:
    - df: Blended DataFrame (after cleaning)
    - site_name: Site identifier
    - plots_path: Path to save plots
    """
    try:
        df_plot = df.copy()
        df_plot['DateTime'] = pd.to_datetime(df_plot['DateTime'])
        
        fig, axes = plt.subplots(2, 1, figsize=(16, 10))
        fig.suptitle(f'{site_name} - NETRAD Diagnostic', fontsize=16, fontweight='bold')
        
        # Time series plot
        ax1 = axes[0]
        ax1.plot(df_plot['DateTime'], df_plot['NETRAD'], 
                color='purple', alpha=0.6, linewidth=0.5)
        ax1.set_ylabel('NETRAD (W m⁻²)')
        ax1.set_xlabel('Date')
        ax1.set_title('NETRAD Time Series (Values > 1500 removed)')
        ax1.grid(True, alpha=0.3)
        
        # Histogram plot
        ax2 = axes[1]
        netrad_clean = df_plot['NETRAD'].dropna()
        if len(netrad_clean) > 0:
            ax2.hist(netrad_clean, bins=50, color='purple', alpha=0.7, edgecolor='black')
            ax2.set_xlabel('NETRAD (W m⁻²)')
            ax2.set_ylabel('Frequency')
            ax2.set_title('NETRAD Distribution (Values > 1500 Removed)')
            ax2.grid(True, alpha=0.3)
            
            # Add statistics
            text = f'n={len(netrad_clean)}\nMean={netrad_clean.mean():.1f}\nStd={netrad_clean.std():.1f}\nMin={netrad_clean.min():.1f}\nMax={netrad_clean.max():.1f}'
            ax2.text(0.95, 0.95, text, transform=ax2.transAxes, 
                    verticalalignment='top', horizontalalignment='right',
                    fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(os.path.join(plots_path, f'{site_name}_NETRAD_diagnostics.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"    Saved NETRAD diagnostic plot")
        
    except Exception as e:
        print(f"  Warning: Could not create NETRAD diagnostic plot: {e}")

def blending_ameri_era(df, site_name):
    """
    Blend AmeriFlux and ERA5 data with site-specific corrections
    Uses direct calibration: AmeriFlux = slope * ERA5 + intercept
    
    Parameters:
    - df: Merged DataFrame with AmeriFlux and ERA5 data
    - site_name: Site identifier (e.g., 'US-Skr')
    
    Returns:
    - df_b: Blended DataFrame
    """
    df_b = df.copy()  # Ensure the original data is not modified
    
    # Clean NETRAD outliers first (values > 1500 removed)
    df_b = clean_netrad(df_b)
    
    # List of sites that need VPD multiplied by 10
    vpd_correction_sites = ['US-Atq','US-Brw', 'US-EvM', 'US-Snd', 'US-S02','US-S03',
                            'US-S04','US-StS','US-MRf']
    
    # Apply VPD correction for specific sites (converting kPa to hPa)
    if site_name in vpd_correction_sites:
        print(f"  Applying VPD correction (multiply by 10) for {site_name}")
        if 'VPD' in df_b.columns:
            df_b['VPD'] = df_b['VPD'] * 10
    
    # Convert DateTime column
    df_b['DateTime'] = pd.to_datetime(df_b['DateTime'])
    target_months = [5, 6, 7, 8]  # Growing season months
    
    # Compute ERA-based variables
    df_b["Rg_era"] = (df_b.ssrd / 3600).round(3)
    df_b['Tair_era'] = (df_b.t2m - 273.15).round(3)
    df_b['dew_p'] = (df_b.d2m - 273.15).round(3)
    
    # Calculate VPD from ERA5 (result in hPa)
    es = 0.6108 * np.exp((17.27 * df_b['Tair_era']) / (df_b['Tair_era'] + 237.3))
    ea = 0.6108 * np.exp((17.27 * df_b['dew_p']) / (df_b['dew_p'] + 237.3))
    df_b['VPD_era'] = ((es - ea) * 10).round(3)  # Convert kPa to hPa
    
    # Handle missing Tair or Rg if they are all NaN
    if df_b['Tair'].isna().all():
        print("  All Tair values are NaN. Using Tair_era for missing Tair values.")
        df_b['Tair'] = df_b['Tair_era']
    
    if df_b['Rg'].isna().all():
        print("  All Rg values are NaN. Using Rg_era for missing Rg values.")
        df_b['Rg'] = df_b['Rg_era']
    
    # CORRECTED: Direct calibration for Tair (AmeriFlux ~ ERA5)
    valid_mask_tair = df_b['Tair'].notna() & df_b['Tair_era'].notna()
    if valid_mask_tair.sum() > 1:
        slope_tair, intercept_tair, r_value, p_value, std_err = linregress(
            df_b.loc[valid_mask_tair, 'Tair_era'],  # Predictor: ERA5
            df_b.loc[valid_mask_tair, 'Tair']       # Target: AmeriFlux
        )
        df_b['Tair_era_c'] = (slope_tair * df_b['Tair_era'] + intercept_tair).round(3)
        print(f"    Tair calibration: slope={slope_tair:.3f}, intercept={intercept_tair:.3f}, R²={r_value**2:.3f}")
    else:
        df_b['Tair_era_c'] = df_b['Tair_era']
        print(f"    Tair: Using ERA5 directly (insufficient data for calibration)")
    
    # CORRECTED: Direct calibration for Rg (AmeriFlux ~ ERA5)
    valid_mask_rg = df_b['Rg'].notna() & df_b['Rg_era'].notna()
    if valid_mask_rg.sum() > 1:
        slope_rg, intercept_rg, r_value, p_value, std_err = linregress(
            df_b.loc[valid_mask_rg, 'Rg_era'],  # Predictor: ERA5
            df_b.loc[valid_mask_rg, 'Rg']       # Target: AmeriFlux
        )
        df_b['Rg_era_c'] = (slope_rg * df_b['Rg_era'] + intercept_rg).round(3)
        print(f"    Rg calibration: slope={slope_rg:.3f}, intercept={intercept_rg:.3f}, R²={r_value**2:.3f}")
    else:
        df_b['Rg_era_c'] = df_b['Rg_era']
        print(f"    Rg: Using ERA5 directly (insufficient data for calibration)")
    
    # Fill missing values based on the blending rules
    for year in df_b['DateTime'].dt.year.unique():
        mask_year = df_b['DateTime'].dt.year == year
        mask_target_months = mask_year & df_b['DateTime'].dt.month.isin(target_months)
        
        if mask_target_months.sum() == 0:
            continue  
        
        # Calculate missing data % for target months
        missing_tair = df_b.loc[mask_target_months, 'Tair'].isna().sum() / mask_target_months.sum()
        missing_rg = df_b.loc[mask_target_months, 'Rg'].isna().sum() / mask_target_months.sum()
        
        # Fill missing values for Tair
        if missing_tair <= 0.5:
            # Use calibrated ERA5 when missingness <= 50%
            df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair'] = df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair_era_c']
        else:
            # Use raw ERA5 when missingness > 50%
            df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair'] = df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair_era']
        
        # Fill missing values for Rg
        if missing_rg <= 0.5:
            df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg'] = df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg_era_c']
        else:
            df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg'] = df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg_era']
    
    # Process VPD for each year
    first_year = df_b['DateTime'].dt.year.min()
    
    for year in df_b['DateTime'].dt.year.unique():
        mask_year = df_b['DateTime'].dt.year == year
        mask_target_months = mask_year & df_b['DateTime'].dt.month.isin(target_months)
        
        if mask_target_months.sum() == 0:
            continue  
        
        missing_vpd = df_b.loc[mask_target_months, 'VPD'].isna().sum() / mask_target_months.sum()
        
        # Use Tair and RH to calculate VPD if both are available
        mask_vpd_formula = mask_year & df_b['VPD'].isna() & df_b['Tair'].notna() & df_b['RH'].notna()
        if mask_vpd_formula.any():
            df_b.loc[mask_vpd_formula, 'VPD'] = (
                0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_formula, 'Tair']) / 
                                (df_b.loc[mask_vpd_formula, 'Tair'] + 237.3)) - 
                (df_b.loc[mask_vpd_formula, 'RH'] / 100 * 
                 0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_formula, 'Tair']) / 
                                  (df_b.loc[mask_vpd_formula, 'Tair'] + 237.3)))
            ) * 10
        
        # Handle completely missing VPD
        if df_b.loc[mask_year, 'VPD'].isna().all():
            if df_b.loc[mask_year, 'RH'].isna().all():
                df_b.loc[mask_year, 'VPD'] = df_b.loc[mask_year, 'VPD_era']
            else:
                mask_vpd_fill = mask_year & df_b['VPD'].isna() & df_b['Tair'].notna() & df_b['RH'].notna()
                if mask_vpd_fill.any():
                    df_b.loc[mask_vpd_fill, 'VPD'] = (
                        0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_fill, 'Tair']) / 
                                        (df_b.loc[mask_vpd_fill, 'Tair'] + 237.3)) - 
                        (df_b.loc[mask_vpd_fill, 'RH'] / 100 * 
                         0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_fill, 'Tair']) / 
                                          (df_b.loc[mask_vpd_fill, 'Tair'] + 237.3)))
                    ) * 10
        
        # CORRECTED: Direct calibration for VPD (AmeriFlux ~ ERA5)
        valid_mask_vpd = mask_year & df_b['VPD'].notna() & df_b['VPD_era'].notna()
        if valid_mask_vpd.sum() > 1:
            slope_vpd, intercept_vpd, r_value, p_value, std_err = linregress(
                df_b.loc[valid_mask_vpd, 'VPD_era'],  # Predictor: ERA5
                df_b.loc[valid_mask_vpd, 'VPD']       # Target: AmeriFlux
            )
            df_b.loc[mask_year, 'VPD_era_c'] = (slope_vpd * df_b.loc[mask_year, 'VPD_era'] + intercept_vpd).round(3)
            if year == first_year:
                print(f"    VPD calibration: slope={slope_vpd:.3f}, intercept={intercept_vpd:.3f}, R²={r_value**2:.3f}")
        else:
            df_b.loc[mask_year, 'VPD_era_c'] = df_b.loc[mask_year, 'VPD_era']
        
        # Final fallback based on missingness threshold
        if missing_vpd <= 0.5:
            # Use calibrated ERA5 when missingness <= 50%
            df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD'] = df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD_era_c']
        else:
            # Use raw ERA5 when missingness > 50%
            df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD'] = df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD_era']
    
    # Cleanup: Remove unrealistic VPD and Rg values
    df_b.loc[(df_b['VPD'] < 0) | (df_b['VPD'] > 70), 'VPD'] = np.nan
    df_b.loc[df_b['Rg'] < 0, 'Rg'] = 0
    
    # Return the final blended dataset
    return df_b

def process_all_blended_files(input_path, output_path, create_plots=True):
    """
    Process all merged files and create blended outputs
    
    Parameters:
    - input_path: Path to merged files (merged_*.csv)
    - output_path: Path to save blended files
    - create_plots: Boolean to create diagnostic plots
    
    Returns:
    - dict: Dictionary with site names and processed DataFrames
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)
    
    # Get all merged files
    merged_files = [f for f in os.listdir(input_path) 
                    if f.startswith('merged_') and f.endswith('.csv')]
    
    if not merged_files:
        print(f"No merged files found in {input_path}")
        return {}
    
    print(f"Found {len(merged_files)} merged files to process")
    print(f"{'='*60}")
    
    processed_results = {}
    successful = 0
    failed = 0
    
    for merged_file in merged_files:
        # Extract site name from filename (merged_US-A03.csv -> US-A03)
        site_name = merged_file.replace('merged_', '').replace('.csv', '')
        print(f"\n{'='*50}")
        print(f"Processing site: {site_name}")
        print(f"{'='*50}")
        
        try:
            # Read merged data
            merged_file_path = os.path.join(input_path, merged_file)
            df = pd.read_csv(merged_file_path)
            
            print(f"  Data loaded: {len(df)} rows, {len(df.columns)} columns")
            
            # Check if required columns exist
            required_cols = ['DateTime', 'Tair', 'Rg', 'VPD', 'RH', 'ssrd', 't2m', 'd2m']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                print(f"  Warning: Missing columns: {missing_cols}")
                # Continue anyway if essential columns exist
                if 'DateTime' not in df.columns:
                    print(f"  ✗ Error: DateTime column missing, skipping")
                    failed += 1
                    continue
            
            # Apply blending (includes NETRAD cleaning)
            df_blended = blending_ameri_era(df, site_name)
            
            # Save blended data
            blended_save(df_blended, site_name, output_path)
            
            # Create diagnostic plots if requested
            if create_plots:
                create_diagnostic_plots(df_blended, site_name, output_path)
            
            # Store in results
            processed_results[site_name] = df_blended
            successful += 1
            
            # Print data completeness summary after blending
            print(f"  Blending complete - Data completeness after blending:")
            key_vars = ['Tair', 'Rg', 'VPD', 'NETRAD']
            for var in key_vars:
                if var in df_blended.columns:
                    pct = df_blended[var].count() / len(df_blended) * 100
                    print(f"    {var}: {pct:.1f}%")
            
        except Exception as e:
            print(f"  ✗ Error processing site {site_name}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
            continue
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Processing Complete!")
    print(f"Successfully processed: {successful} sites")
    print(f"Failed: {failed} sites")
    print(f"Output saved to: {output_path}")
    print(f"{'='*60}")
    
    return processed_results

# Main execution
if __name__ == "__main__":
    # Define paths
    merged_path = r'M:\Research\WUE_CUE\ameri_data\reddy_gaps\merged_ameri_era5'
    blended_output_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps'
    
    print("="*60)
    print("STARTING BLENDING PROCESS (Direct Calibration Method)")
    print("="*60)
    print("Using direct calibration: AmeriFlux = slope * ERA5 + intercept")
    print("NETRAD cleaning: values > 1500 removed")
    print(f"Input path: {merged_path}")
    print(f"Output path: {blended_output_path}")
    print("="*60)
    
    # Process all blended files
    results = process_all_blended_files(merged_path, blended_output_path, create_plots=True)
    
    print("\n" + "="*60)
    print("BLENDING PROCESS COMPLETE!")
    print(f"Successfully processed: {len(results)} sites")
    print(f"Output saved to: {blended_output_path}")
    print("="*60)