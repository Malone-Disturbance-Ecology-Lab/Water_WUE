# -*- coding: utf-8 -*-
"""
Created on Thu May  8 16:04:18 2025

@author: ammar
"""

##

############################################################################################

import pandas as pd
import os

def process_modis_site_data(input_folder, output_folder):
    """
    Separately processes ET and GPP data per site from CSV files in the input folder
    and saves one file for ET and one for GPP per site. Does NOT combine ET and GPP together.

    Parameters:
        input_folder (str): Path to the folder containing the input CSV files.
        output_folder (str): Path to the folder where output files will be saved.
    """
    os.makedirs(output_folder, exist_ok=True)

    all_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]

    # Identify unique site names
    site_names = set()
    for f in all_files:
        parts = f.split('-')
        if len(parts) >= 2:
            site_names.add(f"{parts[0]}-{parts[1]}")

    # Process ET files only
    for site in sorted(site_names):
        et_files = [f for f in all_files if f.startswith(site) and 'ET' in f]
        if len(et_files) != 3:
            print(f"Skipping ET for {site}: Expected 3 ET files, found {len(et_files)}.")
            continue

        try:
            et_df = pd.concat([
                pd.read_csv(os.path.join(input_folder, f))[['aid', 'Date', 'Mean']]
                for f in et_files
            ], ignore_index=True).rename(columns={'Mean': 'ET'})

            et_df['site_name'] = site
            et_df['DoY'] = pd.to_datetime(et_df['Date']).dt.dayofyear
            et_df = et_df.sort_values(by=['aid', 'Date'])

            et_outfile = os.path.join(output_folder, f"{site}_ET.csv")
            et_df.to_csv(et_outfile, index=False)
        except Exception as e:
            print(f"Error processing ET for {site}: {e}")

    # Process GPP files only
    for site in sorted(site_names):
        gpp_files = [f for f in all_files if f.startswith(site) and 'GPP' in f]
        if len(gpp_files) != 3:
            print(f"Skipping GPP for {site}: Expected 3 GPP files, found {len(gpp_files)}.")
            continue

        try:
            gpp_df = pd.concat([
                pd.read_csv(os.path.join(input_folder, f))[['aid', 'Date', 'Mean']]
                for f in gpp_files
            ], ignore_index=True).rename(columns={'Mean': 'GPP'})

            gpp_df['site_name'] = site
            gpp_df['DoY'] = pd.to_datetime(gpp_df['Date']).dt.dayofyear
            gpp_df = gpp_df.sort_values(by=['aid', 'Date'])

            gpp_outfile = os.path.join(output_folder, f"{site}_GPP.csv")
            gpp_df.to_csv(gpp_outfile, index=False)
        except Exception as e:
            print(f"Error processing GPP for {site}: {e}")




input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data"
output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\combo"

process_modis_site_data(input_folder, output_folder)

##############################################################################################################

import pandas as pd
import os

def modis_growing_batch(df1_folder, df2_path, save_folder):
    """
    Merge multiple MODIS ET/GPP files with average phenological growing season metrics by site.

    Parameters:
    df1_folder (str): Folder path containing multiple MODIS CSVs (each includes Date, ET, GPP, DoY, site_id)
    df2_path (str): Path to phenofit growing season CSV (includes site_name, year, avg_sos, avg_eos)
    save_folder (str): Folder path to save the merged output CSVs

    Returns:
    None
    """
    # Load phenofit growing season data
    df2 = pd.read_csv(df2_path)
    df2['site_name'] = df2['site_name'].str.strip()

    # Compute average SOS/EOS per site
    df2_avg = (
        df2[['site_name', 'avg_sos', 'avg_eos']]
        .dropna(subset=['avg_sos', 'avg_eos'])
        .groupby('site_name', as_index=False)
        .mean()
    )

    # Ensure save folder exists
    os.makedirs(save_folder, exist_ok=True)

    # Process each MODIS file
    for filename in os.listdir(df1_folder):
        if filename.endswith(".csv"):
            file_path = os.path.join(df1_folder, filename)
            df1 = pd.read_csv(file_path)

            # Standardize site name
            df1 = df1.rename(columns={'site_id': 'site_name'}) if 'site_id' in df1.columns else df1
            df1['site_name'] = df1['site_name'].str.strip()

            # Extract year from Date
            df1['Date'] = pd.to_datetime(df1['Date'], errors='coerce')
            df1['year'] = df1['Date'].dt.year

            # Merge with site-level average SOS/EOS
            merged = pd.merge(df1, df2_avg, on='site_name', how='left')

            # Save output
            output_path = os.path.join(save_folder, filename)
            merged.to_csv(output_path, index=False)
            print(f"✅ Saved merged file: {output_path}")




######################################################################################################################


# Set paths
df1_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\combo"
df2_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\phenofit_growing_season.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\modis_grow"

# Run the batch merge
modis_growing_batch(df1_folder, df2_path, save_folder)

#################################################################################################################

import pandas as pd
import os
import numpy as np
from glob import glob

def filter_growing_season_modis(input_folder, output_folder):
    """
    Filters MODIS ET and GPP files for each site to growing season using avg_sos and avg_eos.
    Keeps only years with >=2 unique months. Saves one CSV per site with ET + GPP merged.

    Parameters:
        input_folder (str): Folder with *_ET.csv and *_GPP.csv files containing avg_sos and avg_eos.
        output_folder (str): Folder to save filtered growing season data per site.
    """
    os.makedirs(output_folder, exist_ok=True)

    et_files = glob(os.path.join(input_folder, "*_ET.csv"))
    site_names = [os.path.basename(f).split('_ET.csv')[0] for f in et_files]

    for site in site_names:
        try:
            et_path = os.path.join(input_folder, f"{site}_ET.csv")
            gpp_path = os.path.join(input_folder, f"{site}_GPP.csv")

            if not os.path.exists(gpp_path):
                print(f"Missing GPP file for {site}, skipping.")
                continue

            et_df = pd.read_csv(et_path)
            gpp_df = pd.read_csv(gpp_path)

            # Merge on aid and Date
            merged = pd.merge(et_df, gpp_df[['aid', 'Date', 'GPP']], on=['aid', 'Date'], how='inner')
            merged['Date'] = pd.to_datetime(merged['Date'])
            merged['year'] = merged['Date'].dt.year
            merged['month'] = merged['Date'].dt.month

            # Filter growing season data per year
            final_list = []
            for (site_name, year), group in merged.groupby(['site_name', 'year']):
                sos = group['avg_sos'].iloc[0]
                eos = group['avg_eos'].iloc[0]

                sos_idx = (group['DoY'] - sos).abs().idxmin()
                eos_idx = (group['DoY'] - eos).abs().idxmin()

                doy_start = group.loc[sos_idx, 'DoY']
                doy_end = group.loc[eos_idx, 'DoY']

                season_df = group[(group['DoY'] >= doy_start) & (group['DoY'] <= doy_end)]

                if season_df['month'].nunique() >= 2:
                    final_list.append(season_df)

            # Save output
            if final_list:
                final_df = pd.concat(final_list, ignore_index=True)
                final_df = final_df[['aid', 'Date', 'ET', 'GPP', 'site_name', 'DoY', 'year', 'avg_sos', 'avg_eos']]
                out_path = os.path.join(output_folder, f"{site}.csv")
                final_df.to_csv(out_path, index=False)
                print(f"Saved growing season data for {site}")
            else:
                print(f"No valid growing season data for {site}, skipping.")

        except Exception as e:
            print(f"Error processing {site}: {e}")


input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\modis_grow"
output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP"

filter_growing_season_modis(input_folder, output_folder)
#################################################################################################################3

## check if WUE changing in lon term data for modis
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import theilslopes, kendalltau

def analyze_wue_trends(input_folder, output_folder):
    """
    Analyzes long-term WUE trends for all sites in the input folder.
    Processes CSV files, calculates annual WUE, computes Theil-Sen slopes and Mann-Kendall p-values,
    and saves high-resolution trend plots for each site.
    
    Args:
        input_folder (str): Path to folder containing CSV files
        output_folder (str): Path to save output plots
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Process each CSV file in the input folder
    for filename in os.listdir(input_folder):
        if not filename.endswith('.csv'):
            continue
            
        filepath = os.path.join(input_folder, filename)
        site_name = os.path.splitext(filename)[0]
        
        try:
            # Read and process data
            df = pd.read_csv(filepath)
            df['Date'] = pd.to_datetime(df['Date'])
            
            # Handle multiple AIDs by averaging values for each date
            daily_avg = df.groupby(['Date', 'site_name']).agg({
                'ET': 'mean',
                'GPP': 'mean'
            }).reset_index()
            
            # Extract year and filter date range
            daily_avg['Year'] = daily_avg['Date'].dt.year
            annual = daily_avg[(daily_avg['Year'] >= 2003) & (daily_avg['Year'] <= 2023)]
            
            # Aggregate annual sums
            annual = annual.groupby('Year', as_index=False).agg(
                ET=('ET', 'sum'),
                GPP=('GPP', 'sum')
            )
            
            # Calculate WUE and handle zero/invalid cases
            annual['WUE'] = annual['GPP'] / annual['ET']
            annual = annual.replace([np.inf, -np.inf], np.nan).dropna()
            
            if len(annual) < 2:
                print(f"Skipping {site_name}: Insufficient valid data points")
                continue
                
            # Calculate trend statistics
            years = annual['Year'].values
            wue_values = annual['WUE'].values
            
            # CORRECTED: Theil-Sen slope calculation
            slope = theilslopes(wue_values, years)[0]  # Extract first value from tuple
            
            # Mann-Kendall test
            tau, p_value = kendalltau(years, wue_values)
            
            # Generate plot
            plt.figure(figsize=(10, 6))
            bars = plt.bar(annual['Year'], annual['WUE'], color='skyblue', width=0.7)
            
            # Add value labels on bars
            for bar in bars:
                height = bar.get_height()
                plt.annotate(f'{height:.4f}',  # More decimal places for WUE values
                             xy=(bar.get_x() + bar.get_width() / 2, height),
                             xytext=(0, 3),
                             textcoords='offset points',
                             ha='center', 
                             va='bottom',
                             fontsize=8)
            
            # Configure plot appearance
            plt.title(f"Annual Water Use Efficiency: {site_name}", fontsize=14)
            plt.xlabel('Year', fontsize=12)
            plt.ylabel('WUE (gC mm⁻¹H₂O)', fontsize=12)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            
            # Add trend information
            trend_text = f"Theil-Sen Slope: {slope:.6f}\nMann-Kendall p-value: {p_value:.4f}"
            plt.figtext(0.5, 0.01, trend_text, 
                        ha='center', 
                        fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Save high-resolution plot
            output_path = os.path.join(output_folder, f"{site_name}_WUE_trend.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Processed {site_name}: slope={slope:.6f}, p={p_value:.4f}")
            
        except Exception as e:
            print(f"Error processing {site_name}: {str(e)}")

# Example usage
if __name__ == "__main__":
    input_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP\long_term_plot"
    analyze_wue_trends(input_dir, output_dir)
##################################################################################################################
#3 this code merge modis and ameri data

import os
import pandas as pd
from datetime import datetime

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import theilslopes, kendalltau

def analyze_wue_trends(input_folder, output_folder):
    """
    Analyzes long-term WUE trends for all sites in the input folder.
    Processes CSV files, calculates annual WUE, computes Theil-Sen slopes and Mann-Kendall p-values,
    and saves high-resolution trend plots for each site.
    
    Args:
        input_folder (str): Path to folder containing CSV files
        output_folder (str): Path to save output plots
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Process each CSV file in the input folder
    for filename in os.listdir(input_folder):
        if not filename.endswith('.csv'):
            continue
            
        filepath = os.path.join(input_folder, filename)
        site_name = os.path.splitext(filename)[0]
        
        try:
            # Read and process data
            df = pd.read_csv(filepath)
            df['Date'] = pd.to_datetime(df['Date'])
            
            # Multiply GPP by 1000 BEFORE any aggregation
            df['GPP'] = df['GPP'] * 1000
            
            # Handle multiple AIDs by averaging values for each date
            daily_avg = df.groupby(['Date', 'site_name']).agg({
                'ET': 'mean',
                'GPP': 'mean'
            }).reset_index()
            
            # Extract year and filter date range
            daily_avg['Year'] = daily_avg['Date'].dt.year
            annual = daily_avg[(daily_avg['Year'] >= 2003) & (daily_avg['Year'] <= 2023)]
            
            # Aggregate annual sums
            annual = annual.groupby('Year', as_index=False).agg(
                ET=('ET', 'sum'),
                GPP=('GPP', 'sum')
            )
            
            # Calculate WUE and handle zero/invalid cases
            annual['WUE'] = annual['GPP'] / annual['ET']
            annual = annual.replace([np.inf, -np.inf], np.nan).dropna()
            
            if len(annual) < 2:
                print(f"Skipping {site_name}: Insufficient valid data points")
                continue
                
            # Calculate trend statistics
            years = annual['Year'].values
            wue_values = annual['WUE'].values
            
            # Theil-Sen slope calculation
            slope, intercept, _, _ = theilslopes(wue_values, years)
            
            # Mann-Kendall test
            tau, p_value = kendalltau(years, wue_values)
            
            # Generate plot - LINE GRAPH instead of bar plot
            plt.figure(figsize=(10, 6))
            
            # Plot annual WUE as a line with markers
            plt.plot(annual['Year'], annual['WUE'], 'o-', color='blue', 
                     linewidth=2, markersize=8, label='Annual WUE')
            
            # Add value labels on points
            for year, wue in zip(annual['Year'], annual['WUE']):
                plt.annotate(f'{wue:.4f}',
                             xy=(year, wue),
                             xytext=(0, 8),
                             textcoords='offset points',
                             ha='center', 
                             va='bottom',
                             fontsize=8)
            
            # Plot trend line
            trend_line = intercept + slope * years
            plt.plot(annual['Year'], trend_line, 'r--', 
                     linewidth=2, label='Theil-Sen Trend')
            
            # Configure plot appearance
            plt.title(f"Annual Water Use Efficiency: {site_name}", fontsize=14)
            plt.xlabel('Year', fontsize=12)
            plt.ylabel('WUE (gC mm⁻¹H₂O)', fontsize=12)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.legend()
            
            # Add trend information
            trend_text = f"Theil-Sen Slope: {slope:.6f}\nMann-Kendall p-value: {p_value:.4f}"
            plt.figtext(0.5, 0.01, trend_text, 
                        ha='center', 
                        fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Save high-resolution plot
            output_path = os.path.join(output_folder, f"{site_name}_WUE_trend.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Processed {site_name}: slope={slope:.6f}, p={p_value:.4f}")
            
        except Exception as e:
            print(f"Error processing {site_name}: {str(e)}")

# Example usage
if __name__ == "__main__":
    input_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP\long_term_plot"
    analyze_wue_trends(input_dir, output_dir)
##################################################################################################################

# ameriflux and modis comparison plots


    

# ameriflux and modis comparison plots with unified legend

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from sklearn.metrics import r2_score, mean_squared_error

def plot_ameri_modis_comparison(input_folder, save_folder):
    """
    Generate scatter plots with:
    - Single unified legend showing grid colors and stats
    - Simplified grid labels (grid1, grid2, etc.)
    - grid5 in red with triangle markers
    - Other grids with circle markers
    - Increased font sizes for all elements
    """
    # Load data
    file_path = os.path.join(input_folder, "ameri_modis_monthly_gs.csv")
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find file: {file_path}")
    
    # Validate columns
    required_cols = ["site_name", "ET_obs", "ET_mod", "GPP_obs", "GPP_mod", "WUE_obs", "WUE_mod", "aid"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns: {missing_cols}")
    
    # Convert all aids to compact grid labels (grid1, grid2, etc.)
    df["aid"] = df["aid"].astype(str).str.extract('(\d+)')[0]
    df = df[df["aid"].notna()]  # Remove rows with NaN in aid
    df["aid"] = "grid" + (df["aid"].astype(int).astype(str))  # Ensure consistent numbering
    
    # Get unique sites
    sites = df["site_name"].unique()
    
    # Create plots for each variable pair
    variables = [("ET_obs", "ET_mod"), ("GPP_obs", "GPP_mod"), ("WUE_obs", "WUE_mod")]
    os.makedirs(save_folder, exist_ok=True)
    
    # Enhanced font size parameters
    font_sizes = {
        'title': 18,
        'axis_labels': 16,
        'tick_labels': 14,
        'legend_title': 14,
        'legend_text': 12,
        'metrics_text': 11
    }
    
    # Marker and color settings
    grid5_style = {
        'color': 'red',      # Red color
        'marker': '^',       # Triangle marker
        'size': 80,          # Larger size
        'label': 'grid5'     # Simplified label
    }
    other_grids_style = {
        'color': plt.cm.tab10,
        'marker': 'o',       # Circle marker
        'size': 60,
        'label_prefix': 'grid'
    }
    
    for site in sites:
        site_df = df[df["site_name"] == site]
        print(f"\nProcessing site: {site}")
        
        for x_var, y_var in variables:
            plt.figure(figsize=(12, 8))
            
            # Get unique grids and sort them properly
            grids = sorted(site_df["aid"].unique(), 
                         key=lambda x: int(x.replace('grid','')))
            
            # Store handles and labels for legend
            legend_handles = []
            legend_labels = []
            
            for grid in grids:                
                subset = site_df[(site_df["aid"] == grid) & 
                               (site_df[x_var].notna()) & 
                               (site_df[y_var].notna())]
                
                if len(subset) == 0:
                    print(f"Skipping {grid} - no valid data")
                    continue
                    
                x, y = subset[x_var], subset[y_var]
                grid_num = int(grid.replace('grid',''))
                
                # Apply different styling for grid5
                if grid == 'grid5':
                    style = grid5_style
                    color = style['color']
                    marker = style['marker']
                    size = style['size']
                else:
                    color = other_grids_style['color']((grid_num-1) % 10)
                    marker = other_grids_style['marker']
                    size = other_grids_style['size']
                
                try:
                    # Scatter plot with custom styling
                    scatter = plt.scatter(x, y, color=color, marker=marker, 
                                        s=size, alpha=0.7, edgecolors='w', linewidths=0.5)
                    
                    # Linear fit
                    coeffs = np.polyfit(x, y, 1)
                    fit_line = np.poly1d(coeffs)
                    plt.plot(x, fit_line(x), color=color, linestyle="--", lw=1.5)
                    
                    # Calculate grid metrics
                    r2 = r2_score(y, fit_line(x))
                    rmse = np.sqrt(mean_squared_error(y, fit_line(x)))
                    
                    # Create legend entry with stats
                    legend_handles.append(scatter)
                    legend_labels.append(f"{grid} (R²={r2:.2f}, RMSE={rmse:.2f})")
                        
                except Exception as e:
                    print(f"Error processing {site}, {grid}: {str(e)}")
                    continue
            
            # Calculate overall metrics
            valid_data = site_df[[x_var, y_var]].dropna()
            if len(valid_data) > 1:
                try:
                    coeffs_mean = np.polyfit(valid_data[x_var], valid_data[y_var], 1)
                    fit_line_mean = np.poly1d(coeffs_mean)
                    x_range = np.linspace(valid_data[x_var].min(), valid_data[x_var].max(), 100)
                    overall_line = plt.plot(x_range, fit_line_mean(x_range), 
                                         color="black", lw=2, label="Overall Fit")[0]
                    
                    # Calculate overall metrics
                    y_pred = fit_line_mean(valid_data[x_var])
                    overall_r2 = r2_score(valid_data[y_var], y_pred)
                    overall_rmse = np.sqrt(mean_squared_error(valid_data[y_var], y_pred))
                    
                    # Add overall to legend
                    legend_handles.append(overall_line)
                    legend_labels.append(f"Overall (R²={overall_r2:.2f}, RMSE={overall_rmse:.2f})")
                    
                except Exception as e:
                    print(f"Error calculating overall fit for {site}: {str(e)}")
            
            # Create unified legend
            if legend_handles:
                legend = plt.legend(legend_handles, legend_labels,
                                  loc='upper left',
                                  bbox_to_anchor=(1, 1),
                                  title="Grid Metrics",
                                  fontsize=font_sizes['legend_text'],
                                  handletextpad=0.5,
                                  borderaxespad=0.5)
                legend.get_title().set_fontsize(font_sizes['legend_title'])
            
            # Plot aesthetics with increased font sizes
            plt.xlabel(f"AMERI {x_var.split('_')[0]}", fontsize=font_sizes['axis_labels'])
            plt.ylabel(f"MODIS {y_var.split('_')[0]}", fontsize=font_sizes['axis_labels'])
            plt.title(f"{site}: {x_var} vs {y_var}", fontsize=font_sizes['title'])
            plt.grid(True, alpha=0.3)
            
            # Increase tick label size
            plt.xticks(fontsize=font_sizes['tick_labels'])
            plt.yticks(fontsize=font_sizes['tick_labels'])
            
            plt.tight_layout()
            
            # Save plot
            plot_name = f"{site}_{x_var}_vs_{y_var}.png".replace(":", "_")
            save_path = os.path.join(save_folder, plot_name)
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Saved plot to: {save_path}")

# Run the function
if __name__ == "__main__":
    input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
    save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\compare_ameri_modis"
    plot_ameri_modis_comparison(input_folder, save_folder)






























#############################################################################################################
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from scipy.stats import linregress
import pymannkendall as mk

def calculate_and_plot_wue(input_folder, output_folder):
    """
    Calculates yearly summed ET and GPP per aid, computes WUE, and plots WUE trend over years.
    Mann-Kendall test is used for Sen's slope and p-value. Linear regression is used to show R².
    Each year's 9 aid-points are shown in different colors. One average point per year is plotted
    in black and used for the trend line.
    """
    os.makedirs(output_folder, exist_ok=True)
    site_files = glob(os.path.join(input_folder, "*.csv"))

    for site_file in site_files:
        try:
            df = pd.read_csv(site_file)
            site_name = os.path.splitext(os.path.basename(site_file))[0]

            # Convert GPP from kg C/m² to g C/m²
            df["GPP"] = df["GPP"] * 1000

            # Convert Date to datetime and extract year
            df["Date"] = pd.to_datetime(df["Date"])
            df["year"] = df["Date"].dt.year

            # Sum ET and GPP per site, aid, year
            yearly_df = df.groupby(["site_name", "aid", "year"]).agg(
                ET_sum=("ET", "sum"),
                GPP_sum=("GPP", "sum")
            ).reset_index()

            # Calculate WUE
            yearly_df["WUE"] = yearly_df["GPP_sum"] / yearly_df["ET_sum"]

            # Calculate average WUE per site per year (1 point per year for trend line)
            avg_df = yearly_df.groupby(["site_name", "year"]).agg(
                WUE_avg=("WUE", "mean")
            ).reset_index()

            # Linear regression on average WUE
            slope, intercept, r_value, _, _ = linregress(avg_df["year"], avg_df["WUE_avg"])

            # Mann-Kendall test on yearly averages
            mk_result = mk.original_test(avg_df["WUE_avg"])
            sen_slope = mk_result.slope
            p_value = mk_result.p

            # Plotting
            plt.figure(figsize=(10, 6))
            sns.scatterplot(data=yearly_df, x="year", y="WUE", hue="aid", palette="tab10", s=60, edgecolor="black", legend=False)
            sns.scatterplot(data=avg_df, x="year", y="WUE_avg", color="black", s=80, marker="X", label="Yearly Avg WUE")
            sns.lineplot(data=avg_df, x="year", y="WUE_avg", color="black", linewidth=2, label="Trend Line")

            plt.title(f"{site_name} - WUE Trend\nSen Slope: {sen_slope:.3f}, p: {p_value:.3g}, R²: {r_value**2:.3f}")
            plt.xlabel("Year")
            plt.ylabel("WUE (g C / mm H₂O)")
            plt.tight_layout()

            # Save plot
            plot_path = os.path.join(output_folder, f"{site_name}_WUE_trend.png")
            plt.savefig(plot_path, dpi=300)
            plt.close()

        except Exception as e:
            print(f"Error processing {site_file}: {e}")


input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP"
output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\growing_season_ET_GPP\plots"

calculate_and_plot_wue(input_folder, output_folder)






















#################################################################################################################

#point data 

##################################################################################################################
import pandas as pd
import os
import glob

def modis_process(input_folder, save_folder):
    """
    Process MODIS ET and GPP data from CSV files, merge by date and site,
    and save the combined result as a single CSV file.

    Parameters:
    input_folder (str): Path to folder containing input CSV files.
    save_folder (str): Path to folder where output CSV file will be saved.

    Returns:
    modis_df (pd.DataFrame): Combined DataFrame with all sites stacked.
    """
    os.makedirs(save_folder, exist_ok=True)

    all_files = glob.glob(os.path.join(input_folder, "*.csv"))

    et_files = [f for f in all_files if "MYD16A2GF" in os.path.basename(f)]
    gpp_files = [f for f in all_files if "MYD17A2HGF-061" in os.path.basename(f)]

    all_site_data = []

    for et_file in et_files:
        base_name = os.path.basename(et_file)
        site_id = base_name.split("-")[0] + "-" + base_name.split("-")[1]

        gpp_file = next((g for g in gpp_files if site_id in os.path.basename(g)), None)
        if not gpp_file:
            print(f"No matching GPP file found for {site_id}, skipping.")
            continue

        try:
            et_df = pd.read_csv(et_file, usecols=["Date", "MYD16A2GF_061_ET_500m"])
            et_df = et_df.rename(columns={"MYD16A2GF_061_ET_500m": "ET"})
        except Exception as e:
            print(f"Error reading ET file {et_file}: {e}")
            continue

        try:
            gpp_df = pd.read_csv(gpp_file, usecols=["Date", "MYD17A2HGF_061_Gpp_500m"])
            gpp_df = gpp_df.rename(columns={"MYD17A2HGF_061_Gpp_500m": "GPP"})
        except Exception as e:
            print(f"Error reading GPP file {gpp_file}: {e}")
            continue

        merged_df = pd.merge(et_df, gpp_df, on="Date", how="outer")
        merged_df["Date"] = pd.to_datetime(merged_df["Date"], errors='coerce')
        merged_df["DoY"] = merged_df["Date"].dt.dayofyear
        merged_df["site_id"] = site_id

        all_site_data.append(merged_df)

    # Combine all site DataFrames
    modis_df = pd.concat(all_site_data, ignore_index=True)

    # Save the combined file
    output_file = os.path.join(save_folder, "modis_all_sites_combined.csv")
    modis_df.to_csv(output_file, index=False)

    print(f"Saved combined file for all sites to {output_file}")
    
    return modis_df






# Define input and output paths
input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"

# Run the function and capture the output for Spyder variable explorer
modis_df = modis_process(input_folder=input_folder, save_folder=save_folder)



#######################################################################################################
## merge modis data with phenofit sos and eos dates



import pandas as pd
import os

def modis_growing(df1_path, df2_path, save_folder):
    """
    Merge MODIS ET/GPP data with phenological growing season metrics by site and year, and save the result.

    Parameters:
    df1_path (str): Path to MODIS combined data CSV (includes Date, ET, GPP, DoY, site_id)
    df2_path (str): Path to phenofit growing season CSV (includes site_name, year, avg_sos, avg_eos)
    save_folder (str): Folder path to save the merged output CSV

    Returns:
    pd.DataFrame: Merged DataFrame (modis_pheno_fit) with site_name, avg_sos, avg_eos
    """
    # Load both datasets
    df1 = pd.read_csv(df1_path)
    df2 = pd.read_csv(df2_path)

    # Standardize and clean site name column
    df1 = df1.rename(columns={'site_id': 'site_name'}) if 'site_id' in df1.columns else df1
    df1['site_name'] = df1['site_name'].str.strip()
    df2['site_name'] = df2['site_name'].str.strip()

    # Extract year from Date
    df1['Date'] = pd.to_datetime(df1['Date'], errors='coerce')
    df1['year'] = df1['Date'].dt.year

    # Drop rows with missing avg_sos or avg_eos
    df2_valid = df2[['site_name', 'year', 'avg_sos', 'avg_eos']].dropna(subset=['avg_sos', 'avg_eos'])

    # Merge MODIS with phenofit data
    modis_pheno_fit = pd.merge(df1, df2_valid, on=['site_name', 'year'], how='left')

    # Save the result
    output_path = os.path.join(save_folder, 'modis_pheno_fit.csv')
    modis_pheno_fit.to_csv(output_path, index=False)

    print(f"✅ Merged data saved to: {output_path}")
    return modis_pheno_fit



df1_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\modis_all_sites_combined.csv"
df2_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\phenofit_growing_season.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"

modis_pheno_fit = modis_growing(df1_path, df2_path, save_folder)




change this code
for new code df1_path is path to several csv files:
    

df1_path = \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\combo

change save folder to this

save_folder =\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\modis_data\spatial_grid_data\modis_grow
#######################################################################################################

## this function merge all the available modis data with all the available salinity data, making sure
# mean salinity is also calculated based on sos and eos after merging salinity with phnofit data

# product of this code has all modis data for all years when salinity daat is present, but only for growing season


import pandas as pd
import numpy as np
import os

import pandas as pd
import numpy as np
import os

def process_modis_salinity_growing_season(modis_csv, salinity_folder, growing_season_csv, save_folder):
    # Load MODIS data
    modis_df = pd.read_csv(modis_csv)
    modis_df['Date'] = pd.to_datetime(modis_df['Date'], errors='coerce')
    modis_df['year'] = modis_df['Date'].dt.year
    modis_df['DoY'] = modis_df['Date'].dt.dayofyear
    modis_df['month'] = modis_df['Date'].dt.month
    modis_df['site_name'] = modis_df['site_name'].str.strip()

    # Convert GPP from kg C/m² to g C/m²
    modis_df['GPP'] = modis_df['GPP'] * 1000

    # Check if SOS/EOS columns exist
    if 'avg_sos' not in modis_df.columns or 'avg_eos' not in modis_df.columns:
        raise ValueError("MODIS CSV must contain 'avg_sos' and 'avg_eos' columns.")

    # Use site-level mean SOS/EOS for filtering
    gs_mean = modis_df.groupby('site_name')[['avg_sos', 'avg_eos']].mean().reset_index()
    modis_df = pd.merge(modis_df, gs_mean, on='site_name', suffixes=('', '_site_mean'), how='left')
    modis_df = modis_df.dropna(subset=['avg_sos_site_mean', 'avg_eos_site_mean'])

    # Filter MODIS data to growing season
    def filter_gs_modis(group):
        sos = int(round(group['avg_sos_site_mean'].iloc[0]))
        eos = int(round(group['avg_eos_site_mean'].iloc[0]))
        doys = group['DoY'].values
        if len(doys) == 0:
            return pd.DataFrame()
        closest_sos = doys[np.abs(doys - sos).argmin()]
        closest_eos = doys[np.abs(doys - eos).argmin()]
        return group[(group['DoY'] >= closest_sos) & (group['DoY'] <= closest_eos)]

    modis_gs = modis_df.groupby(['site_name', 'year'], group_keys=False).apply(filter_gs_modis)

    # Monthly sum of ET and GPP
    monthly_modis = modis_gs.groupby(['site_name', 'year', 'month'])[['ET', 'GPP']].sum().reset_index()
    monthly_modis['WUE'] = monthly_modis['GPP'] / monthly_modis['ET']

    # Load and process salinity data
    salinity_dfs = []
    for file in os.listdir(salinity_folder):
        if file.endswith('.csv'):
            site = file.replace('.csv', '')
            df = pd.read_csv(os.path.join(salinity_folder, file))
            df['DateTime'] = pd.to_datetime(df['DateTime'], errors='coerce')
            df['DoY'] = df['DateTime'].dt.dayofyear
            df['year'] = df['DateTime'].dt.year
            df['month'] = df['DateTime'].dt.month
            df['salinity'] = pd.to_numeric(df['salinity'], errors='coerce')
            df = df[df['salinity'] >= 0]
            df['site_name'] = site
            salinity_dfs.append(df[['site_name', 'year', 'month', 'DoY', 'salinity']])

    salinity_all = pd.concat(salinity_dfs, ignore_index=True)
    salinity_all['site_name'] = salinity_all['site_name'].str.strip()

    # Load phenofit data to get site-level avg_sos/eos for salinity filtering
    gs_info = pd.read_csv(growing_season_csv)
    gs_info['site_name'] = gs_info['site_name'].str.strip()
    gs_mean_sal = gs_info.groupby('site_name')[['avg_sos', 'avg_eos']].mean().reset_index()

    # Merge salinity with SOS/EOS
    sal_merge = pd.merge(salinity_all, gs_mean_sal, on='site_name', how='inner')

    # Filter salinity data to growing season
    def filter_gs_salinity(group):
        sos = int(round(group['avg_sos'].iloc[0]))
        eos = int(round(group['avg_eos'].iloc[0]))
        doys = group['DoY'].values
        if len(doys) == 0:
            return pd.DataFrame()
        closest_sos = doys[np.abs(doys - sos).argmin()]
        closest_eos = doys[np.abs(doys - eos).argmin()]
        return group[(group['DoY'] >= closest_sos) & (group['DoY'] <= closest_eos)]

    salinity_gs = sal_merge.groupby(['site_name', 'year'], group_keys=False).apply(filter_gs_salinity)

    # Monthly mean salinity
    monthly_salinity = salinity_gs.groupby(['site_name', 'year', 'month'])['salinity'].mean().reset_index()
    monthly_salinity.rename(columns={'salinity': 'mean_salinity'}, inplace=True)

    # Merge salinity with MODIS monthly data
    final_df = pd.merge(monthly_modis, monthly_salinity, on=['site_name', 'year', 'month'], how='inner')

    # Save output
    output_file = os.path.join(save_folder, 'modis_salinity_gs.csv')
    final_df.to_csv(output_file, index=False)

    return final_df
    

modis_sali_gs= process_modis_salinity_growing_season(
    modis_csv=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\modis_pheno_fit.csv",
    salinity_folder=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\continous_salinity\salinity",
    growing_season_csv=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\phenofit_growing_season.csv",
    save_folder=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
)


############################################################################################


import pandas as pd
import os

def modis_ameri(ameri_csv, modis_csv, save_folder):
    """
    Merge AmeriFlux and MODIS datasets with renamed MODIS variables.
    Uses 'site_name', 'year', 'month' as join keys.

    Parameters:
    ameri_csv (str): Path to AmeriFlux + salinity data CSV.
    modis_csv (str): Path to MODIS + salinity data CSV.
    save_folder (str): Folder where merged output will be saved.

    Returns:
    pd.DataFrame: Merged DataFrame.
    """
    # Load datasets
    ameri_df = pd.read_csv(ameri_csv)
    modis_df = pd.read_csv(modis_csv)

    # Standardize site name
    ameri_df['site_name'] = ameri_df['site_name'].astype(str).str.strip()
    modis_df['site_name'] = modis_df['site_name'].astype(str).str.strip()

    # Rename 'Year' to 'year' in ameri_df for consistency
    if 'Year' in ameri_df.columns:
        ameri_df.rename(columns={'Year': 'year'}, inplace=True)

    # Rename MODIS columns to avoid conflicts
    modis_df = modis_df.rename(columns={
        'ET': 'ET_mod',
        'GPP': 'GPP_mod',
        'WUE': 'WUE_mod',
        'mean_salinity': 'mean_salinity_mod'
    })

    # Merge on site_name, year, month
    final_df = pd.merge(
        ameri_df,
        modis_df,
        on=['site_name', 'year', 'month'],
        how='left'
    )

    # Save output
    output_path = os.path.join(save_folder, 'amer_modis_sal_gs.csv')
    final_df.to_csv(output_path, index=False)

    return final_df


processed_df = modis_ameri(
    ameri_csv=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\ameri_salinity_gs_monthly.csv",
    modis_csv=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\modis_salinity_gs.csv",
    save_folder=r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
)


## little bit difference in salinity because not exact sos and eos due to 8 day resolution of modis data

####################################################################################################



import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress

def plot_modis_ameri(input_path_csv, save_folder):
    df = pd.read_csv(input_path_csv)

    ylabel_map = {
        'ET': 'ET (mm)',
        'GPP': 'GPP (g C m⁻² month⁻¹)',
        'WUE': 'WUE (g C per Kg H₂O)',
    }

    target_vars = ['ET', 'GPP', 'WUE']
    mod_vars = ['ET_mod', 'GPP_mod', 'WUE_mod']

    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'H', 'X']
    palette = sns.color_palette("tab10", n_colors=len(df['site_name'].unique()))
    site_names = sorted(df['site_name'].dropna().unique())
    site_styles = {site: (markers[i % len(markers)], palette[i % len(palette)]) for i, site in enumerate(site_names)}

    # Normalize by site
    df_norm = df.copy()
    for col in target_vars + mod_vars:
        if col in df.columns:
            def min_max_norm(x):
                if x.max() == x.min():
                    return pd.Series([0] * len(x), index=x.index)
                return 2 * (x - x.min()) / (x.max() - x.min()) - 1
            df_norm[col] = df.groupby('site_name')[col].transform(min_max_norm)

    for target, mod_target in zip(target_vars, mod_vars):
        if target not in df_norm.columns or mod_target not in df_norm.columns:
            print(f"⚠️ Skipping {target} — required columns not found.")
            continue

        df_clean = df_norm.dropna(subset=['mean_salinity', target, mod_target])
        df_clean = df_clean[df_clean['mean_salinity'] >= 0]

        if df_clean.empty:
            print(f"❌ No valid data with non-negative salinity and {target} found.")
            continue

        for site in site_names:
            site_df = df_clean[df_clean['site_name'] == site]
            if len(site_df) < 2:
                continue

            x = site_df['mean_salinity'].values.reshape(-1, 1)
            y_am = site_df[target].values
            y_mod = site_df[mod_target].values

            model_am = LinearRegression().fit(x, y_am)
            model_mod = LinearRegression().fit(x, y_mod)

            y_am_pred = model_am.predict(x)
            y_mod_pred = model_mod.predict(x)

            slope_am, _, r_value_am, p_value_am, _ = linregress(site_df['mean_salinity'], site_df[target])
            slope_mod, _, r_value_mod, p_value_mod, _ = linregress(site_df['mean_salinity'], site_df[mod_target])

            fig, ax = plt.subplots(figsize=(7, 5))
            marker, color = site_styles[site]

            # Ameri data
            ax.scatter(site_df['mean_salinity'], site_df[target], label='Ameri', marker='o', s=100, color='blue')
            ax.plot(site_df['mean_salinity'], y_am_pred, linestyle='--', linewidth=2.5, color='blue')

            # MODIS data
            ax.scatter(site_df['mean_salinity'], site_df[mod_target], label='MODIS', marker='s', s=100, color='red')
            ax.plot(site_df['mean_salinity'], y_mod_pred, linestyle='--', linewidth=2.5, color='red')

            ax.set_xlabel("Salinity (ppt)", fontsize=18)
            ax.set_ylabel(ylabel_map[target], fontsize=18)

            title = (f"{site} | {target}\n"
                     f"Ameri: slope = {slope_am:.2f}, $R^2$ = {r_value_am**2:.2f}, p = {'<0.01' if p_value_am < 0.01 else f'{p_value_am:.2f}'}\n"
                     f"MODIS: slope = {slope_mod:.2f}, $R^2$ = {r_value_mod**2:.2f}, p = {'<0.01' if p_value_mod < 0.01 else f'{p_value_mod:.2f}'}")

            ax.set_title(title, fontdict={'fontsize': 15})
            ax.tick_params(axis='both', labelsize=16)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.legend(fontsize=14)
            fig.tight_layout()

            filename = f"{site}_{target}_salinity_plot.png"
            fig.savefig(os.path.join(save_folder, filename), dpi=300)
            plt.close()



########################################################################################

input_path_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\amer_modis_sal_gs.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\continous_salinity\ameri_modis"

plot_modis_ameri(input_path_csv, save_folder)






####################################################################################################
## comparison between modis data point and flux tower data point 

import math
import matplotlib.pyplot as plt

# Define coordinates
lat1, lon1 = 25.1908, -80.6391 # Point 1
lat2, lon2 = 25.3, -80.7    # Point 2




# Haversine function to calculate distance in meters
def haversine(lat1, lon1, lat2, lon2):
    R = 6371000  # Earth's radius in meters
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    
    a = math.sin(delta_phi / 2)**2 + \
        math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    return R * c

# Calculate distance
distance_m = haversine(lat1, lon1, lat2, lon2)
print(f"Distance between points: {distance_m:.2f} meters")

# Plotting
plt.figure(figsize=(6, 6))
plt.plot([lon1, lon2], [lat1, lat2], 'ro-')  # red dots connected
plt.text(lon1, lat1, ' US-TaS', verticalalignment='bottom', horizontalalignment='right')
plt.text(lon2, lat2, ' Modis Point', verticalalignment='bottom', horizontalalignment='left')
plt.title(f"Distance: {distance_m:.2f} meters")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.grid(True)
plt.show()

################################################################################################
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# Define coordinates (Florida Keys area)
lat1, lon1 = 25.1908, -80.6391  # Point 1
lat2, lon2 = 25.3, -80.7        # Point 2

# Haversine distance function (in kilometers)
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # Radius of Earth in km
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

# Compute distance
distance_km = haversine(lat1, lon1, lat2, lon2)
print(f"Distance between points: {distance_km:.4f} km")

# Plot
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Zoom out to show more of Southern Florida
ax.set_extent([-83, -79.5, 24.5, 27.5], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=1)

# Plot the two points
ax.plot([lon1, lon2], [lat1, lat2], 'ro', markersize=8, transform=ccrs.PlateCarree())
ax.text(lon1 + 0.05, lat1, 'US-TaS', fontsize=9, transform=ccrs.PlateCarree())
ax.text(lon2 + 0.05, lat2, 'Modis point', fontsize=9, transform=ccrs.PlateCarree())

# Annotate distance between points
mid_lat = (lat1 + lat2) / 2
mid_lon = (lon1 + lon2) / 2
ax.text(mid_lon, mid_lat, f'{distance_km:.2f} km', color='blue', fontsize=10, transform=ccrs.PlateCarree())

# Gridlines with labels
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

plt.title('Southern Florida with Two Points', fontsize=14)
plt.show()



import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from scipy.stats import linregress
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point

import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from scipy.stats import linregress
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point

## write a function

# first find how many years have 4 or more than 4 yrears of data 

first find mankendall slope by site for sites where data is 4 or > 4 years

include slope in final df


    # Calculate slopes
    results = []
    for site in df['site_name'].unique():
        sub = df[df['site_name'] == site]
        if len(sub) < 2:
            continue
        try:
            et_slope, _, _, _, _ = linregress(sub['ET'], sub['WUE'])
            gpp_slope, _, _, _, _ = linregress(sub['GPP'], sub['WUE'])
 
    
 also calculate 


def plot_slope_maps(csv_file, save_folder):
    # Load and filter data
    df = pd.read_csv(csv_file)
    keep_cols = ['site_name', 'ET', 'GPP', 'WUE', 'lat', 'long']
    df = df[keep_cols].dropna()

    # Calculate slopes
    results = []
    for site in df['site_name'].unique():
        sub = df[df['site_name'] == site]
        if len(sub) < 2:
            continue
        try:
            et_slope, _, _, _, _ = linregress(sub['ET'], sub['WUE'])
            gpp_slope, _, _, _, _ = linregress(sub['GPP'], sub['WUE'])
            lat = sub['lat'].iloc[0]
            lon = sub['long'].iloc[0]
            results.append({
                'site_name': site,
                'ET_slope': et_slope,
                'GPP_slope': gpp_slope,
                'lat': lat,
                'lon': lon

        
            })
        except:
            continue

    slope_df = pd.DataFrame(results)
    gpp_path = os.path.join(save_folder, 'GPP_WUE_Slope_Signs_Map.png')
    plt.tight_layout(pad=0)  # Remove extra white space around the map
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Remove the borders at the four corners
    plt.savefig(gpp_path, dpi=300)
    plt.close()

    print(f"Maps saved:\n - {et_path}\n - {gpp_path}")

# Example usage
csv_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\maps"

plot_slope_maps(csv_file, save_folder)