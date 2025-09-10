# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 17:27:52 2025

@author: ammar
"""
## tested on 28 ameriflux sites
## see paths to find path to csv





##############################################################################################################

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import re
from datetime import datetime

def robust_datetime_parsing(series):
    """Robust datetime parsing with multiple format support"""
    if series.dtype == object:
        # Try multiple common formats
        formats = [
            '%m/%d/%Y %H:%M',   # "1/1/2014 0:00"
            '%Y-%m-%d %H:%M:%S', # ISO format
            '%Y%m%d%H%M',        # Compact format
            '%Y-%m-%d %H:%M',    # Without seconds
            '%m/%d/%Y %H:%M:%S', # With seconds
            '%d/%m/%Y %H:%M',    # European format
            '%d/%m/%Y %H:%M:%S'  # European with seconds
        ]
        
        for fmt in formats:
            try:
                parsed = pd.to_datetime(series, format=fmt, errors='raise')
                print(f"    Parsed DateTime with format: {fmt}")
                return parsed
            except (ValueError, TypeError):
                continue
        
        # If no format works, use flexible parsing
        print("    Using flexible datetime parsing")
        return pd.to_datetime(series, errors='coerce')
    else:
        return pd.to_datetime(series, errors='coerce')

def create_diagnostic_plots(df, site_name, output_folder):
    """
    Create diagnostic plots for key variables and save as PNG
    
    Args:
        df: DataFrame with processed data
        site_name: Site identifier
        output_folder: Folder to save the plot
    """
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Create figure with subplots
    fig, axes = plt.subplots(3, 3, figsize=(18, 15), sharex=True)
    fig.suptitle(f'{site_name} Driver Variables', fontsize=16, fontweight='bold')
    
    # Flatten axes array for easy iteration
    axes = axes.flatten()
    
    # Define variables to plot with their labels
    plot_vars = [
        ('PAR', 'PAR (μmol m⁻² s⁻¹)'),
        ('PA', 'PA (kPa)'),
        ('NETRAD', 'NETRAD (W m⁻²)'),
        ('VPD_f', 'VPD (hPa)'),
        ('Tair_f', 'Air Temperature (°C)'),
        ('LE_f', 'Latent Heat Flux (W m⁻²)'),
        ('RH', 'Relative Humidity (%)'),
        ('H', 'Sensible Heat Flux (W m⁻²)'),
        ('Co2', 'CO₂ (ppm)')
    ]
    
    # Create plots for each variable
    for i, (var, ylabel) in enumerate(plot_vars):
        ax = axes[i]
        
        # Skip if variable doesn't exist
        if var not in df.columns:
            ax.text(0.5, 0.5, f'{var} not available', 
                    ha='center', va='center', fontsize=12)
            ax.set_title(f'{var} - missing')
            continue
        
        # Filter out PA values >150 kPa for plotting
        if var == 'PA':
            plot_data = df[df['PA'] <= 150].copy()
        else:
            plot_data = df.copy()
        
        # Calculate data availability
        pct_available = plot_data[var].notna().mean() * 100
        title = f"{var} ({pct_available:.1f}% available)"
        
        # Create scatter plot with reduced marker size
        ax.scatter(plot_data['DateTime'], plot_data[var], s=1, alpha=0.7, color='blue')
        
        # Add monthly mean line
        monthly = plot_data.set_index('DateTime')[var].resample('M').mean()
        ax.plot(monthly.index, monthly, color='red', linewidth=1.5, label='Monthly Mean')
        
        # Set labels and title
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Format x-axis for the bottom row
        if i >= 6:
            ax.tick_params(axis='x', rotation=45)
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))
        else:
            ax.tick_params(labelbottom=False)
    
    # Hide any unused axes
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    # Set common xlabel on the last subplot
    axes[-1].set_xlabel('Date', fontsize=12)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Save plot
    plot_path = os.path.join(output_folder, f"{site_name}_diagnostics.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved diagnostic plot to: {plot_path}")

def process_and_merge_ameriflux_drivers(input_folder, final_save_folder, merge_folder1, merge_folder2):
    """
    Process AmeriFlux data, merge with additional sources, and save final drivers
    
    Args:
        input_folder: Folder with raw AmeriFlux gap files
        final_save_folder: Folder for final merged drivers
        merge_folder1: Folder with Reddy gaps data
        merge_folder2: Folder with AmeriFlux filled data
    """
    # Create final save folder
    os.makedirs(final_save_folder, exist_ok=True)
    
    # Get list of CSV files
    csv_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.csv')]
    if not csv_files:
        print(f"No CSV files found in: {input_folder}")
        return
    
    print(f"Found {len(csv_files)} CSV files to process")
    
    for filename in csv_files:
        file_path = os.path.join(input_folder, filename)
        site_name = filename.split('_')[1]
        print(f"\nProcessing site: {site_name} ({filename})")
        
        try:
            # Step 1: Process raw AmeriFlux data
            df = process_raw_ameriflux(file_path)
            if df is None:
                continue
            
            # Step 2: Extract variables and apply filters
            df = extract_and_filter_variables(df)
            
            # Step 3: Merge with additional data sources
            df = merge_additional_data(df, site_name, merge_folder1, merge_folder2)
            
            # Step 4: Handle PAR and RH
            df = handle_par(df)
            df = handle_rh(df)
            
            # Step 5: Apply H filter
            df = apply_h_filter(df)
            
            # Step 6: Save final data
            save_final_data(df, site_name, final_save_folder)
            
            # Step 7: Create diagnostic plots
            create_diagnostic_plots(df, site_name, final_save_folder)
            
        except Exception as e:
            print(f"  Error processing file: {e}")
            import traceback
            traceback.print_exc()
    
    # Step 8: Fix CO2 for specific sites
    print("\nFixing CO2 data for specific sites...")
    fix_co2_for_sites(final_save_folder)

def process_raw_ameriflux(file_path):
    """Process raw AmeriFlux CSV file"""
    # Header detection
    skip_rows, header_source = detect_header(file_path)
    if skip_rows is None:
        return None
    
    # Read and preprocess data
    df = read_and_preprocess(file_path, skip_rows)
    if df is None:
        return None
    
    # Create time columns
    df = create_time_columns(df)
    
    return df

def detect_header(file_path):
    """Detect header position in CSV file"""
    with open(file_path, 'r', newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = []
        for _ in range(3):
            try:
                rows.append(next(reader))
            except StopIteration:
                break
    
    skip_rows = None
    header_source = "Unknown"
    
    if len(rows) > 0:
        non_empty_count = sum(1 for cell in rows[0] if cell and cell.strip() != '')
        if non_empty_count >= 2:
            skip_rows = 0
            header_source = "First Row"
    
    if skip_rows is None and len(rows) >= 3:
        non_empty_count = sum(1 for cell in rows[2] if cell and cell.strip() != '')
        if non_empty_count >= 2:
            skip_rows = 2
            header_source = "Third Row"
    
    if skip_rows is None:
        print(f"  Error: Could not determine header row")
        return None, None
    
    print(f"  Header detected in: {header_source} (skip_rows={skip_rows})")
    return skip_rows, header_source

def read_and_preprocess(file_path, skip_rows):
    """Read and preprocess CSV file"""
    # First try reading with standard parameters
    try:
        df = pd.read_csv(file_path, skiprows=skip_rows, low_memory=False)
    except:
        # Try with different encodings if needed
        try:
            df = pd.read_csv(file_path, skiprows=skip_rows, encoding='latin1', low_memory=False)
        except Exception as e:
            print(f"  Error reading file: {e}")
            return None
    
    # Handle tab-separated files
    if df.shape[1] == 1:
        try:
            df = pd.read_csv(file_path, skiprows=skip_rows, sep='\t', low_memory=False)
        except:
            try:
                df = pd.read_csv(file_path, skiprows=skip_rows, sep='\t', encoding='latin1', low_memory=False)
            except Exception as e:
                print(f"  Error reading tab-separated file: {e}")
                return None
    
    # Parse timestamp
    timestamp_col = None
    for col in ['TIMESTAMP_START', 'TIMESTAMP', 'datetime', 'DateTime']:
        if col in df.columns:
            timestamp_col = col
            break
    
    if not timestamp_col:
        print("  No timestamp column found")
        return None
    
    # Parse timestamp column with robust method
    if timestamp_col == 'TIMESTAMP_START':
        # Handle numeric timestamp format (YYYYMMDDHHMM)
        try:
            df['TIMESTAMP'] = pd.to_datetime(
                df['TIMESTAMP_START'].astype(str), 
                format='%Y%m%d%H%M',
                errors='coerce'
            )
        except:
            df['TIMESTAMP'] = robust_datetime_parsing(df['TIMESTAMP_START'])
    else:
        df['TIMESTAMP'] = robust_datetime_parsing(df[timestamp_col])
    
    # Check for successful datetime conversion
    if df['TIMESTAMP'].isna().all():
        print("  Failed to parse datetime")
        return None
    
    # Process data
    df = df.drop_duplicates(subset='TIMESTAMP').copy()
    df.set_index('TIMESTAMP', inplace=True)
    
    # Generate complete 30-minute range
    start = df.index.min().floor('D')
    end = df.index.max().ceil('D') - pd.Timedelta(minutes=30)
    full_index = pd.date_range(start=start, end=end, freq='30min')
    
    # Reindex and fill missing values
    df = df.reindex(full_index)
    df.index.name = 'TIMESTAMP'
    df = df.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
    
    return df

def create_time_columns(df):
    """Create time-related columns"""
    df = df.reset_index()
    df['DateTime'] = df['TIMESTAMP']
    df['Year'] = df['TIMESTAMP'].dt.year
    df['DoY'] = df['TIMESTAMP'].dt.dayofyear
    df['Hour'] = df['TIMESTAMP'].dt.hour + df['TIMESTAMP'].dt.minute / 60 + 0.5
    df['Hour'] = df['Hour'].apply(lambda x: round(x * 2) / 2)
    return df

def extract_and_filter_variables(df):
    """Extract and calculate required variables with filtering"""
    # Use case-insensitive column matching
    # CO2 extraction
    co2_cols = [col for col in df.columns if re.match(r'co2', str(col), re.IGNORECASE)]
    if co2_cols:
        # Extract the first CO2 column and ensure it's numeric
        co2_data = df[co2_cols].bfill(axis=1).iloc[:, 0]
        
        # Convert to numeric, forcing any non-numeric values to NaN
        df['Co2'] = pd.to_numeric(co2_data, errors='coerce')
        
        print(f"  Using CO2 columns: {', '.join(co2_cols)}")
        print(f"  CO2 data type: {df['Co2'].dtype}")
        print(f"  CO2 non-NaN values: {df['Co2'].notna().sum()}")
        
        # Check if we have any non-numeric values that were converted to NaN
        if co2_data.dtype == 'object' and df['Co2'].isna().any():
            non_numeric_count = df['Co2'].isna().sum()
            print(f"  Warning: {non_numeric_count} non-numeric values found in CO2 data and converted to NaN")
    else:
        df['Co2'] = np.nan
        print("  No CO2 columns found")
    
    # PAR extraction
    par_cols = [col for col in df.columns if re.match(r'ppfd_in', str(col), re.IGNORECASE)]
    if par_cols:
        # Convert to numeric, handling errors
        par_data = df[par_cols].bfill(axis=1).iloc[:, 0]
        df['PAR'] = pd.to_numeric(par_data, errors='coerce')
        print(f"  Using PAR columns: {', '.join(par_cols)}")
    else:
        df['PAR'] = np.nan
        print("  No PAR columns found")
    
    # WS/WD extraction
    ws_cols = [col for col in df.columns if re.match(r'ws', str(col), re.IGNORECASE)]
    if ws_cols:
        # Convert to numeric, handling errors
        ws_data = df[ws_cols].bfill(axis=1).iloc[:, 0]
        df['WS'] = pd.to_numeric(ws_data, errors='coerce')
    else:
        df['WS'] = np.nan
        
    wd_cols = [col for col in df.columns if re.match(r'wd', str(col), re.IGNORECASE)]
    if wd_cols:
        # Convert to numeric, handling errors
        wd_data = df[wd_cols].bfill(axis=1).iloc[:, 0]
        df['WD'] = pd.to_numeric(wd_data, errors='coerce')
    else:
        df['WD'] = np.nan
    
    # NETRAD extraction and calculation
    netrad_cols = [col for col in df.columns if re.match(r'netrad|rn', str(col), re.IGNORECASE)]
    if netrad_cols:
        # Convert to numeric, handling errors
        netrad_data = df[netrad_cols].bfill(axis=1).iloc[:, 0]
        df['NETRAD'] = pd.to_numeric(netrad_data, errors='coerce')
        print(f"  Using NETRAD columns: {', '.join(netrad_cols)}")
    else:
        print("  Calculating NETRAD from components")
        # Find component columns with case-insensitive matching
        sw_in_cols = [col for col in df.columns if re.match(r'sw_in', str(col), re.IGNORECASE)]
        sw_out_cols = [col for col in df.columns if re.match(r'sw_out', str(col), re.IGNORECASE)]
        lw_in_cols = [col for col in df.columns if re.match(r'lw_in', str(col), re.IGNORECASE)]
        lw_out_cols = [col for col in df.columns if re.match(r'lw_out', str(col), re.IGNORECASE)]
        
        # Create combined columns (convert to numeric)
        if sw_in_cols:
            sw_in_data = df[sw_in_cols].bfill(axis=1).iloc[:, 0]
            df['SW_IN_combined'] = pd.to_numeric(sw_in_data, errors='coerce')
        else:
            df['SW_IN_combined'] = np.nan
        
        if sw_out_cols:
            sw_out_data = df[sw_out_cols].bfill(axis=1).iloc[:, 0]
            df['SW_OUT_combined'] = pd.to_numeric(sw_out_data, errors='coerce')
        else:
            df['SW_OUT_combined'] = np.nan
        
        if lw_in_cols:
            lw_in_data = df[lw_in_cols].bfill(axis=1).iloc[:, 0]
            df['LW_IN_combined'] = pd.to_numeric(lw_in_data, errors='coerce')
        else:
            df['LW_IN_combined'] = np.nan
        
        if lw_out_cols:
            lw_out_data = df[lw_out_cols].bfill(axis=1).iloc[:, 0]
            df['LW_OUT_combined'] = pd.to_numeric(lw_out_data, errors='coerce')
        else:
            df['LW_OUT_combined'] = np.nan
        
        # Calculate NETRAD only when all components are available
        condition = df[['SW_IN_combined', 'SW_OUT_combined', 
                        'LW_IN_combined', 'LW_OUT_combined']].notna().all(axis=1)
        df.loc[condition, 'NETRAD'] = (df['SW_IN_combined'] - df['SW_OUT_combined']) + \
                                      (df['LW_IN_combined'] - df['LW_OUT_combined'])
        
        # Clean up temporary columns
        df = df.drop(columns=['SW_IN_combined', 'SW_OUT_combined', 
                              'LW_IN_combined', 'LW_OUT_combined'])
    
    # PA extraction
    pa_cols = [col for col in df.columns if re.match(r'pa', str(col), re.IGNORECASE)]
    if pa_cols:
        # Convert to numeric, handling errors
        pa_data = df[pa_cols].bfill(axis=1).iloc[:, 0]
        df['PA'] = pd.to_numeric(pa_data, errors='coerce')
        print(f"  Using PA columns: {', '.join(pa_cols)}")
    else:
        df['PA'] = np.nan
        print("  No PA columns found")
    
    # Apply data quality filters (only if columns exist and are numeric)
    # CO2 filter
    if 'Co2' in df.columns and pd.api.types.is_numeric_dtype(df['Co2']):
        co2_count_before = df['Co2'].count()
        df['Co2'] = np.where((df['Co2'] > 1200) | (df['Co2'] < 0), np.nan, df['Co2'])
        co2_count_after = df['Co2'].count()
        print(f"  CO2 data filtered: {co2_count_before - co2_count_after} values removed")
    elif 'Co2' in df.columns:
        print(f"  CO2 column is not numeric, skipping filter")
    
    # WS filter
    if 'WS' in df.columns and pd.api.types.is_numeric_dtype(df['WS']):
        ws_count_before = df['WS'].count()
        df['WS'] = np.where((df['WS'] > 20) | (df['WS'] < 0), np.nan, df['WS'])
        ws_count_after = df['WS'].count()
        print(f"  WS data filtered: {ws_count_before - ws_count_after} values removed")
    elif 'WS' in df.columns:
        print(f"  WS column is not numeric, skipping filter")
    
    # PAR filter
    if 'PAR' in df.columns and pd.api.types.is_numeric_dtype(df['PAR']):
        par_count_before = df['PAR'].count()
        df['PAR'] = np.where(df['PAR'] < 0, np.nan, df['PAR'])
        par_count_after = df['PAR'].count()
        print(f"  PAR data filtered: {par_count_before - par_count_after} values removed")
    elif 'PAR' in df.columns:
        print(f"  PAR column is not numeric, skipping filter")
    
    return df

def merge_additional_data(df, site_name, merge_folder1, merge_folder2):
    """Merge with additional data sources with proper datetime handling"""
    # 1. First merge AmeriFlux filled data
    print(f"  Searching for AmeriFlux filled files in: {merge_folder2}")
    ameriflux_files = [f for f in os.listdir(merge_folder2) 
                      if f.lower().endswith('.csv') and site_name.lower() in f.lower()]
    
    if ameriflux_files:
        # Prefer files with "_fill" in name
        fill_files = [f for f in ameriflux_files if "_fill" in f.lower()]
        selected_file = fill_files[0] if fill_files else ameriflux_files[0]
        
        ameriflux_path = os.path.join(merge_folder2, selected_file)
        try:
            print(f"  Found AmeriFlux filled file: {selected_file}")
            
            # Read AmeriFlux filled data
            ameriflux_df = pd.read_csv(ameriflux_path, low_memory=False)
            
            # Print columns for debugging
            print(f"    Columns in AmeriFlux file: {list(ameriflux_df.columns)}")
            
            # Find datetime column
            dt_col = None
            for col in ['DateTime', 'TIMESTAMP', 'TIMESTAMP_START', 'datetime']:
                if col in ameriflux_df.columns:
                    dt_col = col
                    break
            
            if not dt_col:
                print("    No datetime column found in AmeriFlux file")
                return df
            
            # Apply robust datetime parsing
            ameriflux_df['DateTime'] = robust_datetime_parsing(ameriflux_df[dt_col])
            
            # Drop rows with invalid dates
            ameriflux_df = ameriflux_df.dropna(subset=['DateTime'])
            
            # Convert to 30-minute frequency
            ameriflux_df['DateTime'] = ameriflux_df['DateTime'].dt.round('30min')
            
            # Print sample data with actual values
            print("    Sample AmeriFlux data with values:")
            sample_cols = ['DateTime']
            for col in ['Tair_f', 'VPD_f', 'LE_f', 'Rg_f']:
                if col in ameriflux_df.columns:
                    sample_cols.append(col)
            
            sample_data = ameriflux_df[sample_cols].dropna(how='all').head(3)
            print(sample_data if not sample_data.empty else "    No valid sample data available")
            
            # Verify index alignment
            common_dates = set(df['DateTime']).intersection(set(ameriflux_df['DateTime']))
            print(f"    DateTime overlap: {len(common_dates)} timesteps")
            
            # Create a mapping from DateTime to columns
            ameriflux_map = ameriflux_df.set_index('DateTime')
            
            # Create a temporary index for the main DataFrame
            df_temp = df.set_index('DateTime')
            
            # Merge AmeriFlux columns using direct assignment
            for col in ['Tair_f', 'VPD_f', 'Rg_f', 'LE_f']:
                if col in ameriflux_map.columns:
                    # Add or update column in main DataFrame
                    if col in df_temp.columns:
                        # Fill existing NaNs with values from AmeriFlux
                        mask = df_temp[col].isna() & ameriflux_map[col].notna()
                        df_temp.loc[mask, col] = ameriflux_map.loc[mask, col]
                        merged_count = mask.sum()
                        print(f"    Merged {col} from AmeriFlux data - {merged_count} values added")
                    else:
                        # Create new column from AmeriFlux data
                        df_temp[col] = ameriflux_map[col]
                        merged_count = df_temp[col].notna().sum()
                        print(f"    Added {col} from AmeriFlux data - {merged_count} values")
                else:
                    print(f"    Column {col} not found in AmeriFlux file")
            
            # Reset index for main DataFrame
            df = df_temp.reset_index()
        except Exception as e:
            print(f"    Error loading AmeriFlux file: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"  No matching AmeriFlux filled file found for site: {site_name}")
    
    # 2. Merge Reddy gaps data
    print(f"  Searching for Reddy files in: {merge_folder1}")
    reddy_files = [f for f in os.listdir(merge_folder1) 
                  if f.lower().endswith('.csv') and site_name.lower() in f.lower() and 'gaps_blend' in f.lower()]
    
    if reddy_files:
        selected_file = reddy_files[0]
        reddy_path = os.path.join(merge_folder1, selected_file)
        try:
            print(f"  Found Reddy file: {selected_file}")
            reddy_df = pd.read_csv(reddy_path, low_memory=False)
            
            # Find datetime column
            dt_col = None
            for col in ['DateTime', 'TIMESTAMP', 'TIMESTAMP_START', 'datetime']:
                if col in reddy_df.columns:
                    dt_col = col
                    break
            
            if not dt_col:
                print("    No datetime column found in Reddy file")
                return df
            
            # Apply robust datetime parsing
            reddy_df['DateTime'] = robust_datetime_parsing(reddy_df[dt_col])
            
            # Drop rows with invalid dates
            reddy_df = reddy_df.dropna(subset=['DateTime'])
            
            # Convert to 30-minute frequency
            reddy_df['DateTime'] = reddy_df['DateTime'].dt.round('30min')
            
            # Print sample data with actual values
            print("    Sample Reddy data with values:")
            sample_cols = ['DateTime']
            for col in ['PA', 'RH', 'NETRAD', 'H']:
                if col in reddy_df.columns:
                    sample_cols.append(col)
            
            sample_data = reddy_df[sample_cols].dropna(how='all').head(3)
            print(sample_data if not sample_data.empty else "    No valid sample data available")
            
            # Verify index alignment
            common_dates = set(df['DateTime']).intersection(set(reddy_df['DateTime']))
            print(f"    DateTime overlap: {len(common_dates)} timesteps")
            
            # Create a mapping from DateTime to columns
            reddy_map = reddy_df.set_index('DateTime')
            
            # Create a temporary index for the main DataFrame
            df_temp = df.set_index('DateTime')
            
            # Merge Reddy columns
            for col in ['PA', 'RH', 'NETRAD', 'H']:
                if col in reddy_map.columns:
                    # Add or update column in main DataFrame
                    if col in df_temp.columns:
                        # Fill existing NaNs with values from Reddy
                        mask = df_temp[col].isna() & reddy_map[col].notna()
                        df_temp.loc[mask, col] = reddy_map.loc[mask, col]
                        merged_count = mask.sum()
                        print(f"    Merged {col} from Reddy data - {merged_count} values added")
                    else:
                        # Create new column from Reddy data
                        df_temp[col] = reddy_map[col]
                        merged_count = df_temp[col].notna().sum()
                        print(f"    Added {col} from Reddy data - {merged_count} values")
                else:
                    print(f"    Column {col} not found in Reddy file")
            
            # Reset index for main DataFrame
            df = df_temp.reset_index()
        except Exception as e:
            print(f"    Error loading Reddy file: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  No matching Reddy file found")
    
    return df

def handle_par(df):
    """Handle PAR data, filling with Rg_f if necessary"""
    # Conversion factor from Rg (W/m²) to PAR (μmol m⁻² s⁻¹)
    CONVERSION_FACTOR = 2.3
    
    # Create PAR column if it doesn't exist
    if 'PAR' not in df.columns:
        df['PAR'] = np.nan
    
    # Fill any missing PAR values with converted Rg_f
    if 'Rg_f' in df.columns and not df['Rg_f'].isna().all():
        # Calculate how many values we're filling
        missing_before = df['PAR'].isna().sum()
        
        # Fill NaNs in PAR with converted Rg_f
        df.loc[df['PAR'].isna(), 'PAR'] = df.loc[df['PAR'].isna(), 'Rg_f'] * CONVERSION_FACTOR
        
        # Calculate how many values we filled
        missing_after = df['PAR'].isna().sum()
        filled_count = missing_before - missing_after
        
        print(f"  Filled {filled_count} missing PAR values using Rg_f conversion")
    else:
        print("  No Rg_f available to fill missing PAR values")
    
    return df

def handle_rh(df):
    """Calculate and fill RH using Tair_f and VPD_f where possible"""
    # Check if we have the required columns for RH calculation
    if 'Tair_f' in df.columns and 'VPD_f' in df.columns:
        print("  Calculating RH from Tair_f and VPD_f where possible")
        
        # Calculate saturation vapor pressure (es) in hPa
        df['es'] = 6.1094 * np.exp((17.625 * df['Tair_f']) / (df['Tair_f'] + 243.04))
        # Calculate actual vapor pressure (ea) in hPa
        df['ea'] = df['es'] - df['VPD_f']
        
        # Calculate relative humidity (RH) in %
        df['RH_calculated'] = (df['ea'] / df['es']) * 100
        
        # Cap at 100%
        df['RH_calculated'] = np.minimum(df['RH_calculated'], 100)
        
        # Handle different cases for existing RH data
        if 'RH' in df.columns:
            # Fill gaps in existing RH with calculated values
            rh_count_before = df['RH'].count()
            df['RH'] = df['RH'].fillna(df['RH_calculated'])
            rh_count_after = df['RH'].count()
            gaps_filled = rh_count_after - rh_count_before
            print(f"    Filled {gaps_filled} gaps in RH using calculated values")
        else:
            # Create RH column from calculated values
            df['RH'] = df['RH_calculated']
            print("    Created new RH column from Tair_f and VPD_f")
        
        # Clean up temporary columns
        df = df.drop(columns=['es', 'ea', 'RH_calculated'])
    else:
        # Handle case where we can't calculate RH
        if 'RH' not in df.columns:
            df['RH'] = np.nan
            print("  No Tair_f or VPD_f available to calculate RH")
        else:
            print("  Using existing RH data (no calculation possible)")
    
    return df


def apply_h_filter(df):
    """Apply filter to H data"""
    if 'H' in df.columns:
        h_count_before = df['H'].count()
        df.loc[(df['H'] < -200) | (df['H'] > 800), 'H'] = np.nan
        h_count_after = df['H'].count()
        filtered_count = h_count_before - h_count_after
        print(f"  H data filtered: {filtered_count} values removed")
    return df

def save_final_data(df, site_name, final_save_folder):
    """Save final merged data with data validation"""
    # Define required columns
    required_columns = [
        "DateTime", "Year", "DoY", "Hour", "Co2", "PAR", "WS", "WD", 
        "NETRAD", "PA", "RH", "H", "Tair_f", "VPD_f", "Rg_f", "LE_f"
    ]
    
    # Create a new DataFrame with all required columns
    final_df = pd.DataFrame(index=range(len(df)))
    
    # Add columns if they exist in df, otherwise create as NaN
    for col in required_columns:
        if col in df.columns:
            final_df[col] = df[col]
        else:
            final_df[col] = np.nan
            print(f"  Warning: {col} not found in final dataframe - creating as NaN")
    
    # Apply PA filter: Set values >150 kPa to NaN
    if 'PA' in final_df.columns:
        pa_count_before = final_df['PA'].count()
        final_df.loc[final_df['PA'] > 150, 'PA'] = np.nan
        pa_count_after = final_df['PA'].count()
        pa_filtered = pa_count_before - pa_count_after
        print(f"  PA data filtered: {pa_filtered} values >150 kPa removed")
    
    # Add data validation report
    print("  Final data completeness (% non-missing):")
    completeness = {}
    for col in required_columns:
        pct = final_df[col].count() / len(final_df) * 100
        completeness[col] = pct
        print(f"    {col}: {pct:.1f}%")
    
    # Save to final location
    output_filename = f"{site_name}.csv"
    output_path = os.path.join(final_save_folder, output_filename)
    final_df.to_csv(output_path, index=False)
    print(f"  Saved final data to: {output_path}")
    
    # Print sample of saved data
    print("  Sample of saved data:")
    print(final_df.head(3))

def fix_co2_for_sites(filled_driver_folder):
    """Fix CO2 data for specific sites"""
    site_replacements = {'US-EDN': 'US-Myb', 'US-LA1': 'US-LA2'}
    
    for site, donor_site in site_replacements.items():
        site_file = os.path.join(filled_driver_folder, f"{site}.csv")
        donor_file = os.path.join(filled_driver_folder, f"{donor_site}.csv")
        
        if not os.path.exists(site_file) or not os.path.exists(donor_file):
            print(f"  Files not found for {site} or {donor_site}")
            continue
            
        try:
            site_df = pd.read_csv(site_file, parse_dates=['DateTime'])
            donor_df = pd.read_csv(donor_file, parse_dates=['DateTime'])
            
            # Check if CO2 is completely missing
            if 'Co2' not in site_df.columns or not site_df['Co2'].isna().all():
                print(f"  Site {site} doesn't have all NaN CO2 - skipping")
                continue
                
            if 'Co2' not in donor_df.columns:
                print(f"  Donor site {donor_site} doesn't have CO2 data")
                continue
                
            # Create mapping and merge
            co2_map = donor_df[['DateTime', 'Co2']].rename(columns={'Co2': 'Co2_donor'})
            merged_df = site_df.merge(co2_map, on='DateTime', how='left')
            merged_df['Co2'] = merged_df['Co2_donor']
            merged_df = merged_df.drop(columns=['Co2_donor'])
            
            # Save updated file
            merged_df.to_csv(site_file, index=False)
            print(f"  Replaced CO2 in {site} with data from {donor_site}")
            
            # Report completeness
            pct_coverage = merged_df['Co2'].count() / len(merged_df) * 100
            print(f"    CO2 coverage: {pct_coverage:.1f}%")
            
        except Exception as e:
            print(f"  Error fixing CO2 for {site}: {e}")

# Configuration
input_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps'
final_save_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\fill_driver'
merge_folder1 = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps2'
merge_folder2 = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_fill'

# Run the complete processing pipeline
process_and_merge_ameriflux_drivers(
    input_folder, 
    final_save_folder,
    merge_folder1,
    merge_folder2
)


###############################################################################################################

import os
import pandas as pd
import xarray as xr
import pytz
import numpy as np
from timezonefinder import TimezoneFinder
from pathlib import Path
import traceback
import gc
from collections import defaultdict
import matplotlib.pyplot as plt

def plot_pressure_rnet(processed_data, site_name):
    """
    Plot pressure and Rnet from processed ERA5 data
    
    Parameters:
        processed_data (pd.DataFrame): Processed ERA5 data
        site_name (str): Name of the site for plot title
    """
    if 'PA_era' in processed_data.columns and 'NETRAD_era' in processed_data.columns:
        print(f"Creating plots for {site_name}...")
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
        
        # Plot pressure
        ax1.plot(processed_data['TIMESTAMP'], processed_data['PA_era'], 
                color='blue', linewidth=1, alpha=0.8)
        ax1.set_ylabel('Pressure (kPa)', fontsize=12)
        ax1.set_title(f'{site_name} - Surface Pressure (PA_era)', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot Rnet
        ax2.plot(processed_data['TIMESTAMP'], processed_data['NETRAD_era'], 
                color='red', linewidth=1, alpha=0.8)
        ax2.set_ylabel('Net Radiation (W/m²)', fontsize=12)
        ax2.set_title(f'{site_name} - Net Radiation (NETRAD_era)', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Time', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.tick_params(axis='x', rotation=45)
        
        # Format x-axis for better readability
        for ax in [ax1, ax2]:
            ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M'))
        
        plt.tight_layout()
        plt.show()
        
    else:
        print(f"\n⚠️ Could not create plots for {site_name} - required variables not found")
        if 'PA_era' not in processed_data.columns:
            print("  Missing: PA_era (pressure)")
        if 'NETRAD_era' not in processed_data.columns:
            print("  Missing: NETRAD_era")

def process_era5_data(era5_path):
    """
    Processes ERA5 data with automatic variable grouping:
    1. Groups files by their variable sets
    2. Merges files within each variable group
    3. Combines different variable groups
    4. Converts UTC to local standard time
    5. Resamples to 30-minute intervals
    6. Converts pressure from Pa to kPa
    7. Calculates Rnet (net radiation) from ssr and str and converts from J/m² to W/m²
    8. Keeps only TIMESTAMP, latitude, longitude, PA_era, and NETRAD_era variables
    
    Parameters:
        era5_path (str): Path to the ERA5 data directory
    
    Returns:
        pd.DataFrame: Processed ERA5 DataFrame with only selected variables
        
    """
    # Convert to Path object
    era5_path = Path(era5_path)
    
    # Verify directory exists
    if not era5_path.is_dir():
        raise FileNotFoundError(f"ERA5 directory not found: {era5_path}")
    
    print(f"Processing ERA5 data from: {era5_path}")
    
    # Get all netCDF files
    all_files = sorted(list(era5_path.glob('*.nc')))
    print(f"Found {len(all_files)} netCDF files")
    
    if not all_files:
        raise ValueError("No netCDF files found in the specified directory")
    
    # Group files by their variable sets
    variable_groups = defaultdict(list)
    file_info = {}
    
    print("\nGrouping files by their variable sets:")
    for file in all_files:
        try:
            with xr.open_dataset(file) as ds:
                variables = tuple(sorted(ds.data_vars.keys()))
                variable_groups[variables].append(file)
                file_info[file] = {
                    'variables': variables,
                    'dims': dict(ds.dims),
                    'coords': list(ds.coords)
                }
                print(f"  {file.name}: Variables: {', '.join(variables)}")
        except Exception as e:
            print(f"⚠️ Error inspecting {file.name}: {str(e)}")
            # Add to a special group for problematic files
            variable_groups[('error',)].append(file)
    
    # Process each variable group
    group_datasets = {}
    
    for var_set, files in variable_groups.items():
        if var_set == ('error',):
            print(f"\nSkipping {len(files)} problematic files")
            continue
            
        print(f"\nProcessing variable group: {', '.join(var_set)}")
        print(f"  Contains {len(files)} files")
        
        if len(files) > 1:
            # Try to merge files in the group
            try:
                group_ds = xr.open_mfdataset(
                    files, 
                    combine='nested', 
                    concat_dim='valid_time',
                    parallel=False  # Disable parallel processing to avoid crashes
                )
                print(f"  Successfully merged {len(files)} files")
            except Exception as e:
                print(f"  Error merging files: {str(e)}")
                print("  Processing files individually...")
                ds_list = []
                for file in files:
                    try:
                        with xr.open_dataset(file) as ds:
                            ds_list.append(ds)
                    except Exception as e:
                        print(f"    ⚠️ Error opening {file.name}: {str(e)}")
                if ds_list:
                    group_ds = xr.concat(ds_list, dim='valid_time')
                else:
                    group_ds = None
        else:
            # Single file in group
            try:
                with xr.open_dataset(files[0]) as ds:
                    group_ds = ds
            except Exception as e:
                print(f"  ⚠️ Error opening {files[0].name}: {str(e)}")
                group_ds = None
        
        if group_ds is not None:
            # Create a unique group name based on variables
            group_name = "_".join(var_set)
            group_datasets[group_name] = group_ds
    
    if not group_datasets:
        raise ValueError("No valid datasets created from any files")
    
    # Combine all groups
    print("\nCombining variable groups...")
    try:
        # First try simple merge
        combined = xr.merge(group_datasets.values(), compat='override')
        print("Merged all groups using simple merge")
    except Exception as e:
        print(f"Error merging groups: {str(e)}")
        print("Combining groups with alignment...")
        # Align datasets to common time dimension
        combined = xr.align(*group_datasets.values(), join='outer')[0]
        combined = xr.merge(group_datasets.values(), compat='no_conflicts')
    
    # Find time dimension
    time_dim = 'valid_time' if 'valid_time' in combined.dims else 'time'
    print(f"Using time dimension: {time_dim}")
    
    # Convert to pandas DataFrame in chunks to save memory
    print("\nConverting to DataFrame (chunked processing)...")
    if time_dim not in combined.dims:
        raise ValueError(f"Time dimension '{time_dim}' not found in combined dataset")
    
    time_steps = combined.dims[time_dim]
    chunk_size = min(10000, time_steps)  # Number of time steps per chunk
    chunks = list(range(0, time_steps, chunk_size))
    df_chunks = []
    
    for i, start in enumerate(chunks):
        end = min(start + chunk_size, time_steps)
        print(f"Processing chunk {i+1}/{len(chunks)}: time steps {start}-{end}")
        
        # Select chunk along time dimension
        chunk = combined.isel({time_dim: slice(start, end)})
        
        # Convert chunk to DataFrame
        chunk_df = chunk.to_dataframe().reset_index()
        
        # Handle time conversion for this chunk
        time_col = time_dim
        chunk_df[time_col] = pd.to_datetime(chunk_df[time_col])
        
        # Determine timezone offset using first row coordinates
        if i == 0:  # Only need to do this once
            first_row = chunk_df.iloc[0]
            tf = TimezoneFinder()
            timezone_str = tf.timezone_at(lng=first_row['longitude'], lat=first_row['latitude'])
            
            if not timezone_str:
                # Fallback to average coordinates
                avg_lat = chunk_df['latitude'].mean()
                avg_lon = chunk_df['longitude'].mean()
                timezone_str = tf.timezone_at(lng=avg_lon, lat=avg_lat)
                if not timezone_str:
                    print("⚠️ Could not determine timezone from coordinates. Using Alaska Standard Time (UTC-9)")
                    timezone_str = 'America/Anchorage'
            
            print(f"Determined timezone: {timezone_str}")
            
            # Get UTC offset without DST
            local_tz = pytz.timezone(timezone_str)
            sample_time = pd.Timestamp('2020-01-01 00:00:00')  # Winter time to avoid DST
            offset = local_tz.utcoffset(sample_time)
            fixed_offset_hours = offset.total_seconds() / 3600
            print(f"Fixed UTC offset (without DST): {fixed_offset_hours} hours")
        
        # Apply time offset
        chunk_df['TIMESTAMP'] = chunk_df[time_col] + pd.Timedelta(hours=fixed_offset_hours)
        chunk_df.drop(columns=[time_col], inplace=True)
        
        df_chunks.append(chunk_df)
        
        # Clean up memory
        del chunk, chunk_df
        gc.collect()
    
    # Combine chunks
    print("Combining chunks...")
    df = pd.concat(df_chunks, ignore_index=True)
    
    # Clean up memory
    del combined, group_datasets
    gc.collect()
    
    # Resample to 30-minute intervals
    print("\nResampling to 30-minute intervals...")
    resampled = (df.set_index('TIMESTAMP')
                .resample('30T')
                .interpolate(method='linear')
                .reset_index())
    
    # Convert pressure from Pa to kPa
    print("\nConverting pressure from Pa to kPa...")
    if 'sp' in resampled.columns:
        resampled['PA_era'] = resampled['sp'] / 1000  # Convert Pa to kPa
        print("Pressure conversion completed")
    else:
        print("⚠️ 'sp' (surface pressure) variable not found in dataset")
    
    # Calculate Rnet (net radiation) and convert from J/m² to W/m²
    print("\nCalculating Rnet (net radiation) and converting from J/m² to W/m²...")
    if all(var in resampled.columns for var in ['ssr', 'str']):
        # ERA5 radiation variables are cumulative over the time period (J/m²)
        # Convert to average flux (W/m²) by dividing by time period in seconds
        
        # Determine time resolution from the data
        time_diff = resampled['TIMESTAMP'].diff().mean()
        time_resolution_seconds = time_diff.total_seconds()
        
        print(f"Detected time resolution: {time_resolution_seconds} seconds")
        
        # Convert from cumulative energy to average power
        ssr_wm2 = resampled['ssr'] / time_resolution_seconds  # J/m² to W/m²
        str_wm2 = resampled['str'] / time_resolution_seconds  # J/m² to W/m²
        
        # Rnet = Net shortwave radiation + Net longwave radiation
        resampled['NETRAD_era'] = ssr_wm2 + str_wm2
        
        print("Rnet calculation and unit conversion completed")
        
    else:
        missing_vars = [var for var in ['ssr', 'str'] if var not in resampled.columns]
        print(f"⚠️ Missing variables for Rnet calculation: {missing_vars}")
    
    # Select only the required variables
    print("\nSelecting only required variables: TIMESTAMP, latitude, longitude, PA_era, NETRAD_era")
    required_vars = ['TIMESTAMP', 'latitude', 'longitude', 'PA_era', 'NETRAD_era']
    
    # Check which required variables are available
    available_vars = [var for var in required_vars if var in resampled.columns]
    missing_vars = [var for var in required_vars if var not in resampled.columns]
    
    if missing_vars:
        print(f"⚠️ Missing required variables: {missing_vars}")
        print(f"Available variables: {list(resampled.columns)}")
    
    # Keep only the available required variables
    final_df = resampled[available_vars].copy()
    
    print("Processing complete!")
    print(f"Final dataset size: {final_df.shape}")
    print(f"Time range: {final_df['TIMESTAMP'].min()} to {final_df['TIMESTAMP'].max()}")
    print(f"Final variables: {list(final_df.columns)}")
    
    return final_df

def process_all_era5_folders(main_folder, output_folder):
    """
    Process all US-* subfolders in the main ERA5 directory
    
    Parameters:
        main_folder (str): Path to the main ERA5 directory
        output_folder (str): Path to save processed files
    """
    main_folder = Path(main_folder)
    output_folder = Path(output_folder)
    
    # Create output directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Find all US-* subfolders
    subfolders = sorted([f for f in main_folder.iterdir() if f.is_dir() and f.name.startswith('US-')])
    
    print(f"Found {len(subfolders)} US-* subfolders:")
    for folder in subfolders:
        print(f"  {folder.name}")
    
    # Process each subfolder
    processed_datasets = {}
    
    for subfolder in subfolders:
        try:
            print(f"\n{'='*80}")
            print(f"Processing: {subfolder.name}")
            print(f"{'='*80}")
            
            # Process the data
            processed_data = process_era5_data(subfolder)
            
            # Create plots for this site
            plot_pressure_rnet(processed_data, subfolder.name)
            
            # Save to CSV
            output_file = output_folder / f"{subfolder.name}.csv"
            processed_data.to_csv(output_file, index=False)
            print(f"Saved processed data to: {output_file}")
            
            # Store in dictionary
            processed_datasets[subfolder.name] = processed_data
            
            # Clean up memory
            del processed_data
            gc.collect()
            
        except Exception as e:
            print(f"❌ Error processing {subfolder.name}: {str(e)}")
            traceback.print_exc()
            continue
    
    print(f"\n{'='*80}")
    print("Processing complete!")
    print(f"Processed {len(processed_datasets)} out of {len(subfolders)} folders")
    print(f"Results saved to: {output_folder}")
    print(f"{'='*80}")
    
    return processed_datasets

# Main execution
if __name__ == "__main__":
    # Define paths
    main_era5_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\ERA5'
    output_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\ERA5\2-processed'
    
    # Process all folders
    all_processed_data = process_all_era5_folders(main_era5_folder, output_folder)
    
    # Print summary
    print("\nSummary of processed datasets:")
    for site_name, df in all_processed_data.items():
        print(f"{site_name}: {len(df)} rows, {len(df.columns)} columns")
        print(f"  Time range: {df['TIMESTAMP'].min()} to {df['TIMESTAMP'].max()}")
        print(f"  Variables: {list(df.columns)}")
        print()

#############################################################################################################
## next function blend and merge era data

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import linregress
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import seaborn as sns
import traceback


def driver_blend_merge_ameri_data():
    """
    Blend and merge AMERI and ERA5 data for matching sites with proper blending rules
    and control significant digits to 3 decimal places
    """
    # Define folder paths
    ameri_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\fill_driver')
    era_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\ERA5\2-processed')
    output_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\driver_ameri_era_merge_blend')
    
    # Verify folders exist
    if not ameri_folder.is_dir():
        raise FileNotFoundError(f"AMERI folder not found: {ameri_folder}")
    if not era_folder.is_dir():
        raise FileNotFoundError(f"ERA5 folder not found: {era_folder}")
    
    # Create output directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Find all CSV files
    ameri_files = sorted(list(ameri_folder.glob('US-*.csv')))
    era_files = sorted(list(era_folder.glob('US-*.csv')))
    
    print(f"Found {len(ameri_files)} AMERI files and {len(era_files)} ERA5 files")
    
    # Create dictionary of files by site name
    ameri_dict = {f.stem: f for f in ameri_files}
    era_dict = {f.stem: f for f in era_files}
    
    # Find common sites
    common_sites = set(ameri_dict.keys()) & set(era_dict.keys())
    print(f"Found {len(common_sites)} common sites to process")
    
    # Process each common site
    processed_sites = {}
    
    for site in sorted(common_sites):
        try:
            print(f"\n{'='*80}")
            print(f"Processing: {site}")
            print(f"{'='*80}")
            
            # Read AMERI data
            ameri_df = pd.read_csv(ameri_dict[site])
            print(f"AMERI data: {ameri_df.shape}")
            
            # Read ERA5 data  
            era_df = pd.read_csv(era_dict[site])
            print(f"ERA5 data: {era_df.shape}")
            
            # Convert DateTime columns
            ameri_df['DateTime'] = pd.to_datetime(ameri_df['DateTime'])
            era_df['TIMESTAMP'] = pd.to_datetime(era_df['TIMESTAMP'])
            
            # Merge data on DateTime/TIMESTAMP - KEEP TIMESTAMP FROM ERA
            df_b = pd.merge(
                ameri_df, 
                era_df[['TIMESTAMP', 'latitude', 'longitude', 'PA_era', 'NETRAD_era']],
                left_on='DateTime', 
                right_on='TIMESTAMP',
                how='left'
            )
            
            print(f"Merged data: {df_b.shape}")
            
            # Define target months (growing season)
            target_months = [5, 6, 7, 8]
            
            # NETRAD BLENDING - APPLY TO WHOLE YEAR
            print("\nProcessing NETRAD blending (applied to whole year)...")
            
            # Create copy of original NETRAD before blending
            df_b['NETRAD_original'] = df_b['NETRAD'].copy()
            
            if df_b['NETRAD'].isna().all():
                print("All NETRAD values are NaN. Using NETRAD_era for all NETRAD values.")
                df_b['NETRAD_f'] = df_b['NETRAD_era'].round(3)  # Final blended column
            else:
                # Linear regression for NETRAD bias correction
                valid_mask_netrad = df_b['NETRAD'].notna() & df_b['NETRAD_era'].notna()
                if valid_mask_netrad.sum() > 10:
                    slope_netrad, intercept_netrad, r_value, _, _ = linregress(
                        df_b.loc[valid_mask_netrad, 'NETRAD'], 
                        df_b.loc[valid_mask_netrad, 'NETRAD_era']
                    )
                    print(f"NETRAD regression: slope={slope_netrad:.3f}, intercept={intercept_netrad:.3f}, R²={r_value**2:.3f}")
                    df_b['NETRAD_era_c'] = ((df_b['NETRAD_era'] - intercept_netrad) / slope_netrad).round(3)
                else:
                    print(f"Not enough data for NETRAD regression ({valid_mask_netrad.sum()} points)")
                    df_b['NETRAD_era_c'] = df_b['NETRAD_era'].round(3)
                
                # Apply blending rules - FILL GAPS FOR WHOLE YEAR
                for year in df_b['DateTime'].dt.year.unique():
                    mask_year = df_b['DateTime'].dt.year == year
                    mask_target_months = mask_year & df_b['DateTime'].dt.month.isin(target_months)
                    
                    if mask_target_months.sum() == 0:
                        print(f"No growing season data for year {year}. Using raw ERA5 for gaps.")
                        # Fill all missing NETRAD for this WHOLE YEAR with raw ERA5 (rounded)
                        df_b.loc[mask_year & df_b['NETRAD'].isna(), 'NETRAD'] = df_b.loc[mask_year & df_b['NETRAD'].isna(), 'NETRAD_era'].round(3)
                        continue
                    
                    # Calculate missing percentage for GROWING SEASON ONLY
                    missing_netrad_growing = df_b.loc[mask_target_months, 'NETRAD'].isna().sum() / mask_target_months.sum()
                    print(f"Year {year}: {missing_netrad_growing:.1%} NETRAD missing in growing season")
                    
                    # Fill ALL missing NETRAD for this WHOLE YEAR based on growing season quality
                    if missing_netrad_growing <= 0.5:  # ≤50% missing in growing season
                        df_b.loc[mask_year & df_b['NETRAD'].isna(), 'NETRAD'] = df_b.loc[mask_year & df_b['NETRAD'].isna(), 'NETRAD_era_c']
                        print(f"  → Filled gaps with bias-corrected ERA5")
                    else:  # >50% missing in growing season
                        df_b.loc[mask_year & df_b['NETRAD'].isna(), 'NETRAD'] = df_b.loc[mask_year & df_b['NETRAD'].isna(), 'NETRAD_era'].round(3)
                        print(f"  → Filled gaps with raw ERA5")
                
                # Create final blended column for the whole year
                df_b['NETRAD_f'] = df_b['NETRAD']
            
            # PA PROCESSING (NO BLENDING - ONLY RENAME TO PA_f)
            print("\nProcessing PA (no blending - only renaming to PA_f)...")
            
            if df_b['PA'].isna().all():
                print("All PA values are NaN. Using PA_era for all PA values.")
                df_b['PA_f'] = df_b['PA_era'].round(3)  # Round to 3 decimal places
            else:
                # Just rename PA to PA_f (no blending, keep original AMERI values)
                df_b['PA_f'] = df_b['PA']
                print("Keeping original AMERI PA values as PA_f (no blending)")
            
            # Set negative PA_f values to NaN
            negative_pa_count = sum(df_b['PA_f'] < 0)
            if negative_pa_count > 0:
                df_b.loc[df_b['PA_f'] < 0, 'PA_f'] = np.nan
                print(f"Set {negative_pa_count} negative PA_f values to NaN")
            
            # REMOVE OUTLIERS FROM PA_f (5th and 95th percentiles)
            print("\nRemoving outliers from PA_f (5th and 95th percentiles)...")
            if df_b['PA_f'].notna().sum() > 0:
                # Calculate 5th and 95th percentiles
                pa_5th_percentile = df_b['PA_f'].quantile(0.05)
                pa_95th_percentile = df_b['PA_f'].quantile(0.95)
                
                print(f"PA_f 5th percentile: {pa_5th_percentile:.3f}")
                print(f"PA_f 95th percentile: {pa_95th_percentile:.3f}")
                
                # Count outliers before removal
                outliers_before = ((df_b['PA_f'] < pa_5th_percentile) | (df_b['PA_f'] > pa_95th_percentile)).sum()
                print(f"Outliers detected in PA_f: {outliers_before}")
                
                # Remove outliers by setting them to NaN
                df_b.loc[(df_b['PA_f'] < pa_5th_percentile) | (df_b['PA_f'] > pa_95th_percentile), 'PA_f'] = np.nan
                
                # Count outliers after removal
                outliers_after = ((df_b['PA_f'] < pa_5th_percentile) | (df_b['PA_f'] > pa_95th_percentile)).sum()
                print(f"Outliers remaining after removal: {outliers_after}")
                print(f"PA_f values after outlier removal: {df_b['PA_f'].notna().sum()}")
            else:
                print("No PA_f values available for outlier removal")
            
            # PLOTTING - Focus on requested plots
            print("\nGenerating regression plots...")
            
            def plot_regression(var1, var2, label1, label2, title):
                """Helper function to plot regression comparisons"""
                try:
                    # Only plot where both variables have data
                    mask = df_b[var1].notna() & df_b[var2].notna()
                    if mask.sum() == 0:
                        print(f"No overlapping data for {title}")
                        return
                    
                    r2 = r2_score(df_b.loc[mask, var1], df_b.loc[mask, var2])
                    plt.figure(figsize=(8, 6))
                    sns.scatterplot(x=df_b.loc[mask, var1], y=df_b.loc[mask, var2], alpha=0.7, label="Data points")
                    slope, intercept, _, _, _ = linregress(df_b.loc[mask, var1], df_b.loc[mask, var2])
                    plt.plot(df_b.loc[mask, var1], slope * df_b.loc[mask, var1] + intercept, 
                            color='red', label=f'Linear Fit (R²={r2:.4f})')
                    plt.xlabel(label1)
                    plt.ylabel(label2)
                    plt.title(f"{site} - {title}")
                    plt.legend()
                    plt.grid(True)
                    plt.show()
                except Exception as e:
                    print(f"Error plotting {title}: {e}")
            
            # Plot 1: NETRAD (from AMERI) vs NETRAD_era (from ERA5)
            plot_regression("NETRAD_original", "NETRAD_era", "AMERI NETRAD", "ERA5 NETRAD_era", "NETRAD vs. NETRAD_era")
            
            # Plot 2: NETRAD (from AMERI) vs NETRAD_era_c (bias-corrected ERA5)
            plot_regression("NETRAD_original", "NETRAD_era_c", "AMERI NETRAD", "NETRAD_era_c (Bias-corrected ERA5)", "NETRAD vs. NETRAD_era_c")
            
            # Time series plot of PA_f after outlier removal
            def plot_pa_timeseries():
                """Time series plot of PA_f after outlier removal"""
                try:
                    plt.figure(figsize=(14, 6))
                    
                    # Plot PA_f after outlier removal
                    plt.plot(df_b['DateTime'], df_b['PA_f'], 
                            label='PA_f (After Outlier Removal)', color='green', linewidth=1)
                    
                    plt.xlabel('DateTime')
                    plt.ylabel('Pressure (kPa)')
                    plt.title(f"{site} - PA_f Time Series (After Outlier Removal)")
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    plt.show()
                    
                    # Print PA statistics after outlier removal
                    print(f"PA_f statistics after outlier removal:")
                    print(f"  Total values: {len(df_b['PA_f'])}")
                    print(f"  Non-NaN values: {df_b['PA_f'].notna().sum()}")
                    print(f"  NaN values: {df_b['PA_f'].isna().sum()}")
                    if df_b['PA_f'].notna().sum() > 0:
                        print(f"  Min: {df_b['PA_f'].min():.3f}")
                        print(f"  Max: {df_b['PA_f'].max():.3f}")
                        print(f"  Mean: {df_b['PA_f'].mean():.3f}")
                    
                except Exception as e:
                    print(f"Error creating PA time series plot: {e}")
            
            # Create PA time series plot
            plot_pa_timeseries()
            
            # Final cleanup and precision control
            print("\nApplying precision control (3 decimal places)...")
            
            # Define which columns to round to 3 decimal places
            numeric_cols_to_round = [
                'Co2', 'PAR', 'WS', 'WD', 'NETRAD', 'PA', 'RH', 'H', 
                'Tair_f', 'VPD_f', 'Rg_f', 'LE_f', 'latitude', 'longitude',
                'NETRAD_era', 'NETRAD_era_c', 'PA_era',
                'NETRAD_f', 'PA_f'
            ]
            
            # Apply rounding only to existing numeric columns
            for col in numeric_cols_to_round:
                if col in df_b.columns:
                    df_b[col] = df_b[col].round(3)
            
            # Ensure we have all required columns including NETRAD_era_c
            required_columns = [
                'DateTime', 'Year', 'DoY', 'Hour', 'Co2', 'PAR', 'WS', 'WD', 
                'NETRAD', 'NETRAD_f', 'PA', 'PA_f', 'RH', 'H', 
                'Tair_f', 'VPD_f', 'Rg_f', 'LE_f', 'latitude', 'longitude',
                'NETRAD_era', 'NETRAD_era_c', 'PA_era'
            ]
            
            # Add any missing columns with NaN values
            for col in required_columns:
                if col not in df_b.columns:
                    df_b[col] = np.nan
                    print(f"Added missing column: {col}")
            
            # Keep all columns in the final DataFrame
            df_b = df_b[required_columns]
            
            # Save processed data with original site name
            output_file = output_folder / f"{site}.csv"
            df_b.to_csv(output_file, index=False)
            print(f"Saved processed data to: {output_file}")
            
            processed_sites[site] = df_b
            
        except Exception as e:
            print(f"❌ Error processing {site}: {str(e)}")
            traceback.print_exc()
            continue
    
    print(f"\n{'='*80}")
    print("Processing complete!")
    print(f"Processed {len(processed_sites)} out of {len(common_sites)} sites")
    print(f"Results saved to: {output_folder}")
    print(f"{'='*80}")
    
    return processed_sites

# Execute the function
if __name__ == "__main__":
    blended_data = driver_blend_merge_ameri_data()

###########################################################################################################

## bring NEE and Ustar column from blended2 data 

import pandas as pd
import numpy as np
import os
from pathlib import Path


# Alternative version with more robust filename matching
def add_ustar_nee_columns_robust():
    """
    More robust version that handles different filename patterns
    """
    # Define folder paths
    driver_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\driver_ameri_era_merge_blend')
    reddy_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps2')
    output_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\reddy_proc')
    
    # Verify folders exist
    if not driver_folder.is_dir():
        raise FileNotFoundError(f"Driver folder not found: {driver_folder}")
    if not reddy_folder.is_dir():
        raise FileNotFoundError(f"Reddy gaps folder not found: {reddy_folder}")
    
    # Create output directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Get list of CSV files
    driver_files = sorted(list(driver_folder.glob('US-*.csv')))
    reddy_files = sorted(list(reddy_folder.glob('*.csv')))
    
    print(f"Found {len(driver_files)} driver files and {len(reddy_files)} reddy files")
    
    # Process each driver file
    processed_count = 0
    for driver_file in driver_files:
        site_name = driver_file.stem  # e.g., "US-A03"
        print(f"\nProcessing site: {site_name}")
        
        try:
            # Read driver data
            driver_df = pd.read_csv(driver_file)
            print(f"  Driver data shape: {driver_df.shape}")
            
            # Find matching reddy file
            reddy_file = None
            for rf in reddy_files:
                rf_name = rf.stem
                # Check various filename patterns
                if (f"gaps_blend_{site_name}" in rf_name or 
                    f"{site_name}_gaps_blend" in rf_name or
                    site_name in rf_name):
                    reddy_file = rf
                    break
            
            if reddy_file is None:
                print(f"  No matching Reddy file found for {site_name}")
                continue
            
            print(f"  Found Reddy file: {reddy_file.name}")
            
            # Read reddy data
            reddy_df = pd.read_csv(reddy_file)
            print(f"  Reddy data shape: {reddy_df.shape}")
            
            # Convert DateTime columns
            driver_df['DateTime'] = pd.to_datetime(driver_df['DateTime'])
            reddy_df['DateTime'] = pd.to_datetime(reddy_df['DateTime'])
            
            # Check for required columns in reddy data
            required_cols = ['Ustar', 'NEE']
            available_cols = [col for col in required_cols if col in reddy_df.columns]
            
            if not available_cols:
                print(f"  None of the required columns {required_cols} found in Reddy file")
                continue
            
            print(f"  Available columns from Reddy: {available_cols}")
            
            # Create mapping from DateTime to available columns
            reddy_map = reddy_df.set_index('DateTime')[available_cols]
            
            # Create temporary index for driver data
            driver_temp = driver_df.set_index('DateTime')
            
            # Merge available columns
            for col in available_cols:
                if col in driver_temp.columns:
                    # Fill existing NaNs with values from Reddy
                    mask = driver_temp[col].isna() & reddy_map[col].notna()
                    driver_temp.loc[mask, col] = reddy_map.loc[mask, col]
                    merged_count = mask.sum()
                    print(f"  Merged {col}: {merged_count} values added")
                else:
                    # Create new column from Reddy data
                    driver_temp[col] = reddy_map[col]
                    merged_count = driver_temp[col].notna().sum()
                    print(f"  Added {col}: {merged_count} values")
            
            # Reset index
            driver_df = driver_temp.reset_index()
            
            # Save to output folder with same filename
            output_file = output_folder / driver_file.name
            driver_df.to_csv(output_file, index=False)
            print(f"  Saved to: {output_file}")
            
            processed_count += 1
            
        except Exception as e:
            print(f"  Error processing {site_name}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\nProcessing complete!")
    print(f"Processed {processed_count} out of {len(driver_files)} files")
    print(f"Results saved to: {output_folder}")

# Execute the function
if __name__ == "__main__":
    add_ustar_nee_columns_robust()
#################################################################################################################
import pandas as pd
from pathlib import Path
import os
import glob

def add_ustar_nee_columns_robust():
    """
    More robust version that handles different filename patterns and then fixes column names
    """
    # Define folder paths
    driver_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\driver_ameri_era_merge_blend')
    reddy_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps2')
    output_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\reddy_proc')
    
    # Verify folders exist
    if not driver_folder.is_dir():
        raise FileNotFoundError(f"Driver folder not found: {driver_folder}")
    if not reddy_folder.is_dir():
        raise FileNotFoundError(f"Reddy gaps folder not found: {reddy_folder}")
    
    # Create output directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Get list of CSV files
    driver_files = sorted(list(driver_folder.glob('US-*.csv')))
    reddy_files = sorted(list(reddy_folder.glob('*.csv')))
    
    print(f"Found {len(driver_files)} driver files and {len(reddy_files)} reddy files")
    
    # Process each driver file
    processed_count = 0
    for driver_file in driver_files:
        site_name = driver_file.stem  # e.g., "US-A03"
        print(f"\nProcessing site: {site_name}")
        
        try:
            # Read driver data
            driver_df = pd.read_csv(driver_file)
            print(f"  Driver data shape: {driver_df.shape}")
            
            # Find matching reddy file
            reddy_file = None
            for rf in reddy_files:
                rf_name = rf.stem
                # Check various filename patterns
                if (f"gaps_blend_{site_name}" in rf_name or 
                    f"{site_name}_gaps_blend" in rf_name or
                    site_name in rf_name):
                    reddy_file = rf
                    break
            
            if reddy_file is None:
                print(f"  No matching Reddy file found for {site_name}")
                continue
            
            print(f"  Found Reddy file: {reddy_file.name}")
            
            # Read reddy data
            reddy_df = pd.read_csv(reddy_file)
            print(f"  Reddy data shape: {reddy_df.shape}")
            
            # Convert DateTime columns
            driver_df['DateTime'] = pd.to_datetime(driver_df['DateTime'])
            reddy_df['DateTime'] = pd.to_datetime(reddy_df['DateTime'])
            
            # Check for required columns in reddy data
            required_cols = ['Ustar', 'NEE']
            available_cols = [col for col in required_cols if col in reddy_df.columns]
            
            if not available_cols:
                print(f"  None of the required columns {required_cols} found in Reddy file")
                continue
            
            print(f"  Available columns from Reddy: {available_cols}")
            
            # Create mapping from DateTime to available columns
            reddy_map = reddy_df.set_index('DateTime')[available_cols]
            
            # Create temporary index for driver data
            driver_temp = driver_df.set_index('DateTime')
            
            # Merge available columns
            for col in available_cols:
                if col in driver_temp.columns:
                    # Fill existing NaNs with values from Reddy
                    mask = driver_temp[col].isna() & reddy_map[col].notna()
                    driver_temp.loc[mask, col] = reddy_map.loc[mask, col]
                    merged_count = mask.sum()
                    print(f"  Merged {col}: {merged_count} values added")
                else:
                    # Create new column from Reddy data
                    driver_temp[col] = reddy_map[col]
                    merged_count = driver_temp[col].notna().sum()
                    print(f"  Added {col}: {merged_count} values")
            
            # Reset index
            driver_df = driver_temp.reset_index()
            
            # Save to output folder with same filename
            output_file = output_folder / driver_file.name
            driver_df.to_csv(output_file, index=False)
            print(f"  Saved to: {output_file}")
            
            processed_count += 1
            
        except Exception as e:
            print(f"  Error processing {site_name}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\nProcessing complete!")
    print(f"Processed {processed_count} out of {len(driver_files)} files")
    print(f"Results saved to: {output_folder}")
    
    # Now call the function to fix column names
    print("\n" + "="*50)
    print("Starting column name fixing process...")
    print("="*50)
    
    fix_reddy_proc_names_robust()

def fix_reddy_proc_names_robust():
    """
    Fix column names for Reddy processed CSV files with error handling.
    This function will process the files created by the main function.
    Only keep the specified columns and remove all others.
    """
    # Define folder paths
    input_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\reddy_proc')
    output_folder = Path(r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\reddy_proc\fix_names')
    
    # Verify input folder exists
    if not input_folder.is_dir():
        raise FileNotFoundError(f"Input folder not found: {input_folder}")
    
    # Create output directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Get all CSV files in input folder
    csv_files = list(input_folder.glob("*.csv"))
    
    if not csv_files:
        print("No CSV files found in input folder")
        return
    
    print(f"Found {len(csv_files)} CSV files to process")
    
    # Define EXACTLY the columns we want in the final DataFrame
    # These are the only columns that will be kept - all others will be removed
    final_columns = [
        'DateTime', 'Year', 'DoY', 'Hour', 'Co2', 'PAR', 'WS', 'WD',
        'NETRAD', 'PA', 'RH', 'H', 'Tair', 'VPD', 'Rg', 'LE',
        'latitude', 'longitude', 'Ustar', 'NEE'
    ]
    
    # Mapping from original column names to final column names
    # This removes the _f suffix from the specified columns
    column_mapping = {
        'NETRAD_f': 'NETRAD',
        'PA_f': 'PA', 
        'Tair_f': 'Tair',
        'VPD_f': 'VPD',
        'Rg_f': 'Rg',
        'LE_f': 'LE'
    }
    
    # Process each CSV file
    processed_count = 0
    for csv_file in csv_files:
        try:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            original_columns = list(df.columns)
            print(f"\nProcessing: {csv_file.name}")
            print(f"  Original columns: {len(original_columns)}")
            
            # First, rename columns that need renaming (remove _f suffix)
            df = df.rename(columns=column_mapping)
            
            # Get available columns that match our final desired columns
            available_columns = [col for col in final_columns if col in df.columns]
            missing_columns = [col for col in final_columns if col not in df.columns]
            extra_columns = [col for col in df.columns if col not in final_columns]
            
            print(f"  Available desired columns: {len(available_columns)}/{len(final_columns)}")
            print(f"  Missing columns: {missing_columns}")
            print(f"  Extra columns that will be removed: {extra_columns}")
            
            # Keep ONLY the available desired columns (remove all extra columns)
            df = df[available_columns]
            
            # Add missing columns as NaN (with proper data type where possible)
            for col in missing_columns:
                if col in ['Year', 'DoY', 'Hour']:
                    df[col] = pd.NA  # Integer columns
                else:
                    df[col] = pd.NA  # Float columns
            
            # Reorder columns to match final desired order
            df = df[final_columns]
            
            # Get the filename for output
            output_path = output_folder / csv_file.name
            
            # Save the processed dataframe
            df.to_csv(output_path, index=False)
            
            print(f"  Final shape: {df.shape}")
            print(f"  Saved to: {output_path.name}")
            
            processed_count += 1
            
        except Exception as e:
            print(f"Error processing {csv_file.name}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    print(f"\nColumn fixing complete!")
    print(f"Processed {processed_count} out of {len(csv_files)} files")
    print(f"Results saved to: {output_folder}")
    
    # Show final column structure
    if processed_count > 0:
        print(f"\nFinal DataFrame columns ({len(final_columns)} total):")
        for i, col in enumerate(final_columns, 1):
            print(f"  {i:2d}. {col}")

# Execute the function
if __name__ == "__main__":
    add_ustar_nee_columns_robust()



















###################################################################################################################
## check if required columns are present 

import os
import csv

def main():
    input_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps'
    
    if not os.path.isdir(input_folder):
        print(f"Directory not found: {input_folder}")
        return
    
    for filename in os.listdir(input_folder):
        if filename.lower().endswith('.csv'):
            file_path = os.path.join(input_folder, filename)
            try:
                with open(file_path, 'r', newline='', encoding='utf-8-sig') as csvfile:
                    reader = csv.reader(csvfile)
                    rows = []
                    for _ in range(3):  # Read up to first 3 rows
                        try:
                            rows.append(next(reader))
                        except StopIteration:
                            break
                
                header = None
                header_source = "Unknown"
                
                # Check first row for header
                if len(rows) > 0:
                    non_empty_count = sum(1 for cell in rows[0] if cell.strip() != '')
                    if non_empty_count >= 2:
                        header = [cell.strip() for cell in rows[0]]
                        header_source = "First Row"
                
                # Check third row if first row not valid
                if header is None and len(rows) >= 3:
                    non_empty_count = sum(1 for cell in rows[2] if cell.strip() != '')
                    if non_empty_count >= 2:
                        header = [cell.strip() for cell in rows[2]]
                        header_source = "Third Row"
                
                if header is None:
                    print(f"File: {filename}")
                    print("  Error: Could not determine header row (no suitable row with ≥2 non-empty columns)")
                    print()
                    continue
                
                # Define column patterns to check
                patterns = {
                    "PPFD_IN": "ppfd_in",
                    "CO2": "co2",
                    "WD": "wd",
                    "WS": "ws"
                }
                
                results = {}
                variations = {}
                
                # Check each pattern
                for col_name, pattern in patterns.items():
                    col_variations = [h for h in header if h.lower().startswith(pattern)]
                    results[col_name] = bool(col_variations)
                    variations[col_name] = col_variations
                
                print(f"File: {filename}")
                print(f"  Header source: {header_source}")
                
                # Print results for all columns
                all_present = True
                for col_name in patterns:
                    present = results[col_name]
                    print(f"  {col_name}: {'Present' if present else 'Not Present'}", end='')
                    if present:
                        print(f" (Variations: {', '.join(variations[col_name])})")
                    else:
                        print()
                    all_present = all_present and present
                
                # Print all columns if any required column is missing
                if not all_present:
                    print("  Available columns:")
                    for col in header:
                        print(f"    - {col}")
                print()  # Blank line between files
                
            except Exception as e:
                print(f"{filename}: Error processing file - {str(e)}")

if __name__ == "__main__":
    main()


############################################################################################################

# in the next column following columns are filled but _f not used because reddy proc does not accept
# columns with _f. soe filled columns are following
#LE, VPD, Rg, Tair


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import matplotlib

# Use a non-interactive backend to avoid memory issues with large figures
matplotlib.use('Agg')

def driver_growing_season():
    # Define paths
    input_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\fill_driver"
    wue_cue_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
    save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\fill_driver\growing_season"
    
    # Create save folder if it doesn't exist
    os.makedirs(save_folder, exist_ok=True)
    
    # Read growing season data
    wue_cue_df = pd.read_csv(wue_cue_path)
    
    # Define variables to plot - explicitly including CO2
    target_vars = ['CO2', 'PAR', 'WS', 'WD', 'NETRAD', 'PA', 'RH', 'H', 
                   'Tair_f', 'VPD_f', 'Rg_f', 'LE_f']
    
    # Process each driver CSV file
    for file in os.listdir(input_folder):
        if not file.endswith('.csv'):
            continue
            
        file_path = os.path.join(input_folder, file)
        site = os.path.splitext(file)[0]
        print(f"Processing site: {site}")
        
        try:
            # Read driver data
            driver_df = pd.read_csv(file_path, parse_dates=['DateTime'])
            driver_df['year'] = driver_df['DateTime'].dt.year
            driver_df['DOY'] = driver_df['DateTime'].dt.dayofyear + driver_df['DateTime'].dt.hour / 24.0
            
            # Get SOS/EOS for this site
            site_years = wue_cue_df[wue_cue_df['site_name'] == site]
            if site_years.empty:
                print(f"  Skipping: no growing season data found")
                continue
                
            # Filter data for growing season
            growing_data = []
            for year, group in driver_df.groupby('year'):
                year_data = site_years[site_years['year'] == year]
                if year_data.empty:
                    continue
                sos = year_data['avg_sos'].iloc[0]
                eos = year_data['avg_eos'].iloc[0]
                year_growing = group[(group['DOY'] >= sos) & (group['DOY'] <= eos)].copy()
                
                # Add season markers
                year_growing['avg_sos'] = sos
                year_growing['avg_eos'] = eos
                
                growing_data.append(year_growing)
            
            if not growing_data:
                print(f"  No growing season data found for available years")
                continue
            growing_df = pd.concat(growing_data)
            
            # Data cleaning
            if 'RH' in growing_df:
                growing_df.loc[(growing_df['RH'] < 0) | (growing_df['RH'] > 100), 'RH'] = np.nan
            if 'PA' in growing_df:
                growing_df.loc[growing_df['PA'] > 120, 'PA'] = np.nan
            if 'VPD_f' in growing_df:
                growing_df.loc[growing_df['VPD_f'] < 0, 'VPD_f'] = np.nan
            
            # 1. Data availability plot - LARGE SIZE (18x14 inches)
            availability = {}
            # Check CO2 explicitly
            co2_col = next((col for col in growing_df.columns if col.upper() in ['CO2', 'CO_2']), 'CO2')
            for var in target_vars:
                # Handle CO2 column name variations
                data_var = co2_col if var == 'CO2' else var
                
                if data_var in growing_df:
                    avail_pct = (1 - growing_df[data_var].isna().mean()) * 100
                    availability[var] = avail_pct
                else:
                    print(f"  Warning: {var} not found in data")
                    availability[var] = 0
            
            # Create large figure for data availability
            plt.figure(figsize=(18, 14))
            sorted_vars = sorted(availability.items(), key=lambda x: x[1])
            bars = plt.barh([v[0] for v in sorted_vars], [v[1] for v in sorted_vars])
            plt.xlabel('Data Availability (%)', fontsize=18)
            plt.ylabel('Variable', fontsize=18)
            plt.title(f'Growing Season Data Availability: {site}', fontsize=22)
            plt.xlim(0, 100)
            
            # Add percentage labels to bars
            for bar in bars:
                width = bar.get_width()
                plt.text(width + 1, bar.get_y() + bar.get_height()/2, 
                         f'{width:.1f}%', 
                         ha='left', va='center', fontsize=14)
            
            plt.tick_params(axis='both', which='major', labelsize=16)
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            avail_plot_path = os.path.join(save_folder, f"{site}_growing_season_availability.png")
            plt.savefig(avail_plot_path, dpi=120)
            plt.close()
            
            # 2. Time series plots - OPTIMIZED 3-PANEL LAYOUT WITH YEAR-ONLY X-AXIS
            # Get list of variables actually present in the data
            present_vars = []
            for var in target_vars:
                data_var = co2_col if var == 'CO2' else var
                if data_var in growing_df:
                    present_vars.append((var, data_var))  # Store (display name, column name)
            
            if not present_vars:
                print("  No target variables found - skipping time series plot")
                continue
                
            # Create figure with subplots (3 panels per row)
            num_vars = len(present_vars)
            ncols = 3
            nrows = (num_vars + ncols - 1) // ncols
            
            # Set dimensions for full-screen viewing (1920x1080 equivalent)
            fig_width = 24
            fig_height = 14
            
            # Create figure with 3 panels per row
            fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
            if nrows == 1:
                axes = axes.reshape(1, -1)  # Ensure 2D array even for single row
            
            # Flatten axes array for easy iteration
            flat_axes = axes.flatten()
            
            # Get unique years for coloring
            years = sorted(growing_df['year'].unique())
            year_colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
            
            # Plot each variable in its own subplot
            for i, (display_var, data_var) in enumerate(present_vars):
                ax = flat_axes[i]
                
                # Plot each year separately
                for j, year in enumerate(years):
                    year_data = growing_df[growing_df['year'] == year]
                    if not year_data.empty:
                        # Downsample to daily means for better visibility
                        daily = year_data.set_index('DateTime').resample('D').mean()
                        ax.plot(daily.index, daily[data_var], 
                                color=year_colors[j], alpha=0.8, linewidth=1.5, label=str(year))
                
                ax.set_title(display_var, fontsize=16)
                ax.set_ylabel(display_var, fontsize=14)
                ax.grid(alpha=0.3)
                ax.tick_params(axis='both', labelsize=12)
                
                # Add SOS/EOS markers
                for j, year in enumerate(years):
                    year_data = growing_df[growing_df['year'] == year]
                    if not year_data.empty:
                        sos_date = year_data[year_data['DOY'] >= year_data['avg_sos'].iloc[0]].iloc[0]['DateTime']
                        eos_date = year_data[year_data['DOY'] <= year_data['avg_eos'].iloc[0]].iloc[-1]['DateTime']
                        ax.axvline(sos_date, color=year_colors[j], linestyle='--', alpha=0.7, linewidth=1.5)
                        ax.axvline(eos_date, color=year_colors[j], linestyle='--', alpha=0.7, linewidth=1.5)
            
            # SIMPLIFIED X-AXIS: ONLY SHOW YEARS
            # Create year boundaries for x-ticks
            year_boundaries = [pd.Timestamp(f'{year}-01-01') for year in years]
            
            # Apply to all subplots
            for ax in flat_axes:
                ax.set_xticks(year_boundaries)
                ax.set_xticklabels(years, fontsize=12)
                ax.set_xlim(year_boundaries[0], year_boundaries[-1] + pd.Timedelta(days=365))
            
            # Hide unused axes
            for i in range(len(present_vars), len(flat_axes)):
                flat_axes[i].axis('off')
            
            # Create custom legend for years
            legend_elements = [Line2D([0], [0], color=year_colors[i], lw=3, label=str(year)) 
                              for i, year in enumerate(years)]
            fig.legend(handles=legend_elements, loc='upper center', 
                      bbox_to_anchor=(0.5, 1.03), ncol=min(10, len(years)), 
                      title='Year', fontsize=14, title_fontsize=16)
            
            plt.suptitle(f'{site} Growing Season Time Series', fontsize=20, y=0.99)
            plt.tight_layout(rect=[0, 0, 1, 0.97])  # Make space for top legend
            
            # Save high-resolution plot optimized for full-screen viewing
            ts_plot_path = os.path.join(save_folder, f"{site}_growing_season_timeseries.png")
            plt.savefig(ts_plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            # Save processed data
            data_path = os.path.join(save_folder, f"{site}_growing_season.csv")
            growing_df.to_csv(data_path, index=False)
            print(f"  Saved data and plots for {site}")
            
        except Exception as e:
            print(f"  Error processing {site}: {str(e)}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    driver_growing_season()










#################################################################################################################








