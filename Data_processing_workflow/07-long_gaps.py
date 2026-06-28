# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 12:36:14 2025
@author: ammar
Function to fill long gaps (>=7 days) in LE and NEE during April and September
using climatology from the same DoY and Hour across all years
Includes comprehensive data filtering and PA replacement with ERA5 sp
"""

import os
import pandas as pd
import numpy as np

def apply_data_filters(df):
    """
    Apply comprehensive data filters to remove unrealistic values
    Replaces extreme values with NaN (does not remove rows)
    
    Parameters:
    - df: DataFrame with AmeriFlux data
    
    Returns:
    - df: DataFrame with filtered values
    """
    print("  Applying data filters...")
    filtered_counts = {}
    
    # WS filter: <0 or >20 m/s
    if 'WS' in df.columns:
        ws_before = df['WS'].count()
        df.loc[(df['WS'] < 0) | (df['WS'] > 20), 'WS'] = np.nan
        ws_after = df['WS'].count()
        filtered_counts['WS'] = ws_before - ws_after
        if filtered_counts['WS'] > 0:
            print(f"    WS: removed {filtered_counts['WS']} values outside [0, 20]")
    
    # PAR filter: <0
    if 'PAR' in df.columns:
        par_before = df['PAR'].count()
        df.loc[df['PAR'] < 0, 'PAR'] = np.nan
        par_after = df['PAR'].count()
        filtered_counts['PAR'] = par_before - par_after
        if filtered_counts['PAR'] > 0:
            print(f"    PAR: removed {filtered_counts['PAR']} negative values")
    
    # Rg filter: <0
    if 'Rg' in df.columns:
        rg_before = df['Rg'].count()
        df.loc[df['Rg'] < 0, 'Rg'] = np.nan
        rg_after = df['Rg'].count()
        filtered_counts['Rg'] = rg_before - rg_after
        if filtered_counts['Rg'] > 0:
            print(f"    Rg: removed {filtered_counts['Rg']} negative values")
    
    # RH filter: <0 or >100%
    if 'RH' in df.columns:
        rh_before = df['RH'].count()
        df.loc[(df['RH'] < 0) | (df['RH'] > 100), 'RH'] = np.nan
        rh_after = df['RH'].count()
        filtered_counts['RH'] = rh_before - rh_after
        if filtered_counts['RH'] > 0:
            print(f"    RH: removed {filtered_counts['RH']} values outside [0, 100]")
    
    # VPD filter: <0 or >70 hPa
    if 'VPD' in df.columns:
        vpd_before = df['VPD'].count()
        df.loc[(df['VPD'] < 0) | (df['VPD'] > 70), 'VPD'] = np.nan
        vpd_after = df['VPD'].count()
        filtered_counts['VPD'] = vpd_before - vpd_after
        if filtered_counts['VPD'] > 0:
            print(f"    VPD: removed {filtered_counts['VPD']} values outside [0, 70]")
    
    # H filter: <-200 or >800
    if 'H' in df.columns:
        h_before = df['H'].count()
        df.loc[(df['H'] < -200) | (df['H'] > 800), 'H'] = np.nan
        h_after = df['H'].count()
        filtered_counts['H'] = h_before - h_after
        if filtered_counts['H'] > 0:
            print(f"    H: removed {filtered_counts['H']} values outside [-200, 800]")
    
    # NETRAD filter: >1500
    if 'NETRAD' in df.columns:
        netrad_before = df['NETRAD'].count()
        df.loc[df['NETRAD'] > 1500, 'NETRAD'] = np.nan
        netrad_after = df['NETRAD'].count()
        filtered_counts['NETRAD'] = netrad_before - netrad_after
        if filtered_counts['NETRAD'] > 0:
            print(f"    NETRAD: removed {filtered_counts['NETRAD']} values > 1500")
    
    return df

def handle_pa_column(df):
    """
    Handle PA column:
    1. Remove values > 150 (set to NaN)
    2. If entire PA column is NaN, replace with ERA5 sp (converted from Pa to kPa)
    
    Parameters:
    - df: DataFrame with PA and sp columns
    
    Returns:
    - df: DataFrame with processed PA
    - bool: Whether PA was replaced with ERA5
    """
    # First filter PA values > 150
    if 'PA' in df.columns:
        pa_before = df['PA'].count()
        df.loc[df['PA'] > 150, 'PA'] = np.nan
        pa_after = df['PA'].count()
        if pa_before - pa_after > 0:
            print(f"    PA: removed {pa_before - pa_after} values > 150")
        
        # Check if entire PA column is NaN
        if df['PA'].isna().all():
            print(f"    PA: Entire column is NaN - replacing with ERA5 sp (converted from Pa to kPa)")
            if 'sp' in df.columns:
                # sp is in Pa, convert to kPa by dividing by 1000
                df['PA'] = df['sp'] / 1000
                # Round to 3 decimal places
                df['PA'] = df['PA'].round(3)
                print(f"    PA replaced with sp/1000 (kPa)")
                return df, True
            else:
                print(f"    Warning: sp column not found, cannot replace PA")
                return df, False
    else:
        print(f"    PA column not found")
    
    return df, False

def fill_long_nans_batch(path, save_path):
    """
    Batch process all CSV files to fill long gaps in LE and NEE
    Includes data filtering and PA handling
    
    Parameters:
    - path: Input directory containing CSV files
    - save_path: Output directory for processed files
    
    Returns:
    - dict: Statistics for each processed file
    """
    # Create output directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)
    
    # Get all CSV files in the input directory
    csv_files = [f for f in os.listdir(path) if f.endswith('.csv')]
    
    if not csv_files:
        print(f"No CSV files found in {path}")
        return {}
    
    print(f"Found {len(csv_files)} CSV files to process")
    print(f"{'='*60}")
    
    results = {}
    successful = 0
    failed = 0
    pa_replaced_sites = []
    
    for file in csv_files:
        print(f"\n{'='*50}")
        print(f"Processing: {file}")
        print(f"{'='*50}")
        
        try:
            # Read the CSV file
            input_path = os.path.join(path, file)
            
            # Determine which timestamp column exists
            df_test = pd.read_csv(input_path, nrows=5)
            
            if 'TIMESTAMP' in df_test.columns:
                timestamp_col = 'TIMESTAMP'
            elif 'DateTime' in df_test.columns:
                timestamp_col = 'DateTime'
            elif 'datetime' in df_test.columns:
                timestamp_col = 'datetime'
            else:
                print(f"  ✗ Error: No timestamp column found")
                failed += 1
                continue
            
            # Read full data with the identified timestamp column
            df = pd.read_csv(input_path, parse_dates=[timestamp_col])
            
            # Rename to standard TIMESTAMP for consistency
            if timestamp_col != 'TIMESTAMP':
                df.rename(columns={timestamp_col: 'TIMESTAMP'}, inplace=True)
            
            # Sort by TIMESTAMP to ensure correct order
            df = df.sort_values('TIMESTAMP').reset_index(drop=True)
            
            # Create time components if they don't exist
            if 'Year' not in df.columns:
                df['Year'] = df['TIMESTAMP'].dt.year
            if 'Month' not in df.columns:
                df['Month'] = df['TIMESTAMP'].dt.month
            if 'DoY' not in df.columns:
                df['DoY'] = df['TIMESTAMP'].dt.dayofyear
            if 'Hour' not in df.columns:
                df['Hour'] = df['TIMESTAMP'].dt.hour + df['TIMESTAMP'].dt.minute / 60
                df['Hour'] = df['Hour'].apply(lambda x: round(x * 2) / 2)
            
            # Apply data filters first
            df = apply_data_filters(df)
            
            # Handle PA column (replace with ERA5 sp if entire column is NaN)
            df, pa_replaced = handle_pa_column(df)
            if pa_replaced:
                pa_replaced_sites.append(file)
            
            # Make a copy for processing
            df_b = df.copy()
            
            # Create climatology from all years
            print(f"  Creating climatology from all years...")
            climatology = df.groupby(['DoY', 'Hour'])[['LE', 'NEE']].mean().reset_index()
            
            # Track statistics for this file
            stats = {'LE_filled': 0, 'NEE_filled': 0, 'total_gaps_filled': 0}
            
            # Process each year separately
            years = sorted(df['Year'].unique())
            
            for year in years:
                df_year = df[df['Year'] == year]
                growing_season = df_year[df_year['Month'].isin([5, 6, 7, 8])]
                n_total = len(growing_season)
                
                for var in ['LE', 'NEE']:
                    if var not in df.columns:
                        continue
                    
                    n_available = growing_season[var].notna().sum()
                    
                    # Check if we have sufficient growing season data (>=50%)
                    if n_total > 0 and (n_available / n_total) >= 0.5:
                        # Process April (month 4) and September (month 9)
                        for month in [4, 9]:
                            # Get data for this year and month
                            month_mask = (df_b['Year'] == year) & (df_b['Month'] == month)
                            month_data = df_b.loc[month_mask]
                            
                            if len(month_data) == 0:
                                continue
                            
                            # Identify long gaps (>= 7 days = 336 half-hour periods)
                            values = df_b.loc[month_mask, var]
                            is_nan = values.isna()
                            
                            # Find consecutive NaN sequences
                            group_id = (is_nan != is_nan.shift()).cumsum()
                            nan_groups = is_nan.groupby(group_id).transform('sum')
                            long_nan_mask = (is_nan) & (nan_groups >= 336)  # 7 days * 48 half-hours/day
                            
                            # Get indices where long gaps occur
                            long_gap_indices = month_mask & long_nan_mask
                            
                            if long_gap_indices.sum() > 0:
                                print(f"    Year {year}, {var}, month {month}: Found {long_gap_indices.sum()} long gap values")
                                
                                # Fill each long gap using climatology
                                for idx in df_b[long_gap_indices].index:
                                    row = df_b.loc[idx]
                                    doy = row['DoY']
                                    hour = row['Hour']
                                    
                                    # Find matching climatology value
                                    match = climatology[
                                        (climatology['DoY'] == doy) & 
                                        (climatology['Hour'] == hour)
                                    ]
                                    
                                    if not match.empty and pd.notna(match[var].values[0]):
                                        fill_value = round(match[var].values[0], 3)
                                        df_b.at[idx, var] = fill_value
                                        stats[f'{var}_filled'] += 1
                                        stats['total_gaps_filled'] += 1
            
            # Round all values to 3 decimal places
            for var in ['LE', 'NEE']:
                if var in df_b.columns:
                    df_b[var] = df_b[var].apply(lambda x: round(x, 3) if pd.notna(x) else x)
            
            # Save the processed file
            output_path = os.path.join(save_path, file)
            df_b.to_csv(output_path, index=False)
            
            print(f"  ✓ Processed successfully")
            print(f"    Total gaps filled: {stats['total_gaps_filled']} (LE: {stats['LE_filled']}, NEE: {stats['NEE_filled']})")
            
            results[file] = stats
            successful += 1
            
        except Exception as e:
            print(f"  ✗ Error processing {file}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
            continue
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Processing Complete!")
    print(f"Successfully processed: {successful} files")
    print(f"Failed: {failed} files")
    print(f"Output saved to: {save_path}")
    
    if pa_replaced_sites:
        print(f"\nPA replaced with ERA5 sp for {len(pa_replaced_sites)} sites:")
        for site in pa_replaced_sites:
            print(f"  - {site}")
    
    print(f"{'='*60}")
    
    return results

# Main execution
if __name__ == "__main__":
    # Define paths
    path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps'
    save_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps2'
    
    print("="*60)
    print("STARTING LONG GAP FILLING PROCESS WITH DATA FILTERING")
    print("="*60)
    print(f"Input path: {path}")
    print(f"Output path: {save_path}")
    print("\nData Filters Applied:")
    print("  - WS: <0 or >20 m/s → NaN")
    print("  - PAR: <0 → NaN")
    print("  - Rg: <0 → NaN")
    print("  - RH: <0 or >100% → NaN")
    print("  - VPD: <0 or >70 hPa → NaN")
    print("  - H: <-200 or >800 → NaN")
    print("  - NETRAD: >1500 → NaN")
    print("  - PA: >150 → NaN, if entire column NaN → replace with ERA5 sp/1000")
    print("\nGap Filling:")
    print("  - Filling long gaps (>=7 days) in LE and NEE during April and September")
    print("  - Using climatology from same DoY and Hour across all years")
    print("="*60)
    
    # Process all files
    results = fill_long_nans_batch(path, save_path)
    
    print("\nAll processing completed!")