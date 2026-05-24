# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 17:27:52 2025

@author: ammar
"""
## tested on 28 ameriflux sites
## see paths to find path to csv


# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 17:27:52 2025
@author: ammar
Automated AmeriFlux data processing with full variable extraction
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
from datetime import datetime

# Define folder paths
input_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps'
output_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps'

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

def detect_header_and_sep(file_path):
    """
    Detect header row and separator for AmeriFlux files
    Returns: (skip_rows, separator)
    """
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        lines = []
        for i in range(10):  # Read first 10 lines
            try:
                lines.append(f.readline())
            except:
                break
    
    # Check for tab-separated vs comma-separated
    for line in lines:
        if '\t' in line and ',' not in line:
            sep = '\t'
            break
        elif ',' in line:
            sep = ','
            break
    else:
        sep = ','
    
    # Find the first line that looks like a header (has TIMESTAMP or similar)
    for i, line in enumerate(lines):
        if 'TIMESTAMP' in line.upper() or 'DATETIME' in line.upper():
            return i, sep
    
    # If no TIMESTAMP found, try to find line with numeric values
    for i, line in enumerate(lines):
        # Skip lines that start with '#' or are comments
        if line.strip().startswith('#') or line.strip().startswith('Site:'):
            continue
        # Check if line contains what looks like column headers
        if any(word in line.upper() for word in ['TIMESTAMP', 'DATE', 'TIME', 'NEE', 'LE', 'H']):
            return i, sep
    
    # Default: assume header is at row 0
    return 0, sep

def robust_datetime_parsing(series):
    """Robust datetime parsing for multiple formats"""
    # Try different formats
    formats = [
        '%Y%m%d%H%M',      # 201401010000
        '%Y%m%d %H%M',     # 20140101 0000
        '%Y-%m-%d %H:%M:%S', # 2014-01-01 00:00:00
        '%Y-%m-%d %H:%M',    # 2014-01-01 00:00
        '%m/%d/%Y %H:%M',    # 1/1/2014 0:00
        '%Y%m%d%H%M%S',      # 20140101000000
    ]
    
    # Convert scientific notation to string if needed
    if hasattr(series, 'dtype') and series.dtype in ['float64', 'int64']:
        series = series.astype(str)
        # Remove decimal points and scientific notation
        series = series.str.replace(r'\.0+$', '', regex=True)
        series = series.str.replace(r'e\+', '', regex=True)
    
    # Try each format
    for fmt in formats:
        try:
            parsed = pd.to_datetime(series, format=fmt, errors='coerce')
            if parsed.notna().any():
                print(f"    Parsed datetime with format: {fmt}")
                return parsed
        except:
            continue
    
    # If all fail, use pandas flexible parsing
    print("    Using flexible datetime parsing")
    return pd.to_datetime(series, errors='coerce')

def select_best_column(df, column_patterns, prefer_f=True):
    """
    Select the best column from a list of patterns
    Prefers columns with _F, _PI_F, etc. if prefer_f is True
    """
    available_cols = []
    for pattern in column_patterns:
        # Match exact column names or patterns
        matching_cols = [col for col in df.columns if re.match(pattern, str(col), re.IGNORECASE)]
        available_cols.extend(matching_cols)
    
    # Remove duplicates while preserving order
    available_cols = list(dict.fromkeys(available_cols))
    
    if not available_cols:
        return None
    
    # If prefer_f is True, try to find columns with _F, _PI_F, etc.
    if prefer_f:
        f_patterns = [r'.*_F$', r'.*_PI_F$', r'.*_F_', r'.*_PI_F_']
        for pattern in f_patterns:
            f_cols = [col for col in available_cols if re.match(pattern, col, re.IGNORECASE)]
            if f_cols:
                return f_cols[0]
    
    # Otherwise return the first available column
    return available_cols[0]

def check_column_validity(df, col_name):
    """
    Check if a column has any valid data (not all NaN or -9999)
    Returns True if column has valid data, False otherwise
    """
    if col_name is None or col_name not in df.columns:
        return False
    
    # Get column values and replace -9999 with NaN
    col_values = pd.to_numeric(df[col_name], errors='coerce')
    col_values = col_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
    
    # Check if there are any valid values
    has_valid = col_values.notna().any()
    
    if not has_valid:
        print(f"    Column {col_name} has no valid data (all NaN or -9999)")
    
    return has_valid

def process_nee_with_fallback(df):
    """
    Process NEE with comprehensive fallback logic:
    1. Direct NEE columns (NEE_PI_F, NEE_F, NEE_PI, NEE) - only if they have valid data
    2. FC + SC (any version including _1_1_1)
    3. FC + SFC (any version including _1_1_1)
    4. FC_PI_F
    5. FC alone (last resort)
    6. NaN if nothing exists
    """
    print("  Processing NEE with fallback logic...")
    
    # Step 1: Try direct NEE columns first - only if they have valid data
    nee_patterns = [
        r'^NEE_PI_F$', r'^NEE_PI_F_', 
        r'^NEE_F$', r'^NEE_F_',
        r'^NEE_PI$', r'^NEE_PI_', 
        r'^NEE$', r'^NEE_'
    ]
    
    # Check each NEE pattern in priority order
    for pattern in nee_patterns:
        matching_cols = [col for col in df.columns if re.match(pattern, str(col), re.IGNORECASE)]
        for col in matching_cols:
            if check_column_validity(df, col):
                # Found a valid NEE column
                df["NEE"] = pd.to_numeric(df[col], errors='coerce')
                # Replace -9999 and -6999 with NaN
                df["NEE"] = df["NEE"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
                print(f"    Step 1: Using direct NEE column: {col} (has valid data)")
                return df
            else:
                print(f"    Skipping {col} - no valid data")
    
    print("    Step 1 failed: No valid direct NEE column found")
    df["NEE"] = np.nan
    
    # Define component patterns
    fc_patterns = [r'^FC$', r'^FC_', r'^FC_1_1_1$']
    sc_patterns = [r'^SC$', r'^SC_', r'^SC_1_1_1$']
    sfc_patterns = [r'^SFC$', r'^SFC_', r'^SFC_1_1_1$']
    fc_pi_f_patterns = [r'^FC_PI_F$', r'^FC_PI_F_']
    
    # Find available components (check if they have valid data)
    fc_col = None
    for pattern in fc_patterns:
        matching_cols = [col for col in df.columns if re.match(pattern, str(col), re.IGNORECASE)]
        for col in matching_cols:
            if check_column_validity(df, col):
                fc_col = col
                break
        if fc_col:
            break
    
    sc_col = None
    for pattern in sc_patterns:
        matching_cols = [col for col in df.columns if re.match(pattern, str(col), re.IGNORECASE)]
        for col in matching_cols:
            if check_column_validity(df, col):
                sc_col = col
                break
        if sc_col:
            break
    
    sfc_col = None
    for pattern in sfc_patterns:
        matching_cols = [col for col in df.columns if re.match(pattern, str(col), re.IGNORECASE)]
        for col in matching_cols:
            if check_column_validity(df, col):
                sfc_col = col
                break
        if sfc_col:
            break
    
    fc_pi_f_col = None
    for pattern in fc_pi_f_patterns:
        matching_cols = [col for col in df.columns if re.match(pattern, str(col), re.IGNORECASE)]
        for col in matching_cols:
            if check_column_validity(df, col):
                fc_pi_f_col = col
                break
        if fc_pi_f_col:
            break
    
    # Step 2A: Try FC + SC
    if fc_col and sc_col:
        fc_values = pd.to_numeric(df[fc_col], errors='coerce')
        sc_values = pd.to_numeric(df[sc_col], errors='coerce')
        # Replace missing value codes
        fc_values = fc_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        sc_values = sc_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        df["NEE"] = fc_values + sc_values
        print(f"    Step 2A: Using {fc_col} + {sc_col}")
        return df
    
    # Step 2B: Try FC + SFC
    if fc_col and sfc_col:
        fc_values = pd.to_numeric(df[fc_col], errors='coerce')
        sfc_values = pd.to_numeric(df[sfc_col], errors='coerce')
        # Replace missing value codes
        fc_values = fc_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        sfc_values = sfc_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        df["NEE"] = fc_values + sfc_values
        print(f"    Step 2B: Using {fc_col} + {sfc_col}")
        return df
    
    # Step 3: Try FC_PI_F
    if fc_pi_f_col:
        fc_pi_f_values = pd.to_numeric(df[fc_pi_f_col], errors='coerce')
        fc_pi_f_values = fc_pi_f_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        df["NEE"] = fc_pi_f_values
        print(f"    Step 3: Using {fc_pi_f_col}")
        return df
    
    # Step 4: Last resort - FC alone
    if fc_col:
        fc_values = pd.to_numeric(df[fc_col], errors='coerce')
        fc_values = fc_values.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        df["NEE"] = fc_values
        print(f"    Step 4 (Last Resort): Using {fc_col} alone")
        return df
    
    # Step 5: Nothing available
    print("    Step 5: No NEE data available - keeping as NaN")
    df["NEE"] = np.nan
    return df

def extract_and_filter_variables(df):
    """Extract and calculate required variables with filtering (excluding CO2)"""
    print("  Extracting additional variables (PAR, WS, WD)...")
    
    # PAR extraction (PPFD_IN or similar)
    par_patterns = [r'^PPFD_IN$', r'^PPFD_IN_', r'^PAR$', r'^PAR_']
    par_col = select_best_column(df, par_patterns, prefer_f=False)
    if par_col and check_column_validity(df, par_col):
        df['PAR'] = pd.to_numeric(df[par_col], errors='coerce')
        df['PAR'] = df['PAR'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        print(f"    Using PAR column: {par_col}")
        
        # Apply PAR filter (no negative values)
        if pd.api.types.is_numeric_dtype(df['PAR']):
            par_count_before = df['PAR'].count()
            df['PAR'] = np.where(df['PAR'] < 0, np.nan, df['PAR'])
            par_count_after = df['PAR'].count()
            if par_count_before - par_count_after > 0:
                print(f"    PAR filtered: {par_count_before - par_count_after} values removed")
    else:
        df['PAR'] = np.nan
        print("    No valid PAR columns found")
    
    # WS (Wind Speed) extraction
    ws_patterns = [r'^WS$', r'^WS_', r'^WIND_SPEED$', r'^WIND_SPEED_']
    ws_col = select_best_column(df, ws_patterns, prefer_f=False)
    if ws_col and check_column_validity(df, ws_col):
        df['WS'] = pd.to_numeric(df[ws_col], errors='coerce')
        df['WS'] = df['WS'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        print(f"    Using WS column: {ws_col}")
        
        # Apply WS filter (0-20 m/s range)
        if pd.api.types.is_numeric_dtype(df['WS']):
            ws_count_before = df['WS'].count()
            df['WS'] = np.where((df['WS'] > 20) | (df['WS'] < 0), np.nan, df['WS'])
            ws_count_after = df['WS'].count()
            if ws_count_before - ws_count_after > 0:
                print(f"    WS filtered: {ws_count_before - ws_count_after} values removed")
    else:
        df['WS'] = np.nan
        print("    No valid WS columns found")
    
    # WD (Wind Direction) extraction
    wd_patterns = [r'^WD$', r'^WD_', r'^WIND_DIR$', r'^WIND_DIR_']
    wd_col = select_best_column(df, wd_patterns, prefer_f=False)
    if wd_col and check_column_validity(df, wd_col):
        df['WD'] = pd.to_numeric(df[wd_col], errors='coerce')
        df['WD'] = df['WD'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        print(f"    Using WD column: {wd_col}")
        
        # Apply WD filter (0-360 degrees range)
        if pd.api.types.is_numeric_dtype(df['WD']):
            wd_count_before = df['WD'].count()
            df['WD'] = np.where((df['WD'] > 360) | (df['WD'] < 0), np.nan, df['WD'])
            wd_count_after = df['WD'].count()
            if wd_count_before - wd_count_after > 0:
                print(f"    WD filtered: {wd_count_before - wd_count_after} values removed")
    else:
        df['WD'] = np.nan
        print("    No valid WD columns found")
    
    return df

def process_ameriflux_data(csv_file_path, save_folder):
    """
    Process a single AmeriFlux CSV file and save the processed data
    """
    try:
        # Extract site name from filename
        filename = os.path.basename(csv_file_path)
        match = re.search(r'AMF_([^_]+)_BASE', filename)
        if match:
            site_name = match.group(1)
        else:
            parts = filename.split('_')
            if len(parts) >= 2:
                site_name = parts[1]
            else:
                site_name = filename.replace('.csv', '')
        
        print(f"\n{'='*60}")
        print(f"Processing site: {site_name}")
        print(f"File: {filename}")
        print(f"{'='*60}")
        
        # Detect header row and separator
        skip_rows, sep = detect_header_and_sep(csv_file_path)
        print(f"  Detected: skip_rows={skip_rows}, separator='{sep}'")
        
        # Read CSV file
        try:
            # Read with low_memory=False to avoid dtype warnings
            df = pd.read_csv(csv_file_path, skiprows=skip_rows, sep=sep, low_memory=False)
            print(f"  Data loaded successfully")
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns)[:10]}...")
        except Exception as e:
            print(f"  Error reading file: {e}")
            return None
        
        # Find timestamp column
        timestamp_col = None
        for col in ['TIMESTAMP_START', 'TIMESTAMP', 'datetime', 'DateTime', 'DATE']:
            if col in df.columns:
                timestamp_col = col
                break
        
        if not timestamp_col:
            print(f"  Error: No timestamp column found. Available columns: {list(df.columns)[:10]}")
            return None
        
        print(f"  Using timestamp column: {timestamp_col}")
        
        # Handle scientific notation in timestamps
        if df[timestamp_col].dtype in ['float64', 'int64']:
            # Convert to string and handle scientific notation
            df[timestamp_col] = df[timestamp_col].astype(str)
            df[timestamp_col] = df[timestamp_col].str.replace(r'\.0+$', '', regex=True)
            df[timestamp_col] = df[timestamp_col].str.replace(r'e\+', '', regex=True)
        
        # Parse timestamp
        df['TIMESTAMP'] = robust_datetime_parsing(df[timestamp_col])
        
        # Remove rows with invalid timestamps
        invalid_timestamps = df['TIMESTAMP'].isna().sum()
        if invalid_timestamps > 0:
            print(f"  Warning: {invalid_timestamps} rows with invalid timestamps removed")
            df = df.dropna(subset=['TIMESTAMP'])
        
        if df.empty:
            print("  Error: No valid timestamps after parsing")
            return None
        
        # Drop duplicates and set index
        df = df.drop_duplicates(subset='TIMESTAMP').copy()
        df.set_index('TIMESTAMP', inplace=True)
        
        # Generate complete 30-minute range
        start = df.index.min().floor('D')
        end = df.index.max().ceil('D') - pd.Timedelta(minutes=30)
        full_index = pd.date_range(start=start, end=end, freq='30min')
        
        # Reindex and fill missing values
        df = df.reindex(full_index)
        df.index.name = 'TIMESTAMP'
        
        # Replace missing value codes with NaN
        df = df.replace([-9999, -6999, -9999.0, -6999.0], np.nan)
        
        # Create required time columns
        df = df.reset_index()
        df['DateTime'] = df['TIMESTAMP']
        df['Year'] = df['TIMESTAMP'].dt.year
        df['DoY'] = df['TIMESTAMP'].dt.dayofyear
        df['Hour'] = df['TIMESTAMP'].dt.hour + df['TIMESTAMP'].dt.minute / 60 + 0.5
        df['Hour'] = df['Hour'].apply(lambda x: round(x * 2) / 2)
        
        # Process NEE with comprehensive fallback logic
        df = process_nee_with_fallback(df)
        
        # Filter NEE
        if pd.api.types.is_numeric_dtype(df["NEE"]):
            df["NEE"] = np.where(df["NEE"].between(-50, 50), df["NEE"], np.nan)
        
        # Process LE - include LE_1_1_1
        le_patterns = [
            r'^LE_PI_F$', r'^LE_PI_F_', r'^LE_F$', r'^LE_F_',
            r'^LE_PI$', r'^LE_PI_', r'^LE$', r'^LE_',
            r'^LE_1_1_1$'
        ]
        le_col = select_best_column(df, le_patterns, prefer_f=True)
        if le_col and check_column_validity(df, le_col):
            df["LE"] = pd.to_numeric(df[le_col], errors='coerce')
            df["LE"] = df["LE"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using LE column: {le_col}")
        else:
            df["LE"] = np.nan
            print("    No valid LE column found")
        
        if pd.api.types.is_numeric_dtype(df["LE"]):
            df["LE"] = df["LE"].apply(lambda x: np.nan if x < -200 or x > 800 else x)
        
        # Process H - include H_1_1_1
        h_patterns = [
            r'^H_PI_F$', r'^H_PI_F_', r'^H_F$', r'^H_F_',
            r'^H_PI$', r'^H_PI_', r'^H$', r'^H_',
            r'^H_1_1_1$'
        ]
        h_col = select_best_column(df, h_patterns, prefer_f=True)
        if h_col and check_column_validity(df, h_col):
            df["H"] = pd.to_numeric(df[h_col], errors='coerce')
            df["H"] = df["H"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using H column: {h_col}")
        else:
            df["H"] = np.nan
            print("    No valid H column found")
        
        if pd.api.types.is_numeric_dtype(df["H"]):
            df["H"] = df["H"].apply(lambda x: np.nan if x < -200 or x > 800 else x)
        
        # Process Rg (shortwave radiation) - include SW_IN_1_1_1
        rg_patterns = [
            r'^SW_IN_PI_F$', r'^SW_IN_F$', r'^SW_IN$', r'^Rg$',
            r'^SW_IN_1_1_1$', r'^SW_IN_'
        ]
        rg_col = select_best_column(df, rg_patterns, prefer_f=True)
        if rg_col and check_column_validity(df, rg_col):
            df["Rg"] = pd.to_numeric(df[rg_col], errors='coerce')
            df["Rg"] = df["Rg"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using Rg column: {rg_col}")
        else:
            df["Rg"] = np.nan
            print("    No valid Rg column found")
        
        if pd.api.types.is_numeric_dtype(df["Rg"]):
            df["Rg"] = df["Rg"].apply(lambda x: 0 if x < 0 else x)
        
        # Process Tair - restricted to specific patterns only
        tair_patterns = [
            r'^TA_PI_F$', r'^TA_F$', r'^TA$', r'^TA_1_1_1$', r'^TA_1_2_1$'
        ]
        tair_col = select_best_column(df, tair_patterns, prefer_f=True)
        if tair_col and check_column_validity(df, tair_col):
            df["Tair"] = pd.to_numeric(df[tair_col], errors='coerce')
            df["Tair"] = df["Tair"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using Tair column: {tair_col}")
        else:
            df["Tair"] = np.nan
            print("    No valid Tair column found")
        
        # Process VPD
        vpd_patterns = [
            r'^VPD_PI_F$', r'^VPD_F$', r'^VPD_PI$', r'^VPD$', r'^VPD_'
        ]
        vpd_col = select_best_column(df, vpd_patterns, prefer_f=True)
        if vpd_col and check_column_validity(df, vpd_col):
            df["VPD"] = pd.to_numeric(df[vpd_col], errors='coerce')
            df["VPD"] = df["VPD"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using VPD column: {vpd_col}")
        else:
            df["VPD"] = np.nan
            print("    No valid VPD column found")
        
        # Process RH - include RH_1_1_1 and RH_PI_F
        rh_patterns = [
            r'^RH_PI_F$', r'^RH_PI_F_', r'^RH_F$', r'^RH_F_',
            r'^RH_PI$', r'^RH_PI_', r'^RH$', r'^RH_',
            r'^RH_1_1_1$'
        ]
        rh_col = select_best_column(df, rh_patterns, prefer_f=True)
        if rh_col and check_column_validity(df, rh_col):
            df["RH"] = pd.to_numeric(df[rh_col], errors='coerce')
            df["RH"] = df["RH"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using RH column: {rh_col}")
        else:
            df["RH"] = np.nan
            print("    No valid RH column found")
        
        # Calculate VPD from Tair and RH if needed
        if "VPD" in df.columns and df["VPD"].isna().all() and "Tair" in df.columns and "RH" in df.columns:
            if not df["RH"].isna().all():
                es = 0.6108 * np.exp((17.27 * df["Tair"]) / (df["Tair"] + 237.3))
                ea = df["RH"] / 100 * es
                df["VPD"] = (es - ea) * 10
                print("    Calculated VPD from Tair and RH")
        
        if pd.api.types.is_numeric_dtype(df["VPD"]):
            df["VPD"] = df["VPD"].apply(lambda x: np.nan if x < 0 else x)
        
        # Process USTAR - include USTAR_1_1_1
        ustar_patterns = [
            r'^USTAR_PI_F$', r'^USTAR_F$', r'^USTAR$', r'^UST$',
            r'^USTAR_1_1_1$', r'^UST_'
        ]
        ustar_col = select_best_column(df, ustar_patterns, prefer_f=True)
        if ustar_col and check_column_validity(df, ustar_col):
            df["USTAR"] = pd.to_numeric(df[ustar_col], errors='coerce')
            df["USTAR"] = df["USTAR"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using USTAR column: {ustar_col}")
        else:
            df["USTAR"] = np.nan
            print("    No valid USTAR column found")
        df["Ustar"] = df["USTAR"]
        
        # Process NETRAD - include LW_IN, LW_OUT, SW_IN, SW_OUT variants
        netrad_patterns = [
            r'^NETRAD$', r'^NETRAD_', r'^Rn$', r'^Rn_'
        ]
        netrad_col = select_best_column(df, netrad_patterns, prefer_f=True)
        if netrad_col and check_column_validity(df, netrad_col):
            df["NETRAD"] = pd.to_numeric(df[netrad_col], errors='coerce')
            df["NETRAD"] = df["NETRAD"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using NETRAD column: {netrad_col}")
        else:
            print("    Calculating NETRAD from components")
            # Look for component columns with various naming patterns
            sw_in_patterns = [r'^SW_IN$', r'^SW_IN_1_1_1$', r'^SW_IN_']
            sw_out_patterns = [r'^SW_OUT$', r'^SW_OUT_1_1_1$', r'^SW_OUT_']
            lw_in_patterns = [r'^LW_IN$', r'^LW_IN_1_1_1$', r'^LW_IN_']
            lw_out_patterns = [r'^LW_OUT$', r'^LW_OUT_1_1_1$', r'^LW_OUT_']
            
            sw_in_col = select_best_column(df, sw_in_patterns, prefer_f=False)
            sw_out_col = select_best_column(df, sw_out_patterns, prefer_f=False)
            lw_in_col = select_best_column(df, lw_in_patterns, prefer_f=False)
            lw_out_col = select_best_column(df, lw_out_patterns, prefer_f=False)
            
            if sw_in_col and check_column_validity(df, sw_in_col):
                df['SW_IN'] = pd.to_numeric(df[sw_in_col], errors='coerce')
                df['SW_IN'] = df['SW_IN'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            if sw_out_col and check_column_validity(df, sw_out_col):
                df['SW_OUT'] = pd.to_numeric(df[sw_out_col], errors='coerce')
                df['SW_OUT'] = df['SW_OUT'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            if lw_in_col and check_column_validity(df, lw_in_col):
                df['LW_IN'] = pd.to_numeric(df[lw_in_col], errors='coerce')
                df['LW_IN'] = df['LW_IN'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            if lw_out_col and check_column_validity(df, lw_out_col):
                df['LW_OUT'] = pd.to_numeric(df[lw_out_col], errors='coerce')
                df['LW_OUT'] = df['LW_OUT'].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            
            # Calculate NETRAD when components are available
            if all(col in df.columns for col in ['SW_IN', 'SW_OUT', 'LW_IN', 'LW_OUT']):
                condition = df[['SW_IN', 'SW_OUT', 'LW_IN', 'LW_OUT']].notna().all(axis=1)
                df.loc[condition, 'NETRAD'] = (df['SW_IN'] - df['SW_OUT']) + \
                                              (df['LW_IN'] - df['LW_OUT'])
                print("    Calculated NETRAD from components")
            else:
                df["NETRAD"] = np.nan
        
        # Process PA
        pa_patterns = [r'^PA$', r'^PA_1_1_1$', r'^PA_']
        pa_col = select_best_column(df, pa_patterns, prefer_f=False)
        if pa_col and check_column_validity(df, pa_col):
            df["PA"] = pd.to_numeric(df[pa_col], errors='coerce')
            df["PA"] = df["PA"].replace([-9999, -6999, -9999.0, -6999.0], np.nan)
            print(f"    Using PA column: {pa_col}")
        else:
            df["PA"] = np.nan
            print("    No valid PA column found")
        
        # Extract additional variables (PAR, WS, WD)
        df = extract_and_filter_variables(df)
        
        # Define final columns to save
        final_columns = ["DateTime", "Year", "DoY", "Hour", "NEE", "LE", "H", 
                        "Rg", "PAR", "WS", "WD", "Tair", "VPD", "Ustar", "PA", "NETRAD", "RH"]
        
        # Select only columns that exist
        existing_columns = [col for col in final_columns if col in df.columns]
        reframed = df[existing_columns].copy()
        
        # Apply PA filter
        if 'PA' in reframed.columns and pd.api.types.is_numeric_dtype(reframed['PA']):
            pa_count_before = reframed['PA'].count()
            reframed.loc[reframed['PA'] > 150, 'PA'] = np.nan
            pa_count_after = reframed['PA'].count()
            pa_filtered = pa_count_before - pa_count_after
            if pa_filtered > 0:
                print(f"  PA filtered: {pa_filtered} values >150 kPa removed")
        
        # Save to CSV
        new_filename = f"gaps_{site_name}.csv"
        csv_save_path = os.path.join(save_folder, new_filename)
        reframed.to_csv(csv_save_path, index=False)
        print(f"  Data saved to: {csv_save_path}")
        
        # Print data completeness summary
        print(f"  Data completeness (% non-missing):")
        for col in final_columns:
            if col in reframed.columns:
                pct = reframed[col].count() / len(reframed) * 100
                print(f"    {col}: {pct:.1f}%")
        
        # Create diagnostic plots
        create_diagnostic_plots(reframed, site_name, save_folder)
        
        return reframed
        
    except Exception as e:
        print(f"  Error processing file: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_diagnostic_plots(df, site_name, save_folder):
    """Create diagnostic plots for key variables"""
    try:
        os.makedirs(save_folder, exist_ok=True)
        
        df_plot = df.copy()
        df_plot.set_index("DateTime", inplace=True)
        
        columns_to_plot = [col for col in df_plot.columns if col not in ["Year", "DoY", "Hour"]]
        
        n_cols = min(3, len(columns_to_plot))
        n_rows = (len(columns_to_plot) + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
        fig.suptitle(f'{site_name} - Processed Variables', fontsize=16, fontweight='bold')
        
        if n_rows == 1 and n_cols == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for idx, col in enumerate(columns_to_plot):
            if idx < len(axes):
                ax = axes[idx]
                ax.plot(df_plot.index, df_plot[col], color='b', alpha=0.7, linewidth=0.5)
                ax.set_xlabel('Date')
                ax.set_ylabel(col)
                ax.set_title(f'{col} - {df_plot[col].count() / len(df_plot) * 100:.1f}% available')
                ax.grid(True, alpha=0.3)
                ax.tick_params(axis='x', rotation=45)
        
        for idx in range(len(columns_to_plot), len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plot_path = os.path.join(save_folder, f"{site_name}_diagnostics.png")
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Diagnostic plot saved to: {plot_path}")
        
    except Exception as e:
        print(f"  Error creating diagnostic plot: {e}")

def process_all_sites(input_folder, output_folder):
    """Process all CSV files in the input folder"""
    csv_files = [f for f in os.listdir(input_folder) 
                 if f.lower().endswith('.csv') and os.path.isfile(os.path.join(input_folder, f))]
    
    if not csv_files:
        print(f"No CSV files found in: {input_folder}")
        return
    
    print(f"\n{'='*60}")
    print(f"Found {len(csv_files)} CSV files to process")
    print(f"Input folder: {input_folder}")
    print(f"Output folder: {output_folder}")
    print(f"{'='*60}")
    
    successful = 0
    failed = 0
    
    for filename in csv_files:
        file_path = os.path.join(input_folder, filename)
        result = process_ameriflux_data(file_path, output_folder)
        
        if result is not None:
            successful += 1
        else:
            failed += 1
    
    print(f"\n{'='*60}")
    print(f"Processing Complete!")
    print(f"Successfully processed: {successful} sites")
    print(f"Failed: {failed} sites")
    print(f"{'='*60}")

# Main execution
if __name__ == "__main__":
    process_all_sites(input_folder, output_folder)
    print("\nAll processing completed!")



####################################################################
# check already 
# US-Brw,US-EvM, US-StS , US-Atq, US-Snd, US-S04,  VPD fix  
#US-Brw Rg, P  # yes missing
#HB4, PAR, Tair # true
#HPY PAR # yes
#KS1 Rg # yes
#KS2 Tair
#KS4  PAR
#LA1, Rg, PA, Netrad
#MRM PAR
#SCs pAR ,PA
#Stj Rg



##############################################################33
# next rounf check PHM NEE


#EKH , EKP, netrad no, check if different varibale name , check NEE

#GCE, no ustar LE, NEE,H
#NC4 NEE, LE ,H

