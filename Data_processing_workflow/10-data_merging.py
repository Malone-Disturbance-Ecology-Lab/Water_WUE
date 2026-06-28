"""
@author: ammar
"""

import os
import pandas as pd
import numpy as np
###############################################################################################

## code for processing each site separately and saving individual files


def process_all_data(fill_folder, save_folder): 
    os.makedirs(save_folder, exist_ok=True)
    log_lines = []

    def write_log(msg):
        log_lines.append(msg)
        print(msg)

    processed_sites = []  # To keep track of processed sites
    csv_files = [f for f in os.listdir(fill_folder) if f.endswith('.csv')]

    for file in csv_files:
        file_path = os.path.join(fill_folder, file)
        try:
            df = pd.read_csv(file_path)
            
            # Extract site ID from filename (remove .csv extension)
            site_id = file.replace('.csv', '')
            write_log(f"\n{'='*60}")
            write_log(f"Processing site: {site_id}")
            write_log(f"{'='*60}")

            # Handle time
            if 'DateTime' in df.columns:
                df['DateTime'] = pd.to_datetime(df['DateTime'], errors='coerce')
                df['month'] = df['DateTime'].dt.month
                df['day'] = df['DateTime'].dt.day
                df['Hour'] = df['DateTime'].dt.hour + df['DateTime'].dt.minute / 60
                df['Year'] = df['DateTime'].dt.year
            else:
                write_log(f"Missing 'DateTime' column in {file}")
                continue

            # List of all variables to preserve
            preserve_vars = ['PA', 'RH', 'WS', 'WD', 'NEE_f', 'LE_f', 'H_f', 
                            'Tair_f', 'VPD_f', 'Rg_f', 'PAR_f', 'NETRAD_f']
            
            # Check and preserve existing variables
            for var in preserve_vars:
                if var not in df.columns:
                    df[var] = np.nan
                    write_log(f"Missing column '{var}' in {file}, filled with NaNs")
            
            # Required columns for unit conversion
            required_cols = ['NEE_f', 'GPP_DT', 'Reco_DT', 'GPP_nt', 'Reco_nt', 'Tair_f', 'LE_f']
            for col in required_cols:
                if col not in df.columns:
                    df[col] = np.nan
                    write_log(f"Missing column '{col}' in {file}, filled with NaNs")

            # Unit conversion (convert all to same units first)
            df['NEE'] = df['NEE_f'] * ((12 / 10**6) * 1800)  # g of C
            df['GPP_DT_converted'] = df['GPP_DT'] * ((12 / 10**6) * 1800)
            df['Reco_DT_converted'] = df['Reco_DT'] * ((12 / 10**6) * 1800)
            df['GPP_nt_converted'] = df['GPP_nt'] * ((12 / 10**6) * 1800)
            df['Reco_nt_converted'] = df['Reco_nt'] * ((12 / 10**6) * 1800)
            
            # Calculate lambda and ET
            df['lambda'] = (3149000 - 2370 * (df['Tair_f'] + 273.16)) * 1e-6
            df['ET'] = (df['LE_f'] / df['lambda']) * (1 / 1e6) * 1800  # in mm

            # ===================================================================
            # MODIFIED LOGIC: NEVER MIX DAYTIME AND NIGHTTIME DATA
            # Either use daytime OR nighttime, but never fill one with the other
            # ===================================================================
            
            # Step 1: Check if site has ANY valid daytime data (GPP_DT or Reco_DT)
            has_daytime_gpp = df['GPP_DT_converted'].notna().any()
            has_daytime_reco = df['Reco_DT_converted'].notna().any()
            
            # Step 2: Check if site has ANY valid nighttime data (GPP_nt or Reco_nt)
            has_nighttime_gpp = df['GPP_nt_converted'].notna().any()
            has_nighttime_reco = df['Reco_nt_converted'].notna().any()
            
            # Step 3: Initialize GPP column - PURE daytime OR PURE nighttime, no mixing
            if has_daytime_gpp:
                # Use daytime GPP exclusively
                df['GPP'] = df['GPP_DT_converted']
                write_log(f"Using ONLY daytime GPP (GPP_DT) - no nighttime mixing")
            elif has_nighttime_gpp:
                # Only use nighttime GPP if NO daytime data exists
                df['GPP'] = df['GPP_nt_converted']
                write_log(f"WARNING: No daytime GPP found, using ONLY nighttime GPP (GPP_nt) - no daytime mixing")
            else:
                # No GPP data at all
                df['GPP'] = np.nan
                write_log(f"ERROR: No GPP data (daytime or nighttime) found")
            
            # Step 4: Initialize Reco column - PURE daytime OR PURE nighttime, no mixing
            if has_daytime_reco:
                # Use daytime Reco exclusively
                df['Reco'] = df['Reco_DT_converted']
                write_log(f"Using ONLY daytime Reco (Reco_DT) - no nighttime mixing")
            elif has_nighttime_reco:
                # Only use nighttime Reco if NO daytime data exists
                df['Reco'] = df['Reco_nt_converted']
                write_log(f"WARNING: No daytime Reco found, using ONLY nighttime Reco (Reco_nt) - no daytime mixing")
            else:
                # No Reco data at all
                df['Reco'] = np.nan
                write_log(f"ERROR: No Reco data (daytime or nighttime) found")
            
            # Step 5: NO CROSS-FILLING between daytime and nighttime
            # The GPP and Reco columns remain as they are - no interpolation between data types
            
            # Step 6: Final gap-filling with ONE-WEEK max gap limit using interpolation
            # This ONLY interpolates within the same data type (daytime-only or nighttime-only)
            gap_limit = 48 * 7  # 7 days of half-hourly data (48 half-hours per day * 7 days = 336 records)
            write_log(f"Using gap-filling limit of {gap_limit} records (7 days)")
            
            # Interpolate GPP gaps (within same source type)
            if df['GPP'].notna().any():
                mask = df['GPP'].isna()
                filled = df['GPP'].interpolate(method='linear', limit=gap_limit, limit_direction='both')
                df['GPP'] = np.where(mask & filled.notna(), filled, df['GPP'])
                if mask.any() and filled.notna().any():
                    write_log(f"Interpolated GPP gaps with max {gap_limit} record gap (7 days)")
            else:
                write_log(f"WARNING: No valid GPP data found, skipping interpolation")
            
            # Interpolate Reco gaps (within same source type)
            if df['Reco'].notna().any():
                mask = df['Reco'].isna()
                filled = df['Reco'].interpolate(method='linear', limit=gap_limit, limit_direction='both')
                df['Reco'] = np.where(mask & filled.notna(), filled, df['Reco'])
                if mask.any() and filled.notna().any():
                    write_log(f"Interpolated Reco gaps with max {gap_limit} record gap (7 days)")
            else:
                write_log(f"WARNING: No valid Reco data found, skipping interpolation")
            
            # Interpolate ET gaps
            if df['ET'].notna().any():
                mask = df['ET'].isna()
                filled = df['ET'].interpolate(method='linear', limit=gap_limit, limit_direction='both')
                df['ET'] = np.where(mask & filled.notna(), filled, df['ET'])
                if mask.any() and filled.notna().any():
                    write_log(f"Interpolated ET gaps with max {gap_limit} record gap (7 days)")
            else:
                write_log(f"WARNING: No valid ET data found, skipping interpolation")
            
            # Interpolate meteorological variables (EXCLUDING WD - Wind Direction)
            # WD is excluded because it's circular/angular data and should not be linearly interpolated
            met_vars = ['PA', 'RH', 'WS', 'H_f', 'VPD_f', 'Rg_f', 'PAR_f', 'NETRAD_f']
            for col in met_vars:
                if col in df.columns and df[col].notna().any():
                    mask = df[col].isna()
                    filled = df[col].interpolate(method='linear', limit=gap_limit, limit_direction='both')
                    df[col] = np.where(mask & filled.notna(), filled, df[col])
                    if mask.any() and filled.notna().any():
                        write_log(f"Interpolated {col} gaps with max {gap_limit} record gap (7 days)")
            
            # Note: WD (Wind Direction) is NOT interpolated - left with original NaN values

            # Calculate NEP (Net Ecosystem Production)
            # NEP = GPP - Reco (positive = net carbon uptake, negative = net carbon loss)
            df['NEP'] = df['GPP'] - df['Reco']
            
            # Log NEP availability
            nep_valid = df['NEP'].notna().sum()
            total_rows = len(df)
            if nep_valid > 0:
                write_log(f"Calculated NEP for {nep_valid}/{total_rows} rows")
            else:
                write_log(f"WARNING: Could not calculate NEP - no overlapping GPP and Reco data")

            # Remove intermediate conversion columns to keep output clean
            cols_to_drop = ['GPP_DT_converted', 'Reco_DT_converted', 'GPP_nt_converted', 'Reco_nt_converted']
            df = df.drop(columns=[col for col in cols_to_drop if col in df.columns])

            # Define the final columns in desired order (without source tracking and site_name)
            final_columns = [
                # Time variables
                'DateTime', 'Year', 'month', 'day', 'Hour',
                # Processed flux variables
                'NEE', 'GPP', 'Reco', 'NEP', 'ET',
                # Original flux variables
                'NEE_f', 'LE_f', 'H_f', 'GPP_DT', 'Reco_DT', 'GPP_nt', 'Reco_nt',
                # Meteorological variables
                'Tair_f', 'PA', 'RH', 'WS', 'WD', 'VPD_f', 'Rg_f', 'PAR_f', 'NETRAD_f',
                # Intermediate variables
                'lambda'
            ]
            
            # Reorder columns if they exist in the DataFrame
            existing_columns = [col for col in final_columns if col in df.columns]
            df = df[existing_columns]
            
            # ===================================================================
            # SAVE INDIVIDUAL SITE FILE
            # ===================================================================
            # Create output filename: [original_name]_gpp_et.csv
            output_filename = f"{site_id}_gpp_et.csv"
            output_path = os.path.join(save_folder, output_filename)
            
            # Save the processed data for this site
            df.to_csv(output_path, index=False)
            write_log(f"[OK] Saved processed data to: {output_path}")
            write_log(f"     File contains {len(df)} rows and {len(df.columns)} columns")
            
            processed_sites.append(site_id)

        except Exception as e:
            write_log(f"[ERROR] Error processing {file}: {e}")
            import traceback
            write_log(traceback.format_exc())

    # Summary of all processed sites
    write_log(f"\n{'='*60}")
    write_log(f"PROCESSING SUMMARY")
    write_log(f"{'='*60}")
    write_log(f"Total sites processed: {len(processed_sites)}")
    write_log(f"Sites: {', '.join(processed_sites)}")
    write_log(f"\n[OK] All individual site files saved in: {save_folder}")
    write_log(f"Files are named as: [site_id]_gpp_et.csv")
    write_log(f"Gap-filling limit: 7 days (336 half-hourly records)")
    
    # Save master log file with UTF-8 encoding to handle special characters
    log_path = os.path.join(save_folder, 'processing_log.txt')
    with open(log_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(log_lines))
    write_log(f"\n[OK] Master processing log saved at: {log_path}")
    
    write_log(f"\nNOTE: Wind Direction (WD) was NOT interpolated due to its circular nature")
    
    return processed_sites

# ===================================================================
# MAIN EXECUTION
# ===================================================================

# Define paths
fill_folder = r"M:\Research\WUE_CUE\ameri_data\ameri_fill"
save_folder = r"M:\Research\WUE_CUE\ameri_data\ameri_ET_GPP"

# Run the processing
processed_sites = process_all_data(fill_folder, save_folder)

# Print final summary
print(f"\n{'='*60}")
print(f"PROCESSING COMPLETE")
print(f"{'='*60}")
print(f"Processed {len(processed_sites)} sites")
print(f"Input folder: {fill_folder}")
print(f"Output folder: {save_folder}")
print(f"Each file named: [site_id]_gpp_et.csv")
print(f"Gap-filling limit: 7 days (336 half-hourly records)")