# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 20:51:07 2025

@author: ammar
"""

import os
import pandas as pd
import numpy as np
###############################################################################################

## code for merging all data


def process_all_data(fill_folder, info_csv, save_folder): 
    os.makedirs(save_folder, exist_ok=True)
    log_lines = []

    def write_log(msg):
        log_lines.append(msg)
        print(msg)

    # Read site metadata
    info_df = pd.read_csv(info_csv)

    all_data = []  # To collect all processed data from each site
    csv_files = [f for f in os.listdir(fill_folder) if f.endswith('.csv')]

    for file in csv_files:
        file_path = os.path.join(fill_folder, file)
        try:
            df = pd.read_csv(file_path)

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

            # Columns needed for unit conversion
            required_cols = ['NEE_f', 'GPP_DT', 'Reco_DT', 'GPP_nt', 'Reco_nt', 'Tair_f', 'LE_f']
            for col in required_cols:
                if col not in df.columns:
                    df[col] = np.nan
                    write_log(f"Missing column '{col}' in {file}, filled with NaNs")

            # Unit conversion
            df['NEE'] = df['NEE_f'] * ((12 / 10**6) * 1800) # g of C
            df['GPP'] = df['GPP_DT'] * ((12 / 10**6) * 1800)
            df['Reco'] = df['Reco_DT'] * ((12 / 10**6) * 1800)
            df['GPP_nt'] = df['GPP_nt'] * ((12 / 10**6) * 1800)
            df['Reco_nt'] = df['Reco_nt'] * ((12 / 10**6) * 1800)
            df['lambda'] = (3149000 - 2370 * (df['Tair_f'] + 273.16)) * 1e-6
            df['ET'] = (df['LE_f'] / df['lambda']) * (1 / 1e6) * 1800 # in mm

            # Gap-filling GPP and Reco
            df['GPP'] = df['GPP'].fillna(df['GPP_nt'])
            df['Reco'] = df['Reco'].fillna(df['Reco_nt'])

            # Final gap-filling with two-week max gap limit
            gap_limit = 48 * 14
            for col in ['GPP', 'Reco', 'ET']:
                if col in df.columns:
                    mask = df[col].isna()
                    filled = df[col].interpolate(limit=gap_limit, limit_direction='both')
                    df[col] = np.where(mask & filled.notna(), filled, df[col])

            # Calculate NPP
            df['NPP'] = df['GPP'] - df['Reco']

            # Merge metadata from info_csv
            site_id = file.split('_')[0]
            match = info_df[info_df['site_name'].str.contains(site_id, case=False, na=False)]
            if not match.empty:
                for col in match.columns:
                    df[col] = match.iloc[0][col]
            else:
                write_log(f"No match found in site info for {file}")

            all_data.append(df)

        except Exception as e:
            write_log(f"❌ Error processing {file}: {e}")

    # Merge all processed data
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        final_save_path = os.path.join(save_folder, 'all_merged_data.csv')
        final_df.to_csv(final_save_path, index=False)
        write_log(f"\n✅ Final merged file saved at: {final_save_path}")
        return final_df
    else:
        write_log("⚠️ No valid data processed.")
        return pd.DataFrame()  # Return empty DataFrame if no data processed

# Example call to the function
fill_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_fill"


save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"




final_df = process_all_data(fill_folder, info_csv, save_folder)

# To check the result
print(final_df.head())  # Preview the first few rows of the final merged dataframe



####################################################################################################

## check names 
csv_file= r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\all_merged_data.csv"
df = pd.read_csv(csv_file)
print(df['site_name'].unique())


############################################################################################

# code to calculate monthly WUE by reading all merged data

# Example call
#csv_file= r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\all_merged_data.csv"
#save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"

#df = pd.read_csv(csv_file)
#print(df['site_name'].unique())

#csv_file= r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
#df = pd.read_csv(csv_file)




import pandas as pd
import numpy as np
import os

from tqdm import tqdm

def merge_growing_season_with_all(grow_csv_path, all_csv_path, output_csv_path=None):
    # Read CSVs
    grow_df = pd.read_csv(grow_csv_path)
    all_df = pd.read_csv(all_csv_path)

    # Standardize column names
    grow_df = grow_df.rename(columns={"Year": "year"}) if "Year" in grow_df.columns else grow_df
    all_df = all_df.rename(columns={"Year": "year"}) if "Year" in all_df.columns else all_df

    # Ensure datetime and extract year (keep half-hourly resolution)
    all_df['DateTime'] = pd.to_datetime(all_df['DateTime'])
    all_df['year'] = all_df['DateTime'].dt.year

    # Drop duplicates in grow_df
    grow_df = grow_df.drop_duplicates(subset=["site_name", "year"])

    # Growing season columns (everything except ID columns)
    grow_columns = grow_df.columns.difference(['site_name', 'year']).tolist()

    # Create lookup dictionary: site -> DataFrame of yearly data
    grow_dict = {
        site: group.set_index("year").sort_index()
        for site, group in grow_df.groupby("site_name")
    }

    # Counters
    exact_match = 0
    nearby_avg = 0
    fallback_avg = 0
    no_data = 0

    # Helper function to find growing season row
    def get_grow_row(site, year):
        nonlocal exact_match, nearby_avg, fallback_avg, no_data
        site_data = grow_dict.get(site)
        used_sos = pd.NA
        used_eos = pd.NA

        if site_data is None:
            no_data += 1
            return pd.Series([pd.NA] * (len(grow_columns) + 2), index=grow_columns + ['sos', 'eos'])

        if year in site_data.index:
            exact_match += 1
            row = site_data.loc[year]
            used_sos = row.get('sos', pd.NA)
            used_eos = row.get('eos', pd.NA)
            return row.reindex(grow_columns).append(pd.Series({'sos': used_sos, 'eos': used_eos}))

        # Try ±1 to ±3 years
        candidates = []
        sos_list = []
        eos_list = []
        for offset in [1, 2, 3]:
            for y in [year - offset, year + offset]:
                if y in site_data.index:
                    row = site_data.loc[y]
                    candidates.append(row)
                    sos_list.append(row.get('sos', pd.NA))
                    eos_list.append(row.get('eos', pd.NA))

        if candidates:
            nearby_avg += 1
            mean_row = pd.DataFrame(candidates).mean(numeric_only=True).reindex(grow_columns)
            used_sos = pd.to_numeric(sos_list, errors='coerce').mean()
            used_eos = pd.to_numeric(eos_list, errors='coerce').mean()
            return mean_row.append(pd.Series({'sos': used_sos, 'eos': used_eos}))

        # Fallback: site-wide average
        fallback_avg += 1
        mean_row = site_data.mean(numeric_only=True).reindex(grow_columns)
        used_sos = pd.to_numeric(site_data['sos'], errors='coerce').mean() if 'sos' in site_data else pd.NA
        used_eos = pd.to_numeric(site_data['eos'], errors='coerce').mean() if 'eos' in site_data else pd.NA
        return mean_row.append(pd.Series({'sos': used_sos, 'eos': used_eos}))

    # Apply grow values per row with tqdm progress bar
    tqdm.pandas(desc="Merging growing season data")
    grow_values = all_df.progress_apply(lambda row: get_grow_row(row['site_name'], row['year']), axis=1)

    # Concatenate and return/save
    final_df = pd.concat([all_df.reset_index(drop=True), grow_values.reset_index(drop=True)], axis=1)

    if output_csv_path:
        final_df.to_csv(output_csv_path, index=False)
        print(f"\n✅ Saved merged data to: {output_csv_path}")

    # Summary
    total = len(all_df)
    print("\n📊 Summary of growing season value sources:")
    print(f"   • Exact year match:     {exact_match}")
    print(f"   • Nearby year average:  {nearby_avg}")
    print(f"   • Site-wide fallback:   {fallback_avg}")
    print(f"   • No data at all:       {no_data}")
    print(f"   • Total rows processed: {total}")

    return final_df



grow_csv_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\phenofit_growing_season.csv"
all_csv_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\all_merged_data.csv"
#output_csv_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
output_csv_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"

merge_growing_season_with_all(grow_csv_path, all_csv_path, output_csv_path)


######################################################################################################
#sos and eos calcualted from phenofit 
# same sos and eos (taking average of all years) for a particular site.
# remove months >9 <3. so growing season march-sep
#WUE_CUE_monthly on merged data with growing season phenofit information

import os
import pandas as pd
import numpy as np



def WUE_CUE_monthly(input_csv, save_folder, info_csv):
    # === Test if save_folder is writable ===
    test_file = os.path.join(save_folder, 'test_write_permission.txt')
    try:
        with open(test_file, 'w') as f:
            f.write('permission test')
        os.remove(test_file)
    except PermissionError:
        print(f"❌ Permission denied for: {save_folder}")
        save_folder = os.path.expanduser("~/Desktop")
        print(f"✅ Falling back to local save_folder: {save_folder}")

    df = pd.read_csv(input_csv)
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['DoY'] = df['DateTime'].dt.dayofyear
    df['month'] = df['DateTime'].dt.month
    df['year'] = df['DateTime'].dt.year

    info = pd.read_csv(info_csv)

    metadata_cols = ['State', 'site_name', 'biome', 'salinity_ppt', 'salinity_fine', 'salini_coarse']
    for col in metadata_cols:
        if col not in df.columns:
            df[col] = None

    all_monthly = []
    all_filtered = []
    short_gs_issues = []

    for site in df['site_name'].unique():
        site_df = df[df['site_name'] == site].copy()
        years = site_df['year'].unique()
        site_state = site_df['State'].iloc[0]

        # --- NEW: compute average sos and eos for this site across all years ---
        avg_sos_site = site_df['avg_sos'].mean()
        avg_eos_site = site_df['avg_eos'].mean()

        valid_years = []

        for yr in years:
            yr_df = site_df[site_df['year'] == yr].copy()

            # USE SAME sos/eos FOR ALL YEARS IN THIS SITE
            sos = avg_sos_site
            eos = avg_eos_site

            year_gs = yr_df[(yr_df['DoY'] >= sos) & (yr_df['DoY'] <= eos)].copy()
            n_months = year_gs['month'].nunique()

            min_months = 3 if site_state == 'CA' else 2

            if n_months < min_months:
                # If only one year, record issue and skip
                if len(years) == 1:
                    short_gs_issues.append((site, yr))
                    continue
                else:
                    # Using site-wide average already, no fallback needed
                    # Just record issue if months still too few
                    short_gs_issues.append((site, yr))
                    continue

            # Store same avg_sos and avg_eos for this site-year
            year_gs['avg_sos'] = sos
            year_gs['avg_eos'] = eos

            valid_years.append(year_gs)

        if not valid_years:
            continue

        site_gs_df = pd.concat(valid_years)

        final_years = []
        for yr in site_gs_df['year'].unique():
            sub = site_gs_df[site_gs_df['year'] == yr]
            if not sub['NEE'].isna().all() and not sub['ET'].isna().all():
                final_years.append(yr)
        site_gs_df = site_gs_df[site_gs_df['year'].isin(final_years)]

        def is_month_valid(sub):
            return sub['NEE'].notna().mean() >= 0.5 and sub['ET'].notna().mean() >= 0.5

        valid_mask = site_gs_df.groupby(['year', 'month']).apply(is_month_valid).reset_index()
        valid_mask.columns = ['year', 'month', 'valid']
        site_gs_df = site_gs_df.merge(valid_mask, on=['year', 'month'])
        site_gs_df = site_gs_df[site_gs_df['valid']]
        site_gs_df.drop(columns=['valid'], inplace=True)

        gap_limit = 48 * 14
        for col in ['GPP', 'Reco', 'ET']:
            if col in site_gs_df.columns:
                mask = site_gs_df[col].isna()
                filled = site_gs_df[col].interpolate(limit=gap_limit, limit_direction='both')
                site_gs_df[col] = np.where(mask & filled.notna(), filled, site_gs_df[col])

        site_gs_df['NEP'] = site_gs_df['GPP'] - site_gs_df['Reco']
        all_filtered.append(site_gs_df.copy())

        numeric_cols = ['NEE', 'GPP', 'Reco', 'ET', 'NEP']
        monthly = site_gs_df.groupby(['site_name', 'year', 'month'])[numeric_cols].sum(numeric_only=True).reset_index()

        for col in metadata_cols:
            monthly[col] = site_df[col].iloc[0]

        monthly['WUE'] = monthly['GPP'] / monthly['ET']
        monthly['CUE'] = monthly['NEP'] / monthly['GPP']
        monthly = monthly[(monthly['WUE'] >= 0) & (monthly['WUE'] <= 6)]
        monthly = monthly[monthly['CUE'] >= -10]

        monthly['WUE'] = monthly['WUE'].apply(lambda x: round(x, 3) if isinstance(x, (int, float, np.float64)) else x)
        monthly['CUE'] = monthly['CUE'].apply(lambda x: round(x, 3) if isinstance(x, (int, float, np.float64)) else x)

        # === Add climate monthly averages ===
        climate_columns = ['Tair_f', 'VPD_f', 'Rg_f']
        climate_monthly_avg = site_gs_df.groupby(['site_name', 'year', 'month'])[climate_columns].mean(numeric_only=True).reset_index()
        monthly = pd.merge(monthly, climate_monthly_avg, on=['site_name', 'year', 'month'], how='left')

        # === Add avg_length unchanged per site-year ===
        monthly['avg_length'] = site_df['avg_length'].iloc[0]

        # === Add avg_sos and avg_eos as constant values per year-site ===
        # Now same for all years, so can use the site-wide averages directly
        monthly['avg_sos'] = avg_sos_site
        monthly['avg_eos'] = avg_eos_site

        # === Add site-level metadata from info_csv ===
        site_info = info[info['site_name'] == site]
        if not site_info.empty:
            for col in ['lat', 'long', 'climate', 'climate_2', 'salinity_con']:
                monthly[col] = site_info[col].iloc[0]

        all_monthly.append(monthly)

    final_df = pd.concat(all_monthly, ignore_index=True)
    all_filtered_df = pd.concat(all_filtered, ignore_index=True)

    final_df = final_df.apply(lambda x: x.map(lambda val: round(val, 3) if isinstance(val, (int, float, np.float64)) else val))
    all_filtered_df = all_filtered_df.apply(lambda x: x.map(lambda val: round(val, 3) if isinstance(val, (int, float, np.float64)) else val))

    # === REMOVE months outside the growing season (March to September) from final_df ===
    final_df = final_df[(final_df['month'] >= 3) & (final_df['month'] <= 9)].reset_index(drop=True)

    final_file = os.path.join(save_folder, 'WUE_CUE_monthly.csv')
    final_df.to_csv(final_file, index=False)

    if short_gs_issues:
        print("Sites with less than required growing season months after adjustment:")
        for site, yr in short_gs_issues:
            print(f"Site: {site}, Year: {yr}")

    return final_df, all_filtered_df


# Example usage:
input_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"

final_df, filtered_df = WUE_CUE_monthly(input_csv, save_folder, info_csv)

###################################################################################################
## yearly data set sum
import os
import pandas as pd
import numpy as np

def WUE_CUE_yearly(input_csv, save_folder, info_csv):
    # === Test save folder permissions ===
    test_file = os.path.join(save_folder, 'test_write_permission.txt')
    try:
        with open(test_file, 'w') as f:
            f.write('permission test')
        os.remove(test_file)
    except PermissionError:
        print(f"❌ Permission denied for: {save_folder}")
        save_folder = os.path.expanduser("~/Desktop")
        print(f"✅ Falling back to local save_folder: {save_folder}")

    # Read and preprocess data
    df = pd.read_csv(input_csv)
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['DoY'] = df['DateTime'].dt.dayofyear
    df['month'] = df['DateTime'].dt.month
    df['year'] = df['DateTime'].dt.year

    # Read site info
    info = pd.read_csv(info_csv)

    # Ensure metadata columns exist
    metadata_cols = ['State', 'site_name', 'biome', 'salinity_ppt', 'salinity_fine', 'salini_coarse']
    for col in metadata_cols:
        if col not in df.columns:
            df[col] = None

    # Initialize storage
    all_yearly = []
    short_gs_issues = []

    for site in df['site_name'].unique():
        site_df = df[df['site_name'] == site].copy()
        years = site_df['year'].unique()
        site_state = site_df['State'].iloc[0]

        # Calculate site-wide average SOS and EOS
        avg_sos_site = site_df['avg_sos'].mean()
        avg_eos_site = site_df['avg_eos'].mean()

        valid_years = []

        for yr in years:
            yr_df = site_df[site_df['year'] == yr].copy()
            
            # Filter to growing season
            gs_mask = (yr_df['DoY'] >= avg_sos_site) & (yr_df['DoY'] <= avg_eos_site)
            year_gs = yr_df[gs_mask].copy()
            
            # Check for sufficient data
            min_months = 3 if site_state == 'CA' else 2
            n_months = year_gs['month'].nunique()
            
            if n_months < min_months:
                short_gs_issues.append((site, yr))
                continue  # Skip this year
            
            # Add site-wide averages
            year_gs['avg_sos'] = avg_sos_site
            year_gs['avg_eos'] = avg_eos_site
            valid_years.append(year_gs)

        if not valid_years:
            continue  # Skip site if no valid years

        site_gs_df = pd.concat(valid_years)
        
        # Gap filling for key variables
        gap_limit = 48 * 14  # 14 days of half-hourly data
        for col in ['GPP', 'Reco', 'ET']:
            if col in site_gs_df.columns:
                mask = site_gs_df[col].isna()
                filled = site_gs_df[col].interpolate(limit=gap_limit, limit_direction='both')
                site_gs_df[col] = np.where(mask & filled.notna(), filled, site_gs_df[col])
        
        # Calculate NEP
        site_gs_df['NEP'] = site_gs_df['GPP'] - site_gs_df['Reco']
        
        # Group by year for yearly aggregation
        numeric_cols = ['GPP', 'Reco', 'ET', 'NEP']
        climate_cols = ['Tair_f', 'VPD_f', 'Rg_f']
        
        # Aggregation methods
        agg_methods = {col: 'sum' for col in numeric_cols}
        agg_methods.update({col: 'mean' for col in climate_cols})
        
        yearly = site_gs_df.groupby(['site_name', 'year']).agg(agg_methods).reset_index()
        
        # Calculate WUE and CUE
        yearly['WUE'] = yearly['GPP'] / yearly['ET']
        yearly['CUE'] = yearly['NEP'] / yearly['GPP']
        
        # Filter unreasonable values
        yearly = yearly[(yearly['WUE'] >= 10) & (yearly['WUE'] <= 0.1)]
        yearly = yearly[yearly['CUE'] >= -10]
        
        # Round values
        for col in numeric_cols + ['WUE', 'CUE']:
            yearly[col] = yearly[col].round(3)
        
        # Add metadata
        for col in metadata_cols:
            yearly[col] = site_df[col].iloc[0]
        
        # Add growing season characteristics
        yearly['avg_length'] = site_df['avg_length'].iloc[0]
        yearly['avg_sos'] = avg_sos_site
        yearly['avg_eos'] = avg_eos_site
        
        # Add site info
        site_info = info[info['site_name'] == site]
        if not site_info.empty:
            for col in ['lat', 'long', 'climate', 'climate_2', 'salinity_con']:
                yearly[col] = site_info[col].iloc[0]
        
        all_yearly.append(yearly)

    # Combine all sites
    if not all_yearly:
        print("⚠️ No valid yearly data created!")
        return pd.DataFrame(), pd.DataFrame()
    
    final_df = pd.concat(all_yearly, ignore_index=True)
    
    # Save results
    final_file = os.path.join(save_folder, 'WUE_CUE_yearly.csv')
    final_df.to_csv(final_file, index=False)
    
    # Report issues
    if short_gs_issues:
        print("Sites with insufficient growing season data:")
        for site, yr in short_gs_issues:
            print(f"  - {site} (Year: {yr})")
    
    return final_df

# Example usage
input_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"

yearly_df = WUE_CUE_yearly(input_csv, save_folder, info_csv)

################################################################################################################

import os
import pandas as pd
import numpy as np

def WUE_CUE_yearly(input_csv, save_folder, info_csv):
    # === Test save folder permissions ===
    test_file = os.path.join(save_folder, 'test_write_permission.txt')
    try:
        with open(test_file, 'w') as f:
            f.write('permission test')
        os.remove(test_file)
    except PermissionError:
        print(f"❌ Permission denied for: {save_folder}")
        save_folder = os.path.expanduser("~/Desktop")
        print(f"✅ Falling back to local save_folder: {save_folder}")

    # Read and preprocess data
    df = pd.read_csv(input_csv)
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['DoY'] = df['DateTime'].dt.dayofyear
    df['month'] = df['DateTime'].dt.month
    df['year'] = df['DateTime'].dt.year

    # Read site info
    info = pd.read_csv(info_csv)

    # Ensure metadata columns exist
    metadata_cols = ['State', 'site_name', 'biome', 'salinity_ppt', 'salinity_fine', 'salini_coarse']
    for col in metadata_cols:
        if col not in df.columns:
            df[col] = None

    # Initialize storage
    all_yearly = []
    short_gs_issues = []

    for site in df['site_name'].unique():
        site_df = df[df['site_name'] == site].copy()
        years = site_df['year'].unique()
        site_state = site_df['State'].iloc[0]

        # Calculate site-wide average SOS and EOS
        avg_sos_site = site_df['avg_sos'].mean()
        avg_eos_site = site_df['avg_eos'].mean()

        valid_years = []

        for yr in years:
            yr_df = site_df[site_df['year'] == yr].copy()
            
            # Filter to growing season
            gs_mask = (yr_df['DoY'] >= avg_sos_site) & (yr_df['DoY'] <= avg_eos_site)
            year_gs = yr_df[gs_mask].copy()
            
            # Check for sufficient data
            min_months = 3 if site_state == 'CA' else 2
            n_months = year_gs['month'].nunique()
            
            if n_months < min_months:
                short_gs_issues.append((site, yr))
                continue  # Skip this year
            
            # Add site-wide averages
            year_gs['avg_sos'] = avg_sos_site
            year_gs['avg_eos'] = avg_eos_site
            valid_years.append(year_gs)

        if not valid_years:
            continue  # Skip site if no valid years

        site_gs_df = pd.concat(valid_years)
        
        # CRITICAL FILTER: Remove years with no valid NEE or ET data
        final_years = []
        for yr in site_gs_df['year'].unique():
            sub = site_gs_df[site_gs_df['year'] == yr]
            if not sub['NEE'].isna().all() and not sub['ET'].isna().all():
                final_years.append(yr)
        site_gs_df = site_gs_df[site_gs_df['year'].isin(final_years)]
        
        if site_gs_df.empty:
            continue
        
        # Gap filling for key variables
        gap_limit = 48 * 14  # 14 days of half-hourly data
        for col in ['GPP', 'Reco', 'ET']:
            if col in site_gs_df.columns:
                mask = site_gs_df[col].isna()
                filled = site_gs_df[col].interpolate(limit=gap_limit, limit_direction='both')
                site_gs_df[col] = np.where(mask & filled.notna(), filled, site_gs_df[col])
        
        # Calculate NEP
        site_gs_df['NEP'] = site_gs_df['GPP'] - site_gs_df['Reco']
        
        # Group by year for yearly aggregation
        numeric_cols = ['GPP', 'Reco', 'ET', 'NEP']
        climate_cols = ['Tair_f', 'VPD_f', 'Rg_f']
        
        # Aggregation methods
        agg_methods = {col: 'sum' for col in numeric_cols}
        agg_methods.update({col: 'mean' for col in climate_cols})
        
        yearly = site_gs_df.groupby(['site_name', 'year']).agg(agg_methods).reset_index()
        
        # Calculate WUE and CUE
        yearly['WUE'] = yearly['GPP'] / yearly['ET']
        yearly['CUE'] = yearly['NEP'] / yearly['GPP']
        
        # Handle division issues
        yearly['WUE'] = yearly['WUE'].replace([np.inf, -np.inf], np.nan)
        yearly['CUE'] = yearly['CUE'].replace([np.inf, -np.inf], np.nan)
        
        # CORRECTED: Filter unreasonable values
        yearly = yearly[(yearly['WUE'] >= 0) & (yearly['WUE'] <= 6)]  # Reasonable WUE range
        yearly = yearly[yearly['CUE'] >= -10]  # CUE filter
        
        # Round values
        for col in numeric_cols + ['WUE', 'CUE']:
            yearly[col] = yearly[col].round(3)
        
        # Add metadata
        for col in metadata_cols:
            yearly[col] = site_df[col].iloc[0]
        
        # Add growing season characteristics
        yearly['avg_length'] = site_df['avg_length'].iloc[0]
        yearly['avg_sos'] = avg_sos_site
        yearly['avg_eos'] = avg_eos_site
        
        # Add site info
        site_info = info[info['site_name'] == site]
        if not site_info.empty:
            for col in ['lat', 'long', 'climate', 'climate_2', 'salinity_con']:
                yearly[col] = site_info[col].iloc[0]
        
        all_yearly.append(yearly)

    # Combine all sites
    if not all_yearly:
        print("⚠️ No valid yearly data created!")
        return pd.DataFrame()
    
    final_df = pd.concat(all_yearly, ignore_index=True)
    
    # Save results
    final_file = os.path.join(save_folder, 'WUE_CUE_yearly.csv')
    final_df.to_csv(final_file, index=False)
    
    # Report issues
    if short_gs_issues:
        print("Sites with insufficient growing season data:")
        for site, yr in short_gs_issues:
            print(f"  - {site} (Year: {yr})")
    
    return final_df

# Example usage
input_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"

yearly_df = WUE_CUE_yearly(input_csv, save_folder, info_csv)
##############################################################################################################




# function to plot yearly WUE
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Load the data
file_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_yearly.csv"
df = pd.read_csv(file_path)

# Calculate mean WUE per site and create ordered site list
site_order = df.groupby('site_name')['WUE'].mean().sort_values().index.tolist()

# Create output directory if needed
output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots"
os.makedirs(output_dir, exist_ok=True)

# Create a 1x3 grid of subplots with larger size
plt.figure(figsize=(30, 10))

# 1. WUE Box Plot
plt.subplot(1, 3, 1)
sns.boxplot(data=df, x='site_name', y='WUE', order=site_order, showmeans=True,
            meanprops={"marker":"D", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"10"})
plt.title('Water Use Efficiency (WUE) by Site', fontsize=16, fontweight='bold')
plt.xlabel('Site Name', fontsize=14)
plt.ylabel('WUE [µmol CO₂/mmol H₂O]', fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 2. GPP Box Plot - Using same site order
plt.subplot(1, 3, 2)
sns.boxplot(data=df, x='site_name', y='GPP', order=site_order, showmeans=True,
            meanprops={"marker":"D", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"10"})
plt.title('Gross Primary Production (GPP) by Site', fontsize=16, fontweight='bold')
plt.xlabel('Site Name', fontsize=14)
plt.ylabel('GPP [g C m⁻² year⁻¹]', fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 3. ET Box Plot - Using same site order
plt.subplot(1, 3, 3)
sns.boxplot(data=df, x='site_name', y='ET', order=site_order, showmeans=True,
            meanprops={"marker":"D", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"10"})
plt.title('Evapotranspiration (ET) by Site', fontsize=16, fontweight='bold')
plt.xlabel('Site Name', fontsize=14)
plt.ylabel('ET [mm H₂O year⁻¹]', fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Adjust layout and save high-res version
plt.tight_layout()
output_path = os.path.join(output_dir, "WUE_GPP_ET_boxplots_ordered.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"High-resolution plots saved to: {output_path}")

# Show plot
plt.show()

###########################################################################################################



# function for yearly max WUE 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Load the data
#file_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_yearly.csv"
file_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_yearly_max.csv"
df = pd.read_csv(file_path)

# Calculate mean WUE per site and create ordered site list
site_order = df.groupby('site_name')['WUE'].mean().sort_values().index.tolist()

# Create output directory if needed
output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots"
os.makedirs(output_dir, exist_ok=True)

# Create a 1x3 grid of subplots with larger size
plt.figure(figsize=(30, 10))

# 1. WUE Box Plot
plt.subplot(1, 3, 1)
sns.boxplot(data=df, x='site_name', y='WUE', order=site_order, showmeans=True,
            meanprops={"marker":"D", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"10"})
plt.title('Water Use Efficiency (WUE) by Site', fontsize=16, fontweight='bold')
plt.xlabel('Site Name', fontsize=14)
plt.ylabel('WUE [µmol CO₂/mmol H₂O]', fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 2. GPP Box Plot - Using same site order
plt.subplot(1, 3, 2)
sns.boxplot(data=df, x='site_name', y='GPP', order=site_order, showmeans=True,
            meanprops={"marker":"D", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"10"})
plt.title('Gross Primary Production (GPP) by Site', fontsize=16, fontweight='bold')
plt.xlabel('Site Name', fontsize=14)
plt.ylabel('GPP [g C m⁻² year⁻¹]', fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 3. ET Box Plot - Using same site order
plt.subplot(1, 3, 3)
sns.boxplot(data=df, x='site_name', y='ET', order=site_order, showmeans=True,
            meanprops={"marker":"D", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"10"})
plt.title('Evapotranspiration (ET) by Site', fontsize=16, fontweight='bold')
plt.xlabel('Site Name', fontsize=14)
plt.ylabel('ET [mm H₂O year⁻¹]', fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Adjust layout and save high-res version
plt.tight_layout()
output_path = os.path.join(output_dir, "WUE_GPP_ET_boxplots_ordered.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"High-resolution plots saved to: {output_path}")

# Show plot
plt.show()



##############################################################################################################

import os
import pandas as pd
import numpy as np

def WUE_CUE_yearly_max(input_csv, save_folder, info_csv):
    # === Test save folder permissions ===
    test_file = os.path.join(save_folder, 'test_write_permission.txt')
    try:
        with open(test_file, 'w') as f:
            f.write('permission test')
        os.remove(test_file)
    except PermissionError:
        print(f"❌ Permission denied for: {save_folder}")
        save_folder = os.path.expanduser("~/Desktop")
        print(f"✅ Falling back to local save_folder: {save_folder}")

    # Read and preprocess data
    df = pd.read_csv(input_csv)
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['DoY'] = df['DateTime'].dt.dayofyear
    df['month'] = df['DateTime'].dt.month
    df['year'] = df['DateTime'].dt.year
    df['date'] = df['DateTime'].dt.date

    # Read site info
    info = pd.read_csv(info_csv)

    # Ensure metadata columns exist
    metadata_cols = ['State', 'site_name', 'biome', 'salinity_ppt', 'salinity_fine', 'salini_coarse']
    for col in metadata_cols:
        if col not in df.columns:
            df[col] = None

    # Initialize storage
    all_max_wue = []
    short_gs_issues = []

    for site in df['site_name'].unique():
        site_df = df[df['site_name'] == site].copy()
        years = site_df['year'].unique()
        site_state = site_df['State'].iloc[0]

        # Calculate site-wide average SOS and EOS
        avg_sos_site = site_df['avg_sos'].mean()
        avg_eos_site = site_df['avg_eos'].mean()

        valid_years = []

        for yr in years:
            yr_df = site_df[site_df['year'] == yr].copy()
            
            # Filter to growing season
            gs_mask = (yr_df['DoY'] >= avg_sos_site) & (yr_df['DoY'] <= avg_eos_site)
            year_gs = yr_df[gs_mask].copy()
            
            # Check for sufficient data
            min_months = 3 if site_state == 'CA' else 2
            n_months = year_gs['month'].nunique()
            
            if n_months < min_months:
                short_gs_issues.append((site, yr))
                continue
            
            # Add site-wide averages
            year_gs['avg_sos'] = avg_sos_site
            year_gs['avg_eos'] = avg_eos_site
            valid_years.append(year_gs)

        if not valid_years:
            continue

        site_gs_df = pd.concat(valid_years)
        
        # CRITICAL FILTER: Remove years with no valid NEE or ET data
        final_years = []
        for yr in site_gs_df['year'].unique():
            sub = site_gs_df[site_gs_df['year'] == yr]
            if not sub['NEE'].isna().all() and not sub['ET'].isna().all():
                final_years.append(yr)
        site_gs_df = site_gs_df[site_gs_df['year'].isin(final_years)]
        
        if site_gs_df.empty:
            continue
        
        # Gap filling for key variables
        gap_limit = 48 * 14  # 14 days of half-hourly data
        for col in ['GPP', 'ET']:
            if col in site_gs_df.columns:
                mask = site_gs_df[col].isna()
                filled = site_gs_df[col].interpolate(limit=gap_limit, limit_direction='both')
                site_gs_df[col] = np.where(mask & filled.notna(), filled, site_gs_df[col])
        
        # Calculate daily sums for GPP and ET
        daily_sums = site_gs_df.groupby(['site_name', 'year', 'date']).agg({
            'GPP': 'sum',
            'ET': 'sum'
        }).reset_index()
        
        # Calculate daily WUE and handle division by zero
        daily_sums['WUE'] = daily_sums['GPP'] / daily_sums['ET']
        daily_sums['WUE'] = daily_sums['WUE'].replace([np.inf, -np.inf], np.nan)
        
        # Find maximum WUE for each year
        daily_sums['max_wue'] = daily_sums.groupby(['site_name', 'year'])['WUE'].transform('max')
        yearly_max = daily_sums[daily_sums['WUE'] == daily_sums['max_wue']]
        
        # For years with multiple max values, take the first occurrence
        yearly_max = yearly_max.drop_duplicates(subset=['site_name', 'year'], keep='first')
        
        # Add day of year
        yearly_max['DoY'] = pd.to_datetime(yearly_max['date']).dt.dayofyear
        
        # Round WUE values to 3 decimal places
        yearly_max['WUE'] = yearly_max['WUE'].round(3)
        
        # Add metadata
        for col in metadata_cols:
            yearly_max[col] = site_df[col].iloc[0]
        
        # Add growing season characteristics
        yearly_max['avg_length'] = site_df['avg_length'].iloc[0]
        yearly_max['avg_sos'] = avg_sos_site
        yearly_max['avg_eos'] = avg_eos_site
        
        # Add site info
        site_info = info[info['site_name'] == site]
        if not site_info.empty:
            for col in ['lat', 'long', 'climate', 'climate_2', 'salinity_con']:
                yearly_max[col] = site_info[col].iloc[0]
        
        # Drop temporary column
        yearly_max = yearly_max.drop(columns=['max_wue'])
        
        all_max_wue.append(yearly_max)

    # Combine all sites
    if not all_max_wue:
        print("⚠️ No valid yearly max WUE data created!")
        return pd.DataFrame()
    
    final_df = pd.concat(all_max_wue, ignore_index=True)
    
    # NEW FILTER: Remove unreasonable WUE values
    # Remove zero, infinity, and values outside reasonable range
    final_df = final_df[
        (final_df['WUE'] > 0) &  # Remove zero values
        (final_df['WUE'] >= 0.06) &  # Minimum reasonable WUE
        (final_df['WUE'] <= 90) &    # Maximum reasonable WUE
        (final_df['WUE'].notna())    # Remove NaN values
    ]
    
    # Save results
    final_file = os.path.join(save_folder, 'WUE_CUE_yearly_max.csv')
    final_df.to_csv(final_file, index=False)
    
    # Report issues
    if short_gs_issues:
        print("Sites with insufficient growing season data:")
        for site, yr in short_gs_issues:
            print(f"  - {site} (Year: {yr})")
    
    return final_df

# Example usage
input_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"

max_wue_df = WUE_CUE_yearly_max(input_csv, save_folder, info_csv)

################################################################################################################
# different growing season for each year in same site
#WUE_CUE_monthly on merged data with growing season phenofit information

input_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"

import pandas as pd
import numpy as np
import os



import pandas as pd
import numpy as np
import os

def WUE_CUE_monthly(input_csv, save_folder, info_csv):
    # === Test if save_folder is writable ===
    test_file = os.path.join(save_folder, 'test_write_permission.txt')
    try:
        with open(test_file, 'w') as f:
            f.write('permission test')
        os.remove(test_file)
    except PermissionError:
        print(f"❌ Permission denied for: {save_folder}")
        save_folder = os.path.expanduser("~/Desktop")
        print(f"✅ Falling back to local save_folder: {save_folder}")

    df = pd.read_csv(input_csv)
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['DoY'] = df['DateTime'].dt.dayofyear
    df['month'] = df['DateTime'].dt.month
    df['year'] = df['DateTime'].dt.year

    info = pd.read_csv(info_csv)

    metadata_cols = ['State', 'site_name', 'biome', 'salinity_ppt', 'salinity_fine', 'salini_coarse']
    for col in metadata_cols:
        if col not in df.columns:
            df[col] = None

    all_monthly = []
    all_filtered = []
    short_gs_issues = []

    for site in df['site_name'].unique():
        site_df = df[df['site_name'] == site].copy()
        years = site_df['year'].unique()
        site_state = site_df['State'].iloc[0]
        valid_years = []

        for yr in years:
            yr_df = site_df[site_df['year'] == yr].copy()
            sos = yr_df['avg_sos'].iloc[0]
            eos = yr_df['avg_eos'].iloc[0]
            year_gs = yr_df[(yr_df['DoY'] >= sos) & (yr_df['DoY'] <= eos)].copy()
            n_months = year_gs['month'].nunique()

            min_months = 3 if site_state == 'CA' else 2

            if n_months < min_months:
                if len(years) == 1:
                    short_gs_issues.append((site, yr))
                    continue
                else:
                    other_df = site_df[site_df['year'] != yr]
                    avg_sos = other_df['avg_sos'].mean()
                    avg_eos = other_df['avg_eos'].mean()
                    sos = avg_sos  # updated sos
                    eos = avg_eos  # updated eos
                    year_gs = yr_df[(yr_df['DoY'] >= sos) & (yr_df['DoY'] <= eos)].copy()
                    n_months = year_gs['month'].nunique()

                    if n_months < min_months:
                        short_gs_issues.append((site, yr))
                        continue

            # Store updated sos/eos for later use
            year_gs['avg_sos'] = sos
            year_gs['avg_eos'] = eos

            valid_years.append(year_gs)

        if not valid_years:
            continue

        site_gs_df = pd.concat(valid_years)

        final_years = []
        for yr in site_gs_df['year'].unique():
            sub = site_gs_df[site_gs_df['year'] == yr]
            if not sub['NEE'].isna().all() and not sub['ET'].isna().all():
                final_years.append(yr)
        site_gs_df = site_gs_df[site_gs_df['year'].isin(final_years)]

        def is_month_valid(sub):
            return sub['NEE'].notna().mean() >= 0.5 and sub['ET'].notna().mean() >= 0.5

        valid_mask = site_gs_df.groupby(['year', 'month']).apply(is_month_valid).reset_index()
        valid_mask.columns = ['year', 'month', 'valid']
        site_gs_df = site_gs_df.merge(valid_mask, on=['year', 'month'])
        site_gs_df = site_gs_df[site_gs_df['valid']]
        site_gs_df.drop(columns=['valid'], inplace=True)

        gap_limit = 48 * 14
        for col in ['GPP', 'Reco', 'ET']:
            if col in site_gs_df.columns:
                mask = site_gs_df[col].isna()
                filled = site_gs_df[col].interpolate(limit=gap_limit, limit_direction='both')
                site_gs_df[col] = np.where(mask & filled.notna(), filled, site_gs_df[col])

        site_gs_df['NEP'] = site_gs_df['GPP'] - site_gs_df['Reco']
        all_filtered.append(site_gs_df.copy())

        numeric_cols = ['NEE', 'GPP', 'Reco', 'ET', 'NEP']
        monthly = site_gs_df.groupby(['site_name', 'year', 'month'])[numeric_cols].sum(numeric_only=True).reset_index()

        for col in metadata_cols:
            monthly[col] = site_df[col].iloc[0]

        monthly['WUE'] = monthly['GPP'] / monthly['ET']
        monthly['CUE'] = monthly['NEP'] / monthly['GPP']
        monthly = monthly[(monthly['WUE'] >= 0) & (monthly['WUE'] <= 6)]
        monthly = monthly[monthly['CUE'] >= -10]

        monthly['WUE'] = monthly['WUE'].apply(lambda x: round(x, 3) if isinstance(x, (int, float, np.float64)) else x)
        monthly['CUE'] = monthly['CUE'].apply(lambda x: round(x, 3) if isinstance(x, (int, float, np.float64)) else x)

        # === Add climate monthly averages ===
        climate_columns = ['Tair_f', 'VPD_f', 'Rg_f']
        climate_monthly_avg = site_gs_df.groupby(['site_name', 'year', 'month'])[climate_columns].mean(numeric_only=True).reset_index()
        monthly = pd.merge(monthly, climate_monthly_avg, on=['site_name', 'year', 'month'], how='left')

        # === Add avg_length unchanged per site-year ===
        monthly['avg_length'] = site_df['avg_length'].iloc[0]

        # === Add avg_sos and avg_eos as constant values per year-site ===
        sos_eos_df = site_gs_df.groupby(['site_name', 'year'])[['avg_sos', 'avg_eos']].first().reset_index()
        monthly = pd.merge(monthly, sos_eos_df, on=['site_name', 'year'], how='left')

        # === Add site-level metadata from info_csv ===
        site_info = info[info['site_name'] == site]
        if not site_info.empty:
            for col in ['lat', 'long', 'climate', 'climate_2', 'salinity_con']:
                monthly[col] = site_info[col].iloc[0]

        all_monthly.append(monthly)

    final_df = pd.concat(all_monthly, ignore_index=True)
    all_filtered_df = pd.concat(all_filtered, ignore_index=True)

    final_df = final_df.apply(lambda x: x.map(lambda val: round(val, 3) if isinstance(val, (int, float, np.float64)) else val))
    all_filtered_df = all_filtered_df.apply(lambda x: x.map(lambda val: round(val, 3) if isinstance(val, (int, float, np.float64)) else val))

    final_file = os.path.join(save_folder, 'WUE_CUE_monthly.csv')
    final_df.to_csv(final_file, index=False)

    if short_gs_issues:
        print("Sites with less than required growing season months after adjustment:")
        for site, yr in short_gs_issues:
            print(f"Site: {site}, Year: {yr}")

    return final_df, all_filtered_df


input_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\merged_with_grow_fill.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"

final_df, filtered_df=WUE_CUE_monthly(input_csv, save_folder, info_csv)


## end of automation

##########################################################################################################
# code to calculate monthly WUE by reading one csv by one. not merged data
################################################################################################################3

import os
import pandas as pd
import numpy as np

def WUE_CUE_monthly(fill_folder, info_csv, save_folder):
    os.makedirs(save_folder, exist_ok=True)
    log_lines = []

    def write_log(msg):
        log_lines.append(msg)
        print(msg)

    # Read site metadata
    info_df = pd.read_csv(info_csv)

    all_monthly = []  # To collect monthly summaries from each site
    csv_files = [f for f in os.listdir(fill_folder) if f.endswith('.csv')]

    for file in csv_files:
        file_path = os.path.join(fill_folder, file)
        try:
            df = pd.read_csv(file_path)

            # Handle time
            if 'DateTime' in df.columns:
                df['DateTime'] = pd.to_datetime(df['DateTime'], errors='coerce')
                df['month'] = df['DateTime'].dt.month
                df['day'] = df['DateTime'].dt.day
                df['Hour'] = df['DateTime'].dt.hour + df['DateTime'].dt.minute / 60
            else:
                write_log(f"Missing 'DateTime' column in {file}")
                continue

            # Columns needed for unit conversion
            required_cols = ['NEE_f', 'GPP_DT', 'Reco_DT', 'GPP_nt', 'Reco_nt', 'Tair_f', 'LE_f']
            for col in required_cols:
                if col not in df.columns:
                    df[col] = np.nan
                    write_log(f"Missing column '{col}' in {file}, filled with NaNs")

            # Unit conversion
            df['NEE'] = df['NEE_f'] * ((12 / 10**6) * 1800)
            df['GPP'] = df['GPP_DT'] * ((12 / 10**6) * 1800)
            df['Reco'] = df['Reco_DT'] * ((12 / 10**6) * 1800)
            df['GPP_nt'] = df['GPP_nt'] * ((12 / 10**6) * 1800)
            df['Reco_nt'] = df['Reco_nt'] * ((12 / 10**6) * 1800)
            df['lambda'] = (3149000 - 2370 * (df['Tair_f'] + 273.16)) * 1e-6
            df['ET'] = (df['LE_f'] / df['lambda']) * (1 / 1e6) * 1800

            # Gap-filling GPP and Reco
            df['GPP'] = df['GPP'].fillna(df['GPP_nt'])
            df['Reco'] = df['Reco'].fillna(df['Reco_nt'])

            # Growing season filtering (May–Aug)
            growing_df = df[df['month'].isin([5, 6, 7, 8])]

            # Remove full-year missing data
            valid_years = []
            for yr in growing_df['Year'].unique():
                subset = growing_df[growing_df['Year'] == yr]
                if not subset['NEE'].isna().all() and not subset['LE_f'].isna().all():
                    valid_years.append(yr)
            growing_df = growing_df[growing_df['Year'].isin(valid_years)]

            # Remove months with >90% missing data
            def is_month_valid(sub):
                return sub['NEE'].notna().mean() >= 0.1 and sub['LE_f'].notna().mean() >= 0.1

            valid_mask = growing_df.groupby(['Year', 'month']).apply(is_month_valid).reset_index()
            valid_mask.columns = ['Year', 'month', 'valid']
            growing_df = growing_df.merge(valid_mask, on=['Year', 'month'])
            growing_df = growing_df[growing_df['valid']]
            growing_df.drop(columns=['valid'], inplace=True)

            # Final gap-filling with two-week max gap limit
            gap_limit = 48 * 14
            for col in ['GPP', 'Reco', 'ET']:
                if col in growing_df.columns:
                    mask = growing_df[col].isna()
                    filled = growing_df[col].interpolate(limit=gap_limit, limit_direction='both')
                    growing_df[col] = np.where(mask & filled.notna(), filled, growing_df[col])

            # Calculate NPP
            growing_df['NPP'] = growing_df['GPP'] - growing_df['Reco']

            # Monthly group and calculation
            monthly = growing_df.groupby(['Year', 'month']).sum(numeric_only=True).reset_index()
            monthly['WUE'] = monthly['GPP'] / monthly['ET']
            monthly['CUE'] = monthly['NPP'] / monthly['GPP']

            # Merge metadata from info_csv
            site_id = file.split('_')[0]
            match = info_df[info_df['site_name'].str.contains(site_id, case=False, na=False)]
            if not match.empty:
                for col in match.columns:
                    monthly[col] = match.iloc[0][col]
            else:
                write_log(f"No match found in site info for {file}")

            # Drop unnecessary columns
            monthly.drop(columns=[col for col in ['day', 'DoY', 'Hour'] if col in monthly.columns], inplace=True)

            # Round numeric values to 3 significant digits
            monthly = monthly.apply(lambda x: x.map(lambda val: round(val, 3) if isinstance(val, (int, float, np.float64)) else val))

            all_monthly.append(monthly)

        except Exception as e:
            write_log(f"❌ Error processing {file}: {e}")

    # Merge all monthly results
    if all_monthly:
        final_df = pd.concat(all_monthly, ignore_index=True)
        final_save_path = os.path.join(save_folder, 'WUE_CUE_monthly.csv')
        final_df.to_csv(final_save_path, index=False)
        write_log(f"\n✅ Final merged file saved at: {final_save_path}")
        return final_df
    else:
        write_log("⚠️ No valid data processed.")
        return pd.DataFrame()  # Return empty DataFrame if no data processed


# Example call
fill_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_fill"
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"

df = WUE_CUE_monthly(fill_folder, info_csv, save_folder)
