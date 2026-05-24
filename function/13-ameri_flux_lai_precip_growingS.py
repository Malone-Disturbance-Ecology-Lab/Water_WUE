import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta

# ============================================================================
# CONFIGURATION
# ============================================================================
# Input paths
AMERI_DATA_PATH = r"M:\Research\WUE_CUE\ameri_data\ameri_ET_GPP"
PRECIP_PATH = r"M:\Research\WUE_CUE\drivers\ameri_drivers\precip"
LAI_PATH = r"M:\Research\WUE_CUE\drivers\ameri_drivers\lai\growing_season_LAI"
OUTPUT_PATH = r"M:\Research\WUE_CUE\ameri_data\ameri_fill__lai_precip"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_PATH, exist_ok=True)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_ameriflux_data(file_path):
    """
    Load Ameriflux half-hourly data
    
    Units for key variables:
    - ET: mm per 30 minutes
    - GPP, Reco, NEP: g C m⁻² per 30 minutes
    - NEE: µmol CO₂ m⁻² s⁻¹ (typical, but keep as is)
    """
    df = pd.read_csv(file_path)
    
    # Create datetime column from DateTime string
    # Format: "1/1/2014 0:00" or "1/1/2014 0:30"
    df['datetime'] = pd.to_datetime(df['DateTime'])
    
    # Calculate DOY (Day of Year) from datetime (1-366)
    df['DOY'] = df['datetime'].dt.dayofyear
    
    # Extract year from datetime
    df['Year'] = df['datetime'].dt.year
    
    # Keep original day column (day of month) but don't use for filtering
    # The 'day' column exists but is day of month (1-31), not DOY
    
    return df

def load_precip_data(site_name, precip_path):
    """
    Load daily precipitation data for a specific site
    
    Precipitation units: mm per day
    """
    # Find precip file for this site
    precip_files = glob(os.path.join(precip_path, f"*{site_name}*.csv"))
    
    if not precip_files:
        # Try without wildcards - exact site name
        precip_files = glob(os.path.join(precip_path, f"{site_name}.csv"))
    
    if not precip_files:
        print(f"Warning: No precip file found for site {site_name}")
        return None
    
    precip_df = pd.read_csv(precip_files[0])
    
    # Standardize column names
    if 'date' in precip_df.columns:
        precip_df['date'] = pd.to_datetime(precip_df['date'])
    elif 'Date' in precip_df.columns:
        precip_df['date'] = pd.to_datetime(precip_df['Date'])
    
    # Ensure precip_mm column exists
    precip_col = [col for col in precip_df.columns if 'precip' in col.lower() or 'ppt' in col.lower()]
    if precip_col:
        precip_df['precip_mm'] = precip_df[precip_col[0]]
    else:
        print(f"Warning: No precipitation column found for site {site_name}")
        return None
    
    return precip_df[['date', 'precip_mm']]

def resample_precip_to_halfhourly(precip_daily, ameri_dates):
    """
    Resample daily precipitation to half-hourly values
    
    Daily precipitation (mm/day) is distributed evenly across half-hour periods
    Each half-hour gets precip_mm / 48 (since 24 hours * 2 half-hours = 48)
    
    Units: mm per 30 minutes
    
    Note: This assumes uniform precipitation distribution throughout the day.
    """
    if precip_daily is None:
        return None
    
    # Create date range for half-hourly timestamps matching Ameriflux
    start_date = ameri_dates.min()
    end_date = ameri_dates.max()
    
    # Generate half-hourly timestamps
    half_hourly_dates = pd.date_range(start=start_date, end=end_date, freq='30T')
    
    # Create DataFrame with half-hourly timestamps
    precip_halfhourly = pd.DataFrame({'datetime': half_hourly_dates})
    
    # Add date column for merging
    precip_halfhourly['date'] = precip_halfhourly['datetime'].dt.date
    
    # Prepare daily precip for merging
    precip_daily['date'] = precip_daily['date'].dt.date
    
    # Merge daily precip to half-hourly
    precip_halfhourly = precip_halfhourly.merge(
        precip_daily[['date', 'precip_mm']], 
        on='date', 
        how='left'
    )
    
    # Distribute daily precip across half-hours (divide by 48 half-hours per day)
    # Units: precip_mm per half-hour = (mm/day) / 48
    precip_halfhourly['precip_mm_per_30min'] = precip_halfhourly['precip_mm'] / 48
    
    # Fill NaN with 0 (no precipitation)
    precip_halfhourly['precip_mm_per_30min'] = precip_halfhourly['precip_mm_per_30min'].fillna(0)
    
    return precip_halfhourly[['datetime', 'precip_mm_per_30min']]

def load_lai_data(site_name, lai_path):
    """
    Load LAI data for a specific site
    
    LAI units: m² m⁻² (leaf area index)
    Data available weekly/biweekly during growing season only
    
    Returns:
    - lai_df: DataFrame with LAI measurements
    - sos_eos_df: DataFrame with start and end of growing season per year
    """
    lai_files = glob(os.path.join(lai_path, f"*{site_name}*.csv"))
    
    if not lai_files:
        print(f"Warning: No LAI file found for site {site_name}")
        return None, None
    
    lai_df = pd.read_csv(lai_files[0])
    
    # Standardize column names
    lai_df['Date'] = pd.to_datetime(lai_df['Date'])
    lai_df['DOY'] = lai_df['Date'].dt.dayofyear
    
    # Extract avg_sos and avg_eos (start and end of growing season in DOY)
    if 'avg_sos' in lai_df.columns and 'avg_eos' in lai_df.columns:
        # Get unique SOS/EOS per year
        lai_df['Year'] = lai_df['Date'].dt.year
        sos_eos_df = lai_df.groupby('Year')[['avg_sos', 'avg_eos']].first().reset_index()
        sos_eos_df = sos_eos_df.sort_values('Year')
    else:
        print(f"Warning: SOS/EOS columns not found for site {site_name}")
        return None, None
    
    return lai_df[['Date', 'DOY', 'lai', 'avg_sos', 'avg_eos']], sos_eos_df

def get_sos_eos_for_year(year, sos_eos_df):
    """
    Get SOS and EOS for a given year.
    If year not available, use next available year.
    If next not available, use previous available year.
    
    Returns:
    - sos, eos: Start and end of growing season in DOY
    - used_year: The year that was actually used (for logging)
    """
    # Check if year exists in SOS/EOS data
    if year in sos_eos_df['Year'].values:
        year_data = sos_eos_df[sos_eos_df['Year'] == year].iloc[0]
        return int(year_data['avg_sos']), int(year_data['avg_eos']), year
    
    # If not found, try to find next available year
    future_years = sos_eos_df[sos_eos_df['Year'] > year].sort_values('Year')
    if len(future_years) > 0:
        used_year = future_years.iloc[0]['Year']
        sos = int(future_years.iloc[0]['avg_sos'])
        eos = int(future_years.iloc[0]['avg_eos'])
        return sos, eos, used_year
    
    # If no future year, try previous available year
    past_years = sos_eos_df[sos_eos_df['Year'] < year].sort_values('Year', ascending=False)
    if len(past_years) > 0:
        used_year = past_years.iloc[0]['Year']
        sos = int(past_years.iloc[0]['avg_sos'])
        eos = int(past_years.iloc[0]['avg_eos'])
        return sos, eos, used_year
    
    # No SOS/EOS data at all
    return None, None, None

def filter_growing_season(ameri_df, sos_eos):
    """
    Filter Ameriflux data to only include days within growing season
    based on SOS (Start of Season) and EOS (End of Season) from LAI data
    
    For years without SOS/EOS data, uses next available year.
    If next not available, uses previous available year.
    
    CRITICAL: Uses DOY (Day of Year, 1-366) from datetime, NOT the 'day' column
    which is day of month (1-31)
    
    Returns:
    - filtered_df: DataFrame with growing season data including avg_sos and avg_eos columns
    """
    if sos_eos is None or len(sos_eos) == 0:
        print("Warning: No SOS/EOS data available, returning empty DataFrame")
        return pd.DataFrame()
    
    # Create a boolean mask for growing season
    growing_season_mask = pd.Series(False, index=ameri_df.index)
    
    # For each year in the Ameriflux data, find appropriate SOS/EOS
    for year in ameri_df['Year'].unique():
        sos, eos, used_year = get_sos_eos_for_year(year, sos_eos)
        
        if sos is None:
            print(f"  WARNING: No SOS/EOS data found for year {year} and no nearby years available")
            continue
        
        # Apply filter using DOY
        year_mask = (ameri_df['Year'] == year) & (ameri_df['DOY'] >= sos) & (ameri_df['DOY'] <= eos)
        growing_season_mask = growing_season_mask | year_mask
    
    # Filter to growing season only
    ameri_growing = ameri_df[growing_season_mask].copy()
    
    # Add avg_sos and avg_eos columns to the filtered DataFrame
    # Create arrays for SOS and EOS that match the filtered rows
    avg_sos_values = []
    avg_eos_values = []
    
    for idx in ameri_growing.index:
        year = ameri_growing.loc[idx, 'Year']
        sos, eos, used_year = get_sos_eos_for_year(year, sos_eos)
        avg_sos_values.append(sos)
        avg_eos_values.append(eos)
    
    ameri_growing['avg_sos'] = avg_sos_values
    ameri_growing['avg_eos'] = avg_eos_values
    
    # Report filtering results
    if len(ameri_growing) > 0:
        print(f"  - Filtered to growing season using DOY: {len(ameri_growing)} records "
              f"({len(ameri_growing)/len(ameri_df)*100:.1f}% of original data)")
        
        # Report unique SOS/EOS pairs used
        sos_eos_pairs = ameri_growing[['Year', 'avg_sos', 'avg_eos']].drop_duplicates().sort_values('Year')
        print(f"  - SOS/EOS values used by year:")
        for _, row in sos_eos_pairs.iterrows():
            print(f"      Year {int(row['Year'])}: SOS={int(row['avg_sos'])}, EOS={int(row['avg_eos'])}")
    else:
        print(f"  - No records match growing season criteria")
    
    return ameri_growing

def interpolate_lai(ameri_df, lai_df, sos_eos):
    """
    Interpolate LAI values to half-hourly resolution for growing season only
    
    Uses linear interpolation between available LAI measurements
    For years without LAI data, uses nearest available year's LAI pattern
    Assumes ameri_df is already filtered to growing season
    """
    if lai_df is None or sos_eos is None:
        return None
    
    # Create a copy to avoid warnings
    ameri_with_lai = ameri_df.copy()
    ameri_with_lai['lai'] = np.nan
    
    # Get all available years with LAI data
    available_lai_years = sorted(lai_df['Date'].dt.year.unique())
    
    # For each year in the dataset
    for year in ameri_with_lai['Year'].unique():
        # Get SOS and EOS for this year (using the same mapping logic)
        sos, eos, used_sos_year = get_sos_eos_for_year(year, sos_eos)
        
        if sos is None:
            continue
        
        # Find LAI data for this year
        lai_year = lai_df[lai_df['Date'].dt.year == year].copy()
        
        # If no LAI data for this year, use nearest year with LAI data
        if len(lai_year) == 0:
            # Find closest year with LAI data
            if len(available_lai_years) > 0:
                closest_year = min(available_lai_years, key=lambda x: abs(x - year))
                lai_year = lai_df[lai_df['Date'].dt.year == closest_year].copy()
                if len(lai_year) > 0:
                    print(f"    Year {year}: using LAI data from year {closest_year}")
        
        if len(lai_year) > 0:
            # Create daily LAI interpolated values
            lai_year = lai_year.sort_values('DOY')
            
            # Create a continuous DOY range for this year's growing season
            doy_range = np.arange(sos, eos + 1)
            
            # Interpolate LAI for each DOY
            if len(lai_year) >= 2:
                # Linear interpolation between available measurements
                lai_interp = np.interp(
                    doy_range, 
                    lai_year['DOY'].values, 
                    lai_year['lai'].values
                )
                lai_daily = pd.DataFrame({
                    'DOY': doy_range,
                    'lai_daily': lai_interp
                })
            elif len(lai_year) == 1:
                # Constant LAI if only one measurement
                lai_daily = pd.DataFrame({
                    'DOY': doy_range,
                    'lai_daily': lai_year['lai'].values[0]
                })
            else:
                continue
            
            # Merge daily LAI to half-hourly
            for idx, row in lai_daily.iterrows():
                mask = (ameri_with_lai['DOY'] == row['DOY']) & (ameri_with_lai['Year'] == year)
                ameri_with_lai.loc[mask, 'lai'] = row['lai_daily']
    
    # Fill any remaining NaN values (should be minimal within growing season)
    ameri_with_lai['lai'] = ameri_with_lai['lai'].fillna(method='ffill').fillna(method='bfill').fillna(0)
    
    return ameri_with_lai

# ============================================================================
# MAIN PROCESSING LOOP
# ============================================================================

def process_all_sites():
    """Process all Ameriflux sites and attach precip and LAI data, keeping only growing season"""
    
    # Get all Ameriflux CSV files
    ameri_files = glob(os.path.join(AMERI_DATA_PATH, "*.csv"))
    
    print(f"Found {len(ameri_files)} Ameriflux files")
    
    for ameri_file in ameri_files:
        # Extract site name from filename
        file_basename = os.path.basename(ameri_file)
        # Pattern: US-A03_fill_gpp_et.csv -> US-A03
        site_name = file_basename.split('_')[0]
        
        print(f"\n{'='*60}")
        print(f"Processing site: {site_name}")
        print(f"{'='*60}")
        
        # Step 1: Load Ameriflux data and calculate DOY properly
        ameri_df = load_ameriflux_data(ameri_file)
        print(f"  - Loaded Ameriflux data: {len(ameri_df)} half-hourly records")
        print(f"  - Date range: {ameri_df['datetime'].min()} to {ameri_df['datetime'].max()}")
        print(f"  - Years present: {sorted(ameri_df['Year'].unique())}")
        
        # Step 2: Load LAI data to get SOS/EOS for filtering
        lai_df, sos_eos = load_lai_data(site_name, LAI_PATH)
        
        if sos_eos is None or len(sos_eos) == 0:
            print(f"  - No SOS/EOS data for site {site_name}, skipping...")
            continue
        
        print(f"  - Available SOS/EOS years: {sorted(sos_eos['Year'].unique())}")
        
        # Step 3: Filter to growing season using DOY with fallback logic
        # This now returns the DataFrame WITH avg_sos and avg_eos columns added
        ameri_growing = filter_growing_season(ameri_df, sos_eos)
        
        if len(ameri_growing) == 0:
            print(f"  - No growing season data found for site {site_name}, skipping...")
            continue
        
        # Step 4: Load and attach precipitation (only for growing season dates)
        precip_daily = load_precip_data(site_name, PRECIP_PATH)
        if precip_daily is not None:
            precip_halfhourly = resample_precip_to_halfhourly(precip_daily, ameri_growing['datetime'])
            # Merge precip with Ameriflux data
            ameri_growing = ameri_growing.merge(precip_halfhourly, on='datetime', how='left')
            print(f"  - Attached precipitation data (units: mm per 30 minutes)")
            print(f"    Note: Daily precip (mm/day) divided by 48 to get mm/30min")
        else:
            ameri_growing['precip_mm_per_30min'] = 0
            print(f"  - No precipitation data found, setting to 0")
        
        # Step 5: Attach LAI data with interpolation
        if lai_df is not None:
            ameri_growing = interpolate_lai(ameri_growing, lai_df, sos_eos)
            print(f"  - Attached LAI data (units: m² m⁻²) with linear interpolation")
        else:
            ameri_growing['lai'] = np.nan
            print(f"  - No LAI data found")
        
        # Step 6: Save the enriched dataset (growing season only) - SIMPLE SITE NAME
        output_file = os.path.join(OUTPUT_PATH, f"{site_name}.csv")
        
        # Select columns to save (keep all original plus new ones)
        # Now includes avg_sos and avg_eos (original names from LAI data)
        columns_to_save = ['datetime', 'Year', 'month', 'day', 'DOY', 'Hour', 
                          'NEE', 'GPP', 'Reco', 'NEP', 'ET',
                          'NEE_f', 'LE_f', 'H_f', 'Tair_f', 'PA', 
                          'RH', 'WS', 'WD', 'VPD_f', 'Rg_f', 'PAR_f', 
                          'NETRAD_f', 'lambda',
                          'precip_mm_per_30min', 'lai',
                          'avg_sos', 'avg_eos']
        
        # Check which columns actually exist
        existing_cols = [col for col in columns_to_save if col in ameri_growing.columns]
        
        # Save to CSV
        ameri_growing[existing_cols].to_csv(output_file, index=False)
        
        print(f"  - Saved to: {output_file}")
        print(f"  - Records: {len(ameri_growing)} (growing season only)")
        print(f"  - Variables: {len(existing_cols)}")
        print(f"  - Growing season DOY range: {ameri_growing['DOY'].min()} to {ameri_growing['DOY'].max()}")
        print(f"  - SOS/EOS columns added: avg_sos, avg_eos (original LAI file column names)")

# ============================================================================
# RUN THE PROCESSING
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("AMERIFLUX DATA ENRICHMENT WITH PRECIPITATION AND LAI")
    print("GROWING SEASON ONLY OUTPUT WITH FALLBACK YEAR MAPPING")
    print("=" * 80)
    print("\nUNITS NOTE:")
    print("  - ET: mm per 30 minutes")
    print("  - GPP, Reco, NEP: g C m⁻² per 30 minutes")
    print("  - Precipitation: mm per 30 minutes (resampled from daily mm/day ÷ 48)")
    print("  - LAI: m² m⁻² (interpolated to half-hourly)")
    print("\nFILTERING NOTE:")
    print("  - Output includes ONLY growing season data")
    print("  - Growing season defined by DOY (Day of Year) between avg_sos and avg_eos")
    print("  - DOY is calculated from DateTime column, NOT from the 'day' column")
    print("  - For years without SOS/EOS: uses next available year, then previous year")
    print("  - For years without LAI data: uses nearest available year's LAI pattern")
    print("\nCOLUMNS ADDED:")
    print("  - avg_sos: Start of growing season (DOY) used for filtering (original LAI file column name)")
    print("  - avg_eos: End of growing season (DOY) used for filtering (original LAI file column name)")
    print("=" * 80)
    
    process_all_sites()
    
    print("\n" + "=" * 80)
    print("PROCESSING COMPLETE!")
    print(f"Output saved to: {OUTPUT_PATH}")
    print("All files contain only growing season data filtered by DOY")
    print("Each file includes avg_sos and avg_eos columns (original LAI file names)")
    print("=" * 80)