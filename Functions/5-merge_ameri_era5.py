################################################################################################################


import os
import pandas as pd
import xarray as xr
import pytz
from timezonefinder import TimezoneFinder

def merge_ameri_era5(ameri_path, era5_path, ameri_file):
    """
    Merges AmeriFlux data with ERA5 data after converting UTC timestamps to local time.
    
    Parameters:
        ameri_path (str): Path to the AmeriFlux data directory.
        era5_path (str): Path to the ERA5 data directory.
        ameri_file (str): Name of the AmeriFlux CSV file.
    
    Returns:
        pd.DataFrame: Merged DataFrame with AmeriFlux and ERA5 data.
    """
    # Construct full path for AmeriFlux file
    ameri_file_path = os.path.join(ameri_path, ameri_file)
    
    # Read AmeriFlux data
    ameri = pd.read_csv(ameri_file_path)
    ameri["TIMESTAMP"] = pd.to_datetime(ameri["DateTime"])
    
    # Load ERA5 data
    accum_files = sorted([os.path.join(era5_path, f) for f in os.listdir(era5_path) if 'data_stream-oper_stepType-accum' in f])
    instant_files = sorted([os.path.join(era5_path, f) for f in os.listdir(era5_path) if 'data_stream-oper_stepType-instant' in f])
    
    accum_data = xr.open_mfdataset(accum_files, combine='nested', concat_dim='valid_time')
    instant_data = xr.open_mfdataset(instant_files, combine='nested', concat_dim='valid_time')
    
    merged_data = xr.merge([accum_data, instant_data], compat='override')
    df = merged_data.to_dataframe().reset_index()
    df['valid_time'] = pd.to_datetime(df['valid_time'])
    
    # Convert UTC to Local Time
    def convert_utc_to_local(df):
        tf = TimezoneFinder()

        def get_fixed_offset(row):
            utc_time, latitude, longitude = row.valid_time, row.latitude, row.longitude
            if not isinstance(utc_time, pd.Timestamp):
                utc_time = pd.to_datetime(utc_time)
            if utc_time.tzinfo is None:
                utc_time = utc_time.tz_localize('UTC')
            timezone_str = tf.timezone_at(lng=longitude, lat=latitude)
            if timezone_str is None:
                return None
            local_zone = pytz.timezone(timezone_str)
            local_dt = utc_time.astimezone(local_zone)
            standard_time = local_dt - local_dt.dst()
            return standard_time.utcoffset().total_seconds() / 3600

        fixed_offset_hours = df.apply(get_fixed_offset, axis=1).iloc[0]
        print(f"Fixed offset hours: {fixed_offset_hours}")  # Print the fixed offset hours
        def convert_time(row):
            utc_time = pd.to_datetime(row.valid_time)
            if utc_time.tzinfo is None:
                utc_time = utc_time.tz_localize('UTC')
            new_local_time = utc_time + pd.Timedelta(hours=fixed_offset_hours)
            return new_local_time.replace(tzinfo=None)

        return df.apply(convert_time, axis=1)
    
    df["TIMESTAMP"] = convert_utc_to_local(df)
    df.drop(columns=['valid_time'], inplace=True)
    df.set_index('TIMESTAMP', inplace=True)
    df = df.resample('30T').interpolate(method='linear')
    df.reset_index(inplace=True)
    
    # Merge AmeriFlux and ERA5 data
    era5 = df.copy()
    era5["TIMESTAMP"] = pd.to_datetime(era5["TIMESTAMP"])
    merged_df = pd.merge(ameri, era5, on='TIMESTAMP', how='left')
    
    return merged_df

###########################################################################################################################

ameri_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps'
era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_alaska\A03'
ameri_file='gaps_US-A03.csv'
merged_data = merge_ameri_era5(ameri_path,era5_path,ameri_file)