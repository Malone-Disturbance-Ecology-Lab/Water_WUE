# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 17:27:52 2025

@author: ammar
"""
## tested on 28 ameriflux sites
## see paths to find path to csv


import os
import datetime
import pandas as pd
import numpy as np
import zipfile
from pathlib import Path
import matplotlib.pyplot as plt


save_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps'



###########################################################################################################################
#Florida =5
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-EvM_BASE_HH_2-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-KS3_BASE_HH_1-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-KS4_BASE_HH_3-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Skr_BASE_HH_2-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-TaS_BASE_HH_1-5.csv'  # Replace with the actual path
########################################################################################################################
#california=7

csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Dmg_BASE_HH_4-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-EDN_BASE_HH_3-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-EKP_BASE_HH_1-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Myb_BASE_HH_14-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Srr_BASE_HH_1-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Tw1_BASE_HH_11-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-EKH_BASE_HH_1-5.csv'  # Replace with the actual path

#########################################################################################################################

# alaska=4

csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-A03_BASE_HH_5-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-A10_BASE_HH_4-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Atq_BASE_HH_1-1.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-NGB_BASE_HH_5-5.csv'  # Replace with the actual path

###########################################################################################################################
#delaware 1

csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-StJ_BASE_HH_2-5.csv'  # Replace with the actual path
#####################################################################################################################

# lousiana 3
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-LA1_BASE_HH_2-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-LA2_BASE_HH_4-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-LA3_BASE_HH_2-5.csv'  # Replace with the actual path

############################################################################################################################
# massachusetts  1

csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-PHM_BASE_HH_3-5.csv'  # Replace with the actual path
############################################################################################################################
# new jersey 2

csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Hpy_BASE-HH_2-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-Mrm_BASE_HH_2-5.csv'  # Replace with the actual path

#########################################################################################################################

# south carolina 4

csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-HB1_BASE_HH_3-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-HB2_BASE_HH_1-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-HB3_BASE_HH_2-5.csv'  # Replace with the actual path
csv_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps\AMF_US-StS_BASE_HH_2-4.csv'  # Replace with the actual path

########################################################################################################################


df=process_ameriflux_data(csv_file_path, save_folder)



def process_ameriflux_data(csv_file_path, save_folder):
    os.makedirs(save_folder, exist_ok=True)

    def check_header(file_path):
        try:
            first_row = pd.read_csv(file_path, nrows=5)
            first_line = first_row.columns.tolist()
            if all(isinstance(col, str) and col.strip() for col in first_line):
                return 0
            else:
                return 2
        except Exception:
            return 2

    skip_rows = check_header(csv_file_path)

    try:
        df = pd.read_csv(csv_file_path, skiprows=skip_rows)
        if df.shape[1] == 1:
            df = pd.read_csv(csv_file_path, skiprows=skip_rows, sep='\t')
        print(f"Data from {os.path.basename(csv_file_path)} (skip_rows={skip_rows}):\n", df.head())
    except Exception as e:
        print(f"Error reading {csv_file_path}: {e}")
        return

    # Parse timestamp
    if 'TIMESTAMP_START' in df.columns:
        df['TIMESTAMP'] = pd.to_datetime(df['TIMESTAMP_START'], format='%Y%m%d%H%M')
    elif 'TIMESTAMP' in df.columns:
        df['TIMESTAMP'] = pd.to_datetime(df['TIMESTAMP'])
    elif 'datetime' in df.columns:
        df['TIMESTAMP'] = pd.to_datetime(df['datetime'])
    else:
        raise ValueError("No suitable timestamp column found in the dataset.")

    df = df.drop_duplicates(subset='TIMESTAMP').copy()
    df.set_index('TIMESTAMP', inplace=True)

    # Generate complete 30-minute range covering full days
    start = df.index.min().floor('D')
    end = df.index.max().ceil('D') - pd.Timedelta(minutes=30)
    full_index = pd.date_range(start=start, end=end, freq='30min')

    # Reindex and fill missing values
    df = df.reindex(full_index)
    df.index.name = 'TIMESTAMP'
    df = df.replace([-9999, -6999], np.nan)

    # Create required time columns
    df = df.reset_index()
    df['DateTime'] = df['TIMESTAMP']
    df['Year'] = df['TIMESTAMP'].dt.year
    df['DoY'] = df['TIMESTAMP'].dt.dayofyear
    df['Hour'] = df['TIMESTAMP'].dt.hour + df['TIMESTAMP'].dt.minute / 60 + 0.5
    df['Hour'] = df['Hour'].apply(lambda x: round(x * 2) / 2)

    cols_to_check = ["NEE_PI_F", "NEE_PI", "NEE"]
    df["NEE"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan
    df["NEE"].fillna(df.get("FC", 0) + df.get("SC", 0), inplace=True)
    df["NEE"].fillna(df.get("FC", 0) + df.get("SFC", 0), inplace=True)
    df["NEE"].fillna(df.get("FC", df.get("FC_PI_F")), inplace=True)
    df["NEE"] = np.where(df["NEE"].between(-50, 50), df["NEE"], np.nan)

    cols_to_check = ["LE_PI_F", "LE_PI", "LE"]
    df["LE"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan
    df["LE"] = df["LE"].apply(lambda x: np.nan if x < -200 or x > 800 else x)

    cols_to_check = ["H_PI_F", "H"]
    df["H"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan
    df["H"] = df["H"].apply(lambda x: np.nan if x < -200 or x > 800 else x)

    cols_to_check = ["SW_IN_PI_F", "SW_IN_F", "SW_IN", "SW_IN_1_1_1", 'Rg']
    df["Rg"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan
    df["Rg"] = df["Rg"].apply(lambda x: 0 if x < 0 else x)

    cols_to_check = ["TA_PI_F", "TA", "TA_1_1_1"]
    df["Tair"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan

    cols_to_check = ["VPD_PI_F", "VPD_PI", "VPD"]
    df["VPD"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan

    if "RH" not in df.columns:
        df["RH"] = np.nan
    
    if df["VPD"].isna().all():  
        cols_to_check = ["RH", "RH_1_1_1"]
        if any(col in df for col in cols_to_check):
            df["RH"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0]
        else:
            df["RH"] = np.nan
        
        if "Tair" in df.columns and "RH" in df.columns and not df["RH"].isna().all():
            es = 0.6108 * np.exp((17.27 * df["Tair"]) / (df["Tair"] + 237.3))
            ea = df["RH"] / 100 * es
            df["VPD"] = (es - ea) * 10

    df["VPD"] = df["VPD"].apply(lambda x: np.nan if x < 0 else x)

    cols_to_check = ["USTAR", "UST"]
    df["USTAR"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan
    df["Ustar"] = df["USTAR"]

    cols_to_check = ["NETRAD", "NETRAD_1_1_1", "Rn"]
    df["NETRAD"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan

    cols_to_check = ["PA", "PA_1_1_1"]
    df["PA"] = df[[col for col in cols_to_check if col in df]].bfill(axis=1).iloc[:, 0] if any(col in df for col in cols_to_check) else np.nan

    reframed_columns = ["DateTime", "Year", "DoY", "Hour", "NEE", "LE", "H", "Rg", "Tair", "VPD", "Ustar", "PA", "NETRAD", "RH"]
    reframed = pd.concat([df[col] for col in reframed_columns if col in df.columns], axis=1)

    site_name = os.path.basename(csv_file_path).split('_')[1]  
    new_filename = f"gaps_{site_name}.csv"
    csv_save_path = os.path.join(save_folder, new_filename)
    reframed.to_csv(csv_save_path, index=False)
    print(f"DataFrame saved as '{csv_save_path}'")

    reframed.set_index("DateTime", inplace=True)
    
    columns_to_plot = [col for col in reframed.columns if col not in ["Year", "DoY", "Hour"]]

    for col in columns_to_plot:
        plt.figure(figsize=(12, 5))
        plt.plot(reframed.index, reframed[col], color='b', alpha=0.7)
        plt.xlabel('Date')
        plt.ylabel(col)
        plt.title(f'{site_name} {col}')
        plt.grid(True)
        plt.show()

    return df



        





