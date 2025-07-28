# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 12:15:51 2025

@author: ammar
"""



import os 

# function to properly save blended data

def blended_save(df_b, ameri_file, save_path):
    """
    Save the blended DataFrame as a CSV file with a dynamic name based on the ameri_file name.

    Parameters:
    - df_b: DataFrame to be saved (blended data).
    - ameri_file: The name of the original Ameri data file used for the merge.
    - save_path: The path where the blended CSV file will be saved.

    Returns:
    - None
    """
    # Extract the identifier from ameri_file name (e.g., 'US-Skr' from 'gaps_US-Skr.csv')
    file_identifier = ameri_file.split('_')[1].replace('.csv', '')

    # Create the new file name based on the identifier
    output_file_name = f'gaps_blend_{file_identifier}.csv'

    # Full path to save the file
    output_file_path = os.path.join(save_path, output_file_name)

    # Save the blended DataFrame as CSV
    df_b.to_csv(output_file_path, index=False)

    print(f"Blended data saved to: {output_file_path}")
    



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from sklearn.metrics import r2_score


### function to blend and plot ameri era5 data

def blending_ameri_era(df, ameri_file):
    df_b = df.copy()  # Ensure the original data is not modified
    site_name = ameri_file.replace("gaps_", "").replace(".csv", "")

    # Convert DateTime column
    df_b['DateTime'] = pd.to_datetime(df_b['DateTime'])
    target_months = [5, 6, 7, 8]

    # Compute ERA-based variables
    df_b["Rg_era"] = (df_b.ssrd / 3600).round(3)
    df_b['Tair_era'] = (df_b.t2m - 273.15).round(3)
    df_b['dew_p'] = (df_b.d2m - 273.15).round(3)
    es = 0.6108 * np.exp((17.27 * df_b['Tair_era']) / (df_b['Tair_era'] + 237.3))
    ea = 0.6108 * np.exp((17.27 * df_b['dew_p']) / (df_b['dew_p'] + 237.3))
    df_b['VPD_era'] = ((es - ea) * 10).round(3)

    # Handle missing Tair or Rg if they are all NaN
    if df_b['Tair'].isna().all():
        print("All Tair values are NaN. Using Tair_era for missing Tair values.")
        df_b['Tair'] = df_b['Tair_era']
    
    if df_b['Rg'].isna().all():
        print("All Rg values are NaN. Using Rg_era for missing Rg values.")
        df_b['Rg'] = df_b['Rg_era']

    # Linear regression for Tair (only on rows with no NaNs)
    valid_mask_tair = df_b['Tair'].notna() & df_b['Tair_era'].notna()
    if valid_mask_tair.sum() > 1:
        slope_tair, intercept_tair, _, _, _ = linregress(df_b.loc[valid_mask_tair, 'Tair'], df_b.loc[valid_mask_tair, 'Tair_era'])
        df_b['Tair_era_c'] = ((df_b['Tair_era'] - intercept_tair) / slope_tair).round(3)

    # Linear regression for Rg
    valid_mask_rg = df_b['Rg'].notna() & df_b['Rg_era'].notna()
    if valid_mask_rg.sum() > 1:
        slope_rg, intercept_rg, _, _, _ = linregress(df_b.loc[valid_mask_rg, 'Rg'], df_b.loc[valid_mask_rg, 'Rg_era'])
        df_b['Rg_era_c'] = ((df_b['Rg_era'] - intercept_rg) / slope_rg).round(3)

    # Fill missing values based on the blending rules
    for year in df_b['DateTime'].dt.year.unique():
        mask_year = df_b['DateTime'].dt.year == year
        mask_target_months = mask_year & df_b['DateTime'].dt.month.isin(target_months)

        if mask_target_months.sum() == 0:
            print(f"No data for target months in year {year}. Skipping.")
            continue  

        # Calculate missing data % for target months for all variables
        missing_tair = df_b.loc[mask_target_months, 'Tair'].isna().sum() / mask_target_months.sum()
        missing_rg = df_b.loc[mask_target_months, 'Rg'].isna().sum() / mask_target_months.sum()
 
        # Fill missing values for Tair, Rg for all months in this year
        if missing_tair <= 0.5:
            df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair'] = df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair_era_c']
        else:
            df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair'] = df_b.loc[mask_year & df_b['Tair'].isna(), 'Tair_era']

        if missing_rg <= 0.5:
            df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg'] = df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg_era_c']
        else:
            df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg'] = df_b.loc[mask_year & df_b['Rg'].isna(), 'Rg_era']

    # Corrected indentation for the next loop
    for year in df_b['DateTime'].dt.year.unique():
        mask_year = df_b['DateTime'].dt.year == year
        mask_target_months = mask_year & df_b['DateTime'].dt.month.isin(target_months)

        # If no target months data available for this year, skip it (it means year ended before growing season)
        if mask_target_months.sum() == 0:
            print(f"No data for target months in year {year}. Skipping.")
            continue  

        missing_vpd = df_b.loc[mask_target_months, 'VPD'].isna().sum() / mask_target_months.sum()

        # --- 2. Use Tair and RH to calculate VPD if both are available ---
        mask_vpd_formula = mask_year & df_b['VPD'].isna() & df_b['Tair'].notna() & df_b['RH'].notna()
        df_b.loc[mask_vpd_formula, 'VPD'] = (
            0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_formula, 'Tair']) / 
                            (df_b.loc[mask_vpd_formula, 'Tair'] + 237.3)) - 
            (df_b.loc[mask_vpd_formula, 'RH'] / 100 * 
             0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_formula, 'Tair']) / 
                              (df_b.loc[mask_vpd_formula, 'Tair'] + 237.3)))
        ) * 10  # convert from kPa to hPa

        # --- 3. Handle completely missing VPD ---
        if df_b.loc[mask_year, 'VPD'].isna().all():
            if df_b.loc[mask_year, 'RH'].isna().all():
                print(f"[{year}] All VPD and RH values are NaN. Using VPD_era to fill missing VPD.")
                df_b.loc[mask_year, 'VPD'] = df_b.loc[mask_year, 'VPD_era']
            else:
                print(f"[{year}] All VPD values are NaN, but RH is available. Recalculating VPD using Tair and RH.")
                mask_vpd_fill = mask_year & df_b['VPD'].isna() & df_b['Tair'].notna() & df_b['RH'].notna()
                df_b.loc[mask_vpd_fill, 'VPD'] = (
                    0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_fill, 'Tair']) / 
                                    (df_b.loc[mask_vpd_fill, 'Tair'] + 237.3)) - 
                    (df_b.loc[mask_vpd_fill, 'RH'] / 100 * 
                     0.6108 * np.exp((17.27 * df_b.loc[mask_vpd_fill, 'Tair']) / 
                                      (df_b.loc[mask_vpd_fill, 'Tair'] + 237.3)))
                ) * 10

        # --- 4. Linear regression to correct VPD_era → VPD_era_c ---
        valid_mask_vpd = mask_year & df_b['VPD'].notna() & df_b['VPD_era'].notna()
        if valid_mask_vpd.sum() > 1:
            slope_vpd, intercept_vpd, *_ = linregress(
                df_b.loc[valid_mask_vpd, 'VPD'],
                df_b.loc[valid_mask_vpd, 'VPD_era']
            )
            df_b.loc[mask_year, 'VPD_era_c'] = ((df_b.loc[mask_year, 'VPD_era'] - intercept_vpd) / slope_vpd).round(3)
        else:
            df_b.loc[mask_year, 'VPD_era_c'] = df_b.loc[mask_year, 'VPD_era']

        # --- 5. Final fallback based on missingness threshold ---
        if missing_vpd <= 0.5:
            # Prioritize VPD_era_c if missing_vpd <= 0.5, else use VPD_era
            df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD'] = df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD_era_c']
        else:
            # If missing_vpd > 0.5, fallback to VPD_era
            df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD'] = df_b.loc[mask_year & df_b['VPD'].isna(), 'VPD_era']


    # 5. Cleanup: Remove unrealistic VPD and Rg values
    df_b.loc[(df_b['VPD'] < 0) | (df_b['VPD'] > 70), 'VPD'] = np.nan
    df_b.loc[df_b['Rg'] < 0, 'Rg'] = 0

    # Function to calculate R² and plot
    def plot_regression(var1, var2, label1, label2, title):
        mask = df[var1].notna() & df_b[var2].notna()
        if mask.sum() < 2:  # Ensure there are at least 2 points for regression
            print(f"Skipping plot for {title} due to insufficient data or missing variable.")
            return  # Skip this plot and continue with the next

        try:
            r2 = r2_score(df.loc[mask, var1], df_b.loc[mask, var2])
            plt.figure(figsize=(8, 6))
            sns.scatterplot(x=df.loc[mask, var1], y=df_b.loc[mask, var2], alpha=0.7, label="Data points")
            slope, intercept, _, _, _ = linregress(df.loc[mask, var1], df_b.loc[mask, var2])
            plt.plot(df.loc[mask, var1], slope * df.loc[mask, var1] + intercept, color='red', label=f'Linear Fit (R²={r2:.4f})')
            plt.xlabel(label1)
            plt.ylabel(label2)
            plt.title(f"{site_name} - {title}")
            plt.legend()
            plt.grid(True)
            plt.show()
        except Exception as e:
            print(f"Error occurred during regression for {title}: {e}. Skipping this plot.")

    # Plot R² comparisons
    plot_regression("Tair", "Tair_era", "Tair", "Tair_era", "Tair vs. Tair_era")
    plot_regression("Tair", "Tair_era_c", "Tair", "Tair_era_c", "Tair vs. Tair_era_c")
    plot_regression("Rg", "Rg_era", "Rg", "Rg_era", "Rg vs. Rg_era")
    plot_regression("Rg", "Rg_era_c", "Rg", "Rg_era_c", "Rg vs. Rg_era_c")
    plot_regression("VPD", "VPD_era", "VPD", "VPD_era", "VPD vs. VPD_era")
    plot_regression("VPD", "VPD_era_c", "VPD", "VPD_era_c", "VPD vs. VPD_era_c")

    return df_b  # Return the final blended dataset

save_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps'
ameri_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps'

# change these two links 
#########################################################################################################
# Florida: 5 sites

# make sure to mutiply VPD to 10 for proper unit at df step
era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_florida\evm'
ameri_file='gaps_US-EvM.csv'  

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_florida\ks3'
ameri_file='gaps_US-KS3.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_florida\ks4'
ameri_file='gaps_US-KS4.csv'


era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_florida\skr'
ameri_file='gaps_US-Skr.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_florida\tas'
ameri_file='gaps_US-TaS.csv'


############################################################################################################
##california: 7 sites

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\dmg'
ameri_file='gaps_US-Dmg.csv'


era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\ekh'
ameri_file='gaps_US-EKH.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\edn'
ameri_file='gaps_US-EDN.csv'


era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\ekp'
ameri_file='gaps_US-EKP.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\myb'
ameri_file='gaps_US-Myb.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\srr'
ameri_file='gaps_US-Srr.csv'


era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_california\tw1'
ameri_file='gaps_US-Tw1.csv'


##############################################################################################################

# alaska: 4 sites

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_alaska\A03'
ameri_file='gaps_US-A03.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_alaska\A10'
ameri_file='gaps_US-A10.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_alaska\ngb'
ameri_file='gaps_US-NGB.csv'

# make sure to mutiply VPD to 10 for proper unit at df step
era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_alaska\atq'
ameri_file='gaps_US-Atq.csv'  

####################################################################################################################

## south carolina: 4 sites

# make sure to mutiply VPD to 10 for proper unit at df step
era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_s_carolina\sts'
ameri_file='gaps_US-StS.csv'  

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_s_carolina\HB1'
ameri_file='gaps_US-HB1.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_s_carolina\HB2'
ameri_file='gaps_US-HB2.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_s_carolina\HB3'
ameri_file='gaps_US-HB3.csv'


######################################################################################################################
## delaware: 1 site 


era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_delaware\stj'
ameri_file='gaps_US-StJ.csv'

##################################################################################################################
# louisina: 3 sites
era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_louisiana\LA1'
ameri_file='gaps_US-LA1.csv'


era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_louisiana\LA2'
ameri_file='gaps_US-LA2.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_louisiana\LA3'
ameri_file='gaps_US-LA3.csv'

##################################################################################################################
# New jersey: 2 sites

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_New_Jersey\hpy'
ameri_file='gaps_US-HPY.csv'

era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_data_New_Jersey\mrm'
ameri_file='gaps_US-MRM.csv'

#################################################################################################################
# massachusetts: 1 site
era5_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\Era5_data\ERAS_massachusetts\phm'
ameri_file='gaps_US-PHM.csv'

# total sites 27
####################################################################################################################
df=merge_ameri_era5(ameri_path,era5_path,ameri_file)
# use VPD*10 step only on three sites  (US-EvM, US-StS, US-AtQ)
#df["VPD"]=df["VPD"]*10  # uncomment it soon so that you dont mutiply it twice
# Perform blending of data
df_b=blending_ameri_era(df,ameri_file)
# Call the blended_save function to save the data
blended_save(df_b, ameri_file, save_path)


#21840

# check 99 percentile of a variable 
df_b['VPD'].quantile(0.99)
# %age of missing data

df['VPD'].isna().mean() * 100