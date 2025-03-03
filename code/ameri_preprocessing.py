# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:22:26 2024

@author: ammar
"""

# es jobs google sheet
#eco log listserv
#indeed
# chronical of higher education 

#NEE µmolCO2 m-2 s-1

#NEE= FC+SC
from sklearn.metrics import r2_score
def stat (obs,pre):
    nse=NSE(pre,obs)
   # pear=pearsonr(obs,pre)
    wil=will(pre,obs)   
    MAE=mean_absolute_error(obs,pre)
    RMSE=sqrt(mean_squared_error(obs,pre))
   # reg_=reg.score(x1, y1)
    pbias=(np.sum(pre-obs)/np.sum(obs))*100
    return nse,pear,MAE,reg_,pbias


#from sklearn.linear_model import LinearRegression
#from sklearn.metrics import mean_absolute_error
#import matplotlib as mpl

#from matplotlib import scale as mscale
#from matplotlib import transforms as mtransforms
#from math import sqrt
#import matplotlib.dates as mdates
#from sklearn.metrics import r2_score
#import scipy
#from numpy.polynomial.polynomial import polyfit
#import matplotlib.lines as mlines
#import seaborn as sns
#from matplotlib.legend import Legend
#import matplotlib.dates as mdates
#from matplotlib.ticker import FormatStrFormatter
#import pymannkendall as mk
#import math
#import statsmodels.formula.api as smf
#from scipy.stats.stats import pearsonr
##################################################################################
import os
import datetime
import pandas as pd
import numpy as np
import zipfile
from pathlib import Path
import matplotlib.pyplot as plt
# this code will unzip the file
###########################################################################################################
# take desired data from zipped folder and save it in a different folder

# change the path to zipped file 
zip_file_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\AMF_US-EvM_BASE-BADM_2-5.zip'  # Replace with the actual path

# Directory to extract the contents
extracted_folder = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\ameri_gaps'  # Replace with your desired directory to extract

# Create the directory if it doesn't exist
os.makedirs(extracted_folder, exist_ok=True)

# Open the zip file and extract only the files with 'BASE_HH' in their name
with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
    # Iterate through each file in the zip archive
    for file_name in zip_ref.namelist():
        # Check if 'BASE_HH' is in the file name
        if 'HH' in file_name or 'BASE_HH' in file_name:
            # Extract the file to the specified directory
            zip_ref.extract(file_name, extracted_folder)
            
            # Construct the full path of the extracted file
            extracted_file_path = os.path.join(extracted_folder, file_name)
            
            # Read the file using pandas, skipping the first two rows
            try:
                df = pd.read_csv(extracted_file_path, skiprows=2)
                print(f"Data from {file_name}:\n", df.head())  # Print first few rows for confirmation
            except Exception as e:
                print(f"Error reading {file_name}: {e}")


# Create the TIMESTAMP column using timestamp_start, in the desired formatn
df['TIMESTAMP_START'] = pd.to_datetime(df['TIMESTAMP_START'], format='%Y%m%d%H%M')
df['TIMESTAMP_END'] = pd.to_datetime(df['TIMESTAMP_END'], format='%Y%m%d%H%M')
df['TIMESTAMP'] = df['TIMESTAMP_START'].dt.strftime('%Y-%m-%d %H:%M:%S')
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
df["datetime"]=df['TIMESTAMP']

# find the number of dates in the timeseries, that will help to produce index needed for creating half hour column 
start_date = str(df['datetime'].iloc[0])
end_date = str(df['datetime'].iloc[-1])
# confirm that not a single half hour timestamp is missing 
idx = pd.date_range(start=start_date, end=end_date, freq="30T")
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
# create other columns needed for reddy proc gap filling 
df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
k=pd.concat([s] * int(idx.size / (24 * 2)), axis=0)
k.reset_index(drop=True, inplace=True)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114] # remove extra columns 
df.head(10)
df.tail(10)
back=df
#######################################################################################################

df=back

###############################################################################################################
df=back

# Assign 'NEE' column based on availability of 'NEE_PI' or 'NEE_PI_F'
df["NEE_"] = df.get("NEE_PI_F", df.get("NEE_PI", df.get("NEE")))

# Fill gaps in 'NEE' using fallback values
df["NEE_"].fillna(df.get("FC", 0) + df.get("SC", 0), inplace=True)  # Use 'FC' + 'SC' if available
df["NEE_"].fillna(df.get("FC", df.get("FC_PI_F")), inplace=True)  # Use 'FC' if available, else 'FC_PI_F'
df["NEE_"].fillna(df.get("NEE_PI_F", df.get("NEE_PI", df.get("NEE"))), inplace=True)  # As a last resort, use 'NEE_PI_F' or 'NEE_PI'or 'NEE' for gap filling again
# FINAL STEP: If 'NEE' is still missing, use 'FC' or 'SC' as a last resort
df["NEE_"].fillna(df.get("FC", df.get("SC")), inplace=True)  
df["NEE"]=df["NEE_"]

plt.plot(df.NEE)

df["LE_"] = df.get("LE_PI_F", df.get("LE_PI", df.get("LE")))
# Fill gaps in 'LE' using available columns
df["LE_"].fillna(df.get("LE_PI_F", df.get("LE_PI", df.get("LE"))), inplace=True)
df["LE"]=df["LE_"]

plt.plot(df.LE)

df["Rg_"] = df.get("SW_IN_PI_F", df.get("SW_IN", df.get("SW_IN_1_1_1")))
df["Rg_"].fillna(df.get("SW_IN_PI_F", df.get("SW_IN", df.get("SW_IN_1_1_1"))), inplace=True) # fill in gaps 
df['Rg']=df['Rg_']
plt.plot(df.Rg)



#try different names 
df["NEE"]=df["NEE_PI"]

df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]
df["NEE_check"]=df["FC"]+df["SC"]
np.nanpercentile(df["NEE_check"],[1, 99])
df["NEE1"] = np.where(df["NEE_check"]<-50, np.NaN, df["NEE_check"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE_check"]=df["NEE2"]

np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]
df["Rg"]=df["SW_IN"]
df.loc[df['Rg'] < 0, 'Rg'] = 0
df["Tair"]=df["TA_1_1_1"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["RH"]=df["RH_1_1_1"]
df['NEE'] = df['NEE'].fillna(df['NEE_check'])
ameri=df

#####################################################################################################

merged_df = pd.merge(ameri, era, on='TIMESTAMP', how='left')

df=merged_df 
back3=df


df['Rg'] = df['Rg'].fillna(df['ssrd_'])
df['Tair'] = df['Tair'].fillna(df['Tair_er'])
df['VPD'] = df['VPD'].fillna(df['VPD_er'])

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"],df["LE"],
   df["H"],df['Rg'],df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#only fill columns you need
### save data in folder
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
reframed.to_csv('gaps_US-Skr.csv', index=False, header=True)

















# Step 1: Unzip the file
with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
    zip_ref.extractall(extracted_dir)  # Extract all files to the specified directory

# Step 2: List files in the extracted directory
extracted_files = os.listdir(extracted_dir)
print(f"Extracted Files: {extracted_files}")
################################################################################################################



base_dir = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data")

df = pd.read_csv(base_dir / "florida" / "AMF_US-Skr_BASE_HH_2-5.csv")

df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
df["datetime"]=df['TIMESTAMP']
start_date = str(df['datetime'].iloc[0])
end_date = str(df['datetime'].iloc[-1])
# Generate half-hour intervals for all years
idx = pd.date_range(start=start_date, end=end_date, freq="30T")
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################
#change this based on data size (number of days of data), 

df=back2
k=pd.concat([s] * int(idx.size / (24 * 2)), axis=0)
k.reset_index(drop=True, inplace=True)

df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.head(10)
df.tail(10)
#try different names 
df["NEE"]=df["NEE_PI"]

df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]
df["NEE_check"]=df["FC"]+df["SC"]
np.nanpercentile(df["NEE_check"],[1, 99])
df["NEE1"] = np.where(df["NEE_check"]<-50, np.NaN, df["NEE_check"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE_check"]=df["NEE2"]

np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]
df["Rg"]=df["SW_IN"]
df.loc[df['Rg'] < 0, 'Rg'] = 0
df["Tair"]=df["TA_1_1_1"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["RH"]=df["RH_1_1_1"]
df['NEE'] = df['NEE'].fillna(df['NEE_check'])
ameri=df

#####################################################################################################

merged_df = pd.merge(ameri, era, on='TIMESTAMP', how='left')

df=merged_df 
back3=df


df['Rg'] = df['Rg'].fillna(df['ssrd_'])
df['Tair'] = df['Tair'].fillna(df['Tair_er'])
df['VPD'] = df['VPD'].fillna(df['VPD_er'])

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"],df["LE"],
   df["H"],df['Rg'],df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#only fill columns you need
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
reframed.to_csv('gaps_US-Skr.csv', index=False, header=True)

######################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
df=pd.read_csv("AMF_US-TaS_BASE_HH_1-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df

#2016-2023  # all years

df=back
df["datetime"]=df['TIMESTAMP']
df.head(10)
df.tail(10)

start_date = "2016-01-01 00:00:00"
end_date = "2023-12-31 23:59:59"

# Generate half-hour intervals for all years
idx = pd.date_range(start=start_date, end=end_date, freq="30T")
#idx=pd.date_range(start='01/01/2016 00:00:00', periods=140256, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df


df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################
#change this based on data size (number of days of data), 
#([s]*(#number of days)
df=back2
k=pd.concat([s] * int(idx.size / (24 * 2)), axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.head(10)
df.tail(10)
#try different names 
#df["NEE"]=df["NEE_PI"]
df["NEE"]=df["FC"]+df["SC"]
np.nanpercentile(df["NEE"],[1, 99]) # -21.6, 10

df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]


df["Rg"]=df["SW_IN"]
df["Tair"]=df["TA_1_1_1"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["RH"]=df["RH_1_1_1"]
df["NETRAD"]=(df["SW_IN"]+df["LW_IN"])-(df["SW_OUT"] + df["LW_OUT"])
ameri=df



years_to_remove = [2016]
# Remove rows based on years
df = df[~df['DateTime'].dt.year.isin(years_to_remove)]
df.index=np.arange(0, len(df))

########################################################################################################
reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"],df["LE"],df["H"],
                    df['Rg'],df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
#only fill columns you need
reframed.to_csv('gaps_US-TaS.csv', index=False, header=True)
##########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
df=pd.read_csv("AMF_US-EvM_BASE_HH_2-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2020 00:00:00', periods=70128, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df


df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################
#change this based on data size (number of days of data), 
#([s]*(#number of days)
df=back2
k=pd.concat([s]*(1461),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.head(10)
df.tail(10)
#try different names 
df["NEE"]=df["NEE_PI_F"]
np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]


df["NEE_check"]=df["FC"]+df["SC"]
np.nanpercentile(df["NEE_check"],[1, 99])
df["NEE1"] = np.where(df["NEE_check"]<-50, np.NaN, df["NEE_check"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE_check"]=df["NEE2"]

df["LE"]=df["LE_PI_F"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]


df["Rg"]=df["SW_IN"]
df["Tair"]=df["TA_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]

df['NEE'] = df['NEE'].fillna(df['NEE_check'])
########################################################################################################
reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],
                    df['Rg'],df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
#only fill columns you need
reframed.to_csv('gaps_US-EvM.csv', index=False, header=True)

#######################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
df=pd.read_csv("AMF_US-KS3_BASE_HH_1-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df

df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2018 00:00:00', periods=21840, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2

k=pd.concat([s]*(455),axis=0)
k.reset_index(drop=True, inplace=True)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

df["NEE"]=df["FC"]+df["SC"]
np.nanpercentile(df["NEE"],[1, 99])

df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]
df["Rg"]=df["SW_IN"]
df["Tair"]=df["TA"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["LE"]=df["LE"]


########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
reframed.to_csv('gaps_US-KS3.csv', index=False, header=True)

#####################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
df=pd.read_csv("AMF_US-KS4_BASE_HH_3-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2017 00:00:00', periods=92208, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(1921),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

df["NEE"]=df["FC"]+df["SC"]
np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]



df["Rg"]=df["SW_IN"]

#df["Rg"]=df["PPFD_IN_PI_F_1_1_1"]
df['Rg'] = df['Rg'].apply(lambda x: max(x, 0))
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]

########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\Florida")
#only fill columns you need
reframed.to_csv('gaps_US-KS4.csv', index=False, header=True)

#######################################################################################################

# california sites

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-Dmg_BASE_HH_4-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2021 00:00:00', periods=47808, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(996),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question

np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

df["LE"]=df["LE_PI_F"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]


df["H"]=df["H_PI_F"]
df["Rg"]=df["SW_IN_PI_F"]
df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["VPD"]=df["VPD_PI_F"]
df["Ustar"]=df["USTAR"]

df["reco"]=df["RECO_PI_F"]
df["GPP"]=df["GPP_PI_F"]


########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-Dmg.csv', index=False, header=True)

# fill all and see how your GPP compare 
##########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-EDN_BASE_HH_3-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2018 00:00:00', periods=60624, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(1263),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question

np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]


#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
df["Rg"]=df["SW_IN"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]

#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]


########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-EDN.csv', index=False, header=True)

########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-EDN_BASE_HH_3-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2018 00:00:00', periods=60624, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(1263),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question

np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]


#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
df["Rg"]=df["SW_IN"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]

#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]


########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-EDN.csv', index=False, header=True)

########################################################################################################
#######################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-EKP_BASE_HH_1-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2021 00:00:00', periods=63696, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(1327),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question

np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
VPD=(es-ea)*10
df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
#df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]

#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]
########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-EKP.csv', index=False, header=True)
########################################################################################################

# check why Rg did not turn out to Rg_f

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-Myb_BASE_HH_14-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2010 00:00:00', periods=245424, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(5113),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

df["LE"]=df["LE_PI_F"]
#df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

df["H"]=df["H_PI_F"]
#df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["reco"]=df["RECO_PI_F"]
df["GPP"]=df["GPP_PI_F"]

years_to_remove = [2010]
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]


###################################################################################################

## ustar issue so can not fill

df["LE"]=df.LE.interpolate() 
df.index = np.arange(0, len(df))

########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-Myb.csv', index=False, header=True)
########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-Srr_BASE_HH_1-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2014 00:00:00', periods=65808, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(1371),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

df["LE"]=df["LE_PI_F"]
#df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

df["H"]=df["H_PI_F"]
#df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["reco"]=df["RECO_PI_F"]
df["GPP"]=df["GPP_PI_F"]
########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-Srr.csv', index=False, header=True)

########################################################################################################



os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-Tw1_BASE_HH_11-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2011 00:00:00', periods=227904, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(4748),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

df["LE"]=df["LE_PI_F"]
#df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

df["H"]=df["H_PI_F"]
#df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["reco"]=df["RECO_PI_F"]
df["GPP"]=df["GPP_PI_F"]
########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-Tw1.csv', index=False, header=True)

########################################################################################################
## massachusetts sites
#######################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\massachusetts")
df=pd.read_csv("AMF_US-PHM_BASE_HH_3-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2013 00:00:00', periods=140256, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(2922),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question #NEE_PI_F is missing in this site even if name of the column is 
#written
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Rg"]=df["SW_IN"]

df["Rg"]=df["SW_IN_1_1_1"]

#df["Rg"]=df["PPFD_IN_PI_F_1_1_1"]
df['Rg'] = df['Rg'].apply(lambda x: max(x, 0))

df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
df["GPP"]=df["GPP_PI_F"]

########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-PHM.csv', index=False, header=True)

#########################################################################################################

##### south carolina

#####################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\south carolina")
df=pd.read_csv("AMF_US-HB1_BASE_HH_3-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2019 00:00:00', periods=70128, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(1461),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

# NETRAD is missing on this site 
########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-HB1.csv', index=False, header=True)

########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\south carolina")
df=pd.read_csv("AMF_US-HB2_BASE_HH_1-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2019 00:00:00', periods=17520, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(365),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Rg"]=df["SW_IN"]
df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

# NETRAD is missing on this site 
########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-HB2.csv', index=False, header=True)

##########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\south carolina")
df=pd.read_csv("AMF_US-HB3_BASE_HH_2-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2019 00:00:00', periods=70128, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(1461),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

df["PA"]=df["PA_1_1_1"]

# NETRAD is missing on this site 
########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-HB3.csv', index=False, header=True)
##########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\new jersey")
df=pd.read_csv("US-Hpy_HH_201501010000_201801010000.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2015 00:00:00', periods=52608, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(1096),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df = df.drop('Hour', axis=1) # drop hour column is already exist 
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
#df["VPD"]=df["VPD_PI"]
df["VPD"]=df["VPD"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

#df["PA"]=df["PA_1_1_1"]

# NETRAD is missing on this site 
########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-Hpy.csv', index=False, header=True)
##########################################################################################################
from datetime import datetime, timedelta

def count_half_hours(start_date, end_date ):
    # Convert string input to datetime objects if necessary
    if isinstance(start_date, str):
        start_date = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
    if isinstance(end_date, str):
        end_date = datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")
    
    # Calculate total duration in minutes
    total_minutes = (end_date - start_date).total_seconds() / 60
    
    # Calculate half-hour intervals and include the start & end
    half_hours = (total_minutes / 30) + 1
    
    return int(half_hours)
######################################################################################################

### alaska 

#################################################################################################

start_date = "2014-01-01 00:00:00"
end_date = "2021-12-31 23:30:00"

count_half_hours(start_date, end_date )

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\alaska")
df=pd.read_csv("AMF_US-A03_BASE_HH_5-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2014 00:00:00', periods=140256, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(2922),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
#df = df.drop('Hour', axis=1) # drop hour column is already exist 
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question
df["NEE"]=df["NEE_PI"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Rg"]=df["SW_IN"]
df["Rg"]=df["SW_IN"]
plt.plot(df.Rg)



#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]

#df["Tair"]=df["TA_1_1_1"]
ameri=df

#merged_df = pd.merge(ameri, Rg_LA1, on='TIMESTAMP', how='inner')

merged_df = pd.merge(ameri, ta_A03, left_on="TIMESTAMP", right_on="TIMESTAMP2")

# Drop one of the date columns
merged_df = merged_df.drop(columns=["TIMESTAMP"])

df=merged_df
## rh = (ea/es)*100
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
#df["VPD"]=df["VPD_PI"]
df["VPD"]=df["VPD"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

#df["PA"]=df["PA_1_1_1"]

# NETRAD is missing on this site 



####################################################################################################
df=merged_df
########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-Mrm.csv', index=False, header=True)
#############################################################################################################


#55 link forum
#https://forum.ecmwf.int/latest






















































df=pd.read_csv("AMF_US-Atq_BASE_HH_1-1.csv") 


df=pd.read_csv("AMF_US-A10_BASE_HH_4-5.csv") 











df=pd.read_csv("AMF_US-A03_BASE_HH_5-5.csv") 






















########################################################################################################

#delaware

############################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\delaware")

df=pd.read_csv("AMF_US-StJ_BASE_HH_2-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2014 00:00:00', periods=70128, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(1461),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
#df["NEE"]=df["FC"] # no SC term as question
df["NEE"]=df["FC_PI_F"] # no SC term as question
df["NEE"]=df["FC_PI"] 

np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

df["PA"]=df["PA_1_1_1"]

########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-HB3.csv', index=False, header=True)
######################################################################################################

#Louisiana
###################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")

df=pd.read_csv("AMF_US-LA1_BASE_HH_2-5.csv") 

df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)
#####################################################################################################3
# remove 2011 gaps
# don't need that in automation code 
years_to_remove = [2011]
# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
df.head(10)
df.tail(10)
######################################################################################################3
df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2012 00:00:00', periods=17568, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(366),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Rg"]=df["SW_IN"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

#df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

#df["PA"]=df["PA_1_1_1"]

# netrad and p is missing for LA1
ameri=df

#merged_df = pd.merge(ameri, Rg_LA1, on='TIMESTAMP', how='inner')

merged_df = pd.merge(ameri, Rg_LA1, left_on="TIMESTAMP", right_on="TIMESTAMP2")

# Drop one of the date columns
merged_df = merged_df.drop(columns=["TIMESTAMP"])

####################################################################################################
df=merged_df

########################################################################################################

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["VPD"],df["Ustar"]),axis=1)

#only fill columns you need
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")
reframed.to_csv('gaps_US-LA1.csv', index=False, header=True)

###########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")

df=pd.read_csv("AMF_US-LA2_BASE_HH_3-5.csv") 

df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)
#####################################################################################################3
# remove 2011 gaps
# don't need that in automation code 
years_to_remove = [2011]
# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
df.head(10)
df.tail(10)
######################################################################################################3

df=back
df["datetime"]=df['TIMESTAMP']


idx=pd.date_range(start='01/01/2012 00:00:00', periods=192864, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(4018),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
df["Rg"]=df["SW_IN"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
df["Tair"]=df["TA"]
#df["Tair"]=df["TA_1_1_1"]
#es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

#ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
#VPD=(es-ea)*10
#df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

df["NETRAD"]=df["NETRAD"]
#df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

#df["PA"]=df["PA_1_1_1"]
df["PA"]=df["PA"]


ameri=df


merged_df = pd.merge(ameri, Rg_LA2, left_on="TIMESTAMP", right_on="TIMESTAMP2")

df=merged_df

df['Rg']=df['Rg'].fillna(df['Rg_'])

####################################################################################################


reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"]),axis=1)

#only fill columns you need
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")
reframed.to_csv('gaps_US-LA2.csv', index=False, header=True)
# 2012, 2013 
###########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")

df=pd.read_csv("AMF_US-LA3_BASE_HH_1-5.csv") 

df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)
#####################################################################################################3
#
df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2019 00:00:00', periods=70128, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

df=back2
k=pd.concat([s]*(1461),axis=0) # number of days 
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question
#df["NEE"]=df["FC_PI_F"] # no SC term as question


np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]

#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
#df["Rg"]=df["SW_IN_1_1_1"]
#df["Rg"]=df["SW_IN"]
df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
VPD=(es-ea)*10
df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
#df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
#df["NETRAD"]=df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]
#df["NETRAD"]=df["SW_IN"]-df["SW_OUT"]+df["LW_IN"]-df["LW_OUT"]

#df["NETRAD"]=df["NETRAD"]
df["NETRAD"]=df["NETRAD_1_1_1"]
#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]

#df["PA"]=df["PA_1_1_1"]
df["PA"]=df["PA"]

ameri=df

merged_df = pd.merge(ameri, Rg_LA3, left_on="TIMESTAMP", right_on="TIMESTAMP2")

df=merged_df

df['Rg']=df['Rg'].fillna(df['Rg_'])

####################################################################################################


reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
               df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 # df["Tair"],df["VPD"],df["Ustar"],df["PA"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["GPP"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"]),axis=1)

#only fill columns you need
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")
reframed.to_csv('gaps_US-LA3.csv', index=False, header=True)

##########################################################################################################

# site removed



os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("AMF_US-EKH_BASE_HH_1-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
df.head(10)
df.tail(10)

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2022 00:00:00', periods=46176, freq='30T')  # half hourly. 
df.index = pd.DatetimeIndex(df.datetime)
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
#df= df.replace(-9999, np.nan) 
df.isnull().values.any()
df=df.drop('TIMESTAMP', axis=1)  
df=df.rename_axis('TIMESTAMP').reset_index()  
back1=df

df["DateTime"]=df['TIMESTAMP']
df['DoY'] = pd.to_datetime (df['DateTime']).dt.dayofyear 
df['Year'] = df['TIMESTAMP'].dt.year
df['day'] = df['TIMESTAMP'].dt.dayofyear

s=pd.Series(np.arange(0.5,24.5,0.5))
back2=df

#######################################################################################################

df=back2
k=pd.concat([s]*(962),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.tail(10)
back3=df

#df["NEE"]=df["FC"]+df["SC"]
df["NEE"]=df["FC"] # no SC term as question

np.nanpercentile(df["NEE"],[1, 99])
df["NEE1"] = np.where(df["NEE"]<-50, np.NaN, df["NEE"])
df["NEE2"] = np.where(df["NEE1"]>50, np.NaN, df["NEE1"])
df["NEE"]=df["NEE2"]

#df["LE"]=df["LE_PI_F"]
df["LE"]=df["LE"]
np.nanpercentile(df["LE"],[1, 99])
df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
df["LE"]=df["LE2"]


#df["H"]=df["H_PI_F"]
df["H"]=df["H"]
#df["Rg"]=df["SW_IN_PI_F"]
df["Rg"]=df["SW_IN_1_1_1"]
#df["Tair"]=df["TA_PI_F"]
#df["Tair"]=df["TA"]
df["Tair"]=df["TA_1_1_1"]
es=0.6108*np.exp((17.27*df["Tair"])/(df["Tair"]+237.3))   # in kpa
#ea is actual vapoure pressure in millibar and es is saturated vapor pressure in kpa

## rh = (ea/es)*100

ea=df["RH_1_1_1"]/100*es    ### same as multiplying with es
VPD=(es-ea)*10
df['VPD']=VPD   # in kpa

#df["VPD"]=df["VPD_PI_F"]
#df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["NETRAD"]=(df["SW_IN_1_1_1"]-df["SW_OUT_1_1_1"]+df["LW_IN_1_1_1"]-df["LW_OUT_1_1_1"]).astype(float)

#df["reco"]=df["RECO_PI_F"]
#df["GPP"]=df["GPP_PI_F"]
########################################################################################################

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
                df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"], df["LE"],df["H"],df['Rg'],
 #                   df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"],df["reco"],df["GPP"]),axis=1)

#only fill columns you need
reframed.to_csv('gaps_US-EKH.csv', index=False, header=True)




















######################################################################################################
