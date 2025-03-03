# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:51:24 2025

@author: ammar
"""


# ERA data is in utc 
import xarray as xr
import pandas as pd
from zoneinfo import ZoneInfo
import pandas as pd
from netCDF4 import Dataset




#os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida\ERAS_data\srk")

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida\ERAS_data\tas")
ds1 = xr.open_dataset("data_stream-oper_stepType-accum1.nc")
ds2 = xr.open_dataset("data_stream-oper_stepType-accum2.nc")
#ds3 = xr.open_dataset("data_stream-oper_stepType-accum3.nc")
#ds4 = xr.open_dataset("data_stream-oper_stepType-accum4.nc")
#ds5 = xr.open_dataset("data_stream-oper_stepType-accum5.nc")

#ds6 = xr.open_dataset("data_stream-oper_stepType-instant1.nc")
#ds7 = xr.open_dataset("data_stream-oper_stepType-instant2.nc")
#ds8 = xr.open_dataset("data_stream-oper_stepType-instant3.nc")
#ds9 = xr.open_dataset("data_stream-oper_stepType-instant4.nc")
#ds10 = xr.open_dataset("data_stream-oper_stepType-instant5.nc")


# Merge along an existing dimension (e.g., 'time', 'depth', or another dimension)
merged_ds1 = xr.concat([ds1, ds2, ds3, ds4, ds5], dim="valid_time")  # Change "time" to the desired dimension
merged_ds2 = xr.concat([ds6, ds7, ds8, ds9, ds10], dim="valid_time")
# Save the merged dataset to a new file
merged_ds1.to_netcdf("merged_file1.nc")
merged_ds2.to_netcdf("merged_file2.nc")


merged_ds = xr.merge([merged_ds1, merged_ds2])


merged_ds.to_netcdf("final_merged_file.nc")

# Open the NetCDF file
dataset = Dataset("final_merged_file.nc", mode='r')


# Print the variable names
print("Variables in the NetCDF file:")
print(dataset.variables.keys())

dataset = xr.open_dataset("final_merged_file.nc")
###### if denided permission use below two line code
#dataset = xr.open_dataset("final_merged_file.nc")
dataset.close()


# Select specific variables (e.g., 'temperature', 'humidity')
selected_data = dataset[['valid_time', 'latitude','longitude', 'ssrd','d2m', 't2m']]

# Convert to a pandas DataFrame
df = selected_data.to_dataframe().reset_index()

# Convert the 'utc_time' column to datetime
df['utc_time'] = pd.to_datetime(df['valid_time'])

back=df

# Convert UTC to CST
df['est_time'] = df['utc_time'].dt.tz_localize('UTC').dt.tz_convert('America/New_York')
df['est_time'] = df['est_time'].dt.tz_localize(None, ambiguous='NaT')

df['TIMESTAMP'] = df['est_time'].dt.floor('T')
df.drop(['valid_time','utc_time','expver','number','est_time'], axis=1, inplace=True)

chec=df[df.duplicated(subset=['TIMESTAMP'], keep=False)] # check for duplicated values 

df = df.drop_duplicates(subset=['TIMESTAMP'], keep='first') # remove duplicated values

check=df[df.duplicated(subset=['TIMESTAMP'], keep=False)]# check if removing duplicates worked 

back1=df

df['TIMESTAMP'] = pd.to_datetime(df['TIMESTAMP'])

# Set TIMESTAMP as the index
df.set_index('TIMESTAMP', inplace=True)
back2=df

# Now resample every 30 minutes and interpolate missing values
df = df.resample('30T').interpolate(method='linear')
df = df.reset_index()

back3=df

#df = df.drop_duplicates(subset='TIMESTAMP')
df.head()
df.tail()
df['ssrd_']=df.ssrd/3600
df['Tair_er']=df.t2m -273.15
df['d2m_']=df.d2m -273.15

#############################################################################
## calculate VPD
svpt=0.6108 * np.exp((17.27 * df['Tair_er']) / (df['Tair_er'] + 237.3))
svptd=0.6108 * np.exp((17.27 * df['d2m_']) / (df['d2m_'] + 237.3))
df['VPD_er'] = (svpt - svptd)*10 ## hecta pascal unit
era=df

#########################################################################################################

















































#########################################################################################################
## comparison between ERA and ameriflux data 

df=back3
exclude_months = [1,2,3,4,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]

df.index=np.arange(0, len(df))

check=df

df=check
df_clean = df.dropna(subset=['Rg', 'ssrd_'])
df=df_clean
df.index=np.arange(0, len(df))


correlation_matrix = np.corrcoef(df['Rg'], df['ssrd_'])
r_squared = correlation_matrix[0, 1] ** 2

r_squared






plt.plot(df.Rg,'r--')
plt.plot(df.ssrd_,'b--')
back4=df


plt.plot(df.Rg, color='blue', linestyle='--',label='Rg (ameri)')  # Set line color and label
plt.plot( df.ssrd_, color='black', linestyle='--',label='Rg (ERA5')  # Set line color and label
plt.xlabel('No. of obs')  # Label for x-axis
plt.ylabel('Rg (wm-2)')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()  

## see below for scatter plot
######################

df=check
df_clean = df.dropna(subset=['Tair', 'Tair_er'])
df=df_clean
df.index=np.arange(0, len(df))
                   
plt.plot(df.Tair, color='blue', linestyle='--',label='TA (ameri)')  # Set line color and label
plt.plot( df.Tair_er, color='black', linestyle='--',label='TA(ERA5')  # Set line color and label
plt.xlabel('No. of obs')  # Label for x-axis
plt.ylabel('Tair(C)')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()  

df=check
df_clean = df.dropna(subset=['VPD', 'VPD_er'])
df=df_clean
df.index=np.arange(0, len(df))
                   
                   
plt.plot(df.VPD, color='blue', linestyle='--',label='VPD (ameri)')  # Set line color and label
plt.plot( df.VPD_er, color='black', linestyle='--',label='VPD (ERA5)')  # Set line color and label
plt.xlabel('No. of obs')  # Label for x-axis
plt.ylabel('VPD (hPa)')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()

######################################################################################################
## scatter plot 

correlation_matrix = np.corrcoef(df['Rg'], df['ssrd_'])
r_squared = correlation_matrix[0, 1] ** 2

r_squared


# Create scatter plot
plt.figure(figsize=(6, 4))
plt.scatter(df['Rg'], df['ssrd_'], color='blue')

# Plot the line of best fit (optional, for visual reference)
m, b = np.polyfit(df['Rg'], df['ssrd_'], 1)
#plt.plot(df['Rg'], m*df['ssrd_'] + b, color='red', linestyle='--', label=f'Line of Best Fit')
z1 = np.polyfit(df['Rg'], df['ssrd_'], 1)
p1 = np.poly1d(z1)
plt.plot(df['Rg'], p1(df['Rg']),color='black')

# Add R² label to the plot
#plt.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')

# Customize plot
plt.title('R2=0.56')
plt.xlabel('Rg_ameri')
plt.ylabel('Rg_ERA5')
plt.legend()

# Show plot
plt.show()

#####################################################################################

correlation_matrix = np.corrcoef(df['Tair'], df['Tair_er'])
r_squared = correlation_matrix[0, 1] ** 2

r_squared


# Create scatter plot
plt.figure(figsize=(6, 4))
plt.scatter(df['Tair'], df['Tair_er'], color='blue')

# Plot the line of best fit (optional, for visual reference)
m, b = np.polyfit(df['Tair'], df['Tair_er'], 1)
#plt.plot(df['Rg'], m*df['ssrd_'] + b, color='red', linestyle='--', label=f'Line of Best Fit')
z1 = np.polyfit(df['Tair'], df['Tair_er'], 1)
p1 = np.poly1d(z1)
plt.plot(df['Tair'], p1(df['Tair']),color='black')

# Add R² label to the plot
#plt.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')

# Customize plot
plt.title('R2=0.56')
plt.xlabel('TA_ameri')
plt.ylabel('TA_ERA5')
plt.legend()

# Show plot
plt.show()

#####################################################################################################

correlation_matrix = np.corrcoef(df['VPD'], df['VPD_er'])
r_squared = correlation_matrix[0, 1] ** 2
r_squared

# Create scatter plot
plt.figure(figsize=(6, 4))
plt.scatter(df['VPD'], df['VPD_er'], color='blue')

# Plot the line of best fit (optional, for visual reference)
m, b = np.polyfit(df['VPD'], df['VPD_er'], 1)
#plt.plot(df['Rg'], m*df['ssrd_'] + b, color='red', linestyle='--', label=f'Line of Best Fit')
z1 = np.polyfit(df['VPD'], df['VPD_er'], 1)
p1 = np.poly1d(z1)
plt.plot(df['VPD'], p1(df['VPD_er']),color='black')

# Add R² label to the plot
#plt.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')

# Customize plot
plt.title('R2=0.17')
plt.xlabel('VPD_ameri (hpa)')
plt.ylabel('VPD_ERA5')
plt.legend()

# Show plot
plt.show()


###############################################################################################



plt.plot(df1.WUE, color='blue', linestyle='--',label='WUE (prev)')  # Set line color and label
plt.plot(  df2.WUE, color='black', linestyle='--',label='WUE (ERA5)')  # Set line color and label
plt.xlabel('year')  # Label for x-axis
plt.ylabel('WUE')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()

plt.plot(df1.CUE, color='blue', linestyle='--',label='CUE (prev)')  # Set line color and label
plt.plot( df2.CUE, color='black', linestyle='--',label='CUE (ERA5)')  # Set line color and label
plt.xlabel('year')  # Label for x-axis
plt.ylabel('CUE')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()


plt.plot(df1.GPP, color='blue', linestyle='--',label='GPP (prev)')  # Set line color and label
plt.plot( df2.GPP, color='black', linestyle='--',label='GPP (ERA5)')  # Set line color and label
plt.xlabel('year_month')  # Label for x-axis
plt.ylabel('GPP')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()

plt.plot(df1.reco, color='blue', linestyle='--',label='Reco (prev)')  # Set line color and label
plt.plot( df2.reco, color='black', linestyle='--',label='Reco (ERA5)')  # Set line color and label
plt.xlabel('year')  # Label for x-axis
plt.ylabel('Reco')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()


plt.plot(df1.NEE-df2.NEE, color='blue', linestyle='--',label='NEE diff')  # Set line color and label
#plt.plot( df2.GPP, color='black', linestyle='--',label='GPP (ERA5)')  # Set line color and label
plt.xlabel('year_month')  # Label for x-axis
plt.ylabel('NEE difference (gCm-2)')  # Label for y-axis
plt.title('US-Skr')  # Title of the plot
plt.legend()  # Display legend
plt.show()


#############################################################################################################
# Define the start and end dates for the entire range
start_date = "2004-01-01 00:00:00"
end_date = "2021-12-31 23:59:59"

# Generate half-hour intervals for all years
idx = pd.date_range(start=start_date, end=end_date, freq="30T")

#################################################################################################
back=df

df=back

## fix timestamp 
df["datetime"]=df['cst_time']
back=df
df=back
df.index = pd.DatetimeIndex(df.datetime)
#df["datetime"]=df['TIMESTAMP']
df = df[~df.index.duplicated(keep='first')] 
df=df.drop('datetime', axis=1)  
df = df.reindex(idx, fill_value=-9999) #add missing dates and add  nan for missing values
df= df.replace(-9999, np.nan, regex=True) 
df.isnull().values.any()

df=df.rename_axis('TIMESTAMP').reset_index()  

df["TIMESTAMP2"]=df['TIMESTAMP']
df=df.drop(columns=["TIMESTAMP"])


#################################################################################################
df['Rg_']=df.ssrd/3600

df['Tair']=df.t2m -273.15

plt.plot(df.Tair)

ta_A03=df

ta_A03.to_csv('ta_A03.csv', index=False, header=True)

