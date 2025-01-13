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

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("AMF_US-Skr_BASE_HH_2-5.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
back=df
# 2004-2023

df=back
df["datetime"]=df['TIMESTAMP']
idx=pd.date_range(start='01/01/2004 00:00:00', periods=350640, freq='30T')  # half hourly. 
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

k=pd.concat([s]*(7305),axis=0)
k.reset_index(drop=True, inplace=True)
#k.index = np.arange(1, len(k)+1)
df.insert(3,'Hour',k)
df=df.iloc[:,0:114]
df.head(10)
df.tail(10)
#try different names 
df["NEE"]=df["NEE_PI"]
#np.nanpercentile(df["NEE"],[1, 99]) # -21.6, 10
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
df["Tair"]=df["TA_1_1_1"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["RH"]=df["RH_1_1_1"]

df['NEE'] = df['NEE'].fillna(df['NEE_check'])

back3=df

#######################################################################################################
reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"],df["LE"],
   df["H"],df['Rg'],df["Tair"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

#only fill columns you need

reframed.to_csv('gaps_US-Skr.csv', index=False, header=True)

######################################################################################################


