# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:22:26 2024

@author: ammar
"""


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
df=pd.read_csv("fill_US-Skr.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
#years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2023]
years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2023, 2005,2011,2020,2021,2022]
# complete years: 2007,2008,2009,2010,2018,2019

# missing months 2005/8/15 (missing 15 days of aug)   
# missing months 2011/8/19 (missing 15 days of aug)  
# missing months 2018/4/01 (missing 20 days of april)  
# missing months 2020 (missing july, aug)  
# missing months 2021 (missing july, aug)  
# missing months 2022 (missing april,may,june)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]

df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
# start
plt.plot(df['NEE_f'],color='red')
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
df=df.reset_index()
df_Skr=df

###############################################################################################
df_Skr=df
# remove specific month of a year
#df = df[~((df['date'].dt.year == 2005) & (df['date'].dt.month == 8))]
#################################################################################################
# seasonal
df=df_Skr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df.to_csv('check.csv', index=False, header=True)


df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df['ID'] = 'US-Skr'
df['salin'] = 'med'

#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
df_Skr_m=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Skr_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_Skr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df['ID'] = 'US-Skr'
df['salin'] = 'med'
df_Skr_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Skr_y.csv', index=False, header=True)


#########################################################################################

#plots
df=df_Skr_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Skr:Shark River Slough, Growing season")

##################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
df=pd.read_csv("fill_US-TaS.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
years_to_remove = [2023]
# complete years: 2017-2022
# missing months 2023 (missing jun-aug)  
# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)

df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
df['ET_']=(df['LE_f']/2.43)*(1/1e6)*1800
df['ET'] = df['ET'].fillna(df['ET_'])
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
# start
plt.plot(df['NEE_f'],color='red')
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_Tas=df

#################################################################################################
# seasonal
df=df_Tas
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df['ID'] = 'US-Tas'
df['salin'] = 'high'

#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])

df_Tas_m=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Tas_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_Tas
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df['ID'] = 'US-Tas'
df['salin'] = 'high'
df_Tas_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Tas_y.csv', index=False, header=True)


###############################################################################################
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
df=pd.read_csv("fill_US-EvM.csv")  

df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

#years_to_remove = [2020]
# complete years: 2021-2023

# missing months 2020 (missing april)  

# Remove rows based on years
#df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
# start
#df.to_csv('check_Tas.csv', index=False, header=True)
#plt.plot(df['NEE_f'],color='red')
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_EvM=df

#################################################################################################
# seasonal
df=df_EvM
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-EvM'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
df_EvM_m=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_EvM_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_EvM
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df['ID'] = 'US-EvM'
df['salin'] = 'low'
df_EvM_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_EvM_y.csv', index=False, header=True)

############################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
df=pd.read_csv("fill_US-KS3.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2019]

# missing months 2018 (missing 10 days of april)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)

df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ks3=df

#################################################################################################
# seasonal
df=df_ks3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

df['ID'] = 'US-ks3'
df['salin'] = 'med'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])

df_ks3_m=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks3_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_ks3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

df['ID'] = 'US-ks3'
df['salin'] = 'med'
df_ks3_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks3_y.csv', index=False, header=True)

#######################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
df=pd.read_csv("fill_US-KS4.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2022]

# missing months 2017 (missing 17 days of april)  
# missing months 2019 (missing 3 days of aug)  
#missing months 2022 (only 5 days of data for april is present. missing everything else)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)

############################################################################################

## 3 days of aug missing in 2019, use interpolation
df['LE_f']=df.LE_f.interpolate() 
df['Tair_f']=df.Tair_f.interpolate() 
df['GPP_DT']=df.GPP_DT.interpolate() 
df['Reco_DT']=df.Reco_DT.interpolate() 
df['NEE_f']=df.NEE_f.interpolate() 
#df.to_csv('check.csv', index=False, header=True)


##########################################################################################
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ks4=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ks4=df
#df.to_csv('check.csv', index=False, header=True)
#################################################################################################
# seasonal
df=df_ks4
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df['ID'] = 'US-ks4'
df['salin'] = 'med'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])

df_ks4_m=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks4_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_ks4
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df['ID'] = 'US-ks4'
df['salin'] = 'med'

df_ks4_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks4_y.csv', index=False, header=True)



##############################################################################################

## california sites 

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-Tw1.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2011,2012]
# whole 2011 missing
# missing months 2012 (missing april, may june, 2 weeks of july)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_tw1=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_tw1=df

#################################################################################################
# seasonal
df=df_tw1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-tw1'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
df_tw1_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_tw1_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_tw1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df['ID'] = 'US-tw1'
df['salin'] = 'low'
df_tw1_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_tw1_y.csv', index=False, header=True)

############################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-Myb.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# Remove rows based on years
#df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_Myb=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_Myb=df

#################################################################################################
# seasonal
df=df_Myb
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-Myb'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
df_Myb_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Myb_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_Myb
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-Myb'
df['salin'] = 'low'
df_Myb_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Myb_y.csv', index=False, header=True)


#############################################################################################



os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-Srr.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# Remove rows based on years
#df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_srr=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_srr=df

#################################################################################################
# seasonal
df=df_srr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-srr'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
df_srr_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_srr_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_srr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-srr'
df['salin'] = 'low'
df_srr_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_srr_y.csv', index=False, header=True)

############################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-Dmg.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# Remove rows based on years

years_to_remove = [2021]

df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_dmg=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_dmg=df

#################################################################################################
# seasonal
df=df_dmg
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-dmg'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_dmg_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_dmg_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_dmg
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-dmg'
df['salin'] = 'low'
df_dmg_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_dmg_y.csv', index=False, header=True)
########################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-EDN.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# Remove rows based on years

years_to_remove = [2021]

df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_edn=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_edn=df

#################################################################################################
# seasonal
df=df_edn
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-edn'
df['salin'] = 'high'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_edn_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_edn_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_edn
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-edn'
df['salin'] = 'high'
df_edn_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_edn_y.csv', index=False, header=True)
##############################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-EKH.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# Remove rows based on years
# 2022 missing may and half june
# 2024 missing last 11 days of aug

years_to_remove = [2022,2024]

df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ekh=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ekh=df

#################################################################################################
# seasonal
df=df_ekh
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-ekh'
df['salin'] = 'high'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_ekh_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ekh_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_ekh
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-ekh'
df['salin'] = 'high'
df_ekh_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ekh_y.csv', index=False, header=True)
#########################################################################################3

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-EKP.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# whole 2021 missing 
# 2022 missing may and half june
# 2024 missing last 11 days of aug

years_to_remove = [2021, 2022,2024]

df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,4,9,10,11,12]

# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ekp=df

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_ekp=df

#################################################################################################
# seasonal
df=df_ekp
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-ekp'
df['salin'] = 'high'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_ekp_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ekp_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_ekp
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-ekp'
df['salin'] = 'high'
df_ekp_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ekp_y.csv', index=False, header=True)
#########################################################################################3

## compare sites
### yearly

df = pd.read_excel('analysis.xlsx', sheet_name='yearly')


plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.05       # the width of the bars

#### high
val1 = df["WUE_tas"]
rec1 = ax.bar(ind, val1, width, color='lightcoral')
val2 = df["WUE_edn"]
rec2 = ax.bar(ind+width, val2, width, color='indianred')
val3= df["WUE_ekh"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='red')
val4 = df["WUE_ekp"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='tomato')
###medium
val5 = df["WUE_ks4"]
rec5 = ax.bar(ind+width+0.5, val7, width, color='blue')
val6 = df["WUE_skr"]
rec6 = ax.bar(ind+width+0.6, val8, width, color='royalblue')
val7 = df["WUE_ks3"]
rec7 = ax.bar(ind+width+0.7, val9, width, color='cornflowerblue')

###low
val8 = df["WUE_evm"]
rec8 = ax.bar(ind+width+0.3, val5, width, color='forestgreen')
val9 = df["WUE_dmg"]
rec9 = ax.bar(ind+width+0.4, val6, width, color='limegreen')
val10 = df["WUE_myb"]
rec10 = ax.bar(ind+width+0.8, val10, width, color='lime')
val11 = df["WUE_srr"]
rec11 = ax.bar(ind+width+0.9, val11, width, color='seagreen')
val12 = df["WUE_tw1"]
rec12 = ax.bar(ind+width+1, val12, width, color='darkgreen')

ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010','2011','2012','2013','2014','2015',
                     '2016','2017','2018','2019','2020','2021','2022','2023'))
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0],rec6[0],rec7[0],rec8[0],rec9[0],rec10[0],
rec11[0] ,rec12[0]), ('US-Tas,H','US-EDN:H','US-EKH:H','US-EKP:H','US-KS4:M','US-Skr:M, For',
 'US-KS3:M', 'US-EvM:L','US-Dmg:L', 'US-Myb:L', 'US-Srr:L', 'US-Tw1:L',) ,ncol=3,loc='upper left', fontsize = 'large')
#ax.set_ylim(0,5)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida sites: May-Aug (total)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

###############################################################################################

plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.05       # the width of the bars

#### high
val1 = df["WUE_tas"]
rec1 = ax.bar(ind, val1, width, color='red')
val2 = df["WUE_edn"]
rec2 = ax.bar(ind+width, val2, width, color='red')
val3= df["WUE_ekh"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='red')
val4 = df["WUE_ekp"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='red')
###medium
val5 = df["WUE_ks4"]
rec5 = ax.bar(ind+width+0.5, val5, width, color='blue')
val6 = df["WUE_skr"]
rec6 = ax.bar(ind+width+0.6, val6, width, color='blue')
val7 = df["WUE_ks3"]
rec7 = ax.bar(ind+width+0.7, val7, width, color='blue')

###low
val8 = df["WUE_evm"]
rec8 = ax.bar(ind+width+0.3, val8, width, color='green')
val9 = df["WUE_dmg"]
rec9 = ax.bar(ind+width+0.4, val9, width, color='green')
val10 = df["WUE_myb"]
rec10 = ax.bar(ind+width+0.8, val10, width, color='green')
val11 = df["WUE_srr"]
rec11 = ax.bar(ind+width+0.9, val11, width, color='green')
val12 = df["WUE_tw1"]
rec12 = ax.bar(ind+width+1, val12, width, color='green')

ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010','2011','2012','2013','2014','2015',
                     '2016','2017','2018','2019','2020','2021','2022','2023'))
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0],rec6[0],rec7[0],rec8[0],rec9[0],rec10[0],
rec11[0] ,rec12[0]), ('US-Tas,H','US-EDN:H','US-EKH:H','US-EKP:H','US-KS4:M','US-Skr:M, For',
 'US-KS3:M', 'US-EvM:L','US-Dmg:L', 'US-Myb:L', 'US-Srr:L', 'US-Tw1:L',) ,ncol=3,loc='upper left', fontsize = 'large')
#ax.set_ylim(0,5)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Salinity effect on WUE: May-Aug (total)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()





###########################################################################################3



plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.05       # the width of the bars

#### florida
val1 = df["WUE_tas"]
rec1 = ax.bar(ind, val1, width, color='firebrick')
val2 = df["WUE_ks4"]
rec2 = ax.bar(ind+width+0.5, val2, width, color='red')
val3 = df["WUE_skr"]
rec3 = ax.bar(ind+width+0.6, val3, width, color='red')
val4 = df["WUE_ks3"]
rec4 = ax.bar(ind+width+0.7, val4, width, color='red')
val5 = df["WUE_evm"]
rec5 = ax.bar(ind+width+0.3, val5, width, color='coral')
###califor

val6 = df["WUE_edn"]
rec6 = ax.bar(ind+width, val6, width, color='blue')
val7= df["WUE_ekh"]
rec7 = ax.bar(ind+width+0.1, val7, width, color='blue')
val8 = df["WUE_ekp"]
rec8 = ax.bar(ind+width+0.2, val8, width, color='blue')
###low
val9 = df["WUE_dmg"]
rec9 = ax.bar(ind+width+0.4, val9, width, color='lightsteelblue')
val10 = df["WUE_myb"]
rec10 = ax.bar(ind+width+0.8, val10, width, color='lightsteelblue')
val11 = df["WUE_srr"]
rec11 = ax.bar(ind+width+0.9, val11, width, color='lightsteelblue')
val12 = df["WUE_tw1"]
rec12 = ax.bar(ind+width+1, val12, width, color='lightsteelblue')

ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010','2011','2012','2013','2014','2015',
                     '2016','2017','2018','2019','2020','2021','2022','2023'))
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0],rec6[0],rec7[0],rec8[0],rec9[0],rec10[0],
rec11[0] ,rec12[0]), ('US-Tas,H', 'US-KS4:M','US-Skr:M, For','US-KS3:M', 'US-EvM:L',
 'US-EDN:H','US-EKH:H','US-EKP:H',
'US-Dmg:L', 'US-Myb:L', 'US-Srr:L', 'US-Tw1:L',) ,ncol=3,loc='upper left', fontsize = 'large')
#ax.set_ylim(0,3)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida and California sites: May-Aug (total)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

#########################################################################3


plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.05       # the width of the bars

#### florida
val1 = df["WUE_myb"]
rec1 = ax.bar(ind, val1, width, color='firebrick')
val2 = df["WUE_tw1"]
rec2 = ax.bar(ind+width+0.5, val2, width, color='red')

ax.set_xticks(ind+width)
ax.set_xticklabels( ('2011','2012','2013','2014','2015',
                     '2016','2017','2018','2019','2020','2021','2022','2023'))
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0],rec6[0],rec7[0],rec8[0],rec9[0],rec10[0],
rec11[0] ,rec12[0]), ( 'US-Myb:L', 'US-Tw1:L',) ,ncol=3,loc='upper left', fontsize = 'large')
#ax.set_ylim(0,3)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida and California sites: May-Aug (total)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

################################################################################################

plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total CUE', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.1       # the width of the bars
val1 = df["CUE_evm"]
rec1 = ax.bar(ind, val1, width, color='g')
val2 = df["CUE_tas"]
rec2 = ax.bar(ind+width, val2, width, color='orange')
val3= df["CUE_ks3"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='blue')
val4 = df["CUE_ks4"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='red')
val5 = df["CUE_skr"]
rec5 = ax.bar(ind+width+0.3, val5, width, color='black')


ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010',
                     '2017','2018','2019','2020','2021','2022','2023'
) )
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0]), ('US-EVM:For, 9 ppt','US-Tas:34 ppt','US-KS3: 23 ppt','US-KS4: 23 ppt','US-Skr: 17 ppt') ,loc='upper left', fontsize = 'large')
#ax.set_ylim(0,2.3)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida sites: April-Aug",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

#################################################################################################



df = pd.read_excel('analysis.xlsx', sheet_name='max')



plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.1       # the width of the bars

#### high
val1 = df["WUE_tas"]
rec1 = ax.bar(ind, val1, width, color='red')
val2 = df["WUE_edn"]
rec2 = ax.bar(ind+width, val2, width, color='red')
val3= df["WUE_ekh"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='red')
val4 = df["WUE_ekp"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='red')
val5 = df["WUE_evm"]
rec5 = ax.bar(ind+width+0.3, val5, width, color='red')
val6 = df["WUE_dmg"]
rec6 = ax.bar(ind+width+0.4, val6, width, color='red')

###medium
val7 = df["WUE_myb"]
rec7 = ax.bar(ind+width+0.5, val7, width, color='blue')
val8 = df["WUE_srr"]
rec8 = ax.bar(ind+width+0.6, val8, width, color='blue')
val9 = df["WUE_tw1"]
rec9 = ax.bar(ind+width+0.7, val9, width, color='blue')

###low
val10 = df["WUE_ks4"]
rec10 = ax.bar(ind+width+0.8, val10, width, color='green')
val11 = df["WUE_skr"]
rec11 = ax.bar(ind+width+0.9, val11, width, color='green')
val12 = df["WUE_ks3"]
rec12 = ax.bar(ind+width+1, val12, width, color='green')


ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010','2011','2012','2013','2014','2015',
                     '2016','2017','2018','2019','2020','2021','2022','2023'))
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0],rec6[0],rec7[0],rec8[0],rec9[0],rec10[0],
rec11[0] ,rec12[0]), ('US-Tas,H','US-EDN:H','US-EKH:H','US-EKP:H','US-KS4:M','US-Skr:M, For',
 'US-KS3:M', 'US-EvM:L','US-Dmg:L', 'US-Myb:L', 'US-Srr:L', 'US-Tw1:L',) ,loc='upper left', fontsize = 'large')
ax.set_ylim(0,5)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida sites: April-Aug (max)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

########################################################################################3

plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('maximum CUE', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.1       # the width of the bars
val1 = df["CUE_evm"]
rec1 = ax.bar(ind, val1, width, color='g')
val2 = df["CUE_tas"]
rec2 = ax.bar(ind+width, val2, width, color='orange')
val3= df["CUE_ks3"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='blue')
val4 = df["CUE_ks4"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='red')
val5 = df["CUE_skr"]
rec5 = ax.bar(ind+width+0.3, val5, width, color='black')


ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010',
                     '2017','2018','2019','2020','2021','2022','2023'
) )
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0]), ('US-EVM:For, 9 ppt','US-Tas:34 ppt','US-KS3: 23 ppt','US-KS4: 23 ppt','US-Skr: 17 ppt') ,loc='upper left', fontsize = 'large')
ax.set_ylim(0,1.4)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida sites: April-Aug (max)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

############################################################################################



























####################################################################################################
plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-EDN:Eden Landing Ecological Reserve, wetland,CA, 30–35 ppt",fontsize=8)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=100))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()


















































































######################################################################################################
## fluxnet 

###################################################################################################
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-EDN.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

#np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -1.82) & (df['LE_F_MDS'] <= 202)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_NT_VUT_REF'],color='green')
df['GPP_NT_VUT_REF'].isnull().sum()

#np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -5) & (df['GPP_NT_VUT_REF'] <= 15)
#mask = (df['GPP_DT_VUT_REF'] >= -5)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_DT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_DT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum
df=df.reset_index()

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()

plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_EDN=df


# seasonal
df=US_EDN
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df=df.reset_index()
df["WUE"]=df["GPP_g"]/df["ET"]
US_EDN_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-EDN:Eden Landing Ecological Reserve, wetland,CA, 30–35 ppt",fontsize=8)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=100))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

#####################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-NGB.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

#np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -1.82) & (df['LE_F_MDS'] <= 202)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
plt.plot(df['LE'],color='green')
df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_NGB=df

# seasonal
df=US_NGB
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
df=df.reset_index()

#np.percentile(df["WUE"], [1,99])
#mask = (df['WUE'] >= -28) & (df['WUE'] <= 74)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["WUE"]=df['WUE'].interpolate()
plt.plot(df['WUE'],color='green')
df=df.reset_index()
US_NGB_m=df

plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Year')  # we already handled the x-label with ax1
ax.set_title("US-NGB: NGEE Arctic Barrow, Tundra",fontsize=8)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=500))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()


#####################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-ASH.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])

df['LE_F_MDS']=np.where(df['LE_F_MDS']<=-2 , np.nan, df['LE_F_MDS'])
df['LE_F_MDS']=np.where(df['LE_F_MDS']>=1000 , np.nan, df['LE_F_MDS'])
df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_ASH=df

#np.percentile(df["WUE"], [1,99])
#mask = (df['WUE'] >= -28) & (df['WUE'] <= 74)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))


# seasonal
df=US_ASH
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
US_ASH_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Year')  # we already handled the x-label with ax1
ax.set_title("US-ASH: USSL San Joaquin Valley Almond High Salinity,CA, 0.65-0.96 ppt",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=100))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

####################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-Srr.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)


df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -2) & (df['LE_F_MDS'] <= 250)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
#plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_Srr=df

#np.percentile(df["WUE"], [1,99])
#mask = (df['WUE'] >= -0.1) & (df['WUE'] <= 5)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['WUE'],color='green')
#US_Srr=df


# seasonal
df=US_Srr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
US_Srr_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-Srr: Suisun marsh - Rush Ranch, 4.6 (0.25–12.3) PSU, CA",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=200))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

#################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-ASM.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
mask = (df['LE_F_MDS'] >= -2) & (df['LE_F_MDS'] <= 250)
df = df.loc[mask]
df.index = np.arange(0, len(df))
df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
#plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
#plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_ASM=df


# seasonal
df=US_ASM
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
np.percentile(df["WUE"], [1,99])
#mask = (df['WUE'] >= -0.1) & (df['WUE'] <= 5)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
plt.plot(df['WUE'],color='green')
df=df.reset_index()
US_ASM_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-ASM: USSL San Joaquin Valley Almond Medium Salinity, (0.91- 1.84), Deciduous Broadleaf Forests",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=200))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

########################################################################################################
#ammaramd

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-PSH.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -2) & (df['LE_F_MDS'] <= 250)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
#plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df["WUE"]=df["GPP_g"]/df["ET"]
df=df.reset_index()
US_PSH=df

# seasonal
df=US_PSH
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
US_PSH_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-PSH: USSL San Joaquin Valley Pistachio High, 2016-present (1.18-2.5), Deciduous Broadleaf Forests",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=200))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

#########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-PSL.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -2) & (df['LE_F_MDS'] <= 250)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
#plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_PSL=df

#np.percentile(df["WUE"], [1,99])
#mask = (df['WUE'] >= -0.1) & (df['WUE'] <= 5)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))


df=US_PSL
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
plt.plot(df['WUE'],color='green')
US_PSL_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-PSL: USSL San Joaquin Valley Pistachio Low,(0.545-1.07), Deciduous Broadleaf Forests, fluxnet",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=200))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

#####################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-HB1.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -2) & (df['LE_F_MDS'] <= 250)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
#plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_HB1=df


df=US_HB1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
plt.plot(df['WUE'],color='green')
US_HB1_m=df


plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-HB1: North Inlet Crab Haul Creek, wetland, 32 PSU",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=50))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

#################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-KS3.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
#mask = (df['LE_F_MDS'] >= -2) & (df['LE_F_MDS'] <= 250)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
#plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_KS3=df


df=US_KS3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
plt.plot(df['WUE'],color='green')
US_KS3_m=df


np.percentile(df["WUE"], [1,99])
#mask = (df['WUE'] >= -0.1) & (df['WUE'] <= 10)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
plt.plot(df['WUE'],color='green')

plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-KS3: Kennedy Space Center (salt marsh), wetlands, >6ppt",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=50))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()

########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("US-StJ.csv")  
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['TA_F']+273.16))*1e-6
plt.plot(df['lambda'],color='green')
plt.plot(df['LE_F_MDS'],color='green') #LE

np.percentile(df["LE_F_MDS"], [1,99])
mask = (df['LE_F_MDS'] >= -100) & (df['LE_F_MDS'] <= 550)
df = df.loc[mask]
df.index = np.arange(0, len(df))
df["LE"]=df['LE_F_MDS'].interpolate()

df["LE"]=df['LE_F_MDS']
plt.plot(df['LE'],color='green')

df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()

plt.plot(df['GPP_DT_VUT_REF'],color='green')
df['GPP_DT_VUT_REF'].isnull().sum()

np.percentile(df["GPP_DT_VUT_REF"], [1,99])
#mask = (df['GPP_DT_VUT_REF'] >= -3) & (df['GPP_NT_VUT_REF'] <= 15)
#df = df.loc[mask]
#df.index = np.arange(0, len(df))
#plt.plot(df['GPP_NT_VUT_REF'],color='green')
#df["GPP"]=df['GPP_NT_VUT_REF'].interpolate()

df["GPP"]=df['GPP_DT_VUT_REF']
plt.plot(df['GPP'],color='green')
df["GPP_g"]=(df["GPP"])*((12/10**6)*1800)
plt.plot(df['GPP_g'],color='green')
## daily sum

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
plt.plot(df['GPP_g'],color='green')
df=df.reset_index()
US_StJ=df


df=US_StJ
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['month']).sum()
df["WUE"]=df["GPP_g"]/df["ET"]
plt.plot(df['WUE'],color='green')
df=df.reset_index()
plt.plot(df['WUE'],color='green')
US_StJ_m=df


np.percentile(df["WUE"], [1,99])
mask = (df['WUE'] >= -10) & (df['WUE'] <= 35)
df = df.loc[mask]
df.index = np.arange(0, len(df))
plt.plot(df['WUE'],color='green')

plt.rcParams['figure.figsize'] = (4, 4)
fig, ax = plt.subplots()
plt.plot(df["TIMESTAMP"], df['WUE'],color='green')
plt.ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
#plt.plot(df['WUE'], color='crimson', ls='-',marker='o',markersize=5 )
ax.tick_params(axis='y', labelcolor='black')
#ax.legend()
#ax.legend(loc='upper left', fontsize = '9')
ax.xaxis.label.set_size(10)
ax.yaxis.label.set_size(9) 
ax.set_ylabel('Total Annual WUE (g C per kg $H_{2}$O)')  # we already handled the x-label with ax1
ax.set_xlabel('Water Year')  # we already handled the x-label with ax1
ax.set_title("US-StJ: St Jones Reserve, tidal, 10 (0.3–23), wetlands",fontsize=7)
ax.tick_params(axis = 'x', which = 'major', labelsize = 10)
ax.tick_params(axis = 'y', which = 'major', labelsize =10)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=200))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
fig.autofmt_xdate() 
#plt.ticklabel_format(useOffset=False)
plt.show()
fig.savefig('testfig.jpg',dpi=600,bbox_inches="tight")
plt.tight_layout()






























#######################################################################################################
back=df

df=back
X=df["X"]
Y=df["Y"]

mask = (df['X'] >= 0.968) & (df['X'] <= 50.2627)
df = df.loc[mask]
df.index = np.arange(0, len(df))

df['ET']=(df['LE_F_MDS']/df['lambda'])*(1/1e6)*1800

plt.plot(df['ET'],color='green')

back=df

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
df.tail(10)

df["NEE"]=df["NEE_PI"]
df["Rg"]=df["SW_IN"]
df["Tair"]=df["TA_1_1_1"]
df["VPD"]=df["VPD_PI"]
df["Ustar"]=df["USTAR"]
df["RH"]=df["RH_1_1_1"]

reframed=pd.concat((df["DateTime"],df['Year'],df['DoY'],df['Hour'],df["NEE"],df["LE"],df["H"],df['Rg'],
                    df["Tair"],df["RH"],df["VPD"],df["Ustar"],df["PA"],df["NETRAD"]),axis=1)

reframed.to_csv('gaps.csv', index=False, header=True)

######################################################################################################


df['lambda']=(3149000-2370*(df['Tair']+273.16))*1e-6

plt.plot(df['lambda'],color='green')
df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')


s=df.groupby([df['TIMESTAMP'].dt.date]).sum()


plt.plot(s['ET'],color='red')
plt.ylabel("Daily ET")
plt.title("US-Skr:Shark River Slough")



df=s
df=df.reset_index()
df['Year'] = pd.to_datetime(df["TIMESTAMP"]).dt.year
EMS_daily_ET=df
df.isnull().values.any()

df=EMS_daily_ET
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.resample('m', on='TIMESTAMP').sum()
df=df.reset_index()
df['Year'] = df['TIMESTAMP'].dt.year
df['Year'] = df['TIMESTAMP'].dt.month
EMS_total_month_ET=df

df=EMS_daily_ET
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.resample('Y', on='TIMESTAMP').sum()
df=df.reset_index()
df['Year'] = df['TIMESTAMP'].dt.year
EMS_total_annual_ET=df

df=EMS_daily_ET
df=df.resample('Y', on='TIMESTAMP').mean()
df=df.reset_index()
df['Year'] = df['TIMESTAMP'].dt.year
EMS_mean_annual_ET=df

df=df.rename_axis('TIMESTAMP').reset_index()  
plt.plot(df["Year"], df['ET'],color='red')
plt.ylabel("Annual ET")
plt.title("US-Skr:Shark River Slough")











#######################################################################################################
df['TIMESTAMP1'] = pd.to_datetime(df["TIMESTAMP"], format='%Y-%m-%d') # 1964-2021
df.head(10)
df.tail(10)


df['TIMESTAMP']=back_up_HEM["TIMESTAMP"]
fill_hem_2022=df


plt.plot(df['LE'],color='green')
reddy=fill_hem
#df["LE1"] = np.where(df["LE"]<-200, np.NaN, df["LE"])
#df["LE2"] = np.where(df["LE1"]>800, np.NaN, df["LE1"])
#df["LE"]=df["LE2"]
plt.plot(df['LE'],color='green')
df['lambda']=(3149000-2370*(df['Tair']+273.16))*1e-6
plt.plot(df['lambda'],color='green')


df['ET']=(df['LE']/df['lambda'])*(1/1e6)*1800
