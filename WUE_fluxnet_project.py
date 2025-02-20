# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:22:26 2024

@author: ammar
"""


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("fill_US-Skr.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
#years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2023]
years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2023, 2005,2011,2018,2020,2021,2022]
# complete years: 2007,2008,2009,2010,2019

# missing months 2005/8/15 (missing 15 days of aug)   
# missing months 2011/8/19 (missing 15 days of aug)  
# missing months 2018/4/01 (missing 20 days of april)  
# missing months 2020 (missing july, aug)  
# missing months 2021 (missing july, aug)  
# missing months 2022 (missing april,may,june)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,9,10,11,12]
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
##########################################################################
df_Skr=df
# remove specific month of a year
#df = df[~((df['date'].dt.year == 2005) & (df['date'].dt.month == 8))]
#################################################################################################
# seasonal
df=df_Skr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

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
df_Skr_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Skr_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_Skr
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()
df_Skr_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Skr_ym.csv', index=False, header=True)
#########################################################################################

#plots
df=df_Skr_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Skr:Shark River Slough, Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-Skr:Shark River Slough, Growing season")


df=df_Skr_y
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Skr:Shark River Slough, yearly")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-Skr:Shark River Slough, yearly")


df=df_Skr_ym
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Skr:Shark River Slough, max")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-Skr:Shark River Slough, max")



##################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("fill_US-TaS.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
years_to_remove = [2023]
# complete years: 2017-2022
# missing months 2023 (missing jun-aug)  
# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,9,10,11,12]
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
df_Tas_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Tas_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_Tas
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_Tas_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Tas_ym.csv', index=False, header=True)

###################################################################################################
#plots
##############################################################################################

#plots
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_Tas_ym=df


df=df_Tas_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Tas:Taylor Slough/Panhandle, Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-Tas:Taylor Slough/Panhandle, Growing season")


df=df_Tas_y
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Tas:Taylor Slough/Panhandle, yearly")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-TasTaylor Slough/Panhandle, yearly")


df=df_Tas_ym
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Tas:Taylor Slough/Panhandle, max")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-Tas:Taylor Slough/Panhandle, max")

###############################################################################################
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("fill_US-EvM.csv")  

df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2020]
# complete years: 2021-2023

# missing months 2020 (missing april)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,9,10,11,12]
# remove rows 
df = df[~df['TIMESTAMP'].dt.month.isin(exclude_months)]
df.index=np.arange(0, len(df))
back=df
df.head(10)
df.tail(10)
df['lambda']=(3149000-2370*(df['Tair_f']+273.16))*1e-6
plt.plot(df['lambda'],color='green')

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
#df['ET_']=(df['LE_f']/2.43)*(1/1e6)*1800
#df['ET'] = df['ET'].fillna(df['ET_'])
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
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_EvM_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_EvM_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_EvM
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_EvM_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_EvM_ym.csv', index=False, header=True)

################################################################################################

df=df_EvM_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-EvM:Everglades Saltwater intrusion marsh, Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-EvM:Everglades Saltwater intrusion marsh, Growing season")


########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("fill_US-KS3.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2019]

# missing months 2018 (missing 10 days of april)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,9,10,11,12]
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

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

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
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_ks3_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks3_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_ks3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_ks3_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks3_ym.csv', index=False, header=True)

################################################################################################

df=df_ks3_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-KS3:US-KS3: Kennedy Space Center (salt marsh), Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-KS3:US-KS3: Kennedy Space Center (salt marsh), Growing season")


#######################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data")
df=pd.read_csv("fill_US-KS4.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2017,2022]

# missing months 2017 (missing 17 days of april)  
# missing months 2019 (missing 3 days of aug)  
#missing months 2022 (only 5 days of data for april. missing everything else)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
exclude_months = [1,2,3,9,10,11,12]
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

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

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
plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_ks4_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks4_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_ks4
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_ks4_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_ks4_ym.csv', index=False, header=True)

################################################################################################

df=df_ks4_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-ks4: Kennedy Space Center (Spartina marsh), Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-ks4: Kennedy Space Center (Spartina marsh), Growing season")

##############################################################################################

## california sites 

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-Tw1.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

years_to_remove = [2011]

# missing months 2017 (missing 17 days of april)  
# missing months 2019 (missing 3 days of aug)  
#missing months 2022 (only 5 days of data for april. missing everything else)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
#exclude_months = [1,2,3,9,10,11,12]

exclude_months = [1,9,10,11,12]
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
df_tw1_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_tw1_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_tw1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_tw1_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_tw1_ym.csv', index=False, header=True)

################################################################################################

df=df_tw1_ym
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-tw1: Kennedy Space Center (Spartina marsh), Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-tw1: Kennedy Space Center (Spartina marsh), Growing season")
###############################################################################################

##################################################################################################


plt.plot(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-tw1: Twitchell Wetland West Pond, Growing season")


plt.plot(df["year"], df['CUE'],color='blue')
plt.ylabel("CUE")
plt.title("US-tw1: Twitchell Wetland West Pond, Growing season")

############################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-Myb.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

#years_to_remove = [2011]

# missing months 2017 (missing 17 days of april)  
# missing months 2019 (missing 3 days of aug)  
#missing months 2022 (only 5 days of data for april. missing everything else)  

# Remove rows based on years
df = df[~df['TIMESTAMP'].dt.year.isin(years_to_remove)]
#exclude_months = [1,2,3,9,10,11,12]

exclude_months = [1,2,3,9,10,11,12]
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
df_Myb_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Myb_y.csv', index=False, header=True)


##########################################################################
## max by year
df=df_Myb
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]
df= df.groupby(['year']).max()
df=df.reset_index()

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df_Myb_ym=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Myb_ym.csv', index=False, header=True)

################################################################################################

df=df_Myb_y
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Myb: Kennedy Space Center (Spartina marsh), Growing season")


plt.bar(df["year"], df['CUE'],color='green')
plt.ylabel("CUE")
plt.title("US-Myb: Kennedy Space Center (Spartina marsh), Growing season")
###############################################################################################

##################################################################################################
df=df_Myb_ym

plt.plot(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Myb: Mayberry Wetland, Growing season")


plt.plot(df["year"], df['CUE'],color='blue')
plt.ylabel("CUE")
plt.title("US-Myb: Mayberry Wetland, Growing season")

x1=df['ET']
y1=df['GPP'].astype(float)
#(x1, y1)
mk.original_test(y1)

############################################################################################



























################################################################################################

## compare sites
### yearly

df = pd.read_excel('analysis.xlsx', sheet_name='yearly')


plt.bar(df["year"], df['WUE_evm'],color='blue')
plt.ylabel("WUE")
plt.title("US-ks4: Kennedy Space Center (Spartina marsh), Growing season")


plt.rcParams['figure.figsize'] = (12, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.1       # the width of the bars
val1 = df["WUE_evm"]
rec1 = ax.bar(ind, val1, width, color='g')
val2 = df["WUE_tas"]
rec2 = ax.bar(ind+width, val2, width, color='orange')
val3= df["WUE_ks3"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='blue')
val4 = df["WUE_ks4"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='red')
val5 = df["WUE_skr"]
rec5 = ax.bar(ind+width+0.3, val5, width, color='black')


ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010',
                     '2017','2018','2019','2020','2021','2022','2023'
) )
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0]), ('US-EVM:For, 9 ppt','US-Tas:34 ppt','US-KS3: 23 ppt','US-KS4: 23 ppt','US-Skr: 17 ppt') ,loc='upper left', fontsize = 'large')
ax.set_ylim(0,2.3)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Florida sites: April-Aug",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()

###############################################################################################


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
ax.set_ylabel('Maximum WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.1       # the width of the bars
val1 = df["WUE_evm"]
rec1 = ax.bar(ind, val1, width, color='g')
val2 = df["WUE_tas"]
rec2 = ax.bar(ind+width, val2, width, color='orange')
val3= df["WUE_ks3"]
rec3 = ax.bar(ind+width+0.1, val3, width, color='blue')
val4 = df["WUE_ks4"]
rec4 = ax.bar(ind+width+0.2, val4, width, color='red')
val5 = df["WUE_skr"]
rec5 = ax.bar(ind+width+0.3, val5, width, color='black')


ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010',
                     '2017','2018','2019','2020','2021','2022','2023'
) )
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0]), ('US-EVM:For, 9 ppt','US-Tas:34 ppt','US-KS3: 23 ppt','US-KS4: 23 ppt','US-Skr: 17 ppt') ,loc='upper left', fontsize = 'large')
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
