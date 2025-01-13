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
######################################################################################################################################################################################
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


######################################################################################################################################################################################

#plots
df=df_Skr_m
plt.bar(df["year"], df['WUE'],color='blue')
plt.ylabel("WUE")
plt.title("US-Skr:Shark River Slough, Growing season")

########################################################################################################################################################################################


