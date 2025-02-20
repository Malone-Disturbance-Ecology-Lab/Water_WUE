# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:22:26 2024

@author: ammar
"""
######## function to find negative values in GPP

def percentage_below_zero(df, column_name):
   
    valid_data = df[column_name].dropna()  # Remove NaN values
    total_count = len(valid_data)
    
    if total_count == 0:
        return 0
    
    count_below_zero = (valid_data < 0).sum()
    return (count_below_zero / total_count) * 100


result = percentage_below_zero(df, 'values')
print(f"Percentage of values below zero: {result:.2f}%")




#os.chdir(r"C:\ammara_MD\a_yale\Research Project\previous\data\florida")

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
df=pd.read_csv("fill_US-Skr.csv") 
#df=pd.read_csv("gaps_US-Skr.csv")  
#df=pd.read_csv("fill_US-Skr.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])
years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2020,2021, 2022, 2023] # remove 8 years
#years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2023, 2005,2011]# remove 11 years 
#years_to_remove = [2004,2006,2012,2013,2014,2015,2016,2017,2023, 2005,2011,2020,2021,2022]
# complete years: 2007,2008,2009,2010,2018,2019

# missing months 2005/8/15 (missing 15 days of aug)   
# missing months 2011/8/19 (missing 15 days of aug)  
 
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

df_Skr_h=df
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)

# count %age of data with negative GPP values
percentage_below_zero(df, 'GPP_DT')


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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)



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
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')
df['ID'] = 'US-Skr'
df['salin'] = 'med'

df1=df
df.to_csv('check2.csv', index=False, header=True)
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df_Skr_m=df
df_Skr_m_pr=df
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_Skr_m_pr.csv', index=False, header=True)
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


plt.plot(df.DateTime,df.VPD_PI)
plt.ylabel("VPD (kpa)")
plt.title("US-EvM: Everglades Saltwater intrusion marsh")
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
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))

df_Tas_h=df
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)

# count %age of data with negative GPP values
percentage_below_zero(df, 'GPP_DT')

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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)

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

df=df_Tas_m

#plt.plot(df.year, df.GPP, color='blue', marker='o', linestyle='--',label='GPP_Rg (-ve)')  # Set line color and label
plt.plot(df.month, df.GPP, color='blue', marker='o', linestyle='--',label='GPP_Rg (-ve)')  # Set line color and label

plt.xlabel('Year')  # Label for x-axis
plt.ylabel('GPP (gCm-2)')  # Label for y-axis
plt.title('GPP Over the Years')  # Title of the plot
plt.legend()  # Display legend
plt.show()  


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

df_EvM_h=df
df.to_csv('check.csv', index=False, header=True)


plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)

# count %age of data with negative GPP values
percentage_below_zero(df, 'GPP_DT')


# start

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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)

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

plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)
df_KS3_h=df
df.to_csv('check.csv', index=False, header=True)

# count %age of data with negative GPP values
percentage_below_zero(df, 'GPP_DT')
percentage_below_zero(df, 'Rg_f')
percentage_below_zero(df, 'Reco_DT')

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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)


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

years_to_remove = [2020, 2022]

# dont remove 2020 for monthly data 

# missing months 2017 (missing 17 days of april)  
# missing months 2019 (missing 3 days of aug)  
# missing months 2020 whole month of aug for 2020
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
df.to_csv('check.csv', index=False, header=True)
############################################################################################

## 3 days of aug missing in 2019, use interpolation
#2017- 2022
df['LE_f']=df.LE_f.interpolate() 
df['Tair_f']=df.Tair_f.interpolate() 
df['GPP_DT']=df.GPP_DT.interpolate() 
df['Reco_DT']=df.Reco_DT.interpolate() 
df['NEE_f']=df.NEE_f.interpolate() 
df['Rg_f']=df.NEE_f.interpolate() 


plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)
df_KS4_h=df
# count %age of data with negative GPP values
percentage_below_zero(df, 'GPP_DT')
percentage_below_zero(df, 'Rg_f')
percentage_below_zero(df, 'Reco_DT')



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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)

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

#df_ks4_y=df



df_ks4_rgp=df
#df_ks4_par=df
#df_ks4_y=df


plt.plot(df_ks4_y.year, df_ks4_y.GPP, color='blue', marker='o', linestyle='--',label='GPP_Rg-')  # Set line color and label
plt.plot(df_ks4_rgp.year, df_ks4_rgp.GPP, color='black', marker='o', linestyle='--',label='GPP_Rg+')  # Set line color and label

#plt.plot(df_ks4_par.year, df_ks4_par.GPP, color='green', marker='o', linestyle='--',label='GPP_PAR+')  # Set line color and label

#plt.plot(df.month, df.GPP, color='blue', marker='o', linestyle='--',label='GPP_Rg (-ve)')  # Set line color and label

plt.xlabel('Year')  # Label for x-axis
plt.ylabel('GPP (gCm-2)')  # Label for y-axis
plt.title('GPP Over the Years')  # Title of the plot
plt.legend()  # Display legend
plt.show()  


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


plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)
df_Tw1_h=df
df.to_csv('check.csv', index=False, header=True)


percentage_below_zero(df, 'GPP_DT')
percentage_below_zero(df, 'Rg_f')
percentage_below_zero(df, 'Reco_DT')

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
plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)


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

df["Rg_f"]=df["Rg"]
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)
df_Myb_h=df

percentage_below_zero(df, 'GPP_DT')
percentage_below_zero(df, 'Rg_f')
percentage_below_zero(df, 'Reco_DT')



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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)

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
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)
df_srr_h=df

percentage_below_zero(df, 'GPP_DT')
percentage_below_zero(df, 'Rg_f')
percentage_below_zero(df, 'Reco_DT')


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
plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)
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
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)

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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)
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
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)


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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)
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

#########################################################################################3

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\california")
df=pd.read_csv("fill_US-EKP.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# whole 2021 missing 
# 2022 missing may and half june
# 2023 first 13 days of may Rg is missing, used PI GPP recalculate it because first 13 days of may GPP 
# is zero but reco is present 
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
df.to_csv('check.csv', index=False, header=True)


plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)

df.to_csv('check.csv', index=False, header=True)

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

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)
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
#####################################################################################################


### include EKH




########################################################################################33
os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\massachusetts")
df=pd.read_csv("fill_US-PHM.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


# 2017,2020 few weeks of GPP missing 

#years_to_remove = [2017,2022] # if don't want to replace GPP_dt with GPP_nt

#2020, 2017 missing rg 
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

df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f)
plt.plot(df.Rg_f)
df_Phm_h=df


percentage_below_zero(df, 'GPP_DT')
percentage_below_zero(df, 'Rg_f')
percentage_below_zero(df, 'Reco_DT')



df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_phm=df

df['GPP_DT']=df['GPP_DT'].fillna(df['GPP_nt'])
df['GPP_DT'].isnull().any()




#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_phm=df


plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)
#################################################################################################
# seasonal
df=df_phm
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-phm'
df['salin'] = 'med'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_phm_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_phm_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_phm
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-phm'
df['salin'] = 'med'

#df_phm_par=df
df_phm_rgp=df
#df_phm_y=df


plt.plot(df_phm_y.year, df_phm_y.GPP, color='blue', marker='o', linestyle='--',label='GPP_Rg-')  # Set line color and label
plt.plot(df_phm_rgp.year, df_phm_rgp.GPP, color='black', marker='o', linestyle='--',label='GPP_Rg+')  # Set line color and label

plt.plot(df_phm_par.year, df_phm_par.GPP, color='green', marker='o', linestyle='--',label='GPP_PAR+')  # Set line color and label

#plt.plot(df.month, df.GPP, color='blue', marker='o', linestyle='--',label='GPP_Rg (-ve)')  # Set line color and label

plt.xlabel('Year')  # Label for x-axis
plt.ylabel('GPP (gCm-2)')  # Label for y-axis
plt.title('GPP Over the Years')  # Title of the plot
plt.legend()  # Display legend
plt.show()  


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_phm_y.csv', index=False, header=True)
#####################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\south carolina")
df=pd.read_csv("fill_US-HB1.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


# 2017,2020 few weeks of GPP missing 

#years_to_remove = [2017,2022] # if don't want to replace GPP_dt with GPP_nt

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
df.to_csv('check.csv', index=False, header=True)


plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)

df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hb1=df

#df['GPP_DT']=df['GPP_DT'].fillna(df['GPP_nt'])
#df['GPP_DT'].isnull().any()

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hb1=df

plt.plot(df.NEE)
plt.plot(df.GPP)
plt.plot(df.reco)
#################################################################################################
# seasonal
df=df_hb1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hb1'
df['salin'] = 'med'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_hb1_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hb1_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_hb1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hb1'
df['salin'] = 'med'
df_hb1_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hb1_y.csv', index=False, header=True)

####################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\south carolina")
df=pd.read_csv("fill_US-HB2.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


# 2017,2020 few weeks of GPP missing 

#years_to_remove = [2017,2022] # if don't want to replace GPP_dt with GPP_nt

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
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)


df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hb2=df

#df['GPP_DT']=df['GPP_DT'].fillna(df['GPP_nt'])
#df['GPP_DT'].isnull().any()

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hb2=df

#################################################################################################
# seasonal
df=df_hb2
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hb2'
df['salin'] = 'med'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_hb2_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hb2_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_hb2
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hb2'
df['salin'] = 'med'
df_hb2_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hb2_y.csv', index=False, header=True)

########################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\south carolina")
df=pd.read_csv("fill_US-HB3.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


# 2017,2020 few weeks of GPP missing 

#years_to_remove = [2017,2022] # if don't want to replace GPP_dt with GPP_nt

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
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)



df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hb3=df

#df['GPP_DT']=df['GPP_DT'].fillna(df['GPP_nt'])
#df['GPP_DT'].isnull().any()

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hb3=df

#################################################################################################
# seasonal
df=df_hb3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hb3'
df['salin'] = 'med'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_hb3_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hb3_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_hb3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hb3'
df['salin'] = 'med'
df_hb3_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hb3_y.csv', index=False, header=True)

#####################################################################################################

### new jersey 


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\new jersey")
df=pd.read_csv("fill_US-Hpy.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


# 2017,2020 few weeks of GPP missing 

years_to_remove = [2017] # NEE and LE is missing for aug in 2017 

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
df.to_csv('check.csv', index=False, header=True)



plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)



df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hpy=df

#df['GPP_DT']=df['GPP_DT'].fillna(df['GPP_nt'])
#df['GPP_DT'].isnull().any()

#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_hpy=df

#################################################################################################
# seasonal
df=df_hpy
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hpy'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_hpy_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hpy_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_hpy
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-hpy'
df['salin'] = 'low'
df_hpy_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_hpy_y.csv', index=False, header=True)
#######################################################################################################3

##########################################################################################################

### louisiana 

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")
df=pd.read_csv("fill_US-LA1.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])

# 2017,2020 few weeks of GPP missing 

#years_to_remove = [2017,2022] # if don't want to replace GPP_dt with GPP_nt

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
df.to_csv('check.csv', index=False, header=True)


plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)



df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_la1=df
#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_la1=df

df.to_csv('check.csv', index=False, header=True)
#################################################################################################
# seasonal
df=df_la1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-la1'
df['salin'] = 'low'
#df=df.sort_values(['ET','year','month'], ascending=[False, False, False])
#df=df.sort_values(['GPP','year','month'], ascending=[False, False, False])

df_la1_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_la1_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_la1
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-la1'
df['salin'] = 'low'
df_la1_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_la1_y.csv', index=False, header=True)

##########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")
df=pd.read_csv("fill_US-LA2.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


years_to_remove = [2014,2015,2016,2017,2018,2019,2020] # for monthly 
#years_to_remove = [2014,2015,2016,2017,2018,2019,2020,2021,2022] # for yearly


# 2021 may is missing
# 2022 2 week of may missing 

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
df.to_csv('check.csv', index=False, header=True)


plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)



df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_la2=df
#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_la2=df

df.to_csv('check.csv', index=False, header=True)
#################################################################################################
# seasonal
df=df_la2
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-la2'
df['salin'] = 'low'

df_la2_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_la2_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_la2
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-la2'
df['salin'] = 'low'
df_la2_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_la2_y.csv', index=False, header=True)
####################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana")
df=pd.read_csv("fill_US-LA3.csv")  
df["TIMESTAMP"]=df["DateTime"]
df['TIMESTAMP'] = pd.to_datetime(df["TIMESTAMP"])


#years_to_remove = [2019] # for monthly 
#years_to_remove = [2019,2020,2021] # for yearly


# 2021 may is missing
# 2022 2 week of may missing 

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
df.to_csv('check.csv', index=False, header=True)

plt.plot(df.Tair_f)
plt.plot(df.VPD_f) 
plt.plot(df.Rg_f)



df['ET']=(df['LE_f']/df['lambda'])*(1/1e6)*1800
plt.plot(df['ET'],color='red')
df['ET'].isnull().sum()
df.index = np.arange(0, len(df))
df["NEE"]=df["NEE_f"]*((12/10**6)*1800)
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_la3=df
#################################################################################################
df["GPP"]=df["GPP_DT"]*((12/10**6)*1800)
df["reco"]=df["Reco_DT"]*((12/10**6)*1800)
df["NPP"]=df["GPP"]-df["reco"]

df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])

df=df.groupby([df['TIMESTAMP'].dt.date]).sum()
plt.plot(df['ET'],color='green')
df=df.reset_index()
df_la3=df

df.to_csv('check.csv', index=False, header=True)
#################################################################################################
# seasonal
df=df_la3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year','month']).sum()
df=df.reset_index()

df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-la3'
df['salin'] = 'imed'

df_la3_m=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_la3_m.csv', index=False, header=True)
# monthly -13 to -55
# highest 
######################################################################################################3
## yearly
df=df_la3
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df['year'] = df['TIMESTAMP'].dt.year
df['month'] = df['TIMESTAMP'].dt.month
df= df.groupby(['year']).sum()
df=df.reset_index()
df["WUE"]=df["GPP"]/df["ET"]
df["CUE"]=df["NPP"]/df["GPP"]

plt.plot(df['WUE'],color='green')
plt.plot(df['CUE'],color='green')

df['ID'] = 'US-la3'
df['salin'] = 'imed'
df_la3_y=df

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\results_all")
df.to_csv('df_la3_y.csv', index=False, header=True)
##########################################################################################################

## Alaska Tair data from ERA


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\alaska")
# ERA data is in utc 
import xarray as xr
import pandas as pd
from zoneinfo import ZoneInfo
import pandas as pd


# Path to your NetCDF file
file_path = 'A03_tair.nc'

from netCDF4 import Dataset

# Open the NetCDF file
dataset = Dataset(file_path, mode='r')

# Print the variable names
print("Variables in the NetCDF file:")
print(dataset.variables.keys())

# Open the NetCDF file
dataset = xr.open_dataset(file_path)

# Print the dataset to inspect available variables
print(dataset.variables)

# Select specific variables (e.g., 'temperature', 'humidity')
selected_data = dataset[['valid_time', 'latitude','longitude','expver', 't2m']]

# Convert to a pandas DataFrame
df = selected_data.to_dataframe().reset_index()

# Convert the 'utc_time' column to datetime
df['utc_time'] = pd.to_datetime(df['valid_time'])

# Convert UTC to CST
df['cst_time'] = df['utc_time'].dt.tz_localize('UTC').dt.tz_convert('America/Anchorage')
df['cst_time']=pd.to_datetime(df['cst_time'])
df.drop(['valid_time','utc_time'], axis=1, inplace=True)

df.set_index('cst_time', inplace=True)

df_resampled = df.resample('30T').asfreq()

# Resample to half-hourly frequency and interpolate
df = df.resample('30T').interpolate(method='linear')

df.reset_index(inplace=True)

df['cst_time'] = df['cst_time'].dt.tz_localize(None)

# Define the start and end dates for the entire range
start_date = "2014-01-01 00:00:00"
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


#convert jm-2 to Wm-2

df['Tair']=df.t2m -273.15

plt.plot(df.Tair)

ta_A03=df

ta_A03.to_csv('ta_A03.csv', index=False, header=True)
























##########################################################################################################
### read netcdf file

#################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana\Rg_data")
# ERA data is in utc 
import xarray as xr
import pandas as pd
from zoneinfo import ZoneInfo
import pandas as pd


# Path to your NetCDF file
file_path = 'LA_1data_stream-oper_stepType-accum.nc'

from netCDF4 import Dataset

# Open the NetCDF file
dataset = Dataset(file_path, mode='r')

# Print the variable names
print("Variables in the NetCDF file:")
print(dataset.variables.keys())

# Open the NetCDF file
dataset = xr.open_dataset(file_path)

# Print the dataset to inspect available variables
print(dataset.variables)

# Select specific variables (e.g., 'temperature', 'humidity')
selected_data = dataset[['valid_time', 'latitude','longitude','expver', 'ssrd']]

# Convert to a pandas DataFrame
df = selected_data.to_dataframe().reset_index()

# Convert the 'utc_time' column to datetime
df['utc_time'] = pd.to_datetime(df['valid_time'])

# Convert UTC to CST
df['cst_time'] = df['utc_time'].dt.tz_localize('UTC').dt.tz_convert('America/Chicago')
df['cst_time']=pd.to_datetime(df['cst_time'])
df.drop(['valid_time','utc_time'], axis=1, inplace=True)

df.set_index('cst_time', inplace=True)

df_resampled = df.resample('30T').asfreq()

# Resample to half-hourly frequency and interpolate
df = df.resample('30T').interpolate(method='linear')

df.reset_index(inplace=True)

df['cst_time'] = df['cst_time'].dt.tz_localize(None)

######################################################################################################

# Define the start and end dates for the entire range
start_date = "2012-01-01 00:00:00"
end_date = "2022-12-31 23:59:59"

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


#convert jm-2 to Wm-2
df['Rg']=df.ssrd/3600

plt.plot(df.Rg)

Rg_LA1=df

Rg_LA1.to_csv('Rg_LA1.csv', index=False, header=True)

##########################################################################################################


os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana\Rg_data")
# ERA data is in utc 
# Path to your NetCDF file
file_path = 'LA_2data_stream-oper_stepType-accum.nc'

# Open the NetCDF file
dataset = Dataset(file_path, mode='r')

# Print the variable names
print("Variables in the NetCDF file:")
print(dataset.variables.keys())

# Open the NetCDF file
dataset = xr.open_dataset(file_path)

# Print the dataset to inspect available variables
print(dataset.variables)

# Select specific variables (e.g., 'temperature', 'humidity')
selected_data = dataset[['valid_time', 'latitude','longitude','expver', 'ssrd']]

# Convert to a pandas DataFrame
df = selected_data.to_dataframe().reset_index()

# Convert the 'utc_time' column to datetime
df['utc_time'] = pd.to_datetime(df['valid_time'])

# Convert UTC to CST
df['cst_time'] = df['utc_time'].dt.tz_localize('UTC').dt.tz_convert('America/Chicago')
df['cst_time']=pd.to_datetime(df['cst_time'])
df.drop(['valid_time','utc_time'], axis=1, inplace=True)

df.set_index('cst_time', inplace=True)

df_resampled = df.resample('30T').asfreq()

# Resample to half-hourly frequency and interpolate
df = df.resample('30T').interpolate(method='linear')

df.reset_index(inplace=True)

df['cst_time'] = df['cst_time'].dt.tz_localize(None)
# Reset index if needed


######################################################################################################

# Define the start and end dates for the entire range
start_date = "2012-01-01 00:00:00"
end_date = "2022-12-31 23:59:59"

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



##########################################################################################################
#convert jm-2 to Wm-2
df['Rg_']=df.ssrd/3600

plt.plot(df.Rg_)

Rg_LA2=df

Rg_LA2.to_csv('Rg_LA2.csv', index=False, header=True)

##################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\lousiana\Rg_data")
# ERA data is in utc 
# Path to your NetCDF file
file_path = 'LA_3data_stream-oper_stepType-accum.nc'

# Open the NetCDF file
dataset = Dataset(file_path, mode='r')

# Print the variable names
print("Variables in the NetCDF file:")
print(dataset.variables.keys())

# Open the NetCDF file
dataset = xr.open_dataset(file_path)

# Print the dataset to inspect available variables
print(dataset.variables)

# Select specific variables (e.g., 'temperature', 'humidity')
selected_data = dataset[['valid_time', 'latitude','longitude','expver', 'ssrd']]

# Convert to a pandas DataFrame
df = selected_data.to_dataframe().reset_index()

# Convert the 'utc_time' column to datetime
df['utc_time'] = pd.to_datetime(df['valid_time'])

# Convert UTC to CST
df['cst_time'] = df['utc_time'].dt.tz_localize('UTC').dt.tz_convert('America/Chicago')
df['cst_time']=pd.to_datetime(df['cst_time'])
df.drop(['valid_time','utc_time'], axis=1, inplace=True)

df.set_index('cst_time', inplace=True)

df_resampled = df.resample('30T').asfreq()

# Resample to half-hourly frequency and interpolate
df = df.resample('30T').interpolate(method='linear')

df.reset_index(inplace=True)

df['cst_time'] = df['cst_time'].dt.tz_localize(None)
# Reset index if needed



######################################################################################################

# Define the start and end dates for the entire range
start_date = "2019-01-01 00:00:00"
end_date = "2022-12-31 23:59:59"

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


##########################################################################################################
#convert jm-2 to Wm-2
df['Rg_']=df.ssrd/3600

plt.plot(df.Rg_)

Rg_LA3=df

Rg_LA3.to_csv('Rg_LA3.csv', index=False, header=True)



######################################################################################################

os.chdir(r"C:\ammara_MD\a_yale\Research Project\data\florida")
# ERA data is in utc 
import xarray as xr
import pandas as pd
from zoneinfo import ZoneInfo
import pandas as pd
from netCDF4 import Dataset

# Path to your NetCDF file
#file_path = 'data_stream-oper_stepType-instant.nc'
file_path = 'data_stream-oper_stepType-accum1.nc'



# Open the NetCDF file
dataset = Dataset(file_path, mode='r')

# Print the variable names
print("Variables in the NetCDF file:")
print(dataset.variables.keys())

# Open the NetCDF file
dataset = xr.open_dataset(file_path)

# Print the dataset to inspect available variables
print(dataset.variables)

# Select specific variables (e.g., 'temperature', 'humidity')
selected_data = dataset[['valid_time', 'latitude','longitude','expver', 't2m']]

# Convert to a pandas DataFrame
df = selected_data.to_dataframe().reset_index()

# Convert the 'utc_time' column to datetime
df['utc_time'] = pd.to_datetime(df['valid_time'])

# Convert UTC to CST
df['cst_time'] = df['utc_time'].dt.tz_localize('UTC').dt.tz_convert('America/Anchorage')
df['cst_time']=pd.to_datetime(df['cst_time'])
df.drop(['valid_time','utc_time'], axis=1, inplace=True)

df.set_index('cst_time', inplace=True)

df_resampled = df.resample('30T').asfreq()

# Resample to half-hourly frequency and interpolate
df = df.resample('30T').interpolate(method='linear')

df.reset_index(inplace=True)

df['cst_time'] = df['cst_time'].dt.tz_localize(None)

# Define the start and end dates for the entire range
start_date = "2017-01-01 00:00:00"
end_date = "2020-12-31 23:59:59"

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


#convert jm-2 to Wm-2

df['Tair']=df.t2m -273.15

plt.plot(df.Tair)

ta_A03=df


































#######################################################################################################
## compare salinity
### yearly

df = pd.read_excel('analysis.xlsx', sheet_name='yearly')



plt.rcParams['figure.figsize'] = (16, 7)
fig, ax = plt.subplots()
ax.set_ylabel('Total WUE (g C per kg $H_{2}$O)', color='black')
N = len(df["year"])
ind = np.arange(N)  # the x locations for the groups
width = 0.05       # the width of the bars

#### high
val1 = df["WUE_edn"]
rec1 = ax.bar(ind, val1, width, color='red')
val2 = df["WUE_tas"]
rec2 = ax.bar(ind+width, val2, width, color='red')
val3 = df["WUE_ekp"]
rec3 = ax.bar(ind+width+0.2, val3, width, color='red')
### upper medium

val4 = df["WUE_hb1"]
rec4 = ax.bar(ind+width+0.3, val4, width, color='blue')
val5 = df["WUE_hb2"]
rec5 = ax.bar(ind+width+0.4, val5, width, color='blue')
val6 = df["WUE_hb3"]
rec6 = ax.bar(ind+width+0.5, val6, width, color='blue')

### lower medium
val7 = df["WUE_ks4"]
rec7 = ax.bar(ind+width+0.6, val7, width, color='cornflowerblue')
val8 = df["WUE_ks3"]
rec8 = ax.bar(ind+width+0.7, val8, width, color='cornflowerblue')
val9 = df["WUE_phm"]
rec9 = ax.bar(ind+width+0.8, val9, width, color='cornflowerblue')
val10 = df["WUE_skr"]
rec10 = ax.bar(ind+width+0.9, val10, width, color='cornflowerblue')
val11 = df["WUE_la3"]
rec11 = ax.bar(ind+width+1, val11, width, color='cornflowerblue')

####low 
val12 = df["WUE_evm"]
rec12 = ax.bar(ind+width+1.1, val12, width, color='green')
val13 = df["WUE_srr"]
rec13 = ax.bar(ind+width+1.2, val13, width, color='green')
val14 = df["WUE_myb"]
rec14 = ax.bar(ind+width+1.3, val14, width, color='green')
val15 = df["WUE_la1"]
rec15 = ax.bar(ind+width+1.4, val15, width, color='green')


###very low
val16 = df["WUE_hpy"]
rec16 = ax.bar(ind+width+1.5, val16, width, color='lime')
val17 = df["WUE_dmg"]
rec17 = ax.bar(ind+width+1.6, val17, width, color='lime')
val18 = df["WUE_la2"]
rec18 = ax.bar(ind+width+1.7, val18, width, color='lime')
val19 = df["WUE_tw1"]
rec19 = ax.bar(ind+width+1.8, val19, width, color='lime')

ax.set_xticks(ind+width)
ax.set_xticklabels( ('2007','2008','2009','2010','2011','2012','2013','2014','2015',
                     '2016','2017','2018','2019','2020','2021','2022','2023'))
#ax.legend( (rects1[0], rects2[0],rects3[0]), ('Pines','Potatoes','Precip+Irrigation') ,loc='upper left', fontsize = 'large')
ax.legend( (rec1[0], rec2[0],rec3[0],rec4[0],rec5[0],rec6[0],rec7[0],rec8[0],rec9[0],rec10[0],
rec11[0] ,rec12[0],rec13[0],rec14[0],rec15[0],rec16[0],rec17[0],rec18[0],rec19[0]), ('US-EDN,H','US-TAS:H','US-EKP:H','US-HB1:HM','US-HB2:HM',
 'US-HB3:HM', 'US-KS4:LM', 'US-KS3:LM', 'US-PhM:LM', 'US-SKR:LM','US-LA3:LM','US-EVM:L', 'US-Myb:L', 
 'US-LA1:L', 'US-Hpy:LL','US-Dmg:LL','US-LA2:LL','US-Tw1:LL') ,ncol=5,loc='upper left', fontsize = 'large')
ax.set_ylim(0,5)

ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18) #there is no label 
ax.tick_params(axis = 'y', which = 'major', labelsize = 18)
ax.tick_params(axis = 'x', which = 'major', labelsize = 14)
ax.set_title("Coastal sites: May-Aug (total)",fontsize=15)
fig.autofmt_xdate() 
plt.show()
fig.savefig('testfig1.png',dpi=300, bbox_inches = "tight")
plt.tight_layout()














#########################################################################################################

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

o