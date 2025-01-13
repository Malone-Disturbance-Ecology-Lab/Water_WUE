# Water_WUE

 The purpose of my project is to understand, how much WUE, CUE changes within ecosystems. The specific objective of my research includes, 1) whether differences in WUE, CUE across ecosystems resemble ecosystem variability and 2) determine if WUE, CUE in coastal ecosystems denotes salinity stress. As part of this project, I downloaded data from Ameri flux website.

 There are three steps to execute the analysis for this project:

 1) Preprocess Ameriflux raw data from coastal sites in the United States and prepare a dataframe to be used for the reddyproc package for further processing
 2) Gap fill pre proceeded ameriflux dataset using reddyproc package
 3) Calculate GPP and Reco using the light response curve method
 4) Calulate carbon use efficiency and water use efficiency 
    

**gaps_US-site.py** Python code that reads in AmeriFlux data from coastal sites in the United States works through the data preprocessing step, prepares a proper time stamp, adds columns such as day of the year (DoY), hour, and year, calculates VPD from RH if VPD is not already calculated. Replace -9999 with nans to be recognized as missing value in the Reddyproc package, and create a dataset with columns such as 'DateTime,' 'Year,' 'DoY,' 'Hour,' 'NEE,' 'LE,' 'H,' 'Rg.' 

**fill_reddyproc.R** This code fills Ameriflux preprocessed data using the Ustar filtering method (reddyproc) for fluxes, calculates night time and daytime "gross primary productivity (GPP)," and respiration using "night time fluxes (Reichstein 2005) and "light response curve method (Lasslop, light response curve). 

**WUE_CUE.py** This code calculates ET from LE and calculates "net primary productivity (NPP)," "Carbon use efficiency (CUE)," and "water use efficiency (WUE) on monthly and yearly time steps."

