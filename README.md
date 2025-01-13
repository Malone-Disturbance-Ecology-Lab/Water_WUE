# Water_WUE

**gaps_US-site.py** Python code that reads in AmeriFlux data from coastal sites in the United States works through the data preprocessing step, prepares a proper time stamp, adds columns such as day of the year (DoY), hour, and year, calculates VPD from RH if VPD is not already calculated. Replace -9999 with nans to be recognized as missing value in the Reddyproc package, and create a dataset with columns such as 'DateTime,' 'Year,' 'DoY,' 'Hour,' 'NEE,' 'LE,' 'H,' 'Rg.' 

**fill_reddyproc.R** This code fills Ameriflux preprocessed data using the Ustar filtering method (reddyproc) for fluxes, calculates night time and daytime "gross primary productivity (GPP)," and respiration using "night time fluxes (Reichstein 2005) and "light response curve method (Lasslop, light response curve). 

**WUE_CUE.py** This code calculates ET from LE and calculates "net primary productivity (NPP)," "Carbon use efficiency (CUE)," and "water use efficiency (WUE) on monthly and yearly time steps."

