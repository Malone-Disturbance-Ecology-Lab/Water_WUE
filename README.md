# Water_WUE

 The purpose of my project is to understand, how much WUE, CUE changes within ecosystems. The specific objective of my research includes, 1) whether differences in WUE, CUE across ecosystems resemble ecosystem variability and 2) determine if WUE, CUE in coastal ecosystems denotes salinity stress. As part of this project, I downloaded data from Ameri flux website.


 There are three steps to execute the analysis for this project:

 1) Preprocess Ameriflux raw data from coastal sites in the United States and prepare a data frame to be used for the reddyproc package for further processing
 2) Gap-fill pre proceeded ameriflux dataset using reddyproc package
 3) Calculate GPP and Reco using the light response curve method
 4) Calculate carbon use efficiency and water use efficiency 


    

**1-ameri_api**  This script uses the amerifluxr R package to programmatically download BASE-BADM metadata files for selected AmeriFlux sites. The amf_download_base() function retrieves data for sites related to water use efficiency (WUE) and salinity impact studies. Metadata is saved locally to a shared network directory. All downloads comply with CCBY4.0 licensing via user agreement.

Function Used: amf_download_base()
Output: BADM .csv files saved to ameri_data directory' 




**fill_reddyproc.R** This code fills Ameriflux preprocessed data using the Ustar filtering method (reddyproc) for fluxes, calculates night time and daytime "gross primary productivity (GPP)," and respiration using "night time fluxes (Reichstein 2005) and "light response curve method (light response curve,Lasslop 2010). 




**WUE_CUE.py** This code calculates ET from LE and calculates "net primary productivity (NPP)," "Carbon use efficiency (CUE)," and "water use efficiency (WUE) on monthly and yearly time steps."

