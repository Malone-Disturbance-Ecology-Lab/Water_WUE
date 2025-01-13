# Water_WUE
AT

**gaps_US-site ** Paython code that reads in AmeriFlux data from coastal sites in United states, works through data preprocessing step, prepares a proper time stamp, add columns such as
day of year (DoY), hour, year, calculate VPD from RH, if VPD is not already calculated. replace -9999 with nans to be recognized as missing value in Reddyproc package. 



estimates light response curve and temperature response curve parameters, models site level carbon exchange based on LAI, radiation and temperature, and then makes figures. This code calls Light_Reponse_Function.m and Temp_Response_Funtion.m.

** Light_Reponse_Function.m ** Estimates light response curve paraments when given radiation and NEE data. This function calls NEP_1PredVar_Model.m.

** Temp_Reponse_Function.m ** Estimates light response curve paraments when given temperature and NEE data.

** NEP_1PredVar_Model.m ** Defines nonrectangular light response function used in Light_Reponse_Function.m.
