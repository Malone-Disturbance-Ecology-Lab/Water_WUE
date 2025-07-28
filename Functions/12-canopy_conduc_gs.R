
rm(list=ls(all=TRUE)) 
library(Rcpp)
library(REddyProc)
library(readr)

# Calculating gs (canopy conductance) from eddy flux data
# Analytical code based on file AmeriFluxHo1gsCalc.ipf from Wehr et al Biogeosci 2021
# on Dryad archive at https://datadryad.org/stash/dataset/doi:10.5061/dryad.h44j0zpgp
# to produce inverted Penman-Monteith (iPM) and Flux Gradient (FG) estimates 
# of canopy conductance (gs). 
# WARNING: still need to filter for not recent rain events after this!
#
# Inputs: 
#   1. flux_dat input dataframe with columns: 
#     datetime, Tair (C), PA (kPa), RH (%), VPD (kPa), NETRAD (W m-2),  
#     LE (W m-2), H (W m-2), PAR (umol m-2 s-1) 
#     optional: lambda (latent heat vap.; J kg-1) or will use constant 2.481
#   2. bounds for month_min and month_max for growing season (numeric)
#   3. bound for LE_min (only keep data where LE is greater than this for gs estimation)
#   4. bound for PAR_min (only keep data where PAR above)
#   5. out output format to return: can be "all" for all values, or "daily" or "monthly" for summarized
#
# Output: options to return full data frame "all" with original flux_dat + gsv added
#   or summarized stats (median, 25/75 %tiles) at "monthly" or "daily" 

growseas_gs = function(flux_dat, month_min, month_max, 
                       LE_min,PAR_min, out){
  
  # Subset flux data to month min/max and LE > LE_min and PAR > PAR_min
  flux_subset = flux_dat 
  flux_subset$datetime= as.Date(ems_ustfilt$datetime,"%m/%d/%Y %H:%M")
  flux_subset$year = as.numeric(format(flux_subset$datetime, "%Y"))
   flux_subset$month = as.numeric(format(flux_subset$datetime, "%m"))
  flux_subset$hour = as.numeric(format(flux_subset$datetime, "%H"))
  flux_subset = flux_subset[flux_subset$month >= month_min & flux_subset$month <= month_max,]
  flux_subset = flux_subset[!is.na(flux_subset$LE) & flux_subset$LE >= LE_min,]
  flux_subset = flux_subset[!is.na(flux_subset$PAR) & flux_subset$PAR >= PAR_min,]

  
  # Define constants
  sidesOfLeafWithStomata = 1 #1: hypostomatous leaves, 2: amphistomatous leaves
  gasConstR = 8.314472 #J/mol/K
  c_p = 1012 #specific heat of air, J/kg/K
  airMolarMass_dry = 0.02897 #molar mass of dry air, kg mol-1
  if(!is.null(flux_subset$lambda)){
    lambda = flux_subset$lambda #latent heat of vaporization of water, J/kg
  } else {
    lambda = 2.45e6 #latent heat of vaporization of water, J/kg (approx)
  }
  Sc_CO2 = 1.05 #1.02 //Schmidt number = 1.02 for CO2 from Ogee et al (2003), 1.0347 based on ozone value of Lamaud and ratio of binary diffusion coefficients, 1.14 from a website, 1.1 by my calculations based on another website
  Sc_H2O = Sc_CO2/1.57278 #1.57278 is ratio of binary diffusion coefficients of H2O and CO2 //Schmidt number = 0.6 for H2O, from Kramm et al (2002)
  Pr = 0.71 #Prandtl number = 0.72 from Ogee et al (2003) but can be between 0.5 and 1 according to Kramm et al 2002
  
  # Finding hourly energy budget gap
  G_guess = flux_subset$NETRAD*15/700 # this is a rough estimate of heat flux to ground, considered roughly propotional to Rnet based on Lindroth et al 2010
  TotalTurbHeatFlux = flux_subset$H + flux_subset$LE
  AvailableEnergy = flux_subset$NETRAD - G_guess
  EnergyImbalance = AvailableEnergy - TotalTurbHeatFlux
  
  # Estimate slope between turbulent flux (H+LE) and available energy (Rnet-G)
  hourlyEbudgratio = lm(TotalTurbHeatFlux ~ AvailableEnergy)
  hourlyEnergyBudgetRatio = coef(summary(hourlyEbudgratio))[2]
  
  # Flux correction to achieve daily budget closure (using annual-average gap)
  TotalTurbHeatFlux_hCorr = TotalTurbHeatFlux/hourlyEnergyBudgetRatio
  H_hCorr = flux_subset$H/hourlyEnergyBudgetRatio
  LE_hCorr = flux_subset$LE/hourlyEnergyBudgetRatio
  
  # Variables to calculate air physics vars & convert H and LE to flux
  P_pa = flux_subset$PA*1000  #Atmospheric pressure [Pa]
  ma = 28.964/1000    #molar mass of dry air [Kg/mol]
  mv = 18/1000        #molar mass of water vapor [Kg/mol]
  R = 8.314          #Universal gas constant dry air [J/(K mol)]
  Cpa_dry = 1004.67  #J Kg-1 K-1 - specific heat of dry air
  
  #Compute wet and dry air concentrations and heat capacities...
  SatVP = 100*6.112*exp(17.62*flux_subset$Tair/(243.12 + flux_subset$Tair)) #Pa
  SatVPslope = 4098*SatVP/(flux_subset$Tair+237.3)^2
  VP = (flux_subset$RH)/100*(SatVP) #1000 is to convert wH2O from mmol/mol to mol/mol
  VPDcheck = SatVP - VP #Pa
  VP = SatVP - flux_subset$VPD*1e3 #Pa
  VP_n = VP #because Re is taken to be zero in this version of the code
  AirConc_dry = (flux_subset$PA*1e3 - VP)/(gasConstR*(flux_subset$Tair+273.15)) #dry air molar concentration, mol/m3
  AirConc_wet = flux_subset$PA*1e3/(gasConstR*(flux_subset$Tair+273.15)) #wet air molar concentration, should be mol/m3
  
  HeatCapacity_dry = 1003 + (1008 - 1003)*((flux_subset$Tair + 23.16)/100) #J/kg/K, linear interpolation of values from http://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
  HeatCapacity_waterVapor = HeatCapacity_dry*(1+0.84*flux_subset$RH/100) # TEST!! Specific heat of moist air [J kg-1 K-1] 
  AirDensity_dry = AirConc_dry*airMolarMass_dry
  WaterVaporDensity = (AirConc_wet - AirConc_dry)*0.018
  AirDensity_wet = AirDensity_dry + WaterVaporDensity
  HeatCapacity_wet = (AirDensity_dry*HeatCapacity_dry + WaterVaporDensity*HeatCapacity_waterVapor)/AirDensity_wet #J/kg/K, heat capacity of wet air
  Gamma = HeatCapacity_wet*flux_subset$PA*1e3/(lambda*0.62198) #~66.1 //psychrometric constant, Pa/K
  
  # Eddy resistance
  R_e = 1e-16 #in absence of CO2 profiles, this line just assumes no turbulent eddy resistance
  G_e = 1/R_e
  R_e = flux_subset$PA/(gasConstR*(flux_subset$Tair+273.15)*G_e)
  
  # Boundary resistance to heat
  R_bH = 10 # s m-1, boundary resistance to heat
  R_b = (2/sidesOfLeafWithStomata)*R_bH*((Sc_CO2/Pr)^(2/3)) #(Sc_CO2/Pr)^(2/3) gives 1.26, not (1.4*0.92), boundary layer resistance (s m-1) for CO2 //1.4 ~ 1.6^(2/3) for quasi-diffusion
  R_bV = (2/sidesOfLeafWithStomata)*R_bH*((Sc_H2O/Pr)^(2/3)) #(Sc_H2O/Pr)^(2/3) gives 0.88, not 0.92 #valid for transpiration only, factor of 2 because leaves are hypostomatic but heat comes from both sides, wR_b/1.4
  
  # Boundary conductance from resistance
  G_bH = flux_subset$PA*1e3/(gasConstR*(flux_subset$Tair+273.15)*R_bH) #but temperature here should be closer to leaf temperature (maybe not significant)?
  G_b = flux_subset$PA*1e3/(gasConstR*(flux_subset$Tair+273.15)*R_b) #but temperature here should be closer to leaf temperature (maybe not significant)?
  G_bV = flux_subset$PA*1e3/(gasConstR*(flux_subset$Tair+273.15)*R_bV) #but temperature here should be closer to leaf temperature (maybe not significant)?
  
  # Water flux, leaf temp & VPD
  E = LE_hCorr/(lambda*0.018) #mol m-2 s-1, water flux
  Tair_n = (H_hCorr*R_e/(AirDensity_wet*HeatCapacity_wet)) + flux_subset$Tair
  Tair_leaf = (H_hCorr*R_bH/(AirDensity_wet*HeatCapacity_wet)) + Tair_n
  SatVP_leaf = 100*6.112*exp(17.62*Tair_leaf/(243.12 + Tair_leaf)) #Pa
  LeafAirVPD = (SatVP_leaf - VP_n)
  LeafAirConcDiff = (1/(gasConstR*(Tair_n+273.15)))*(SatVP_leaf - VP_n)
  
  # Flux gradient method
  LE_hCorr = (lambda*0.018)*E
  
  # Canopy resistance to water vapor
  R_sV = LeafAirConcDiff/E - R_bV #this will be more correct if concentration gradients drive diffusion, rather than partial pressure gradients
  R_s = R_sV*1.57278 #because binary diffusion rate of CO2-N2 is about 1.6 times slower than that of H2O-N2
  
  # Canopy conductance to water vapor
  G_sV_hCorr = flux_subset$PA*1e3/(gasConstR*(Tair_leaf+273.15)*R_sV) #conversion from s/m to mol/m2/s given by Grace et al. (1995)
  G_s = flux_subset$PA*1e3/(gasConstR*(Tair_leaf+273.15)*R_s) #conversion from s/m to mol/m2/s given by Grace et al. (1995)
  
  # Penman-Monteith inversion method
  # PM Parameters
  k = 0.41		# Von Karman constant
  h =10			# Canopy height [m] # 25 for harvard forest
  zm = 29			# Measurement height [m]
  zd = 0.67*h		# Zero plane displacement [m]
  zo = 0.1*h		# Momentum roughness length [m]
  Z_combine = (zm - zd)/zo
  
  # iPM resistance to water vapor
  R_sV_PM = (SatVPslope*(flux_subset$NETRAD-5-flux_subset$LE)*R_bH +
               AirDensity_wet*HeatCapacity_wet*(SatVP - VP))/
    (Gamma*flux_subset$LE) - R_bV
  
  # iPM conductance to water vapor
  G_sV_PM_hCorr = (flux_subset$PA*1e3)/(gasConstR*(flux_subset$Tair+273.15)*R_sV_PM)
  
  # Compile data frame with FG and iPM E-balance-corrected conductance
  flux_subset$gs_iPM = G_sV_PM_hCorr
  flux_subset$gs_FG = G_sV_hCorr
  flux_subset$VPDl = LeafAirVPD
  
  # Calculate median, 25/75 %tiles growing season gs by year
  
  if(out == "all"){
    return(flux_subset)
  } else if(out == "monthly") {
    gs_group = dplyr::group_by(flux_subset, month, year)
    gs_stats = dplyr::summarize(gs_group, 
                                gs_iPM_med = median(gs_iPM, na.rm=T),
                                gs_iPM_25 = quantile(gs_iPM, 0.25, na.rm=T),
                                gs_iPM_75 = quantile(gs_iPM, 0.75, na.rm=T),
                                gs_FG_med = median(gs_FG, na.rm=T),
                                gs_FG_25 = quantile(gs_FG, 0.25, na.rm=T),
                                gs_FG_75 = quantile(gs_FG, 0.75, na.rm=T))
    return(gs_stats)
  } else if (out=="daily"){
    gs_group = dplyr::group_by(flux_subset, date)
    gs_stats = dplyr::summarize(gs_group, 
                                gs_iPM_med = median(gs_iPM, na.rm=T),
                                gs_iPM_25 = quantile(gs_iPM, 0.25, na.rm=T),
                                gs_iPM_75 = quantile(gs_iPM, 0.75, na.rm=T),
                                gs_FG_med = median(gs_FG, na.rm=T),
                                gs_FG_25 = quantile(gs_FG, 0.25, na.rm=T),
                                gs_FG_75 = quantile(gs_FG, 0.75, na.rm=T))
    return(gs_stats)
  } else {
    "No output file format (all, daily, monthly, gs) given."
  }
}


# read in file for ustar filtered, non-gapfilled eddy fluxes
#ems_ustfilt = read_csv("~/Desktop/ems_archive_2023/ems_ustfilt_2023.csv")

#ems_ustfilt <-read_csv(file="C:/ammara_MD/a_harvard/Project/HF_paper/ems_ustfilt_2023.csv")
#ems_ustfilt <-read_csv(file="C:/ammara_MD/a_harvard/Project/HF_paper/ems_ustfilt_2023_filled.csv")

ems_ustfilt <-read_csv(file="C:/ammara_MD/a_harvard/Project/HF_paper/gaps_ML_HEM_2023.csv")



##### CANOPY CONDUCTANCE
# Estimate growing season canopy conductance by year
flux_dat = ems_ustfilt
#flux_dat$lambda = flux_dat$lambda*1e6 # use it if EMS lambda is 2.6 
flux_dat$lambda = flux_dat$lambda  # use this for hemlock

flux_dat$PAR = flux_dat$PPFD

# Calculate hourly whole growing season canopy conductance
ems_gs = growseas_gs(flux_dat, month_min = 1, month_max = 12,
                     LE_min = 50, PAR_min = 500, out="all")


#write.csv(ems_gs,"C:/ammara_MD/a_harvard/Project/HF_paper/ems_gs_all.csv")


write.csv(ems_gs,"C:/ammara_MD/a_harvard/Project/HF_paper/hem_gs_all.csv")



 # Calculate hourly growing season canopy conductance
ems_gs = growseas_gs(flux_dat, month_min = 6, month_max = 8,
                     LE_min = 50, PAR_min = 500, out="daily")


#mol m−2 s−1

