
#rm(list=ls(all=TRUE)) 

phenofit_growing_season <- function(file_path, save_folder) {
  library(tidyverse)
  library(phenofit)
  library(data.table)
  library(lubridate)
  
  average_duplicate_years <- function(df) {
    df %>%
      mutate(year = as.integer(format(origin, "%Y"))) %>%
      group_by(year) %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  }
  
  df <- read.csv(file_path)
  site_list <- unique(str_trim(df$site_name))
  all_sites_results <- list()
  
  for (site in site_list) {
    cat("Processing site:", site, "\n")
    
    sub_df <- df %>%
      filter(str_trim(site_name) == site) %>%
      mutate(Date = as.Date(DateTime)) %>%
      group_by(Date) %>%
      summarise(GPP = mean(GPP, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(w = 1.0)
    
    if (nrow(sub_df) < 50) next
    
    nptperyear <- 365
    wFUN <- wTSM
    wmin <- 0.2
    minExtendMonth <- 0.5
    maxExtendMonth <- 2
    minPercValid <- 0
    south_hemisphere <- FALSE
    
    INPUT <- check_input(
      t = sub_df$Date,
      y = sub_df$GPP,
      w = sub_df$w,
      nptperyear = nptperyear,
      maxgap = nptperyear / 4,
      wmin = wmin,
      mask_spike = TRUE,
      south = south_hemisphere
    )
    
    # Adjust lambda based on site
    if (site %in% c("US-Atq", "US-EvM", "US-StS")) {
      lambda_value <- ifelse(site == "US-Atq", 100, 1000)
      brks <- season_mov(INPUT,
                         list(FUN = "smooth_wWHIT", wFUN = wFUN,
                              maxExtendMonth = 3,
                              wmin = 0.2,
                              r_min = 0.02,
                              ypeak_min = 0.02,
                              lambda = lambda_value,
                              iters = 5,
                              MaxPeaksPerYear = 1,
                              MaxTroughsPerYear = 2))
    } else {
      brks <- season_mov(INPUT,
                         list(FUN = "smooth_wWHIT", wFUN = wFUN,
                              maxExtendMonth = 2,
                              wmin = wmin,
                              r_min = 0.05,
                              ypeak_min = 0.05,
                              lambda = 50000,
                              iters = 5,
                              MaxPeaksPerYear = 1,
                              MaxTroughsPerYear = 2))
    }
    
    fit <- curvefits(INPUT, brks,
                     list(methods = c("Beck", "Elmore", "Gu", "AG", "Zhang"),
                          wFUN = wFUN, iters = 5,
                          wmin = wmin,
                          maxExtendMonth = maxExtendMonth,
                          minExtendMonth = minExtendMonth,
                          minPercValid = minPercValid))
    
    TRS <- c(0.1, 0.25, 0.5)
    l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = FALSE)
    
    doy_Beck   <- average_duplicate_years(l_pheno$doy$Beck)
    doy_Elmore <- average_duplicate_years(l_pheno$doy$Elmore)
    doy_Gu     <- average_duplicate_years(l_pheno$doy$Gu)
    doy_AG     <- average_duplicate_years(l_pheno$doy$AG)
    doy_Zhang  <- average_duplicate_years(l_pheno$doy$Zhang)
    
    common_years <- Reduce(intersect, list(
      doy_Beck$year,
      doy_Elmore$year,
      doy_Gu$year,
      doy_AG$year,
      doy_Zhang$year
    ))
    
    if (length(common_years) == 0) next
    
    doy_Beck   <- filter(doy_Beck, year %in% common_years)
    doy_Elmore <- filter(doy_Elmore, year %in% common_years)
    doy_Gu     <- filter(doy_Gu, year %in% common_years)
    doy_AG     <- filter(doy_AG, year %in% common_years)
    doy_Zhang  <- filter(doy_Zhang, year %in% common_years)
    
    nyears <- length(common_years)
    
    data_grow_season <- tibble(
      site_ID = rep(site, nyears),
      year = common_years,
      avg_DER_sos = rowMeans(cbind(doy_Beck$DER.sos, doy_Elmore$DER.sos, doy_Gu$DER.sos, doy_AG$DER.sos, doy_Zhang$DER.sos), na.rm = TRUE),
      avg_DER_eos = rowMeans(cbind(doy_Beck$DER.eos, doy_Elmore$DER.eos, doy_Gu$DER.eos, doy_AG$DER.eos, doy_Zhang$DER.eos), na.rm = TRUE),
      avg_TRS_sos = rowMeans(cbind(doy_Beck$TRS2.5.sos, doy_Elmore$TRS2.5.sos, doy_Gu$TRS2.5.sos, doy_AG$TRS2.5.sos, doy_Zhang$TRS2.5.sos), na.rm = TRUE),
      avg_TRS_eos = rowMeans(cbind(doy_Beck$TRS2.5.eos, doy_Elmore$TRS2.5.eos, doy_Gu$TRS2.5.eos, doy_AG$TRS2.5.eos, doy_Zhang$TRS2.5.eos), na.rm = TRUE),
      avg_length = rowMeans(cbind(
        doy_Beck$DER.eos - doy_Beck$DER.sos,
        doy_Elmore$DER.eos - doy_Elmore$DER.sos,
        doy_Gu$DER.eos - doy_Gu$DER.sos,
        doy_AG$DER.eos - doy_AG$DER.sos,
        doy_Zhang$DER.eos - doy_Zhang$DER.sos
      ), na.rm = TRUE)
    ) %>%
      mutate(
        avg_sos = (avg_DER_sos + avg_TRS_sos) / 2,
        avg_eos = (avg_DER_eos + avg_TRS_eos) / 2
      ) %>%
      rename(site_name = site_ID)
    
    all_sites_results[[site]] <- data_grow_season
  }
  
  final_df <- bind_rows(all_sites_results)
  out_path <- file.path(save_folder, "phenofit_growing_season.csv")
  write.csv(final_df, out_path, row.names = FALSE)
  cat("All done. Results saved to:", out_path, "\n")
}



##########################################################################
# Define the path to the CSV file
file_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/data_products/all_merged_data.csv"
save_folder <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/data_products"
phenofit_growing_season(file_path, save_folder)
##############################################################################















################################################################################
## how to use this function step by step on one site at a time if debugging 

# Load necessary libraries
library(librarian)
shelf(tidyverse, phenofit, data.table, lubridate, cowplot, ggpubr)
library(stringr)
library(dplyr)

# try this for site EvM
# 100 for Atq
#1000 EvM

# use this for US-EvM
brks <- season_mov(INPUT,
                   list(FUN = "smooth_wWHIT",
                        wFUN = wFUN,
                        maxExtendMonth = 3,
                        wmin = 0.2,            
                        r_min = 0.02,         
                        ypeak_min = 0.02, # make this parameter smaller from 0.05 to 0.02 if cant find growing season
                        lambda = 1000, # make it smaller too if cant find growing       
                        iters = 5,
                        MaxPeaksPerYear = 1,
                        MaxTroughsPerYear = 2))



# Define a function to extract year and average duplicate rows
average_duplicate_years <- function(df) {
  df %>%
    mutate(year = as.integer(format(origin, "%Y"))) %>%
    group_by(year) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}


df <- read.csv(file_path)
full_df=df

unique(df$site_name)
#unique(df$site_name)
#site_to_test <- "US-A03"
#site_to_test <- "US-A10"
#site_to_test <- "US-Atq"

#site_to_test <- "US-EKH"
#site_to_test <- "US-HB1"
#site_to_test <- "US-LA2"
#site_to_test <- "US-Skr" # 
#site_to_test <- "US-MRM"
#site_to_test <- "US-Myb"

#site_to_test <- "US-Tw1"  
#site_to_test <- "US-Dmg"  
#site_to_test <- "US-TaS" 
#site_to_test <- "US-EDN" 
#site_to_test <- "US-EKP" 
#site_to_test <- "US-EvM"
#site_to_test <- "US-HB1"
#site_to_test <- "US-HB2"
#site_to_test <- "US-HB3"
#site_to_test <-"US-HPY"
#site_to_test <-"US-KS3"
#site_to_test <-"US-KS4"
#site_to_test <-"US-LA1"
#site_to_test <-"US-LA2"
#site_to_test <-"US-NGB"
#site_to_test <-"US-PhM"
#site_to_test <-"US-Srr"
#site_to_test <-"US-StJ"
#################################################################
## start from here

site_to_test <- "US-Atq" 
site_to_test <- "US-StS" 

site_ID <- site_to_test
###############################################################################
##########################################################################

#Step 1: Aggregate Half-Hourly Data to Daily GPP
sub_df <- df %>%
  filter(str_trim(site_name) == site_to_test) %>%
  mutate(Date = as.Date(DateTime)) %>%
  group_by(Date) %>%
  summarise(GPP = mean(GPP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(w = 1.0)  # Assign uniform weights since data is quality controlled

unique(format(sub_df$Date, "%Y"))
###############################################################################

#Parameters
nptperyear <- 365
wFUN <- wTSM
wmin <- 0.2
methods_fine <- c("AG", "Zhang", "Beck", "Elmore", "Gu")
minExtendMonth <- 0.5
maxExtendMonth <- 2
minPercValid <- 0
south_hemisphere <- FALSE  # Set to TRUE if the site is in the Southern Hemisphere



#Step 2: Prepare Input for Phenofit

INPUT <- check_input(
  t = sub_df$Date,
  y = sub_df$GPP,
  w = sub_df$w,
  QC_flag = NULL,  # QC_flag is not available
  nptperyear = nptperyear,
  maxgap = nptperyear / 4,
  wmin = wmin,
  mask_spike = TRUE,
  south = south_hemisphere
)

## use this for all sites excep US-Atq
brks <- season_mov(INPUT,
                   list(FUN = "smooth_wWHIT", wFUN = wFUN,
                        #        minpeakdistance = 100,
                        maxExtendMonth = 2,
                        wmin = wmin, r_min = 0.05, ypeak_min = 0.05, 
                        lambda = 50000,
                        iters = 5,
                        MaxPeaksPerYear = 1,
                        MaxTroughsPerYear = 2
                   ))


## lambda values

# 50000 # Dmg, A03 
# 10000 # TW1



#Step 4: Curve Fitting Using Multiple Methods
fit <- curvefits(
  INPUT,
  brks,
  list(
    methods = methods_fine,
    wFUN = wFUN,
    iters = 5,
    wmin = wmin,
    maxExtendMonth = maxExtendMonth,
    minExtendMonth = minExtendMonth,
    minPercValid = minPercValid
  )
)

#Step 5: Extract Phenological Metrics

TRS <- c(0.1, 0.25, 0.5)
l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = FALSE)
plot_season(INPUT, brks)   # Before curve fitting




## fix problem here x
# Extract data frames for each method
doy_Beck <- l_pheno$doy$Beck
doy_Elmore <- l_pheno$doy$Elmore
doy_Gu <- l_pheno$doy$Gu
doy_AG <- l_pheno$doy$AG
doy_Zhang <- l_pheno$doy$Zhang


# Apply aggregate function to merge years when two sos and eod dates each dataset
doy_Beck   <- average_duplicate_years(doy_Beck)
doy_Elmore <- average_duplicate_years(doy_Elmore)
doy_Gu     <- average_duplicate_years(doy_Gu)
doy_AG     <- average_duplicate_years(doy_AG)
doy_Zhang  <- average_duplicate_years(doy_Zhang)


# Step 2: Determine common years across methods
common_years <- Reduce(intersect, list(
  doy_Beck$year,
  doy_Elmore$year,
  doy_Gu$year,
  doy_AG$year,
  doy_Zhang$year
))

# Step 3: Filter data frames to include only common years
doy_Beck <- filter(doy_Beck, year %in% common_years)
doy_Elmore <- filter(doy_Elmore, year %in% common_years)
doy_Gu <- filter(doy_Gu, year %in% common_years)
doy_AG <- filter(doy_AG, year %in% common_years)
doy_Zhang <- filter(doy_Zhang, year %in% common_years)

# Step 4: Build Averaged Phenology Results
nyears <- length(common_years)
data_grow_season <- tibble(
  site_ID = rep(site_ID, nyears),
  year = common_years,
  avg_DER_sos = rowMeans(cbind(doy_Beck$DER.sos, doy_Elmore$DER.sos, doy_Gu$DER.sos, doy_AG$DER.sos, doy_Zhang$DER.sos), na.rm = TRUE),
  avg_DER_eos = rowMeans(cbind(doy_Beck$DER.eos, doy_Elmore$DER.eos, doy_Gu$DER.eos, doy_AG$DER.eos, doy_Zhang$DER.eos), na.rm = TRUE),
  avg_TRS_sos = rowMeans(cbind(doy_Beck$TRS2.5.sos, doy_Elmore$TRS2.5.sos, doy_Gu$TRS2.5.sos, doy_AG$TRS2.5.sos, doy_Zhang$TRS2.5.sos), na.rm = TRUE),
  avg_TRS_eos = rowMeans(cbind(doy_Beck$TRS2.5.eos, doy_Elmore$TRS2.5.eos, doy_Gu$TRS2.5.eos, doy_AG$TRS2.5.eos, doy_Zhang$TRS2.5.eos), na.rm = TRUE),
  avg_length = rowMeans(cbind(
    doy_Beck$DER.eos - doy_Beck$DER.sos,
    doy_Elmore$DER.eos - doy_Elmore$DER.sos,
    doy_Gu$DER.eos - doy_Gu$DER.sos,
    doy_AG$DER.eos - doy_AG$DER.sos,
    doy_Zhang$DER.eos - doy_Zhang$DER.sos
  ), na.rm = TRUE)
)


data_grow_season <- data_grow_season %>%
  mutate(
    avg_sos = (avg_DER_sos + avg_TRS_sos) / 2,
    avg_eos = (avg_DER_eos + avg_TRS_eos) / 2
  )

data_grow_season <- data_grow_season %>%
  rename(site_name = site_ID)



################################################################################









































## if errors because of nans use this 

## filter data 
sub_df <- df %>%
  filter(str_trim(site_name) == site_to_test) %>%
  mutate(Date = as.Date(DateTime),
         Year = year(Date),
         Month = month(Date))

# Calculate the number of valid NEE_f entries per year for months 2 to 10
valid_years <- sub_df %>%
  filter(Month %in% 5:8) %>%
  group_by(Year) %>%
  summarise(
    total_days = n(),
    valid_days = sum(!is.na(NEE_f))
  ) %>%
  mutate(valid_ratio = valid_days / total_days) %>%
  filter(valid_ratio > 0.5) %>%
  pull(Year)

sub_df_filtered <- sub_df %>%
  filter(Year %in% valid_years)

# Identify numeric columns to aggregate
numeric_vars <- sub_df_filtered %>%
  select(-DateTime, -site_name, -Date, -Year, -Month) %>%
  select(where(is.numeric)) %>%
  names()

# Aggregate to daily means
daily_df <- sub_df_filtered %>%
  group_by(Date) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)),
            .groups = 'drop') %>%
  mutate(w = 1.0)  # Assign uniform weights since data is quality controlled


daily_df$DateTime=daily_df$Date

daily_gpp <- daily_df %>%
  select(Date, GPP) %>%
  rename(t = Date, y = GPP) %>%
  mutate(w = 1.0)  # Assign uniform weights since data is quality controlled


colnames(daily_gpp) <- c("Date", "GPP", "w")
unique(format(daily_df$Date, "%Y"))
#sub_df=daily_gpp
# 4289 for Tw1
## filter data step ends 
###########################################################################








