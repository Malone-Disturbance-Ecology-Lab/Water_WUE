rm(list=ls(all=TRUE)) 
library(Rcpp)
library(REddyProc)
#packageVersion("REddyProc") 1.3.3
library(lubridate)
# Get general information for all AmeriFlux sites
library("amerifluxr")
library(lutz)
site_info <- amf_site_info()
#https://rdrr.io/cran/flux/man/gpp.html
#https://rdrr.io/rforge/REddyProc/man/sEddyProc.example.html

## this function will change the parameters if default parameters don't work to calculate ustar

file_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/ameri_data/reddy_gaps/blended_gaps2"
save_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/ameri_data/ameri_fill"



# Loop through all CSV files in the file path
files <- list.files(file_path, pattern = "*.csv", full.names = TRUE)

## return year with problematic year where ustar can not be calculated

estimateUStarThreshold_start <- function(EddyProc.C, EddyDataWithPosix.F) {
  tryCatch({
    # Try estimating UStar threshold
    uStarTh_default <- EddyProc.C$sEstUstarThold()
    message("No problematic year found during uStar threshold estimation.")
    return(NULL)
  }, error = function(e) {
    if (grepl("Expected valid uStar records for year", e$message)) {
      problematic_year <- as.numeric(sub(".*for year (\\d{4}).*", "\\1", e$message))
      message("Problematic year detected: ", problematic_year)
      return(problematic_year)
    } else {
      message("No problematic year found during uStar threshold estimation.")
      return(NULL)
    }
  })
}


## this function calculate default ustar first, if error then it adjust parameters
# if still error tehn use user defined ustar

estimateUStarThreshold <- function(EddyProc.C, EddyDataWithPosix.F) {
  uStarTh_default <- EddyProc.C$sEstUstarThold()
  
  if (is.data.frame(uStarTh_default) && "uStar" %in% names(uStarTh_default)) {
    uStar_values <- uStarTh_default$uStar
  } else {
    stop("Unexpected output format from EddyProc.C$sEstUstarThold()")
  }
  
  if (all(is.na(uStar_values))) {
    message("Estimated UStar threshold is NA. Using custom parameters.")
    
    seasonFactor <- usCreateSeasonFactorMonth(EddyDataWithPosix.F$DateTime, startMonth = 1)
    ctrlUstarSub <- usControlUstarSubsetting(
      taClasses = 7, 
      UstarClasses = 20,
      swThr = 10, 
      minRecordsWithinTemp = 10,
      minRecordsWithinSeason = 100, 
      minRecordsWithinYear = 3000,
      isUsingOneBigSeasonOnFewRecords = TRUE
    )
    
    uStarTh <- EddyProc.C$sEstUstarThold(
      seasonFactor = seasonFactor,  
      ctrlUstarSub = ctrlUstarSub
    )
    
    if (all(is.na(uStarTh$uStar))) {
      message("Fallback also resulted in NA. Setting uStar = 0.1.")
      uStar <- 0.1
      return(uStar)
    }
  } else {
    message("Valid UStar threshold found. Using default estimation.")
    uStarTh <- uStarTh_default
  }
  
  return(uStarTh)
}


reddy_proc <- function(file_path, save_path, files) {
  for (full_file_path in files) {
    tryCatch({
      # Extract file name and SITE_ID
      file <- basename(full_file_path)
      SITE_ID <- sub("^gaps_blend_(.*)\\.csv$", "\\1", file)  # Extract site ID
      message("Processing file: ", file, " | SITE_ID: ", SITE_ID)
      
      # Load data
      data <- read.csv(full_file_path, na.strings = c("NA", "NaN", "", " "), stringsAsFactors = FALSE)
      df <- data[, -1]  # Remove the first column (usually index column)
      
      # Convert time
      EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
      
      # Define and filter required columns
      required_columns <- c('NEE', 'LE', 'H', 'Rg', 'Tair', 'VPD', 'Ustar', 'PA', 'NETRAD', 'WS', 'WD', 'RH', 'PAR')
      available_columns <- required_columns[required_columns %in% colnames(EddyDataWithPosix.F)]
      
      # Ensure PA is numeric
      if (!is.numeric(EddyDataWithPosix.F$PA)) {
        EddyDataWithPosix.F$PA <- as.numeric(EddyDataWithPosix.F$PA)
      }
      
      # Drop columns from available_columns that are all NA
      available_columns <- available_columns[!sapply(EddyDataWithPosix.F[available_columns], function(col) all(is.na(col)))]
      
      # Initialize sEddyProc
      EddyProc.C <- sEddyProc$new('df', EddyDataWithPosix.F, available_columns, DTS = 48)
      
      # Estimate uStar threshold (non-stopping version)
      problematic_info <- estimateUStarThreshold_start(EddyProc.C, EddyDataWithPosix.F)
      
      # Extract unique problematic years from the result (if any)
      if (!is.null(problematic_info)) {
        if (is.numeric(problematic_info)) {
          problematic_years <- unique(na.omit(problematic_info))
        } else if (is.data.frame(problematic_info) && "year" %in% names(problematic_info)) {
          problematic_years <- unique(na.omit(problematic_info$year))
        } else {
          problematic_years <- NULL
        }
        
        if (!is.null(problematic_years) && length(problematic_years) > 0) {
          df <- df[!df$Year %in% problematic_years, ]
          EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
          message("Removing problematic year(s): ", paste(problematic_years, collapse = ", "))
        } else {
          message("No problematic years detected. Proceeding to next step.")
        }
      } else {
        message("No problematic years detected. Proceeding to next step.")
      }
      
      # After removing problematic year, redo the analysis
      EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
      
      # Continue with the rest of the process
      required_columns <- c('NEE', 'LE', 'H', 'Rg', 'Tair', 'VPD', 'Ustar', 'PA', 'NETRAD', 'WS', 'WD', 'RH', 'PAR')
      available_columns <- required_columns[required_columns %in% colnames(EddyDataWithPosix.F)]
      
      # Ensure PA is numeric
      if (!is.numeric(EddyDataWithPosix.F$PA)) {
        EddyDataWithPosix.F$PA <- as.numeric(EddyDataWithPosix.F$PA)
      }
      
      # Store original column info before dropping for filling later
      original_cols_present <- colnames(df)
      
      # Drop columns from available_columns that are all NA
      available_columns <- available_columns[!sapply(EddyDataWithPosix.F[available_columns], function(col) all(is.na(col)))]
      
      # Initialize sEddyProc
      EddyProc.C <- sEddyProc$new('df', EddyDataWithPosix.F, available_columns, DTS = 48)
      
      # Estimate uStar threshold
      uStarThreshold <- estimateUStarThreshold(EddyProc.C, EddyDataWithPosix.F)
      # Perform gap-filling with and without uStar filtering
      
      if (is.numeric(uStarThreshold) && uStarThreshold == 0.1) {
        uStar <- 0.1
        EddyProc.C$sMDSGapFillAfterUstar('NEE', uStarTh = uStar)
        EddyProc.C$sMDSGapFillAfterUstar('LE', uStarTh = uStar)
        EddyProc.C$sMDSGapFillAfterUstar('H', uStarTh = uStar)
        EddyProc.C$sExportResults()
      } else {
        EddyProc.C$sMDSGapFillAfterUstar('NEE')
        EddyProc.C$sMDSGapFillAfterUstar('LE')
        EddyProc.C$sMDSGapFillAfterUstar('H')
        EddyProc.C$sExportResults()
      }
      
      # Fill meteorological variables
      EddyProc.C$sMDSGapFill('Rg')
      EddyProc.C$sMDSGapFill('Tair', FillAll.b = FALSE)
      EddyProc.C$sMDSGapFill('VPD', FillAll.b = FALSE)
      
      # Fill PAR only if it exists AND has at least some non-NA values
      if ("PAR" %in% original_cols_present && "PAR" %in% colnames(EddyDataWithPosix.F) && !all(is.na(EddyDataWithPosix.F$PAR))) {
        EddyProc.C$sMDSGapFill('PAR', FillAll.b = FALSE)
      }
      
      # NETRAD - fill only if it exists AND has at least some non-NA values
      if ("NETRAD" %in% original_cols_present && "NETRAD" %in% colnames(EddyDataWithPosix.F) && !all(is.na(EddyDataWithPosix.F$NETRAD))) {
        EddyProc.C$sMDSGapFill('NETRAD', FillAll.b = FALSE)
      }
      
      filled1 <- EddyProc.C$sExportResults()
      
      # Get site lat/lon and timezone info
      site_info_specific <- site_info[site_info$SITE_ID == SITE_ID, ]
      lat <- site_info_specific$LOCATION_LAT
      lon <- site_info_specific$LOCATION_LONG
      timezone_info <- tz_lookup_coords(lat, lon, method = "accurate")
      standard_time <- as.POSIXct("2025-01-01", tz = timezone_info)
      standard_offset <- as.numeric(format(standard_time, "%z")) / 100
      
      # Set location and perform flux partitioning
      EddyProc.C$sSetLocationInfo(Lat_deg.n = lat, Long_deg.n = lon, TimeZone_h.n = standard_offset)
      EddyProc.C$sMRFluxPartition(Suffix.s = 'uStar')  # night time #Reichstein 2005
      EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar')  # day time  Reco_DT, GPP_DT #Lasslop, light response curve
      filled2 <- EddyProc.C$sExportResults()
      
      # Start with required/base columns
      filled3 <- data.frame(
        DateTime = df$TIMESTAMP,
        Year     = df$Year,
        DoY      = df$DoY,
        Hour     = df$Hour,
        stringsAsFactors = FALSE
      )
      
      # Add original variables without filling (these will have empty cells where data missing, not "NA" text)
      # Ustar - original, no filling
      if ("Ustar" %in% colnames(df)) {
        filled3$Ustar <- df$Ustar
      }
      
      # PA - original, no filling
      if ("PA" %in% colnames(df)) {
        filled3$PA <- df$PA
      }
      
      # RH - original, no filling
      if ("RH" %in% colnames(df)) {
        filled3$RH <- df$RH
      }
      
      # WS - original, no filling
      if ("WS" %in% colnames(df)) {
        filled3$WS <- df$WS
      }
      
      # WD - original, no filling
      if ("WD" %in% colnames(df)) {
        filled3$WD <- df$WD
      }
      
      # Define mapping for filled columns from filled2
      optional_cols <- list(
        NEE_uStar_f    = "NEE_f",
        LE_uStar_f     = "LE_f",
        H_uStar_f      = "H_f",
        Tair_f         = "Tair_f",
        VPD_f          = "VPD_f",
        Rg_f           = "Rg_f",
        GPP_DT_uStar   = "GPP_DT",
        GPP_uStar_f    = "GPP_nt",
        Reco_DT_uStar  = "Reco_DT",
        Reco_uStar     = "Reco_nt"
      )
      
      # Add filled columns to filled3 if they exist in filled2
      for (col_name in names(optional_cols)) {
        if (col_name %in% names(filled2)) {
          filled3[[optional_cols[[col_name]]]] <- filled2[[col_name]]
        }
      }
      
      # Add PAR_f if filled version exists (meaning data had some non-NA values and was filled)
      if ("PAR_f" %in% names(filled2)) {
        filled3$PAR_f <- filled2$PAR_f
      } else if ("PAR" %in% original_cols_present) {
        # If PAR column exists in original but was all NA, add empty column (all NA)
        filled3$PAR <- df$PAR
      }
      
      # Add NETRAD_f if filled version exists (meaning data had some non-NA values and was filled)
      if ("NETRAD_f" %in% names(filled2)) {
        filled3$NETRAD_f <- filled2$NETRAD_f
      } else if ("NETRAD" %in% original_cols_present) {
        # If NETRAD column exists in original but was all NA, add empty column (all NA)
        filled3$NETRAD <- df$NETRAD
      }
      
      # Ensure all NA values are properly represented as empty cells in CSV
      # write.csv will automatically write NA as blank cells when na = "" is used
      output_file <- file.path(save_path, paste0(SITE_ID, "_fill.csv"))
      write.csv(filled3, output_file, row.names = FALSE, na = "")
      cat("filled data saved to:", output_file, "\n")
      print(paste("Standard offset:", standard_offset))
      
      # Print summary of columns added
      cat("\nColumns in output file:\n")
      cat(paste(names(filled3), collapse = ", "), "\n\n")
      
    }, error = function(e) {
      message("Failed to process file: ", full_file_path)
      message("Error: ", e$message)
    })
  }
}

# Run the process with the files from the specified directory
reddy_proc(file_path, save_path, files)