#################################################################################
## without loop

rm(list=ls(all=TRUE)) 
library(Rcpp)
library(REddyProc)
library(lubridate)
library("amerifluxr")
library(lutz)

# Get general information for all AmeriFlux sites
site_info <- amf_site_info()
# Specify file paths

file_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/drivers/ameri_drivers/reddy_proc"
save_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/drivers/ameri_drivers/reddy_proc/gap_filled"

# Specify the single file to process
file_name <- "US-HB3.csv"  # Change this to any site you want to process
full_file_path <- file.path(file_path, file_name)

# Create a local copy of the file to avoid permission issues
local_temp_dir <- tempdir()
local_file_path <- file.path(local_temp_dir, file_name)

# Copy the file to local temporary directory
cat("Copying file to local temporary directory...\n")
file.copy_result <- file.copy(full_file_path, local_file_path)

if (file.copy_result) {
  cat("Successfully copied file to:", local_file_path, "\n")
  cat("Local file exists:", file.exists(local_file_path), "\n")
  cat("Local file read permission:", file.access(local_file_path, 4) == 0, "\n\n")
  
  # Now use the local file path instead
  full_file_path <- local_file_path
} else {
  stop("Failed to copy file to local directory. Please check if you have write permissions in your temp directory.")
}

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

## Robust uStar threshold estimation that handles various edge cases
estimateUStarThreshold <- function(EddyProc.C, EddyDataWithPosix.F) {
  # Function to validate uStar threshold results
  validate_uStar_threshold <- function(uStarTh) {
    if (!is.data.frame(uStarTh) || !"uStar" %in% names(uStarTh)) {
      return(FALSE)
    }
    
    uStar_values <- uStarTh$uStar
    if (all(is.na(uStar_values)) || any(!is.finite(uStar_values)) || length(uStar_values) == 0) {
      return(FALSE)
    }
    
    return(TRUE)
  }
  
  # Try multiple estimation methods with different parameters
  estimation_methods <- list(
    # Method 1: Default estimation
    list(
      name = "Default estimation",
      func = function() EddyProc.C$sEstUstarThold(),
      params = NULL
    ),
    
    # Method 2: Custom parameters 1
    list(
      name = "Custom parameters (aggressive)",
      func = function() {
        seasonFactor <- usCreateSeasonFactorMonth(EddyDataWithPosix.F$DateTime, startMonth = 1)
        ctrlUstarSub <- usControlUstarSubsetting(
          taClasses = 5, 
          UstarClasses = 15,
          swThr = 5, 
          minRecordsWithinTemp = 5,
          minRecordsWithinSeason = 50, 
          minRecordsWithinYear = 1000,
          isUsingOneBigSeasonOnFewRecords = TRUE
        )
        EddyProc.C$sEstUstarThold(seasonFactor = seasonFactor, ctrlUstarSub = ctrlUstarSub)
      },
      params = NULL
    ),
    
    # Method 3: Custom parameters 2 (more lenient)
    list(
      name = "Custom parameters (lenient)",
      func = function() {
        seasonFactor <- usCreateSeasonFactorMonth(EddyDataWithPosix.F$DateTime, startMonth = 1)
        ctrlUstarSub <- usControlUstarSubsetting(
          taClasses = 10, 
          UstarClasses = 25,
          swThr = 15, 
          minRecordsWithinTemp = 3,
          minRecordsWithinSeason = 30, 
          minRecordsWithinYear = 500,
          isUsingOneBigSeasonOnFewRecords = TRUE
        )
        EddyProc.C$sEstUstarThold(seasonFactor = seasonFactor, ctrlUstarSub = ctrlUstarSub)
      },
      params = NULL
    )
  )
  
  # Try each estimation method
  for (method in estimation_methods) {
    message("Trying ", method$name, "...")
    
    uStarTh <- tryCatch({
      result <- method$func()
      if (validate_uStar_threshold(result)) {
        message("✓ ", method$name, " succeeded with valid uStar thresholds")
        return(result)
      } else {
        message("✗ ", method$name, " returned invalid values")
        NULL
      }
    }, error = function(e) {
      message("✗ ", method$name, " failed: ", e$message)
      NULL
    })
    
    if (!is.null(uStarTh)) {
      return(uStarTh)
    }
  }
  
  # If all methods failed, try to calculate a reasonable default based on data statistics
  message("All estimation methods failed. Calculating default uStar from data statistics...")
  
  # Calculate default uStar based on quantiles of available Ustar data
  ustar_data <- EddyDataWithPosix.F$Ustar
  ustar_data <- ustar_data[is.finite(ustar_data) & !is.na(ustar_data)]
  
  if (length(ustar_data) > 0) {
    # Use a reasonable quantile of the available ustar data
    default_uStar <- quantile(ustar_data, probs = 0.25, na.rm = TRUE)
    message("Using data-derived default uStar = ", round(default_uStar, 3))
    return(as.numeric(default_uStar))
  } else {
    # Final fallback - use site-specific default based on common values
    message("No Ustar data available. Using site-agnostic default uStar = 0.25")
    return(0.25)
  }
}

# Modified function to process a single file with comprehensive error handling
reddy_proc_single <- function(full_file_path, save_path) {
  tryCatch({
    # First, try to read the file
    cat("Attempting to read file:", full_file_path, "\n")
    
    data <- try(read.csv(full_file_path), silent = TRUE)
    if (inherits(data, "try-error")) {
      stop("Failed to read CSV file: ", attr(data, "condition")$message)
    }
    
    # Extract file name and SITE_ID
    file <- basename(full_file_path)
    SITE_ID <- sub("\\.csv$", "", file)
    message("Processing file: ", file, " | SITE_ID: ", SITE_ID)
    
    df <- data[, -1]  # Remove the first column
    
    # Convert time
    EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
    
    # Define and filter required columns
    required_columns <- c('NEE', 'LE', 'PAR', 'H', 'Rg', 'Tair', 'VPD', 'Ustar','RH', 'PA', 'NETRAD', 'Co2', 'WS', 'WD')
    available_columns <- required_columns[required_columns %in% colnames(EddyDataWithPosix.F)]
    
    # Ensure critical columns are numeric
    numeric_columns <- c('PA', 'Ustar', 'Rg', 'Tair')
    for (col in numeric_columns) {
      if (col %in% colnames(EddyDataWithPosix.F) && !is.numeric(EddyDataWithPosix.F[[col]])) {
        EddyDataWithPosix.F[[col]] <- as.numeric(EddyDataWithPosix.F[[col]])
      }
    }
    
    # Drop columns from available_columns that are all NA
    available_columns <- available_columns[!sapply(EddyDataWithPosix.F[available_columns], function(col) all(is.na(col)))]
    
    # Check if we have enough data for processing
    if (!'Ustar' %in% available_columns) {
      message("Warning: Ustar column not available or all NA. Using default uStar value.")
    }
    
    # Initialize sEddyProc
    EddyProc.C <- sEddyProc$new('df', EddyDataWithPosix.F, available_columns, DTS = 48)
    
    # Set location information
    site_info_specific <- site_info[site_info$SITE_ID == SITE_ID, ]
    
    message("Site info rows found: ", nrow(site_info_specific))
    if (nrow(site_info_specific) > 0) {
      message("Site found: ", site_info_specific$SITE_ID)
    } else {
      # Try case-insensitive match
      site_info_specific <- site_info[tolower(site_info$SITE_ID) == tolower(SITE_ID), ]
      if (nrow(site_info_specific) > 0) {
        SITE_ID <- site_info_specific$SITE_ID[1]
        message("Using site ID: ", SITE_ID)
      } else {
        stop("Site ", SITE_ID, " not found in site_info. Cannot proceed without coordinates.")
      }
    }
    
    if (nrow(site_info_specific) > 0) {
      lat <- site_info_specific$LOCATION_LAT[1]
      lon <- site_info_specific$LOCATION_LONG[1]
    }
    
    message("Site coordinates: Lat = ", lat, ", Lon = ", lon)
    
    timezone_info <- tz_lookup_coords(lat, lon, method = "accurate")
    standard_time <- as.POSIXct("2025-01-01", tz = timezone_info)
    standard_offset <- as.numeric(format(standard_time, "%z")) / 100
    
    EddyProc.C$sSetLocationInfo(Lat_deg.n = lat, Long_deg.n = lon, TimeZone_h.n = standard_offset)
    message("Location set: Lat = ", lat, ", Lon = ", lon, ", Timezone = ", timezone_info, ", Offset = ", standard_offset)
    
    # Estimate uStar threshold
    message("Estimating uStar threshold...")
    uStarThreshold <- estimateUStarThreshold(EddyProc.C, EddyDataWithPosix.F)
    
    # Perform gap-filling
    if (is.numeric(uStarThreshold)) {
      message("Using numeric uStar threshold: ", uStarThreshold)
      EddyProc.C$sMDSGapFillAfterUstar('H', uStarTh = uStarThreshold)
    } else if (is.data.frame(uStarThreshold) && "uStar" %in% names(uStarThreshold)) {
      message("Using uStar threshold data frame")
      EddyProc.C$sMDSGapFillAfterUstar('H')
    } else {
      message("Invalid uStar threshold. Using safe default uStar = 0.25")
      EddyProc.C$sMDSGapFillAfterUstar('H', uStarTh = 0.25)
    }
    
    # Gap fill other variables
    EddyProc.C$sMDSGapFill('PAR')
    EddyProc.C$sMDSGapFill('RH')
    EddyProc.C$sMDSGapFill('PA')
    
    # Export results
    filled2 <- EddyProc.C$sExportResults()
    
    # Create output dataframe
    filled3 <- data.frame(
      DateTime = EddyDataWithPosix.F$DateTime,
      Year     = df$Year,
      DoY      = df$DoY,
      Hour     = df$Hour,
      Tair     = df$Tair,
      VPD      = df$VPD,
      Rg       = df$Rg,
      Ustar    = df$Ustar,
      NETRAD   = df$NETRAD,
      WS       = df$WS,
      WD       = df$WD,
      Co2      = df$Co2,
      LE       = df$LE,
      NEE      = df$NEE
    )
    
    # Add gap-filled columns
    gap_filled_cols <- list(
      H_uStar_f = 'H_f',
      PAR_f     = 'PAR_f',
      RH_f      = 'RH_f',
      PA_f      = 'PA_f'
    )
    
    for (col_name in names(gap_filled_cols)) {
      if (col_name %in% names(filled2)) {
        filled3[[gap_filled_cols[[col_name]]]] <- filled2[[col_name]]
      }
    }
    
    # Rename original columns to add _f suffix
    rename_mapping <- c(
      'Tair'    = 'Tair_f',
      'VPD'     = 'VPD_f', 
      'Rg'      = 'Rg_f',
      'NETRAD'  = 'NETRAD_f',
      'LE'      = 'LE_f'
    )
    
    for (old_name in names(rename_mapping)) {
      if (old_name %in% colnames(filled3)) {
        names(filled3)[names(filled3) == old_name] <- rename_mapping[old_name]
      }
    }
    
    # Save file
    output_file <- file.path(save_path, file)
    write.csv(filled3, output_file, row.names = FALSE)
    message("Data successfully saved to: ", output_file)
    
  }, error = function(e) {
    message("Failed to process file: ", full_file_path)
    message("Error: ", e$message)
    traceback()
  })
}

# Process the file
cat("Processing file...\n")
reddy_proc_single(full_file_path, save_path)

# Clean up
cat("Cleaning up temporary file...\n")
if (file.exists(local_file_path)) {
  file.remove(local_file_path)
  cat("Temporary file removed:", !file.exists(local_file_path), "\n")
}


