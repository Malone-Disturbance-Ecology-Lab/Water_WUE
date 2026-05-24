rm(list=ls(all=TRUE)) 
library(Rcpp)
library(REddyProc)
#packageVersion("REddyProc") 1.3.3
library(lubridate)
library(amerifluxr)
library(lutz)

# Get general information for all AmeriFlux sites
site_info <- amf_site_info()

file_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/ameri_data/reddy_gaps/blended_gaps2"
save_path <- "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/ameri_data/ameri_fill"

files <- list.files(file_path, pattern = "*.csv", full.names = TRUE)

# ========================================
# Helper Functions
# ========================================

#' Detect problematic year in uStar estimation
estimateUStarThreshold_start <- function(EddyProc.C, EddyDataWithPosix.F) {
  tryCatch({
    EddyProc.C$sEstUstarThold()
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

#' Estimate uStar threshold with fallbacks
estimateUStarThreshold <- function(EddyProc.C, EddyDataWithPosix.F) {
  tryCatch({
    uStarTh_default <- EddyProc.C$sEstUstarThold()
    
    if (is.data.frame(uStarTh_default) && "uStar" %in% names(uStarTh_default)) {
      uStar_values <- uStarTh_default$uStar
    } else {
      stop("Unexpected output format from sEstUstarThold()")
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
        return(0.1)
      } else {
        return(uStarTh)
      }
    } else {
      message("Valid UStar threshold found. Using default estimation.")
      return(uStarTh_default)
    }
  }, error = function(e) {
    message("uStar estimation failed: ", e$message)
    message("Using default uStar = 0.1")
    return(0.1)
  })
}

#' Check data quality before partitioning
check_data_quality <- function(EddyProc.C, EddyDataWithPosix.F, SITE_ID) {
  message("\n--- Data Quality Check for ", SITE_ID, " ---")
  
  filled_results <- EddyProc.C$sExportResults()
  
  # Check NEE quality
  if ("NEE_uStar_f" %in% names(filled_results)) {
    nee_valid <- sum(!is.na(filled_results$NEE_uStar_f))
    nee_total <- length(filled_results$NEE_uStar_f)
    nee_pct <- round(100 * nee_valid / nee_total, 1)
    message("NEE_uStar_f valid records: ", nee_valid, "/", nee_total, " (", nee_pct, "%)")
  }
  
  # Check day/night classification
  valid_night_NEE <- 0
  if ("Rg" %in% colnames(EddyDataWithPosix.F)) {
    daytime_idx <- which(EddyDataWithPosix.F$Rg > 10 & !is.na(EddyDataWithPosix.F$Rg))
    nighttime_idx <- which(EddyDataWithPosix.F$Rg <= 10 & !is.na(EddyDataWithPosix.F$Rg))
    
    message("Daytime records (Rg>10): ", length(daytime_idx), 
            " (", round(100*length(daytime_idx)/nrow(EddyDataWithPosix.F), 1), "%)")
    message("Nighttime records (Rg<=10): ", length(nighttime_idx), 
            " (", round(100*length(nighttime_idx)/nrow(EddyDataWithPosix.F), 1), "%)")
    
    # Check valid NEE during nighttime
    if ("NEE_uStar_f" %in% names(filled_results)) {
      valid_day_NEE <- sum(!is.na(filled_results$NEE_uStar_f[daytime_idx]))
      valid_night_NEE <- sum(!is.na(filled_results$NEE_uStar_f[nighttime_idx]))
      
      message("Valid NEE_uStar_f during daytime: ", valid_day_NEE, "/", length(daytime_idx), 
              " (", round(100*valid_day_NEE/length(daytime_idx), 1), "%)")
      message("Valid NEE_uStar_f during nighttime: ", valid_night_NEE, "/", length(nighttime_idx), 
              " (", round(100*valid_night_NEE/length(nighttime_idx), 1), "%)")
    }
  }
  
  # Check temperature data
  if ("Tair_f" %in% names(filled_results)) {
    tair_valid <- sum(!is.na(filled_results$Tair_f))
    message("Tair_f valid records: ", tair_valid, "/", length(filled_results$Tair_f), 
            " (", round(100*tair_valid/length(filled_results$Tair_f), 1), "%)")
  }
  
  message("-----------------------------------\n")
  return(valid_night_NEE)
}

#' Extract E0 from MR partitioning result
extract_E0_from_MR <- function(EddyProc.C) {
  tryCatch({
    # Try to get E0 from the internal object
    if (exists("E0_est", where = EddyProc.C)) {
      e0_val <- EddyProc.C$E0_est
      if (is.finite(e0_val) && e0_val > 50 && e0_val < 450) {
        message("Extracted E0 from MR object: ", round(e0_val, 2))
        return(e0_val)
      }
    }
    
    # Try from results
    res <- EddyProc.C$sExportResults()
    e0_cols <- grep("^E0_|^E_0_|FP_E0", names(res), value = TRUE)
    if (length(e0_cols) > 0) {
      e0_vals <- res[[e0_cols[1]]]
      e0_val <- median(e0_vals, na.rm = TRUE)
      if (is.finite(e0_val) && e0_val > 50 && e0_val < 450) {
        message("Extracted E0 from MR result: ", round(e0_val, 2))
        return(e0_val)
      }
    }
    
    message("Could not extract valid E0 from MR result. Using default 171.01.")
    return(171.01)
  }, error = function(e) {
    message("Error extracting E0: ", e$message, ". Using default 171.01.")
    return(171.01)
  })
}

# ========================================
# Main Processing Function
# ========================================

reddy_proc <- function(file_path, save_path, files) {
  for (full_file_path in files) {
    tryCatch({
      file <- basename(full_file_path)
      SITE_ID <- sub("^gaps_blend_(.*)\\.csv$", "\\1", file)
      message("\n========================================")
      message("Processing file: ", file, " | SITE_ID: ", SITE_ID)
      message("========================================")
      
      # Load data
      data <- read.csv(full_file_path, na.strings = c("NA", "NaN", "", " "), stringsAsFactors = FALSE)
      df <- data[, -1]  # Remove the first column (usually index column)
      
      # Initial time conversion
      EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
      
      # Define and filter required columns
      required_columns <- c('NEE', 'LE', 'H', 'Rg', 'Tair', 'VPD', 'Ustar', 'PA', 'NETRAD', 'WS', 'WD', 'RH', 'PAR')
      available_columns <- required_columns[required_columns %in% colnames(EddyDataWithPosix.F)]
      
      # Ensure PA is numeric
      if (!is.numeric(EddyDataWithPosix.F$PA)) {
        EddyDataWithPosix.F$PA <- as.numeric(EddyDataWithPosix.F$PA)
      }
      
      # Drop columns that are all NA
      available_columns <- available_columns[!sapply(EddyDataWithPosix.F[available_columns], function(col) all(is.na(col)))]
      
      # Check if we have enough data
      if (nrow(df) < 100) {
        message("Insufficient data rows (<100) for site ", SITE_ID, ". Skipping.")
        next
      }
      
      EddyProc.C <- sEddyProc$new('df', EddyDataWithPosix.F, available_columns, DTS = 48)
      
      # Handle problematic years
      problematic_info <- estimateUStarThreshold_start(EddyProc.C, EddyDataWithPosix.F)
      if (!is.null(problematic_info)) {
        if (is.numeric(problematic_info)) {
          problematic_years <- unique(na.omit(problematic_info))
        } else if (is.data.frame(problematic_info) && "year" %in% names(problematic_info)) {
          problematic_years <- unique(na.omit(problematic_info$year))
        } else {
          problematic_years <- NULL
        }
        
        if (!is.null(problematic_years) && length(problematic_years) > 0) {
          df_remaining <- df[!df$Year %in% problematic_years, ]
          if (nrow(df_remaining) < 100) {
            message("Removing problematic years would leave insufficient data. Using all data.")
            problematic_years <- NULL
          } else {
            df <- df_remaining
            EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
            message("Removing problematic year(s): ", paste(problematic_years, collapse = ", "))
          }
        }
      }
      
      # Recreate EddyProc after potential data removal
      EddyDataWithPosix.F <- fConvertTimeToPosix(df, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')
      available_columns <- required_columns[required_columns %in% colnames(EddyDataWithPosix.F)]
      
      if (!is.numeric(EddyDataWithPosix.F$PA)) {
        EddyDataWithPosix.F$PA <- as.numeric(EddyDataWithPosix.F$PA)
      }
      
      available_columns <- available_columns[!sapply(EddyDataWithPosix.F[available_columns], function(col) all(is.na(col)))]
      EddyProc.C <- sEddyProc$new('df', EddyDataWithPosix.F, available_columns, DTS = 48)
      
      # ========================================
      # uStar Estimation
      # ========================================
      has_valid_ustar <- "Ustar" %in% colnames(EddyDataWithPosix.F) && 
        any(!is.na(EddyDataWithPosix.F$Ustar)) &&
        sum(!is.na(EddyDataWithPosix.F$Ustar)) >= 100
      
      if (!has_valid_ustar) {
        message("Insufficient valid uStar records. Using default uStar = 0.1")
        uStar_result <- 0.1
        use_default_ustar <- TRUE
      } else {
        uStar_result <- tryCatch({
          estimateUStarThreshold(EddyProc.C, EddyDataWithPosix.F)
        }, error = function(e) {
          message("uStar estimation failed: ", e$message)
          return(0.1)
        })
        
        # Determine if we should use default based on result type
        if (is.numeric(uStar_result) && length(uStar_result) == 1 && uStar_result == 0.1) {
          use_default_ustar <- TRUE
        } else if (is.data.frame(uStar_result) && "uStar" %in% names(uStar_result)) {
          if (all(is.na(uStar_result$uStar))) {
            message("uStar threshold values are all NA. Using default 0.1")
            uStar_result <- 0.1
            use_default_ustar <- TRUE
          } else if (all(uStar_result$uStar == 0.1, na.rm = TRUE)) {
            use_default_ustar <- TRUE
          } else {
            use_default_ustar <- FALSE
          }
        } else if (is.numeric(uStar_result) && length(uStar_result) > 1) {
          use_default_ustar <- FALSE
        } else {
          message("Unexpected uStar estimation result. Using default 0.1")
          uStar_result <- 0.1
          use_default_ustar <- TRUE
        }
      }
      
      # ========================================
      # Gap-filling with uStar Filtering
      # ========================================
      if (use_default_ustar) {
        message("Using default uStar = 0.1 for gap-filling")
        tryCatch({
          EddyProc.C$sMDSGapFillAfterUstar('NEE', uStarTh = 0.1)
          EddyProc.C$sExportResults()
          EddyProc.C$sMDSGapFillAfterUstar('LE', uStarTh = 0.1)
          EddyProc.C$sMDSGapFillAfterUstar('H', uStarTh = 0.1)
        }, error = function(e) {
          message("Error in uStar filtering: ", e$message, ". Falling back to simple gap-fill.")
          EddyProc.C$sMDSGapFill('NEE')
          EddyProc.C$sMDSGapFill('LE')
          EddyProc.C$sMDSGapFill('H')
        })
      } else {
        message("Using estimated uStar thresholds")
        tryCatch({
          EddyProc.C$sMDSGapFillAfterUstar('NEE')
          EddyProc.C$sExportResults()
          EddyProc.C$sMDSGapFillAfterUstar('LE')
          EddyProc.C$sMDSGapFillAfterUstar('H')
        }, error = function(e) {
          message("Error in uStar filtering: ", e$message, ". Using simple gap-fill.")
          EddyProc.C$sMDSGapFill('NEE')
          EddyProc.C$sMDSGapFill('LE')
          EddyProc.C$sMDSGapFill('H')
        })
      }
      
      # ========================================
      # Fill Meteorological Variables
      # ========================================
      tryCatch({
        EddyProc.C$sMDSGapFill('Rg')
        EddyProc.C$sMDSGapFill('Tair', FillAll.b = FALSE)
        EddyProc.C$sMDSGapFill('VPD', FillAll.b = FALSE)
        
        if ("PAR" %in% colnames(EddyDataWithPosix.F) && !all(is.na(EddyDataWithPosix.F$PAR))) {
          EddyProc.C$sMDSGapFill('PAR', FillAll.b = FALSE)
        }
        
        if ("NETRAD" %in% colnames(EddyDataWithPosix.F) && !all(is.na(EddyDataWithPosix.F$NETRAD))) {
          EddyProc.C$sMDSGapFill('NETRAD', FillAll.b = FALSE)
        }
      }, error = function(e) {
        message("Warning: Some meteorological gap-filling failed: ", e$message)
      })
      
      # ========================================
      # Get Site Location Info
      # ========================================
      site_info_specific <- site_info[site_info$SITE_ID == SITE_ID, ]
      if (nrow(site_info_specific) == 0) {
        stop_message <- paste0("No location info found for site ", SITE_ID, ". Cannot proceed without coordinates.")
        message("\n!!! ERROR: ", stop_message)
        stop(stop_message)
      } else {
        lat <- site_info_specific$LOCATION_LAT
        lon <- site_info_specific$LOCATION_LONG
        message("Site location found: Lat = ", lat, ", Lon = ", lon)
      }
      
      timezone_info <- tz_lookup_coords(lat, lon, method = "accurate")
      standard_time <- as.POSIXct("2025-01-01", tz = timezone_info)
      standard_offset <- as.numeric(format(standard_time, "%z")) / 100
      
      # Set location
      EddyProc.C$sSetLocationInfo(Lat_deg.n = lat, Long_deg.n = lon, TimeZone_h.n = standard_offset)
      
      # ========================================
      # Diagnostic Check Before Partitioning
      # ========================================
      valid_night_NEE <- check_data_quality(EddyProc.C, EddyDataWithPosix.F, SITE_ID)
      
      # ========================================
      # FLUX PARTITIONING
      # ========================================
      message("\n--- Starting Flux Partitioning ---")
      
      # Initialize flags and storage
      nighttime_success <- FALSE
      daytime_success <- FALSE
      mr_E0 <- NA
      
      # Create vectors to store nighttime results
      GPP_nt <- rep(NA, nrow(EddyDataWithPosix.F))
      Reco_nt <- rep(NA, nrow(EddyDataWithPosix.F))
      
      # ========================================
      # NIGHTTIME PARTITIONING (MR Method)
      # ========================================
      tryCatch({
        message("Attempting nighttime flux partitioning (MR method)...")
        EddyProc.C$sMRFluxPartition(Suffix.s = 'uStar')
        nighttime_success <- TRUE
        message("✓ Nighttime flux partitioning (MR) SUCCESSFUL")
        
        # Extract nighttime results immediately after MR partitioning
        temp_results <- EddyProc.C$sExportResults()
        
        # Check for GPP and Reco columns
        if ("GPP_uStar_f" %in% names(temp_results)) {
          GPP_nt <- temp_results$GPP_uStar_f
          message("  ✓ Found GPP_uStar_f in results")
        } else if ("GPP_NT_uStar" %in% names(temp_results)) {
          GPP_nt <- temp_results$GPP_NT_uStar
          message("  ✓ Found GPP_NT_uStar in results")
        } else {
          message("  WARNING: No GPP column found after MR partitioning")
        }
        
        if ("Reco_uStar" %in% names(temp_results)) {
          Reco_nt <- temp_results$Reco_uStar
          message("  ✓ Found Reco_uStar in results")
        } else if ("Reco_NT_uStar" %in% names(temp_results)) {
          Reco_nt <- temp_results$Reco_NT_uStar
          message("  ✓ Found Reco_NT_uStar in results")
        } else {
          message("  WARNING: No Reco column found after MR partitioning")
        }
        
        # Extract E0 for use in daytime partitioning
        mr_E0 <- extract_E0_from_MR(EddyProc.C)
        
      }, error = function(e) {
        message("✗ Nighttime flux partitioning (MR) FAILED: ", e$message)
        nighttime_success <- FALSE
      })
      
      # ========================================
      # DAYTIME PARTITIONING (GL Method) - ALWAYS ATTEMPT
      # ========================================
      message("\nAttempting daytime flux partitioning (GL method)...")
      
      # Prepare E0 for fixed-temp strategies (use MR E0 if available, otherwise default)
      e0_for_GL <- if (is.finite(mr_E0)) mr_E0 else 171.01
      
      # Multiple strategies in order of preference
      gl_success <- FALSE
      
      # Strategy 1: Lasslop-compatible (requires good nighttime data)
      if (!gl_success && valid_night_NEE >= 5000) {
        tryCatch({
          message("  Strategy 1: Lasslop-compatible settings...")
          ctrlGL <- partGLControlLasslopCompatible()
          EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar', controlGLPart = ctrlGL)
          gl_success <- TRUE
          daytime_success <- TRUE
          message("  ✓ GL succeeded with Lasslop-compatible settings")
        }, error = function(e) {
          message("  Strategy 1 failed: ", e$message)
        })
      }
      
      # Strategy 2: Manual controls with default parameters
      if (!gl_success) {
        tryCatch({
          message("  Strategy 2: Manual GL controls...")
          ctrlGL <- partGLControl(
            smoothTempSensEstimateAcrossTime = FALSE,
            nBootUncertainty = 0,
            minNRecInDayWindow = 10,
            isRefitMissingVPDWithNeglectVPDEffect = FALSE,
            isNeglectVPDEffect = FALSE
          )
          EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar', controlGLPart = ctrlGL)
          gl_success <- TRUE
          daytime_success <- TRUE
          message("  ✓ GL succeeded with manual controls")
        }, error = function(e) {
          message("  Strategy 2 failed: ", e$message)
        })
      }
      
      # Strategy 3: Fixed temperature sensitivity (use E0 from MR or default)
      if (!gl_success) {
        tryCatch({
          message("  Strategy 3: Fixed E0 = ", round(e0_for_GL, 2), "...")
          ctrlGL_fixed <- partGLControl(
            smoothTempSensEstimateAcrossTime = FALSE,
            fixedTempSens = data.frame(E0 = e0_for_GL, sdE0 = 50),
            nBootUncertainty = 0,
            minNRecInDayWindow = 10,
            isRefitMissingVPDWithNeglectVPDEffect = FALSE,
            isNeglectVPDEffect = FALSE
          )
          EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar', controlGLPart = ctrlGL_fixed)
          gl_success <- TRUE
          daytime_success <- TRUE
          message("  ✓ GL succeeded with fixed E0")
        }, error = function(e) {
          message("  Strategy 3 failed: ", e$message)
        })
      }
      
      # Strategy 4: Neglect VPD effect (helps when VPD data is poor)
      if (!gl_success) {
        tryCatch({
          message("  Strategy 4: Neglect VPD effect...")
          ctrlGL_novpd <- partGLControl(
            smoothTempSensEstimateAcrossTime = FALSE,
            fixedTempSens = data.frame(E0 = e0_for_GL, sdE0 = 50),
            nBootUncertainty = 0,
            minNRecInDayWindow = 10,
            isRefitMissingVPDWithNeglectVPDEffect = TRUE,
            isNeglectVPDEffect = TRUE
          )
          EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar', controlGLPart = ctrlGL_novpd)
          gl_success <- TRUE
          daytime_success <- TRUE
          message("  ✓ GL succeeded with VPD effect neglected")
        }, error = function(e) {
          message("  Strategy 4 failed: ", e$message)
        })
      }
      
      # Strategy 5: Use PAR instead of Rg for radiation
      if (!gl_success) {
        tryCatch({
          message("  Strategy 5: Using PAR as radiation variable...")
          ctrlGL_par <- partGLControl(
            smoothTempSensEstimateAcrossTime = FALSE,
            fixedTempSens = data.frame(E0 = e0_for_GL, sdE0 = 50),
            nBootUncertainty = 0,
            minNRecInDayWindow = 10,
            isRefitMissingVPDWithNeglectVPDEffect = TRUE,
            isNeglectVPDEffect = TRUE
          )
          EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar', RadVar = "PAR_f", controlGLPart = ctrlGL_par)
          gl_success <- TRUE
          daytime_success <- TRUE
          message("  ✓ GL succeeded with PAR as radiation")
        }, error = function(e) {
          message("  Strategy 5 failed: ", e$message)
        })
      }
      
      # Strategy 6: Minimal constraints (for sites with zero uncertainty issues)
      if (!gl_success) {
        tryCatch({
          message("  Strategy 6: Minimal constraints (for difficult sites)...")
          ctrlGL_minimal <- partGLControl(
            smoothTempSensEstimateAcrossTime = FALSE,
            fixedTempSens = data.frame(E0 = e0_for_GL, sdE0 = 50),
            nBootUncertainty = 0,
            minNRecInDayWindow = 5,  # Reduced minimum records needed
            isRefitMissingVPDWithNeglectVPDEffect = TRUE,
            isNeglectVPDEffect = TRUE
          )
          EddyProc.C$sGLFluxPartition(Suffix.s = 'uStar', controlGLPart = ctrlGL_minimal)
          gl_success <- TRUE
          daytime_success <- TRUE
          message("  ✓ GL succeeded with minimal constraints")
        }, error = function(e) {
          message("  Strategy 6 failed: ", e$message)
        })
      }
      
      # Final check
      if (!gl_success) {
        message("✗ All GL strategies FAILED - No daytime partitioning results")
        daytime_success <- FALSE
      }
      
      # Export final results after GL partitioning
      final_results <- EddyProc.C$sExportResults()
      
      # ========================================
      # Partitioning Summary
      # ========================================
      message("\n--- Flux Partitioning Summary ---")
      if (nighttime_success && daytime_success) {
        message("✓ Both partitioning methods succeeded - Full outputs available")
      } else if (nighttime_success && !daytime_success) {
        message("⚠ Only nighttime partitioning succeeded - GPP_nt and Reco_nt available")
      } else if (!nighttime_success && daytime_success) {
        message("⚠ Only daytime partitioning succeeded - GPP_DT and Reco_DT available")
      } else {
        message("✗ Both partitioning methods failed - No GPP/Reco outputs")
      }
      
      # ========================================
      # Build Final DataFrame
      # ========================================
      filled3 <- data.frame(
        DateTime = df$TIMESTAMP,
        Year = df$Year,
        DoY = df$DoY,
        Hour = df$Hour,
        stringsAsFactors = FALSE
      )
      
      # Add original variables
      original_vars <- c("Ustar", "PA", "RH", "WS", "WD")
      for (var in original_vars) {
        if (var %in% colnames(df)) {
          filled3[[var]] <- df[[var]]
        }
      }
      
      # Add filled columns
      optional_cols <- list(
        NEE_uStar_f = "NEE_f",
        LE_uStar_f = "LE_f",
        H_uStar_f = "H_f",
        Tair_f = "Tair_f",
        VPD_f = "VPD_f",
        Rg_f = "Rg_f"
      )
      
      for (col_name in names(optional_cols)) {
        if (col_name %in% names(final_results)) {
          filled3[[optional_cols[[col_name]]]] <- final_results[[col_name]]
        }
      }
      
      # Add nighttime partitioning results
      if (nighttime_success) {
        # Check if GPP_nt has any non-NA values
        if (any(!is.na(GPP_nt))) {
          filled3$GPP_nt <- GPP_nt
          message("✓ Added GPP_nt (nighttime method) with ", sum(!is.na(GPP_nt)), " valid values")
        } else {
          message("✗ GPP_nt has no valid values - check MR partitioning")
        }
        
        if (any(!is.na(Reco_nt))) {
          filled3$Reco_nt <- Reco_nt
          message("✓ Added Reco_nt (nighttime method) with ", sum(!is.na(Reco_nt)), " valid values")
        } else {
          message("✗ Reco_nt has no valid values - check MR partitioning")
        }
      }
      
      # Add daytime partitioning results
      if (daytime_success) {
        if ("GPP_DT_uStar" %in% names(final_results)) {
          filled3$GPP_DT <- final_results$GPP_DT_uStar
          message("✓ Added GPP_DT (daytime method)")
        }
        if ("Reco_DT_uStar" %in% names(final_results)) {
          filled3$Reco_DT <- final_results$Reco_DT_uStar
          message("✓ Added Reco_DT (daytime method)")
        }
      }
      
      # Add PAR and NETRAD
      if ("PAR_f" %in% names(final_results)) {
        filled3$PAR_f <- final_results$PAR_f
      } else if ("PAR" %in% colnames(df)) {
        filled3$PAR <- df$PAR
      }
      
      if ("NETRAD_f" %in% names(final_results)) {
        filled3$NETRAD_f <- final_results$NETRAD_f
      } else if ("NETRAD" %in% colnames(df)) {
        filled3$NETRAD <- df$NETRAD
      }
      
      # ========================================
      # Save File - ALWAYS HAPPENS
      # ========================================
      output_file <- file.path(save_path, paste0(SITE_ID, "_fill.csv"))
      write.csv(filled3, output_file, row.names = FALSE, na = "")
      
      cat("\n========================================")
      cat("\n✓ SUCCESSFULLY SAVED: ", output_file)
      cat("\n  Columns saved:", paste(names(filled3), collapse = ", "))
      cat("\n  Partitioning status: Nighttime=", nighttime_success, "| Daytime=", daytime_success)
      if (nighttime_success) {
        cat("\n  Nighttime valid values: GPP_nt=", if("GPP_nt" %in% names(filled3)) sum(!is.na(filled3$GPP_nt)) else 0, 
            "| Reco_nt=", if("Reco_nt" %in% names(filled3)) sum(!is.na(filled3$Reco_nt)) else 0)
      }
      cat("\n========================================\n")
      
    }, error = function(e) {
      message("\n!!! CRITICAL ERROR for file: ", full_file_path)
      message("Error: ", e$message)
      message("This file was skipped. Continuing to next file...\n")
    })
  }
}

# Run the process
reddy_proc(file_path, save_path, files)