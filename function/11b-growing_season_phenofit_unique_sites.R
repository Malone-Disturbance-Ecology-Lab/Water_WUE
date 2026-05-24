# Rescue version for difficult/low-quality sites
# This version processes sites from the rescue_sites folder only

phenofit_growing_season_rescue <- function(input_folder, save_folder) {
  library(tidyverse)
  library(phenofit)
  library(data.table)
  library(lubridate)
  library(zoo)  # For rollmean smoothing
  
  # Create save folder if it doesn't exist
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  
  # Function to average duplicate years
  average_duplicate_years <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      mutate(year = as.integer(format(origin, "%Y"))) %>%
      group_by(year) %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  }
  
  # Function to safely parse DateTime with multiple format attempts (NO UTC)
  safe_parse_datetime <- function(datetime_vec) {
    # Try standard format first (YYYY-MM-DD HH:MM:SS)
    parsed <- as.POSIXct(datetime_vec, format = "%Y-%m-%d %H:%M:%S")
    
    # If that fails, try without seconds (YYYY-MM-DD HH:MM)
    if (all(is.na(parsed))) {
      parsed <- as.POSIXct(datetime_vec, format = "%Y-%m-%d %H:%M")
    }
    
    # If still failing, try with lubridate's ymd_hms (no tz)
    if (all(is.na(parsed))) {
      parsed <- ymd_hms(datetime_vec, quiet = TRUE)
    }
    
    # Try ymd_hm if that fails
    if (all(is.na(parsed))) {
      parsed <- ymd_hm(datetime_vec, quiet = TRUE)
    }
    
    return(parsed)
  }
  
  # Function to screen and keep only good years
  screen_good_years <- function(daily_df, site_name) {
    # Check if daily_df has any data
    if (is.null(daily_df) || nrow(daily_df) == 0) {
      cat("  No data available for screening\n")
      return(NULL)
    }
    
    # Add a small smoothing step to reduce noise
    daily_df <- daily_df %>%
      arrange(Date) %>%
      mutate(GPP_smoothed = rollmean(GPP, k = min(7, n()), fill = NA, align = "center")) %>%
      mutate(GPP = ifelse(!is.na(GPP_smoothed), GPP_smoothed, GPP)) %>%
      select(-GPP_smoothed) %>%
      filter(!is.na(GPP))
    
    # Add year column
    daily_df <- daily_df %>%
      mutate(year = year(Date))
    
    # Get all years
    all_years <- unique(daily_df$year)
    good_years <- c()
    
    for (yr in all_years) {
      year_data <- daily_df %>% filter(year == yr)
      n_days <- n_distinct(year_data$Date)
      
      # For single-year sites, keep them even if incomplete
      if (length(all_years) == 1) {
        cat("    Single-year site, keeping year", yr, "with", n_days, "days\n")
        good_years <- c(good_years, yr)
        next
      }
      
      # Skip incomplete final year only if there are multiple years
      if (yr == max(all_years) && n_days < 300) {
        cat("    Year", yr, "is incomplete (", n_days, "days) - excluding\n")
        next
      }
      
      # Regular thresholds for other sites
      if (n_days < 150) {
        cat("    Year", yr, "has only", n_days, "unique days - excluding\n")
        next
      }
      
      # Check if year has enough GPP amplitude
      gpp_range <- diff(range(year_data$GPP, na.rm = TRUE))
      
      if (is.na(gpp_range) || gpp_range < 0.01) {
        cat("    Year", yr, "has insufficient GPP amplitude (", round(gpp_range, 4), ") - excluding\n")
        next
      }
      
      # If we got here, this is a good year
      good_years <- c(good_years, yr)
      cat("    Year", yr, "accepted - amplitude:", round(gpp_range, 4), ", days:", n_days, "\n")
    }
    
    if (length(good_years) == 0) {
      cat("  No good years found for", site_name, "\n")
      return(NULL)
    }
    
    cat("  Keeping", length(good_years), "good years:", paste(good_years, collapse = ", "), "\n")
    
    # Return filtered data with only good years
    return(daily_df %>% filter(year %in% good_years))
  }
  
  # Get list of CSV files in the rescue_sites folder
  csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    stop("No CSV files found in rescue folder: ", input_folder)
  }
  
  cat("Found", length(csv_files), "rescue site CSV files to process\n")
  
  successful_sites <- list()
  failed_sites <- list()
  
  for (file_path in csv_files) {
    filename <- basename(file_path)
    site_name <- gsub("_fill_gpp_et\\.csv$", "", filename)
    site_name <- gsub("\\.csv$", "", site_name)
    
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("RESCUE MODE - Processing site:", site_name, "\n")
    cat("File:", filename, "\n")
    
    # SPECIAL CASE FOR US-StS - Use manually provided values
    if (site_name == "US-StS") {
      cat("[INFO] Special handling for US-StS - using manually verified phenology values\n")
      
      # Manually entered values based on historical data
      data_grow_season <- tibble(
        site_name = "US-StS",
        year = 2015,
        avg_DER_sos = 158.4,
        avg_DER_eos = 219.0,
        avg_TRS_sos = 155.1,
        avg_TRS_eos = 223.8,
        avg_length = 60.6,
        avg_sos = 156.75,
        avg_eos = 221.4
      )
      
      # SAVE INDIVIDUAL SITE FILE
      output_filename <- paste0(site_name, "_phenofit.csv")
      output_path <- file.path(save_folder, output_filename)
      
      write.csv(data_grow_season, output_path, row.names = FALSE)
      
      successful_sites[[site_name]] <- list(
        file = output_filename,
        years = 1,
        rows = 1,
        manual = TRUE
      )
      
      cat("[OK] RESCUE MODE - Manually processed US-StS with verified values\n")
      cat("[OK] Saved to:", output_path, "\n")
      cat("    Values: SOS = 156.75, EOS = 221.4, Length = 60.6 days\n")
      
      next  # Skip the regular processing for this site
    }
    
    # SPECIAL CASE FOR US-BFS - Use manually provided values
    if (site_name == "US-BFS") {
      cat("[INFO] Special handling for US-BFS - using manually verified phenology values\n")
      
      # Manually entered values based on your data
      data_grow_season <- tibble(
        site_name = "US-BFS",
        year = 2025,
        avg_DER_sos = 23,
        avg_DER_eos = 190,
        avg_TRS_sos = 23,
        avg_TRS_eos = 190,
        avg_length = 167,
        avg_sos = 23,
        avg_eos = 190
      )
      
      # SAVE INDIVIDUAL SITE FILE
      output_filename <- paste0(site_name, "_phenofit.csv")
      output_path <- file.path(save_folder, output_filename)
      
      write.csv(data_grow_season, output_path, row.names = FALSE)
      
      successful_sites[[site_name]] <- list(
        file = output_filename,
        years = 1,
        rows = 1,
        manual = TRUE
      )
      
      cat("[OK] RESCUE MODE - Manually processed US-BFS with verified values\n")
      cat("[OK] Saved to:", output_path, "\n")
      cat("    Values: SOS = 23, EOS = 190, Length = 167 days\n")
      
      next  # Skip the regular processing for this site
    }
    
    # Regular processing for all other rescue sites
    tryCatch({
      # Read the CSV file
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Check required columns
      if (!"DateTime" %in% colnames(df) || !"GPP" %in% colnames(df)) {
        cat("[ERROR] Missing required columns - skipping\n")
        failed_sites[[site_name]] <<- "Missing required columns"
        next
      }
      
      # Check if GPP column has any valid data
      if (all(is.na(df$GPP))) {
        cat("[ERROR] All GPP values are NA for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "All GPP values are NA"
        next
      }
      
      # Print first few DateTime values for debugging
      cat("First few DateTime values:", paste(head(df$DateTime)), "\n")
      
      # SAFE DATETIME PARSING (NO UTC)
      df$DateTime <- safe_parse_datetime(df$DateTime)
      
      # Check if parsing worked
      if (all(is.na(df$DateTime))) {
        cat("[ERROR] Failed to parse DateTime column for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "Failed to parse DateTime column"
        next
      }
      
      cat("Successfully parsed DateTime for", nrow(df), "rows\n")
      
      # Remove rows with NA DateTime
      df <- df %>% filter(!is.na(DateTime))
      
      # Check if we have any data after filtering
      if (nrow(df) == 0) {
        cat("[ERROR] No valid DateTime rows for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "No valid DateTime rows"
        next
      }
      
      # Aggregate to daily
      daily_df <- df %>%
        mutate(Date = as.Date(DateTime)) %>%
        group_by(Date) %>%
        summarise(GPP = mean(GPP, na.rm = TRUE), .groups = "drop") %>%
        filter(!is.na(GPP))
      
      # Check if daily aggregation produced any data
      if (nrow(daily_df) == 0) {
        cat("[ERROR] No valid daily GPP data for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "No valid daily GPP data"
        next
      }
      
      cat("Total days before screening:", nrow(daily_df), "\n")
      cat("Date range:", format(min(daily_df$Date), "%Y-%m-%d"), "to", 
          format(max(daily_df$Date), "%Y-%m-%d"), "\n")
      
      # SCREEN AND REMOVE BAD YEARS BEFORE FITTING
      daily_df <- screen_good_years(daily_df, site_name)
      
      if (is.null(daily_df) || nrow(daily_df) < 30) {
        cat("[ERROR] No good years found for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "No good years after screening"
        next
      }
      
      cat("Total days after screening:", nrow(daily_df), "\n")
      cat("Years after screening:", paste(unique(year(daily_df$Date)), collapse = ", "), "\n")
      
      # Add weights
      daily_df <- daily_df %>% mutate(w = 1.0)
      
      # Parameters for phenofit - CAREFULLY RELAXED for rescue mode
      nptperyear <- 365
      wFUN <- wTSM
      wmin <- 0.1  # Lower threshold
      minExtendMonth <- 1  # More relaxed: extend up to 1 month
      maxExtendMonth <- 3  # More relaxed: extend up to 3 months
      minPercValid <- 0.2  # Only need 20% valid data
      south_hemisphere <- FALSE
      
      # Use reasonable maxgap
      maxgap <- nptperyear / 4  # 91 days
      cat("Using maxgap =", maxgap, "days\n")
      
      # Check input data
      INPUT <- check_input(
        t = daily_df$Date,
        y = daily_df$GPP,
        w = daily_df$w,
        nptperyear = nptperyear,
        maxgap = maxgap,
        wmin = wmin,
        mask_spike = TRUE,
        south = south_hemisphere
      )
      
      # Unified relaxed settings for all sites
      cat("Using relaxed default settings for site", site_name, "\n")
      brks <- season_mov(INPUT,
                         list(FUN = "smooth_wWHIT", wFUN = wFUN,
                              maxExtendMonth = 3,
                              wmin = 0.1,
                              r_min = 0.02,
                              ypeak_min = 0.02,
                              lambda = 25000,
                              iters = 10,
                              MaxPeaksPerYear = 1,
                              MaxTroughsPerYear = 2))
      
      # Curve fitting using multiple methods
      fit <- curvefits(INPUT, brks,
                       list(methods = c("Beck", "Elmore", "Gu", "AG", "Zhang"),
                            wFUN = wFUN, iters = 10,
                            wmin = 0.1,
                            maxExtendMonth = maxExtendMonth,
                            minExtendMonth = minExtendMonth,
                            minPercValid = minPercValid))
      
      # Extract phenological metrics
      TRS <- c(0.1, 0.25, 0.5)
      l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = FALSE)
      
      # Check if we got valid phenology data
      if (is.null(l_pheno$doy) || 
          is.null(l_pheno$doy$Beck) || 
          nrow(l_pheno$doy$Beck) == 0) {
        cat("[WARNING] No phenology data extracted for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "No phenology data extracted"
        next
      }
      
      # Extract data frames for each method - SAFELY HANDLE MISSING METHODS
      doy_Beck <- tryCatch(average_duplicate_years(l_pheno$doy$Beck), error = function(e) NULL)
      doy_Elmore <- tryCatch(average_duplicate_years(l_pheno$doy$Elmore), error = function(e) NULL)
      doy_Gu <- tryCatch(average_duplicate_years(l_pheno$doy$Gu), error = function(e) NULL)
      doy_AG <- tryCatch(average_duplicate_years(l_pheno$doy$AG), error = function(e) NULL)
      doy_Zhang <- tryCatch(average_duplicate_years(l_pheno$doy$Zhang), error = function(e) NULL)
      
      # Create a named list of only valid methods
      valid_methods <- list()
      if (!is.null(doy_Beck)) valid_methods[["Beck"]] <- doy_Beck
      if (!is.null(doy_Elmore)) valid_methods[["Elmore"]] <- doy_Elmore
      if (!is.null(doy_Gu)) valid_methods[["Gu"]] <- doy_Gu
      if (!is.null(doy_AG)) valid_methods[["AG"]] <- doy_AG
      if (!is.null(doy_Zhang)) valid_methods[["Zhang"]] <- doy_Zhang
      
      if (length(valid_methods) < 2) {
        cat("[WARNING] Insufficient valid methods for", site_name, "- only", length(valid_methods), "methods, skipping\n")
        failed_sites[[site_name]] <<- paste("Insufficient methods:", length(valid_methods))
        next
      }
      
      cat("Using", length(valid_methods), "methods:", paste(names(valid_methods), collapse = ", "), "\n")
      
      # Find common years across all valid methods
      common_years <- Reduce(intersect, lapply(valid_methods, function(x) x$year))
      
      if (length(common_years) == 0) {
        cat("[WARNING] No common years across methods for", site_name, "- skipping\n")
        failed_sites[[site_name]] <<- "No common years across methods"
        next
      }
      
      # Filter each method to common years
      for (method_name in names(valid_methods)) {
        valid_methods[[method_name]] <- filter(valid_methods[[method_name]], year %in% common_years)
      }
      
      nyears <- length(common_years)
      cat("Found", nyears, "common years across methods\n")
      
      # Build averaged phenology results
      der_sos_matrix <- do.call(cbind, lapply(valid_methods, function(x) x$DER.sos))
      der_eos_matrix <- do.call(cbind, lapply(valid_methods, function(x) x$DER.eos))
      trs_sos_matrix <- do.call(cbind, lapply(valid_methods, function(x) x$TRS2.5.sos))
      trs_eos_matrix <- do.call(cbind, lapply(valid_methods, function(x) x$TRS2.5.eos))
      length_matrix <- do.call(cbind, lapply(valid_methods, function(x) x$DER.eos - x$DER.sos))
      
      data_grow_season <- tibble(
        site_name = rep(site_name, nyears),
        year = common_years,
        avg_DER_sos = rowMeans(der_sos_matrix, na.rm = TRUE),
        avg_DER_eos = rowMeans(der_eos_matrix, na.rm = TRUE),
        avg_TRS_sos = rowMeans(trs_sos_matrix, na.rm = TRUE),
        avg_TRS_eos = rowMeans(trs_eos_matrix, na.rm = TRUE),
        avg_length = rowMeans(length_matrix, na.rm = TRUE)
      ) %>%
        mutate(
          avg_sos = (avg_DER_sos + avg_TRS_sos) / 2,
          avg_eos = (avg_DER_eos + avg_TRS_eos) / 2
        )
      
      # Check if we got valid numbers (not all NA)
      if (all(is.na(data_grow_season$avg_sos))) {
        cat("[WARNING] All phenology metrics are NA for", site_name, "\n")
        cat("Saving with NA values for reference\n")
      }
      
      # SAVE INDIVIDUAL SITE FILE ONLY
      output_filename <- paste0(site_name, "_phenofit.csv")
      output_path <- file.path(save_folder, output_filename)
      
      write.csv(data_grow_season, output_path, row.names = FALSE)
      
      successful_sites[[site_name]] <- list(
        file = output_filename,
        years = nyears,
        rows = nrow(data_grow_season),
        manual = FALSE,
        all_na = all(is.na(data_grow_season$avg_sos))
      )
      
      if (all(is.na(data_grow_season$avg_sos))) {
        cat("[WARNING] RESCUE MODE - Processed but all NA for", site_name, "\n")
      } else {
        cat("[OK] RESCUE MODE - Successfully processed", site_name, "-", nyears, "years\n")
      }
      cat("[OK] Saved to:", output_path, "\n")
      
    }, error = function(e) {
      error_msg <- paste("Unexpected error:", e$message)
      cat("[ERROR]", error_msg, "- skipping", site_name, "\n")
      failed_sites[[site_name]] <<- error_msg
    })
  }
  
  # Summary
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("RESCUE MODE SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Input folder:", input_folder, "\n")
  cat("Output folder:", save_folder, "\n")
  cat("Total successful sites:", length(successful_sites), "\n")
  cat("Total failed sites:", length(failed_sites), "\n")
  
  if (length(successful_sites) > 0) {
    cat("\nSuccessfully processed:\n")
    for (site in names(successful_sites)) {
      if (successful_sites[[site]]$manual) {
        cat("  -", site, ":", successful_sites[[site]]$years, "years →", 
            successful_sites[[site]]$file, " (MANUALLY ENTERED VALUES)\n")
      } else if (successful_sites[[site]]$all_na) {
        cat("  -", site, ":", successful_sites[[site]]$years, "years →", 
            successful_sites[[site]]$file, " (ALL NA - insufficient data)\n")
      } else {
        cat("  -", site, ":", successful_sites[[site]]$years, "years →", 
            successful_sites[[site]]$file, "\n")
      }
    }
  }
  
  if (length(failed_sites) > 0) {
    cat("\nFailed sites:\n")
    for (site in names(failed_sites)) {
      cat("  -", site, ":", failed_sites[[site]], "\n")
    }
    
    # ONLY SAVE rescue_failed_sites.csv IF THERE ARE ACTUAL FAILED SITES
    failed_df <- data.frame(
      site_name = names(failed_sites),
      reason = unlist(failed_sites),
      stringsAsFactors = FALSE
    )
    failed_path <- file.path(save_folder, "rescue_failed_sites.csv")
    write.csv(failed_df, failed_path, row.names = FALSE)
    cat("\n[INFO] Failed sites list saved to:", failed_path, "\n")
  } else {
    cat("\n[INFO] No failed sites - all rescue sites processed successfully!\n")
  }
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PROCESSING COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  return(list(
    successful = successful_sites,
    failed = failed_sites
  ))
}

# ============================================================================
# RUN RESCUE MODE ON RESCUE SITES ONLY
# ============================================================================

# Define paths for rescue sites only
input_folder <- "M:\\Research\\WUE_CUE\\ameri_data\\ameri_ET_GPP\\rescue_sites"
save_folder <- "M:\\Research\\WUE_CUE\\ameri_data\\ameri_phenofit"

# Run the rescue processing
rescue_results <- phenofit_growing_season_rescue(input_folder, save_folder)

# Print final summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Rescue sites processed:", length(rescue_results$successful), "\n")
cat("Rescue sites failed:", length(rescue_results$failed), "\n")

if (length(rescue_results$successful) > 0) {
  cat("\nOutput files saved to:", save_folder, "\n")
  cat("Files:\n")
  for (site in names(rescue_results$successful)) {
    if (rescue_results$successful[[site]]$manual) {
      cat("  -", site, "_phenofit.csv (MANUALLY ENTERED VALUES)\n")
    } else if (rescue_results$successful[[site]]$all_na) {
      cat("  -", site, "_phenofit.csv (ALL NA - insufficient data)\n")
    } else {
      cat("  -", site, "_phenofit.csv\n")
    }
  }
}

if (length(rescue_results$failed) > 0) {
  cat("\nFailed sites:\n")
  for (site in names(rescue_results$failed)) {
    cat("  -", site, ":", rescue_results$failed[[site]], "\n")
  }
}