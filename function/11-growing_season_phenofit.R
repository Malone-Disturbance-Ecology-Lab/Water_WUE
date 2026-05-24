
# use separate version of this code for unique sites:
# unique sites US-EKN, US-Elm


#rm(list=ls(all=TRUE)) 

phenofit_growing_season <- function(input_folder, save_folder) {
  library(tidyverse)
  library(phenofit)
  library(data.table)
  library(lubridate)
  
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
  
  # Get list of all CSV files in input folder
  csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    stop("No CSV files found in input folder: ", input_folder)
  }
  
  cat("Found", length(csv_files), "CSV files to process\n")
  
  all_sites_results <- list()
  failed_sites <- list()  # Track failed sites and reasons
  successful_sites <- list()  # Track successful sites and file paths
  
  for (file_path in csv_files) {
    # Extract site name from filename
    # Example: US-A03_fill_gpp_et.csv -> US-A03
    filename <- basename(file_path)
    site_name <- gsub("_fill_gpp_et\\.csv$", "", filename)
    site_name <- gsub("\\.csv$", "", site_name)  # Fallback if pattern doesn't match
    
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("Processing site:", site_name, "\n")
    cat("File:", filename, "\n")
    
    # Try to process the site, catch any errors
    tryCatch({
      # Read the CSV file
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Check if required columns exist
      if (!"DateTime" %in% colnames(df)) {
        cat("[ERROR] Missing 'DateTime' column in", filename, "- skipping\n")
        failed_sites[[site_name]] <- "Missing 'DateTime' column"
        next
      }
      
      if (!"GPP" %in% colnames(df)) {
        cat("[ERROR] Missing 'GPP' column in", filename, "- skipping\n")
        failed_sites[[site_name]] <- "Missing 'GPP' column"
        next
      }
      
      # Print first few DateTime values to debug
      cat("First few DateTime values:", paste(head(df$DateTime), collapse = ", "), "\n")
      
      # Convert DateTime to proper format - SIMPLIFIED APPROACH
      # The DateTime is already in format "2014-01-01 00:00:00" or similar
      # Just convert to POSIXct directly
      df$DateTime <- as.POSIXct(df$DateTime, format = "%Y-%m-%d %H:%M:%S")
      
      # If that fails, try without seconds
      if (all(is.na(df$DateTime))) {
        cat("Trying format without seconds...\n")
        df$DateTime <- as.POSIXct(df$DateTime, format = "%Y-%m-%d %H:%M")
      }
      
      # If still failing, try with lubridate
      if (all(is.na(df$DateTime))) {
        cat("Trying lubridate...\n")
        df$DateTime <- ymd_hms(df$DateTime, quiet = TRUE)
        if (all(is.na(df$DateTime))) {
          df$DateTime <- ymd_hm(df$DateTime, quiet = TRUE)
        }
      }
      
      # Check if conversion worked
      if (all(is.na(df$DateTime))) {
        cat("[ERROR] Failed to parse DateTime column for", site_name, "- skipping\n")
        cat("Sample of DateTime values:", paste(head(df$DateTime), collapse = ", "), "\n")
        failed_sites[[site_name]] <- "Failed to parse DateTime column"
        next
      }
      
      # Remove rows with NA DateTime
      df <- df %>% filter(!is.na(DateTime))
      
      cat("Successfully parsed", nrow(df), "rows with valid DateTime\n")
      cat("DateTime range:", range(df$DateTime), "\n")
      
      # Aggregate half-hourly data to daily GPP
      sub_df <- df %>%
        mutate(Date = as.Date(DateTime)) %>%
        group_by(Date) %>%
        summarise(GPP = mean(GPP, na.rm = TRUE), .groups = "drop") %>%
        mutate(w = 1.0)  # Uniform weights
      
      # Remove rows with NA GPP
      sub_df <- sub_df %>% filter(!is.na(GPP))
      
      # Check if enough data
      if (nrow(sub_df) < 50) {
        cat("[WARNING] Insufficient data for", site_name, "- only", nrow(sub_df), "days, skipping\n")
        failed_sites[[site_name]] <- paste("Insufficient data: only", nrow(sub_df), "days")
        next
      }
      
      # Format dates for printing properly
      min_date <- format(min(sub_df$Date), "%Y-%m-%d")
      max_date <- format(max(sub_df$Date), "%Y-%m-%d")
      
      cat("Total days with GPP data:", nrow(sub_df), "\n")
      cat("Date range:", min_date, "to", max_date, "\n")
      cat("Years present:", paste(unique(format(sub_df$Date, "%Y")), collapse = ", "), "\n")
      
      # Parameters for phenofit
      nptperyear <- 365
      wFUN <- wTSM
      wmin <- 0.2
      minExtendMonth <- 0.5
      maxExtendMonth <- 2
      minPercValid <- 0
      south_hemisphere <- FALSE
      
      # Check input data
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
      if (site_name %in% c("US-Atq", "US-EvM", "US-StS")) {
        lambda_value <- ifelse(site_name == "US-Atq", 100, 1000)
        cat("Using special lambda =", lambda_value, "for site", site_name, "\n")
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
        cat("Using default lambda = 50000 for site", site_name, "\n")
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
      
      # Curve fitting using multiple methods
      fit <- curvefits(INPUT, brks,
                       list(methods = c("Beck", "Elmore", "Gu", "AG", "Zhang"),
                            wFUN = wFUN, iters = 5,
                            wmin = wmin,
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
        failed_sites[[site_name]] <- "No phenology data extracted"
        next
      }
      
      # Extract data frames for each method and average duplicate years
      doy_Beck <- tryCatch(average_duplicate_years(l_pheno$doy$Beck), error = function(e) NULL)
      doy_Elmore <- tryCatch(average_duplicate_years(l_pheno$doy$Elmore), error = function(e) NULL)
      doy_Gu <- tryCatch(average_duplicate_years(l_pheno$doy$Gu), error = function(e) NULL)
      doy_AG <- tryCatch(average_duplicate_years(l_pheno$doy$AG), error = function(e) NULL)
      doy_Zhang <- tryCatch(average_duplicate_years(l_pheno$doy$Zhang), error = function(e) NULL)
      
      # Check which methods have valid data
      valid_methods <- list(
        Beck = doy_Beck,
        Elmore = doy_Elmore,
        Gu = doy_Gu,
        AG = doy_AG,
        Zhang = doy_Zhang
      )
      
      valid_methods <- valid_methods[!sapply(valid_methods, is.null)]
      
      if (length(valid_methods) < 2) {
        cat("[WARNING] Insufficient valid methods for", site_name, "- only", length(valid_methods), "methods, skipping\n")
        failed_sites[[site_name]] <- paste("Insufficient methods:", length(valid_methods), "methods only")
        next
      }
      
      # Find common years across all valid methods
      common_years <- Reduce(intersect, lapply(valid_methods, function(x) x$year))
      
      if (length(common_years) == 0) {
        cat("[WARNING] No common years across methods for", site_name, "- skipping\n")
        failed_sites[[site_name]] <- "No common years across methods"
        next
      }
      
      # Filter each method to common years
      for (method_name in names(valid_methods)) {
        valid_methods[[method_name]] <- filter(valid_methods[[method_name]], year %in% common_years)
      }
      
      nyears <- length(common_years)
      cat("Found", nyears, "common years across methods\n")
      
      # Build averaged phenology results
      data_grow_season <- tibble(
        site_name = rep(site_name, nyears),
        year = common_years,
        avg_DER_sos = rowMeans(cbind(
          valid_methods$Beck$DER.sos, 
          valid_methods$Elmore$DER.sos, 
          valid_methods$Gu$DER.sos, 
          valid_methods$AG$DER.sos, 
          valid_methods$Zhang$DER.sos
        ), na.rm = TRUE),
        avg_DER_eos = rowMeans(cbind(
          valid_methods$Beck$DER.eos, 
          valid_methods$Elmore$DER.eos, 
          valid_methods$Gu$DER.eos, 
          valid_methods$AG$DER.eos, 
          valid_methods$Zhang$DER.eos
        ), na.rm = TRUE),
        avg_TRS_sos = rowMeans(cbind(
          valid_methods$Beck$TRS2.5.sos, 
          valid_methods$Elmore$TRS2.5.sos, 
          valid_methods$Gu$TRS2.5.sos, 
          valid_methods$AG$TRS2.5.sos, 
          valid_methods$Zhang$TRS2.5.sos
        ), na.rm = TRUE),
        avg_TRS_eos = rowMeans(cbind(
          valid_methods$Beck$TRS2.5.eos, 
          valid_methods$Elmore$TRS2.5.eos, 
          valid_methods$Gu$TRS2.5.eos, 
          valid_methods$AG$TRS2.5.eos, 
          valid_methods$Zhang$TRS2.5.eos
        ), na.rm = TRUE),
        avg_length = rowMeans(cbind(
          valid_methods$Beck$DER.eos - valid_methods$Beck$DER.sos,
          valid_methods$Elmore$DER.eos - valid_methods$Elmore$DER.sos,
          valid_methods$Gu$DER.eos - valid_methods$Gu$DER.sos,
          valid_methods$AG$DER.eos - valid_methods$AG$DER.sos,
          valid_methods$Zhang$DER.eos - valid_methods$Zhang$DER.sos
        ), na.rm = TRUE)
      ) %>%
        mutate(
          avg_sos = (avg_DER_sos + avg_TRS_sos) / 2,
          avg_eos = (avg_DER_eos + avg_TRS_eos) / 2
        )
      
      # ===================================================================
      # SAVE INDIVIDUAL SITE FILE
      # ===================================================================
      # Create output filename: [site_name]_phenofit.csv
      output_filename <- paste0(site_name, "_phenofit.csv")
      output_path <- file.path(save_folder, output_filename)
      
      # Save the phenology results for this site
      write.csv(data_grow_season, output_path, row.names = FALSE)
      
      # Store successful site info
      successful_sites[[site_name]] <- list(
        file = output_filename,
        years = nyears,
        rows = nrow(data_grow_season),
        date_range = paste(min_date, "to", max_date)
      )
      
      cat("[OK] Successfully processed", site_name, "-", nyears, "years of phenology data\n")
      cat("[OK] Saved to:", output_path, "\n")
      
      # Also store in all_sites_results for potential combined output later if needed
      all_sites_results[[site_name]] <- data_grow_season
      
    }, error = function(e) {
      # Catch any unexpected errors during processing
      error_msg <- paste("Unexpected error:", e$message)
      cat("[ERROR]", error_msg, "- skipping", site_name, "\n")
      failed_sites[[site_name]] <- error_msg
    })
  }
  
  # Print summary of successful sites
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SUCCESSFUL SITES SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  if (length(successful_sites) == 0) {
    cat("No sites were processed successfully!\n")
  } else {
    cat("Total successful sites:", length(successful_sites), "\n")
    cat("\nSuccessful sites and their files:\n")
    for (site in names(successful_sites)) {
      cat("  -", site, ":", successful_sites[[site]]$file, 
          "(", successful_sites[[site]]$years, "years,", 
          successful_sites[[site]]$rows, "rows)\n")
      cat("    Date range:", successful_sites[[site]]$date_range, "\n")
    }
  }
  
  # Print summary of failed sites
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("FAILED SITES SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  if (length(failed_sites) == 0) {
    cat("All sites processed successfully!\n")
  } else {
    cat("Total failed sites:", length(failed_sites), "\n")
    cat("\nFailed sites and reasons:\n")
    for (site in names(failed_sites)) {
      cat("  -", site, ":", failed_sites[[site]], "\n")
    }
  }
  
  # Save master summary file
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SAVING SUMMARY FILES\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Save successful sites summary
  if (length(successful_sites) > 0) {
    successful_df <- data.frame(
      site_name = names(successful_sites),
      file_name = sapply(successful_sites, function(x) x$file),
      years = sapply(successful_sites, function(x) x$years),
      rows = sapply(successful_sites, function(x) x$rows),
      date_range = sapply(successful_sites, function(x) x$date_range),
      stringsAsFactors = FALSE
    )
    successful_path <- file.path(save_folder, "successful_sites.csv")
    write.csv(successful_df, successful_path, row.names = FALSE)
    cat("[OK] Successful sites list saved to:", successful_path, "\n")
  }
  
  # Save failed sites to a file for reference
  if (length(failed_sites) > 0) {
    failed_df <- data.frame(
      site_name = names(failed_sites),
      reason = unlist(failed_sites),
      stringsAsFactors = FALSE
    )
    failed_path <- file.path(save_folder, "failed_sites.csv")
    write.csv(failed_df, failed_path, row.names = FALSE)
    cat("[OK] Failed sites list saved to:", failed_path, "\n")
  }
  
  # Optionally save combined results if needed
  if (length(all_sites_results) > 0) {
    final_df <- bind_rows(all_sites_results)
    combined_path <- file.path(save_folder, "all_sites_phenofit_combined.csv")
    write.csv(final_df, combined_path, row.names = FALSE)
    cat("[OK] Combined results saved to:", combined_path, "\n")
  }
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PROCESSING COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Total sites processed successfully:", length(successful_sites), "\n")
  cat("Total sites failed:", length(failed_sites), "\n")
  cat("Output folder:", save_folder, "\n")
  
  return(list(
    successful = successful_sites,
    failed = failed_sites,
    results = all_sites_results
  ))
}

##########################################################################
# Define the paths
input_folder <- "M:\\Research\\WUE_CUE\\ameri_data\\ameri_ET_GPP"
save_folder <- "M:\\Research\\WUE_CUE\\ameri_data\\ameri_phenofit"

# Run the function
phenology_results <- phenofit_growing_season(input_folder, save_folder)

# Print final summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Input folder:", input_folder, "\n")
cat("Output folder:", save_folder, "\n")
cat("Number of sites processed:", length(phenology_results$successful), "\n")
cat("Number of sites failed:", length(phenology_results$failed), "\n")

if (length(phenology_results$successful) > 0) {
  cat("\nOutput files created:\n")
  for (site in names(phenology_results$successful)) {
    cat("  -", site, "_phenofit.csv\n")
  }
  cat("  - successful_sites.csv\n")
  cat("  - all_sites_phenofit_combined.csv\n")
}

if (length(phenology_results$failed) > 0) {
  cat("\nFailed sites saved to: failed_sites.csv\n")
}