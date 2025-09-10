
library(daymetr)
library(readr)
library(dplyr)

get_daymet_precip_for_sites <- function(input_folder,
                                        csv_file = "site_lat_long.csv",
                                        save_folder = "C:/temp/precip_daymet") {
  # Create output folder if needed
  if (!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
  
  # Read input CSV
  input_path <- file.path(input_folder, csv_file)
  sites <- read_csv(input_path, show_col_types = FALSE)
  
  # Check required columns
  required_cols <- c("site_name", "lat", "long", "start_date", "end_date")
  if (!all(required_cols %in% names(sites))) {
    stop("Missing required columns: site_name, lat, long, start_date, end_date")
  }
  
  for (i in seq_len(nrow(sites))) {
    # Get and explicitly preserve the site name
    site_name <- sites$site_name[i]
    lat <- sites$lat[i]
    lon <- sites$long[i]
    
    # Convert year-only to numeric
    start_year <- as.numeric(sites$start_date[i])
    end_year <- as.numeric(sites$end_date[i])
    
    message(paste("Downloading Daymet for", site_name, "(", lat, lon, ") from", start_year, "to", end_year))
    
    tryCatch({
      result <- download_daymet(
        lat = lat,
        lon = lon,
        start = start_year,
        end = end_year,
        internal = TRUE
      )
      
      df <- result$data %>%
        select(year, yday, prcp..mm.day.) %>%
        mutate(site = site_name,  # Use the original site name
               date = as.Date(yday - 1, origin = paste0(year, "-01-01"))) %>%
        rename(precip_mm = prcp..mm.day.) %>%
        select(site, date, precip_mm)
      
      # Create filename with explicit case preservation
      filename <- paste0(site_name, ".csv")
      filepath <- file.path(save_folder, filename)
      
      # Debug: show what we're about to save
      message(paste("Saving to:", filepath))
      
      write_csv(df, filepath)
      
    }, error = function(e) {
      message(paste("  Skipped", site_name, "-", conditionMessage(e)))
    })
  }
}

# Run the function
get_daymet_precip_for_sites(
  input_folder = "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/data_products/info",
  csv_file = "site_lat_long.csv",
  save_folder = "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/drivers/ameri_drivers/precip"
)