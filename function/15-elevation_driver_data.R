# Load required packages
library(elevatr)
library(sf)
library(dplyr)
library(readr)
library(purrr)

################################################################################

## code to get elevation data
get_site_elevations <- function(input_folder) {
  library(elevatr)
  library(sf)
  library(readr)
  library(dplyr)
  
  input_file <- file.path(input_folder, "site_lat_long.csv")
  site_df <- read_csv(input_file, col_types = cols())
  
  if (!all(c("lat", "long") %in% colnames(site_df))) {
    stop("Input CSV must contain 'lat' and 'long' columns.")
  }
  
  site_sf <- st_as_sf(site_df, coords = c("long", "lat"), crs = 4326, remove = FALSE)
  
  # Use AWS terrain tiles — more robust in remote/coastal areas
  site_elev_sf <- get_elev_point(locations = site_sf, src = "aws", units = "meters")
  
  # Attach elevation
  site_df$elevation_m <- site_elev_sf$elevation
  
  output_file <- file.path(input_folder, "site_lat_long_with_elevation.csv")
  write_csv(site_df, output_file)
  
  message("Elevation data saved to: ", output_file)
}


get_site_elevations("\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\data_products\\info")


# Example usage:
input_folder <- "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\data_products\\info"


###########################################################################
## code to get climate data














