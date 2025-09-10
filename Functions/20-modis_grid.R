library(sf)
library(fs)
library(readr)
library(zip)

### for following sites need to move point a bit

# US-LA3, US-LA2, US-LA1 ET
# US-KS3, US-KS4 ET
#US-TaS GPP


create_modis_grids <- function(input_csv, output_folder) {
  # Read input CSV
  df <- read_csv(input_csv)
  
  # Define MODIS Sinusoidal projection
  modis_crs <- "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"
  
  for (i in 1:nrow(df)) {
    site <- df$site_name[i]
    lat <- df$lat[i]
    lon <- df$long[i]
    
    # Create WGS84 point
    pt_wgs <- st_sfc(st_point(c(lon, lat)), crs = 4326)
    
    # Transform to MODIS projection
    pt_modis <- st_transform(pt_wgs, crs = modis_crs)
    
    # 750 m buffer around point
    bbox_area <- st_buffer(pt_modis, dist = 750)
    
    # Create 3x3 grid of 500m cells
    grid_raw <- st_make_grid(bbox_area, cellsize = 500, what = "polygons")
    grid_sf <- st_sf(id = 1:length(grid_raw), geometry = grid_raw)
    grid_sf$id <- matrix(1:9, nrow = 3, byrow = FALSE)[, 3:1] |> as.vector()
    
    # Reproject to WGS84 for export
    grid_wgs84 <- st_transform(grid_sf, crs = 4326)
    
    # Define folder names
    site_folder <- file.path(tempdir(), site)
    dir_create(site_folder)
    
    # Define shapefile path
    shp_name <- file.path(site_folder, paste0("modis_grid_", site, ".shp"))
    
    # Write shapefile
    st_write(grid_wgs84, shp_name, delete_layer = TRUE, quiet = TRUE)
    
    # Zip the folder
    zip_path <- file.path(output_folder, paste0(site, "_shp.zip"))
    zip(zipfile = zip_path, files = dir_ls(site_folder, glob = "*.*"), mode = "cherry-pick")
  }
  
  message("All shapefiles generated and zipped.")
}


create_modis_grids(
  input_csv = "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/modis_data/modis_lat_long.csv",
  output_folder = "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/modis_data/shape_data"
)

#############################################################################

