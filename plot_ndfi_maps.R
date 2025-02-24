#install.packages(c("sf", "raster", "ggplot2"))
library(sf)
library(raster)
library(ggplot2)
library(terra)
library(lubridate)
library(scales) # For the scales::squish function

get_utm_zone <- function(lon) {
  return(floor((lon + 180) / 6) + 1)
}

# Function to convert sf object to appropriate UTM
transform_to_utm <- function(sf_object) {
  # Ensure the sf object is in a geographic CRS
  if (st_crs(sf_object)$epsg != 4326) {
    sf_object <- st_transform(sf_object, crs = 4326)
  }
  
  # Calculate the centroid of the sf object
  centroid <- st_centroid(st_geometry(sf_object))
  centroid_coords <- st_coordinates(centroid)
  
  # Determine the UTM zone
  centroid_lon <- centroid_coords[1, 1]
  centroid_lat <- centroid_coords[1, 2]
  utm_zone <- get_utm_zone(centroid_lon)
  
  # Construct the UTM CRS
  utm_crs <- if (centroid_lat >= 0) {
    sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", utm_zone)
  } else {
    sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs +south", utm_zone)
  }
  
  # Reproject the sf object
  utm_sf_object <- st_transform(sf_object, crs = utm_crs)
  return(utm_sf_object)
}

setwd('C:\\Projetos\\bioflore\\bgci_II')
#2. Read in the shapefile containing the polygons.

shapefile_path <- "geo\\sites\\SITES+REF\\all_sites.shp"
polygons <- st_read(shapefile_path)
# label_list <- c('La_Pena', 'HUB_CERRADO_CECAP', 'HUB_CERRADO_NA_FLORESTA')  # Replace with your actual list of strings
# # Filter the rows where 'label' is in the list
# polygons <- polygons %>%
#   filter(name %in% label_list)%>%
#   mutate(name= gsub("_", "-", name))
# polygons <- transform_to_utm(polygons)

#3. List all the raster files in your directory.

raster_directory <- "geo\\NDFI"
raster_files <- list.files(raster_directory,
                           pattern = "\\.tif$",
                           full.names = TRUE,
                           recursive = T)

#4. Loop through each raster file and match it with the corresponding polygon.
for (raster_file in raster_files) {
  # raster_file <- raster_files[[2]]
  tryCatch({
    # Read the raster data
    raster_data <- rast(raster_file)
    
    non_na_count <- sum(!is.na(values(raster_data)))
    # If all values are NA, skip further processing for this raster
    if (non_na_count == 0) {
      message(paste("Skipping", raster_file, "as it contains only NA values."))
      next
    }
    # Extract the label and date from the raster file name
    label_date <- gsub(".tif", "", basename(raster_file))
    
    # Split the label and date from the filename
    parts <- strsplit(label_date, "_")[[1]]
    label <- parts[1]
    date <- parts[2]
    
    # Filter the polygons to find the corresponding polygon
    matching_polygon <- polygons[polygons$Name == label, ]
    
    matching_polygon <- st_transform(matching_polygon, crs = as.character(crs(raster_data)))
    
    # Convert raster to data frame and rename columns appropriately
    raster_df <- as.data.frame(terra::as.data.frame(raster_data, xy = TRUE), xy = TRUE)
    names(raster_df) <- c("lon", "lat", "NDFI")
    
    # Create a ggplot
    p <- ggplot() +
      geom_raster(data = raster_df, aes(x = lon, y = lat, fill = NDFI)) +
      geom_sf(data = matching_polygon, fill = NA, color = "red", size = 1) +
      theme_minimal() +
      scale_fill_gradientn(limits = c(-1, 1), colors = scales::viridis_pal()(7)) + 
      labs(title = paste("NDFI for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ymd")))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8))
    
    dir <- paste0("png/maps_ndfi/", label,"/")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(paste0(dir,"map_",label,"_", date, ".png"), plot = p)
    print(p)
  }, error = function(e) {
    message("Error processing ", raster_file, ": ", e$message)
  })
}
