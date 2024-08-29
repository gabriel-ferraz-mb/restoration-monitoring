#install.packages(c("sf", "raster", "ggplot2"))
library(sf)
library(raster)
library(ggplot2)
library(terra)
library(lubridate)
library(scales) # For the scales::squish function

setwd('C:\\Projetos\\bioflore\\ecosia')
#2. Read in the shapefile containing the polygons.

shapefile_path <- "geo\\merge_ipe_espi.shp"
polygons <- st_read(shapefile_path)
polygons <- st_transform(polygons, crs = 4326)
#3. List all the raster files in your directory.

raster_directory <- "NDFI"
raster_files <- list.files(raster_directory, pattern = "\\.tif$", full.names = TRUE)

#4. Loop through each raster file and match it with the corresponding polygon.
for (raster_file in raster_files) {
 # raster_file = raster_files[1]
  # Extract the label and date from the raster file name (assuming the filename format is consistent)
  label_date <- gsub(".tif", "", basename(raster_file))
  
  # Split the label and date from the filename or use a regex if needed
  parts <- strsplit(label_date, "_")[[1]]  # Assuming filename always has a separator like _
  label <- paste0(parts[1],"_",parts[2])
  date <- parts[3]
  
  # Filter the polygons to find the corresponding polygon
  matching_polygon <- polygons[polygons$Name == label, ]
  
  if (nrow(matching_polygon) > 0) {
    # Read the raster data
    raster_data <- rast(raster_file)
    raster_data <- terra::project(raster_data, "EPSG:4326")
    
    # Convert raster to data frame and rename columns appropriately
    raster_df <- as.data.frame(terra::as.data.frame(raster_data, xy = TRUE), xy = TRUE)
    names(raster_df) <- c("x", "y", "value")
    
    # Create a ggplot
    p <- ggplot() +
      geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
      geom_sf(data = matching_polygon, fill = NA, color = "red", size = 1) +
      theme_minimal() +
      scale_fill_gradientn(limits = c(-1, 1), colors = scales::viridis_pal()(7)) + 
      labs(title = paste("NDFI for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ymd")))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8))
    
    # Save or display the plot
    ggsave(paste0("maps_ndfi\\map_", label, "_", date, ".png"), plot = p)
    print(p)
  } else {
    message("No matching polygon found for ", label, " on ", date)
  }
}
