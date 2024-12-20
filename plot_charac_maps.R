library(raster)
library(terra)
library(sf)
library(ggplot2)
library(ggspatial)
library(glue)
library(patchwork)

setwd('C:\\Projetos\\bioflore\\gbs-kenya\\geo')
# Set the paths to your raster and vector files
plot_map <- function(raster_path,
                     polygons,
                     crop,
                     title,
                     legend_title,
                     color,
                     legend_position,
                     discrete,
                     discrete_dict,
                     color_mapping,
                     hide_ticks = F){
  
  # Load the raster data
  raster_data <- raster(raster_path)
  
  # Load the polygon vector data
  # polygons <- st_read(vector_path)
  
  # Ensure the CRSs match
  polygons <- st_transform(st_zm(polygons), crs = crs(raster_data))
  
  if (crop){
    raster_data <- mask(crop(raster_data, extent(polygons)), polygons)
  }
  # Crop the raster using the polygon
  
  raster_df <- as.data.frame(rasterToPoints(raster_data))
  names(raster_df) <- c("x", "y", "value")
  
  if (discrete) {
    raster_df$value <- as.character(raster_df$value)
    raster_df$discrete <- discrete_dict[raster_df$value]
    
   p <-  ggplot() +
      geom_tile(data = raster_df, aes(x = x, y = y, fill = discrete)) +
      scale_fill_manual(values = color_mapping, name = legend_title) +
      geom_sf(data = polygons, fill = NA, color = "red", size = 0.5) +
      coord_sf(xlim = c(st_bbox(polygons)["xmin"], st_bbox(polygons)["xmax"]),
               ylim = c(st_bbox(polygons)["ymin"], st_bbox(polygons)["ymax"]),
               expand = FALSE) +
      ggtitle(title) +
      labs(x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(legend.position = legend_position,
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  else{
    p <-   ggplot() +
      geom_tile(data = raster_df, aes(x = x, y = y, fill = value)) +
      scale_fill_viridis_c(legend_title, option = color) +
      geom_sf(data = polygons, fill = NA, color = "red", size = 0.5) +
      coord_sf(xlim = c(st_bbox(polygons)["xmin"], st_bbox(polygons)["xmax"]),
               ylim = c(st_bbox(polygons)["ymin"], st_bbox(polygons)["ymax"]),
               expand = FALSE) +
      ggtitle(title) +
      labs(x = "Longitude", y = "Latitude") +
      theme_minimal()
    # Conditionally alter the theme based on hide_ticks
    if (hide_ticks) {
      p <- p + theme(
        legend.position = legend_position,
        axis.ticks.x = element_blank(),  # Hide x axis ticks
        axis.ticks.y = element_blank(),  # Hide y axis ticks
        axis.text.x = element_blank(),   # Hide x axis text
        axis.text.y = element_blank()    # Hide y axis text
      )
    } else {
      p <- p + theme(legend.position = legend_position,
                     axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }
  # Plot using ggplot2
  return(p)
}

plot_vector <- function(data_path,
                        polygons,
                        crop,
                        title,
                        legend_title,
                        color,
                        legend_position) {  # Removed extra comma
  
  # Load the vector data
  vector_data <- st_read(data_path)
  
  # Load the polygon vector data
  # polygons <- st_read(vector_path)
  
  # Ensure the CRSs match
  polygons <- st_transform(st_zm(polygons), crs = st_crs(vector_data))  # Use st_crs
  
  if (crop){
    vector_data <- st_intersection(polygons, vector_data)  # Correct variable name
  }
  
  ggplot() +
    geom_sf(data = vector_data, color = "blue", size = 1) +
    geom_sf(data = polygons, fill = NA, color = "red", size = 0.5) +
    coord_sf(xlim = c(st_bbox(polygons)["xmin"], st_bbox(polygons)["xmax"]),
             ylim = c(st_bbox(polygons)["ymin"], st_bbox(polygons)["ymax"]),
             expand = FALSE) +
    ggtitle(title) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = legend_position,
          axis.text.x = element_text(angle = 45, hjust = 1))
}

vector <- st_read('shp/mount_kenya.shp')

# Plot Elevation
plot_map(raster_path = 'elevation-dreinage/001.elevation.tif',
         polygons = vector,
         crop = T,
         title = "Digital Elevation Model",
         legend_title = "Altitude (m)",
         color = "C",
         legend_position = 'right',
         discrete = F)

# Plot Slope
slope_labels <- c(
 '1' = "Flat",
 '2' = "Gentle",
 '3' = "Moderate",
 '4' = "Steep",
 '5'= "Very Steep")

slope_mapping <- c(
  "Flat" = ,
  "Gentle" = ,
  "Moderate" = ,
  "Steep" = ,
  "Very Steep" = )

plot_map(raster_path = 'elevation-dreinage/slope.tif',
         polygons = vector,
         crop = T,
         title = "Slope",
         legend_title = "degrees",
         color = "E",
         legend_position = 'right',
         discrete = T,
         discrete_dict = slope_labels)

# Plot Dreinage
plot_vector(data_path = 'elevation-dreinage/dreinage_lines.shp',
            polygons = vector,
            crop = T,
            title = "Dreinage",
            legend_title = "",
            legend_position = "none")

# Plot Soil Class
mapping <- c(
  '1' = "Clay",
  '2' = "Silty Clay",
  '3' = "Sandy Clay",
  '4' = "Clay Loam",
  '5' = "Silty Clay Loam",
  '6' = "Sandy Clay Loam",
  '7' = "Loam",
  '8' = "Silt Loam",
  '9' = "Sandy Loam",
  '10' = "Silt",
  '11' = "Loamy Sand",
  '12' = "Sand"
)

plot_map('texture_0_20\\texture_class.texture_0_20.tif',
         polygons = vector,
         crop = T,
         "Soil classification 0-20cm",
         "Class",
         "plasma",
         "right",
         discrete = T,
         discrete_dict = mapping)

plot_map('texture_20_50\\texture_class.texture_20_50.tif',
         polygons = vector,
         crop = T,
         "Soil classification 20-50cm",
         "Class",
         "plasma",
         "right",
         discrete = T,
         discrete_dict = mapping)

# Plot Land Cover
lulc <- c(
  '10' = 'Tree cover',
  '20'='Shrubland',
  '30'='Grassland',
  '40'='Cropland',
  '50'='Built-up',
  '60'='Bare / sparse vegetation',
  '70'='Snow and ice',
  '80'='Permanent water bodies',
  '90'='Herbaceous wetland',
  '95'='Mangroves',
  '100'='Moss and lichen'
)

lulc_color <- c(
  'Tree cover'='#006400',
  'Shrubland'='#ffbb22',
  'Grassland'='#ffff4c',
  'Cropland'='#f096ff',
  'Built-up'='#fa0000',
  'Bare / sparse vegetation'='#b4b4b4',
  'Snow and ice'='#f0f0f0',
  'Permanent water bodies'='#0064c8',
  'Herbaceous wetland'='#0096a0',
  'Mangroves'='#00cf75',
  'Moss and lichen'='#fae6a0'
)	

root_directory <- "Map"
raster_files <- list.files(root_directory, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
raster_list <- lapply(raster_files, rast)
merged_raster <- do.call(mosaic, raster_list)
writeRaster(merged_raster, "Map/lulc_merged_raster.tif", overwrite = TRUE)


plot_map(raster_path = 'Map/lulc_merged_raster.tif',
         polygons = vector,
         crop = T,
         title = "Land Cover Classification (2021)",
         legend_title = "Class",
         legend_position = "right",
         discrete = T,
         discrete_dict = lulc,
         color_mapping = lulc_color)

# Plot CHM
restoration_sites = st_read('shp/restoration_sites.shp')
individual_polygons <- st_cast(restoration_sites, "POLYGON")

root_directory <- "cover_code"
raster_files <- list.files(root_directory, pattern = "\\.tif$",
                           full.names = TRUE, recursive = TRUE)
raster_list <- lapply(raster_files, rast)
s <- sprc(raster_list)
merged_raster <- merge(s)
writeRaster(merged_raster, "cover_code/chm_merged_raster.tif", overwrite = TRUE)

plotList = list() 
# Iterate over each polygon and plot
for (i in seq_len(nrow(individual_polygons))) {
  single_polygon <- individual_polygons[i, ]
  
  t <- glue("CHM_{single_polygon$label}")
  
  p <- plot_map(raster_path = "cover_code/chm_merged_raster.tif",
           polygons = single_polygon,
           crop = T,
           title = t,
           legend_title = "Height (m)",
           color = "B",
           legend_position = "none",
           discrete = F,
           hide_ticks = T)
  plotList[[i]] <- p  
  }
combined_plot <- wrap_plots(plotList) 
print(combined_plot)

# Save Plots
lots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="C:\\Projetos\\bioflore\\gbs-kenya\\png")




