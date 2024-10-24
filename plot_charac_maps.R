library(raster)
library(sf)
library(ggplot2)
library(ggspatial)

setwd('C:\\Projetos\\bioflore\\gbs-uganda\\geo')
# Set the paths to your raster and vector files
plot_map <- function(raster_path,
                     vector_path,
                     crop,
                     title,
                     legend_title,
                     color,
                     legend_position,
                     discrete,
                     discrete_dict){
  
  
  # Load the raster data
  raster_data <- raster(raster_path)
  
  # Load the polygon vector data
  polygons <- st_read(vector_path)
  
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
    
    ggplot() +
      geom_tile(data = raster_df, aes(x = x, y = y, fill = discrete)) +
      scale_fill_viridis_d(legend_title, option = color) +
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
    ggplot() +
      geom_tile(data = raster_df, aes(x = x, y = y, fill = value)) +
      scale_fill_viridis_c(legend_title, option = color) +
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
  # Plot using ggplot2
  
}

plot_map(raster_path = 'elevation-dreinage/001.elevation.tif',
         vector_path = 'shp/restoration_sites.shp',
         crop = T,
         title = "Digital Elevation Model",
         legend_title = "Altitude (m)",
         color = "C",
         legend_position = 'right',
         discrete = F)

plot_map(raster_path = 'elevation-dreinage/slope.tif',
         vector_path = 'shp/restoration_sites.shp',
         crop = T,
         title = "Slope",
         legend_title = "degrees",
         color = "E",
         legend_position = 'right',
         discrete = F)

plot_map(raster_path = 'elevation-dreinage/dreinage.tif',
         vector_path = 'shp/restoration_sites.shp',
         crop = F,
         title = "Dreinage",
         legend_title = "",
         color = "A",
         legend_position = "none",
         discrete = F)

plot_map('cover_code\\uganda\\cover_code.tif\\download.cover_code.tif',
         'shp/restoration_sites.shp',
         crop = T,
         "CHM",
         "Height (m)",
         "B",
         "right",
         discrete = F)

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

plot_map('texture_0_20\\uganda\\texture_0_20.tif\\texture_class.texture_0_20.tif',
         'shp/restoration_sites.shp',
         crop = T,
         "Soil classification 0-20cm",
         "Class",
         "plasma",
         "right",
         discrete = T,
         discrete_dict = mapping)

plot_map('texture_20_50\\uganda\\texture_20_50.tif\\texture_class.texture_20_50.tif',
         'shp/restoration_sites.shp',
         crop = T,
         "Soil classification 20-50cm",
         "Class",
         "plasma",
         "right",
         discrete = T,
         discrete_dict = mapping)

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

plot_map('Map\\uganda\\Map.tif\\download.Map.tif',
         'shp/restoration_sites.shp',
         crop = T,
         "Land Cover Classification (2021)",
         "Class",
         "plasma",
         "right",
         discrete =T,
         discrete_dict = lulc)

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="C:\\Projetos\\bioflore\\gbs-uganda\\png")




