library(raster)
library(sf)   # If using sf, use `library(sf)` instead

setwd("C:\\Projetos\\bioflore\\gbs-madagascar\\geo")

polygon_shapefile <- st_read("possible_forest_restauration_areas_WGS84_f38S\\possible_forest_restauration_areas_WGS84_f38S.shp")
polygon_shapefile <- st_transform(polygon_shapefile, crs = 32738)
paths <- grep(list.files('limite-sentinel2', full.name=T),
                             pattern = "limite_", value= T)

for (i in paths){
  # i <- paths[[1]]
  # Load the raster brick (example file path)
  raster_brick <- brick(i)
  # Load the polygon shapefile (example file path)
  
  # Ensure the CRS match between the raster and polygon
  raster_brick <- projectRaster(raster_brick, crs = crs(polygon_shapefile))
  
  # Masking the raster brick with the polygon
  masked_raster_brick <- mask(raster_brick, polygon_shapefile)
  
  # Cropping the raster brick with the polygon's extent
  cropped_raster_brick <- crop(masked_raster_brick, extent(polygon_shapefile))
  
  # Save the result to disk
  output_file_path <- gsub('.tif','_recover.tif', i)
  writeRaster(cropped_raster_brick, filename = output_file_path, format = "GTiff", overwrite = TRUE)
}
