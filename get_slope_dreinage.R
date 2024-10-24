library(raster)
library(rgrass)
library(terra)
library(stringr)
library(sf)

Sys.setlocale("LC_ALL", "C")
Sys.setenv(LANG = "en")

# Load the DEM
dem_path <- "/tmp/elevation/uganda/elevation.tif/001.elevation.tif"  # Specify your DEM file path
dem <- raster(dem_path)
plot(dem)
# Calculate the slope
slope <- terrain(dem, opt = "slope", unit = "degrees")

plot(slope)
# Save the slope to a new file
slope_output_path <- str_replace(dem_path,
                                 '001.elevation.tif', "slope.tif")
writeRaster(slope, slope_output_path, format = "GTiff", overwrite = TRUE)


initGRASS(gisBase = "C:\\Program Files\\GRASS GIS 8.4",  # Path to your GRASS GIS folder
          gisDbase = "C:\\tmp\\grassdb",                # Temporary database directory
          location = "uganda",                          # Location name
          mapset = "PERMANENT",                         # Mapset name
          SG = rast(dem_path),                          # DEM raster
          override = TRUE)

execGRASS("g.region", flags = "p")

# Load DEM into GRASS
execGRASS("r.in.gdal",
          input = dem_path,
          output = "dem_grass",
          flags = "overwrite")

# Set the GRASS region to match the DEM
execGRASS("g.region", raster = "dem_grass")

# Calculate Flow Accumulation and Direction
execGRASS("r.watershed",
          elevation = "dem_grass",
          accumulation = "flow_accumulation",
          drainage = "flow_direction",
          threshold = 100,          # Adjust this threshold according to your data
          flags = "overwrite")

# Extract streams using the flow accumulation
execGRASS("r.stream.extract",
          elevation = "dem_grass",
          accumulation = "flow_accumulation",
          threshold = 100,       # Adjust this threshold based on your specific catchment size
          stream_raster = "streams",
          stream_vector = "streams_vec",
          flags = c("overwrite", "verbose"))

# execGRASS("v.out.ogr",
#           input = "streams_vec",
#           output = "C:/tmp/grassdb/uganda/drainage_lines.shp",
#           format = "ESRI_Shapefile",
#           flags = c("overwrite", "verbose", "quiet"))
# Export the stream vectors to a shapefile

execGRASS("r.out.gdal",
          input = "streams",
          output = "C:/tmp/grassdb/uganda/streams_check.tif",
          format = "GTiff",
          flags = "overwrite")

r <- raster("C:/tmp/grassdb/uganda/streams_check.tif")
r[] <- ifelse(!is.na(r[]), 1, NA)
plot(r, main="Converted Raster", col = c(NA, "blue"), legend = FALSE)
writeRaster(r, "C:/tmp/grassdb/uganda/dreinage.tif", format="GTiff", overwrite=TRUE)

cat("Drainage lines saved as C:/tmp/grassdb/uganda/streams_check.tif\n")