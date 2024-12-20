library(raster)
library(rgrass)
library(terra)
library(stringr)
library(sf)
library(sp)

Sys.setlocale("LC_ALL", "C")
Sys.setenv(LANG = "en")
Sys.setenv('OSGEO4W_ROOT' = file.path('C:', 'OSGeo4W'))

# Load the DEM
dem_path <- "C:\\Projetos\\bioflore\\gbs-kenya\\geo\\elevation-dreinage\\001.elevation.tif"  # Specify your DEM file path
dem <- raster(dem_path)
plot(dem)
# Calculate the slope
slope <- terrain(dem, opt = "slope", unit = "degrees")
plot(slope)

breaks <- c(-Inf, 3, 8, 20, 45, Inf)  # Example breakpoints
discrete_slope_raster <- cut(slope, breaks=breaks)
plot(discrete_slope_raster)

labels <- c("Flat", "Gentle", "Moderate", "Steep", "Very Steep")
raster_levels <- data.frame(ID = 1:(length(breaks)-1), Labels = labels)
levels(discrete_slope_raster) <- raster_levels
# Save the slope to a new file
slope_output_path <- str_replace(dem_path,
                                 '001.elevation.tif', "slope.tif")
writeRaster(discrete_slope_raster, slope_output_path, format = "GTiff", overwrite = TRUE)


GRASS_INSTALLATION <- Sys.getenv("GRASS_INSTALLATION")

crs_string <- crs(dem, describe=T)
execGRASS("g.proj", flags="c", epsg=4326)

initGRASS(gisBase = "C:\\Program Files\\GRASS GIS 8.4",  # Path to your GRASS GIS folder
          gisDbase = "C:\\Projetos\\bioflore\\restoration_monitoring\\grassdb",                # Temporary database directory
          location = "kenya",                          # Location name
          mapset = "PERMANENT",                         # Mapset name
          SG = rast(dem_path),                          # DEM raster
          override = TRUE)

execGRASS("g.proj", flags = "p")
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
          threshold = 1000,          # Adjust this threshold according to your data
          flags = "overwrite")

# Extract streams using the flow accumulation
execGRASS("r.stream.extract",
          elevation = "dem_grass",
          accumulation = "flow_accumulation",
          threshold = 1000,       # Adjust this threshold based on your specific catchment size
          stream_raster = "streams",
          stream_vector = "streams_vec",
          flags = c("overwrite", "verbose"))

# execGRASS("v.out.ogr",
#           input = "streams_vec",
#           output = "C:/Projetos/bioflore/restoration_monitoring/grassdb/kenya/PERMANENT/drainage_lines.shp",
#           format = "ESRI_Shapefile",
#           flags = c("overwrite", "verbose", "quiet"))
# Export the stream vectors to a shapefile

execGRASS("r.out.gdal",
          input = "streams",
          output = "C:/Projetos/bioflore/restoration_monitoring/grassdb/kenya/PERMANENT/streams_check.tif",
          format = "GTiff",
          flags = "overwrite")

r <- raster("C:/Projetos/bioflore/restoration_monitoring/grassdb/kenya/PERMANENT/streams_check.tif")
r[] <- ifelse(!is.na(r[]), 1, NA)
plot(r)

writeRaster(r, "C:/Projetos/bioflore/restoration_monitoring/grassdb/kenya/PERMANENT/dreinage.tif", format="GTiff", overwrite=TRUE)
cat("Drainage lines saved as :/Projetos/bioflore/restoration_monitoring/grassdb/kenya/PERMANENT/\n")
