# install.packages("geojsonio")
library(landscapemetrics)
library(terra)
library(sf)
library(dplyr)

##################################

setwd('C:\\Projetos\\bioflore\\gbs-uganda\\geo')

landscape <- terra::rast("Map\\uganda\\Map.tif\\download.Map.tif")
aoi <- polygons <- st_transform(st_zm(st_read('shp/restoration_sites.shp')),
                                crs = crs(landscape))

abreviations <- lsm_abbreviations_names
options_landscapemetrics(to_disk = NULL)
check_landscape(landscape, verbose = TRUE)
list <- list_lsm()

get_utm_zone <- function(lon) {
  return(floor((lon + 180) / 6) + 1)
}

# Get the central longitude of your raster extent's centroid
raster_extent <- terra::ext(landscape)
centroid_lon <- (raster_extent$xmax + raster_extent$xmin) / 2
centroid_lat <- (raster_extent$ymax + raster_extent$ymin) / 2

utm_zone <- get_utm_zone(centroid_lon)

# Determine the UTM zone
utm_crs <- if (centroid_lat >= 0) {
  sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", utm_zone)
} else {
  sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs +south", utm_zone)
}

utm_raster <- project(landscape, utm_crs)
plot(utm_raster)

r <- calculate_lsm(
  utm_raster,
  # level = 'class',
  # metric = c('area_mn','ca','np','te','ed','shape_mn','tca','enn_mn'),
  name = NULL,
  type = NULL,
  what =c('lsm_c_area_mn', 'lsm_c_enn_mn', 'lsm_c_shape_mn'),
  directions = 8,
  count_boundary = FALSE,
  consider_boundary = FALSE,
  edge_depth = 3,
  cell_center = FALSE,
  classes_max = NULL,
  neighbourhood = 4,
  ordered = TRUE,
  base = "log2",
  full_name = FALSE,
  verbose = TRUE,
  progress = TRUE
)

r_metrics <- calculate_lsm(
  utm_raster,
  level = 'class',
  metric = c('ca','np','te','ed','tca'),
  name = NULL,
  type = NULL,
  # what =c('lsm_c_area_mn', 'lsm_c_enn_mn', 'lsm_c_shape_mn'),
  directions = 8,
  count_boundary = FALSE,
  consider_boundary = FALSE,
  edge_depth = 3,
  cell_center = FALSE,
  classes_max = NULL,
  neighbourhood = 4,
  ordered = TRUE,
  base = "log2",
  full_name = FALSE,
  verbose = TRUE,
  progress = TRUE
)

result = rbind(r,r_metrics)


show_correlation(
  result,
  method = "pearson",
  diag = TRUE,
  labels = FALSE,
  vjust = 0,
  text_size = 15
)

result <- result %>%
  filter(class == "10") %>%
  mutate_if(is.numeric, round, 5)

result <- merge(result, abreviations[, c("metric", "name")], by = "metric", all = F)
result <- result[!duplicated(result), ]

patches <- show_patches(
  landscape,
  class = "10",
  directions = 8,
  labels = F,
  nrow = NULL,
  ncol = NULL
)

cores <- show_cores(
  landscape,
  directions = 8,
  class = "10",
  labels = FALSE,
  nrow = NULL,
  ncol = NULL,
  consider_boundary = F,
  edge_depth = 3
)

options(scipen = 999)
write.table(result , file = "..\\spreadsheets\\lsm.csv")


