# Curso: Engenharia Cartogr?fica - Instituto Militar de Engenharia (IME)
# Disciplina: 06032- PROCESSAMENTO DIGITAL DE IMAGENS II
# Aula 2: Compress?o de imagens - Data: 10/03/2020
# Prof: Dr. Matheus Pinheiro Ferreira - matheus@ime.eb.br

# Primeiramente devemos carregar alguns pacotes, caso os pacotes n?o estejam instalados rode o comando: install.packages()
library(raster) # install.packages('raster')
# library(rgdal)  # install.packages('rgdal')
library(RStoolbox)
library(ggplot2)# install.packages('RStoolbox')
library(data.table)
library(stringr)
library(dplyr)
library(splines)
library(zoo)
library(tidyverse)
library(sf)

suppressWarnings(expr)

setwd('C:\\Projetos\\bioflore\\gbs-uganda')

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

files <- grep(list.files(path="geo\\sentinel2",full.names = T,recursive = T),
               pattern = '.xml', invert =T, value = T)
out <- files[lengths(files)!=0]


endmembers<-read.csv("spreadsheets/em.txt")
colnames(endmembers)<-c("class","Blue","Green","Red","NIR")

melted <- melt(endmembers, id.vars="class")
colnames(melted)<-c("Class","WaveLength","Reflectance")
p <- ggplot(melted, aes(WaveLength, Reflectance, color = Class)) +
  geom_line(aes(group=Class)) + ggtitle("Madagascar Reflectance by class")
p

shp = st_read("geo\\shp\\restoration_sites.shp")
shp <- transform_to_utm(shp)

for (name in out){
  # name = out[[1]]
  label <- str_split(str_split(name, '_')[[1]][1], '/')[[1]][2]
  product <- str_replace_all(name, c(sentinel2= "NDFI", .tif="_NDFI.tif"))
  if (file.exists(product)){
    next
  }
  
  data<-brick(name)
  crs(data)
  bnames <- c('B2','B3','B4','B8','B11',
              'B12','SCL')
  names(data) <- bnames
  
  scl <- data[['SCL']]
  valid_mask <- (scl == 4 | scl == 5)
  data <- mask(data, valid_mask)
  imgBrick <- data[[c('B2','B3','B4','B8')]]
  
  # Altere o nome das bandas
  imgBrick@data@names<-c("Blue","Green","Red","NIR")
    
  polygon <- shp[shp$Name== label,]$geometry
 
  cropped_raster <- crop(imgBrick, as_Spatial(st_zm(polygon)))
  masked_raster <- mask(cropped_raster, as(st_zm(polygon), "Spatial"))
  # 2) Aplique a fun??o 'mesma' para obter as imagens fra??o.
  imgFracao<-mesma(masked_raster, endmembers)
  
  gv_shade <- imgFracao$GV/(1-imgFracao$Shade)
  NDFI <- (gv_shade - (imgFracao$NPV+imgFracao$Soil))/(imgFracao$NPV+imgFracao$Soil+gv_shade)
  plot(NDFI)
  # plot(NDFI,col=gray.colors(256))
  dir.create(paste0('geo\\NDFI\\',label), showWarnings = FALSE)
  writeRaster(NDFI,product)
}

ndfi.files <- list.files(path="geo\\NDFI", full.names = T, recursive = T)
result <- data.frame()

for (name in ndfi.files){
  date <- as.Date(str_split(name, '_')[[1]][2], "%Y-%m-%d")
  talhao <- str_split(str_split(name, '_')[[1]][1], '/')[[1]][2]
  data <- brick(name)[[1]]
  x.stats <- data.frame(talhao=talhao,
    data=date,
    ndfi=cellStats(data, "mean"))
  
  result <- rbind(result,x.stats)
}

colnames(result) <- c("talhao", "date", "ndfi_mean")

result <- result %>%
  filter(date > as.POSIXct("2020-01-01"))

write.csv(result,
          "C:\\Projetos\\bioflore\\gbs-uganda\\spreadsheets\\result_sites_ndfi.csv",
          row.names = FALSE)

result <- read.csv("spreadsheets\\result_sites_ndfi.csv")
result = result[complete.cases(result), ]

RE = 'Reference'
BRE = 'degraded'
labels <- as.list(unique(result$talhao))
labels = labels[! labels %in% c(RE,BRE)] 

for (label in labels){
  # label = labels[[1]]
  vc <- c(label,RE,BRE)
  data<-result[result$talhao %in% vc,]
  data <- data[!duplicated(data[c('talhao','date')]),]
  
  pivoted_data <- data %>% pivot_wider(
    names_from = talhao,
    values_from = ndfi_mean
  )
  # Prepare data for plotting
  pivoted_data$date <- as.Date(pivoted_data$date) # Ensure dates are in Date format
  pivoted_data <- pivoted_data[order(pivoted_data$date),]
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(date)
  #pivoted_data %>% filter(date >= date_non_na)
  # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$date 
  end_date <- pivoted_data[nrow(pivoted_data),1]$date
  
  # Open a new graphics window (specific to Windows)
  index_name = "NDFI"
  #win.graph(10.90625, 16.83333)
  
  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(-1.0, 1.0), main = label, xaxt = 'n',cex.lab = 1.5)
  
  # Manually add x-axis labels at a 45-degree angle, showing labels from Apr 2018 to Apr 2022 every 3 months
  date_ticks <- seq(from = start_date, to = end_date, by = "5 months")
  axis(1, at = as.numeric(date_ticks), labels = FALSE) # Add tick marks without labels
  text(x = as.numeric(date_ticks), y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]),
       labels = format(date_ticks, "%b %Y"), srt = 45, adj = 1.3, xpd = TRUE)
  
  x <- pivoted_data[[label]]
  x_re <- pivoted_data[[RE]]
  x_bre <- pivoted_data[[BRE]]
  # Forward fill interpolation
  x <-na.locf(x, na.rm = FALSE)
  x_re<- na.locf(x_re, na.rm = FALSE)
  x_bre<- na.locf(x_bre, na.rm = FALSE)
  y <- pivoted_data$date
  # y_re <- $date
  degree <- 15 # Degree of the spline
  df <- 5 # Degrees of freedom or roughly the number of knots
  
  # Fit B-spline
  fit_bs <- lm(x ~ bs(y, degree = degree, df = df))
  fit_bs_re <- lm(x_re ~ bs(y, degree = degree, df = df))
  fit_bs_bre <- lm(x_bre ~ bs(y, degree = degree, df = df))
  # Predictions for a smooth curve across the entire date range
  y_pred <- seq(date_non_na, end_date, by = "day")
  x_pred <- predict(fit_bs, newdata = list(y = y_pred))
  
  y_pred_re <- seq(start_date, end_date, by = "day")
  x_pred_re <- predict(fit_bs_re, newdata = list(y = y_pred_re))
  x_pred_bre <- predict(fit_bs_bre, newdata = list(y = y_pred_re))
  
  max_value <- max(x_pred)
  if (max_value > 1){
    x_pred <- x_pred / max_value
  }
  
  min_value <- min(x_pred)
  if (min_value < -1){
    x_pred <- x_pred / -min_value
  }
  
  max_value_re <- max(x_pred_re)
  if (max_value_re > 1){
    x_pred_re <- x_pred_re / max_value_re
  }
  
  min_value_re <- min(x_pred_re)
  if (min_value_re < -1){
    x_pred_re <- x_pred_re / -min_value_re
  }
  
  max_value_bre <- max(x_pred_bre)
  if (max_value_bre > 1){
    x_pred_bre <- x_pred_bre / max_value_bre
  }
  
  min_value_bre <- min(x_pred_bre)
  if (min_value_bre < -1){
    x_pred_bre <- x_pred_bre / -min_value_bre
  }
  
  # Plot the B-spline curve for this typology
  lines(y_pred, x_pred, col = 'black', lwd=3.5)
  lines(y_pred_re, x_pred_re, col = 'green', lwd=3.5)
  lines(y_pred_re, x_pred_bre, col = 'orange', lwd=3.5)
  
  # specific_date <- as.Date(shp[shp$identifier_akvo==label,]$planting_date)
  # 
  # # Add intermittent vertical line
  # n_segments <- 50  # Number of segments for the intermittent line
  # x_pos <- as.numeric(specific_date)
  # ymin <- -1
  # ymax <- 1
  # segment_length <- (ymax - ymin) / (n_segments * 2)
  # 
  # for (i in 0:(n_segments-1)) {
  #   segments(x_pos, ymin + i * 2 * segment_length, x_pos, ymin + (i * 2 + 1) * segment_length, col = "red")
  # }
  
  # legend('topleft', legend=c("Restoration Area", "Reference Area", 'Planting Date'),
  #        col=c("black", "green", 'red'), lty=c(1,1,2), cex=0.8,lwd=c(3.5,3.5,1))
  
  legend('bottomright', legend=c("Restoration Area", "Reference Area", "Degraded Area"),
         col=c("black", "green", 'orange'), lty=c(1,1,1), cex=0.8,lwd=c(3.5,3.5,3.5))
  
  # Initialize vectors to store the minimum points for connecting them later
  min_points_x <- c()
  min_points_y <- c()
  
  # Find and plot minimum points along the spline curve for each year within the range
  years <- unique(format(y_pred, "%Y"))
  for (year in years) {
    in_year <- format(y_pred, "%Y") == year
    if (sum(in_year) > 0) {
      # Find the minimum value of the spline within this year
      min_val_index <- which.min(x_pred[in_year])
      min_date <- y_pred[in_year][min_val_index]
      min_val <- x_pred[in_year][min_val_index]
      
      # Accumulate the minimum points
      min_points_x <- c(min_points_x, min_date)
      min_points_y <- c(min_points_y, min_val)
      
      # Plot the minimum point
      points(min_date, min_val, col = "#0101ff", pch=19)
    }
  }
  
  # Connect the minimum points with lines
  if (length(min_points_x) > 1) {
    lines(min_points_x, min_points_y, col = "#0101ff", lty=1)
  }
}

#dev.off()

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="C:/Projetos/bioflore/gbs-uganda/png/plots_ndfi")



