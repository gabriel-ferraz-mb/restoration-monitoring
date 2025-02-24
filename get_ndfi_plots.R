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

setwd('C:\\Projetos\\bioflore\\bgci_II')

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
  geom_line(aes(group=Class)) + ggtitle('Kenya project "Reflectance by class"')
p

shp = st_read(
  "geo\\sites\\SITES+REF\\all_sites.shp"
  )

# label_list <- c('La_Pena', 'HUB_CERRADO_CECAP', 'HUB_CERRADO_NA_FLORESTA')  # Replace with your actual list of strings
# # Filter the rows where 'label' is in the list
# shp <- shp %>%
#   filter(name %in% label_list)%>%
#   mutate(name= gsub("_", "-", name))

for (name in out){
  # name = out[[1]]
  label <- str_split(name, '/')[[1]][2]
  product <- str_replace_all(name, c(sentinel2= "NDFI", .tif="_NDFI.tif"))
  
  if (file.exists(product)){
    next
  }
  
  data<-brick(name)
  crs(data)
  bnames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL')
  names(data) <- bnames
  
  scl <- data[['SCL']]
  valid_mask <- (scl == 4 | scl == 5)
  data <- mask(data, valid_mask, maskvalue = TRUE,inverse = TRUE)
  imgBrick <- data[[c('B2','B3','B4','B8')]]
  
  # Altere o nome das bandas
  imgBrick@data@names<-c("Blue","Green","Red","NIR")
    
  polygon <- shp[shp$Name== label,]$geometry
  polygon <- st_transform(polygon, crs = as.character(crs(imgBrick)))
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

patterns <- c('HUB-CERRADO', 'La-Pena')  # Add any additional patterns here
# ndfi.files <- ndfi.files[vapply(
#   ndfi.files, function(file) any(sapply(patterns, function(pattern) grepl(pattern, file))), logical(1))]

result <- data.frame()

for (name in ndfi.files){
  # name = ndfi.files[[1]]
  date <- as.Date(
    str_split(str_split(name, '/')[[1]][3],"_")[[1]][2],
  "%Y-%m-%d")
  talhao <- str_split(str_split(name, '/')[[1]][3],"_")[[1]][1]
  data <- brick(name)[[1]]
  x.stats <- data.frame(talhao=talhao,
    data=date,
    ndfi=cellStats(data, "mean"))
  
  result <- rbind(result,x.stats)
}

colnames(result) <- c("talhao", "date", "ndfi_mean")

write.csv(result,
          "spreadsheets\\result_sites_ndfi.csv",
          row.names = FALSE)

# tab <- read.csv("spreadsheets\\result_sites_ndfi.csv")
# tab$ndfi_mean<-as.numeric(tab$ndfi_mean)
# tab$date <- as.Date(tab$date)
tab <- result
tab = tab[complete.cases(tab), ]

result_bspline <- data.frame(date=as.Date(character()),
                             ndfi_mean=numeric(), 
                             talhao=character(), 
                             stringsAsFactors=FALSE)

labels <- as.list(unique(tab$talhao))
for (label in labels){
  #label = "ESP_RE"
  print(label)
  vc <- c(label)
  data<-tab[tab$talhao %in% vc,]
  data <- data[!duplicated(data[c('talhao','date')]),]
  
  pivoted_data <- data %>% pivot_wider(
    names_from = talhao,
    values_from = ndfi_mean
  )
  # Prepare data for plotting
  pivoted_data$date <- as.Date(pivoted_data$date) # Ensure dates are in Date format
  pivoted_data <- pivoted_data[order(pivoted_data$date),]
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(date)
  pivoted_data %>% filter(date >= date_non_na)
  # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$date 
  end_date <- pivoted_data[nrow(pivoted_data),1]$date
  
  # Open a new graphics window (specific to Windows)
  
  x <- pivoted_data[[label]]
  # Forward fill interpolation
  x <- na.locf(x, na.rm = FALSE)
  y <- pivoted_data$date
  #y_re <- IPE_RE$date
  degree <-20 # Degree of the spline
  df <- 15 # Degrees of freedom or roughly the number of knots
  
  # Fit B-spline
  fit_bs <- lm(x ~ bs(y, degree = degree, df = df))
  # Predictions for a smooth curve across the entire date range
  y_pred <- seq(start_date, end_date, by = "day")
  x_pred <- predict(fit_bs, newdata = list(y = y_pred))
  l = rep( label , length( x_pred)  ) 
  
  df <- data.frame(y_pred[(31):(length(y_pred) - 30)],
                   x_pred[(31):(length(x_pred) - 30)],
                   l[(31):(length(l) - 30)])
  names(df) <- c('date','ndfi_mean','talhao')
  #print(df)
  result_bspline <- rbind(result_bspline, df)
}
write.csv(result_bspline, "spreadsheets/bspline_ndfi_lmite.csv", row.names=FALSE)

result_bspline<-read.table("spreadsheets/bspline_ndfi_lmite.csv", sep=",")
colnames(result_bspline)<-result_bspline[1,]
result_bspline<-result_bspline[2:nrow(result_bspline),]
result_bspline$ndvi_mean<-as.numeric(result_bspline$ndfi_mean)
result_bspline$date <- as.Date(result_bspline$date)


pivoted_data <- result_bspline %>% pivot_wider(
  id_cols = date,
  names_from = talhao,
  values_from = ndfi_mean
)
# Prepare data for plotting
pivoted_data$date <- as.Date(pivoted_data$date) # Ensure dates are in Date format
pivoted_data <- pivoted_data[order(pivoted_data$date),]

labels <- list(
  c("ABP","2010-01-01"),
  c("Bigodi","1998-01-01"),
  c("Dilips","2020-01-01"),
  c("IPEECOSIA","2020-01-01"),
  c("Jubiya","2021-01-01"),
  c("KB","2016-01-01"),
  c("KOREN","2021-01-01"),
  c("LUM","2017-01-01"),
  c("Magaca","2021-01-01"),
  c("MARLENE","2019-06-01"),
  c("Nteko","2014-01-01"),
  c("Pandal","2019-06-01"),
  c("RPPNRENOPOLIS","1950-01-01"),
  c("Soraypampa","2015-01-01"),
  c("Wangari","2018-01-01"),
  c("HUB-CERRADO-NA-FLORESTA","2000-01-01"),
  c("HUB-CERRADO-CECAP","2017-01-01"),
  c("La-Pena","2021-01-01")
)

for (i in labels){
  
  label = i[[1]]
  specific_date = i[[2]]
  
  label = str_replace_all(label,"_","-")
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(date)
  pivoted_data %>% filter(date >= date_non_na)
  # # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$date
  end_date <- pivoted_data[nrow(pivoted_data),1]$date
  
  # Open a new graphics window (specific to Windows)
  RE = paste0(label, ".ref")
  index_name = "NDFI"
  #win.graph(10.90625, 16.83333)
  
  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(-1, 1), main = label, xaxt = 'n',cex.lab = 1.5)
  
  # Manually add x-axis labels at a 45-degree angle, showing labels from Apr 2018 to Apr 2022 every 3 months
  date_ticks <- seq(from = start_date, to = end_date, by = "5 months")
  axis(1, at = as.numeric(date_ticks), labels = FALSE) # Add tick marks without labels
  text(x = as.numeric(date_ticks), y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]),
       labels = format(date_ticks, "%b %Y"), srt = 45, adj = 1.3, xpd = TRUE)
  x <- pivoted_data[[label]]
  x_re <- pivoted_data[[RE]]
  
  y <- pivoted_data$date
  
  # Plot the B-spline curve for this typology
  lines(y, x, col = 'black', lwd=3.5)
  lines(y, x_re, col = 'green', lwd=3.5)
  
  
  specific_date_index <- which(y == specific_date)
  if (length(specific_date_index) > 0) {
    # If the specific date is found, get the numeric position for plotting
    x_pos <- as.numeric(y[specific_date_index])
    
    # Add intermittent vertical line at the specific date
    n_segments <- 50
    ymin <- -1  # Ensure this value suits your plot's y-limits
    ymax <- 1
    segment_length <- (ymax - ymin) / (n_segments * 2)
    
    for (j in 0:(n_segments - 1)) {
      segments(x_pos, ymin + j * 2 * segment_length, x_pos, ymin + (j * 2 + 1) * segment_length, col = "red")
    }
    legend('bottomright', legend=c("Restoration Area", "Reference Area","Baseline"),
           col=c("black", "green",'red'), lty=c(1,1,1,2), cex=0.8,lwd=c(3.5,3.5,3.5,1))
  } else {
    warning(paste("Specific date", specific_date, "not found in dates vector for label", label))
    legend('bottomright', legend=c("Restoration Area", "Reference Area"),
           col=c("black", "green"), lty=c(1,1,1,2), cex=0.8,lwd=c(3.5,3.5,3.5,1))
  }
  
  min_points_x <- c()
  min_points_y <- c()
  
  # Find and plot minimum points along the spline curve for each year within the range
  # years <- unique(format(y, "%Y"))
  # for (year in years) {
  #   in_year <- format(y, "%Y") == year
  #   if (sum(in_year) > 0) {
  #     # Find the minimum value of the spline within this year
  #     min_val_index <- which.min(x[in_year])
  #     min_date <- y[in_year][min_val_index]
  #     min_val <- x[in_year][min_val_index]
  #     
  #     # Accumulate the minimum points
  #     min_points_x <- c(min_points_x, min_date)
  #     min_points_y <- c(min_points_y, min_val)
  #     
  #     # Plot the minimum point
  #     points(min_date, min_val, col = "#0101ff", pch=19)
  #   }
  # }
  # 
  # # Connect the minimum points with lines
  # if (length(min_points_x) > 1) {
  #   lines(min_points_x, min_points_y, col = "#0101ff", lty=1)
  # }
}

#dev.off()

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="png/plots_ndfi")



