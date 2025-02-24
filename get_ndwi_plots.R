# Plot NDWI ---------------------------------------------------------------

library(dplyr)
library(splines)
library(zoo)
library(tidyverse)
library(sf)
library(raster)
rm(list = ls())

suppressWarnings(expr)

setwd('C:\\Projetos\\bioflore\\bgci_II')

files <- grep(list.files(path="geo\\sentinel2_water",full.names = T, recursive = T),
              pattern = '.aux', invert = T, value = T)

out <- files[lengths(files)!=0]

# shp = st_read("geo\\Contract129_reference_Polygons_corrected_final.gpkg")
# shp <- st_transform(shp, crs = 32738)

result <- data.frame()
for (name in files){
  # name <- files[[1]]
  date <- as.Date(str_replace(str_split(name, '_')[[1]][-1],'.tif',''),
                  "%Y-%m-%d")
  
  talhao <-  str_split(name, '/')[[1]][2]
  data <- brick(name)
  bnames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL')
  names(data) <- bnames
  
  nir <- data[['B8']]
  green <- data[['B3']]
  
  # Calculate ndwi
  NDWI <- (green - nir) / (green + nir)
  plot(NDWI)
  
  scl <- data[['SCL']]
  plot(scl)
  valid_mask <- scl == 6
  plot(valid_mask)
  NDWI_masked <- mask(NDWI, mask = valid_mask, maskvalue = TRUE,inverse = TRUE)
  plot(NDWI_masked)
  x.stats <- data.frame(talhao=talhao,
                        data=date,
                        ndwi=cellStats(NDWI_masked, "mean"))
  
  result <- rbind(result,x.stats)
}

colnames(result) <- c("talhao", "date", "ndwi_mean")

write.csv(result,
          "spreadsheets\\result_sites_ndwi.csv",
          row.names = FALSE)

tab <- result
tab = tab[complete.cases(tab), ]

result_bspline <- data.frame(date=as.Date(character()),
                             ndwi_mean=numeric(), 
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
    values_from = ndwi_mean
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
  names(df) <- c('date','ndwi_mean','talhao')
  #print(df)
  result_bspline <- rbind(result_bspline, df)
}
write.csv(result_bspline, "spreadsheets/bspline_ndwi_lmite.csv", row.names=FALSE)
pivoted_data <- result_bspline %>% pivot_wider(
  id_cols = date,
  names_from = talhao,
  values_from = ndwi_mean
)
# Prepare data for plotting
pivoted_data$date <- as.Date(pivoted_data$date) # Ensure dates are in Date format
pivoted_data <- pivoted_data[order(pivoted_data$date),]

# RE = 'reference2'
BRE = ''
# 
# labels = labels[! labels %in% c(RE,BRE)]

labels <- list(
  c("LUM","2017-01-01"),
  c("Pandal","2019-06-01"),
  c("Siru","2019-06-01")
)

for (i in labels){
  label = i[[1]]
  specific_date = i[[2]]
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(date)
  pivoted_data %>% filter(date >= date_non_na)
  # # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$date
  end_date <- pivoted_data[nrow(pivoted_data),1]$date
  # label = labels[[1]]
  # Open a new graphics window (specific to Windows)
  # RE = paste0(label, ".ref")
  
  index_name = "NDWI"
  
  #win.graph(10.90625, 16.83333)
  
  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(-0.5, 1), main = label, xaxt = 'n',cex.lab = 1.5)
  
  # Manually add x-axis labels at a 45-degree angle, showing labels from Apr 2018 to Apr 2022 every 3 months
  date_ticks <- seq(from = start_date, to = end_date, by = "5 months")
  axis(1, at = as.numeric(date_ticks), labels = FALSE) # Add tick marks without labels
  text(x = as.numeric(date_ticks), y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]),
       labels = format(date_ticks, "%b %Y"), srt = 45, adj = 1.3, xpd = TRUE)
  x <- pivoted_data[[label]]
  # x_re <- pivoted_data[[RE]]
  
  y <- pivoted_data$date
  
  # Plot the B-spline curve for this typology
  lines(y, x, col = 'navy', lwd=3.5)
  # lines(y, x_re, col = 'green', lwd=3.5)
  
 
  specific_date_index <- which(y == specific_date)
  if (length(specific_date_index) > 0) {
    # If the specific date is found, get the numeric position for plotting
    x_pos <- as.numeric(y[specific_date_index])
    
    # Add intermittent vertical line at the specific date
    n_segments <- 50
    ymin <- -0.5  # Ensure this value suits your plot's y-limits
    ymax <- 1
    segment_length <- (ymax - ymin) / (n_segments * 2)
    
    for (j in 0:(n_segments - 1)) {
      segments(x_pos, ymin + j * 2 * segment_length, x_pos, ymin + (j * 2 + 1) * segment_length, col = "red")
    }
    legend('bottomright', legend=c("Water bodie NDWI values", "Baseline"),
           col=c("navy", 'red'), lty=c(1,1), cex=0.8,lwd=c(3.5,3.5))
  } else {
    warning(paste("Specific date", specific_date, "not found in dates vector for label", label))
    legend('bottomright', legend=c("Water bodie NDWI values"),
           col=c("navy"), lty=c(1,1), cex=0.8,lwd=c(3.5,3.5))
  }
  
}


plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="png\\plots_ndwi")
