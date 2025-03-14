library(dplyr)
library(splines)
library(zoo)
library(sf)
library(raster)
library(tidyverse)
rm(list = ls())

suppressWarnings(expr)

setwd('C:\\Projetos\\bioflore\\bgci_II')

files <- grep(list.files(path="geo\\sentinel2",full.names = T, recursive = T),
                    pattern = '.xml', invert = T, value = T)

out <- files[lengths(files)!=0]

result <- data.frame()
for (name in files){
  # name = files[[1]]
  date <- as.Date(str_replace(str_split(name, '_')[[1]][-1],'.tif',''),
                  "%Y-%m-%d")
  
  talhao <- str_split(name, '/')[[1]][2]
  data <- brick(name)
  
  # bnames <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11',
  #             'B12','AOT','WVP','SCL','TCI_R','TCI_G','TCI_B','MSK_CLDPRB',
  #             'MSK_SNWPRB','QA10','QA20','QA60','MSK_CLASSI_OPAQUE',
  #             'MSK_CLASSI_CIRRUS','MSK_CLASSI_SNOW_ICE')
  bnames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL')
  
  names(data) <- bnames
  
  NIR <- data[['B8']]
  red <- data[['B4']]
  
  # Calculate NDMI
  NDVI <- (NIR - red) / (NIR + red)
  
  scl <- data[['SCL']]
  # plot(scl)
  valid_mask <- (scl == 4 | scl == 5)
  # plot(valid_mask)
  NDVI_masked <- mask(NDVI, valid_mask, maskvalue = TRUE,inverse = TRUE)
  # plot(NDVI_masked)
  x.stats <- data.frame(talhao=talhao,
                        data=date,
                        ndvi=cellStats(NDVI_masked, "mean"))
  
  result <- rbind(result,x.stats)
}

colnames(result) <- c("talhao", "date", "ndvi_mean")

tname =  "spreadsheets\\result_sites_ndvi.csv"
write.csv(result,
          tname,
          row.names = FALSE)

# tab<-read.table(tname, sep=",")
# colnames(tab)<-tab[1,]
# tab<-tab[2:nrow(tab),]
# tab$ndvi_mean<-as.numeric(tab$ndvi_mean)
# tab$date <- as.Date(tab$date)

tab <- result

tab = tab[complete.cases(tab), ]

result_bspline <- data.frame(date=as.Date(character()),
                     ndvi_mean=numeric(), 
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
    values_from = ndvi_mean
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
  names(df) <- c('date','ndvi_mean','talhao')
  #print(df)
  result_bspline <- rbind(result_bspline, df)
}
write.csv(result_bspline, "spreadsheets/bspline_ndvi_limite.csv", row.names=FALSE)

pivoted_data <- result_bspline %>% pivot_wider(
  id_cols = date,
  names_from = talhao,
  values_from = ndvi_mean
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
  RE = paste0(label, ".ref")
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(date)
  pivoted_data %>% filter(date >= date_non_na)
  # # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$date
  end_date <- pivoted_data[nrow(pivoted_data),1]$date
  # Open a new graphics window (specific to Windows)
  index_name = "NDVI"
  #win.graph(10.90625, 16.83333)

  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(-0.1, 1), main = label, xaxt = 'n',cex.lab = 1.5)

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
    ymin <- -0.1  # Ensure this value suits your plot's y-limits
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
}


plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="png\\plots_ndvi")
