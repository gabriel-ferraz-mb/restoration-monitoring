# Plot NDMI ---------------------------------------------------------------

library(dplyr)
library(splines)
library(zoo)
library(tidyverse)
library(sf)
library(raster)
rm(list = ls())

suppressWarnings(expr)

setwd('C:\\Projetos\\bioflore\\gbs-uganda')

files <- grep(list.files(path="geo\\sentinel2",full.names = T, recursive = T),
                    pattern = '.aux', invert = T, value = T)

out <- files[lengths(files)!=0]

# shp = st_read("geo\\Contract129_reference_Polygons_corrected_final.gpkg")
# shp <- st_transform(shp, crs = 32738)

result <- data.frame()
for (name in files){
  # name <- files[[1]]
  date <- as.Date(str_replace(str_split(name, '_')[[1]][-1],'.tif',''),
                              "%Y-%m-%d")
  
  talhao <- str_split(str_split(name, '_')[[1]][1], '/')[[1]][2]
  data <- brick(name)
  bnames <- c('B2','B3','B4','B8','B11',
              'B12','SCL')
  names(data) <- bnames
  
  NIR <- data[['B8']]
  SWIR <- data[['B11']]
  
  # Calculate NDMI
  NDMI <- (NIR - SWIR) / (NIR + SWIR)
  
  scl <- data[['SCL']]
  valid_mask <- (scl == 4 | scl == 5)
  NDMI_masked <- mask(NDMI, valid_mask)
  x.stats <- data.frame(talhao=talhao,
                        data=date,
                        ndmi=cellStats(NDMI_masked, "mean"))
  
  result <- rbind(result,x.stats)
}

colnames(result) <- c("talhao", "date", "ndmi_mean")

result <- result %>%
  filter(date > as.POSIXct("2020-01-01"))

write.csv(result,
          "C:\\Projetos\\bioflore\\gbs-uganda\\spreadsheets\\result_sites_ndmi.csv",
          row.names = FALSE)

tab<-read.table(
  "C:\\Projetos\\bioflore\\gbs-uganda\\spreadsheets\\result_sites_ndmi.csv",
                sep=",")
colnames(tab)<-tab[1,]
tab<-tab[2:nrow(tab),]
tab$ndmi_mean<-as.numeric(tab$ndmi_mean)
tab$date <- as.Date(tab$date)

RE = 'Reference'
BRE = 'degraded'
labels <- as.list(unique(tab$talhao))
labels = labels[! labels %in% c(RE,BRE)]

tab = tab[complete.cases(tab), ]

for (label in labels){
  vc <- c(label,RE,BRE)
  data<-tab[tab$talhao %in% vc,]
  
  pivoted_data <- data %>% pivot_wider(
    id_cols = date,
    names_from = talhao,
    values_from = ndmi_mean
  )
  # Prepare data for plotting
  # pivoted_data$date <- as.Date(pivoted_data$date) # Ensure dates are in Date format
  pivoted_data <- pivoted_data[order(pivoted_data$date),]
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(date)
  pivoted_data %>% filter(date >= date_non_na)
  # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$date 
  end_date <- pivoted_data[nrow(pivoted_data),1]$date
  
  # Open a new graphics window (specific to Windows)
  index_name = "NDMI"
  #win.graph(10.90625, 16.83333)
  
  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(-0.5, 0.5), main = label, xaxt = 'n',cex.lab = 1.5)
  
  # Manually add x-axis labels at a 45-degree angle, showing labels from Apr 2018 to Apr 2022 every 3 months
  date_ticks <- seq(from = start_date, to = end_date, by = "5 months")
  axis(1, at = as.numeric(date_ticks), labels = FALSE) # Add tick marks without labels
  text(x = as.numeric(date_ticks), y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]),
       labels = format(date_ticks, "%b %Y"), srt = 45, adj = 1.3, xpd = TRUE)
  
  x <- pivoted_data[[label]]
  x_re <- pivoted_data[[RE]]
  x_bre <- pivoted_data[[BRE]]
  # Forward fill interpolation
  x <- na.locf(x, na.rm = FALSE)
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
  
  min_value <- min(x_pred)
  if (min_value < -1){
    x_pred <- x_pred / -min_value
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
  # legend('bottomright', legend=c("Restoration Area", "Reference Area","Bad Reference Area", 'Planting Date'),
  #        col=c("black", "green",'orange', 'red'), lty=c(1,1,1,2), cex=0.8,lwd=c(3.5,3.5,3.5,1))
  legend('bottomright', legend=c("Restoration Area", "Reference Area","Degraded Area"),
         col=c("black", "green",'orange'), lty=c(1,1,1), cex=0.8,lwd=c(3.5,3.5,3.5))

}

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="C:\\Projetos\\bioflore\\gbs-uganda\\png\\plots_ndmi")
