

# Plot NDVI ---------------------------------------------------------------

library(dplyr)
library(splines)
library(zoo)
library(tidyverse)
rm(list = ls())

tab<-read.table("C:/Projetos/bioflore/ecosia/resultado_22062024.csv", sep=",")
colnames(tab)<-tab[1,]
tab<-tab[2:nrow(tab),]
tab$ndvi_mean<-as.numeric(tab$ndvi_mean)

labels <- as.list(unique(tab$talhao))
for (label in labels){
  #png(paste0(label, ".png"))
  if (startsWith(label, "IPE")){
    RE = 'IPE_RE'
  } else{
    RE = 'ESP_RE'
  }
  
  vc <- c(label,RE)
  data<-tab[tab$talhao %in% vc,]
  
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
  index_name = "NDVI"
  #win.graph(10.90625, 16.83333)
  
  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(0, 1), main = label, xaxt = 'n',cex.lab = 1.5)
  
  # Manually add x-axis labels at a 45-degree angle, showing labels from Apr 2018 to Apr 2022 every 3 months
  date_ticks <- seq(from = start_date, to = end_date, by = "5 months")
  axis(1, at = as.numeric(date_ticks), labels = FALSE) # Add tick marks without labels
  text(x = as.numeric(date_ticks), y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]),
       labels = format(date_ticks, "%b %Y"), srt = 45, adj = 1.3, xpd = TRUE)
  
  x <- pivoted_data[[label]]
  x_re <- pivoted_data[[RE]]
  # Forward fill interpolation
  x <- na.locf(x, na.rm = FALSE)
  x_re<- na.locf(x_re, na.rm = FALSE)
  y <- pivoted_data$date
  #y_re <- IPE_RE$date
  degree <- 15 # Degree of the spline
  df <- 5 # Degrees of freedom or roughly the number of knots
  
  # Fit B-spline
  fit_bs <- lm(x ~ bs(y, degree = degree, df = df))
  fit_bs_re <- lm(x_re ~ bs(y, degree = degree, df = df))
  # Predictions for a smooth curve across the entire date range
  y_pred <- seq(date_non_na, end_date, by = "day")
  x_pred <- predict(fit_bs, newdata = list(y = y_pred))
  
  y_pred_re <- seq(start_date, end_date, by = "day")
  x_pred_re <- predict(fit_bs_re, newdata = list(y = y_pred_re))
  
  # Plot the B-spline curve for this typology
  lines(y_pred, x_pred, col = 'black', lwd=3.5)
  lines(y_pred_re, x_pred_re, col = 'green', lwd=3.5)
  
  legend('bottomright', legend=c("Restoration Area", "Reference Area"),
         col=c("black", "green"), lty=1:1, cex=0.8,lwd=3.5)
  
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
file.copy(from=plots.png.paths, to="C:/Projetos/bioflore/ecosia/plots_ndvi")
