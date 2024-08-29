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

result <- data.frame(date=as.Date(character()),
                 ndvi_mean=numeric(), 
                 talhao=character(), 
                 stringsAsFactors=FALSE)

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
  degree <- 15 # Degree of the spline
  df <- 5 # Degrees of freedom or roughly the number of knots
  
  # Fit B-spline
  fit_bs <- lm(x ~ bs(y, degree = degree, df = df))
  # Predictions for a smooth curve across the entire date range
  y_pred <- seq(date_non_na, end_date, by = "day")
  x_pred <- predict(fit_bs, newdata = list(y = y_pred))
  l = rep( label , length( x_pred)  ) 
  
  df <- data.frame(y_pred,
                   x_pred,
                   l)
  names(df) <- c('date','ndvi_mean','talhao')
  #print(df)
  result <- rbind(result, df)
}

write.csv(result, "C:\\Projetos\\bioflore\\ecosia\\Ecosia\\bspline_ndvi.csv", row.names=FALSE)
#combined_df <-  do.call(rbind, df_list)
#dev.off()
