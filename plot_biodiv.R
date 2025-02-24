
# script_path <- rstudioapi::getActiveDocumentContext()$path
# funchir::stale_package_check(script_path)

library(utils)
library(terra)
library(ggplot2)
library(lubridate)
library(glue)
library(sf)
library(stringr)

setwd("C:\\Projetos\\bioflore\\bgci_II")

################################ Define Fuctions ###############################

plot_alpha_get_stats <- function(label, date, cropped_alpha, matching_polygon){
  x.stats <- data.frame(talhao=label,
                        data=date,
                        shannon=global(cropped_alpha, fun = "mean", na.rm = TRUE))
  # result_df <- rbind(result_df,x.stats)
  # Convert raster to data frame and rename columns appropriately
  raster_df_alpha <- as.data.frame(terra::as.data.frame(cropped_alpha, xy = TRUE), xy = TRUE)
  names(raster_df_alpha) <- c("lon", "lat", "shannon")
  max_shannon <- max(raster_df_alpha$value, na.rm=T)
  
  # Create a ggplot
  p <- ggplot() +
    geom_raster(data = raster_df_alpha, aes(x = lon, y = lat, fill = shannon)) +
    geom_sf(data = matching_polygon, fill = NA, color = "red", lwd = 1) +
    theme_minimal() +
    scale_fill_gradientn(limits = c(0,2), colors = scales::viridis_pal()(7)) + 
    labs(title = paste("Alpha Diversitiy for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ymd")))) +
    theme( plot.title = element_text(hjust = 0.5, size = 16),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
           axis.text.y = element_text(size = 8))
  
  # Save or display the plot
  ggsave(glue("png\\biodiv\\shannon_{label}_{date}.png"), plot = p)
  
  return(x.stats)
}

plot_beta <- function(label, date, cropped_beta, matching_polygon){
  raster_df_beta <- as.data.frame(cropped_beta, xy = TRUE)
  colnames(raster_df_beta) <- c('x', 'y', 'Red', 'Green', 'Blue') 
  raster_df_beta$Red <- pmax(raster_df_beta$Red, 0) / 255
  raster_df_beta$Green <- pmax(raster_df_beta$Green, 0) / 255
  raster_df_beta$Blue <- pmax(raster_df_beta$Blue, 0) / 255
  
  maxVal <- max(raster_df_beta[c("Red","Green","Blue")])
  
  b <- ggplot() + 
    geom_raster(data = raster_df_beta, aes(x = x, y = y, fill = rgb(Red, Green, Blue,
                                                                    maxColorValue=maxVal))) +
    scale_fill_identity() +
    theme_minimal() +
    geom_sf(data = matching_polygon, fill = NA, color = "red", lwd = 1) +
    # geom_path(data=polygon_fortified, aes(x = x, y = y, group = group), color = "red", size = 1) +
    xlab('Longitude') + ylab('Latitude') +
    ggtitle(paste("Beta Diversitiy for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ymd")))) +
    theme_bw() + 
    coord_equal() +
    coord_sf() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8)
    )
  ggsave(glue("png\\biodiv\\map_{label}_{date}.png"), plot = b)
}

plot_scatter_beta <- function(label, date, cropped_beta){
  band1 <- values(cropped_beta [[1]])
  band2 <- values(cropped_beta [[2]])
  band3 <- values(cropped_beta [[3]])

  # Combine bands into a data frame
  data <- data.frame(Band1 = band1, Band2 = band2, Band3 = band3)

  # Remove any NA values
  data <- na.omit(data)

  ### Creating Scatter Plots
  ##Now, let's create scatter plots comparing each pair of bands using `ggplot2`.

  # Scatter plot for Band 1 vs Band 2
  plot1 <- ggplot(data, aes(x = Band1, y = Band2)) +
    geom_point(alpha = 0.5) +
    labs(
         x = "Band 1",
         y = "Band 2") +
    theme_minimal()

  # Scatter plot for Band 1 vs Band 3
  plot2 <- ggplot(data, aes(x = Band1, y = Band3)) +
    geom_point(alpha = 0.5) +
    labs(
         x = "Band 1",
         y = "Band 3") +
    theme_minimal()

  # Scatter plot for Band 2 vs Band 3
  plot3 <- ggplot(data, aes(x = Band2, y = Band3)) +
    geom_point(alpha = 0.5) +
    labs(
         x = "Band 2",
         y = "Band 3") +
    theme_minimal()

  g <- grid.arrange(
    arrangeGrob(plot1, plot2, plot3, ncol = 3),
    top = textGrob(paste("Scatter plot of PCoA bands for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ym"))),
                   gp = gpar(fontsize = 10))
  )

  ggsave(file=glue("png\\biodiv\\scatter_{label}_{date}.png"), g)

}

################################ Start Processes ###############################

raster_files <-  grep(list.files(path="geo\\biodiv_output",
                                 full.names = T,recursive = T),
                      pattern = 'spectral', invert =T, value = T)
out <- raster_files[lengths(raster_files)!=0]

shapefile_path <- "geo\\sites\\SITES+REF\\all_sites.shp"
polygons <- st_read(shapefile_path)
# polygons <- transform_to_utm(polygons)

result_df <- data.frame()
for (file in raster_files){
  # file <- raster_files[1]
  
  label <- str_split(str_split(file, '/')[[1]][2], '_')[[1]][1]
  date <- str_remove(
    str_split(str_split(file, '/')[[1]][2], '_')[[1]][4],".tif")
  
  matching_polygon <- polygons[polygons$Name == label,]
  
  raster_data <- rast(file)
  
  matching_polygon <- st_transform(matching_polygon, crs = as.character(crs(raster_data)))
  
  cropped_data <- crop(raster_data, as_Spatial(st_zm(matching_polygon$geom)))
  
  if (str_detect(file, "alpha")){
    x.stats <- plot_alpha_get_stats(label, date, cropped_data,matching_polygon)
    result_df <- rbind(result_df,x.stats)
  } else {
    plot_beta(label, date, cropped_data,matching_polygon)
  }
}

write.csv(result_df,
          file = "spreadsheets\\result_shannon_limite.csv",
          row.names = FALSE)


result_df$talhao <- snakecase::to_title_case(result_df$talhao)

# custom_order <- c( "Reference 1", "Reference 2", "Lower Imenti", "Gathiuru",
#                    "Hombe", "Hombe 2", "Hombe 3", "Kangaita", "Untitled",
#                    "Ontulili")
# 
# # Convert the 'talhao' column to a factor with the specified levels
# result_df$talhao <- factor(result_df$talhao, levels = custom_order)
for(label in unique(result_df$talhao)) {
  row <- result_df[result_df$talhao == label,]
  p <- ggplot(row, aes(x = talhao, y = mean, fill = data)) + 
    geom_bar(stat = 'identity', position = position_dodge(width = 0.9), color = "black") +  # Add border with color
    geom_text(aes(label = data), position = position_dodge(width = 0.9), angle = 90,
              vjust = 0.5, hjust = 1.5, color = "black", size = 3) +
    labs(x = 'Date', y = 'Shannon Index Mean Value', title = glue('Shannon index mean by date {label}')) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis tick labels
      axis.ticks.x = element_blank(), # Remove x-axis ticks
      # axis.title.x = element_blank(), # Remove x-axis title if necessary
      plot.title = element_text(hjust = 0.5), # Center title
      legend.position = "none"
    )
  
  # print(p)
  ggsave(file=glue("png\\biodiv\\barplot_shannon_{label}.png"), p)
}

