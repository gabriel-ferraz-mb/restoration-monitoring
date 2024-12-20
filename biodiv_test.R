library(biodivMapR)
library(utils)
library(stars)
library(terra)
library(ggplot2)
library(lubridate)
library(dplyr)
library(grid)
library(tidyverse)
library(glue)
library(sf)
library(gridExtra)

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

setwd("C:\\Projetos\\bioflore\\restoration_monitoring\\BiodivMapR")

shapefile_path <- "C:\\Projetos\\bioflore\\gbs-kenya\\geo\\shp\\restoration_sites.shp"
polygons <- st_read(shapefile_path)
# polygons <- transform_to_utm(polygons)
unique_names <- unique(polygons$label)

result_df <- data.frame()
# for (label in unique_names){
# label = unique_names[[1]]
raster_files = grep(list.files(path="C:\\Projetos\\bioflore\\gbs-kenya\\geo\\biodiv_input",
                               full.names = T,recursive = T),
                      pattern = '.xml', invert =T, value = T)
out <- raster_files[lengths(raster_files)!=0]
# matching_polygon <- polygons[polygons$Name == label,]
# matching_polygon <- polygons
years <- seq(2021, 2024, by = 1)

for (year in years) {
  # sentinel2 <- raster_files[1]
  # year = years[1]
  
  year_pattern <- as.character(year)
  year_raster <- out[grepl(year_pattern, out)]
  
  raster_list <- lapply(year_raster, function(file) {
    # Read the raster file
    original_raster <- rast(file)
    # Reproject to CRS 4326
    reprojected_raster <- project(original_raster, "EPSG:4326", res = c(1/11113.95,1/11113.95))
    # plot(reprojected_raster) # This will plot each raster, be careful with large numbers
    return(reprojected_raster)
  })
  
  # Use a reference raster for alignment purposes
  s <- sprc(raster_list)
  # Use mosaic function to combine the aligned rasters
  imgMosaic <- merge(s)
  
  polygons <- st_transform(polygons, crs = as.character(crs(imgMosaic)))
  # imgBrick <- brick(sentinel2)
  # polygon <- matching_polygon$geom
  # tiff_data  <- crop(imgBrick, as_Spatial(st_zm(polygons)))
  # tiff_data  <- mask(tiff_data, as(st_zm(polygons), "Spatial"))
  # 
  # bnames <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11',
  #             'B12','AOT','WVP','SCL','TCI_R','TCI_G','TCI_B','MSK_CLDPRB',
  #             'MSK_SNWPRB','QA10','QA20','QA60','MSK_CLASSI_OPAQUE',
  #             'MSK_CLASSI_CIRRUS','MSK_CLASSI_SNOW_ICE')
  
  bnames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL')
  names(imgMosaic) <- bnames
  
  tiff_data  <- imgMosaic[[c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')]]
  # tiff_data <- img[[c('B2','B3','B4','B8','B11','B12')]]
  
  NameRaster <- 's2a_Subset.tif'
  destfile <- file.path(getwd(),NameRaster)
  
  writeRaster(tiff_data , filename=destfile,
              overwrite=TRUE)
  
  # BandName <- c('band_02', 'band 03', 'band_04', 'band_08', 'band_11', 'band_12')
  # SpectralBands <- c(496.6, 560.0, 664.5, 835.1, 1613.7, 2202.4)
  # WLunits <- 'Nanometers'
  
  # read ENVI file with stars
  # create_hdr(ImPath = destfile, Sensor = 'MyOwnSensor', 
  #            SpectralBands = SpectralBands, BandName = BandName,
  #            WLunits = WLunits)
  
  create_hdr(ImPath = destfile, Sensor = 'SENTINEL_2A')
  
  Input_Image_File <- destfile
  Output_Dir <- 'RESULTS'
  TypePCA <- 'SPCA'
  
  Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File,
                                                   Mask_Path = FALSE,
                                                   Output_Dir = Output_Dir,
                                                   TypePCA = TypePCA,
                                                   NDVI_Thresh = 0.5,
                                                   Blue_Thresh = 800,
                                                   NIR_Thresh = 1200)
  
  mask <- rast(Input_Mask_File)
  # Step 3: Count the pixels for each unique value
  values <- values(mask)
  value_counts <- as.data.frame(table(values))
  valid_pixels <- value_counts[value_counts$values==1,]$Freq

  if(length(valid_pixels) == 0) {
    unlink(list.files(full.names = TRUE, recursive = TRUE), recursive = TRUE)
    next # skip 3rd iteration and go to next iteration
  }
  
  
  FilterPCA <- FALSE
  window_size <- 5
  nbCPU <- 7
  MaxRAM <- 0.9
  nbclusters <- 20
  
  print("PERFORM DIMENSIONALITY REDUCTION")
  PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File, 
                            Input_Mask_File = Input_Mask_File,
                            Output_Dir = Output_Dir, 
                            TypePCA = TypePCA, 
                            FilterPCA = FilterPCA,
                            nbCPU = nbCPU, 
                            MaxRAM = MaxRAM, 
                            Continuum_Removal = F)
  
  # path of the raster resulting from dimensionality reduction
  PCA_Files <- PCA_Output$PCA_Files
  # path for the updated mask
  Input_Mask_File <- PCA_Output$MaskPath
  
  pca_vars <- PCA_Output$PCA_model$sdev^2
  pca_var_percent <- pca_vars / sum(pca_vars) * 100
  cumulative_explained_variance <- cumsum(pca_var_percent)
  # Step 3: Identify the number of components needed to reach at least 75% cumulative variance
  num_components <- which(cumulative_explained_variance >= 75)[1] # Find the first component where cumulative variance >= 75%
  
  lines = paste( unlist(1:num_components), collapse="\n")
  # Write the lines to a text file
  writeLines(lines, con = "RESULTS\\s2a_Subset\\SPCA\\PCA\\Selected_Components.txt")
  
  skip_to_next <- FALSE
  tryCatch(
    expr = {
      print("MAP SPECTRAL SPECIES")
      Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File, 
                                          Input_Mask_File = PCA_Output$MaskPath,
                                          Output_Dir = Output_Dir,
                                          SpectralSpace_Output = PCA_Output, 
                                          nbclusters = nbclusters, 
                                          nbCPU = nbCPU, MaxRAM = MaxRAM)
    },
    error = function(e){
      print(e)
      skip_to_next <<- TRUE
      unlink(list.files(full.names = TRUE, recursive = TRUE), recursive = TRUE)
    }
  )
  
  if(skip_to_next) next
  
  print("MAP ALPHA DIVERSITY")
  # Index.Alpha   = c('Shannon','Simpson')
  Index_Alpha <- c('Shannon')
  map_alpha_div(Input_Image_File = Input_Image_File, 
                Output_Dir = Output_Dir, 
                TypePCA = TypePCA,
                window_size = window_size, 
                nbCPU = nbCPU, 
                MaxRAM = MaxRAM,
                Index_Alpha = Index_Alpha, 
                nbclusters = nbclusters)
  
  print("MAP BETA DIVERSITY")
  map_beta_div(Input_Image_File = Input_Image_File, 
               Output_Dir = Output_Dir, 
               TypePCA = TypePCA,
               window_size = window_size, 
               nbCPU = nbCPU, 
               MaxRAM = MaxRAM,
               nbclusters = nbclusters)
  
  alpha_r = list.files("RESULTS\\s2a_Subset\\SPCA\\ALPHA", full.names = T)
  pat = c("SD","MeanFilter",".hdr")
  combined_pattern <- paste(pat, collapse = "|")
  # Filter the list to exclude strings matching the combined pattern
  filtered_list <- alpha_r[!grepl(combined_pattern, alpha_r)]
  
  alpha_r <- filtered_list[[1]]
  
  beta_r <- grep(list.files(path="RESULTS\\s2a_Subset\\SPCA\\BETA", full.names = T),
                 pattern='.hdr',
                 invert=TRUE, 
                 value=TRUE)
  
  
  date = substr(gsub('.tif','', strsplit(year_raster[1],'_')[[1]][3]),
                1,7)
  
  alpha_data <- rast(alpha_r)
  crs(alpha_data) <- crs(tiff_data)
  writeRaster(alpha_data,
              glue("C:\\Projetos\\bioflore\\gbs-kenya\\geo\\biodiv_output\\alpha_extent_{date}.tif"),
             overwrite=T)
  
  beta_data <- rast(beta_r)
  crs(beta_data) <- crs(tiff_data)
  writeRaster(beta_data,
              glue("C:\\Projetos\\bioflore\\gbs-kenya\\geo\\biodiv_output\\beta_extent_{date}.tif"),
              overwrite=T
              )
  
  spectral_data <- rast("RESULTS\\s2a_Subset\\SPCA\\SpectralSpecies\\SpectralSpecies")
  crs(spectral_data) <- crs(tiff_data)
  writeRaster(spectral_data,
              glue("C:\\Projetos\\bioflore\\gbs-kenya\\geo\\biodiv_output\\spectral_extent_{date}.tif"),
              overwrite=T
              )
  
  for (label in unique_names){
    # label <- unique_names[4]
    
    matching_polygon <- polygons[polygons$label == label,]
    
    cropped_alpha <- crop(alpha_data, as_Spatial(st_zm(matching_polygon$geom)))
    
    x.stats <- data.frame(talhao=label,
                          data=date,
                          shannon=global(cropped_alpha, fun = "mean", na.rm = TRUE))
    result_df <- rbind(result_df,x.stats)
    # Convert raster to data frame and rename columns appropriately
    raster_df_alpha <- as.data.frame(terra::as.data.frame(cropped_alpha, xy = TRUE), xy = TRUE)
    names(raster_df_alpha) <- c("x", "y", "value")
    max_shannon <- max(raster_df_alpha$value, na.rm=T)
    
    # Create a ggplot
    p <- ggplot() +
      geom_raster(data = raster_df_alpha, aes(x = x, y = y, fill = value)) +
      geom_sf(data = matching_polygon, fill = NA, color = "red", lwd = 1) +
      theme_minimal() +
      scale_fill_gradientn(limits = c(0,2), colors = scales::viridis_pal()(7)) + 
      labs(title = paste("Alpha Diversitiy for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ym")))) +
      theme( plot.title = element_text(hjust = 0.5, size = 16),
             axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
             axis.text.y = element_text(size = 8))
    
    # Save or display the plot
    ggsave(paste0("C:\\Projetos\\bioflore\\gbs-kenya\\png\\biodiv\\shannon_", label, "_", date, ".png"), plot = p)
    print(p)
    
    cropped_beta <- crop(beta_data, as_Spatial(st_zm(matching_polygon$geom)))
    
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
      ggtitle(paste("Beta Diversitiy for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ym")))) +
      theme_bw() + 
      coord_equal() +
      coord_sf() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)
      )
    ggsave(paste0("C:\\Projetos\\bioflore\\gbs-kenya\\png\\biodiv\\map_", label, "_", date, ".png"), plot = b)
    print(b)
    
    # band1 <- values(beta_data [[1]])
    # band2 <- values(beta_data [[2]])
    # band3 <- values(beta_data [[3]])
    # 
    # # Combine bands into a data frame
    # data <- data.frame(Band1 = band1, Band2 = band2, Band3 = band3)
    # 
    # # Remove any NA values
    # data <- na.omit(data)
    # 
    # ### Creating Scatter Plots
    # ##Now, let's create scatter plots comparing each pair of bands using `ggplot2`.
    # 
    # # Scatter plot for Band 1 vs Band 2
    # plot1 <- ggplot(data, aes(x = Band1, y = Band2)) +
    #   geom_point(alpha = 0.5) +
    #   labs(
    #        x = "Band 1",
    #        y = "Band 2") +
    #   theme_minimal()
    # 
    # # Scatter plot for Band 1 vs Band 3
    # plot2 <- ggplot(data, aes(x = Band1, y = Band3)) +
    #   geom_point(alpha = 0.5) +
    #   labs(
    #        x = "Band 1",
    #        y = "Band 3") +
    #   theme_minimal()
    # 
    # # Scatter plot for Band 2 vs Band 3
    # plot3 <- ggplot(data, aes(x = Band2, y = Band3)) +
    #   geom_point(alpha = 0.5) +
    #   labs(
    #        x = "Band 2",
    #        y = "Band 3") +
    #   theme_minimal()
    # 
    # g <- grid.arrange(
    #   arrangeGrob(plot1, plot2, plot3, ncol = 3),
    #   top = textGrob(paste("Scatter plot of PCoA bands for Label:", label, "Date:", gsub(" UTC","",parse_date_time(date, "ym"))),
    #                  gp = gpar(fontsize = 10))
    # )
    # 
    # ggsave(file=paste0("C:\\Projetos\\bioflore\\gbs-kenya\\png\\biodiv\\scatter_", label, "_", date, ".png"), g)
    
  }
  unlink(list.files(full.names = TRUE, recursive = TRUE), recursive = TRUE)
  date_time<-Sys.time()
  while((as.numeric(Sys.time()) - as.numeric(date_time))<60){}
}
# }

write.csv(result_df,
          file = "C:\\Projetos\\bioflore\\gbs-kenya\\spreadsheets\\result_shannon_limite.csv",
          row.names = FALSE)


# result_df <- read.csv(
#   "C:\\Projetos\\bioflore\\gbs-kenya\\spreadsheets\\result_shannon_limite.csv"
# )
# result_df <- result_df[result_df$talhao != 'vpn1-b6ru-baej',]
# result_df$data <- as.character(result_df$data)

result_df$talhao <- snakecase::to_title_case(result_df$talhao)

custom_order <- c( "Reference 1", "Reference 2", "Lower Imenti", "Gathiuru",
                  "Hombe", "Hombe 2", "Hombe 3", "Kangaita", "Untitled",
                  "Ontulili")

# Convert the 'talhao' column to a factor with the specified levels
result_df$talhao <- factor(result_df$talhao, levels = custom_order)

p <- ggplot(result_df, aes(x = talhao, y = mean, fill = data)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) +
  geom_text(aes(label = data), position = position_dodge(width = 0.9), angle = 90,
            vjust = 0.5, hjust = 2, color = "black", size = 3) +
  labs(x = 'Name', y = 'Value', title = 'Shannon index mean by date') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Omit the legend

p
ggsave(file="C:\\Projetos\\bioflore\\gbs-kenya\\png\\result_shannon.png", p)
