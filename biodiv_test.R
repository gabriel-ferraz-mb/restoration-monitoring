# install.packages("remotes")
# remotes::install_github('cran/dissUtils')
# remotes::install_github('jbferet/biodivMapR')

library(biodivMapR)
library(terra)
library(glue)

setwd("C:\\Projetos\\bioflore\\restoration_monitoring\\BiodivMapR")

############################## Start Process ###################################
raster_files <-  grep(list.files(path="C:\\Projetos\\bioflore\\bgci_II\\geo\\biodiv_input",
                               full.names = T,recursive = T),
                  pattern = '.xml', invert =T, value = T)
out <- raster_files[lengths(raster_files)!=0]
# matching_polygon <- polygons[polygons$Name == label,]
# matching_polygon <- polygons
years <- seq(2019, 2024, by = 1)

talhoes <- list(
  "ABP",
  "Bigodi",
  "Dilips",
  "IPEECOSIA",
  "Jubiya",
  "KB",
  "KOREN",
  "LUM",
  # "Magaca",
  "MARLENE",
  # "Nteko",
  "Pandal",
  "RPPNRENOPOLIS",
  "Soraypampa",
  "Wangari"
)

for (talhao in talhoes){
  # talhao = talhoes[1]
  label_raster <- out[grepl(talhao, out)]
  for (year in years) {
    # year = years[1]
    
    year_raster <- label_raster[grepl(as.character(year), label_raster)]
    
    date = gsub('.tif','', strsplit(year_raster[1],'_')[[1]][4])
    
    if (file.exists(
      glue("C:\\Projetos\\bioflore\\bgci_II\\geo\\biodiv_output\\{talhao}_alpha_extent_{date}.tif"))
      ){
      print("Process already completed")
      unlink(list.files(full.names = TRUE, recursive = TRUE), recursive = TRUE)
      next 
    } else {
      print(glue("starting {talhao}-{date}"))
    }
    
    
    if (length(year_raster) > 1){
      
      raster_list <- lapply(year_raster, function(file) {
        original_raster <- rast(file)
        reprojected_raster <- project(original_raster, "EPSG:4326", res = c(1/11113.95, 1/11113.95))
        return(reprojected_raster)
      })
      
      # Use a reference raster for alignment purposes
      s <- sprc(raster_list)
      # Use mosaic function to combine the aligned rasters
      imgMosaic <- merge(s)
      
    } else {
      imgMosaic <- rast(year_raster[1])
    }
    
    # bnames <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11',
    #             'B12','AOT','WVP','SCL','TCI_R','TCI_G','TCI_B','MSK_CLDPRB',
    #             'MSK_SNWPRB','QA10','QA20','QA60','MSK_CLASSI_OPAQUE',
    #             'MSK_CLASSI_CIRRUS','MSK_CLASSI_SNOW_ICE')
    
    bnames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL')
    names(imgMosaic) <- bnames
    
    tiff_data  <- imgMosaic[[c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')]]
    
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
                                                     Blue_Thresh = 500,
                                                     NIR_Thresh = 1500)
    
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
    nbCPU <- 8
    MaxRAM <- 0.5
    nbclusters <- 10
    
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
    
    print("Done!")
    
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
    
    alpha_data <- rast(alpha_r)
    crs(alpha_data) <- crs(tiff_data)
    writeRaster(alpha_data,
                glue("C:\\Projetos\\bioflore\\bgci_II\\geo\\biodiv_output\\{talhao}_alpha_extent_{date}.tif"),
                overwrite=T)
    
    beta_data <- rast(beta_r)
    crs(beta_data) <- crs(tiff_data)
    writeRaster(beta_data,
                glue("C:\\Projetos\\bioflore\\bgci_II\\geo\\biodiv_output\\{talhao}_beta_extent_{date}.tif"),
                overwrite=T
    )
    
    spectral_data <- rast("RESULTS\\s2a_Subset\\SPCA\\SpectralSpecies\\SpectralSpecies")
    crs(spectral_data) <- crs(tiff_data)
    writeRaster(spectral_data,
                glue("C:\\Projetos\\bioflore\\bgci_II\\geo\\biodiv_output\\{talhao}_spectral_extent_{date}.tif"),
                overwrite=T
    )
    
    
    unlink(list.files(full.names = TRUE, recursive = TRUE), recursive = TRUE)
    date_time<-Sys.time()
    while((as.numeric(Sys.time()) - as.numeric(date_time))<60){}
  }
}
