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

suppressWarnings(expr)

# Agora, diga ao R em qual diret?rio a imagem est?. Lembre-se de inserir "\\" entre o nome das pastas
setwd('C:\\Projetos\\bioflore\\ecosia')

files <- list.files(path="PSdata", pattern='ESP', full.names = T)
files <- lapply(files, grep, pattern = "QA.tif", invert = TRUE, value = TRUE)
out <- files[lengths(files)!=0]

# data<-brick(out[[1]])
# m<-brick(str_replace(out[[1]], ".tif", "_QA.tif"))[[1]]
# values(m) <- 0
# imgBrick <- data*m
# contains_non_zero <- any(values(imgBrick) != 0, na.rm = TRUE)  # `values(r)` extracts the raster data as a vector

# 1) Extraia a resposta espectral de pixels "puros" da imagem e construa uma "biblioteca espectral", veja
# o exemplo abaixo para pixels de vegeta??o, ?gua e solo

# Crie biblioteca espectral no formato que a fun??o 'mesma' aceita, ou seja, linhas representam
# um ?nico endmember, as colunas s?o as bandas espectrais.

endmembers<-read.csv("em.txt")
colnames(endmembers)<-c("class","Blue","Green","Red","NIR")

melted <- melt(endmembers, id.vars="class")
colnames(melted)<-c("Class","WaveLength","Reflectance")
p <- ggplot(melted, aes(WaveLength, Reflectance, color = Class)) +
  geom_line(aes(group=Class)) + ggtitle("Ipê Reflectance by class")
p

# out = out[1:2]

for (name in out){
  
  product <- str_replace_all(name, c(PSdata="NDFI", .tif="_NDFI.tif"))
  
  if (file.exists(product)){
    next
  }

# Carregue uma imagem do sat?lite Landsat-8 adqurida sobre a cidade do Rio de Janeiro - RJ
  data<-brick(name)
  m<-brick(str_replace(name, ".tif", "_QA.tif"))[[1]]
  imgBrick <- data*m
  
  contains_non_zero <- any(values(imgBrick) != 0, na.rm = TRUE)  # `values(r)` extracts the raster data as a vector
  
  if(!contains_non_zero){
    next
  }
  
  # Altere o nome das bandas
  imgBrick@data@names<-c("Blue","Green","Red","NIR")
  
  # imgRGBcv<-stack(imgBrick$Red,imgBrick$Green,imgBrick$Blue)
  # 
  # maxBands<-max(getValues(imgRGBcv)) # A fun??o 'plotRGB' necessita saber qual o valor m?ximo dos pixels nas bandas, j? que a imagem n?o ? 8-bits (0-255)
  # 
  # plotRGB(imgRGBcv, scale=maxBands, stretch='lin')
  
  # 2) Aplique a fun??o 'mesma' para obter as imagens fra??o.
  imgFracao<-mesma(imgBrick, endmembers)
  
  gv_shade <- imgFracao$veg/(1-imgFracao$agua)
  NDFI <- (gv_shade - (imgFracao$palha+imgFracao$areia))/(imgFracao$palha+imgFracao$areia+gv_shade)
  # plot(NDFI,col=gray.colors(256))
  writeRaster(NDFI,product)
}

ndfi.files <- list.files(path="NDFI", pattern='ESP', full.names = T)

result <- data.frame()
for (name in ndfi.files){
  #name <- ndfi.files[[i]]
  date <- as.Date(str_split(name, '_')[[1]][3], "%Y%m%d")
  talhao <- substring(name,6,11)
  data <- brick(name)[[1]]
  x.stats <- data.frame(talhao=talhao,
    data=date,
    ndfi=cellStats(data, "mean"))
  
  result <- rbind(result,x.stats)
}

labels <- as.list(unique(result$talhao))
for (label in labels){
  
  #label <- labels[[1]]
  #png(paste0(label, ".png"))
  if (startsWith(label, "IPE")){
    RE = 'IPE_RE'
  } else{
    RE = 'ESP_RE'
  }
  
  vc <- c(label,RE)
  data<-result[result$talhao %in% vc,]
  data <- data[!duplicated(data[c('talhao','data')]),]
  
  pivoted_data <- data %>% pivot_wider(
    names_from = talhao,
    values_from = ndfi
  )
  # Prepare data for plotting
  pivoted_data$data <- as.Date(pivoted_data$data) # Ensure dates are in Date format
  pivoted_data <- pivoted_data[order(pivoted_data$data),]
  
  date_non_na <- pivoted_data %>% filter(!is.na(pivoted_data[[label]])) %>% slice(1) %>% pull(data)
  #pivoted_data %>% filter(date >= date_non_na)
  # Define the specific range for the x-axis
  start_date <- pivoted_data[1,1]$data 
  end_date <- pivoted_data[nrow(pivoted_data),1]$data
  
  # Open a new graphics window (specific to Windows)
  index_name = "NDFI"
  #win.graph(10.90625, 16.83333)
  
  #jpeg(filename=paste0(label,".png"))
  # Create an empty plot with x-axis labels suppressed
  plot(1, type = "n", xlab = "", ylab = paste0(index_name), xlim = c(start_date, end_date),
       ylim = c(-1.5, 1.5), main = label, xaxt = 'n',cex.lab = 1.5)
  
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
  y <- pivoted_data$data
  #y_re <- IPE_RE$date
  degree <- 15 # Degree of the spline
  df <- 5 # Degrees of freedom or roughly the number of knots
  
  # Fit B-spline
  fit_bs <- lm(x ~ bs(y, degree = degree, df = df))
  fit_bs_re <- lm(x_re ~ bs(y, degree = degree, df = df))
  # Predictions for a smooth curve across the entire date range
  y_pred <- seq(date_non_na, end_date, by = "day")
  x_pred <- predict(fit_bs, newdata = list(y = y_pred))
  
  # x_pred[x_pred< -1] <- -1
  # x_pred[x_pred> 1] <- 1
  
  y_pred_re <- seq(start_date, end_date, by = "day")
  x_pred_re <- predict(fit_bs_re, newdata = list(y = y_pred_re))
  
  # x_pred_re[x_pred_re< -1] <- -1
  # x_pred_re[x_pred_re> 1] <- 1
  
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
file.copy(from=plots.png.paths, to="C:/Projetos/bioflore/ecosia/plots_ndfi")
