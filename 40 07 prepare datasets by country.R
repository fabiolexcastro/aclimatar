


# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)


# Load data ---------------------------------------------------------------

## Raster data
fles <- dir_ls('./common_data/climate_indices', regexp = '.tif$')
ccoa <- rast('./common_data/climate_indices/clima_index_cocoa.tif')
cffe <- rast('./common_data/climate_indices/clima_index_coffee.tif')

## Vector data 
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
zne <- vect(zne)

# Function to extract by mask and prepare table ---------------------------
make.table <- function(rst, iso, crp){
  
  # rst <- ccoa
  # iso <- 'COL'
  # crp <- 'cocoa'
  
  ## To extract by mask
  cat('To start the process!\n')
  cnt <- zne[zne$sov_a3 == iso,]
  rst <- terra::crop(rst, cnt)
  rst <- terra::mask(rst, cnt)
  
  ## Raster to table
  tbl <- rst %>% 
    terra::as.data.frame(xy = T) %>% 
    as_tibble()
  
  ## To write the table 
  out <- glue('./common_data/tables_send1')
  write.csv(tbl, paste0(out, '/', crp, '_', iso, '.csv'), row.names = FALSE)
  cat('Done!\n')
  
}

# To apply the function ---------------------------------------------------
isos <- c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI')

## Cocoa
for(i in 1:length(isos)){
  cat(isos[i], '\t')
  make.table(rst = ccoa, iso = isos[i], crp = 'cocoa')
}

## Coffee
for(i in 1:length(isos)){
  cat(isos[i], '\t')
  make.table(rst = cffe, iso = isos[i], crp = 'coffee')
}





