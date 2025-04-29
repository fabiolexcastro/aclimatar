
# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
clma <- rast('./common_data/climate/clima_bsl-ftr_prec-tmin-tavg-tmax-etps-baln.tif')
indx.ccoa <- rast('./common_data/indices_global/ssp370/index_rcl-cocoa.tif')
indx.cffe <- rast('./common_data/indices_global/ssp370/index_rcl-coffee.tif')

# Check baseline vars -----------------------------------------------------
prec <- clma[[grep('prec', names(clma), value = F)]]
tmin <- clma[[grep('tmin', names(clma), value = F)]]
tmax <- clma[[grep('tmax', names(clma), value = F)]]

# Stacking ----------------------------------------------------------------
stck.ccoa <- c(clma, indx.ccoa)
stck.cffe <- c(clma, indx.cffe)

# To check NAs ------------------------------------------------------------
rmve.nas <- function(stck){
  
  ## Raster to table
  cat('To start the analysis!\n')
  tble <- terra::as.data.frame(stck, xy = T)
  tble <- as_tibble(tble)
  
  ## To remove NAs
  nrow(tble)
  dfrm <- drop_na(tble)
  nrow(dfrm)
  
  ## Table to raster 
  rstr <- rast(dfrm, type = 'xyz', crs = 'EPSG:4326')
  
  ## Finish 
  cat('Done!\n')
  return(rstr)
  
  
}

## To apply the function 
rstr.ccoa <- rmve.nas(stck = stck.ccoa)
rstr.cffe <- rmve.nas(stck = stck.cffe)

# To write the rasters ----------------------------------------------------
terra::writeRaster(x = rstr.ccoa, filename = './common_data/climate_indices/clima_index_cocoa.tif', overwrite = TRUE)
terra::writeRaster(x = rstr.cffe, filename = './common_data/climate_indices/clima_index_coffee.tif', overwrite = TRUE)



