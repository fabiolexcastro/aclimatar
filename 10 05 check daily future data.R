

# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

root <- './common_data/indices'
dirs <- as.character(dir_ls(root))
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
ssps <- c('ssp245', 'ssp370')

## 
mskr <- dir_ls('./common_data/chirps_cmip6_america/Prec_ACCESS-ESM1-5_ssp245_2021_2040')[1] %>% rast()
mskr <- mskr * 0 + 1

## 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- filter(wrld, sov_a3 %in% c('SLV', 'GTM', 'MEX', 'NIC', 'HND', 'PAN', 'CRI', 'COL', 'ECU', 'PER'))

# Sample points -----------------------------------------------------------
pnts <- map_dfr(.x = 1:nrow(zone), .f = function(i){
  zne <- zone[i,]
  zne <- vect(zne)
  msk <- terra::crop(mskr, zne)
  msk <- terra::mask(msk, zne)
  pnt <- as.data.frame(raptr::randomPoints(mask = msk, n = 1))
  pnt <- mutate(pnt, iso = zne$sov_a3)
  return(pnt)
})
write.csv(pnts, './tble/points_sample.csv', row.names = FALSE)

# Precipitation -----------------------------------------------------------

check.prec <- function(sspe){
  
  ## 
  sspe <- 'ssp245'
  
  ## To start
  cat('To process: ', sspe, '\n')
  pths <- as.character(dir_ls(grep(sspe, dirs, value = T)))
  pths <- grep(sspe, pths, value = T)
  
  map(.x = 1:length(gcms), .f = function(i){
    
    drs <- grep(gcms[i], pths, value = T)
    fls <- map(drs, dir_ls)
    fls <- unlist(fls)
    fls <- as.character(fls)
    
    ## To reas as raster
    
  })
  
  
}