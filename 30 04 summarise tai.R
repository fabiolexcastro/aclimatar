
# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Vector data 
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'BLZ', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI', 'DOM'))
zne <- vect(zne)

## Directories
dirs <- dir_ls('./common_data/indices')

dirs.bsln <- grep('historical', dirs, value = T) %>% 
  dir_ls() %>% 
  grep('TAI', ., value = T) %>% 
  as.character()

dirs.ftre <- dirs[-grep('historical', dirs, value = F)] %>% 
  map(., dir_ls) %>% 
  unlist() %>% 
  as.character() %>% 
  grep('ssp370', ., value = T) %>% 
  grep('TAI', ., value = T)

gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')

# Summarise baseline ------------------------------------------------------

##
fles.bsln <- dirs.bsln %>% dir_ls(., regexp = '.tif$') %>% as.character()
year.bsln <- basename(fles.bsln) %>% str_split('-') %>% map_chr(2) %>% unique()
year.bsln <- 1995:2014

##
rstr.bsln <- map(.x = 1:length(year.bsln), .f = function(i){
  
  cat('To process: ', year.bsln[i], '\n')
  file <- grep(year.bsln[i], fles.bsln, value = T)
  rstr <- rast(file)
  rstr <- mean(rstr)
  rstr <- terra::crop(rstr, zne)
  rstr <- terra::mask(rstr, zne)  
  names(rstr) <- glue('TAI_{year.bsln[i]}')
  return(rstr)
  
}) %>% 
  reduce(c)

##
rstr.bsln <- mean(rstr.bsln)
names(rstr.bsln) <- glue('TAI_bsl')

##
terra::writeRaster(x = rstr.bsln, filename = glue('./common_data/indices_average/historical/TAI/TAI_raw.tif'), overwrite = TRUE)

# Summarise future --------------------------------------------------------
dirs.ftre
sspe <- 'ssp370'

##
smmr.gcme <- function(gcme){
  
  ##
  cat('To process: ', gcme, '\n')
  fles <- as.character(dir_ls(grep(gcme, dirs.ftre, value = T)))
  # year <- basename(fles) %>% str_split('-') %>% map_chr(2) %>% unique()
  year <- 2025:2055
  
  rstr <- map(.x = 1:length(year), .f = function(i){
    
    cat('To process: ', year[i], '\n')
    file <- grep(year[i], fles, value = T)
    rstr <- rast(file)
    rstr <- mean(rstr)
    rstr <- terra::crop(rstr, zne)
    rstr <- terra::mask(rstr, zne)  
    names(rstr) <- glue('TAI_{year[i]}')
    return(rstr)
    
  }) %>% 
    reduce(c)
  
  ##
  rstr <- mean(rstr)
  names(rstr) <- glue('TAI_{sspe}_{gcme}')
  
  ##
  dout <- glue('./common_data/indices_average/{sspe}_{gcme}/TAI')
  dir_create(dout)
  terra::writeRaster(x = rstr, filename = glue('./common_data/indices_average/{sspe}_{gcme}/TAI/TAI_raw.tif'), overwrite = TRUE)
  cat('Done!\n')
  return(rstr)
  
}

## 
rstr.ftre <- map(gcms, smmr.gcme)
rstr.ftre <- reduce(rstr.ftre, c)

## 
rstr.ftre.avrg <- mean(rstr.ftre)
names(rstr.ftre.avrg) <- glue('TAI_{sspe}')
dir_create(glue('./common_data/indices_global/{sspe}'))
terra::writeRaster(x = rstr.ftre.avrg, filename = glue('./common_data/indices_global/{sspe}/TAI_raw.tif'))



