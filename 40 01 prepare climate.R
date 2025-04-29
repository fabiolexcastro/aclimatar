

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Vector data 
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
zne <- vect(zne)

## Raster data baseline
tmin.bsln <- as.character(dir_ls('./common_data/chirts_hist_america/monthly/Tmin')) %>% rast()
tmax.bsln <- as.character(dir_ls('./common_data/chirts_hist_america/monthly/Tmax')) %>% rast()
prec.bsln <- as.character(dir_ls('./common_data/chirps_hist_america/monthly')) %>% rast()

## Raster data future
fles.ftre <- dir_ls('./common_data/indices',  regexp = 'ssp370') %>% as.character() %>% map(., dir_ls) %>% unlist() %>% as.character()
prec.ftre <- grep('PTOT', fles.ftre, value = T)
tmin.ftre <- grep('TMIN', fles.ftre, value = T)
tmax.ftre <- grep('TMAX', fles.ftre, value = T)
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')

# Summarise baseline dataset -----------------------------------------------

## Function
avg.bsl <- function(r){
  
  avg <- map(.x = 1:12, .f = function(i){
    a <- r[[grep(paste0('-', i, '$'), names(r))]]
    a <- mean(a)
  })
  
  avg <- reduce(avg, c)
  nme <- str_sub(names(r), 1, 4)
  nme <- unique(nme)
  names(avg) <- glue('{nme}_{1:12}')
  cat('Done!\n')
  return(avg)

}

## To calcualte the average
tmin.bsln <- avg.bsl(tmin.bsln)
tmax.bsln <- avg.bsl(tmax.bsln)
prec.bsln <- avg.bsl(prec.bsln)

## Make one stack 
stck.bsln <- c(prec.bsln, tmin.bsln, tmax.bsln)

# Summarise future dataset ------------------------------------------------

## Function 
avg.ftr <- function(fls, gcm){
  
  # fls <- prec.ftre
  # gcm <- gcms[1]
  
  ## To list the files
  fls <- grep(gcm, fls, value = T) %>% 
    dir_ls() %>% 
    grep(paste0(2025:2045, collapse = '|'), ., value = T) %>% 
    as.character()
  
  nme <- dirname(fls) %>% basename() %>% unique()
  
  ## Loop by each month 
  avg <- list()
  
  for(i in 1:12){
    
    cat('>>> Month: ', i, '\n')
    mnt <- ifelse(i < 10, paste0('0', i), as.character(i))
    rst <- grep(paste0('-', mnt, '.'), fls, value = T)
    rst <- rast(rst)
    avg[[i]] <- mean(rst)
    
  }
  
  avg <- reduce(avg, c)
  names(avg) <- glue('{nme}_{1:12}')
  return(avg)
  
}

## To apply the function 
r.tx <- r.tn <- list()
for(gc in gcms){
  cat(gc,'\n')
  r.tn[[gc]] <- avg.ftr(fls = tmin.ftre, gcm = gc)   # r.pp[[gc]] <- avg.ftr(fls = prec.ftre, gcm = gc)
  r.tx[[gc]] <- avg.ftr(fls = tmax.ftre, gcm = gc)
}
r.tn <- reduce(r.tn, c)
r.tx <- reduce(r.tx, c)

## Average 
a.tn <- a.tx <- list()
for(i in 1:12){
  a.tn[[i]] <- mean(r.tn[[grep(paste0('_', i, '$'), names(r.tn))]])
  a.tx[[i]] <- mean(r.tx[[grep(paste0('_', i, '$'), names(r.tx))]])
}
a.tn <- reduce(a.tn, c)
a.tx <- reduce(a.tx, c)

names(a.tn) <- glue('tmin_{c(paste0("0", 1:9), 10:12)}_ftr')
names(a.tx) <- glue('tmax_{c(paste0("0", 1:9), 10:12)}_ftr')

stck.ftre <- c(a.tn, a.tx)

# Min - max (prec) --------------------------------------------------------
prec.ftre <- prec.ftre %>% map(dir_ls) %>% unlist() %>% as.character() %>% grep(paste0(2025:2055, collapse = '|'), ., value = T)
prec.gcms <- map(.x = 1:length(gcms), .f = function(g){
  
  ## Filter the GCMe
  cat('>>> ', gcms[g], '\n')
  gcme <- gcms[g]
  prec <- grep(gcme, prec.ftre, value = T)
  mnts <- c(paste0('0', 1:9), 10:12)
  
  ## To make the average
  rstr <- map(.x = mnts, .f = function(m){
    
    cat(m, '\n')
    rst <- rast(grep(paste0(m, '.tif$'), prec, value = T))
    rst <- mean(rst)
    names(rst) <- glue('prec_{m}')
    return(rst)
    
  })
  
  ## Reduce the raster
  rstr <- reduce(rstr, c)
  return(rstr)
  
})
prec.gcms <- reduce(prec.gcms, c)
prec.sttn <- map(.x = 1:12, .f = function(m){
  mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
  ppt <- prec.gcms[[grep(mnt, names(prec.gcms))]]
  min <- min(ppt)
  max <- max(ppt)
  avg <- mean(ppt)
  names(min) <- glue('prec_{mnt}_min-ftr')
  names(avg) <- glue('prec_{mnt}_avg-ftr')
  names(max) <- glue('prec_{mnt}_max-ftr')
  return(c(min, avg, max)) 
}) 
prec.sttn <- reduce(prec.sttn, c)
prec.sttn
stck.ftre <- c(stck.ftre, prec.sttn)

# To write the raster -----------------------------------------------------
'./common_data/climate'
terra::writeRaster(x = stck.bsln, filename = './common_data/climate/clima_bsl.tif', overwrite = TRUE)
terra::writeRaster(x = stck.ftre, filename = './common_data/climate/clima_ftr.tif', overwrite = TRUE)

stck.bsln <- terra::rast('./common_data/climate/clima_bsl.tif')
stck.ftre <- terra::rast('./common_data/climate/clima_ftr.tif')

## To change the names 
stck.bsln
names(stck.bsln) <- c(
  c(glue('prec_0{1:9}_bsl'), glue('prec_{10:12}_bsl')),
  c(glue('tmin_0{1:9}_bsl'), glue('tmin_{10:12}_bsl')),
  c(glue('tmax_0{1:9}_bsl'), glue('tmax_{10:12}_bsl'))
  
)

stck.clma <- c(stck.bsln, stck.ftre)

# Forecast datasets -------------------------------------------------------

## Southamerica
rstr.sta <- rast('//catalogue/workspace-cluster9/CL/v1/data/tif/results/climate_indx-raw_ssp370_col-ecu-per_v2.tif')
rstr.sta <- rstr.sta[[grep('frc', names(rstr.sta))]]

## Centralamerica
fles.cam <- as.character(dir_ls('//catalogue/workspace-cluster9/CAM/v1/data/tif/climate/forecast'))
rstr.cam <- rast(grep('stack', fles.cam, value = T))
names(rstr.sta) <- sub("_", "-", names(rstr.sta), fixed = TRUE)

# Check the variables -----------------------------------------------------
names(rstr.sta) %>% length()
names(rstr.cam) %>% length()
rstr.cam <- rstr.cam[[grep(paste0(names(rstr.sta), collapse = '|'), names(rstr.cam))]]
nmes <- names(rstr.cam) 

# Mosaic the variables ----------------------------------------------------
rstr <- map(.x = 1:length(nmes), .f = function(i){
  
  cat(i, '\n')
  r.cam <- rstr.cam[[grep(paste0(nmes[i], '$'), names(rstr.cam), value = F)]]
  r.sta <- rstr.sta[[grep(paste0(nmes[i], '$'), names(rstr.sta), value = F)]]
  r.lst <- sprc(r.cam, r.sta)
  r.msc <- mosaic(r.lst)
  return(r.msc)
  
}) %>% reduce(., c)

rstr <- rstr[[1:48]]
mnts <- c(paste0('0', 1:9), 10:12)

## To change the names forecast 
names(rstr) <- c(glue('prec_{mnts}_frc'), glue('tmin_{mnts}_frc'), glue('tavg_{mnts}_frc'), glue('tmax_{mnts}_frc'))

# Join climate  -----------------------------------------------------------
base <- read_csv('./tble/base/clima-index_ken-cff.csv', show_col_types = FALSE)

## Baseline / Forecast / Future
stck.clma <- c(stck.clma, rstr)
names(stck.clma)[1:36] <- names(stck.bsln)

terra::writeRaster(x = stck.clma, filename = './common_data/climate/clima_bsl-ftr_prec-tmin-tavg-tmax.tif', overwrite = TRUE)

stck.clma <- terra::rast('./common_data/climate/clima_bsl-ftr_prec-tmin-tavg-tmax-etps-baln.tif')
