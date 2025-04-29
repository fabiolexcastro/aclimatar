
# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Vector data -------------------------------------------------------------
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
zne <- vect(zne)

base <- read_csv('./tble/base/clima-index_ken-cff.csv', show_col_types = FALSE)
names(base)

# Raster data -------------------------------------------------------------
rstr <- terra::rast('./common_data/climate/clima_bsl-ftr_prec-tmin-tavg-tmax.tif')

rstr.bsln <- rstr[[grep('bsl', names(rstr))]]
rstr.frcs <- rstr[[grep('frc', names(rstr))]]
rstr.ftre <- rstr[[grep('ftr', names(rstr))]]

nlyr(rstr.bsln) + nlyr(rstr.frcs) + nlyr(rstr.ftre) == nlyr(rstr)

rstr.bsln

## Solar radiation 
srad <- dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad') %>% mixedsort() %>% grep('et_solrad', ., value = T) %>% as.character() %>% rast() %>% crop(., zne) %>% mask(., zne)

# To calculate the variables ----------------------------------------------

##
calc.vrs <- function(rst, tme){
  
  # rst <- rstr.ftre
  # tme <- 'ftr'
  
  cat('To start!\n')
  
  if(tme == 'ftr'){
    ppt <- rst[[grep('prec', names(rst))]]  
    ppt <- ppt[[grep('avg-ftr', names(ppt))]]
  } else {
    ppt <- rst[[grep('prec', names(rst))]]
  }
  
  tmn <- rst[[grep('tmin', names(rst))]]
  tmx <- rst[[grep('tmax', names(rst))]]
  
  ## 
  tav <- (tmn + tmx) / 2
  names(tav) <- gsub('tmin', 'tavg', names(tav))
  
  ##
  srd <- terra::resample(srad, tav, method = 'bilinear')
  
  ## 
  etp <- 0.0013 * 0.408 * srd * (tav + 17) * (tmx - tmn - 0.0123 * ppt) ^ 0.76
  etp <- etp * c(31,29,31,30,31,30,31,31,30,31,30,31)
  for(i in 1:12){etp[[i]][which.lyr(is.na(etp[[i]]))] <- 0}
  etp <- terra::crop(etp, zne) %>% terra::mask(., zne)
  names(etp) <- glue('etps_{c(paste0("0", 1:9), 10:12)}_{tme}')
  
  ## 
  bal <- ppt - etp
  names(bal) <- glue('baln_{c(paste0("0", 1:9), 10:12)}_{tme}')
  
  ## 
  mtx <- matrix(c(-Inf, -10, -1, -10, 10, 0, 0, Inf, 1), ncol = 3, nrow = 3, byrow = T)
  rcl <- terra::classify(bal, mtx)
  names(rcl) <- glue('baln-rclf_{c(paste0("0", 1:9), 10:12)}_{tme}')
  
  ##
  if(tme == 'ftr'){
    ppt <- rst[[grep('prec', names(rst))]]  
    ppt.min <- ppt[[grep('min-ftr', names(ppt))]]
    ppt.max <- ppt[[grep('max-ftr', names(ppt))]]
    ppt.avg <- ppt[[grep('avg-ftr', names(ppt))]]
    stk <- c(ppt.min, ppt.avg, ppt.max, tmn, tmx, tav, etp, bal, rcl)
    
  } else {
    stk <- c(ppt, tmn, tmx, tav, etp, bal, rcl)
  }
  
  ##
  cat('Done!\n')
  return(stk)
  
}

##
stck.bsl <- calc.vrs(rst = rstr.bsln, tme = 'bsl')
stck.frc <- calc.vrs(rst = rstr.frcs, tme = 'frc')
stck.ftr <- calc.vrs(rst = rstr.ftre, tme = 'ftr')

## 
stck <- c(stck.bsl, stck.frc, stck.ftr)
terra::writeRaster(x = stck, filename = './common_data/climate/clima_bsl-ftr_prec-tmin-tavg-tmax-etps-baln.tif', overwrite = TRUE)


