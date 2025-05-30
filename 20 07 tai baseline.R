## Thornthwaite's Aridity Index (TAI)
## By: H. Achicanoy
## December, 2022
## Eddited: F Castro-Llanos
## April, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,raster,lubridate,envirem))
#source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- './common_data'
ref <- terra::rast('./common_data/chirps_hist_america/monthly/Prec_1995.tif')[[1]]
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'BLZ', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI', 'DOM'))
zne <- vect(zne)

# sce_climate <- "historical" #"future"

# ET SRAD, load and process only once
srf <- list.dirs(paste0('//catalogue/workspace-cluster9/DATA/ET_SolRad'), full.names = T, recursive = F)
srf <- srf[-length(srf)]
srf <- srf %>% gtools::mixedsort()
srd <- srf %>% rast()
srd <- srd %>% terra::resample(ref)
srd <- srd %>% terra::crop(zne)
srd <- srd %>% terra::mask(zne)
names(srd) <- c(paste0('SRAD_0',1:9),paste0('SRAD_', 10:12))

# Calculate TAI function
calc_tai <- function(yr, dataset="CHIRPS"){
  
  outfile <- paste0(out_dir,'/TAI-',yr,'.tif')
  cat(outfile, "\n")
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
    dts <- paste0(yr, '-', 1:12)
    
    # Files
    pr_fls <- paste0(pr_pth,'/daily/Prec_',dts,'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    
    if (dataset == "CHIRPS") {
      tx_fls <- paste0(tx_pth,'/','Tmax_',gsub(pattern='-', replacement='-', x=dts, fixed=T),'.tif')
      tm_fls <- paste0(tm_pth,'/','Tmin_',gsub(pattern='-', replacement='-', x=dts, fixed=T),'.tif')
    } else if (dataset == "AgERA5") {
      tx_fls <- paste0(ae5tx_pth,'/Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
      tm_fls <- paste0(ae5tm_pth,'/Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
    }
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- tm_fls[file.exists(tm_fls)]
    
    # Read variables, and calculate monthly summaries, do by month to reduce memory consumption
    prc_ls <- tav_ls <- rng_ls <- c()
    
    system.time(
      expr = {
        for (j in 1:12) {
          cat("month=", j, "\n")
          # Load precip
          mth_fls <- paste0(pr_pth, '/daily/Prec_', yr, '-', j, '.tif')
          prc <- terra::rast(pr_fls[pr_fls %in% mth_fls])
          # prc <- prc %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
          # prc[prc < 0] <- NA
          prc_month <- sum(prc, na.rm=TRUE)
          rm(prc)
          gc(verbose=F, full=T, reset=T)
          
          # Load tmax
          if (dataset == "CHIRPS") {
            mth_fls <- paste0(tx_pth, '/', 'Tmax_', yr, '-', j, '.tif')
          } else if (dataset == "AgERA5") {
            mth_fls <- paste0(ae5tx_pth, '/Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_', yr, sprintf("%02.0f", j), sprintf("%02.0f", 1:31), '_final-v1.0.nc')
          }
          tmx <- terra::rast(tx_fls[tx_fls %in% mth_fls])
          if (dataset == "AgERA5") {
            tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
            tmx <- tmx - 273.15
          } else {
            tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
            tmx[tmx == -9999] <- NA
          }
          
          #Load tmin
          if (dataset == "CHIRPS") {
            mth_fls <- paste0(tm_pth, '/', 'Tmin_', yr, '-', j, '.tif')
          } else if (dataset == "AgERA5") {
            mth_fls <- paste0(ae5tm_pth, '/Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_', yr, sprintf("%02.0f", j), sprintf("%02.0f", 1:31), '_final-v1.0.nc')
          }
          tmn <- terra::rast(tm_fls[tm_fls %in% mth_fls])
          if (dataset == "AgERA5") {
            # tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
            # tmn <- tmn - 273.15
          } else {
            # tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
            # tmn[tmn == -9999] <- NA
          }
          
          # Calculate average temperature
          tav <- (tmx + tmn)/2
          tav_month <- mean(tav, na.rm=TRUE)
          rm(tav)
          gc(verbose=F, full=T, reset=T)
          
          # Calculate temperature range
          rnge <- abs(tmx - tmn)
          rng_month <- mean(rnge, na.rm=TRUE)
          rm(rnge, tmx, tmn)
          gc(verbose=F, full=T, reset=T)
          
          # Append
          prc_ls <- c(prc_ls, prc_month)
          tav_ls <- c(tav_ls, tav_month)
          rng_ls <- c(rng_ls, rng_month)
        }
      }
    )
    
    
    # Assign precipitation names in envirem environment
    envirem::assignNames(solrad='SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
    
    TMEAN <- tav_ls %>% terra::rast() #%>% raster::stack()
    PREC  <- prc_ls %>% terra::rast() #%>% raster::stack()
    TRNG  <- rng_ls %>% terra::rast() #%>% raster::stack()
    
    names(TMEAN) <- c(paste0('TMEAN_0',1:9),paste0('TMEAN_', 10:12))
    names(PREC)  <- c(paste0('PREC_0',1:9),paste0('PREC_', 10:12))
    names(TRNG)  <- c(paste0('TRNG_0',1:9),paste0('TRNG_', 10:12))
    
    #clean up
    rm(tav_ls, prc_ls, rng_ls, tav_month, prc_month, rng_month)
    gc(verbose=F, full=T, reset=T)
    
    # Thornthwaite's Aridity Index
    PET <- envirem::monthlyPET(TMEAN, srd, TRNG)# %>% raster::stack()
    names(PET)  <- c(paste0('PET_0',1:9), paste0('PET_', 10:12))
    PET <- terra::resample(PET, PREC[[1]])
    TAI <- envirem::aridityIndexThornthwaite(PREC, PET)
    names(TAI) <- yr
    
    #write output
    terra::writeRaster(TAI, outfile)
    
    #clean up
    rm(PET, TMEAN, PREC, TRNG, TAI)
    gc(verbose=F, full=T, reset=T)
    
  }
}

#runs
if (sce_climate == "historical") {
  # Historical setup
  stp <- data.frame(yrs = 1995:2014) #1995:2014
  pr_pth <- paste0(root,'/chirps_hist_america') # Precipitation
  tm_pth <- paste0(root,'/chirts_hist_america/daily/Tmin') # Minimum temperature
  tx_pth <- paste0(root,'/chirts_hist_america/daily/Tmax') # Maximum temperature
  # ae5tm_pth <- paste0(root,'/ecmwf_agera5/2m_temperature-24_hour_minimum') # Minimum temperature
  # ae5tx_pth <- paste0(root,'/ecmwf_agera5/2m_temperature-24_hour_maximum') # Maximum temperature
  out_dir <- paste0(root,'/indices/historical/TAI')
  1:nrow(stp) %>%
    purrr::map(.f = function(i){
      calc_tai(yr = stp$yrs[i], dataset="CHIRPS")
      gc(verbose=F, full=T, reset=T)
      tmpfls <- list.files(tempdir(), full.names=TRUE)
      1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
    })
} else if (sce_climate == "future") {
  # Future setup
  #gcm <- 'EC-Earth3' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
  # gcm <- 'ACCESS-ESM1-5'
  # ssp <- 'ssp245'
  # prd <- '2021_2040'
  
  for (gcm in c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")) {
    for (ssp in c('ssp370')) {
      for (prd in c('2021_2040')) {
        cmb <- paste0(ssp,'_',gcm,'_',prd)
        prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
        stp <- data.frame(yrs = prd_num[1]:prd_num[2])
        pr_pth <- paste0(root,'/chirps_cmip6_america/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
        tm_pth <- paste0(root,'/chirts_cmip6_america/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
        tx_pth <- paste0(root,'/chirts_cmip6_america/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
        out_dir <- paste0(root,'/indices/',ssp, '_', gcm, '/TAI')
        1:nrow(stp) %>%
          purrr::map(.f = function(i){calc_tai(yr = stp$yrs[i])})
      }
    }
  }
}


# # ----------------------------------------------------------------------
# # Data fixes
# # Get reruns file.
# source("~/Repositories/hazards/R/03_bias_correction/getReruns.R")
# other_bfiles <- c(paste0(root, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2043.03.18.tif"), 
#                   paste0(root, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2054.10.10.tif"),
#                   paste0(root, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2056.03.29.tif"))
# reruns_df <- getReruns(newfiles=other_bfiles) %>%
#   dplyr::select(-varname, -mn) %>%
#   unique(.)
# 
# # Do the reruns
# for (j in 1:nrow(reruns_df)) {
#   gcm <- reruns_df$gcm[j]
#   ssp <- reruns_df$ssp[j]
#   prd <- reruns_df$prd[j]
#   cmb <- paste0(ssp,'_',gcm,'_',prd)
#   prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
#   stp <- data.frame(yrs = prd_num[1]:prd_num[2])
#   
#   #folders
#   pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
#   tm_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
#   tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
#   out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/TAI')
#   
#   #remove and redo file
#   this_file <- paste0(out_dir, "/TAI-", reruns_df$yr[j], ".tif")
#   cat("redoing file", this_file, "\n")
#   system(paste0("rm -f ", this_file))
#   system(paste0("Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ", gcm, " ", ssp, " ", prd, " 1 ", which(stp$yr %in% as.numeric(reruns_df$yr[j]))))
#   #calc_tai(yr = reruns_df$yr[j])
# }
# 