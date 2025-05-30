
## Get daily future data
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,sf,lubridate,furrr,rnaturalearthdata,rnaturalearth))

root <- './common_data' # root <- '/home/jovyan/common_data'

# Worldataset
wrld <- (rnaturalearth::ne_countries(returnclass = 'sf', scale = 50))
zone <- filter(wrld, sov_a3 %in% c('COL', 'ECU', 'BLZ', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI', 'DOM'))

# America reference raster
lons=c(-179.0625, -56.8760); lats=c(-56.8760, 85.6250)
ref <- terra::rast('//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps/chirps-v2.0.1981.01.01.tif') #ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))
ref <- terra::crop(ref, c(lons, lats))
ref <- terra::crop(ref, zone)
ref <- terra::mask(ref, zone)

# Interpolated monthly anomalies directory
anm_pth <- paste0(root,'/esfg_cmip6/intermediate/interpolated_mthly_anomaly')

# Setup
gcms <- c('ACCESS-ESM1-5') # 'ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0'
ssps <- c('ssp245','ssp370')
vrss <- c('pr','tasmax','tasmin')
prds <- c('2021_2040','2041_2060', '2061_2080', '2081_2100')
stp <- base::expand.grid(gcms,ssps,vrss,prds) %>% base::as.data.frame()
names(stp) <- c('gcm','ssp','var','prd'); rm(gcms, ssps, vrss, prds)
stp <- stp %>%
  dplyr::arrange(gcm,ssp,prd,var) %>%
  base::as.data.frame()

# gcm <- gcms[1]
# ssp <- ssps[1]
# var <- vrss[5]
# prd <- prds[1]

# Read monthly deltas
get_daily_future_data <- function(gcm, ssp, var, prd){
  cat(paste0("processing ",var,'_',gcm,'_',ssp,'_',prd,' \n'))
  prd <- as.character(prd)
  file <- paste0('CMIP6_',gcm,'_',ssp,'_r1i1p1f1_',var,'_America_monthly_intp_anomaly_',prd,'.tif')
  # Load deltas
  dlts <- terra::rast(paste0(anm_pth,'/',file))
  dlts <- dlts %>% terra::resample(ref)
  dlts <- dlts %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
  prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
  his_yrs <- 1995:2014
  fut_yrs <- prd_num[1]:prd_num[2]
  
  # Temporal mapping
  Baseline = seq(from = as.Date('1995-01-01'),
                 to   = as.Date('2014-12-31'),
                 by   = 'day')
  Future   = seq(from = as.Date(paste0(prd_num[1],'-01-01')),
                 to   = as.Date(paste0(prd_num[2],'-12-31')),
                 by   = 'day')
  # Remove feb-29 from all leap years (not coincidence between years)
  Baseline <- Baseline[!(format(Baseline,"%m") == "02" & format(Baseline, "%d") == "29"), drop = FALSE]
  Future   <- Future[!(format(Future,"%m") == "02" & format(Future, "%d") == "29"), drop = FALSE]
  mpg <- data.frame(Baseline, Future)
  mpg$year  <- lubridate::year(mpg$Baseline)
  mpg$year_fut <- lubridate::year(mpg$Future)
  mpg$month <- lubridate::month(mpg$Baseline)
  
  if(var %in% c('pr','rsds')){
    # Paths
    if (var == "pr") {
      his_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps' 
      fut_pth <- paste0(root,'/chirps_cmip6_america/Prec_',gcm,'_',ssp,'_',prd); dir.create(fut_pth, F, T)
    } else if (var == "rsds") {
      his_pth <- paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/solar_radiation_flux')
      fut_pth <- paste0(root,'/ecmwf_agera5_cmip6_america/solar_radiation_flux_',gcm,'_',ssp,'_',prd)
      dir.create(fut_pth, F, T)
    }
    if(length(list.files(fut_pth)) < 7300){
      # File structure
      if (var == "pr") {
        his_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else if (var == "rsds") {
        his_str <- paste0('Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=mpg$Baseline, fixed=T),'_final-v1.0.nc')
        fut_str <- paste0('Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=mpg$Future, fixed=T),'_final-v1.0.nc')
      }
      # Split by months
      his_lst <- split(his_str, mpg$month)
      fut_lst <- split(fut_str, mpg$month)
      1:length(his_lst) %>%
        purrr::map(.f = function(j){
          delta <- dlts[[j]]
          his_daily <- his_lst[[j]]
          fut_daily <- fut_lst[[j]]
          plan(multicore, workers = 5)
          1:length(his_daily) %>%
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',fut_daily[k])
              if(!file.exists(outfile)){
                r <- terra::rast(paste0(his_pth,'/',his_daily[k]))
                if (var == "pr") {
                  r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
                  r[r == -9999] <- NA
                } else if (var == "rsds") {
                  r <- r %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
                  r <- r * 1e-6
                }
                r <- r * (1 + delta)
                terra::writeRaster(r, outfile)
              }
            })
          plan(sequential); gc(reset = T)
        })
    }
  }
  if(var %in% c('tasmax','tasmin','hurs')){
    # Paths
    his_pth <- ifelse(var == 'tasmax',
                      paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmax'),
                      ifelse(var == 'tasmin',
                             paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmin'), 
                             paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Rhum')))
    fut_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts_cmip6_america/Tmax_',gcm,'_',ssp,'_',prd),
                      ifelse(var == 'tasmin',
                             paste0(root,'/chirts_cmip6_america/Tmin_',gcm,'_',ssp,'_',prd),
                             paste0(root,'/chirts_cmip6_america/RHum_',gcm,'_',ssp,'_',prd)))
    if(length(list.files(fut_pth)) < 7300){
      # File structure
      if(var == 'tasmax'){
        his_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else if (var == 'tasmin'){
        his_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else {
        his_str <- paste0('RH.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('RH.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      }
      yrs_str   <- mpg$year
      yrs_f_str <- mpg$year_fut
      # Split by months
      his_lst   <- split(his_str, mpg$month)
      fut_lst   <- split(fut_str, mpg$month)
      yrs_lst   <- split(yrs_str, mpg$month)
      yrs_f_lst <- split(yrs_f_str, mpg$month)
      1:length(his_lst) %>%
        purrr::map(.f = function(j){
          delta <- dlts[[j]]
          his_daily   <- his_lst[[j]]
          fut_daily   <- fut_lst[[j]]
          yrs_daily   <- yrs_lst[[j]]
          yrs_f_daily <- yrs_f_lst[[j]]
          plan(multicore, workers = 5)
          1:length(his_daily) %>%
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',yrs_f_daily[k],'/',fut_daily[k]); dir.create(dirname(outfile),F,T)
              if(!file.exists(outfile)){
                r <- terra::rast(paste0(his_pth,'/',his_daily[k])) # yrs_daily[k],'/',
                r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
                r[r == -9999] <- NA
                r <- r + delta
                #for hurs limit min to 0, max to 100
                if (var == "hurs") {
                  r <- min(r, 100)
                  r <- max(r, 0)
                }
                terra::writeRaster(r, outfile)
              }
            })
          plan(sequential); gc(reset = T)
        })
    }
  }
  
  return(cat(paste0(var,'_',gcm,'_',ssp,'_',prd,' ready.\n')))
  
}

1:nrow(stp) %>%
  purrr::map(.f = function(i){get_daily_future_data(gcm = stp$gcm[i],
                                                    ssp = stp$ssp[i],
                                                    var = stp$var[i],
                                                    prd = stp$prd[i])
    gc(verbose=F, full=T, reset=T)
  })
