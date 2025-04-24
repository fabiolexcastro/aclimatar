## Number of waterlogging days (days) at start of saturation (NDWL0)
## By: H. Achicanoy & A. Mendez
## April, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,sf,terra,gtools,lubridate,compiler, raster, ncdf4))
#source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

peest <- function(srad, tmin, tmean, tmax){
  
  #gc()
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  
  # Soil heat flux parameters
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Net radiation
  #rn <- (1-albedo) * srad
  rn <- lapply(srad, function(x){(1-albedo) * x})
  # Soil heat flux
  temp_term <- tmean + c_eslope
  exp_term <- exp(b_eslope * tmean / temp_term)
  eslope <- (a_eslope * b_eslope * c_eslope * exp_term) / (temp_term^2)
  
  
  #eslope <- a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
  
  # Estimate vpd
  esat_min <- 0.61120*exp((17.67*tmin)/(tmin+243.5))
  esat_max <- 0.61120*exp((17.67*tmax)/(tmax+243.5))
  vpd <- vpd_cte*(esat_max-esat_min) #kPa
  
  # Priestley-Taylor
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  pt_coef <- pt_fact*pt_const
  pt_coef <- 1 + (pt_coef-1) * vpd / vpd_ref
  
  #*10^6? To convert fluxes MJ to J
  #rlat_ht? Latent heat flux to water flux
  #100/rho_w? Kg/m^2 to cm
  
  es_psy <- (eslope+psycho)
  div_es_psy <- eslope/es_psy
  et_max <- (pt_coef * rn * div_es_psy * 10E6 / rlat_ht * 100/rho_w)*10 #in mm
  return(et_max)
}

peest2 <- function(srad, tmin, tmean, tmax){
  # Convert rasters to matrices for faster processing
  srad_m <- raster::values(srad)
  tmin_m <- terra::values(tmin)
  tmean_m <- terra::values(tmean)
  tmax_m <- terra::values(tmax)
  gc()
  
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Pre-allocate matrices
  n_cells <- ncol(srad_m)
  n_layers <- nrow(srad_m)
  rn <- matrix(0, nrow = n_layers, ncol = n_cells)
  eslope <- matrix(0, nrow = n_layers, ncol = n_cells)
  et_max <- matrix(0, nrow = n_layers, ncol = n_cells)
  
  # Optimized calculations
  rn <- (1-albedo) * srad_m
  rm(srad_m); gc()
  
  temp_denom <- tmean_m + c_eslope
  eslope <- (a_eslope * b_eslope * c_eslope * exp(b_eslope * tmean_m/temp_denom)) / (temp_denom^2)
  rm(tmean_m, temp_denom); gc()
  
  vpd <- vpd_cte * (0.61120 * (exp((17.67*tmax_m)/(tmax_m+243.5)) - 
                                 exp((17.67*tmin_m)/(tmin_m+243.5))))
  rm(tmin_m, tmax_m); gc()
  
  # Constants
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  # Final calculation
  pt_coef <- 1 + (pt_fact*pt_const-1) * vpd / vpd_ref
  conversion_factor <- 1E6 * 100 / (rlat_ht * rho_w) * 10
  et_max <- pt_coef * rn * eslope/(eslope+psycho) * conversion_factor
  
  # Convert back to raster
  et_max_rast <- tmin
  terra::values(et_max_rast) <- et_max
  
  return(et_max_rast)
  
}

root <- './common_data'
ref <- terra::rast('./common_data/chirps_hist_america/monthly/Prec_1995.tif')[[1]]
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries', returnclas = 'sf') 
zne <- dplyr::filter(wrl, sov_a3 %in% c('COL', 'ECU', 'BLZ', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI', 'DOM'))
zne <- vect(zne)

# Soil variables
scp <- rast('//catalogue/workspace-cluster9/2025/soils_world/soils-world_ssp-ssat/tif/ssat_world.tif')
scp <- scp |> terra::resample(ref) |> terra::mask(zne)
sst <- rast('//catalogue/workspace-cluster9/2025/soils_world/soils-world_ssp-ssat/tif/sscp_world.tif')
sst <- sst |> terra::resample(ref) |> terra::mask(zne)

cpeest <- compiler::cmpfun(peest)

# Calculate NDWS function
yr <- '2021'
mn <- '01'

calc_ndwl0 <- function(yr, mn){
  
  outfile <- paste0(out_dir,'/NDWL0-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  if (!file.exists(outfile)) {
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    # Files
    # t1 <- system.time({
      pr_fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
      pr_fls <- pr_fls[file.exists(pr_fls)]
      tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
      tx_fls <- tx_fls[file.exists(tx_fls)]
      tm_fls <- paste0(tm_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
      tm_fls <- tm_fls[file.exists(tm_fls)]
      if(as.numeric(yr) > 2020){
        dts_chr <- as.character(dts)
        dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
        dts <- as.Date(dts_chr); rm(dts_chr)
      }
      sr_fls <- paste0(sr_pth,'/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
      sr_fls <- sr_fls[file.exists(sr_fls)]
      
      # Read variables
      prc <<- terra::rast(pr_fls)
      #prc <- prc |> terra::crop(terra::ext(ref)) |> terra::mask(ref)
      #prc <- terra::classify(prc, rcl = cbind(-9999, NA)) # prc[prc == -9999] <- NA
      tmx <<- terra::rast(tx_fls)
      #tmx <- tmx |> terra::crop(terra::ext(ref)) |> terra::mask(ref)
      #tmx <- terra::classify(tmx, rcl = cbind(-9999, NA)) # tmx[tmx == -9999] <- NA
      tmn <<- terra::rast(tm_fls)
      #tmn <- tmn |> terra::crop(terra::ext(ref)) |> terra::mask(ref)
      #tmn <- terra::classify(tmn, rcl = cbind(-9999, NA)) # tmn[tmn == -9999] <- NA
      tav <<- (tmx + tmn)/2
      #avg <- function(x,y){(x+y)/2}; cavg <- compiler::cmpfun(avg)
      #tav <- terra::lapp(x = terra::sds(tmn, tmx), fun = cavg)
      srd <- terra::rast(sr_fls)
      srd <- srd/1000000
      srd <- srd |> terra::resample(ref) |> terra::mask(ref)
      srd <- raster::stack(srd)
    # })
    
    #cat("carga de rasters: ", t1[3], " seg \n")
    #cat("Numero de capas:", raster::nlayers(srd), "\n")
    # Maximum evapotranspiration
    t2 <- system.time({
      #sptd <<- terra::sds(srd,tmn,tav,tmx)
    })
    #cat("sds: ", t2[3], " seg \n")
    t2.5 <- system.time({
      #ETMAX <<- terra::lapp(x = sptd, fun = cpeest)
      ETMAX <<- peest2(srad = srd, tmin = tmn, tmean = tav, tmax = tmx)
      gc()
    })
    
    #cat("ETMAX: ", t2.5[3], " seg \n")
    rm(list = c('tmn','tmx','tav','srd','sptd'))
    gc(verbose = F, full = T, reset = T)
    
    # Compute water balance model
    t3 <- system.time({
      # Compute water balance model
      date <- paste0(yr,'-',mn)
      if(date %in% c('1995-01','2021-01','2041-01','2061-01','2081-01')){
        AVAIL <<- ref
        AVAIL[!is.na(AVAIL)] <- 0
      } else {
        AVAIL <<- terra::rast(paste0(dirname(outfile),'/AVAIL.tif'))
      }
    })
    #cat("AVAIL: ", t3[3], " seg \n")
    
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL, rain = prc[[1]], evap = ETMAX[[1]]){
      gc()
      
      soilcp_m  <- terra::values(soilcp)
      avail_m   <- terra::values(avail)
      soilsat_m <- terra::values(soilsat)
      avail_m   <- terra::values(avail)
      rain_m    <- terra::values(rain)
      evap_m    <- terra::values(evap)
      
      avail_m   <- pmin(avail_m, soilcp_m)
      
      # ERATIO
      percwt <- pmin(avail_m/soilcp_m*100, 100)
      percwt <- pmax(percwt, 1)
      eratio <- pmin(percwt/(97-3.868*sqrt(soilcp_m)), 1)
      
      demand  <- eratio * evap_m
      result  <- avail_m + rain_m - demand
      rm(rain_m, evap_m); gc()
      logging <- result - soilcp_m
      logging <- pmax(logging, 0)
      logging <- pmin(logging, soilsat_m)
      # runoff  <- result - logging + soilcp
      avail_m   <- pmin(soilcp_m, result)
      avail_m   <- pmax(avail_m, 0)
      
      #volver a raster
      
      avail_r <- rain
      eratio_r <- rain
      logging_r <- rain
      terra::values(avail_r) <- avail_m
      terra::values(logging_r) <- logging
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail_r),
                      # Demand       = demand,
                      # Eratio       = eratio_r
                      Logging      = logging_r
      )
      return(out)
    }
    
    ceabyep_calc <- compiler::cmpfun(eabyep_calc)
    
    t4 <- system.time({
      watbal <- 1:terra::nlyr(ETMAX) |>
        purrr::map(.f = function(i){
          water_balance <- eabyep_calc(soilcp  = scp,
                                       soilsat = sst,
                                       avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                       rain    = prc[[i]],
                                       evap    = ETMAX[[i]])
          AVAIL <<- water_balance$Availability
          return(water_balance)
        })
    })
    
    #cat("water_balance - eabyep_calc: ", t4[3], " seg \n")
    
    t5 <- system.time({
      LOGGING <- watbal |> purrr::map('Logging') |> terra::rast()
      NDWL0  <- sum(LOGGING > 0)
      terra::writeRaster(NDWL0, outfile, overwrite = T)
      terra::writeRaster(AVAIL[[terra::nlyr(AVAIL)]], paste0(dirname(outfile),'/AVAIL.tif'), overwrite = T)
    })
    #cat("ERATIO y save Rasters: ", t5[3], " seg \n")
    
    #cat("TIEMPO TOTAL: ", sum(c(t1[3], t2[3], t2.5[3], t3[3], t4[3], t5[3] )), "\n \n")
    #clean up
    #rm(list=c('prc','ETMAX','AVAIL','watbal','ERATIO','NDWS'))
    #gc(verbose = F, full = T, reset = T)
  }
}
#rm(list = c('tmn','tmx','tav','srd','sptd'))
#rm(list = c('prc','ETMAX','AVAIL','watbal','ERATIO','NDWS'))
#gc(verbose = F, full = T, reset = T)

# Future setup
gcm <- 'ACCESS-ESM1-5'
ssp <- 'ssp370'
prd <- c('2021_2040')

for (gcm in c('ACCESS-ESM1-5')) { # 'ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0'
  for (ssp in c('ssp370')) {
    for (prd in c('2021_2040','2041_2060')) {
      
      cmb <- paste0(ssp,'_',gcm,'_',prd)
      prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
      yrs <- prd_num[1]:prd_num[2]
      mns <- c(paste0('0',1:9),10:12)
      stp <- base::expand.grid(yrs, mns) |> base::as.data.frame(); rm(yrs,mns)
      names(stp) <- c('yrs','mns')
      stp <- stp |>
        dplyr::arrange(yrs, mns) |>
        base::as.data.frame()
      pr_pth <- paste0(root,'/chirps_cmip6_america/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
      tx_pth <- paste0(root,'/chirts_cmip6_america/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
      tm_pth <- paste0(root,'/chirts_cmip6_america/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
      sr_pth <- paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/solar_radiation_flux') # Solar radiation
      out_dir <- paste0(root,'/indices/',ssp, '_', gcm, '/NDWL0')
      
      yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 4)),
                            Future = as.character(c(2021:2040,2041:2060,2061:2080,2081:2100)))
      
      1:nrow(stp) |>
        purrr::map(.f = function(i){
          calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)
          # if (i%%5 == 0) {
          #   tmpfls <- list.files(tempdir(), full.names=TRUE)
          #   1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
          # }
        })
      
    }
  }
}
