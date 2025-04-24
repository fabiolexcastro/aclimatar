

## Calc number of days with temperature above 35°C
## By: H. Achicanoy / J. Ramirez-Villegas
## December, 2022
## Edited by Fabio Castro-Llanos
## March, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

# Root folder 
root <- './common_data'

# Function
calc_ntx35 <- function(yr, mn, thr = 35){
  
  # yr <- '2021'; mn <- '01'
  
  outfile <- paste0(out_dir, 'NTX35-', yr, '-', mn, '.tif')
  cat(basename(outfile), '\n')
  
  if(!file.exists(outfile)){
    
    ## Create output directory
    dir.create(dirname(outfile), F, T)
    
    ## Tidy the dates
    cat('Processing -------> ', yr, mn, '\n')
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    fls <- paste0(tx_pth, '/', yr, '/', 'Tmax.', gsub('-', '\\.', dts), '.tif')
    fls <- fls[file.exists(fls)]
    
    ## Read the precipitation data
    tmx <- rast(fls)
    
    ## Calculate ntx35
    terra::app(x   = tmx,
               fun = function(x){ ntxval = sum(x >= thr, na.rm = T); return(ntxval)},
               filename = outfile)
    
    rm(tmx, dts, fls)
    
  }
  
}

# Future setup
ssp <- 'ssp245'
gcm <- 'EC-Earth3'
prd <- '2021_2040'

for(ssp in c('ssp245', 'ssp370')){
  
  for(prd in c('2021_2040', '2041_2060', '2061_2080', '2081_2100')){
    
    for(gcm in c('ACCESS-ESM1-5', 'MPI-ESM1-2-HR', 'EC-Earth3', 'INM-CM5-0', 'MRI-ESM2-0')){
      
      ## To start the process
      cat('----------------------', ssp, ' ', gcm, '--------------------------\n')
      
      ## Parameters
      cmb <- paste0(ssp, '_', gcm)
      mnt <- c(paste0('0', 1:9), 10:12)
      yrs <- str_sub(prd, 1, 4):str_sub(prd, 6, 9)
      stp <- base::expand.grid(yrs, mnt) %>% base::as.data.frame() %>% setNames(c('yrs', 'mnt')) %>% arrange(yrs, mnt) %>% as.data.frame(); rm(yrs, mnt)
      
      ## Setup in/out files
      tx_pth  <- paste0(root, '/chirts_cmip6_america/Tmax_', gcm, '_', ssp, '_', prd)
      out_dir <- paste0(root, '/indices/', ssp, '_', gcm, '/NTX35/') 
      
      1:nrow(stp) %>% 
        purrr::map(.f = function(i){calc_ntx35(yr = stp$yrs[i], mn = stp$mnt[i], thr = 35)})
      
      ##
      cat('----Finish----\n') 
      
    }
    
  }
  
}

