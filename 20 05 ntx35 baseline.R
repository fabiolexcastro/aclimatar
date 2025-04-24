

## Calc number of days with temperature above 30°C
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
thr <- 35
calc_ntx35 <- function(yr, mn, thr = 35){
  
  # yr <- '1995'; mn <- '01'
  
  outfile <- paste0(out_dir, '/NTX35-', yr, '-', mn, '.tif')
  cat(basename(outfile), '\n')
  
  if(!file.exists(outfile)){
    
    ## Create output directory
    dir.create(dirname(outfile), F, T)
    
    ## Tidy the dates
    cat('Processing -------> ', yr, mn, '\n')
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    fls <- paste0(tx_pth, '/', 'Tmax_', yr, '-', as.numeric(mn), '.tif')
    fls <- fls[file.exists(fls)]
    
    ## Read the precipitation data
    tmx <- rast(fls)
    
    ## Calculate ntx30
    terra::app(x   = tmx,
               fun = function(x){ ntxval = sum(x >= thr, na.rm = T); return(ntxval)},
               filename = outfile)
    
    rm(tmx, dts, fls)
    
  }
  
}

# Historical setup
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
tx_pth <- paste0(root,'/chirts_hist_america/daily/Tmax') # Daily maximum temperature
out_dir <- paste0(root,'/indices/historical/NTX35')
1:nrow(stp) %>%
  purrr::map(.f = function(i){
    calc_ntx35(yr = stp$yrs[i], mn = stp$mns[i], thr=35)
    tmpfls <- list.files(tempdir(), full.names=TRUE)
    1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
  })
