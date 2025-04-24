
## Total precipitation per month
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
calc_ptot <- function(yr, mn){
  
  # yr <- '2021'; mn <- '01'
  
  outfile <- paste0(out_dir, 'PTOT-', yr, '-', mn, '.tif')
  cat(basename(outfile), '\n')
  
  if(!file.exists(outfile)){
    
    ## Create output directory
    dir.create(dirname(outfile), F, T)
    
    ## Tidy the dates
    cat('Processing -------> ', yr, mn, '\n')
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    fls <- paste0(pr_pth, '/', 'chirps-v2.0.', gsub('-', '\\.', dts), '.tif')
    fls <- fls[file.exists(fls)]
    
    ## Read the precipitation data
    prc <- rast(fls)

    ## To calculate the number of dry days 
    prc <- sum(prc)
    
    ## To write the raster 
    terra::writeRaster(x = prc, filename = outfile, overwrite = T)
    rm(prc, dts, fls)
    
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
      pr_pth  <- paste0(root, '/chirps_cmip6_america/Prec_', gcm, '_', ssp, '_', prd)
      out_dir <- paste0(root, '/indices/', ssp, '_', gcm, '/PTOT/') 
      
      1:nrow(stp) %>% 
        purrr::map(.f = function(i){calc_ptot(yr = stp$yrs[i], mn = stp$mnt[i])})
      
      ##
      cat('----Finish----\n')
      
    }
    
  }

}

lsta <- './common_data/indices/ssp245_EC-Earth3/PTOT' %>% 
  dir_ls()
