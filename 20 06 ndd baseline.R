

## Calc number of dry days
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
calc_ndd <- function(yr, mn){
  
  # yr <- '1995'; mn <- '06'
  
  outfile <- paste0(out_dir, '/NDD-', yr, '-', mn, '.tif')
  cat(basename(outfile), '\n')
  
  system.time(
    expr = {
      if(!file.exists(outfile)){
        
        ## Create output directory
        dir.create(dirname(outfile), F, T)
        
        ## Tidy the dates
        cat('Processing -------> ', yr, mn, '\n')
        # last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
        dts <- paste0(yr, '-', as.numeric(mn))
        fls <- paste0(pr_pth, 'daily/', 'Prec_', dts, '.tif')
        fls <- fls[file.exists(fls)]
        
        ## Read the precipitation data
        prc <- rast(fls)
        
        ## Calculate ndd
        terra::app(x = prc, 
                   fun = function(x){ndd = sum(x < 1, na.rm = T); return(ndd)},
                   filename = './common_data/indices/historical/NTX30/r.tif', overwrite = T)
        
        
        rm(prc, dts, fls)
        
      }
    }
  )

  
}

## Historical setup
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>% dplyr::arrange(yrs, mns) %>% base::as.data.frame()

pr_pth <- paste0(root,'/chirps_hist_america/' ) # Daily precipitation
out_dir <- paste0(root,'/indices/historical/NDD')

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ndd(yr = stp$yrs[i], mn = stp$mns[i])})
