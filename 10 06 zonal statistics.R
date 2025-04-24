
# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
root <- './common_data/indices'
dirs <- as.character(dir_ls(root))
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
ssps <- c('ssp245', 'ssp370')

wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'BLZ', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI', 'DOM'))
zne <- vect(zne)

# Baseline ----------------------------------------------------------------
prec.bsln <- as.character(dir_ls('./common_data/chirps_hist_america/monthly', regexp = '.tif$')) %>% rast()
tmin.bsln <- as.character(dir_ls('./common_data/chirts_hist_america/monthly/Tmin')) %>% rast()
tmax.bsln <- as.character(dir_ls('./common_data/chirts_hist_america/monthly/Tmax')) %>% rast()

# Future ------------------------------------------------------------------
gcms <- c('ACCESS-ESM1-5', 'MPI-ESM1-2-HR', 'EC-Earth3', 'INM-CM5-0', 'MRI-ESM2-0')
prec.ftre <- as.character(unlist(map(grep('ssp370', as.character(dir_ls('./common_data/chirps_cmip6_america', regexp = 'Prec')), value = T), dir_ls)))
tmin.ftre <- as.character(dir_ls(as.character(unlist(map(grep('ssp370', as.character(dir_ls('./common_data/chirts_cmip6_america', regexp = 'Tmin')), value = T), dir_ls)))))
tmax.ftre <- as.character(dir_ls(as.character(unlist(map(grep('ssp370', as.character(dir_ls('./common_data/chirts_cmip6_america', regexp = 'Tmax')), value = T), dir_ls)))))

# Function baseline -------------------------------------------------------

## Summarise
calc.avrg.bsln <- function(rst){
  
  cat('To start!\n')
  avg <- map(.x = 1:12, .f = function(i){
    trr <- rst[[grep(paste0('-', i, '$'), names(rst))]]
    trr <- mean(trr)
    return(trr)
  })
  avg <- reduce(avg, c)
  names(avg) <- glue('rstr_{1:12}')
  return(avg)
  
}
prec.bsln.avrg <- calc.avrg.bsln(prec.bsln) 
tmin.bsln.avrg <- calc.avrg.bsln(tmin.bsln) 
tmax.bsln.avrg <- calc.avrg.bsln(tmax.bsln) 

# To write the rasters
names(prec.bsln.avrg) <- glue('prec_{1:12}')
names(tmin.bsln.avrg) <- glue('tmin_{1:12}')
names(tmax.bsln.avrg) <- glue('tmax_{1:12}')

bsln <- c(prec.bsln.avrg, tmin.bsln.avrg, tmax.bsln.avrg)
terra::writeRaster(x = bsln, filename = glue('./cmip6_check_monthly-avrg/rst_hist.tif'), overwrite = TRUE)

# Function future ---------------------------------------------------------

## Summarise
calc.avrg.ftre <- function(gcm){
  
  gcm <- gcms[4]
  
  ## To list the files
  fls.ppt <- grep(gcm, prec.ftre, value = T)
  fls.tmn <- grep(gcm, tmin.ftre, value = T)
  fls.tmx <- grep(gcm, tmax.ftre, value = T)
  
  ## Periods
  prd <- c('2021_2040', '2041_2060')
  
  ## Loop by periods 
  map(.x = prd, .f = function(p){
    
    outfile <- glue('./cmip6_check_monthly-avrg/rst_ssp370_{p}_{gcm}.tif')
    
    if(!file.exists(outfile)){
      
      cat(p, '\n')
      fl.ppt <- grep(p, fls.ppt, value = T)
      fl.tmn <- grep(p, fls.tmn, value = T)
      fl.tmx <- grep(p, fls.tmx, value = T)
      
      years <- str_sub(p, 1, 4):str_sub(p, 6, nchar(p))
      
      r.year <- map(.x = years, .f = function(y){
        
        r.mnt <- map(.x = 1:12, .f = function(m){
          
          cat('>>>', y, '-', m, '\n')
          mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
          f.ppt <- grep(paste0(y, '.', mnt), fls.ppt, value = T)
          f.tmn <- grep(paste0(y, '.', mnt), fls.tmn, value = T)
          f.tmx <- grep(paste0(y, '.', mnt), fls.tmx, value = T)
          
          r.ppt <- rast(f.ppt)
          r.ppt <- sum(r.ppt)
          
          r.tmn <- rast(f.tmn)
          r.tmn <- mean(r.tmn)
          
          r.tmx <- rast(f.tmx)
          r.tmx <- mean(r.tmx)
          
          names(r.ppt) <- glue('prec_{y}-{mnt}')
          names(r.tmn) <- glue('tmin_{y}-{mnt}')
          names(r.tmx) <- glue('tmax_{y}-{mnt}')
          
          return(c(r.ppt, r.tmn, r.tmx))
          
          
        })
        
        r.mnt <- reduce(r.mnt, c)
        return(r.mnt)
        
      })
      
      gc(reset = T)
      
      r.prec <- map(.x = 1:length(r.year), .f = function(x){
        cat(x, '\n')
        r <- r.year[[x]][[grep('prec', names(r.year[[x]]))]]
        return(r)
      })
      r.prec <- reduce(r.prec, c)
      
      r.tmin <- map(.x = 1:length(r.year), .f = function(x){
        cat(x, '\n')
        r <- r.year[[x]][[grep('tmin', names(r.year[[x]]))]]
        return(r)
      })
      r.tmin <- reduce(r.tmin, c)
      
      r.tmax <- map(.x = 1:length(r.year), .f = function(x){
        cat(x, '\n')
        r <- r.year[[x]][[grep('tmax', names(r.year[[x]]))]]
        return(r)
      })
      r.tmax <- reduce(r.tmax, c)
      
      r.avrg <- map(.x = 1:12, .f = function(m){
        
        cat(m, '\n')
        mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
        
        ppt <- r.prec[[grep(paste0(paste0('prec_', years, '-', mnt), collapse = '|'), names(r.prec))]]
        ppt <- mean(ppt)
        
        tmn <- r.tmin[[grep(paste0(paste0('tmin_', years, '-', mnt), collapse = '|'), names(r.tmin))]]
        tmn <- mean(tmn)
        
        tmx <- r.tmax[[grep(paste0(paste0('tmax_', years, '-', mnt), collapse = '|'), names(r.tmax))]]
        tmx <- mean(tmx)
        
        names(ppt) <- glue('prec_{m}')
        names(tmn) <- glue('tmin_{m}')
        names(tmx) <- glue('tmax_{m}')
        
        return(c(ppt, tmn, tmx))
        
      })
      
      r.avrg <- reduce(r.avrg, c)
      terra::writeRaster(x = r.avrg, filename = outfile, overwrite = TRUE)
      rm(r.avrg, r.year, r.prec, r.tmin, r.tmax)
      gc(reset = T)
      
    }
    
  })
  
  
  
  
}


# To calculate the difference ---------------------------------------------
fles <- as.character(dir_ls('./cmip6_check_monthly-avrg', regexp = 'ssp370'))
gcms <- basename(fles) %>% str_split(., '_') %>% map_chr(5) %>% unlist() %>% unique() %>% gsub('.tif$', '', .)

make.dfrn <- function(gcme){
  
  # gcme <- 'ACCESS-ESM1-5'
  
  cat('To process: ', gcme, '\n')
  fls <- grep(gcme, fles, value = T)
  ft1 <- as.character(grep('2021_2040', fls, value = T))
  ft2 <- as.character(grep('2041_2060', fls, value = T))
  
  ft1 <- rast(ft1)
  ft2 <- rast(ft2)
  
  ftr <- list(ft1, ft2)
  
  a.bsl.ppt <- sum(prec.bsln.avrg)
  a.bsl.tmn <- mean(tmin.bsln.avrg)
  a.bsl.tmx <- mean(tmax.bsln.avrg)
  
  prd <- c('2021_2040', '2041_2060')
  
  dfr <- map(.x = 1:length(ftr), .f = function(i){
    
    r.ftr <- ftr[[i]]
    
    r.ftr.ppt <- r.ftr[[grep('prec', names(r.ftr))]]
    r.ftr.tmn <- r.ftr[[grep('tmin', names(r.ftr))]]
    r.ftr.tmx <- r.ftr[[grep('tmax', names(r.ftr))]]
    
    a.ftr.ppt <- sum(r.ftr.ppt)
    a.ftr.tmn <- mean(r.ftr.tmn)
    a.ftr.tmx <- mean(r.ftr.tmx)
    
    dfr.ppt <- ((a.ftr.ppt - a.bsl.ppt) / a.bsl.ppt) * 100
    dfr.tmn <- a.ftr.tmn - a.bsl.tmn
    dfr.tmx <- a.ftr.tmx - a.bsl.tmx
    
    names(dfr.ppt) <- glue('prec_{prd[i]}')
    names(dfr.tmn) <- glue('tmin_{prd[i]}')
    names(dfr.tmx) <- glue('tmax_{prd[i]}')
    
    return(c(dfr.ppt, dfr.tmn, dfr.tmx))
    
  })
  
  dfr <- reduce(dfr, c)
  names(dfr) <- glue('{names(dfr)}_{gcme}')
  return(dfr)
  
}

dfrn <- map(gcms, make.dfrn)
dfrn <- reduce(dfrn, c)

# To draw the maps --------------------------------------------------------

## 
r.prec <- dfrn[[grep('prec', names(dfrn))]]
r.tmin <- dfrn[[grep('tmin', names(dfrn))]]
r.tmax <- dfrn[[grep('tmax', names(dfrn))]]

##
make.map <- function(rst, prd){
  
  # rst <- r.prec
  # prd <- '2021_2040'
  
  cat('To process: ', prd, '\n')
  rst <- rst[[grep(prd, names(rst))]]
  lbl <- names(rst) %>% str_sub(., 1, 4) %>% unique()
  
  ## Color 
  if(lbl == 'prec'){
    library(RColorBrewer)
    clr <- brewer.pal(n = 9, name = 'BrBG')
  } else {
    library(RColorBrewer)
    clr <- brewer.pal(n = 9, name = 'YlOrRd')
  }
  
  ## Raster to table 
  tbl <- rst %>% 
    terra::as.data.frame(xy = T) %>% 
    as_tibble() %>% 
    gather(var, value, -c(x, y)) %>% 
    mutate(gcme = str_sub(var, 16, nchar(var)))
  
  ## To draw the map 
  g.dfrn <- ggplot() + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = value), alpha = 1) + 
    facet_wrap(~gcme) + 
    scale_fill_gradientn(colors = clr) +
    geom_sf(data = st_as_sf(zne), fill = NA, col = 'grey30') +
    labs(fill = lbl, x = '', y = '') +
    ggtitle(label = glue('{lbl} - {prd} vs Baseline')) +
    coord_sf() +
    theme_minimal() +
    theme(
      legend.position = 'bottom', 
      plot.title = element_text(hjust = 0.5, face = 'bold'), 
      strip.text = element_text(face = 'bold'),
      legend.key.width = unit(3, 'line'),
      axis.text.x = element_text(size = 5), 
      axis.text.y = element_text(size = 5, angle = 90, hjust = 0.5)
    )
  
  ggsave(plot = g.dfrn, filename = glue('./png/maps_dfrn/{lbl}_{prd}.jpg'), units = 'in', width = 9, height = 7, dpi = 300, create.dir = T)
  cat('Done!\n')
  rm(tbl, rst)
  gc(reset = T)
  
}

## 
prds <- c('2021_2040', '2041_2060')

map(.x = 1:length(prds), .f = function(i){
  make.map(rst = r.prec, prd = prds[i])
  make.map(rst = r.tmin, prd = prds[i])
  make.map(rst = r.tmax, prd = prds[i])
})


