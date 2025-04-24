
# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, rnaturalearth, gtools, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

fles <- '//catalogue/WFP_ClimateRiskPr1/1.Data/CMIP6/America_all' %>% 
  dir_ls(., regexp = '.tif$') %>% 
  as.character()

wrld <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
amrc <- filter(wrld, region_un == 'Americas')

# A simple plot -----------------------------------------------------------
make.plot <- function(ssp){
  
  ##
  cat('To process: ', ssp, '\n')
  fls <- grep(ssp, fles, value = T)

  ##
  tmn <- grep('tasmin', fls, value = T)
  tmx <- grep('tasmax', fls, value = T)
  ppt <- grep('pr', fls, value = T)
  
  ## 
  ssps <- fls %>% basename() %>% str_split(., '_') %>% map_chr(2) %>% unique()
  
  ## 
  for(i in 1:length(tmn)){
    
    cat(i, '\n')
    sp <- ssps[i]
    tn <- rast(tmn[i])
    tx <- rast(tmx[i])
    pt <- rast(ppt[i])
    
    png(filename = glue('./png/maps_raw/{sp}.png'), width = 12, height = 4, units = 'in', res = 300)
    
    layout(matrix(1:3, ncol = 3))
    par(mar = c(3, 3, 3, 3), oma = c(0, 0, 2, 0))
    
    plot(tn[[1]], main = 'Minimum Temperature')
    plot(sf::st_geometry(amrc), add = TRUE)
    
    plot(tx[[1]], main = 'Maximum Temperature')
    plot(sf::st_geometry(amrc), add = TRUE)
    
    plot(pt[[1]], main = 'Precipitation')
    plot(sf::st_geometry(amrc), add = TRUE)
    
    mtext(text = glue('Climate Variables - Scenario {sp}'), 
          side = 3, outer = TRUE, line = 0.5, cex = 1.5)
    
    dev.off()
    
  }
  
  ## 
  cat('Done!\n')
  
}
make.plot(ssp = 'ssp370')

