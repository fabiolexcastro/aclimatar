

# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Matrix ------------------------------------------------------------------
rcl.hsh <- tibble(value = 1:4, min = c(0, 25, 30, 35), max = c(25, 30, 35, Inf), class = c('Mild or no stress', 'Caution', 'Danger', 'Extreme danger'))

# Functions ---------------------------------------------------------------
my.rcl <- function(rst, sov, max){
  
  # rst <- n35
  # sov <- 'ECU'
  # max <- 366
  
  ## To make the filtering
  cat('To start the process!\n')
  pnt <- filter(pnts, iso %in% sov)
  vct <- zone[zone$sov_a3 == sov,]
  
  ## To extract by mask 
  rst <- terra::crop(rst, vct)
  rst <- terra::mask(rst, vct)
  
  ## Extract the values for the points 
  vls <- cbind(pnt, terra::extract(rst, pnt[,c('lon', 'lat')]))
  vls <- as_tibble(vls)
  colnames(vls)[5] <- 'bsl'
  
  ## Quantile
  qnt <- quantile(vls$bsl, c(0, 0.05, 0.75, 1), na.rm = T)
  mtr <- matrix(c(0, qnt[2], 1, qnt[2], qnt[3], 2, qnt[3], max, 3), byrow = T, ncol = 3)
  
  ## To classify 
  rcl <- terra::classify(rst, mtr, include.lowest = T)
  
  ## Finish 
  cat('Done!\n')
  return(rcl)
  
}

# Load data ---------------------------------------------------------------

## List the index (baseline versus future)
bsln <- as.character(unlist(map(dir_ls('./common_data/indices_average/historical'), dir_ls)))
ftre <- as.character(dir_ls('./common_data/indices_global/ssp370'))

## Vector data
wrld <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zone <- filter(wrld, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
zone <- vect(zone)

## Points datasets coffee --------------------------------------------------

### Southamerica
pnts.col <- read_csv('//alliancedfs.alliance.cgiar.org/CL9_Coffee_Cocoa2/_coffeeColombia/tbl/points/coffee_v1.csv')

### Centralamerica
pnts.cam <- read_csv('./tble/points/Arabica_cleaned_up_points.csv')

### Spatial intersection 
pnts.cam$sov_a3 <- terra::extract(vect(wrld), pnts.cam[,c('Longitude', 'Latitude')])$sov_a3
pnts.cam <- filter(pnts.cam, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
pnts.cam <- pnts.cam[,-1]
length(unique(pnts.cam$sov_a3))

pnts.cff <- pnts.cam

### To write the rasters
write.csv(pnts.cff, './tble/points/points_coffee.csv', row.names = FALSE)
pnts <- read_csv('./tble/points/points_coffee.csv', show_col_types = FALSE)
colnames(pnts) <- c('lon', 'lat', 'iso')

# To classify HSH ---------------------------------------------------------
hsh <- c(bsln, ftre) %>% grep('HSH', ., value = T)
hsh <- rast(hsh)
mtx.hsh <- as.matrix(dplyr::select(rcl.hsh, min, max, value))
hsh.rcl <- terra::classify(hsh, mtx.hsh)


# To classify NTx30 -------------------------------------------------------
n30 <- c(bsln, ftre) %>% grep('NTX30', ., value = T)
n30 <- rast(n30)

n30.rcl <- list()

for(i in 1:nrow(zone)){
  
  sve <- zone$sov_a3[i]
  n30.rcl[[i]] <- my.rcl(rst = n30, sov = sve, max = 366)
  
}

n30.rcl <- sprc(n30.rcl)
n30.rcl <- mosaic(n30.rcl, fun = 'modal')

# To classify NTx35 -------------------------------------------------------
n35 <- c(bsln, ftre) %>% grep('NTX35', ., value = T)
n35 <- rast(n35)

n35.rcl <- list()

for(i in 1:nrow(zone)){
  
  sve <- zone$sov_a3[i]
  n35.rcl[[i]] <- my.rcl(rst = n35, sov = sve, max = 366)
  
}

n35.rcl <- sprc(n35.rcl)
n35.rcl <- mosaic(n35.rcl, fun = 'modal')
n35.rcl <- terra::mask(n35.rcl, zone)
plot(n35.rcl)

# To classify NDD ---------------------------------------------------------
ndd <- c(bsln, ftre) %>% grep('NDD', ., value = T)
ndd <- rast(ndd)

ndd.rcl <- list()

for(i in 1:nrow(zone)){
  
  sve <- zone$sov_a3[i]
  ndd.rcl[[i]] <- my.rcl(rst = ndd, sov = sve, max = 366)
  
}

ndd.rcl <- sprc(ndd.rcl)
ndd.rcl <- mosaic(ndd.rcl, fun = 'modal')
plot(ndd.rcl)

# To classify NDWL0 -------------------------------------------------------
ndw <- c(bsln, ftre) %>% grep('NDWL', ., value = T)
ndw <- rast(ndw)

ndw.rcl <- list()

for(i in 1:nrow(zone)){
  
  sve <- zone$sov_a3[i]
  ndw.rcl[[i]] <- my.rcl(rst = ndw, sov = sve, max = 365)
  
}

ndw.rcl <- sprc(ndw.rcl)
ndw.rcl <- mosaic(ndw.rcl, fun = 'modal')

# To classify TAI ---------------------------------------------------------
tai <- c(bsln, ftre) %>% grep('TAI', ., value = T)
tai <- rast(tai)

tai.rcl <- list()

for(i in 1:nrow(zone)){
  
  sve <- zone$sov_a3[i]
  tai.rcl[[i]] <- my.rcl(rst = tai, sov = sve, max = 101)
  
}

tai.rcl <- sprc(tai.rcl)
tai.rcl <- mosaic(tai.rcl, fun = 'modal')

# To write the raster -----------------------------------------------------

stck <- list(hsh.rcl, n30.rcl, n35.rcl, ndd.rcl, ndw.rcl, tai.rcl)
stck <- reduce(stck, c)
terra::writeRaster(x = stck, filename = './common_data/indices_global/ssp370/index_rcl-coffee.tif', overwrite = TRUE)

