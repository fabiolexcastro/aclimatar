
# Load libraires ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Vector data 
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
zne <- vect(zne)

## Directories 
as.character(dir_ls('./cmip6_check_monthly-avrg'))
bsln <- dir_ls('./common_data/indices_average/historical') %>% map(dir_ls) %>% unlist() %>% as.character()
ftre <- as.character(dir_ls('./common_data/indices_global/ssp370'))
indx <- c('HSH', 'NDWL0', 'NDD', 'NTX30', 'NTX35', 'TAI')

# Function for making the maps --------------------------------------------
make.maps <- function(idx){
  
  idx <- 'NDD'
  
  ## Filtering
  cat('>>>> To start: ', idx, '\n')
  r.bsln <- rast(grep(idx, bsln, value = T))
  r.ftre <- rast(grep(idx, ftre, value = T))
  
  ## Make just one stack 
  stck <- c(r.bsln, r.ftre)
  
  ## Extract by mask 
  stck <- terra::crop(stck, zne)
  stck <- terra::mask(stck, zne)
  
  ## Raster to table 
  tble <- terra::as.data.frame(stck, xy = T) %>% 
    as_tibble() %>% 
    gather(var, value, -c(x, y)) %>% 
    separate(data = ., col = 'var', into = c('index', 'period'), sep = '_', remove = T) %>% 
    mutate(period = ifelse(period == 'bsl', 'Baseline', 'Future (SSP 370)'), 
           period = factor(period, levels = c('Baseline', 'Future (SSP 370)'))) 
  
  ## To select the colors 
  if(idx == 'HSH'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 9, name = 'Reds')
    ttle <- 'Human Heat Stress'
    sbtt <- 'Index'
  } else if(idx == 'NDWL0'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 9, name = 'BrBG')
    ttle <- 'Number of dry days with soil waterlogging at moisture content at start of saturarion or above'
    sbtt <- 'Days'
  } else if(idx == 'NTX30'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 9, name = 'YlOrRd')
    ttle <- 'Number of heat stress (Tmax > 30º)'
    sbtt <- 'Days'
  } else if(idx == 'TAI'){
    cat(idx, '\t')
    clrs <- rev(brewer.pal(n = 9, name = 'BrBG'))
    ttle <- 'Thornthwaite aridity index'
    sbtt <- '%'
  } else if(idx == 'NTX35'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 9, name = 'YlOrRd')
    ttle <- 'Number of heat stress (Tmax > 35º)'
    sbtt <- 'Days'
  } else if(idx == 'NDD'){
    cat(idx, '\t')
    clrs <- rev(brewer.pal(n = 9, name = 'BrBG'))
    ttle <- 'Number of dry days'
    sbtt <- 'Days'
  }
  
  ## To draw the map (raw values)
  g.vls <- ggplot() + 
    geom_tile(data = tble, aes(x = x, y = y, fill = value)) + 
    facet_wrap(.~period) + 
    scale_fill_gradientn(colors = clrs, na.value = 'white', 
                         guide = guide_colorbar(
                           title.position = 'top',
                           title.hjust = 0.5,
                           barwidth = unit(20, 'line'),
                           barheight = unit(0.5, 'line')
                         )) +
    geom_sf(data = st_as_sf(zne), fill = NA, col = 'grey30') +
    coord_sf() +
    labs(x = '', y = '', fill = idx) +
    ggtitle(label = ttle, 
            subtitle = sbtt) +
    theme_minimal() +
    theme(
      legend.position = 'bottom', # legend.key.width = unit(3, 'line'), 
      axis.text.x = element_text(size = 5), 
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 13),
      plot.subtitle = element_text(face = 'bold', hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 5, angle = 90, hjust = 0.5), 
      strip.text = element_text(face = 'bold', hjust = 0.5)
    )
  
  ### To save the map 
  ggsave(plot = g.vls, filename = glue('./png/maps_index/values_{idx}.jpg'), units = 'in', width = 10, height = 7, dpi = 300, create.dir = T)
  
  ## To make the difference
  dfrn <- tble %>% spread(period, value) %>% mutate(Difference = `Future (SSP 370)` - Baseline)
  # dfrn <- dfrn %>% filter(Difference > -350)
  
  ## To draw the map (difference values)
  g.dfr <- ggplot() + 
    geom_tile(data = dfrn, aes(x = x, y = y, fill = Difference)) + 
    scale_fill_gradientn(colors = clrs, na.value = 'white', 
                         guide = guide_colorbar(
                           title.position = 'top',
                           title.hjust = 0.5,
                           barwidth = unit(20, 'line'),
                           barheight = unit(0.5, 'line')
                         )) +
    geom_sf(data = st_as_sf(zne), fill = NA, col = 'grey30') +
    coord_sf() +
    labs(x = '', y = '', fill = idx) +
    theme_minimal() +
    theme(
      legend.position = 'bottom', # legend.key.width = unit(3, 'line'), 
      axis.text.x = element_text(size = 5), 
      axis.text.y = element_text(size = 5, angle = 90, hjust = 0.5), 
      strip.text = element_text(face = 'bold', hjust = 0.5)
    )
  
  ggsave(plot = g.dfr, filename = glue('./png/maps_index/difference_{idx}.jpg'), units = 'in', width = 7, height = 6, dpi = 300, create.dir = T)
  
  
}

# To make the maps --------------------------------------------------------
map(indx, make.maps)

