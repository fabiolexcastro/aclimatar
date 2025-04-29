
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rnaturalearthdata, raptr, rnaturalearth, gtools, glue, RColorBrewer) 

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Raster data
fles <- dir_ls('./common_data/climate_indices', regexp = '.tif$')
ccoa <- rast('./common_data/climate_indices/clima_index_cocoa.tif')
cffe <- rast('./common_data/climate_indices/clima_index_coffee.tif')

nlyr(ccoa)
nlyr(cffe)

## Vector data 
wrl <- rnaturalearth::ne_countries(scale = 50, type = 'countries') 
zne <- filter(wrl, sov_a3 %in% c('COL', 'ECU', 'PER', 'MEX', 'NIC', 'SLV', 'GTM', 'HND', 'PAN', 'CRI'))
zne <- vect(zne)

# Just climate ------------------------------------------------------------
clma <- ccoa[[1:252]]

# To make the maps for the climate ----------------------------------------

## Function
make.maps.vars <- function(vrb){
  
  # vrb <- 'prec'
  
  ## To filter just the target variable
  cat('To process: ', vrb, '\n')
  stk <- clma[[grep(vrb, names(clma))]]
  
  if(vrb == 'prec'){
    stk <- stk[[-grep('min', names(stk), value = F)]]
    stk <- stk[[-grep('max', names(stk), value = F)]]
    names(stk)[25:36] <- gsub('avg-ftr', 'ftr', names(stk)[25:36])
  }
  
  ## Raster to table 
  tbl <- stk %>% 
    terra::as.data.frame(xy = T) %>% 
    as_tibble() %>% 
    gather(var, value, -c(x, y)) %>% 
    separate(col = var, into = c('variable', 'month', 'period'), sep = '_')
  
  tbl <- tbl %>% 
    mutate(month = as.numeric(month)) %>% 
    inner_join(., tibble(month = 1:12, month_abb = month.abb), by = 'month') %>% 
    mutate(month_abb = factor(month_abb, levels = month.abb), 
           period = ifelse(period == 'bsl', 'Baseline', ifelse(period == 'frc', 'Forecast', 'Future')))
  
  tbl <- tbl %>% 
    mutate(type = paste0(period, ' - ', month_abb))
  
  ## To order 
  tbl <- tbl %>% 
    mutate(type = factor(type, levels = c(
      'Baseline - Jan', 'Forecast - Jan', 'Future - Jan',
      'Baseline - Feb', 'Forecast - Feb', 'Future - Feb',
      'Baseline - Mar', 'Forecast - Mar', 'Future - Mar',
      'Baseline - Apr', 'Forecast - Apr', 'Future - Apr',
      'Baseline - May', 'Forecast - May', 'Future - May',
      'Baseline - Jun', 'Forecast - Jun', 'Future - Jun',
      'Baseline - Jul', 'Forecast - Jul', 'Future - Jul',
      'Baseline - Aug', 'Forecast - Aug', 'Future - Aug',
      'Baseline - Sep', 'Forecast - Sep', 'Future - Sep',
      'Baseline - Oct', 'Forecast - Oct', 'Future - Oct',
      'Baseline - Nov', 'Forecast - Nov', 'Future - Nov',
      'Baseline - Dec', 'Forecast - Dec', 'Future - Dec'
    )))
  
  # tbl %>% 
  #   filter(month_abb == 'Jan') %>% 
  #   dplyr::select(-period, -month_abb, -month) %>% 
  #   spread(type, value)
  
  ## To select the colors
  if(vrb == 'prec'){
    clr <- brewer.pal(n = 9, name = 'BrBG')
  } else if(vrb == 'etps'){
    print('etps')
    clr <-rev(brewer.pal(n = 9, name = 'BrBG'))
  } else {
    clr <- brewer.pal(n = 9, name = 'YlGnBu')
  }
  
  ## To draw the map 
  mnts <- month.abb
  
  map(.x = 1:12, .f = function(m){
    
    cat(mnts[m], '\t')
    tb <- filter(tbl, month_abb == mnts[m])
    
    gm <- ggplot() + 
      geom_tile(data = tb, aes(x = x, y = y, fill = value)) + 
      facet_wrap(.~type, ncol = 3) + 
      scale_fill_gradientn(colors = clr) +
      coord_sf() +
      labs(fill = vrb) +
      ggtitle(label = glue('Raw values - {mnts[m]}'), 
              subtitle = vrb) +
      theme_void() +
      theme(
        legend.position = 'bottom',
        legend.title = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line'),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        legend.title.position = 'top',
        strip.text = element_text(face = 'bold')
      )
    
    ggsave(plot = gm, filename = glue('./png/maps_climate/{vrb}-{mnts[m]}_periods.jpg'), units = 'in', width = 8, height = 5, dpi = 300)
    
    tb <- tb %>% 
      dplyr::select(-month, -type) %>% 
      spread(period, value) %>% 
      mutate(
        dfr_frc = Forecast - Baseline, 
        dfr_ftr = Future - Baseline
      ) %>% 
      dplyr::select(x, y, dfr_frc, dfr_ftr) %>% 
      gather(var, value, -c(x, y)) %>% 
      mutate(period = ifelse(var == 'dfr_frc', 'Difference Forecast', 'Difference Future'), 
             period = factor(period, levels = c('Difference Forecast', 'Difference Future')))
    
    gd <- ggplot() + 
      geom_tile(data = tb, aes(x = x, y = y, fill = value)) + 
      facet_wrap(.~period, ncol = 3) + 
      scale_fill_gradientn(colors = clr) +
      coord_sf() +
      labs(fill = vrb) +
      ggtitle(label = glue('Difference - {mnts[m]}'), 
              subtitle = vrb) +
      theme_void() +
      theme(
        legend.position = 'bottom',
        legend.title = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line'),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold'),
        legend.title.position = 'top',
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        strip.text = element_text(face = 'bold')
      )
    
    ggsave(plot = gd, filename = glue('./png/maps_climate/{vrb}-{mnts[m]}_difference-periods.jpg'), units = 'in', width = 8, height = 5, dpi = 300)
    rm(tb)
    gc(reset = T)
    cat('Saved maps\n')
      
  })
  
 cat('Done!\n')
  
  
}

## To make the maps
make.maps.vars(vrb = 'tmin')
make.maps.vars(vrb = 'tmax')
make.maps.vars(vrb = 'prec')
make.maps.vars(vrb = 'etps')

# To make the maps for the index ------------------------------------------

## Function 
make.maps.indx <- function(stk, idx, crp){
  
  stk <- ccoa
  idx <- 'HSH'
  # crp <- 'Cocoa'
  
  ## To start the analysis
  cat('To starth: ', idx, '\n')
  rst <- stk[[grep(idx, names(stk))]]
  
  ## Raster to table 
  tbl <- rst %>% 
    terra::as.data.frame(xy = T) %>% 
    as_tibble() %>% 
    gather(var, value, -c('x', 'y')) %>% 
    separate(data = ., col = 'var', into = c('Index', 'Period'), sep = '_') %>% 
    mutate(
      Period = ifelse(Period == 'bsl', 'Baseline', 'Future'), 
      Period = factor(Period, levels = c('Baseline', 'Future'))
    )
  
  ## To classify the values
  tbl <- tbl %>% 
      mutate(
        class = case_when(
          value == 1 ~ 'Low',
          value == 2 ~ 'Medium',
          value == 3 ~ 'High',
          TRUE ~ 'Very High'
      )
    )
  
  ## As factor column 
  tbl <- tbl %>% 
    mutate(
      class = factor(class, levels = c('Low', 'Medium', 'High', 'Very high'))
    )
  
  ## To select the colors 
  if(idx == 'HSH'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 4, name = 'Reds')
    ttle <- 'Human Heat Stress'
    sbtt <- 'Index'
  } else if(idx == 'NDWL0'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 3, name = 'BrBG')
    ttle <- 'Number of dry days with soil waterlogging at moisture content at start of saturarion or above'
    sbtt <- 'Days'
  } else if(idx == 'NTX30'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 3, name = 'YlOrRd')
    ttle <- 'Number of heat stress (Tmax > 30º)'
    sbtt <- 'Days'
  } else if(idx == 'TAI'){
    cat(idx, '\t')
    clrs <- rev(brewer.pal(n = 3, name = 'BrBG'))
    ttle <- 'Thornthwaite aridity index'
    sbtt <- '%'
  } else if(idx == 'NTX35'){
    cat(idx, '\t')
    clrs <- brewer.pal(n = 3, name = 'YlOrRd')
    ttle <- 'Number of heat stress (Tmax > 35º)'
    sbtt <- 'Days'
  } else if(idx == 'NDD'){
    cat(idx, '\t')
    clrs <- rev(brewer.pal(n = 3, name = 'BrBG'))
    ttle <- 'Number of dry days'
    sbtt <- 'Days'
  }
  
  ## To draw the map
  g.vls <- ggplot() + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = class)) + 
    facet_wrap(.~Period) + 
    scale_fill_manual(values = clrs, na.value = 'white') +
    geom_sf(data = st_as_sf(zne), fill = NA, col = 'grey30') +
    coord_sf() +
    labs(x = '', y = '', fill = idx) +
    ggtitle(label = ttle, 
            subtitle = crp) +
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
  ggsave(plot = g.vls, filename = glue('./png/maps_index_classified/values_{idx}_{crp}.jpg'), units = 'in', width = 10, height = 7, dpi = 300, create.dir = T)
  cat('Done!\n')
  
}

## To make the maps
indx <- c('NTX30', 'NTX35', 'NDD', 'NDWL0', 'TAI', 'HSH')

### Cocoa 
map(.x = 1:length(indx), .f = function(i) make.maps.indx(stk = ccoa, idx = indx[i], crp = 'Cocoa'))

### Coffee
map(.x = 1:length(indx), .f = function(i) make.maps.indx(stk = cffe, idx = indx[i], crp = 'Coffee'))


# Map for HSH -------------------------------------------------------------

stk <- ccoa
idx <- 'HSH'

## To start the analysis
cat('To starth: ', idx, '\n')
rst <- stk[[grep(idx, names(stk))]]

## Raster to table 
tbl <- rst %>% 
  terra::as.data.frame(xy = T) %>% 
  as_tibble() %>% 
  gather(var, value, -c('x', 'y')) %>% 
  separate(data = ., col = 'var', into = c('Index', 'Period'), sep = '_') %>% 
  mutate(
    Period = ifelse(Period == 'bsl', 'Baseline', 'Future'), 
    Period = factor(Period, levels = c('Baseline', 'Future'))
  )

## To classify the values
tbl <- tbl %>% 
  mutate(
    class = case_when(
      value == 1 ~ 'Low',
      value == 2 ~ 'Medium',
      value == 3 ~ 'High',
      TRUE ~ 'Very High'
    )
  )

## As factor column 
tbl <- tbl %>% 
  mutate(
    class = factor(class, levels = c('Low', 'Medium', 'High', 'Very High'))
  )

## To select the colors 
if(idx == 'HSH'){
  cat(idx, '\t')
  clrs <- brewer.pal(n = 4, name = 'Reds')
  ttle <- 'Human Heat Stress'
  sbtt <- 'Index'
} else if(idx == 'NDWL0'){
  cat(idx, '\t')
  clrs <- brewer.pal(n = 3, name = 'BrBG')
  ttle <- 'Number of dry days with soil waterlogging at moisture content at start of saturarion or above'
  sbtt <- 'Days'
} else if(idx == 'NTX30'){
  cat(idx, '\t')
  clrs <- brewer.pal(n = 3, name = 'YlOrRd')
  ttle <- 'Number of heat stress (Tmax > 30º)'
  sbtt <- 'Days'
} else if(idx == 'TAI'){
  cat(idx, '\t')
  clrs <- rev(brewer.pal(n = 3, name = 'BrBG'))
  ttle <- 'Thornthwaite aridity index'
  sbtt <- '%'
} else if(idx == 'NTX35'){
  cat(idx, '\t')
  clrs <- brewer.pal(n = 3, name = 'YlOrRd')
  ttle <- 'Number of heat stress (Tmax > 35º)'
  sbtt <- 'Days'
} else if(idx == 'NDD'){
  cat(idx, '\t')
  clrs <- rev(brewer.pal(n = 3, name = 'BrBG'))
  ttle <- 'Number of dry days'
  sbtt <- 'Days'
}

## To draw the map
g.vls <- ggplot() + 
  geom_tile(data = tbl, aes(x = x, y = y, fill = class)) + 
  facet_wrap(.~Period) + 
  scale_fill_manual(values = clrs, na.value = 'white') +
  geom_sf(data = st_as_sf(zne), fill = NA, col = 'grey30') +
  coord_sf() +
  labs(x = '', y = '', fill = idx) +
  ggtitle(label = ttle) +
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
ggsave(plot = g.vls, filename = glue('./png/maps_index_classified/values_{idx}.jpg'), units = 'in', width = 10, height = 7, dpi = 300, create.dir = T)
cat('Done!\n')





