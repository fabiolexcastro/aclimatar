

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

## 
mskr <- dir_ls('./common_data/chirps_cmip6_america/Prec_ACCESS-ESM1-5_ssp245_2021_2040')[1] %>% rast()
mskr <- mskr * 0 + 1

## 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- filter(wrld, sov_a3 %in% c('SLV', 'GTM', 'MEX', 'NIC', 'HND', 'PAN', 'CRI', 'COL', 'ECU', 'PER'))

# Sample points -----------------------------------------------------------
pnts <- map_dfr(.x = 1:nrow(zone), .f = function(i){
  zne <- zone[i,]
  zne <- vect(zne)
  msk <- terra::crop(mskr, zne)
  msk <- terra::mask(msk, zne)
  pnt <- as.data.frame(raptr::randomPoints(mask = msk, n = 1))
  pnt <- mutate(pnt, iso = zne$sov_a3)
  return(pnt)
})
# write.csv(pnts, './tble/points_sample.csv', row.names = FALSE)
pnts <- read_csv('./tble/points_sample.csv', show_col_types = FALSE)

gpnt <- ggplot() +
  geom_sf(data = zone, fill = NA, col = 'grey30') +
  geom_point(data = pnts, aes(x = x, y = y), col = 'red') +
  ggtitle(label = 'Puntos para validación de datos') +
  labs(x = '', y = '') +
  theme_minimal() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5)
  )

ggsave(plot = gpnt, filename = './png/maps_raw/points.jpg', units = 'in', width = 9, height = 7, dpi = 300, create.dir = T)

# Get baseline values -----------------------------------------------------

## Precipitation
root.bsln <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps'
fles.bsln <- as.character(dir_ls(root.bsln))
year.bsln <- paste0('v2.0.', 1995:2014, '.')
fles.bsln <- grep(paste0(year.bsln, collapse = '|'), fles.bsln, value = T)

### Function 
get.value_bsl <- function(year){
  
  # year <- year.bsln[1]
  
  cat('>>>>>', year, '\n')
  fles <- grep(year, fles.bsln, value = T)
  yr   <- str_sub(year, 6, 9)
  
  ## Classify by each month
  rstr <- map(.x = 1:12, .f = function(m){
    
    ## To read as a raster
    cat('To process: ', m, '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    fls <- grep(paste0(year, mnt), fles, value = T)
    rst <- rast(fls)
    
    ## To extract by mask 
    rst <- terra::crop(rst, zone)
    rst <- terra::mask(rst, zone)
    
    ## Classify 
    mtx <- matrix(c(-Inf, 0, 0, Inf, NA), ncol = 3, byrow = T)
    rcl <- terra::classify(rst, mtx, include.lowest = T, right = FALSE)
    
    ## To sum
    sma <- sum(rcl)
    names(sma) <- glue('prec_{yr}-{m}')
    
    ## To write daily raster 
    terra::writeRaster(x = rcl, filename = glue('./common_data/chirps_hist_america/daily/Prec_{yr}-{m}.tif'), overwrite = TRUE)
    
    ## Finish 
    return(sma)
    
  })
  
  ### To write raster
  rstr <- reduce(rstr, c)
  dout <- glue('./common_data/chirps_hist_america/monthly')
  # terra::writeRaster(x = rstr, filename = glue('{dout}/Prec_{yr}.tif'), overwrite = TRUE)
  
  ### To extract the values
  vles <- dplyr::select(cbind(pnts, terra::extract(rstr, pnts[,1:2])), -ID)
  
  ### To remove some objects 
  rm(rstr, dout, fles, yr)
  return(vles)
  
}

### Apply the function (prec - baseline)
tble.prec <- map(year.bsln, get.value_bsl)
tble.prec <- reduce(tble.prec, inner_join, by = c('x', 'y', 'iso'))
tble.prec <- as_tibble(tble.prec)
tble.prec <- tble.prec %>% gather(var, value, -c(x, y, iso))
tble.prec <- separate(data = tble.prec, col = 'var', into = c('variable', 'date'), sep = '_', remove = TRUE)

## Temperature 
root.bsln <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global'
dirs.bsln <- paste0(root.bsln, c('Tmin', 'Tmax'))

### Function 
get.value_bsl <- function(year, varb){
  
  # year <- '1995'
  # varb <- 'Tmin'
  
  ## To list the files and grepping
  cat('To process: ', year, ' ', varb, '\n')
  fles <- as.character(dir_ls(glue('{root.bsln}/{varb}'), regexp = '.tif$'))
  fles <- grep(paste0('.', year, '.'), fles, value = T)
  
  ## Classify by each month
  rstr <- map(.x = 1:12, .f = function(m){
    
    ## To read as a raster
    cat('To process: ', m, '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    fls <- grep(paste0(year, '.', mnt), fles, value = T)
    rst <- rast(fls)
    
    ## To extract by mask 
    rst <- terra::crop(rst, zone)
    rst <- terra::mask(rst, zone)
    
    ## Classify 
    rcl <- terra::subst(rst, -9999, NA)
    
    ## To sum
    avg <- mean(rcl)
    names(avg) <- glue('{varb}_{year}-{m}')
    
    ## To write the raster 
    out <- glue('./common_data/chirts_hist_america/daily/{varb}')
    terra::writeRaster(x = rcl, filename = glue('{out}/{varb}_{year}-{m}.tif'), overwrite = TRUE)
    
    ## Finish 
    return(avg)
    
  })
  rstr <- reduce(rstr, c)
  dout <- glue('./common_data/chirts_hist_america/monthly/{varb}/{varb}_{year}-{m}.tif')
  terra::writeRaster(x = rstr, filename = dout, overwrite = TRUE)
  
  ## To extract the values
  vles <- dplyr::select(cbind(pnts, terra::extract(rstr, pnts[,1:2])), -ID)
  cat('Done!\n')
  return(vles)
  
}

### Apply the function (tmin - tmax / baseline)
stpn <- expand.grid(year = 1995:2014, vars = c('Tmin', 'Tmax'))
tble.tasm <- 1:nrow(stpn) %>% map(.f = function(i){get.value_bsl(year = stpn$year[i], varb = stpn$vars[i])})
tble.tasm.2 <- tble.tasm %>% reduce(., inner_join, by = c('x', 'y', 'iso'))
tble.tasm.2 <- as_tibble(tble.tasm.2)
tble.tasm.2 <- tble.tasm.2 %>% gather(var, value, -c(x, y, iso))
tble.tasm.2 <- separate(data = tble.tasm.2, col = 'var', into = c('variable', 'date'), sep = '_', remove = TRUE)
tble.tasm.2 <- tble.tasm.2 %>% spread(variable, value)
write.csv(tble.tasm.2, './tble/points_values-tasm_hist.csv', row.names = FALSE)

## Join prec and temperature for the baseline
tble.prec <- tble.prec %>% spread(variable, value)
tble.bsln <- inner_join(tble.prec, tble.tasm.2, by = c('x', 'y', 'iso', 'date'))
write.csv(tble.bsln, './tble/points_values-prec-tmin-tmax_hist.csv', row.names = FALSE)

tble.bsln <- read_csv('./tble/points_values-prec-tmin-tmax_hist.csv', show_col_types = FALSE)

# Get future values -------------------------------------------------------

## Function
get.values_ftr <- function(sspe){
  
  # sspe <- 'ssp370'
  
  ## To start
  cat('To process: ', sspe, '\n')
  pths <- as.character(dir_ls(grep(sspe, dirs, value = T)))
  pths <- grep(paste0(c('PTOT', 'TMIN', 'TMAX'), collapse = '|'), pths, value = T)
  
  outfile <- glue('./tble/points_values-prec-tmin-tmax_{sspe}.csv')
  
  if(!file.exists(outfile)){
    
    tble.gcms <- map(.x = 1:length(gcms), .f = function(i){
      
      ## To list the files
      cat('>>>>>> ', gcms[i], '\n')
      drs <- grep(gcms[i], pths, value = T)
      fls <- map(drs, dir_ls)
      fls <- unlist(fls)
      fls <- as.character(fls)
      
      ## To read as raster
      ppt <- grep('PTOT', fls, value = T)
      tmn <- grep('TMIN', fls, value = T)
      tmx <- grep('TMAX', fls, value = T)
      
      yrs <- 2021:2100
      
      vls <- map(.x = yrs, .f = function(y){
        
        cat('>>> Year: ', y, '\n')
        pt <- rast(grep(y, ppt, value = T))
        tn <- rast(grep(y, tmn, value = T))
        tx <- rast(grep(y, tmx, value = T))
        
        names(pt) <- glue('prec_{y}-{1:12}')
        names(tn) <- glue('tmin_{y}-{1:12}')
        names(tx) <- glue('tmax_{y}-{1:12}')
        
        vl <- list(
          terra::extract(pt, pnts[,1:2]),
          terra::extract(tn, pnts[,1:2]), 
          terra::extract(tx, pnts[,1:2])
        ) %>% 
          reduce(., inner_join, by = 'ID')
        
        rm(pt, tn, tx) 
        gc(reset = T)
        
        return(vl)
        
      })
      
      ## Tidy the final table
      vls <- vls %>% reduce(., inner_join, by = c('ID'))
      vls <- vls %>% gather(var, value, -c(ID))
      vls <- vls %>% separate(data = ., col = 'var', into = c('variable', 'date'), sep = '_', remove = T)
      vls <- vls %>% mutate(ssp = sspe, gcme = gcms[i])
      vls <- vls %>% as_tibble()
      
      ## Finish 
      cat('Done!\n')
      return(vls)
      
    })
    
    tble.gcms <- bind_rows(tble.gcms)
    write.csv(tble.gcms, outfile, row.names = FALSE)
    cat('Done!\n')
    
  } else {
    
    tble.gcms <- read_csv(outfile, show_col_types = FALSE)
    
  }
  
  return(tble.gcms)
  
}

## To apply the function 
tble.ssp245 <- get.values_ftr(sspe = 'ssp245')
tble.ssp245 <- tble.ssp245 %>% spread(variable, value)
tble.ssp370 <- get.values_ftr(sspe = 'ssp370')
tble.ssp370 <- tble.ssp370 %>% spread(variable, value)


# Tidy the tables ---------------------------------------------------------
pnts <- mutate(pnts, ID = 1:10)[,3:4]
tble.bsln <- inner_join(tble.bsln, pnts, by = 'iso')
tble.bsln <- tble.bsln %>% arrange(ID)

tble.ssp245 <- inner_join(tble.ssp245, pnts, by = 'ID')
tble.ssp370 <- inner_join(tble.ssp370, pnts, by = 'ID')

write.csv(tble.bsln, './tble/values/vls_hist.csv', row.names = FALSE)
write.csv(tble.ssp245, './tble/values/vls_ssp245.csv', row.names = FALSE)
write.csv(tble.ssp370, './tble/values/vls_ssp370.csv', row.names = FALSE)

tble.bsln <- read_csv('./tble/values/vls_hist.csv')
tble.ssp245 <- read_csv('./tble/values/vls_ssp245.csv')
tble.ssp370 <- read_csv('./tble/values/vls_ssp370.csv')

# To compare the values ---------------------------------------------------
tble.bsln <- tble.bsln %>% rename(prec = prec, tmax = Tmax, tmin = Tmin)
tble.s245 <- tble.ssp245
tble.s370 <- tble.ssp370

# To calculate the average ------------------------------------------------
smmr.s245 <- tble.s245 %>% 
  group_by(ID, iso, date, ssp) %>% 
  reframe(
    prec = mean(prec), 
    tmin = mean(tmin), 
    tmax = mean(tmax)
  )
smmr.s370 <- tble.s370 %>% 
  group_by(ID, iso, date, ssp) %>% 
  reframe(
    prec = mean(prec), 
    tmin = mean(tmin), 
    tmax = mean(tmax)
  )

## Function 
calc.dfrn <- function(x, y){
  
  x <- tble.bsln
  y <- smmr.s370
  
  ## Rename
  cat('Start\n')
  x <- x %>% rename(prec_bsln = prec, tmax_bsln = tmax, tmin_bsln = tmin)
  x <- x %>% mutate(month = as.numeric(str_sub(date, 6, nchar(date))))
  y <- y %>% mutate(month = as.numeric(str_sub(date, 6, nchar(date))))
  
  ## Summarise
  x <- x %>% 
    group_by(x, y, ID, iso, month) %>% 
    reframe(prec_bsln = mean(prec_bsln), tmin_bsln = mean(tmin_bsln), tmax_bsln = mean(tmax_bsln)) %>% 
    arrange(ID)
  y <- y %>% 
    group_by(ID, iso, month, ssp) %>% 
    reframe(prec = mean(prec), tmin = mean(tmin), tmax = mean(tmax))
  
  ## To make the join and calculate the difference
  xy <- inner_join(x, y)
  xy <- mutate(xy, prec_dfrn = prec - prec_bsln, tmin_dfrn = tmin - tmin_bsln, tmax_dfrn = tmax - tmax_bsln)
  
  ## Return
  cat('Finish!\n')
  return(xy)
  
}

## To calculate the difference 
dfrn.s245 <- calc.dfrn(x = tble.bsln, y = smmr.s245)
dfrn.s370 <- calc.dfrn(x = tble.bsln, y = smmr.s370)
dfrn.s370

# To make the graphs ------------------------------------------------------

make.graph <- function(tble){
  
  tble <- dfrn.s370
  
  ## 
  cat('To process!\n')
  prec <- tble %>% 
    dplyr::select(ID, iso, month, prec_bsln, prec, prec_dfrn) %>% 
    rename(prec_ftre = prec) %>% 
    gather(var, value, -c(ID, iso, month)) %>% 
    separate(data = ., col = 'var', into = c('variable', 'period'), sep = '_')
  tmin <- tble %>% 
    dplyr::select(ID, iso, month, tmin_bsln, tmin, tmin_dfrn) %>% 
    rename(tmin_ftre = tmin) %>% 
    gather(var, value, -c(ID, iso, month)) %>% 
    separate(data = ., col = 'var', into = c('variable', 'period'), sep = '_')
  tmax <- tble %>% 
    dplyr::select(ID, iso, month, tmax_bsln, tmax, tmax_dfrn) %>% 
    rename(tmax_ftre = tmax) %>% 
    gather(var, value, -c(ID, iso, month)) %>% 
    separate(data = ., col = 'var', into = c('variable', 'period'), sep = '_')
  
  tmin <- tmin %>% mutate(period = factor(period, levels = c('bsln', 'ftre', 'dfrn')))
  tmax <- tmax %>% mutate(period = factor(period, levels = c('bsln', 'ftre', 'dfrn')))
  
  sspe <- unique(tble$ssp)
  
  ## To draw the graphs
  gprec <- ggplot(data = filter(prec, period %in% c('bsln', 'ftre')), aes(x = factor(month), y = value, fill = period)) + 
    geom_col(stat = 'identity', position = 'dodge') +
    facet_wrap(.~iso, scales = 'free_y') +
    scale_fill_brewer(palette = "Set2") +
    scale_x_discrete(labels = month.abb) +
    ggtitle(label = 'Precipitation (mm)') +
    labs(x = '', y = 'Prec (mm)', fill = '') +
    theme_minimal() +
    theme(
      plot.title = element_text(face = 'bold', hjust = 0.5), 
      legend.position = 'bottom', 
      axis.text.x = element_text(angle = 90, hjust = 0.5), 
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    )
  
  colors <- c("bsln" = "#bfb966", "ftre" = "#66bf79", "dfrn" = "#e5632f")
  
  gtmin <- ggplot(tmin, aes(x = month)) +
    geom_col(
      data = filter(tmin, period != "dfrn"),
      aes(y = value, fill = period),
      position = position_dodge(width = 0.7),
      width = 0.6,
      alpha = 0.8
    ) +
    geom_line(
      data = filter(tmin, period == "dfrn"),
      aes(y = value, color = period),
      linewidth = 1.2,
      group = 1
    ) +
    geom_point(
      data = filter(tmin, period == "dfrn"),
      aes(y = value, color = period),
      size = 3
    ) +
    facet_wrap(.~iso) +
    # Escalas y temas
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(
      x = "",
      y = "Temperatura Mínima (°C)",
      title = "Temperatura Mínima Mensual: Línea Base vs Futuro",
      subtitle = "Barras: bsln (presente) y ftre (futuro) | Línea naranja: Diferencia (dfrn)",
      fill = "Periodo",
      color = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "top", 
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 0.5)
    )
  
  gtmax <- ggplot(tmax, aes(x = month)) +
    geom_col(
      data = filter(tmax, period != "dfrn"),
      aes(y = value, fill = period),
      position = position_dodge(width = 0.7),
      width = 0.6,
      alpha = 0.8
    ) +
    geom_line(
      data = filter(tmax, period == "dfrn"),
      aes(y = value, color = period),
      linewidth = 1.2,
      group = 1
    ) +
    geom_point(
      data = filter(tmax, period == "dfrn"),
      aes(y = value, color = period),
      size = 3
    ) +
    facet_wrap(.~iso) +
    # Escalas y temas
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(
      x = "",
      y = "Temperatura Máxima (°C)",
      title = "Temperatura Máxima Mensual: Línea Base vs Futuro",
      subtitle = "Barras: bsln (presente) y ftre (futuro) | Línea naranja: Diferencia (dfrn)",
      fill = "Periodo",
      color = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "top", 
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 0.5)
    )
  
  ## To save the graphs
  ggsave(plot = gprec, filename = glue('./png/graphs_difference/difference_{sspe}_prec.jpg'), units = 'in', width = 12, height = 12, dpi = 300, create.dir = T)
  ggsave(plot = gtmin, filename = glue('./png/graphs_difference/difference_{sspe}_tmin.jpg'), units = 'in', width = 12, height = 12, dpi = 300, create.dir = T)
  ggsave(plot = gtmax, filename = glue('./png/graphs_difference/difference_{sspe}_tmax.jpg'), units = 'in', width = 12, height = 12, dpi = 300, create.dir = T)
  
}







