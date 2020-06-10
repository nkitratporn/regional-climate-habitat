#[ status: on-going]
#[ note: 1. can improve using nest() and map() for functional and parallel processing (purrr + furrr)
#        2. data.table can be use as alternative to c++ function (have not tried)]
#@ This code is for calculating daily KBDI from raster images
#@ 1. convert rasters to data.frame/data.table for faster calculation
#@ 2. perform vectorize calculation over column
#@ 3. create recursive function for KBDI calculation using c++
#@ 5. perform statistical analysis
#@ 6. convert kbdi output / post-calculation results from data.frame/data.table to raster stack (if necessary)
#=====================================# 
# loading library
## data handling
library(tidyverse)

## raster and shapefile processing ##
library(rgeos)
library(rgdal)
library(raster)
library(sp)
#devtools::install_github("Pakillo/rgis")

## generating graph ##
#library(ggplot2)   # already loaded under tidyverse
library(maptools)
library(mapdata)

## c++ communication (lower load with recursive function/looping)
library(Rcpp)

#๒ parallel processing
library(foreach)
library(doParallel)


#=====================================#
## KBDI Calculation
#@ KBDI equation and symbol explained
########
#@ KBDI = KBDIy - (100  * NPR) + DF                                                             (eq.1)
#@ dQ   = (800.0-KBDIy)*(0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.441*prAnnaul)) (eq.2)
#@ dQ   = 0  if tasmax < 50F else dQ = dQ                                                       (eq.3)
########
#@ KBDI    : KBDI value at day t
#@ KBDIy   : KBDI value at day t-1
#@ NPR     : net rainfall at day t
#@ pr      : total 24-h rainfall at t in inch
#@ tasmax  : maximum temperature at t in Fahrenheit
#@ prAnnual: mean annual precipitation
#@ NF      : net precipitation
#@ dQ      : drought factor calculate using eq.2 and 3

#=====================================#
#Register CoreCluster
UseCores <- detectCores() - 10 #Define how many cores you want to use
cl <- makeCluster(UseCores)
registerDoParallel(cl)

## load polygon data
rangeDetail <- readOGR('data/shapefiles/rsStudy_dissolved.shp')
rangeCountry <- readOGR('data/shapefiles/rsStudy.shp')

## check number of life 
pr.files <- list.files('data/ERA5-mask/',pattern = '.*pr-([0-9]+).*')
tasmax.files <- list.files('data/ERA5-mask/',pattern = '.*tasmax-([0-9]+).*')

## load raster data

prAnnual <- stack('data/ERA5-mask/pr-Annual1995-2019.tif')                            # GEE calculation caused some slight mismatch  between raster
prAnnual <- projectRaster(prAnnual,stack(paste0('data/ERA5-mask/',pr.files[1]))[[1]]) # need to project to same extent and resolution

rangeDetail.r <- prAnnual
rangeDetail.r <- rasterize(rangeDetail,rangeDetail.r)

prAnnual <- mask(crop(prAnnual,rangeDetail),rangeDetail)


if(!dir.exists('data/climate-mask/pr-daily')){
  dir.create('data/climate-mask/pr-daily')
}
if(!dir.exists('data/climate-mask/tasmax-daily')){
  dir.create('data/climate-mask/tasmax-daily')
}
#writeRaster(prAnnual,'data/climate-mask/prAnnual.tif',overwrite=T)

prAnnual.df <- as.data.frame(prAnnual,xy=T)
prAnnual.df$location <- seq(1,nrow(prAnnual.df))

## create function for kbdi calculation
## solution 1 using C++ function
cppFunction('NumericVector kbdiC(NumericVector x, NumericVector y, NumericVector z){
	
	int n = x.size();
	double kbdi = x[1];
	NumericVector out(n);
    
    for (int i = 0; i < n; i++){
        kbdi = kbdi + ((800-kbdi)*y[i]) - (100*z[i]);
        if(kbdi < 0.0) { 
          kbdi = 0.0;
        }
        out[i] = kbdi;
    }
    return out;
}')

last.kbdi <- data.frame()

for (i in 2:length(pr.files)) {
  ## mask all raster
  start <- proc.time()[[3]]
  cat('load raster stack for ', sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'\n')
  if(i==1){
    pr <- stack(paste0('data/ERA5-mask/',pr.files[2]))
    tasmax <- stack(paste0('data/ERA5-mask/', tasmax.files[i]))
  }else{
    pr <- stack(paste0('data/ERA5-mask/',pr.files[2]))
    first.date <- as.Date(paste0(gsub(".*(\\d+{4}).*$", "\\1",names(pr[[1]])),'0101'), format = '%Y%m%d')
    
    
    tasmax <- stack(paste0('data/ERA5-mask/',tasmax.files[2]))
    tasmax.date <- first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",names(tasmax)))
  }
  pr <- setZ(pr,pr.date)
  cat('masked pr: ', nlayers(pr), ' COMPLETED.../ ')
  
  tasmax <- setZ(tasmax,tasmax.date)
  cat('masked tasmax: ', nlayers(tasmax), ' COMPLETED.../ ')
  
  
  #========CONVERT TO DATAFRAME===========#
  
  ### 1. convert to dataframe
  pr.df <- data.frame()
  tasmax.df <- data.frame()
  ## extract regional-scale mean and stddev into df
  pr.df <- foreach(i=1:nlayers(pr), .packages = 'raster',
                   .combine = 'cbind',.export = 'pr.df') %dopar% {
                     if (i > 1){
                       temp <- as.data.frame(pr[[i]],xy=T)
                       temp <- as.data.frame(temp[,-c(1,2)])
                       colnames(temp)[1] <- names(pr[[i]])
                     }else{
                       temp <- as.data.frame(pr[[i]],xy=T)
                       colnames(temp)[3] <- names(pr[[i]])
                     }
                     temp
                   }
  
  tasmax.df <- foreach(i=1:nlayers(tasmax), .packages = 'raster',
                       .combine = 'cbind',.export = 'tasmax.df') %dopar% {
                         if (i > 1){
                           temp <- as.data.frame(tasmax[[i]],xy=T)
                           temp <- as.data.frame(temp[,-c(1,2)])
                           colnames(temp)[1] <- names(tasmax[[i]])
                         }else{
                           temp <- as.data.frame(tasmax[[i]],xy=T)
                           colnames(temp)[3] <- names(tasmax[[i]])
                         }
                         temp
                       }
  
  ### 2. melt by column names and extract date from variable name, change column name
  if(i==1){
    pr.df.m <- pr.df %>% 
      reshape2::melt(id.vars = c('location','x','y')) %>%
      mutate(date = as.Date(gsub('.*X([0-9]+).*','\\1',variable), format = '%Y%m%d')) %>%
      rename(pr = value)
    
    tasmax.df.m <- tasmax.df %>%
      reshape2::melt(id.vars = c('location','x','y')) %>%
      mutate(date = as.Date(gsub('.*X([0-9]+).*','\\1',variable), format = '%Y%m%d')) %>%
      rename(tasmax = value)
  }else{
    pr.df.m <- pr.df %>% 
      reshape2::melt(id.vars = c('location','x','y')) %>%
      mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",variable))) %>%
      rename(pr = value)
    
    tasmax.df.m <- tasmax.df %>%
      reshape2::melt(id.vars = c('location','x','y')) %>%
      mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",variable))) %>%
      rename(tasmax = value)
  }
  
  
  prAnnual.df.m <- prAnnual.df %>%
    reshape2::melt(id.vars = c('location','x','y')) %>%
    rename(prAnnual = value)
  
  cat('df: COMPLETED.../ ')
  
  #========CALCULATION STEPS===========#
  # 1. calculate net precipitation (NF) dataframe
  #@ NPR is computed by subtracting 0.2 in. from the value of daily rainfall.
  #@ If there are consecutive wet days, 0.2 in. is subtracted only once on the day
  #@ when cumulative rainfall exceeds 0.2. 
  #@ A wet period ends when two rainy days are separated by one day without measurable rainfall,
  #@ thus 0.2 has to be subtracted again in the next rain period.
  
  ## calculate cumulative rainfall
  pr.df.m <- pr.df.m %>% 
    group_by(location) %>%
    arrange(date, .by_group = T) %>% 
    group_by(temp = cumsum(pr==0), .add=TRUE) %>%
    mutate(cRain = cumsum(pr)) %>%
    ungroup() %>%
    dplyr::select(-temp)
  
  ## set counter to reduce 0.2 in. start at 1 from cRain > 0.2 and reset when 0 again
  pr.df.m <- pr.df.m %>%
    mutate(counter = ifelse(cRain > 0.2,1,0)) %>%            #set 1 to anything cRain > 0.2
    group_by(location) %>%
    group_by(counterTemp = cumsum(counter==0),.add=TRUE) %>%  #group by location first and by counter == 0
    mutate(counter = cumsum(counter)) %>%                     #cumsum over group
    ungroup() %>%
    dplyr::select(-counterTemp)
  
  ## calculate net rainfall (NF)
  pr.df.m <- pr.df.m %>%
    mutate(nf = ifelse(counter==1, cRain-0.2,
                       ifelse(counter>1, pr,0)))
  
  cat('NF: COMPLETED.../ ')
  #head(pr.df.m[!is.na(pr.df.m$pr) & pr.df.m$location==25,],10)
  
  # 2. calculate kbdi
  ## join tasmax and prAnnual to pr
  pr.df.m <- pr.df.m %>%
    inner_join(dplyr::select(tasmax.df.m, c(location,date,tasmax)), by = c('location','date'))
  
  pr.df.m <- pr.df.m %>%
    inner_join(dplyr::select(prAnnual.df.m, c(location,prAnnual)), by = c('location'))
  
  #nrow(pr.df.m[!is.na(pr.df.m$pr),])
  #tail(pr.df.m[!is.na(pr.df.m$pr),],20)
  
  ## calculate PET part inside dQ formula
  #@ dQ   = 10^(-3)*(800Â–Q0) [0.968exp(0.0486T) Â– 8.3] dt / [1+10.88 exp (-0.0441 R)]
  pr.df.m <- pr.df.m %>%
    mutate(pet = (0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.0441*prAnnual)))
  
  if (nrow(last.kbdi)==0) {
    pr.df.m <- pr.df.m %>%
      add_column(kbdi_last = 0)
  }else{
    pr.df.m <- pr.df.m %>%
      inner_join(last.kbdi, by='location')
  }
  
  ## call kbdi function
  
  kbdi.df.m <- pr.df.m %>%
    group_by(location) %>%
    mutate(kbdi = kbdiC(kbdi_last,pet,nf)) %>%
    ungroup()
  
  cat('KBDI: COMPLETED.../ ')
  
  last.kbdi <- kbdi.df.m %>%
    group_by(location) %>%
    slice(n()) %>%
    dplyr::select(location,kbdi) %>%
    rename(kbdi_last = kbdi) %>%
    ungroup()
  
  write_csv(kbdi.df.m,path = paste0('output/kbdi-df/kbdi-',sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'.csv'),
            col_names = T)
  cat('save output: COMPLETED.../ ')
  rm(kbdi.df.m,pr.df.m,tasmax.df.m)
  end <- proc.time()[[3]]
  cat('elapse: ',end-start,'sec \n')
}

stopCluster(cl)
