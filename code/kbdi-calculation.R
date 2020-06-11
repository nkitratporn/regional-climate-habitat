#[ status: on-going]
#[ note: 1. use foreach to perform parallel processing when convert raster to dataframe
#           number of core to be used is usually (available.core-1), 
#           but since other tasks must be process on a shared server, only 85% was allocated (available.core - 5)
#        2. can improve using nest() and future_map() for functional and parallel processing with purrr & furrr pkgs (have not tried)
#        3. data.table can be use as alternative to c++ function (have not tried)]
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
#@ dQ      : drought factor calculate using eq.2
########
#@ @citation: 
#@  Keetch, J.J. and Byram, G.M. (1968) A drought index for forest fire control. USDA Forest Service.
#@  Alexander, M.E., 1990. Computer calculation of the Keetch-Byram Drought Index - programmers beware. Fire Management Notes 51, 23–25.
#=====================================#

#=====================================#
#       ERA5 HISTORICAL PERIOD        # 
#=====================================#
## check number of files
pr.files <- list.files('data/ERA5-mask/',pattern = '.*pr-([0-9]+).*')
tasmax.files <- list.files('data/ERA5-mask/',pattern = '.*tasmax-([0-9]+).*')

pr.files; 
tasmax.files;

## load static pr annual raster
prAnnual <- stack('data/ERA5-mask/pr-Annual1995-2019.tif')                            # GEE calculation caused some slight mismatch  between raster
prAnnual <- projectRaster(prAnnual,
                          stack(paste0('data/ERA5-mask/',pr.files[1]))[[1]]) # need to project to same extent and resolution

## convert to df
prAnnual.df <- as.data.frame(prAnnual,xy=T)
prAnnual.df$location <- seq(1,nrow(prAnnual.df))

## melt df and set relevant name
prAnnual.df.m <- prAnnual.df %>%
  reshape2::melt(id.vars = c('location','x','y')) %>%
  rename(prAnnual = value)

## create function for kbdi calculation
## solution 1 using C++ to speedup recusrive function
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

## declare empty last.kbdi in 1994 assuming kbdi = 0 in rainy season
last.kbdi <- data.frame()

for (i in 1:length(pr.files)) {
  start <- proc.time()[[3]]
  cat('load raster stack for ', sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'\n')
  
  pr <- stack(paste0('data/ERA5-mask/',pr.files[i]))
  tasmax <- stack(paste0('data/ERA5-mask/', tasmax.files[i]))
  
  if(i>1){
    first.date <- as.Date(paste0(gsub(".*(\\d+{4}).*$", "\\1",names(pr[[1]])),'0101'), format = '%Y%m%d')
  }
  
  #========CONVERT TO DATAFRAME===========#
  
  ### 1. convert to dataframe
  pr.df <- data.frame()
  tasmax.df <- data.frame()
  ## extract regional-scale mean and stddev into df
  
  #Register CoreCluster
  #in this code we check first how many core available and use around 85-90 %
  UseCores <- detectCores() - 5 #Define how many cores you want to use
  cl <- makeCluster(UseCores)
  registerDoParallel(cl)
  
  pr.df <- foreach(j=1:nlayers(pr), .packages = 'raster',
                   .combine = 'cbind',.export = 'pr.df') %dopar% {
                     if (j > 1){
                       temp <- as.data.frame(pr[[j]],xy=T)
                       temp <- as.data.frame(temp[,-c(1,2)])
                       colnames(temp)[1] <- names(pr[[j]])
                     }else{
                       temp <- as.data.frame(pr[[j]],xy=T)
                       colnames(temp)[3] <- names(pr[[j]])
                     }
                     temp
                   }
  
  tasmax.df <- foreach(j=1:nlayers(tasmax), .packages = 'raster',
                       .combine = 'cbind',.export = 'tasmax.df') %dopar% {
                         if (j > 1){
                           temp <- as.data.frame(tasmax[[j]],xy=T)
                           temp <- as.data.frame(temp[,-c(1,2)])
                           colnames(temp)[1] <- names(tasmax[[j]])
                         }else{
                           temp <- as.data.frame(tasmax[[j]],xy=T)
                           colnames(temp)[3] <- names(tasmax[[j]])
                         }
                         temp
                       }
  
  
  pr.df$location <- seq(1,nrow(pr.df))
  tasmax.df$location <- seq(1,nrow(tasmax.df))
  
  rm(pr,tasmax)
  
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
      rename(pr = value) %>%
      dplyr::filter(!is.na(pr))
    
    
    tasmax.df.m <- tasmax.df %>%
      reshape2::melt(id.vars = c('location','x','y')) %>%
      mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",variable))) %>%
      rename(tasmax = value) %>%
      dplyr::filter(!is.na(tasmax))
  }
  
  stopCluster(cl)
  
  rm(pr.df,tasmax.df)
  cat('df: COMPLETED... ')
  
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
  
  cat('NF: COMPLETED... ')
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
    cat('first year...')
  }else{
    pr.df.m <- pr.df.m %>%
      inner_join(last.kbdi, by='location')
  }
  
  ## call kbdi function
  
  kbdi.df.m <- pr.df.m %>%
    group_by(location) %>%
    mutate(kbdi = kbdiC(kbdi_last,pet,nf)) %>%
    ungroup() %>%
    dplyr::select(location,x,y,date,pr,tasmax,prAnnual,kbdi)
  
  cat('KBDI: COMPLETED... ')
  
  last.kbdi <- kbdi.df.m %>%
    group_by(location) %>%
    slice(n()) %>%
    dplyr::select(location,kbdi) %>%
    rename(kbdi_last = kbdi) %>%
    ungroup()
  
  rm(pr.df.m)
  
  write_csv(kbdi.df.m,path = paste0('output/kbdi-df/kbdi-',sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'.csv'),
            col_names = T)
  cat('save output: COMPLETED... ')
  #rm(kbdi.df.m,pr.df.m,tasmax.df.m)
  end <- proc.time()[[3]]
  cat('elapse: ',end-start,'sec \n')
}

#=====================================#
#       NEXGDDP FUTURE SCENARIO       # 
#=====================================#

## set scenario and models 
scenario <- c('rcp45','rcp85')
models <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR','CanESM2')

## prep for static variables/files
pr.files <- list.files('data/NEX-mask/',pattern = '.*pr-*')
tasmax.files <- list.files('data/NEX-mask/',pattern = '.*tasmax*')

prAnnual <- stack('data/ERA5-mask/pr-Annual1995-2019.tif')                            # GEE calculation caused some slight mismatch  between raster
prAnnual <- projectRaster(prAnnual,
                          stack(paste0('data/NEX-mask/',pr.files[1]))[[1]]) # need to project to same extent and resolution

prAnnual.df <- as.data.frame(prAnnual,xy=T)
prAnnual.df$location <- seq(1,nrow(prAnnual.df))

prAnnual.df.m <- prAnnual.df %>%
  reshape2::melt(id.vars = c('location','x','y')) %>%
  rename(prAnnual = value)

## set starting kbdi value from last day in historical analysis (2019-12-31)
last.kbdi2015 <- read_csv('output/kbdi-df/kbdi-2015.csv') %>%
  group_by(location) %>%
  slice(n()) %>%
  dplyr::select(location,kbdi) %>%
  rename(kbdi_last = kbdi) %>%
  ungroup()
head(last.kbdi2015)
for (k in seq(scenario)) {
  cat(scenario[k],': ')
  for (m in seq(models)) {
    cat(models[m],'\n')
    pr.select <- pr.files[grep(pr.files, pattern = paste0('.*',scenario[k],'.*m',m,'.*'))]
    tasmax.select <- tasmax.files[grep(tasmax.files, pattern = paste0('.*',scenario[k],'.*m',m,'.*'))]
    
    # test if number of files of pr and tasmax is equal
    if (length(pr.select)!=length(tasmax.select)) {
      cat('\n number of pr and tasmax files is not equal!')
      break
    }
    if (length(pr.select)==0 | length(tasmax.select)==0) {
      cat('\n no files!')
      break
    }
    
    for (n in 1:length(pr.select)) {
      start <- proc.time()[[3]]
      cat('load raster stacks: ', sub(".*(\\d+{4}).*$", "\\1", pr.select[n]),'\n')
      
      # load raster stacks
      pr <- stack(paste0('data/NEX-mask/', pr.select[n]))
      tasmax <- stack(paste0('data/NEX-mask/', tasmax.select[n]))
      
      # setup first date in stack
      first.date <- as.Date(paste0(gsub(".*(\\d+{4}).*$", "\\1",names(pr[[1]])),'0101'), format = '%Y%m%d')
      
      #========CONVERT TO DATAFRAME===========#
      
      ### 1. convert to dataframe
      pr.df <- data.frame()
      tasmax.df <- data.frame()
      
      #Register CoreCluster
      #in this code we check first how many core available and use around 85-90 %
      UseCores <- detectCores() - 5 #Define how many cores you want to use
      cl <- makeCluster(UseCores)
      registerDoParallel(cl)
      
      # convert raster to df in parallel 
      pr.df <- foreach(j=1:nlayers(pr), .packages = 'raster',
                       .combine = 'cbind',.export = 'pr.df') %dopar% {
                         if (j > 1){
                           temp <- as.data.frame(pr[[j]],xy=T)
                           temp <- as.data.frame(temp[,-c(1,2)])
                           colnames(temp)[1] <- names(pr[[j]])
                         }else{
                           temp <- as.data.frame(pr[[j]],xy=T)
                           colnames(temp)[3] <- names(pr[[j]])
                         }
                         temp
                       }
      
      tasmax.df <- foreach(j=1:nlayers(tasmax), .packages = 'raster',
                           .combine = 'cbind',.export = 'tasmax.df') %dopar% {
                             if (j > 1){
                               temp <- as.data.frame(tasmax[[j]],xy=T)
                               temp <- as.data.frame(temp[,-c(1,2)])
                               colnames(temp)[1] <- names(tasmax[[j]])
                             }else{
                               temp <- as.data.frame(tasmax[[i]],xy=T)
                               colnames(temp)[3] <- names(tasmax[[j]])
                             }
                             temp
                           }
      
      # add location id
      pr.df$location <- seq(1,nrow(pr.df))
      tasmax.df$location <- seq(1,nrow(tasmax.df))
      
      # remove raster stack (reduce memory occupency)
      rm(pr,tasmax)
      
      ### 2. melt by column names and extract date from variable name, change column name
      pr.df.m <- pr.df %>% 
        reshape2::melt(id.vars = c('location','x','y')) %>%
        mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){4}([^.]+).*", "\\1",variable))) %>%
        rename(pr = value) %>%
        dplyr::filter(!is.na(pr))
      
      
      tasmax.df.m <- tasmax.df %>%
        reshape2::melt(id.vars = c('location','x','y')) %>%
        mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){4}([^.]+).*", "\\1",variable))) %>%
        rename(tasmax = value) %>%
        dplyr::filter(!is.na(tasmax))
      
      stopCluster(cl)
      
      rm(pr.df,tasmax.df)
      cat('df: COMPLETED... ')
      
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
      
      cat('NF: COMPLETED... ')
      
      # 2. calculate kbdi
      ## join tasmax and prAnnual to pr
      pr.df.m <- pr.df.m %>%
        inner_join(dplyr::select(tasmax.df.m, c(location,date,tasmax)), by = c('location','date'))
      
      pr.df.m <- pr.df.m %>%
        inner_join(dplyr::select(prAnnual.df.m, c(location,prAnnual)), by = c('location'))
      
      
      ## calculate PET part inside dQ formula
      #@ dQ   = 10^(-3)*(800Ã‚Â–Q0) [0.968exp(0.0486T) Ã‚Â– 8.3] dt / [1+10.88 exp (-0.0441 R)]
      pr.df.m <- pr.df.m %>%
        mutate(pet = (0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.0441*prAnnual)))
      
      if(n==1){
        last.kbdi <- last.kbdi2015
      }
      
      pr.df.m <- pr.df.m %>%
        inner_join(last.kbdi, by='location')
      
      ## call kbdi function
      
      kbdi.df.m <- pr.df.m %>%
        group_by(location) %>%
        mutate(kbdi = kbdiC(kbdi_last,pet,nf)) %>%
        ungroup() %>%
        dplyr::select(location,x,y,date,pr,tasmax,prAnnual,kbdi)
      
      cat('KBDI: COMPLETED... ')
      
      last.kbdi <- kbdi.df.m %>%
        group_by(location) %>%
        slice(n()) %>%
        dplyr::select(location,kbdi) %>%
        rename(kbdi_last = kbdi) %>%
        ungroup()
      
      rm(pr.df.m)
      
      write_csv(kbdi.df.m,path = paste0('output/kbdi-df/kbdi-',scenario[k],'-m',m,'-',sub(".*(\\d+{4}).*$", "\\1", pr.select[n]),'.csv'),
                col_names = T)
      cat('save output: COMPLETED... ')
      #rm(kbdi.df.m,pr.df.m,tasmax.df.m)
      end <- proc.time()[[3]]
      cat('elapse: ',end-start,'sec \n')
    }
  }
}
