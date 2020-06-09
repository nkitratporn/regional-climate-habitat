#[ status: on-going]
#[ note: 1. can improve using nest() and map() for functional programming
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

## generating graph ##
#library(ggplot2)   # already loaded under tidyverse
library(maptools)
library(mapdata)

## c++ communication (lower load with recursive function/looping)
library(Rcpp)

#=====================================#
## KBDI Calculation
#@ KBDI equation and symbol explained
########
#@ KBDI = KBDIy - (100  * NPR) + DF                                                             (eq.1)
#@ dQ   = (800.0-KBDIy)*(0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.441*prAnnaul)) (eq.2)
#@ dQ   = 0  if tasmax < 50F else dQ = dQ                                                       (eq.3)
########
#@ KBDI   : KBDI value at day t
#@ KBDIy  : KBDI value at day t-1
#@ NPR    : net rainfall at day t
#@ pr     : total 24-h rainfall at t in inch
#@ tasmax : maximum temperature at t in Fahrenheit
#@ NF     : net precipitation
#@ dQ     : drought factor calculate using eq.2 and 3

#=====================================#
## load polygon data
rangeDetail <- readOGR('data/shapefiles/rsStudy_dissolved.shp')
rangeCountry <- readOGR('data/shapefiles/rsStudy.shp')

## load raster data
pr1994 <- stack('data/pr-1994.tif')
tasmax1994 <- stack('data/tasmax-1994.tif')
prAnnual <- stack('data/prAnnual.tif')
## mask all raster
pr1994 <- mask(crop(pr1994,rangeDetail),rangeDetail)
pr.date <- as.Date(gsub('.*X([0-9]+).*','\\1',names(pr1994)), format = '%Y%m%d')
pr1994 <- setZ(pr1994,pr.date)

tasmax1994 <- mask(crop(tasmax1994,rangeDetail),rangeDetail)
tasmax.date <- as.Date(gsub('.*X([0-9]+).*','\\1',names(tasmax1994)), format = '%Y%m%d')
tasmax1994 <- setZ(tasmax1994,tasmax.date)

prAnnual <- mask(crop(prAnnual,rangeDetail),rangeDetail)

## convert to dataframe
### A. use directly for further calculation
pr.df <- as.data.frame(pr1994,xy=T)
pr.df$location <- seq(1,nrow(pr.df))

tasmax.df <- as.data.frame(tasmax1994,xy=T)
tasmax.df$location <- seq(1,nrow(tasmax.df))

prAnnual.df <- as.data.frame(prAnnual,xy=T)
prAnnual.df$location <- seq(q,nrow(pr.df))

### B. melt by column names and change julian day to date and rbind to existing data
###    combine all to create daily time series graph later
pr.df.m <- pr.df %>% 
  reshape2::melt(id.vars = c('location','x','y')) %>%
  mutate(date = as.Date(gsub('.*X([0-9]+).*','\\1',variable), format = '%Y%m%d')) %>%
  rename(pr = value)

tasmax.df.m <- tasmax.df %>%
  reshape2::melt(id.vars = c('location','x','y')) %>%
  mutate(date = as.Date(gsub('.*X([0-9]+).*','\\1',variable), format = '%Y%m%d')) %>%
  rename(tasmax = value)

prAnnual.df.m <- prAnnual.df %>%
  reshape2::melt(id.vars = c('location','x','y')) %>%
  rename(prAnnual = value)

# CALCULATION STEPS
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
  select(-temp)

## set counter to reduce 0.2 in. start at 1 from cRain > 0.2 and reset when 0 again
pr.df.m <- pr.df.m %>%
  mutate(counter = if_else(cRain > 0.2,1,0)) %>%            #set 1 to anything cRain > 0.2
  group_by(location) %>%
  group_by(counterTemp = cumsum(counter==0),.add=TRUE) %>% #group by location first and by counter == 0
  mutate(counter = cumsum(counter)) %>%                    #cumsum over group
  ungroup() %>%
  select(-counterTemp)

## calculate net rainfall (NF)
pr.df.m <- pr.df.m %>%
  mutate(nf = if_else(counter==1, cRain-0.2,
                     if_else(counter>1, pr,0)))

tail(pr.df.m[!is.na(pr.df.m$pr) & pr.df.m$location==30,],10)

# 2. calculate kbdi

## join tasmax and prAnnual to pr
pr.df.m <- pr.df.m %>%
  full_join(tasmax.df.m, by = c('location','date')) %>%
  select(-x.y, -y.y, -variable.y)

pr.df.m <- pr.df.m %>%
  full_join(prAnnual.df.m, by = c('location')) %>%
  select(-x.y, -y.y, -variable.y)

head(pr.df.m[!is.na(pr.df.m$pr) & pr.df.m$location == 25,],20)

## calculate PET part inside dQ formula
#@ dQ   = (800.0-KBDIy)*PET; PET = (0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.441*prAnnaul))
pr.df.m <- pr.df.m %>%
  group_by(location) %>%
  mutate(pet = (0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.441*prAnnaul))) %>%
  ungroup()

pr.df.m <- pr.df.m %>%
  add_column(tmp = 0)

first(pr.df.m$tmp) <- last.kbdi

## create function for kbdi calculation
## solution 1 using C++ function
cppFunction('NumericVector kbdiC(NumericVector x, NumericVector y, NumericVector z){
	
	int n = x.size();
	double kbdi = x[1];
	NumericVector out(n);
    
    for (i = 0; i < n; i++){
        kbdi = kbdi + ((800-kbdi)*y[i]) - (100*z[i]);
        if(kbdi < = ) { kbdi = 0}
        out[i] = kbdi;
    }
    return out;
}')

## call kbdi function

kbdi.df.m <- pr.df.m %>%
  group_by(location) %>%
  mutate(kbdi = kbdiC(pet,nf,last.kbdi)) %>%
  ungroup()

last.kbdi = last(kbdi.df.m$kbdi)