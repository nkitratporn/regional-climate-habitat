# loading library
## data handling ##
library(tidyverse)
library(dplyr)

## raster and shapefile processing ##
library(rgeos)
library(rgdal)
library(raster)
library(sp)

## generating graph ##
library(ggplot2)
library(reshape2)
library(maptools)
library(lattice)
library(rasterVis)

# ========= ERA5 =========== #
## loading range 
rangeDetail <- readOGR('data/shapefiles/rsStudy_dissolved.shp')
rangeCountry <- readOGR('data/shapefiles/rsStudy.shp')

#Register CoreCluster
UseCores <- detectCores() -1 #Define how many cores you want to use
cl <- makeCluster(UseCores)
registerDoParallel(cl)

## loading raw daily variables
### annual pr
prQuarter.era <- stack(paste0('data/ERA-quarter-stats/',
                           list.files('data/ERA-quarter-stats/',pattern = 'ERA-pr-*')))
prQuarter.nex <- stack(paste0('data/NEX-quarter-stats/',
                           list.files('data/NEX-quarter-stats/',pattern = '.*-pr-*')))

prQuarter.era <- mask(crop(prQuarter.era,rangeDetail),rangeDetail)
prQuarter.nex <- mask(crop(prQuarter.nex,rangeDetail),rangeDetail)

names(prQuarter.era) <- c('ERA5_DJF','ERA5_JJA','ERA5_MAM','ERA5_SON')
names(prQuarter.nex) <- c('NEX_rcp45m1p1q1_DJF','NEX_rcp45m1p1q3_JJA','NEX_rcp45m1p1q2_MAM','NEX_rcp45m1p1q4_SON',
                          'NEX_rcp45m1p0q1_DJF','NEX_rcp45m1p0q3_JJA','NEX_rcp45m1p0q2_MAM','NEX_rcp45m1p0q4_SON',
                          'NEX_rcp45m2p1q1_DJF','NEX_rcp45m2p1q3_JJA','NEX_rcp45m2p1q2_MAM','NEX_rcp45m2p1q4_SON',
                          'NEX_rcp45m2p0q1_DJF','NEX_rcp45m2p0q3_JJA','NEX_rcp45m2p0q2_MAM','NEX_rcp45m2p0q4_SON',
                          'NEX_rcp45m3p1q1_DJF','NEX_rcp45m3p1q3_JJA','NEX_rcp45m3p1q2_MAM','NEX_rcp45m3p1q4_SON',
                          'NEX_rcp45m3p0q1_DJF','NEX_rcp45m3p0q3_JJA','NEX_rcp45m3p0q2_MAM','NEX_rcp45m3p0q4_SON',
                          'NEX_rcp45m4p1q1_DJF','NEX_rcp45m4p1q3_JJA','NEX_rcp45m4p1q2_MAM','NEX_rcp45m4p1q4_SON',
                          'NEX_rcp45m4p0q1_DJF','NEX_rcp45m4p0q3_JJA','NEX_rcp45m4p0q2_MAM','NEX_rcp45m4p0q4_SON',
                          'NEX_rcp45m5p1q1_DJF','NEX_rcp45m5p1q3_JJA','NEX_rcp45m5p1q2_MAM','NEX_rcp45m5p1q4_SON',
                          'NEX_rcp45m5p0q1_DJF','NEX_rcp45m5p0q3_JJA','NEX_rcp45m5p0q2_MAM','NEX_rcp45m5p0q4_SON',
                          'NEX_rcp45m6p1q1_DJF','NEX_rcp45m6p1q3_JJA','NEX_rcp45m6p1q2_MAM','NEX_rcp45m6p1q4_SON',
                          'NEX_rcp45m6p0q1_DJF','NEX_rcp45m6p0q3_JJA','NEX_rcp45m6p0q2_MAM','NEX_rcp45m6p0q4_SON',                          
                          'NEX_rcp85m1p1q1_DJF','NEX_rcp85m1p1q3_JJA','NEX_rcp85m1p1q2_MAM','NEX_rcp85m1p1q4_SON',
                          'NEX_rcp85m1p0q1_DJF','NEX_rcp85m1p0q3_JJA','NEX_rcp85m1p0q2_MAM','NEX_rcp85m1p0q4_SON',
                          'NEX_rcp85m2p1q1_DJF','NEX_rcp85m2p1q3_JJA','NEX_rcp85m2p1q2_MAM','NEX_rcp85m2p1q4_SON',
                          'NEX_rcp85m2p0q1_DJF','NEX_rcp85m2p0q3_JJA','NEX_rcp85m2p0q2_MAM','NEX_rcp85m2p0q4_SON',
                          'NEX_rcp85m3p1q1_DJF','NEX_rcp85m3p1q3_JJA','NEX_rcp85m3p1q2_MAM','NEX_rcp85m3p1q4_SON',
                          'NEX_rcp85m3p0q1_DJF','NEX_rcp85m3p0q3_JJA','NEX_rcp85m3p0q2_MAM','NEX_rcp85m3p0q4_SON',
                          'NEX_rcp85m4p1q1_DJF','NEX_rcp85m4p1q3_JJA','NEX_rcp85m4p1q2_MAM','NEX_rcp85m4p1q4_SON',
                          'NEX_rcp85m4p0q1_DJF','NEX_rcp85m4p0q3_JJA','NEX_rcp85m4p0q2_MAM','NEX_rcp85m4p0q4_SON',
                          'NEX_rcp85m5p1q1_DJF','NEX_rcp85m5p1q3_JJA','NEX_rcp85m5p1q2_MAM','NEX_rcp85m5p1q4_SON',
                          'NEX_rcp85m5p0q1_DJF','NEX_rcp85m5p0q3_JJA','NEX_rcp85m5p0q2_MAM','NEX_rcp85m5p0q4_SON',
                          'NEX_rcp85m6p1q1_DJF','NEX_rcp85m6p1q3_JJA','NEX_rcp85m6p1q2_MAM','NEX_rcp85m6p1q4_SON',
                          'NEX_rcp85m6p0q1_DJF','NEX_rcp85m6p0q3_JJA','NEX_rcp85m6p0q2_MAM','NEX_rcp85m6p0q4_SON' )


pr.djf <- subset(prQuarter.nex,grep('DJF',names(prQuarter.nex),value=T)) - 
  subset(prQuarter.era,grep('DJF',names(prQuarter.era),value=T))
pr.jja <- subset(prQuarter.nex,grep('JJA',names(prQuarter.nex),value=T)) - 
  subset(prQuarter.era,grep('JJA',names(prQuarter.era),value=T))
pr.mam <- subset(prQuarter.nex,grep('MAM',names(prQuarter.nex),value=T)) - 
  subset(prQuarter.era,grep('MAM',names(prQuarter.era),value=T))
pr.son <- subset(prQuarter.nex,grep('SON',names(prQuarter.nex),value=T)) - 
  subset(prQuarter.era,grep('SON',names(prQuarter.era),value=T))

names(pr.djf) <- names(subset(prQuarter.nex,grep('DJF',names(prQuarter.nex),value=T)))
names(pr.jja) <- names(subset(prQuarter.nex,grep('JJA',names(prQuarter.nex),value=T)))
names(pr.mam) <- names(subset(prQuarter.nex,grep('MAM',names(prQuarter.nex),value=T)))
names(pr.son) <- names(subset(prQuarter.nex,grep('SON',names(prQuarter.nex),value=T)))
prQuarter.change <- stack(pr.djf,pr.jja,pr.mam,pr.son)
prQuarter.change <- subset(prQuarter.change,order(names(prQuarter.change)))

dir.create('data/climate-mask/pr-quarter-change');
writeRaster(prQuarter.change,'data/climate-mask/pr-quarter-change/',bylayer=T,names(prQuarter.change),format='GTiff',overwrite=T)


## select models for plot
prQuarter.change.select <- prQuarter.change[[grep('m5',names(prQuarter.change),value=T,invert=T)]]

modelNames <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR', 'ACCESS-1.0','CanESM2')

# setup color and break points
prbreaks <- seq(-2500, 2500, by = 500)
prcol <- colorRampPalette(c("darkred",'red','orange','yellow','lightblue','skyblue','blue','darkblue'))(length(prbreaks)-1)

levelplot(subset(prQuarter.change.select,grep('rcp45',names(prQuarter.change.select), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = prbreaks,col.regions=prcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(subset(prQuarter.change.select,grep('rcp85',names(prQuarter.change.select), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = prbreaks,col.regions=prcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

#-------------------
### annual tasmax
tasmaxQuarter.era <- stack(paste0('data/ERA-quarter-stats/',
                              list.files('data/ERA-quarter-stats/',pattern = 'ERA-tasmax-*')))
tasmaxQuarter.nex <- stack(paste0('data/NEX-quarter-stats/',
                              list.files('data/NEX-quarter-stats/',pattern = '.*-tasmax-*')))

tasmaxQuarter.era <- mask(crop(tasmaxQuarter.era,rangeDetail),rangeDetail)
tasmaxQuarter.nex <- mask(crop(tasmaxQuarter.nex,rangeDetail),rangeDetail)

names(tasmaxQuarter.era) <- c('ERA5_DJF','ERA5_JJA','ERA5_MAM','ERA5_SON')
names(tasmaxQuarter.nex) <- c('NEX_rcp45m1p1q1_DJF','NEX_rcp45m1p1q3_JJA','NEX_rcp45m1p1q2_MAM','NEX_rcp45m1p1q4_SON',
                          'NEX_rcp45m1p0q1_DJF','NEX_rcp45m1p0q3_JJA','NEX_rcp45m1p0q2_MAM','NEX_rcp45m1p0q4_SON',
                          'NEX_rcp45m2p1q1_DJF','NEX_rcp45m2p1q3_JJA','NEX_rcp45m2p1q2_MAM','NEX_rcp45m2p1q4_SON',
                          'NEX_rcp45m2p0q1_DJF','NEX_rcp45m2p0q3_JJA','NEX_rcp45m2p0q2_MAM','NEX_rcp45m2p0q4_SON',
                          'NEX_rcp45m3p1q1_DJF','NEX_rcp45m3p1q3_JJA','NEX_rcp45m3p1q2_MAM','NEX_rcp45m3p1q4_SON',
                          'NEX_rcp45m3p0q1_DJF','NEX_rcp45m3p0q3_JJA','NEX_rcp45m3p0q2_MAM','NEX_rcp45m3p0q4_SON',
                          'NEX_rcp45m4p1q1_DJF','NEX_rcp45m4p1q3_JJA','NEX_rcp45m4p1q2_MAM','NEX_rcp45m4p1q4_SON',
                          'NEX_rcp45m4p0q1_DJF','NEX_rcp45m4p0q3_JJA','NEX_rcp45m4p0q2_MAM','NEX_rcp45m4p0q4_SON',
                          'NEX_rcp45m5p1q1_DJF','NEX_rcp45m5p1q3_JJA','NEX_rcp45m5p1q2_MAM','NEX_rcp45m5p1q4_SON',
                          'NEX_rcp45m5p0q1_DJF','NEX_rcp45m5p0q3_JJA','NEX_rcp45m5p0q2_MAM','NEX_rcp45m5p0q4_SON',
                          'NEX_rcp45m6p1q1_DJF','NEX_rcp45m6p1q3_JJA','NEX_rcp45m6p1q2_MAM','NEX_rcp45m6p1q4_SON',
                          'NEX_rcp45m6p0q1_DJF','NEX_rcp45m6p0q3_JJA','NEX_rcp45m6p0q2_MAM','NEX_rcp45m6p0q4_SON',                          
                          'NEX_rcp85m1p1q1_DJF','NEX_rcp85m1p1q3_JJA','NEX_rcp85m1p1q2_MAM','NEX_rcp85m1p1q4_SON',
                          'NEX_rcp85m1p0q1_DJF','NEX_rcp85m1p0q3_JJA','NEX_rcp85m1p0q2_MAM','NEX_rcp85m1p0q4_SON',
                          'NEX_rcp85m2p1q1_DJF','NEX_rcp85m2p1q3_JJA','NEX_rcp85m2p1q2_MAM','NEX_rcp85m2p1q4_SON',
                          'NEX_rcp85m2p0q1_DJF','NEX_rcp85m2p0q3_JJA','NEX_rcp85m2p0q2_MAM','NEX_rcp85m2p0q4_SON',
                          'NEX_rcp85m3p1q1_DJF','NEX_rcp85m3p1q3_JJA','NEX_rcp85m3p1q2_MAM','NEX_rcp85m3p1q4_SON',
                          'NEX_rcp85m3p0q1_DJF','NEX_rcp85m3p0q3_JJA','NEX_rcp85m3p0q2_MAM','NEX_rcp85m3p0q4_SON',
                          'NEX_rcp85m4p1q1_DJF','NEX_rcp85m4p1q3_JJA','NEX_rcp85m4p1q2_MAM','NEX_rcp85m4p1q4_SON',
                          'NEX_rcp85m4p0q1_DJF','NEX_rcp85m4p0q3_JJA','NEX_rcp85m4p0q2_MAM','NEX_rcp85m4p0q4_SON',
                          'NEX_rcp85m5p1q1_DJF','NEX_rcp85m5p1q3_JJA','NEX_rcp85m5p1q2_MAM','NEX_rcp85m5p1q4_SON',
                          'NEX_rcp85m5p0q1_DJF','NEX_rcp85m5p0q3_JJA','NEX_rcp85m5p0q2_MAM','NEX_rcp85m5p0q4_SON',
                          'NEX_rcp85m6p1q1_DJF','NEX_rcp85m6p1q3_JJA','NEX_rcp85m6p1q2_MAM','NEX_rcp85m6p1q4_SON',
                          'NEX_rcp85m6p0q1_DJF','NEX_rcp85m6p0q3_JJA','NEX_rcp85m6p0q2_MAM','NEX_rcp85m6p0q4_SON' )


tasmax.djf <- subset(tasmaxQuarter.nex,grep('DJF',names(tasmaxQuarter.nex),value=T)) - 
  subset(tasmaxQuarter.era,grep('DJF',names(tasmaxQuarter.era),value=T))
tasmax.jja <- subset(tasmaxQuarter.nex,grep('JJA',names(tasmaxQuarter.nex),value=T)) - 
  subset(tasmaxQuarter.era,grep('JJA',names(tasmaxQuarter.era),value=T))
tasmax.mam <- subset(tasmaxQuarter.nex,grep('MAM',names(tasmaxQuarter.nex),value=T)) - 
  subset(tasmaxQuarter.era,grep('MAM',names(tasmaxQuarter.era),value=T))
tasmax.son <- subset(tasmaxQuarter.nex,grep('SON',names(tasmaxQuarter.nex),value=T)) - 
  subset(tasmaxQuarter.era,grep('SON',names(tasmaxQuarter.era),value=T))

names(tasmax.djf) <- names(subset(tasmaxQuarter.nex,grep('DJF',names(tasmaxQuarter.nex),value=T)))
names(tasmax.jja) <- names(subset(tasmaxQuarter.nex,grep('JJA',names(tasmaxQuarter.nex),value=T)))
names(tasmax.mam) <- names(subset(tasmaxQuarter.nex,grep('MAM',names(tasmaxQuarter.nex),value=T)))
names(tasmax.son) <- names(subset(tasmaxQuarter.nex,grep('SON',names(tasmaxQuarter.nex),value=T)))
tasmaxQuarter.change <- stack(tasmax.djf,tasmax.jja,tasmax.mam,tasmax.son)
tasmaxQuarter.change <- subset(tasmaxQuarter.change,order(names(tasmaxQuarter.change)))

dir.create('data/climate-mask/tasmax-quarter-change');
writeRaster(tasmaxQuarter.change,'data/climate-mask/tasmax-quarter-change/',bylayer=T,names(tasmaxQuarter.change),format='GTiff',overwrite=T)


## select models for plot
tasmaxQuarter.change.select <- tasmaxQuarter.change[[grep('m5',names(tasmaxQuarter.change),value=T,invert=T)]]

modelNames <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR', 'ACCESS-1.0','CanESM2')

# setup color and break points
tasmaxbreaks <- seq(-5, 15, by = 2.5)
tasmaxcol <- colorRampPalette(c('blue',"lightblue", "yellow", "orange", "red", "darkred","#32174D"))(length(tasmaxbreaks)-1)

levelplot(subset(tasmaxQuarter.change.select,grep('rcp45',names(tasmaxQuarter.change.select), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = tasmaxbreaks,col.regions=tasmaxcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))
#names.attr=c(rep(c(modelNames[1],rep('',3)),2),rep(c(modelNames[2],rep('',3)),2),rep(c(modelNames[3],rep('',3)),2),
#rep(c(modelNames[4],rep('',3)),2),rep(c(modelNames[6],rep('',3)),2)),

levelplot(subset(tasmaxQuarter.change.select,grep('rcp85',names(tasmaxQuarter.change.select), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = tasmaxbreaks,col.regions=tasmaxcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

#------------------ THAI
tha <- readOGR('data/shapefiles/rsStudy_NAME_0_Thailand.shp')
prQuarter.change.tha <- mask(crop(prQuarter.change.select,tha),tha)
tasmaxQuarter.change.tha <- mask(crop(tasmaxQuarter.change.select,tha),tha)

# setup color and break points
prbreaks <- seq(-2500, 2500, by = 500)
prcol <- colorRampPalette(c("darkred",'red','orange','yellow','lightblue','skyblue','blue','darkblue'))(length(prbreaks)-1)

levelplot(subset(prQuarter.change.tha,grep('rcp45',names(prQuarter.change.tha), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = prbreaks,col.regions=prcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(subset(prQuarter.change.tha,grep('rcp85',names(prQuarter.change.tha), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = prbreaks,col.regions=prcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

# setup color and break points
tasmaxbreaks <- seq(-5, 15, by = 2.5)
tasmaxcol <- colorRampPalette(c('blue',"lightblue", "yellow", "orange", "red", "darkred","#32174D"))(length(tasmaxbreaks)-1)

levelplot(subset(tasmaxQuarter.change.tha,grep('rcp45',names(tasmaxQuarter.change.tha), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = tasmaxbreaks,col.regions=tasmaxcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(subset(tasmaxQuarter.change.tha,grep('rcp85',names(tasmaxQuarter.change.tha), value=T)),layout=c(8,5),
          names.attr=c(rep(c('DJF','MAM','JJA','SON'),2),rep('',32)),
          at = tasmaxbreaks,col.regions=tasmaxcol,
          colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))
