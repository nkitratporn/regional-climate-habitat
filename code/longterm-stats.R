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
### annual pr and tasmax 
prMean.era <- stack(paste0('data/ERA5-longterm-stats/',
                           list.files('data/ERA5-longterm-stats/',pattern = '*pr-annual-mean.tif$')))
prMean.nex <- stack(paste0('data/NEX-longterm-stats/',
                           list.files('data/NEX-longterm-stats/',pattern = '*pr-annual-mean.tif$')))
prStd.era <- stack(paste0('data/ERA5-longterm-stats/',
                          list.files('data/ERA5-longterm-stats/',pattern = '*pr-annual-std.tif$')))
prStd.nex <- stack(paste0('data/NEX-longterm-stats/',
                          list.files('data/NEX-longterm-stats/',pattern = '*pr-annual-std.tif$')))

tasmaxMean.era <- stack(paste0('data/ERA5-longterm-stats/',
                               list.files('data/ERA5-longterm-stats/',pattern = '*tasmax-annual-mean.tif$')))
tasmaxMean.nex <- stack(paste0('data/NEX-longterm-stats/',
                               list.files('data/NEX-longterm-stats/',pattern = '*tasmax-annual-mean.tif$')))
tasmaxStd.era <- stack(paste0('data/ERA5-longterm-stats/',
                              list.files('data/ERA5-longterm-stats/',pattern = '*tasmax-annual-std.tif$')))
tasmaxStd.nex <- stack(paste0('data/NEX-longterm-stats/',
                              list.files('data/NEX-longterm-stats/',pattern = '*tasmax-annual-std.tif$')))

prMean <- stack(prMean.era,prMean.nex)
tasmaxMean <- stack(tasmaxMean.era,tasmaxMean.nex)
prStd <- stack(prStd.era,prStd.nex)
tasmaxStd <- stack(tasmaxStd.era,tasmaxStd.nex)

layernames <- c('ERA5',           'NEX_rcp45m1p1' , 'NEX_rcp45m1p0',
                'NEX_rcp45m2p1' , 'NEX_rcp45m2p0' , 'NEX_rcp45m3p1',
                'NEX_rcp45m3p0' , 'NEX_rcp45m4p1' , 'NEX_rcp45m4p0', 
                'NEX_rcp45m5p1' , 'NEX_rcp45m5p0' , 'NEX_rcp45m6p1',
                'NEX_rcp45m6p0' , 'NEX_rcp85m1p1' , 'NEX_rcp85m1p0',
                'NEX_rcp85m2p1' , 'NEX_rcp85m2p0' , 'NEX_rcp85m3p1',
                'NEX_rcp85m3p0' , 'NEX_rcp85m4p1' , 'NEX_rcp85m4p0',
                'NEX_rcp85m5p1' , 'NEX_rcp85m5p0' , 'NEX_rcp85m6p1', 'NEX_rcp85m6p0')

names(prMean) <- layernames
names(tasmaxMean) <- layernames
names(prStd) <- layernames
names(tasmaxStd) <- layernames

prMean <- subset(prMean,order(names(prMean)))
tasmaxMean <- subset(tasmaxMean,order(names(tasmaxMean)))
prStd <- subset(prStd,order(names(prStd)))
tasmaxStd <- subset(tasmaxStd,order(names(tasmaxStd)))

#-------------------------------------#
# mask
prMean <- mask(prMean,rangeDetail)
tasmaxMean <- mask(tasmaxMean,rangeDetail)
prStd <- mask(prStd,rangeDetail)
tasmaxStd <- mask(tasmaxStd,rangeDetail)

writeRaster(prMean,'data/climate-mask/prMean/',bylayer=T,names(prMean),format='GTiff',overwrite=T)
writeRaster(prStd,'data/climate-mask/prStd/',bylayer=T,names(prStd),format='GTiff',overwrite=T)
writeRaster(tasmaxMean,'data/climate-mask/tasmaxMean/',bylayer=T,names(tasmaxMean),format='GTiff',overwrite=T)
writeRaster(tasmaxStd,'data/climate-mask/tasmaxStd/',bylayer=T,names(tasmaxStd),format='GTiff',overwrite=T)
#-------------------------------------#

regionalAnnual <-list(prMean,prStd,tasmaxMean,tasmaxStd)
annual.df <- data.frame()
var.name <- var.name <- c('pr_mean','pr_stddev','tasmax_mean','tasmax_stddev')

## extract regional-scale mean and stddev into df
annual.df <- foreach(i=1:length(regionalAnnual), .packages = 'raster',
                     .combine = 'rbind',.export = 'annual.df') %dopar% {
                       temp <- raster::extract(regionalAnnual[[i]],rangeDetail,fun='mean',na.rm=T,df=T)[1,]
                       temp[1,1] <- var.name[i]
                       temp
                     }

annual.df

## extract by states mean and stddev into df
country.df <- data.frame()
country.df <- foreach(i=1:length(regionalAnnual), .packages = c('raster','dplyr','reshape2'), 
                      .combine = 'cbind', .export = 'country.df') %dopar%{
                        temp <- raster::extract(regionalAnnual[[i]],rangeCountry,fun='mean',na.rm=T,df=T,factors=T)
                        temp$ID <- rangeCountry$NAME_0
                        temp <- melt(temp,id.vars='ID')
                        names(temp)[names(temp) == 'ID'] <- 'Country'
                        names(temp)[names(temp) == 'value'] <- var.name[i]
                        temp
                      }
head(country.df,15)



regionAnnual.df <- data.frame(t(annual.df[-1])) #transpose 
colnames(regionAnnual.df) <- annual.df[,1]      #rename column
regionAnnual.df$Country <- 'Regional'
regionAnnual.df$variable <- rownames(regionAnnual.df)
rownames(regionAnnual.df) <- NULL
head(regionAnnual.df)

country.df <- subset(country.df,select= -c(4,5,7,8,10,11))
names(country.df)[names(country.df) == 'ID'] <- 'Country'
country.df <- country.df[,c(3,4,5,6,1,2)]
head(country.df,15)

annual.df <- rbind(regionAnnual.df,country.df)
str(annual.df)
summary(annual.df)
head(annual.df,25)

rm(regionAnnual.df,country.df)

#--------
## Change
### Mean
prChange <- prMean - prMean[[1]]
tasmaxChange <- tasmaxMean - tasmaxMean[[1]]
names(prChange) <- names(prMean)
names(tasmaxChange) <- names(tasmaxMean)
dir.create('data/climate-mask/pr-change'); dir.create('data/climate-mask/tasmax-change')
writeRaster(prChange,'data/climate-mask/pr-change/',bylayer=T,names(prChange),format='GTiff',overwrite=T)
writeRaster(tasmaxChange,'data/climate-mask/tasmax-change/',bylayer=T,names(prStd),format='GTiff',overwrite=T)

regionalAnnual <-list(prChange,tasmaxChange)
change.df <- data.frame()
var.name <- var.name <- c('pr_change','tasmax_change')

## extract regional-scale mean and stddev into df
change.df <- foreach(i=1:length(regionalAnnual), .packages = 'raster',
                     .combine = 'rbind',.export = 'change.df') %dopar% {
                       temp <- raster::extract(regionalAnnual[[i]],rangeDetail,fun='mean',na.rm=T,df=T)[1,]
                       temp[1,1] <- var.name[i]
                       temp
                     }

change.df

change.c.df <- data.frame()
change.c.df <- foreach(i=1:length(regionalAnnual), .packages = c('raster','dplyr','reshape2'), 
                      .combine = 'cbind', .export = 'change.c.df') %dopar%{
                        temp <- raster::extract(regionalAnnual[[i]],rangeCountry,fun='mean',na.rm=T,df=T,factors=T)
                        temp$ID <- rangeCountry$NAME_0
                        temp <- melt(temp,id.vars='ID')
                        names(temp)[names(temp) == 'ID'] <- 'Country'
                        names(temp)[names(temp) == 'value'] <- var.name[i]
                        temp
                      }
head(change.c.df,15)

regionAnnual.df <- data.frame(t(change.df[-1])) #transpose 
colnames(regionAnnual.df) <- change.df[,1]    #rename column
regionAnnual.df$Country <- 'Regional'
regionAnnual.df$variable <- rownames(regionAnnual.df)
rownames(regionAnnual.df) <- NULL
head(regionAnnual.df)

change.c.df <- subset(change.c.df,select= -c(4,5))
names(change.c.df)[names(change.c.df) == 'ID'] <- 'Country'
change.c.df <- change.c.df[,c(3,4,1,2)]
head(change.c.df,15)

change.df <- rbind(regionAnnual.df,change.c.df)
head(change.df,25)

rm(change.c.df,regionAnnual.df)

### Std
prStdChange <- prStd - prStd[[1]]
tasmaxStdChange <- tasmaxStd - tasmaxStd[[1]]
names(prStdChange) <- names(prStd)
names(tasmaxStdChange) <- names(tasmaxStd)
dir.create('data/climate-mask/prStddev-change'); dir.create('data/climate-mask/tasmaxStddev-change')
writeRaster(prStdChange,'data/climate-mask/prStddev-change/',bylayer=T,names(prStdChange),format='GTiff',overwrite=T)
writeRaster(tasmaxStdChange,'data/climate-mask/tasmaxStddev-change/',bylayer=T,names(tasmaxStdChange),format='GTiff',overwrite=T)


regionalAnnual <-list(prStdChange,tasmaxStdChange)
stdChange.df <- data.frame()
var.name <- var.name <- c('prStddev_change','tasmaxStddev_change')

## extract regional-scale mean and stddev into df
stdChange.df <- foreach(i=1:length(regionalAnnual), .packages = 'raster',
                     .combine = 'rbind',.export = 'stdChange.df') %dopar% {
                       temp <- raster::extract(regionalAnnual[[i]],rangeDetail,fun='mean',na.rm=T,df=T)[1,]
                       temp[1,1] <- var.name[i]
                       temp
                     }

stdChange.df

stdChange.c.df <- data.frame()
stdChange.c.df <- foreach(i=1:length(regionalAnnual), .packages = c('raster','dplyr','reshape2'), 
                       .combine = 'cbind', .export = 'stdChange.c.df') %dopar%{
                         temp <- raster::extract(regionalAnnual[[i]],rangeCountry,fun='mean',na.rm=T,df=T,factors=T)
                         temp$ID <- rangeCountry$NAME_0
                         temp <- melt(temp,id.vars='ID')
                         names(temp)[names(temp) == 'ID'] <- 'Country'
                         names(temp)[names(temp) == 'value'] <- var.name[i]
                         temp
                       }
head(stdChange.c.df,15)

regionAnnual.df <- data.frame(t(stdChange.df[-1])) #transpose 
colnames(regionAnnual.df) <- stdChange.df[,1]    #rename column
regionAnnual.df$Country <- 'Regional'
regionAnnual.df$variable <- rownames(regionAnnual.df)
rownames(regionAnnual.df) <- NULL
head(regionAnnual.df)

stdChange.c.df <- subset(stdChange.c.df,select= -c(4,5))
names(stdChange.c.df)[names(stdChange.c.df) == 'ID'] <- 'Country'
stdChange.c.df <- stdChange.c.df[,c(3,4,1,2)]
head(stdChange.c.df,15)

stdChange.df <- rbind(regionAnnual.df,stdChange.c.df)
head(stdChange.df,30)

rm(stdChange.c.df,regionAnnual.df)

#---
# combine
longterm.stats <- cbind(annual.df,change.df,stdChange.df)
head(longterm.stats)
longterm.stats[,6] == longterm.stats[,10]
longterm.stats[,6] == longterm.stats[,14]
longterm.stats[,5] == longterm.stats[,9];
longterm.stats[,5] == longterm.stats[,13]
longterm.stats <- subset(longterm.stats,select= -c(9,10,13,14))
head(longterm.stats)

rm(annual.df,change.df,stdChange.df)

#----
# add relavent column
longterm.stats$Scenario <- ifelse(grepl('rcp45',longterm.stats$variable),'RCP45',
                             ifelse(grepl('rcp85',longterm.stats$variable),'RCP85','History'))
longterm.stats$Model <- ifelse(grepl('m1',longterm.stats$variable),'CESM1-BGC',
                          ifelse(grepl('m2',longterm.stats$variable),'MPI-ESM-MR',
                                 ifelse(grepl('m3',longterm.stats$variable),'MIROC5', 
                                        ifelse(grepl('m4',longterm.stats$variable),'IPSL-CM5A-MR',
                                               ifelse(grepl('m5',longterm.stats$variable),'ACCESS1-0',
                                                      ifelse(grepl('m6',longterm.stats$variable),'CanESM2','ERA5'))))))
longterm.stats$Period <- ifelse(grepl('p0',longterm.stats$variable),'Near Future',
                           ifelse(grepl('p1',longterm.stats$variable),'Far Future','History'))
longterm.stats$Model <- as.factor(longterm.stats$Model)
longterm.stats$Scenario <- as.factor(longterm.stats$Scenario)
longterm.stats$Period <- as.factor(longterm.stats$Period)
longterm.stats$Country <- as.factor(longterm.stats$Country)
summary(longterm.stats)

longterm.stats$Period <- factor(longterm.stats$Period, levels=c('History','Near Future','Far Future'))
longterm.stats$Country <- factor(longterm.stats$Country,levels = c('Regional','Bangladesh','Bhutan','Cambodia','China','India',
                                                                   'Indonesia','Laos','Malaysia','Myanmar','Nepal','Sri Lanka',
                                                                   'Thailand','Vietnam'))


#save
write.csv(longterm.stats,'output/annual-longterm-stats-final.csv')


#--------
# dryDays
dryDays.era <- stack(paste0('data/ERA5-longterm-stats/',
                            list.files('data/ERA5-longterm-stats/',pattern = '*-drydays-*')))
dryDays.nex <- stack(paste0('data/NEX-longterm-stats/',
                            list.files('data/NEX-longterm-stats/',pattern = '*-drydays-*')))

stopCluster(cl)


#-------
# graph
modelNames <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR', 'ACCESS-1.0','CanESM2')

## REGIONAL SCALE ##
## change in pr and tasmax

#selected models
prChange.select <- dropLayer(prChange,c(10,11,22,23))
tasmaxChange.select <- dropLayer(tasmaxChange,c(10,11,22,23))

# setup color and break points
prbreaks <- seq(-3000, 3000, by = 500)
tasbreaks <- seq(-5,15,by=2.5)
tascol <- colorRampPalette(c('blue',"lightblue", "yellow", "orange", "red", "darkred","#32174D"))(length(tasbreaks)-1) #
prcol <- colorRampPalette(c("darkred",'red','orange','yellow','lightblue','skyblue','blue','darkblue'))(length(prbreaks)-1)

levelplot(prChange.select[[2:nlayers(prChange.select)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = prbreaks,
          col.regions=prcol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(tasmaxChange.select[[2:nlayers(tasmaxChange.select)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = tasbreaks,
          col.regions=tascol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

# change in pr and tasmax
#selected models
prStdChange.select <- dropLayer(prStdChange,c(10,11,22,23))
tasmaxStdChange.select <- dropLayer(tasmaxStdChange,c(10,11,22,23))

# setup color and break points
prbreaks <- seq(-1000, 1500, by = 250)
tasbreaks <- seq(-0.6,0.6,by=0.20)
tascol <- colorRampPalette(c('blue',"lightblue", "yellow", "orange", "red"))(length(tasbreaks)-1) #"darkred", "#32174D"
prcol <- colorRampPalette(c("blue",'skyblue', "yellow", "orange", "red",'darkred'))(length(prbreaks)-1)

levelplot(prStdChange.select[[2:nlayers(prStdChange.select)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = prbreaks,
          col.regions=prcol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(tasmaxStdChange.select[[2:nlayers(tasmaxStdChange.select)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = tasbreaks,
          col.regions=tascol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))


## THAILAND ##
tha <- readOGR('data/shapefiles/rsStudy_NAME_0_Thailand.shp')
prChange.tha <- mask(crop(prChange.select,tha),tha)
tasmaxChange.tha <- mask(crop(tasmaxChange.select,tha),tha)
prStdChange.tha <- mask(crop(prStdChange.select,tha),tha)
tasmaxStdChange.tha <- mask(crop(tasmaxStdChange.select,tha),tha)

# setup color and break points
prbreaks <- seq(-3000, 3000, by = 500)
tasbreaks <- seq(-5,15,by=2.5)
tascol <- colorRampPalette(c('blue',"lightblue", "yellow", "orange", "red", "darkred","#32174D"))(length(tasbreaks)-1) #
prcol <- colorRampPalette(c("darkred",'red','orange','yellow','lightblue','skyblue','blue','darkblue'))(length(prbreaks)-1)

levelplot(prChange.tha[[2:nlayers(prChange.tha)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = prbreaks,
          col.regions=prcol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(tasmaxChange.tha[[2:nlayers(tasmaxChange.tha)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = tasbreaks,
          col.regions=tascol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

# setup color and break points
prbreaks <- seq(-1000, 1500, by = 250)
tasbreaks <- seq(-0.6,0.6,by=0.20)
tascol <- colorRampPalette(c('blue',"lightblue", "yellow", "orange", "red"))(length(tasbreaks)-1) #"darkred", "#32174D"
prcol <- colorRampPalette(c("blue",'skyblue', "yellow", "orange", "red",'darkred'))(length(prbreaks)-1)

levelplot(prStdChange.tha[[2:nlayers(prStdChange.tha)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = prbreaks,
          col.regions=prcol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))

levelplot(tasmaxStdChange.tha[[2:nlayers(tasmaxStdChange.tha)]],
          layout = c(4,5), names.attr=c(rep("",20)),
          at = tasbreaks,
          col.regions=tascol, colorkey = list(space = "bottom", height = 1, width = 1),
          xlab = list(label = "", vjust = -.2), ylab = list(label=""))
