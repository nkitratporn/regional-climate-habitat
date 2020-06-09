#=====================================#
# loading library
## download packages not yet in CRAN (e.g. github)
library(devtools)

## raster and shapefile processing ##
library(rgeos)
library(rgdal)
library(raster)

## data handling ##
library(tidyverse)

## generating graph ##
library(rasterVis)
library(gridExtra)
library(ggpubr)

## rgee package allowing R to talk with EE
#remotes::install_github("r-spatial/rgee")
library(rgee)

## setup GDrive
library(googledrive) #install.packages("googledrive")

#=====================================#
# setting up environment

## specify the Python version in conda if it's setup manually. Otherwise,
## run ee_install() to automatically create an isolated Python virtual environment with all rgee dependencies 
## use_condaenv("ee", conda = "auto",required = TRUE) # existing environment
ee_install()

# connect and authenticate GEE and Gdrive
ee_Initialize(email = 'tita@g.ecc.u-tokyo.ac.jp') # even after authenticated successfully, Gdrive is not connected

# checking the setup
ee_check() # somehow always get stuck when authenticating Gdrive, but seems to work if directly use googledrive::drive_auth()
           # install missing packages ('oauth2client') manually using reticulate::conda_install() 

#=====================================#
# pre-processing climatic data from ERA5 and NEX-GDDP
#=====================================#

## Initialize common variables
### range extent
range <- ee$Geometry$Polygon(
  list(
    c(68.06467647557102, 30.625434988022814),
    c(68.06467647557102, -9.853266214563462),
    c(119.04123897557102, -9.853266214563462),
    c(119.04123897557102, 30.625434988022814)
  )
)
rangeExtent <- ee$Geometry(range,{},FALSE)        #planner geometry
rangeExtent <- ee$FeatureCollection(rangeExtent)  #as feature collection

Map$addLayer(rangeExtent,name='Range')            # check extent of range

## functions
### rename ERA5 bands
eraRename = function(img){
  newName <- list('pr','tasmax')
  bands <- list('total_precipitation','maximum_2m_air_temperature')
  return(img$select(bands)$rename(newName))
}

### convert ERA5 values from meter (m) to inch (in) / Kelvin (K) to Fehrenhin (F)
eraUnitConversion <- function(img){
  orig <- img;
  pr <- img$select('pr')$multiply(39.37)
  tasmax <- img$select('tasmax')$subtract(273.15)
  tasmax <- tasmax$multiply(9)$divide(5)
  tasmax <- tasmax$add(32)
  new <- pr$float()$addBands(tasmax$float())$copyProperties(orig,orig$propertyNames())
  
  return(new)
}

### convert NEX-GDDP values from flux of precipiation -> kg/(m^2*s) to in. / K to F and 
nexUnitConversion <- function(img){
    orig <- img;
    pr <- img$select('pr')$multiply(86400)$multiply(39.37)$rename('pr');
    tasmax <- img$select('tasmax')$subtract(273.15)$multiply(9)$divide(5)$add(32)$rename('tasmax');
    new <- pr$float()$addBands(tasmax$float())$copyProperties(orig,orig$propertyNames())
    return(new)
  }

### clip image in collection using rangeExtent
clipImage <- function(img){
  orig = img;
  return(img$clip(rangeExtent)$copyProperties(orig,orig$propertyNames()))
}

#=====================================#
## date for filtering
start <- c('1994-07-20')
end <- c('2019-12-31')

## load ERA5 image collection and filter
era <- ee$ImageCollection('ECMWF/ERA5/DAILY')$
  filterDate(start,end)$
  map(clipImage)$
  map(eraRename)$
  map(eraUnitConversion)

## study period (historical)
historyPeriod <- seq(1995,2019,5)

## download parameter
downConfig = list(region= rangeExtent$geometry(), scale = 28000, maxPixels = 1.0E13, driveFolder = 'ERA-stack')

## loop through each year, create year stack, download to drive (1 stack/yr)
for (y in historyPeriod) {
  pr <- era$select('pr')$filter(ee$Filter$calendarRange(y,y+4,'year'))$toBands()
  tasmax <- era$select('tasmax')$filter(ee$Filter$calendarRange(y,y+4,'year'))$toBands()
  exportpr <- paste0('pr-',y)
  exporttasmax <- paste0('tasmax-',y)
  task1 <- ee$batch$Export$image(pr, exportpr, downConfig)
  task2 <- ee$batch$Export$image(tasmax,exporttasmax, downConfig)
  task1$start()
  task2$start()
  cat('execute: ',y)
}
Map$addLayer(era$first()$select('pr'),list(min=0, max=6),'pr first day')
#=====================================#
## study period (future projection)
futurePeriod <- seq(2020,2074,5)

## scenario and models
scenario <- list('rcp45','rcp85')
models <- list('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR','CanESM2')

## download parameter
downConfig = list(region= rangeExtent$geometry(), scale = 28000, maxPixels = 1.0E13, driveFolder = 'NEX-stack')

for (y in futurePeriod) {
  for (m in 1:length(models)) {
    for (s in 1:length(scenario)) {
      col <- ee$ImageCollection("NASA/NEX-GDDP")$
        filter(ee$Filter$calendarRange(y,y+4,'year'))$
        filterMetadata('scenario','equals',scenario[[s]])$
        filterMetadata('model','equals',models[[m]])$
        map(clipImage)$
        map(nexUnitConversion)
      pr <- col$select('pr')$toBands()
      tasmax <- col$select('tasmax')$toBands()
      
      exportpr <- paste0('pr-',scenario[[s]],'-m',m,'-',y,'-',y+4)
      exporttasmax <- paste0('tasmax-',scenario[[s]],'-m',m,'-',y,'-',y+4)
      task1 <- ee$batch$Export$image(pr, exportpr, downConfig)
      task2 <- ee$batch$Export$image(tasmax,exporttasmax,downConfig)
      task1$start()
      task2$start()
      cat('execute: ',scenario[[s]],'-',models[[m]],' ')
    }
    cat('\n')
  }
  cat('Completed: ',y,'-',y+4,'\n')
}
