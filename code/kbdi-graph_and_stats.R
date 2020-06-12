#[ status: on-going]
#[ note: ]
#@ This code is for visualizing time-series of KBDI from historical to far future scenario, as well as
#@ calculating descriptive statistics of the result within each scenario
#@ A. time-series  
#@ 1. structure of csv file -> 165 pixel locations, each with daily values over 5-year period
#@ 2. the data of the first 5 years (1994-1999) was used as a calibration period and excluded from this analysis
#@ 2. set a boundary and filter only pixels within that areas
#@ 3. perform daily mean within the boundary
#@ 4. combine all into single files
#@ 5. plot time series
#@
#@ B. calculate descriptive statistics over 3 period: 1) 2000-2019, 2) 2025-2044 and 3) 2055-2074
#@ 1. average annual kbdi and 3-month period (DJF, MAM, JJA, SON)
#@ 2. annual drought day duration (consecutive days with KBDI >  600)
#@ 3. repeated anomaly drought years (KBDI > 20-year average)

#=====================================#

# loading library
## data handling ##
library(tidyverse)

## generating graph ##
library(ggplot2)
library(maptools)
library(rasterVis)
library(gridExtra)

# parallel processing
library(foreach)
library(doParallel)