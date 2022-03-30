
### CODE SECTION CONTENTS ###
### setup ###
### compile pollen and occurrence data for Fraxinus ###
### create spatial polygon encompassing study extent ###
### mask Dalton et al 2020 ice sheet layers onto land mass from Lorenz et al 2016 ###
### generate elevation raster for study region ###
### calculate "biotic" velocity of exposed land as null model for real biotic velocity ###
### create publication-ready figure of study region ###

#############
### setup ###
#############
rm(list=ls())

library(sp)
library(raster)
library(dismo)
library(rgeos)
library(neotoma)
library(omnibus)
library(scales)
library(BIEN)
# library(ncdf4)
library(spatialEco)

library(enmSdm) # Adam's custom library on GitHub: adamlilith/enmSdm

rasterOptions(format='GTiff', overwrite=TRUE)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM/')

# get rasters representing land from Lorenz et al
lorenzDir <- paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ')

years <- seq(21000, 0, by = -1000)
rastInterpFx <- 'linear'

if (exists('lorenz')) rm(lorenz)
for (year in years) {
  lorenz <- if (exists('lorenz')) {
    stack(lorenz, raster(paste0(lorenzDir, '/', year, 'BP/an_avg_TMAX.tif')))
  } else {
    raster(paste0(lorenzDir, '/', year, 'BP/an_avg_TMAX.tif'))
  }
}

lorenz <- lorenz * 0 + 1
lorenzYears <- seq(21000, 0, by = -1000)
names(lorenz) <- paste0('yr', lorenzYears, 'bp')

# get carbon and calendar years for each ice layer for Dalton et al ice sheet data
	# note that the calendar years are not exactly the same as listed on Table 1 in Dalton et al 2020

daltonYears <- read.csv('/Volumes/lj_mac_22/Dalton et al 2020 QSR Ice Layers/Dalton et al 2020 QSR Dates from Shapefile Names.csv')

### interpolate land rasters to time periods matching Dalton et al ice sheet representations
interpTo <- -1000 * daltonYears$calKiloYear
lorenzInterp <- interpolateRasters(lorenz, interpFrom = -1 * lorenzYears,
                                   interpTo = interpTo, type = rastInterpFx)

names(lorenzInterp) <- paste0('yr', 1000 * daltonYears$calKiloYear, 'bp')

### for each cell, assign a value from 0 to 1 indicating proportion covered by ice sheet

daltonDir <- '/Volumes/lj_mac_22/Dalton et al 2020 QSR Ice Layers/RDA Files/'

for (countDalton in 1:nrow(daltonYears)) {
  
  # # ice shapefile
  daltonCalYear <- daltonYears$calKiloYear[countDalton]
  load(paste0(daltonDir, 'daltonEtAl2020_', sprintf('%02.2f', daltonCalYear), '_kiloCalYBP.rda'))
  
  
  # extract, remember proportion of cell covered by ice
  iceOnRast <- raster::extract(lorenzInterp[[countDalton]], daltonIce, cellnumbers=TRUE, weights=TRUE, normalizeWeights=FALSE)
  
  # transfer values back to raster
  vals <- values(lorenzInterp[[countDalton]])
  vals <- vals * 0
  
  for (countIce in seq_along(iceOnRast)) {
    
    vals[iceOnRast[[countIce]][ , 'cell']] <- vals[iceOnRast[[countIce]][ , 'cell']] + iceOnRast[[countIce]][ , 'weight']
    
  }
  
  lorenzInterp[[countDalton]] <- setValues(lorenzInterp[[countDalton]], values=vals)
  
}

### add "year 0" to Lorenz terrestrial

year0 <- lorenz[[nlayers(lorenz)]] * 0
lorenzInterp <- stack(lorenzInterp, year0)

lorenzInterp <- stack(lorenzInterp$yr20500bp, lorenzInterp$yr20500bp, lorenzInterp$yr19300bp, 
                      lorenzInterp$yr18000bp,lorenzInterp$yr16800bp, lorenzInterp$yr16100bp, 
                      lorenzInterp$yr14900bp, lorenzInterp$yr14200bp, lorenzInterp$yr12800bp, 
                      lorenzInterp$yr12100bp, lorenzInterp$yr11000bp, lorenzInterp$yr10300bp, 
                      lorenzInterp$yr9000bp, lorenzInterp$yr8100bp, lorenzInterp$yr7300bp, 
                      lorenzInterp$yr6300bp, lorenzInterp$yr5710bp, lorenzInterp$yr4500bp,
                      lorenzInterp$yr3200bp, lorenzInterp$yr2000bp, lorenzInterp$yr910bp,
                      lorenzInterp$yr0bp)

writeRaster(lorenzInterp, 
            './regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')

