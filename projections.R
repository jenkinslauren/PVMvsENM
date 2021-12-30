rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)

library(raster)
library(rgdal)
library(enmSdm)
library(terra)
library(maxnet)
library(sf)

library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

setwd('/Volumes/LJ MacBook Backup/MOBOT/PVMvsENM')
load('./PCA_Beyer')

getClimRasts <- function(pc, climYear) {
  fileName <- paste0('./data_and_analyses/env_data/Beyer/envDataClipped_',
                     climYear, 'KYBP_pc', pc, '.tif')
  if (file.exists(fileName)) {
    envData <- brick(fileName)
    # rename raster layers to pc's
    names(envData) <- paste0('pca', 1:pc)
  } else {
    pcPrediction <- list()
    # label each env layer by variable and year
    for (i in 1:length(clim)) {
      names(clim[[i]]) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                            "bio17", "bio18", "bio19", "bio4", "bio5", "bio6", "bio7",
                            "bio8", "bio9", "cloudiness", "relative_humidity")
      pcPrediction[i] <- raster::predict(clim[[i]], pca, index = 1:pc)
      names(pcPrediction[[i]]) <- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
    }
    
    envDataPca <- stack(pcPrediction)
    
    # keep only rasters (first five rasters) for climate year
    envYr <- pcPrediction[[(climYear/1000) + 1]]
    # plot(envYr) 
    names(envYr) <- paste0('pca', 1:pc)
    
    # check projections = wgs84
    print("Ensure that the projection of these rasters is WGS84:")
    print(paste0("Projection of envYr = ", projection(envYr)))
    
    # define study region & extent, load if already defined
    if (file.exists('./data_and_analyses/study_region/Study Region.Rdata') & 
        file.exists('./data_and_analyses/study_region/Study Region Extent.Rdata')) {
      load('./data_and_analyses/study_region/Study Region.Rdata')
      load('./data_and_analyses/study_region/Study Region Extent.Rdata')
    } else {
      studyRegion <- rgdal::readOGR(dsn = './data_and_analyses/study_region', 'study_region')
      studyRegion <- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
      projection(studyRegion) <- getCRS("WGS84")
      studyExtent <- extent(studyRegion)
      save(studyRegion, 
           file='./data_and_analyses/study_region/Study Region.Rdata', 
           compress=TRUE)
      save(studyExtent, 
           file='./data_and_analyses/study_region/Study Region Extent.Rdata', 
           compress=TRUE)
    }
    
    # clip environmental PCAs to study extent for given species, visualize, and save:
    envDataClipped <- list()
    for (i in 1:nlayers(envYr)) {
      x <- envYr[[i]]
      x <- crop(x, studyExtent)
      x <- mask(x, studyRegion)
      projection(x) <- getCRS("WGS84")
      envDataClipped[[i]] <- x
      envData <- stack(envDataClipped)
      # plot(envData)
      writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
    }
  }  
  return(envData)
}

# load model
speciesAb_ <- 'Frax_Penn'
load(paste0('./Models/', speciesAb_, '_Maxent/Model.Rdata'))

climYears <- seq(21000, 0, by=-1000)
predictors <- c(paste0('pca', 1:5))

if (exists('preds')) rm(preds)
for (climYear in climYears) {
  climate <- getClimRasts(pc = 5, climYear = climYear)
  thisPred <- predict(climate, envModel)
  preds <- if (exists('preds')) {
    stack(preds, thisPred)
  } else {
    thisPred
  }
}

# interpFrom <- -1 * climYears
# interpTo <- seq(-21000, 0, by=30)
# preds <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')

# # project
# preds <- projectRaster(preds, studyRegionRasts)
# preds <- calc(preds, fun=function(x) ifelse(x < 0, 0, x))
# preds <- calc(preds, fun=function(x) ifelse(x > 1, 1, x))

# # mask by study region and force values to be within [0, 1] (can get pushed outside this during re-projection)
# for (i in 1:nlayers(preds)) {

# landMask <- (1 - studyRegionRasts[[i]])
# preds[[i]] <- preds[[i]] * landMask

# }

# names(preds) <- paste0('ybp', seq(21000, 0, by=-30))
# writeRaster(preds, paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo))



