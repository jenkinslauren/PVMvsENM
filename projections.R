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

# constants
gcm <- 'ecbilt'
pc <- 5
workingFolder <- if(gcm = 'Beyer') {
  './data_and_analyses/env_data/Beyer/tifs'
} else {
  paste0('./data_and_analyses/env_data/Lorenz/V2/',
                          gcm, '_21-0k_all_tifs_LJ')
}

load('./PCA_Beyer')

speciesList <- c('Fraxinus americana', 'Fraxinus caroliniana', 'Fraxinus cuspidata'
                 ,'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# set constants
climYears <- seq(21000, 0, by=-1000)
pc <- 5

studyRegionFileName <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/0ByjNJEf91IW5SUlEOUJFVGN0Y28/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

getClimRasts <- function(pc, climYear) {
  fileName <<- paste0('./data_and_analyses/env_data/Beyer/envDataClipped_',
                     climYear, 'KYBP_pc', pc, '.tif')
  if (file.exists(fileName)) {
    envData <<- brick(fileName)
    # rename raster layers to pc's
    names(envData) <<- paste0('pca', 1:pc)
  } else {
    pcPrediction <<- list()
    # label each env layer by variable and year
    for (i in 1:length(clim)) {
      names(clim[[i]]) <<- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                            "bio17", "bio18", "bio19", "bio4", "bio5", "bio6", "bio7",
                            "bio8", "bio9", "cloudiness", "relative_humidity")
      pcPrediction[i] <<- raster::predict(clim[[i]], pca, index = 1:pc)
      names(pcPrediction[[i]]) <<- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
    }
    
    envDataPca <<- stack(pcPrediction)
    
    # keep only rasters (first five rasters) for climate year

    envYr <<- pcPrediction[[(climYear/1000) + 1]]
    # plot(envYr) 
    names(envYr) <<- paste0('pca', 1:pc)
    
    # check projections = wgs84
    print("Ensure that the projection of these rasters is WGS84:")
    print(paste0("Projection of envYr = ", projection(envYr)))
    
    # define study region & extent, load if already defined
    if (file.exists('./data_and_analyses/study_region/Study Region.Rdata') & 
        file.exists('./data_and_analyses/study_region/Study Region Extent.Rdata')) {
      load('./data_and_analyses/study_region/Study Region.Rdata')
      load('./data_and_analyses/study_region/Study Region Extent.Rdata')
    } else {
      studyRegion <<- rgdal::readOGR(dsn = './data_and_analyses/study_region', 'study_region')
      studyRegion <<- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
      projection(studyRegion) <<- getCRS("WGS84")
      studyExtent <<- extent(studyRegion)
      save(studyRegion, 
           file='./data_and_analyses/study_region/Study Region.Rdata', 
           compress=TRUE)
      save(studyExtent, 
           file='./data_and_analyses/study_region/Study Region Extent.Rdata', 
           compress=TRUE)
    }
    
    # clip environmental PCAs to study extent for given species, visualize, and save:
    envDataClipped <<- list()
    for (i in 1:nlayers(envYr)) {
      x <<- envYr[[i]]
      x <<- crop(x, studyExtent)
      # x <<- mask(x, studyRegion)
      projection(x) <<- getCRS("WGS84")
      envDataClipped[[i]] <<- x
      envData <<- stack(envDataClipped)
      # plot(envData)
      writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
    }
  }  
  return(envData)
}

predictions <- function(speciesAb_, pc) {
  
  predictors <<- c(paste0('pca', 1:pc))
  
  for (climYear in climYears) {
    climate <<- getClimRasts(pc = 5, climYear = climYear)
    # plot(climate)
    climate <<- climate[[predictors]]
    thisPred <<- predict(climate, envModel, clamp = F, type='cloglog')
    names(thisPred) <<- (paste0(climYear, ' ybp'))
    preds <<- if (exists('preds')) {
      stack(preds, thisPred)
    } else {
      thisPred
    }
  }
  return(preds)
}


for(sp in speciesList) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  
  fileName <- list.files(path = paste0('./Models/', speciesAb_, '_Maxent'),
                         pattern = 'Model.Rdata',
                         full.names = T)
  load(fileName)
  
  if(exists('preds')) rm(preds)
  preds <- predictions(speciesAb_, pc)
  
  # interpFrom <- -1 * climYears
  # interpTo <- seq(-21000, 0, by=1000)
  # predsInterp <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')
  
  # project
  preds <- projectRaster(preds, studyRegionRasts)
  
  # why isn't this working ? because the max value is -8! so everything just becomes 0
  preds <- raster::calc(preds, fun = function(x) ifelse(x < 0, 0, x))
  preds <- raster::calc(preds, fun = function(x) ifelse(x > 1, 1, x))
  
  r <- c(-Inf, 0, 0, 1, Inf, 1)
  r <- matrix(r, ncol = 3, byrow = TRUE)
  preds <- raster::reclassify(preds, r)
  
  # mask by study region and force values to be within [0, 1] 
  # (can get pushed outside this during re-projection)
  for (i in 1:nlayers(preds)) {
    landMask <- (1 - studyRegionRasts[[i]])
    preds[[i]] <- preds[[i]] * landMask
  }
  
  names(preds) <- paste0('ybp', seq(21000, 0, by=-1000))
  
  writeRaster(stack(preds), paste0('./predictions/Beyer_pca5_', speciesAb_),
              format = 'GTiff', overwrite = T)
  
  save.image(paste0('./workspaces/06 - Projections ENM ', speciesAb_))
}

fileName <- list.files(path = './predictions/',
                       pattern = '.tif',
                       full.names = T)
tmp <- list(brick(fileName[[1]]), brick(fileName[[2]]), brick(fileName[[3]]), 
             brick(fileName[[4]]), brick(fileName[[5]]), brick(fileName[[6]]), 
             brick(fileName[[7]]), brick(fileName[[8]]))
meansList <- list()
maxList <- list()
for (i in 1:length(climYears)) {
  climYear <- climYears[i]
  n <- stack(tmp[[1]][[i]], tmp[[2]][[i]], tmp[[3]][[i]], tmp[[4]][[i]], 
             tmp[[5]][[i]], tmp[[6]][[i]], tmp[[7]][[i]], tmp[[8]][[i]])
  mn <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]], 
                         n[[7]], n[[8]], fun = mean)
  names(mn) <- paste0(climYear, ' ybp')
  meansList <- append(meansList, mn)
  
  mx <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]], 
                        n[[7]], n[[8]], fun = max)
  names(mx) <- paste0(climYear, ' ybp')
  maxList <- append(maxList, mx)
  
}

plot(stack(meansList))
plot(stack(maxList))
  


