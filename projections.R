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

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

# constants
gcm <- 'Lorenz_ecbilt'
pc <- 5

load(paste0('./PCA_', gcm, '_PC', pc))

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda','Fraxinus quadrangulata')

# set constants
climYears <- seq(21000, 0, by=-1000)

studyRegionFileName <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/0ByjNJEf91IW5SUlEOUJFVGN0Y28/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

getClimRasts <- function(pc, climYear) {
  
  if (gcm == 'Beyer') {
    fileName <- paste0('./data_and_analyses/env_data/Beyer/envDataClipped_',
                       climYear, 'KYBP_pc', pc, '.tif')
    vars <- c("BIO1", paste0('BIO', 4:19), "cloudiness", "relative_humidity")
  } else if (gcm == 'Lorenz_ccsm') {
    fileName <- paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/envDataClipped_',
                       climYear, 'KYBP_pc', pc, '.tif')
    clim <- lorenz
  } else {
    fileName <- paste0('./data_and_analyses/env_data/Lorenz/V2/ecbilt_21-0k_all_tifs_LJ/envDataClipped_',
                       climYear, 'KYBP_pc', pc, '.tif')
    clim <- lorenz
  }
  
  if (file.exists(fileName)) {
    envData <- brick(fileName)
    # rename raster layers to pc's
    names(envData) <- paste0('pca', 1:pc)
  } else {
    
    envDataPca <- stack(pcPrediction)
    
    # keep only rasters (first five rasters) for climate year
    # climYear <- as.numeric(sub("KYBP.*", "", sub(".*_", "", names(envDataPca)[k*5])))
    #print(paste0("clim year = ", climYear))
    envYr <- pcPrediction[[(climYear/1000) + 1]]
    # plot(envYr) 
    names(envYr) <- paste0('pca', 1:pc)
    
    # check projections = wgs84
    # print("Ensure that the projection of these rasters is WGS84:")
    # print(paste0("Projection of envYr = ", projection(envYr)))
    # 
    # define study region & extent, load if already defined
    if (file.exists('./data_and_analyses/study_region/Study Region.Rdata') & 
        file.exists('./data_and_analyses/study_region/Study Region Extent.Rdata')) {
      load('./data_and_analyses/study_region/Study Region.Rdata')
      load('./data_and_analyses/study_region/Study Region Extent.Rdata')
    } else {
      studyRegion <- rgdal::readOGR('/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region',
                                    'study_region')
      studyRegion <- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
      projection(studyRegion) <- getCRS("WGS84")
      studyExtent <- extent(studyRegion)
      save(studyRegion, 
           file='/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region/Study Region.Rdata',
           compress=TRUE)
      save(studyExtent, 
           file='/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region/Study Region Extent.Rdata', compress=TRUE)
    }
    
    # clip environmental PCAs to study extent for given species, visualize, and save:
    envDataClipped <- list()
    for (i in 1:nlayers(envYr)) {
      x <- envYr[[i]]
      x <- crop(x, studyExtent)
      # x <<- mask(x, studyRegion)
      projection(x) <- getCRS("WGS84")
      envDataClipped[[i]] <- x
      envData <- stack(envDataClipped)
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
    # plot(thisPred, main = paste0("prediction at ", climYear, " ybp"))
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
  
  fileName <- list.files(path = paste0('./Models/Maxent/', speciesAb_, '_Maxent'),
                         pattern = paste0(gcm, '.Rdata'),
                         full.names = T)
  load(fileName)
  
  if(exists('preds')) rm(preds)
  preds <- predictions(speciesAb_, pc)
  
  # interpFrom <- -1 * climYears
  # interpTo <- seq(-21000, 0, by=1000)
  # predsInterp <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')
  
  # project
  preds <- projectRaster(preds, studyRegionRasts)
  
  # mask by study region and force values to be within [0, 1] 
  # (can get pushed outside this during re-projection)
  preds <- raster::calc(preds, fun = function(x) ifelse(x < 0, 0, x))
  preds <- raster::calc(preds, fun = function(x) ifelse(x > 1, 1, x))
  for (i in 1:nlayers(preds)) {
    landMask <- (1 - studyRegionRasts[[i]])
    preds[[i]] <- preds[[i]] * landMask
  }
  
  # names(preds) <- paste0('ybp', seq(21000, 0, by=-1000))
  
  writeRaster(stack(preds), paste0('./predictions/', speciesAb_, '_GCM_', gcm, '_PC', pc),
              format = 'GTiff', overwrite = T)
  
  # save.image(paste0('./workspaces/06 - Projections ENM ', speciesAb_))
}

fileName <- list.files(path = './predictions/',
                       pattern = paste0(gcm, '_PC', pc,'.tif'),
                       full.names = T)
tmp <- list(brick(fileName[[1]]), brick(fileName[[2]]), brick(fileName[[3]]), 
            brick(fileName[[4]]), brick(fileName[[5]]), brick(fileName[[6]]), 
            brick(fileName[[7]]), brick(fileName[[8]]))

meansList <- list()
maxList <- list()
sumList <- list()

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
  
  mx <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]], 
                       n[[7]], n[[8]], fun = sum)
  names(mx) <- paste0(climYear, ' ybp')
  sumList <- append(sumList, mx)
  
}

pdf(file = paste0('./predictions/pdf/', gcm, '_meansList.pdf'), width = 11, height = 8.5)
for (i in 1:length(meansList)) {
  plot(meansList[[i]], main = names(meansList[[i]]))
}
# plot(stack(meansList))
dev.off()

pdf(file = paste0('./predictions/pdf/', gcm, '_maxList.pdf'), width = 11, height = 8.5)
for (i in 1:length(maxList)) {
  plot(maxList[[i]], main = names(maxList[[i]]))
}
# plot(stack(maxList))
dev.off()

pdf(file = paste0('./predictions/pdf/', gcm, '_sumList.pdf'), width = 11, height = 8.5)
for (i in 1:length(sumList)) {
  plot(sumList[[i]], main = names(sumList[[i]]))
}
dev.off()


