# SDM Template to run on cluster
# Author: Lauren Jenkins
# 27 April 2022

rm(list = ls())
# Load required packages.
# Load packages
library(enmSdm)

# for handling rasters
library(raster)
library(rgdal)
library(dismo)
# library(rgeos) # not used, but may be helpful when debugging

# library(terra)
# library(maxnet)

# visualization tools
library(sp)
library(sf)
library(maps)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# additional tools
library(tools)
library(units)
library(dplyr)
library(tidyr)

setwd('/mnt/research/TIMBER/PVMvsENM')

args <- commandArgs(TRUE)
gcm <- args[1]

# nad27
default_crs = sf::st_crs(4267)

# set constants
ll <- c('longitude', 'latitude')
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

pc <- 5
climYear <- 0

lorenzRast <- raster::raster('./in/env_data/Lorenz/ccsm/0BP/an_avg_ETR.tif')

for(sp in speciesList) {
  sp <- sp
  print(paste0("SPECIES = ", sp))
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  load(paste0('./in/little_range_map/', 
              rangeName, '.Rdata'))
  gcm <- gcm
  # print(paste0('GCM = ', gcm))
  
  # Load environmental PCA rasters. 
  # If already clipped, load that. 
  # Otherwise, clip the data to present.
  studyRegionFileName <- './in/study_region/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
  studyRegionRasts <- brick(studyRegionFileName)
  
  if (gcm == 'Beyer') { # Beyer
    load('./in/env_data/Beyer/PCA_clim.Rdata')
    load(paste0('./in/env_data/Beyer/pca_pc', pc, '.Rdata'))
    fileName <- paste0('./in/env_data/Beyer/envDataClipped_',
                       climYear, 'YBP_pc', pc, '.tif')
    vars <- c("BIO1", paste0('BIO', 4:19), "cloudiness", "relative_humidity")
  } else if (gcm == 'Lorenz_ccsm') { # CCSM
    load(paste0('./in/env_data/Lorenz/ccsm/PCA_', gcm, '_clim.Rdata')) 
    load(paste0('./in/env_data/Lorenz/ccsm/pca_pc', pc, '.Rdata')) 
    fileName <- paste0('./in/env_data/Lorenz/ccsm/envDataClipped_',
                       climYear, 'YBP_pc', pc, '.tif')
    clim <- lorenz
    workingFolder <- paste0('./in/env_data/Lorenz/ccsm/', 
                            climYear, 'BP')
    vars <- names(clim[[1]])
  } else { # ECBilt
    load(paste0('./in/env_data/Lorenz/ecbilt/PCA_', gcm, '_clim.Rdata')) 
    load(paste0('./in/env_data/Lorenz/ecbilt/pca_pc', pc, '.Rdata'))
    fileName <- paste0('./in/env_data/Lorenz/ecbilt/envDataClipped_',
                       climYear, 'YBP_pc', pc, '.tif')
    clim <- lorenz
    workingFolder <- paste0('./in/env_data/Lorenz/ecbilt/', climYear, 'BP')
    vars <- names(clim[[1]])
  }
  
  if (file.exists(fileName)) {
    envData <- brick(fileName)
    # rename raster layers to pc's
    names(envData) <- paste0('pca', 1:pc)
  } else {
    pcPrediction <<- list()
    # label each env layer by variable and year
    for (i in 1:length(clim)) {
      names(clim[[i]]) <- vars
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
    
    # clip environmental PCAs to study extent for given species, visualize, and save:
    envDataClipped <- list()
    for (n in 1:nlayers(envYr)) {
      x <- envYr[[n]]
      x <- crop(x, extent(studyRegionRasts[[yrs[(climYear/1000) + 1]]]))
      # x <<- mask(x, studyRegion)
      projection(x) <- getCRS("WGS84")
      envDataClipped[[n]] <- x
    }
    envData <- stack(envDataClipped)
    writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
  } 
  
  fileName <- paste0('./in/cleaned_records/', 
                     species, '_finalRecords.rData')
  load(fileName)
  records <- data.frame(speciesSf_filtered_final)
  records$geometry <- gsub("[c()]", "", records$geometry)
  records <- separate(data = records, 
                      col = 'geometry', 
                      into = ll, 
                      sep = "\\,")
  records$longitude <- as.double(records$longitude)
  records$latitude <- as.double(records$latitude)
  
  occsEnv <- raster::extract(envData, 
                             cbind(records$longitude, 
                                   records$latitude))
  occsEnvDf <- as.data.frame(occsEnv)
  records <- cbind(records, occsEnvDf)
  
  # records in the water:
  if (any(is.na(rowSums(occsEnvDf)))) {
    water <- records[which(is.na(rowSums(occsEnvDf))), ] 
    water <- SpatialPointsDataFrame(water[,ll], data = water,
                                    proj4 = getCRS('wgs84', TRUE))
  }
  
  # remove records in water from dataset:
  if (any(is.na(rowSums(occsEnvDf)))) records <-
    records[-which(is.na(rowSums(occsEnvDf))), ]
  
  # visualize the points that fall in the water
  # convert to sp object
  recordsSp <- SpatialPointsDataFrame(records[, ll], data = records,
                                      proj4 = getCRS('wgs84', TRUE))
  
  # visualize
  plot(recordsSp, pch = 16, cex = 0.5, col = "red", 
       main = paste0(sp, ' occurrences (BIEN) thinned'))
  if (exists("water")) {
    plot(water, col = 'blue', add = TRUE)
  }
  map("state", add = TRUE)
  map("world", add = TRUE)
  
  save.image(paste0('./in/workspaces/04 - Modeling Workspace - Clipping ', 
                    sp, '_PC_', pc, '_GCM_', gcm))
  # load(paste0('./workspaces/04 - Modeling Workspace - Clipping ', sp, '_PC_', pc, '_GCM_', gcm))
  
  bufferFileName <- paste0('./in/cleaned_records/', species, '_buffer.rData')
  load(bufferFileName)
  
  # calculate calibration region buffer at 160-km to extract bg sites
  calibBuffer <- st_buffer(st_transform(st_as_sf(x = recordsSp), getCRS('albersNA')), 
                           dist = as_units(160, 'km'))
  calibBuffer <- st_union(calibBuffer)
  
  calibRegionSpAlb <- sp::spTransform(as(calibBuffer, 'Spatial'), getCRS('albersNA', TRUE))
  calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
  
  bgFileName <- paste0('./in/bg_sites/Background Sites/Random Background Sites across Study Region - ', 
                       speciesAb, '.Rdata')
  
  # load bg sites in calibration region for a species if they have already been defined
  # loaded in: bgTestSp, bgCalib, bgEnv, bg
  # otherwise, get 20,000 random background sites from calibration region
  # we will only keep 10,000 points...this accounts for points that may fall in water
  if (file.exists(bgFileName)) load(bgFileName) else {
    bgTestSpAlb <- suppressWarnings(sp::spsample(calibRegionSpAlb, n=20000, 
                                                 type='random', iter = 10))
    bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
    bgCalib <- as.data.frame(coordinates(bgTestSp))
    names(bgCalib) <- ll
    save(bgTestSp, bgCalib, file = bgFileName, compress = T)
  }
  
  # extract environment at random background sites
  climate <- envData
  bgEnv <- raster::extract(climate, bgCalib)
  bgEnv <- as.data.frame(bgEnv)
  
  # remove any sites with NA for at least one variable
  isNa <- is.na(rowSums(bgEnv))
  if (any(isNa)) {
    bgCalib <- bgCalib[-which(isNa), ]
    bgEnv <- bgEnv[-which(isNa), ]
  }
  
  # only keep 10,000 background sites
  bgCalib <- bgCalib[1:10000,]
  bgEnv <- bgEnv[1:10000,]
  
  # combine with coordinates and rename coordinate fields
  bg <- cbind(bgCalib, bgEnv)
  names(bg)[1:2] <- ll
  
  presBg <- c(rep(1, nrow(records)), rep(0, nrow(bg)))
  occsEnv <- occsEnv[complete.cases(occsEnv), ]
  
  env <- rbind(occsEnv, bgEnv)
  env <- cbind(presBg, env)
  env <- as.data.frame(env)
  
  env <- env[complete.cases(env), ]
  
  # model species
  # envModel <- enmSdm::trainMaxNet(data = env, resp = 'presBg')
  envModel_tune <- enmSdm::trainMaxNet(data = env, resp = 'presBg',
                                       classes = 'lpq', out = c('models', 'tuning'))
  envModel <- envModel_tune$models[[95]]
  
  predictors <- c(paste0('pca', 1:pc))
  # prediction
  # envMap <- predict(climate[[predictors]], envModel, clamp = F, type = 'cloglog')
  envMap <- predict(
    climate[[predictors]],
    envModel,
    filename = paste0('./in/models/maxent/', speciesAb_, '_Maxent/GCM_', gcm,
                      '_PC', pc, '_', climYear, 'ybp'),
    clamp = F,
    format='GTiff',
    overwrite = T,
    type='cloglog')
  
  envMapSp <- rasterToPolygons(envMap)
  
  plot(rangeMap, border = 'blue', main = paste0('Maxent output, ', sp,
                                                ' occurrences'))
  plot(envMap, add = TRUE)
  plot(rangeMap, border = 'blue', add = TRUE)
  map("state", add = TRUE)
  map("world", add = TRUE)
  points(records$longitude, records$latitude, pch = 16, cex = 0.6, col = 'red')
  
  plot(envMap, main = paste0('Maxent output, ', 
                             sp,
                             ' occurrences'))
  plot(rangeMap, border = 'blue', add = TRUE)
  
  outputFileName <<- paste0('./in/models/maxent/', speciesAb_, 
                            '_Maxent/Model_PC', pc, '_GCM_', gcm, '.rData')
  save(rangeMap, envModel_tune, envModel, records, file = outputFileName, overwrite = T)
  
  outputFileName <<- paste0('./in/models/maxent/all_model_outputs/', speciesAb_,
                            '_GCM', gcm, '_PC', pc, '.Rdata')
  save(rangeMap, envModel_tune, envModel, records, file = outputFileName, overwrite = T)
  
  save.image(paste0('./in/workspaces/05 - Modeling Workspace - ', speciesAb_,
                    ' Model Output - PC', pc, '_GCM_', gcm))
  
}