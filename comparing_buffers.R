# How does buffer size impact model predictions?
# Does a small buffer size overpredict?

# Setup

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

# set constants
ll <- c('longitude', 'latitude')
# all final Fraxinus species to model (8 total)
speciesList <- c('Fraxinus americana', 'Fraxinus caroliniana', 'Fraxinus cuspidata'
                 ,'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# load PCA rasters
fileName <- paste0('./data_and_analyses/env_data/Beyer/envDataClipped_0kybp_pca5.tif')
envData <- brick(fileName)
# rename raster layers to pc's
names(envData) <- paste0('pca', 1:5)
# check rasters
# plot(envData)

calcBuffers <- function(sp) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  
  # Load records for given species.
  recordsRaw <- paste0('./cleaning_records/', species, '_finalRecords.rData')
  load(recordsRaw)
  recordsRaw <- data.frame(speciesSf_filtered_final)
  recordsRaw$geometry <- gsub("[c()]", "", recordsRaw$geometry)
  recordsRaw <- separate(data = recordsRaw, 
                         col = 'geometry', 
                         into = ll, 
                         sep = "\\,")
  recordsRaw$longitude <- as.double(recordsRaw$longitude)
  recordsRaw$latitude <- as.double(recordsRaw$latitude)
  
  # Match environmental data to records
  occsEnv <- raster::extract(envData, 
                             cbind(recordsRaw$longitude, 
                                   recordsRaw$latitude))
  occsEnvDf <- as.data.frame(occsEnv)
  records <- cbind(recordsRaw, occsEnvDf)
  
  # remove records that fall in the water & visualize 
  # records in the water:
  if (any(is.na(rowSums(occsEnvDf)))) {
    water <- records[which(is.na(rowSums(occsEnvDf))), ] 
    water <- SpatialPointsDataFrame(water[,ll], data = water,
                                    proj4 = getCRS('wgs84', TRUE))
  }
  
  # remove records in water from dataset:
  if (any(is.na(rowSums(occsEnvDf)))) records <-
    records[-which(is.na(rowSums(occsEnvDf))), ]
  
  load(paste0('./workspaces/03 - Modeling Workspace - ', speciesAb_, ' Cleaning'))
  rangeMap <- get(rangeName)
  
  bufferFileName <- paste0('./cleaning_records/', 
                           gsub(tolower(species), pattern = ' ', 
                                replacement = '_'), '_buffer.rData')
  load(bufferFileName)
  
  buffer_distance_1 <- buffer_distance
  buffer_distance_2 <- buffer_distance * 2
  buffer_distance_3 <- buffer_distance * 3
  
  # buffer 1
  rangeMapSf <- st_transform(st_as_sf(x = range), getCRS('albersNA'))
  range_buffer <- st_buffer(rangeMapSf, dist = buffer_distance_1)
  buffer_1 <- st_union(range_buffer)
  
  # use the buffer/calibration region for extracting background sites
  calibRegionSpAlb <- as(buffer_1, 'Spatial')
  calibRegionSpAlb <- sp::spTransform(calibRegionSpAlb, getCRS('albersNA', TRUE))
  calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
  
  # get 20,000 random background sites from calibration region
  bgTestSpAlb <- spsample(calibRegionSpAlb, n=20000, type='random')
  bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
  randomBgSites_1 <- as.data.frame(coordinates(bgTestSp))
  names(randomBgSites_1) <- ll
  
  # extract environment at background sites
  climate <- envData
  randomBgEnv_1 <- raster::extract(climate, randomBgSites_1)
  randomBgEnv_1 <- as.data.frame(randomBgEnv_1)
  
  # remove any sites with NA for at least one variable
  isNa <- is.na(rowSums(randomBgEnv_1))
  if (any(isNa)) {
    randomBgSites_1 <- randomBgSites_1[-which(isNa), ]
    randomBgEnv_1 <- randomBgEnv_1[-which(isNa), ]
  }
  
  # only keep 10,000 random background sites
  randomBgSites_1 <- randomBgSites_1[1:10000,]
  randomBgEnv_1 <- randomBgEnv_1[1:10000,]
  # combine with coordinates and rename coordinate fields
  randomBg_1 <- cbind(randomBgSites_1, randomBgEnv_1)
  names(randomBg_1)[1:2] <- ll
  
  presBg <- rep(c(1, 0), c(nrow(recordsRaw), nrow(randomBg_1)))
  
  env <- rbind(occsEnv, randomBgEnv_1)
  env <- cbind(presBg, env)
  env <- as.data.frame(env)
  env <- na.omit(env)
  
  # create output directory for model object and rasters
  dir.create(paste0('./Models/', speciesAb_, '_Maxent'),
             recursive=TRUE, showWarnings=FALSE)
  
  # model species
  envModel_1 <- enmSdm::trainMaxNet(data = env, resp = "presBg")
  # envModel <- maxnet(p = presBg, data = trainData)
  modelFileName <- paste0('./Models/', speciesAb_, '_Maxent/Model_0KYBP_pc5_buffer1_',
                          buffer_distance_1, 'km.Rdata')
  save(envModel_1, file = modelFileName, compress=TRUE)
  
  predictors <- c('pca1', 'pca2', 'pca3', 'pca4', 'pca5')
  # prediction
  envMap_1 <- predict(
    climate[[predictors]],
    envModel_1,
    filename = paste0('./Models/', speciesAb_, '_Maxent/maxentPredictionBeyer_0KYBP_pc5_buffer1_',
                      buffer_distance_1, 'km'), 
    clamp = F, 
    format='GTiff', 
    overwrite = TRUE, 
    type='cloglog')
  
  # buffer 2
  rangeMapSf <- st_transform(st_as_sf(x = range), getCRS('albersNA'))
  range_buffer <- st_buffer(rangeMapSf, dist = buffer_distance_2)
  buffer_2 <- st_union(range_buffer)
  
  # use the buffer/calibration region for extracting background sites
  calibRegionSpAlb <- as(buffer_2, 'Spatial')
  calibRegionSpAlb <- sp::spTransform(calibRegionSpAlb, getCRS('albersNA', TRUE))
  calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
  
  # get 20,000 random background sites from calibration region
  bgTestSpAlb <- spsample(calibRegionSpAlb, n=20000, type='random')
  bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
  randomBgSites_2 <- as.data.frame(coordinates(bgTestSp))
  names(randomBgSites_2) <- ll
  
  # extract environment at background sites
  climate <- envData
  randomBgEnv_2 <- raster::extract(climate, randomBgSites_2)
  randomBgEnv_2 <- as.data.frame(randomBgEnv_2)
  
  # remove any sites with NA for at least one variable
  isNa <- is.na(rowSums(randomBgEnv_2))
  if (any(isNa)) {
    randomBgSites_2 <- randomBgSites_2[-which(isNa), ]
    randomBgEnv_2 <- randomBgEnv_2[-which(isNa), ]
  }
  
  # only keep 10,000 random background sites
  randomBgSites_2 <- randomBgSites_2[1:10000,]
  randomBgEnv_2 <- randomBgEnv_2[1:10000,]
  # combine with coordinates and rename coordinate fields
  randomBg_2 <- cbind(randomBgSites_2, randomBgEnv_2)
  names(randomBg_2)[1:2] <- ll
  
  presBg <- rep(c(1, 0), c(nrow(recordsRaw), nrow(randomBg_2)))
  
  env <- rbind(occsEnv, randomBgEnv_2)
  env <- cbind(presBg, env)
  env <- as.data.frame(env)
  env <- na.omit(env)
  
  # create output directory for model object and rasters
  dir.create(paste0('./Models/', speciesAb_, '_Maxent'),
             recursive=TRUE, showWarnings=FALSE)
  
  # model species
  envModel_2 <- enmSdm::trainMaxNet(data = env, resp = "presBg")
  # envModel <- maxnet(p = presBg, data = trainData)
  modelFileName <- paste0('./Models/', speciesAb_, '_Maxent/Model_0KYBP_pc5_buffer2_',
                          buffer_distance_2, 'km.Rdata')
  save(envModel_2, file = modelFileName, compress=TRUE)
  
  predictors <- c('pca1', 'pca2', 'pca3', 'pca4', 'pca5')
  # prediction
  envMap_2 <- predict(
    climate[[predictors]],
    envModel_2,
    filename = paste0('./Models/', speciesAb_, '_Maxent/maxentPredictionBeyer_0KYBP_pc5_buffer2_',
                      buffer_distance_2, 'km'), 
    clamp = F, 
    format='GTiff', 
    overwrite = TRUE, 
    type='cloglog')
  
  # buffer 3
  range_buffer <- st_buffer(rangeMapSf, dist = buffer_distance_3)
  buffer_3 <- st_union(range_buffer)
  
  # use the buffer/calibration region for extracting background sites
  calibRegionSpAlb <- as(buffer_3, 'Spatial')
  calibRegionSpAlb <- sp::spTransform(calibRegionSpAlb, getCRS('albersNA', TRUE))
  calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
  
  # get 20,000 random background sites from calibration region
  bgTestSpAlb <- spsample(calibRegionSpAlb, n=20000, type='random')
  bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
  randomBgSites_3 <- as.data.frame(coordinates(bgTestSp))
  names(randomBgSites_3) <- ll
  
  # extract environment at background sites
  climate <- envData
  randomBgEnv_3 <- raster::extract(climate, randomBgSites_3)
  randomBgEnv_3 <- as.data.frame(randomBgEnv_3)
  
  # remove any sites with NA for at least one variable
  isNa <- is.na(rowSums(randomBgEnv_3))
  if (any(isNa)) {
    randomBgSites_3 <- randomBgSites_3[-which(isNa), ]
    randomBgEnv_3 <- randomBgEnv_3[-which(isNa), ]
  }
  
  # only keep 10,000 random background sites
  randomBgSites_3 <- randomBgSites_3[1:10000,]
  randomBgEnv_3 <- randomBgEnv_3[1:10000,]
  # combine with coordinates and rename coordinate fields
  randomBg_3 <- cbind(randomBgSites_3, randomBgEnv_3)
  names(randomBg_3)[1:2] <- ll
  
  presBg <- rep(c(1, 0), c(nrow(recordsRaw), nrow(randomBg_3)))
  
  env <- rbind(occsEnv, randomBgEnv_3)
  env <- cbind(presBg, env)
  env <- as.data.frame(env)
  env <- na.omit(env)
  
  # create output directory for model object and rasters
  dir.create(paste0('./Models/', speciesAb_, '_Maxent'),
             recursive=TRUE, showWarnings=FALSE)
  
  # model species
  envModel_3 <- enmSdm::trainMaxNet(data = env, resp = "presBg")
  # envModel <- maxnet(p = presBg, data = trainData)
  modelFileName <- paste0('./Models/', speciesAb_, '_Maxent/Model_0KYBP_pc5_buffer3_',
                          buffer_distance_3, 'km.Rdata')
  save(envModel_3, file = modelFileName, compress=TRUE)
  
  predictors <- c('pca1', 'pca2', 'pca3', 'pca4', 'pca5')
  # prediction
  envMap_3 <- predict(
    climate[[predictors]],
    envModel_3,
    filename = paste0('./Models/', speciesAb_, '_Maxent/maxentPredictionBeyer_0KYBP_pc5_buffer3_',
                      buffer_distance_3, 'km'), 
    clamp = F, 
    format='GTiff', 
    overwrite = TRUE, 
    type='cloglog')
  
  plot(envMap_1, main = substitute(paste('Maxent output, ', 
                                         italic(x),
                                         ' occurrences'), list(x=sp)),
       sub = paste0('Buffer = ', buffer_distance_1, ' km'))
  
  plot(envMap_2, main = substitute(paste('Maxent output, ', 
                                         italic(x),
                                         ' occurrences'), list(x=sp)),
       sub = paste0('Buffer = ', buffer_distance_2, ' km'))
  
  plot(envMap_3, main = substitute(paste('Maxent output, ', 
                                         italic(x),
                                         ' occurrences'), list(x=sp)),
       sub = paste0('Buffer = ', buffer_distance_3, ' km'))
  
  brick <- brick(envMap_1, envMap_2, envMap_3)
  fileName <- paste0('./data_and_analyses/buffer_output/', speciesAb_, '.tif')
  writeRaster(brick, fileName, format = 'GTiff', bylayer = T, 
              suffix = names(brick),overwrite = T)
}

lapply(speciesList, calcBuffers)

loadBufferRasters <- function(sp) {
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  fileNames <- list.files(path = './data_and_analyses/buffer_output/', 
                          pattern = speciesAb_, full.names = T)
  stack <- stack(fileNames)
  names(stack) <- c(paste0(substr(names(stack[[1]]), start = 1, stop = 9), 
                         " ", substr(names(stack[[1]]), start = 33, stop = 37), 
                         " Buffer = ", substr(names(stack[[1]]), start = 51, stop = 55)),
                    paste0(substr(names(stack[[2]]), start = 1, stop = 9), 
                           " ", substr(names(stack[[2]]), start = 33, stop = 37), 
                           " Buffer = ", substr(names(stack[[2]]), start = 51, stop = 55)),
                    paste0(substr(names(stack[[3]]), start = 1, stop = 9), 
                           " ", substr(names(stack[[3]]), start = 33, stop = 37), 
                           " Buffer = ", substr(names(stack[[3]]), start = 51, stop = 55)))
  plot(stack)
}

lapply(speciesList, loadBufferRasters)



