# Fraxinus Models
# Author: Lauren Jenkins
# Date: 20 December 2021

# 111 = full species name with first letter capital (Fraxinus pennsylvanica)
# 222 = full species name, all lowercase, _ inbetween (fraxinus_pennsylvanica)
# 333 = abreviated species name (ex. FraxPenn)
# 444 = abreviated with _ inbetween 

# Setup environment:
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

# load PCA environment
load('./PCA_Beyer')

# set constants
ll <- c('longitude', 'latitude')
speciesList <- c('Fraxinus americana', 'Fraxinus caroliniana', 'Fraxinus cuspidata'
                 ,'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# Load environmental PCA rasters. 
# If already clipped, load that. 
# Otherwise, clip the data to present.
fileName <- paste0('./data_and_analyses/env_data/Beyer/envDataClipped.tif')
if (file.exists(fileName)) {
  envData <- brick(fileName)
  # rename raster layers to pc's
  names(envData) <- paste0('pca', 1:5)
} else {
  pcPrediction <- list()
  # label each env layer by variable and year
  for (i in 1:length(clim)) {
    names(clim[[i]]) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                          "bio17", "bio18", "bio19", "bio4", "bio5", "bio6", "bio7",
                          "bio8", "bio9", "cloudiness", "relative_humidity")
    pcPrediction[i] <- raster::predict(clim[[i]], pca, index = 1:5)
    names(pcPrediction[[i]]) <- paste0("pc", 1:5, "_", (i-1)*1000, "KYBP")
  }
  
  envDataPca <- stack(pcPrediction)
  
  # keep only rasters (first three rasters) from present (0 KYBP)
  envPresences <- stack(envDataPca[[1:5]])
  names(envPresences) <- paste0('pca', 1:5)
  
  # check projections = wgs84
  print("Ensure that the projection of these rasters is WGS84:")
  print(paste0("Projection of envPresences = ", projection(envPresences)))
  
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
  for (i in 1:nlayers(envPresences)) {
    x <- envPresences[[i]]
    x <- crop(x, studyExtent)
    x <- mask(x, studyRegion)
    projection(x) <- getCRS("WGS84")
    envDataClipped[[i]] <- x
    envData <- stack(envDataClipped)
    # plot(envData)
    writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
  }
}


generateModel <- function(sp) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  
  print(paste0("SpeciesAb = ", speciesAb))
  print(paste0("SpeciesAb_ = ", speciesAb_))
  print(paste0("Range Name = ", rangeName))
  
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
  
  # visualize the points that fall in the water
  # convert to sp object
  recordsSp <- SpatialPointsDataFrame(records[, ll], data = records,
                                      proj4 = getCRS('wgs84', TRUE))
  
  # visualize
  plot(recordsSp, pch = 16, cex = 0.5, col = "red", 
       main = substitute(paste(italic(speciesList[i]), ' occurrences (BIEN) thinned')))
  if (exists("water")) {
    plot(water, col = 'blue', add = TRUE)
  }
  map("state", add = TRUE)
  map("world", add = TRUE)
  save.image(paste0('./04 - Modeling Workspace - Clipping ', sp))
  
  load(paste0('./03 - Modeling Workspace - ', speciesAb_, ' Cleaning'))
  rangeMap <- get(rangeName)
  
  # use the buffer/calibration region for extracting background sites
  calibRegionSpAlb <- as(range_buffer_final, 'Spatial')
  calibRegionSpAlb <- sp::spTransform(calibRegionSpAlb, getCRS('albersNA', TRUE))
  calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
  
  # get 10,000 random background sites from calibration region
  bgTestSpAlb <- spsample(calibRegionSpAlb, n=20000, type='random')
  bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
  randomBgSites <- as.data.frame(coordinates(bgTestSp))
  names(randomBgSites) <- ll
  
  # sanity check: plot the background sites
  plot(bgTestSp, pch = 16, cex = 0.5, col = "red", 
       main = substitute(paste(italic(speciesList[i]), ' background sites')))
  plot(calibRegionSpWgs, add = TRUE, border = 'blue')
  map("state", add = TRUE)
  map("world", add = TRUE)
  
  # extract environment at background sites
  climate <- envData
  randomBgEnv <- raster::extract(climate, randomBgSites)
  randomBgEnv <- as.data.frame(randomBgEnv)
  
  # remove any sites with NA for at least one variable
  isNa <- is.na(rowSums(randomBgEnv))
  if (any(isNa)) {
    randomBgSites <- randomBgSites[-which(isNa), ]
    randomBgEnv <- randomBgEnv[-which(isNa), ]
  }
  
  # only keep 10,000 random background sites
  randomBgSites <- randomBgSites[1:10000,]
  randomBgEnv <- randomBgEnv[1:10000,]
  
  # combine with coordinates and rename coordinate fields
  randomBg <- cbind(randomBgSites, randomBgEnv)
  names(randomBg)[1:2] <- ll
  
  dir.create('./Background Sites', recursive=TRUE, showWarnings=FALSE)
  save(randomBg, 
       file = paste0('./Background Sites/Random Background Sites across Study Region - ', 
                     speciesAb, '.Rdata'),
       compress=TRUE)
  
  presBg <- rep(c(1, 0), c(nrow(recordsRaw), nrow(randomBg)))
  
  env <- rbind(occsEnv, randomBgEnv)
  env <- cbind(presBg, env)
  env <- as.data.frame(env)
  
  env <- env[complete.cases(env), ]
  # View(env)
  
  # create output directory for model object and rasters
  dir.create(paste0('./Models/', speciesAb_, '_Maxent'),
             recursive=TRUE, showWarnings=FALSE)
  
  # model species
  envModel <- enmSdm::trainMaxNet(data = env, resp = "presBg")
  # envModel <- maxnet(p = presBg, data = trainData)
  modelFileName <- paste0('./Models/', speciesAb_, '_Maxent/Model.Rdata')
  save(envModel, file = modelFileName, compress=TRUE)
  
  predictors <- c('pca1', 'pca2', 'pca3', 'pca4', 'pca5')
  # prediction
  envMap <- predict(
    climate[[predictors]],
    envModel,
    filename=paste0('./Models/', speciesAb_, '_Maxent/maxentPredictionBeyer0KYBP'), 
    clamp = F, 
    format='GTiff', 
    overwrite = TRUE, 
    type='cloglog')
  
  envMapSp <- rasterToPolygons(envMap)
  
  # visualize 
  plot(rangeMap, border = 'blue', main = substitute(paste('Maxent output, ', 
                                                          italic(speciesList[i]),
                                                          ' occurrences')))
  plot(envMap, add = TRUE)
  plot(rangeMap, border = 'blue', add = TRUE)
  map("state", add = TRUE)
  map("world", add = TRUE)
  points(records$longitude, records$latitude, pch = 16, cex = 0.6, col = 'red')
  
  # zoomed out
  plot(envMap, main = substitute(paste('Maxent output, ', 
                                       italic(speciesList[i]),
                                       ' occurrences')))
  plot(rangeMap, border = 'blue', add = TRUE)
  
  i <- i + 1
  
}

i <- 1
lapply(speciesList, generateModel)

