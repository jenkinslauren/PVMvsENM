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
library(units)

library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

# set constants
pc <- 5
gcmList <- c('Beyer', 'Lorenz_ccsm', 'Lorenz_ecbilt')
gcm <- 'Beyer'
climYear <- 0

# set constants
ll <- c('longitude', 'latitude')
sp <- 'Fraxinus quadrangulata'


# load PCA environment
load(paste0('./PCA_', gcm, '_PC', pc))
load('./workspaces/01 - Modeling Workspace - Fraxinus Range Maps (BIEN + Little)')

species <<- gsub(tolower(sp), pattern=' ', replacement='_')
speciesAb <<- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
speciesAb_ <<- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
rangeName <<- paste0('littleRange_', speciesAb)

# Load environmental PCA rasters. 
# If already clipped, load that. 
# Otherwise, clip the data to present.
studyRegionFileName <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/0ByjNJEf91IW5SUlEOUJFVGN0Y28/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

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
  pcPrediction <<- list()
  # label each env layer by variable and year
  pcPrediction <- list()
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

occsEnv <- raster::extract(envData, 
                           cbind(recordsRaw$longitude, 
                                 recordsRaw$latitude))
occsEnvDf <- as.data.frame(occsEnv)
records <- cbind(recordsRaw, occsEnvDf)

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

# save.image(paste0('./workspaces/04 - Modeling Workspace - Clipping ', sp, '_PC_', pc, '_GCM_', gcm))
# load(paste0('./04 - Modeling Workspace - Clipping ', sp))

# load(paste0('./workspaces/03 - Modeling Workspace - ', speciesAb, ' Cleaning'))
rangeMap <- get(rangeName)

bufferFileName <- paste0('./cleaning_records/', species, '_buffer.rData')
load(bufferFileName)
load(paste0("./workspaces/03 - Modeling Workspace - ", speciesAb, " Cleaning"))

# use the buffer/calibration region for extracting background sites
print(paste0('buffer size = ', buffer_distance))
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
     main = substitute(paste(italic('Fraxinus quadrangulata'), ' background sites')))
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
dir.create(paste0('./Models/', speciesAb_, '_Maxent_PC', pc, '_GCM_', gcm),
           recursive=TRUE, showWarnings=FALSE)

# model species
# try: trainMaxNet(data, resp='presBg', regMult=c(0.5, 1, 2, 3, 4, 5, 7.5, 10))
envModel <- enmSdm::trainMaxNet(data = env, resp = "presBg", regMult=c(0.5, 1, 2, 3, 4, 5, 7.5, 10))
# envModel <- maxnet(p = presBg, data = trainData)
modelFileName <- paste0('./Models/', speciesAb_, '_Maxent/Model_PC', pc, '_GCM_', gcm, '.Rdata')
save(envModel, file = modelFileName, compress=TRUE)

predictors <- c(paste0('pca', 1:pc))
# prediction
envMap <- predict(
  climate[[predictors]],
  envModel,
  filename = paste0('./Models/', speciesAb_, '_Maxent/maxentPredictionBeyer0KYBP_PC', 
                    pc, '_GCM', gcm), 
  clamp = F,
  format='GTiff', 
  overwrite = TRUE, 
  type='cloglog')

envMapSp <- rasterToPolygons(envMap)

plot(rangeMap, border = 'blue', main = substitute(paste('Maxent output, ', 
                                                        italic('Fraxinus quadrangulata '),
                                                        'occurrences')))
plot(envMap, add = TRUE)
plot(rangeMap, border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)
points(records$longitude, records$latitude, pch = 16, cex = 0.6, col = 'red')

plot(envMap, main = substitute(paste('Maxent output, ', 
                                     italic('Fraxinus quadrangulata '),
                                     'occurrences')))
plot(rangeMap, border = 'blue', add = TRUE)

# save.image(paste0('./workspaces/05 - Modeling Workspace - ', speciesAb_, 
#                   ' Model Output - PC', pc, '_GCM_', gcm))



lapply(speciesToModel, buildModel)
