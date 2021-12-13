rm(list = ls())
library(raster)
library(rgdal)
library(enmSdm)
library(dplyr)
library(terra)
setwd('/Volumes/LJ MacBook Backup/MOBOT/PVMvsENM')

envDataPca <- brick('./data_and_analyses/env_data/Beyer/tifs/pca_output.tif')

if (file.exists('./data_and_analyses/study_region/Study Region.Rdata') & 
    file.exists('./data_and_analyses/study_region/Study Region Extent.Rdata')) {
  load('./data_and_analyses/study_region/Study Region.Rdata')
  load('./data_and_analyses/study_region/Study Region Extent.Rdata')
} else {
  studyRegion <- rgdal::readOGR(dsn = './data_and_analyses/study_region', 'study_region')
  studyRegion <- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
  projection(studyRegion) <- getCRS("WGS84")
  studyExtent <- extent(studyRegion)
  save(studyRegion, file='./data_and_analyses/study_region/Study Region.Rdata', compress=TRUE)
  save(studyExtent, file='./data_and_analyses/study_region/Study Region Extent.Rdata', compress=TRUE)
}

envDataClipped <- list()
for (i in 1:nlayers(envDataPca)) {
  x <- envDataPca[[i]]
  x <- crop(x, studyExtent)
  x <- mask(x, studyRegion)
  projection(x) <- getCRS("WGS84")
  envDataClipped[[i]] <- x
}

envData <- stack(envDataClipped)
plot(envData)
writeRaster(envData, './data_and_analyses/env_data/Beyer/envDataClipped.tif', 
            format = 'GTiff', overwrite = T)

# load data for focal species
ll <- c('longitude', 'latitude')
species <- 'fraxinus_pennsylvanica'
records <- paste0('./species_records/02_', species, '_thinned_records.rData')
load(records)
records <- finalThinned # 515 observations

# match environmental data to records
envSpecies <- extract(envData, cbind(records$longitude, records$latitude))
envSpecies <- as.data.frame(envSpecies)
records <- cbind(records, envSpecies)

# remove records that fall in the water
if (any(is.na(rowSums(envSpecies)))) records <-
  records[-which(is.na(rowSums(envSpecies))), ]

# visualize the points that fall in the water
# define these observations:
if (any(is.na(rowSums(envSpecies)))) water <- 
  records[which(is.na(rowSums(envSpecies))), ] 

# convert to sp object
finalThinned <- SpatialPointsDataFrame(finalThinned[, ll], data = finalThinned, 
                                       proj4 = getCRS('wgs84', TRUE))

# visualize
plot(finalThinned, pch = 16, cex = 0.5, col = "red", 
     main = 'BIEN occurrences thinned')
plot(water, col = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# save workspace
save.image('./04 - Modeling Workspace - Clipping Environmental Data')
