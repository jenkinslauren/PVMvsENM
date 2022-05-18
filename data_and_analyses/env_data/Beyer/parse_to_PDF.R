# Write each variable to pdf 

rm(list=ls())
setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
library(ncdf4)
library(raster)

# File path on Lauren's computer:
file <- "./data_and_analyses/env_data/Beyer/LateQuaternary_Environment.nc"
studyRegion <- rgdal::readOGR('/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region',
                              'study_region')
studyRegion <- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
projection(studyRegion) <- getCRS("WGS84")

plotVar <- function(varName) {
  if (grepl('bio', varName)) {
    upperVar <- toupper(varName)
  } else {
    upperVar <- varName
  }
  
  clim <- brick(file, varname = upperVar, level = 1)
  
  # pdf(file = paste0('./data_and_analyses/env_data/Beyer/climRasts_', upperVar, '.pdf'), 
  #     width = 11, height = 8.5)
  for (i in 1:nlayers(clim)) {
    x <- crop(clim[[i]], studyRegion)
    plot(x, main = names(x))
  }
  # dev.off()
}

vars <- c('bio1', paste0('bio', 4:19), 'cloudiness', 'relative_humidity')
lapply(vars, plotVar)

file <- list.files(path = './data_and_analyses/env_data/Beyer/tifs', pattern = '.tif', 
                   full.names = T)

for (i in 1:length(file)) {
  var <- gsub(".tif.*", "", gsub(".*tifs/", "", file[[i]]))
  pdf(file = paste0('./data_and_analyses/env_data/Beyer/climRasts_', var, '.pdf'), 
      width = 11, height = 8.5)
  x <- brick(file[[i]])
  for (j in 1:nlayers(x)) {
    xPlot <- crop(x[[j]], studyRegion)
    plot(xPlot, main = names(xPlot))
  }
  dev.off()
}

fileList <- list.files(path = './data_and_analyses/env_data/Beyer/tifs', 
                   pattern = '*.tif', full.names = T)
fileList <- mixedsort(fileList)
for (var in fileList) {
  pdf(file = paste0('./data_and_analyses/env_data/Beyer/pdf/', var, '.pdf'), width = 11, height = 8.5)
  x <- brick(var)
  for (i in 1:nlayers(x)) {
    xPlot <- crop(x[[i]], studyRegion)
    plot(xPlot, main = names(xPlot))
  }
  dev.off()
}

for (i in 1:length(file)) {
  x <- crop(brick(file[i]), studyRegion)
  plot(x, main = names(x))
}
dev.off()

pdf(file = './data_and_analyses/env_data/Beyer/climRastTest.pdf', width = 11, height = 8.5)
for (i in 1:length(means)) {
  x <- crop(means[[i]], studyRegion)
  plot(x, main = names(means[[i]]))
}
dev.off()
