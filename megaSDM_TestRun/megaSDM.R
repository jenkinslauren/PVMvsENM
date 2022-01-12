rm(list = ls())
setwd('/Volumes/LJ MacBook Backup/MOBOT/PVMvsENM')

library(megaSDM)
library(sf)
library(raster)
library(rgdal)
library(dplyr)
library(ggplot2)
library(tidyr)

ll <- c('longitude', 'latitude')

speciesList <- c('Fraxinus americana', 'Fraxinus caroliniana', 'Fraxinus cuspidata'
                 ,'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

input_TA <- list.files('./data_and_analyses/env_data/Beyer/tifs', 
                       pattern = '.tif', full.names = TRUE)
input_TA <- input_TA[!grepl('.*pca_output.*\\.tif', input_TA)]

# where training & study rasters will be printed out to
envoutput <- "megaSDM_TestRun"

# define extent of training & study region
TSEnv <- TrainStudyEnv(input_TA = input_TA,
                       output = envoutput,
                       clipTrain = c(-178.2166, 83.7759, -55.90223, 83.6236),
                       clipStudy = c(-178.2166, 83.7759, -55.90223, 83.6236))

# PredictEnv clips, resamples, and reprojects the past environmental layers to 
## match the study region rasters created with TrainStudyEnv.
folders <- list()
for (i in (1:19)[-c(2,3)]) {
  folders <- append(folders, 
                    file.path('data_and_analyses','env_data','Beyer','tifs', paste0('BIO', i)))
}
folders <- append(folders, 
                  file.path('data_and_analyses','env_data','Beyer','tifs', 'cloudiness'))
folders <- append(folders, 
                  file.path('data_and_analyses','env_data','Beyer','tifs', 'relative_humidity'))

climYears <- seq(1000, 21000, by=1000)
envList <- list()
# create list of env rasters at each year backwards from 21000 ybp
# iterate through each folder in list and return all files
# unlist those lists of files into a single vector

envYears <- function(yr) {
  return(unlist(sapply(folders, function(folder) {
    paste0('./', list.files(folder, pattern = paste0("_", yr, "ybp.tif"), full.names=TRUE))
  })))
}
  
envList <- append(envList, lapply(climYears, envYears))

names <- list()
for (n in length(names(TSEnv$study))){
  for(i in (1:21)[-c(2,3)]) {
    for (yr in c(0,climYears)) {
      if (i == 20) {
        names <- append(names, paste0("cloudiness", "_", yr, "ybp"))
      } else if (i == 21) {
        names <- append(names, paste0("relative_humidity", "_", yr, "ybp"))
      } else {
        names <- append(names, paste0("BIO", i, "_", yr, "ybp"))
      }
    }
  }
}
names(TSEnv$study) <- c(names)

prediction <- PredictEnv(studylayers = TSEnv$study,
                         futurelayers = envList,
                         time_periods = c(0, climYears),
                         output = envoutput,
                         scenario_name = "Beyer")

toCSV <- function(sp) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  
  recordsRaw <- paste0('./cleaning_records/', species, '_finalRecords.rData')
  load(recordsRaw)
  fileName <- paste0('./megaSDM_TestRun/', species, '.csv')
  st_write(speciesSf_filtered_final, 
           fileName, 
           layer_options = "GEOMETRY=AS_XY",
           append = F)
  temp <- read.csv(fileName)
  colnames(temp)[1:3] <- c('longitude', 'latitude', 'species')
  write.csv(temp, fileName, row.names = FALSE)
}

lapply(speciesList, toCSV)

occlist <- list.files('./megaSDM_TestRun', pattern = '.csv', full.names = T)
occ_output <- "./megaSDM_TestRun/occurrences"

OccurrenceManagement(occlist = occlist,
                     output = occ_output,
                     envextract = TRUE,
                     envsample = TRUE,
                     nbins = 25,
                     envdata = TSEnv$training)

occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)
buff_output <- "megaSDM_TestRun/buffers"
BackgroundBuffers(occlist = occlist,
                  envdata = TSEnv$training,
                  buff_output,
                  ncores = 1)

# 10,000 background points
nbg <- 10000
# proportion of the background points sampled from within buffers
spatial_weights <- 0.5

# random or environmentally subsampled (Varela)?
sampleMethod <- "Varela"

bufflist <- list.files(buff_output, pattern = ".shp$", full.names = TRUE)
bg_output <- "megaSDM_TestRun/backgrounds"

BackgroundPoints(spplist = speciesList,
                 envdata = TSEnv$training,
                 output = bg_output,
                 nbg = nbg,
                 spatial_weights = spatial_weights,
                 buffers = bufflist,
                 method = sampleMethod,
                 ncores = 1)

occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)
model_output <- "megaSDM_TestRun/models"

MaxEntModel(occlist = occlist,
            bglist = bglist,
            model_output = model_output,
            ncores = 1,
            nrep = 4,
            alloutputs = FALSE)
