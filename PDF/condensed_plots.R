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
gcm <- 'Beyer'
pc <- 5

load(paste0('./workspaces/06 - ', gcm, ' Projections'))

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# set constants
climYears <- seq(21000, 0, by=-1000)

studyRegionFileName <- './regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

i <- 1
title <- paste0('Fraxinus, \nGCM = ', gcm)
pdf(file = '/Users/laurenjenkins/Downloads/test.pdf', width = 11, height = 8.5)
for(i in 1:22) {
  par(mfrow=c(2,5))
  plot(meansList[[i]], main = paste0('MEANS, ', title, ', ', names(meansList[[i]])), 
       col = colors, axes = F)
  plot(maxList[[i]], main = paste0('MAX, ', title, ', ', names(maxList[[i]])),  
       col = colors, axes = F)
  for(f in fileName) {
    s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
    speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
    load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
                '_PC', pc, '.rData'))
    b <- brick(f)
    names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
    title <- str_replace_all(gsub('./*GCM', '\nGCM = ', gsub('.*/', '', s)), '_', ' ')
    
    plot(b[[i]], main = paste0(gsub('\\X*', '', names(b[[i]])),'\n ', title),  col = colors, axes = F)
    
  }
}
dev.off()
