rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

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
climYears <- seq(21000, 0, by = -1000)

studyRegionFileName <- './regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

title <- paste0('Fraxinus, \nGCM = ', gcm)
pdf(file = paste0('./PDF/', gcm, '_predictions_allSp.pdf'), width = 11, height = 8.5)
for(i in 1:22) {
  par(mfrow=c(2,5), mar=c(2,1,5,1)+0.1)
  plot(meansList[[i]], 
       main = paste0(sub('\\.', ' ', 
                         gsub('\\X*', '', names(meansList[[i]]))), '\nMEANS, ', title), 
       col = colors, axes = F, legend.mar = 10, box = F)
  plot(maxList[[i]], 
       main = paste0(sub('\\.', ' ', 
                         gsub('\\X*', '', names(maxList[[i]]))), '\nMAX, ', title),  
       col = colors, axes = F, legend.mar = 10, box = F)
  for(f in fileName) {
    s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
    speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
    load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
                '_PC', pc, '.rData'))
    b <- brick(f)
    names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
    title <- str_replace_all(gsub('./*GCM', '\nGCM = ', gsub('.*/', '', s)), '_', ' ')
    
    plot(b[[i]], main = paste0(sub('\\.', ' ', gsub('\\X*', '', names(b[[i]]))),'\n ', title),  
         col = colors, axes = F, legend.mar = 10, box = F)
    
  }
}
dev.off()

fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)
pollenRast <- brick('/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred.tif')

tmp <- list(brick(fileName[[1]]), brick(fileName[[2]]), brick(fileName[[3]]), 
            brick(fileName[[4]]), brick(fileName[[5]]), brick(fileName[[6]]), 
            brick(fileName[[7]]), brick(fileName[[8]]))
skip <- 1
means <- maxes <- list()
for (j in 1:length(tmp)){
  meansSkipList <- list()
  maxSkipList <- list()
  for(k in 1:length(climYears)) {
    nList <- list()
    for(m in 1:length(tmp)) {
      if(skip != m) { # species to include
        # print(names(tmp[[m]][[k]]))
        # thisRast <- raster::resample(tmp[[m]][[k]], pollenRast[[k]], method = 'bilinear')
        thisRast <- tmp[[m]][[k]]
        nList <- append(nList, brick(thisRast))
      } else { # species to remove
        skipped <- speciesList[m]
      }
    }
    
    climYear <- climYears[k]
    n <- stack(nList)
    mnSkip <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]],
                             n[[7]], fun = mean)
    names(mnSkip) <- paste0(climYear, ' ybp, Species skipped = ', skipped)
    meansSkipList <- append(meansSkipList, mnSkip)
    
    mxSkip <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]],
                             n[[7]], fun = max)
    names(mxSkip) <- paste0(climYear, ' ybp, Species skipped = ', skipped)
    maxSkipList <- append(maxSkipList, mxSkip)

  }
  
  means <- append(means, stack(meansSkipList))
  maxes <- append(maxes, stack(maxSkipList))
  
  # title <- paste0('Species skipped = ', skipped, ', GCM = ', gcm)
  # sp <- paste0(substr(skipped,1,4), toupper(substr(skipped,10,10)), 
  #              substr(skipped,11,13))
  # pdf(file = './PDF/predictions_removedSp.pdf', width = 11, height = 8.5)
  # for (i in 1:length(meansSkipList)) {
  #   par(mfrow=c(1,2))
  #   plot(meansSkipList[[i]], main = paste0('MEANS, ', names(meansSkipList[[i]])), 
  #        sub = title, col = colors, axes = F)
  #   plot(meansList[[i]], main = paste0('MEANS without sp removed, ', names(meansList[[i]])), 
  #        col = colors, axes = F)
  # }
  # dev.off()
  # 
  # pdf(file = paste0('./predictions/pdf/', gcm, '/No_', sp, '_maxList.pdf'), width = 11, height = 8.5)
  # for (i in 1:length(maxSkipList)) {
  #   par(mfrow=c(1,2))
  #   plot(maxSkipList[[i]], main = paste0('MAX, ', names(maxSkipList[[i]])), 
  #        sub = title, col = colors, axes = F)
  #   plot(maxList[[i]], main = paste0('MAX without sp removed, ', names(maxList[[i]])), 
  #        col = colors, axes = F)
  # }
  # dev.off()
  
  skip <- skip + 1
}

title <- paste0('Fraxinus, \nGCM = ', gcm)
pdf(file = paste0('./PDF/', gcm, '_predictions_removedSp_mean.pdf'), width = 11, height = 8.5)
for (i in 1:22) {
  par(mfrow=c(3,3), mar=c(2,1,5,1)+0.1)
  plot(meansList[[i]], main = paste0(sub('\\.', ' ', 
                                         gsub('\\X*', '', names(meansList[[i]]))), '\nMEANS, ', title), 
       col = colors, axes = F, legend.mar = 10, box = F)
  for(j in 1:8) {
    t <- gsub('\\.', ' ', gsub('X', '', names(means[[j]][[i]])))
    plot(means[[j]][[i]], main = paste0(sub('skipped ', 'skipped =\n', sub('ybp ', 'ybp\n', t))), 
         col = colors, axes = F, legend.mar = 10, box = F)
  }
}
dev.off()

pdf(file = paste0('./PDF/', gcm, '_predictions_removedSp_max.pdf'), width = 11, height = 8.5)
for (i in 1:22) {
  par(mfrow=c(3,3), mar=c(2,1,5,1)+0.1)
  plot(maxList[[i]], main = paste0(sub('\\.', ' ', 
                                       gsub('\\X*', '', names(maxList[[i]]))), '\nMAX, ', title), 
       col = colors, axes = F, legend.mar = 10, box = F)
  for(j in 1:8) {
    t <- gsub('\\.', ' ', gsub('X', '', names(maxes[[j]][[i]])))
    plot(maxes[[j]][[i]], main = paste0(sub('skipped ', 'skipped =\n', sub('ybp ', 'ybp\n', t))), 
         col = colors, axes = F, legend.mar = 10, box = F)
  }
}
dev.off()

##############################################################################
##############################################################################
rm(list = ls())

# constants
gcmList_ <- c('Lorenz_ccsm', 'ecbilt', 'Beyer', 'pollen')
pc <- 5

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# set constants
climYears <- seq(21000, 0, by=-1000)

studyRegionFileName <- './regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- as(world, "Spatial")

pollenRast <- brick('/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred.tif')

pdf(file = './PDF/all_gcm_predictions_means.pdf', 
    width = 11, height = 8.5)
par(mfrow=c(2,2), mar=c(2,1,5,1)+0.1)
for(a in 1:22) {
  for(gcm in gcmList_) {
    if(gcm == 'pollen') {
      t <- 'Fraxinus, \nGCM = Pollen'
      # scale values between 0 and 1
      mnv <- cellStats(pollenRast[[a]],'min')
      mxv <- cellStats(pollenRast[[a]],'max')
      pollenRast[[a]] <- (pollenRast[[a]] - mnv) / (mxv - mnv)
      
      plot(pollenRast[[a]], main = paste0(climYears[a], ' ybp,\n', t), 
           col = colors, axes = F, legend.mar = 10, box = F)
      plot(sp::spTransform(world, CRS(projection(pollenRast[[a]]))), border = 'black', add = T)
    } else {
      load(paste0('./workspaces/06 - ', gcm, ' Projections'))
      
      if(gcm == 'Lorenz_ccsm') t <- 'Fraxinus, \nGCM = CCSM'
      if(gcm == 'ecbilt') t <- 'Fraxinus, \nGCM = ECBilt'
      if(gcm == 'Beyer') t <- 'Fraxinus, \nGCM = HadAM3H'
      
      # scale values between 0 and 1
      mnv <- cellStats(meansList[[a]],'min')
      mxv <- cellStats(meansList[[a]],'max')
      meansList[[a]] <- (meansList[[a]] - mnv) / (mxv - mnv)
      
      plot(meansList[[a]], main = paste0(sub('\\.', ' ', 
                                             gsub('\\X*', '', names(meansList[[a]]))), '\nMEANS, ', t), 
           col = colors, axes = F, legend.mar = 10, box = F)
      plot(sp::spTransform(world, CRS(projection(meansList[[a]]))), border = 'black', add = T)
    }
  }
}
dev.off()

pdf(file = './PDF/all_gcm_predictions_max.pdf', 
    width = 11, height = 8.5)
par(mfrow=c(2,2), mar=c(2,1,5,1)+0.1)
for(a in 1:22) {
  for(gcm in gcmList_) {
    t <- paste0('Fraxinus, \nGCM = ', gcm)
    if(gcm == 'pollen') {
      # scale values between 0 and 1
      mnv <- cellStats(pollenRast[[a]],'min')
      mxv <- cellStats(pollenRast[[a]],'max')
      pollenRast[[a]] <- (pollenRast[[a]] - mnv) / (mxv - mnv)
      
      plot(pollenRast[[a]], main = paste0(climYears[a], ' ybp,\n', t), 
           col = colors, axes = F, legend.mar = 10, box = F)
      plot(sp::spTransform(world, CRS(projection(pollenRast[[a]]))), border = 'black', add = T)
    } else {
      load(paste0('./workspaces/06 - ', gcm, ' Projections'))
      
      # scale values between 0 and 1
      mnv <- cellStats(maxList[[a]],'min')
      mxv <- cellStats(maxList[[a]],'max')
      maxList[[a]] <- (maxList[[a]] - mnv) / (mxv - mnv)
      
      plot(maxList[[a]], main = paste0(sub('\\.', ' ', gsub('\\X*', '', names(maxList[[a]]))), '\nMAX, ', t), 
           col = colors, axes = F, legend.mar = 10, box = F)
      plot(sp::spTransform(world, CRS(projection(maxList[[a]]))), border = 'black', add = T)
    }
  }
}
dev.off()

