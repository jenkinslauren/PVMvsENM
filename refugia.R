rm(list = ls())
library(dismo)
library(sp)
library(enmSdm)
library(xlsx)
library(geosphere)
library(raster)
library(data.table)
library(dplyr)
library(viridis)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
# setwd('/mnt/research/TIMBER/PVMvsENM')

ll <- c('longitude', 'latitude')

gcmList <- c('Beyer','Lorenz_ccsm', 'ecbilt')
pc <- 5
predictors <- c(paste0('pca', 1:pc))
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')
thresholds <- list()

for(a in 1:length(speciesList)) {
  sp <- speciesList[a]
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  print(paste0("Species = ", species))
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  
  folderName <- paste0('./Models/Maxent/', speciesAb_,
                       '_Maxent/Model Evaluation - Random K-Folds - ', gcm)
  
  # load bg sites, records, and rangeMap
  load(paste0('./Background Sites/Random Background Sites across Study Region - ', 
              speciesAb, '.Rdata'))
  load(paste0('./Models/Maxent/all_model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  load(paste0('./data_and_analyses/study_region/regions/little_range_map/', 
              rangeName, '.Rdata'))
  # load(paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
  
  # load k-folds for presences and background sites from model evaluations
  t <- list()
  for(i in 1:5) {
    print(paste0('K-fold ', i, ':'))
    
    load(paste0(folderName, '/Model ', i, '.Rdata'))
    
    # temp <- enmSdm::thresholdWeighted(predPres, predBg, na.rm = T)
    temp <- enmSdm::thresholdWeighted(predPres, predBg, na.rm = T)
    if(temp['msss'] == 0) t <- append(t, NA)
    else t <- append(t, temp['msss'])
  }
  thresholds[[a]] <- t
  
}

thresholds <- data.frame(t(rbindlist(thresholds, fill = T)))
colnames(thresholds) <- speciesList
rownames(thresholds) <- paste("K-fold", 1:5) 

thresholds <- rbind(thresholds, mean = summarize_all(thresholds, mean))

fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)
t <- c(thresholds['mean',])

for(z in 1:length(t)) {
  threshold <- t[[z]]
  f <- fileName[z]
  s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
  speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
  load(paste0('./Models/Maxent/all_model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  # pdf(file = paste0('./predictions/PDF_output/', s, '.pdf'), width = 11, height = 8.5)
  b <- brick(f)
  
  b <- b[[1]]
  names(b) <- paste0(21, ' Kybp')
  
  title <- gsub('.*/', '', s)
  # par(mfrow=c(1,2))
  refugia <- b >= threshold
  refugiaId <- raster::clump(refugia, directions = 8, gaps = F)
  names(refugiaId) <- 'refugiaId'
  # plot(refugiaId, main = paste0(names(b),' ', title), axes = F)
  abund <- b * refugia
  names(abund) <- 'refugiaAbund'
  
  nrows <- nrow(b)
  ncols <- ncol(b)
  ncells <- raster::ncell(b)
  
  v <- rep(seq(nrows * (ncols - 1) - 1, 1, by=-ncols), each=ncols) + 0:(ncols - 1)
  cellNum <- matrix(v, nrow=nrows, ncol=ncols, byrow=TRUE)
  cellNum <- raster::raster(cellNum, template=b)
  cellNum <- as.vector(cellNum)
  simRefugiaBinary <- as.vector(refugia)
  refugeCellNum <- cellNum[simRefugiaBinary]
  if (any(is.na(refugeCellNum))) refugeCellNum <- refugeCellNum[!is.na(refugeCellNum)]
  
  # mean refuge abundance
  meanRefugeAbund <- raster::cellStats(abund, 'sum') / length(refugeCellNum)
  
  out <- list(
    simulationScale = raster::stack(refugiaId, abund),
    refugeCellNum = refugeCellNum,
    meanRefugeAbund = meanRefugeAbund
  )
  
  par(mfrow=c(1,2))
  plot(out$simulationScale[[1]], main = paste0('Refugia\n', speciesAb_, '\n', gcm), 
       col = viridis(256), axes = F)
  plot(out$simulationScale[[2]], main = paste0('Refugia abundance\n', speciesAb_, '\n', gcm), 
       col = viridis(256), axes = F)
  
  # dev.off()
}
