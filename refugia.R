rm(list = ls())
library(dismo)
library(sp)
library(enmSdm)
library(xlsx)
library(geosphere)
library(raster)
library(data.table)

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
              speciesAb_, '_', gcm, '.Rdata'))
  load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  load(paste0('./data_and_analyses/study_region/regions/little_range_map/', 
              rangeName, '.Rdata'))
  # load(paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
  
  # load k-folds for presences and background sites from model evaluations
  t <- list()
  for(i in 1:5) {
    print(paste0('K-fold ', i, ':'))
    
    load(paste0(folderName, '/Model ', i, '.Rdata'))
    
    temp <- enmSdm::thresholdWeighted(predPres, predBg)
    t <- append(t, temp['msss'])
  }
  thresholds[[a]] <- t
  
}

thresholds <- data.frame(t(rbindlist(thresholds, fill = T)))
colnames(thresholds) <- speciesList
rownames(thresholds) <- paste("K-fold", 1:5) 
