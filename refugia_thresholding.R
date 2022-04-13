rm(list = ls())
library(dismo)
library(sp)
library(enmSdm)
library(xlsx)
library(geosphere)
library(raster)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
# setwd('/mnt/research/TIMBER/PVMvsENM')

ll <- c('longitude', 'latitude')

gcmList <- c('Beyer','Lorenz_ccsm', 'ecbilt')
pc <- 5
predictors <- c(paste0('pca', 1:pc))
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

folderName <- paste0('./Models/Maxent/', speciesAb_,
                     '_Maxent/Model Evaluation - Random K-Folds - ', gcm)

thresholds <- data.frame()

sp <- speciesList[a]
species <- gsub(tolower(sp), pattern=' ', replacement='_')
print(paste0("Species = ", species))
speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
rangeName <- paste0('littleRange_', speciesAb)

# load bg sites, records, and rangeMap
load(paste0('./Background Sites/Random Background Sites across Study Region - ', 
            speciesAb, '_', gcm, '.Rdata'))
load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
            '_PC', pc, '.rData'))
load(paste0('./regions/little_range_map/', rangeName, '.Rdata'))

# calculate k-folds for presences and background sites
kPres <- kfold(records, k = 5)
kBg <- kfold(randomBg, k = 5)

for(i in 1:5) {
  print(paste0('K-fold ', i, ':'))
  
  # create training data, with presences/absences vector of 0/1 with all points 
  # EXCEPT the ones in the fold
  envData <- rbind(records[kPres != i, predictors], randomBg[kBg != i, predictors])
  presBg <- c(rep(1, sum(kPres != i)), rep(0, sum(kBg != i)))
  trainData <- cbind(presBg, envData)
  
  # load(paste0(folderName, '/Model ', i, '.Rdata'))
  model <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg')
  save(model, file = paste0(folderName, '/Model ', i, '.Rdata'))
  
  # predict presences & background sites
  predPres <- raster::predict(model, 
                              newdata = records[kPres == i,],
                              clamp = F,
                              type = 'cloglog')
  predBg <- raster::predict(model, 
                            newdata = randomBg[kPres == i,],
                            clamp = F,
                            type = 'cloglog')
  
  enmSdm::thresholdWeighted(predPres, predBg)
  
}