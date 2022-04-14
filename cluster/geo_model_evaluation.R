# Geographic Model Evaluation
# Date: 17 March 2022
# Author: Lauren Jenkins

rm(list = ls())
library(dismo)
library(sp)
library(enmSdm)
library(xlsx)
library(geosphere)
library(rgeos)

args <- commandArgs(TRUE)
gcm <- as.numeric(args[1])

# setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
setwd('/mnt/research/TIMBER/PVMvsENM')

gcmList <- c('Beyer','Lorenz_ccsm', 'ecbilt')
pc <- 5
predictors <- c(paste0('pca', 1:pc))
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')


a <- 1
for(i in 1:length(speciesList)) {
  rm(list= ls()[!(ls() %in% c('gcm','pc', 'predictors', 'speciesList', 'a'))])
  sp <- speciesList[a]
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  print(paste0("Species = ", species))
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  
  # load bg sites, records, and rangeMap
  load(paste0('./in/bg_sites/Random Background Sites across Study Region - ', 
              speciesAb, '_', gcm, '.Rdata'))
  load(paste0('./in/models/maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  # load(paste0('./regions/little_range_map/', rangeName, '.Rdata'))
  folderName <- paste0('./in/models/maxent/', speciesAb_,
                       '_Maxent/Model Evaluation - Geographic K-Folds - ', gcm)
  # create output directory for geographically distributed model evaluation
  dir.create(folderName)
  
  # for storing model evaluation metrics for each k-fold
  aucGeog <- cbiGeog <- rep(NA, 5)
  
  # create g-folds
  gPres <- geoFold(x = records, k = 5, minIn = 5, minOut = 10, longLat = c('longitude', 'latitude'))
  
  # now, we have our folds, but we want to divide the bg sites into fold based
  # on where they are in relation to records
  
  # initialize vectors to store g-fold assignments
  # gTrainBg <- rep(NA, nrow(randomBg))
  gTestBg <- rep(NA, nrow(randomBg))
  
  # convert records to sp object for gDistance function
  sp.records <- records
  coordinates(sp.records) <- ~longitude+latitude
  sp.randomBg <- randomBg
  coordinates(sp.randomBg) <- ~longitude+latitude
  
  # divide bg sites between training & test
  nearest <- apply(gDistance(sp.records, sp.randomBg, byid = T), 1, which.min)
  for (i in 1:nrow(randomBg)) {
    gTestBg[i] <- gPres[nearest[i]]
  }
  
  for (i in 1:5) {
    # make training data frame with predictors and vector of 1/0
    # for presence/background
    envData <- rbind(
      records[gPres!=i, predictors],
      randomBg[gTestBg!=i, predictors]
    )
    
    presBg <- c(rep(1, sum(gPres!=i)), rep(0, sum(gTestBg!=i)))
    
    trainData <- cbind(presBg, envData)
    
    # Maxent model
    model <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg')
    save(model, file = paste0(folderName, '/Model ', i, '.Rdata'), compress = T)
    
    # predict presences & background sites
    predPres <- raster::predict(model, 
                                newdata = records[gPres == i,],
                                clamp = F,
                                type = 'cloglog')
    predBg <- raster::predict(model, 
                              newdata = randomBg[gTestBg == i,],
                              clamp = F,
                              type = 'cloglog')
    
    # evaluate
    thisEval <- evaluate(p = as.vector(predPres), a = as.vector(predBg))
    thisAuc <- thisEval@auc
    thisCbi <- contBoyce(pres = predPres, bg = predBg)
    
    # print(paste('AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2)))
    
    aucGeog[i] <- thisAuc
    cbiGeog[i] <- thisCbi
  }
  save(aucGeog, cbiGeog, file = paste0(folderName, '/auc_cbi_vals.Rdata'))
  a <- a + 1
}


speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')


a <- data.frame(c(seq(1:5)))
c <- data.frame(c(seq(1:5)))
colnames(a)[1] <- colnames(c)[1] <- 'fold #'
for(sp in speciesList) {
  sp <- sp
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  
  folderName <- paste0('./Models/Maxent/', speciesAb_, 
                       '_Maxent/Model Evaluation - Geographic K-Folds - ', gcm)
  load(paste0(folderName, '/auc_cbi_vals.Rdata'))
  
  a <- cbind(a, aucGeog)
  c <- cbind(c, cbiGeog)
  n <- ncol(a)
  colnames(a)[n] <- colnames(c)[n] <- sp
}
# save(a, c, file = paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
save(a, c, file = paste0('./Models/Maxent/', gcm, '_geoEvals.Rdata'))


# for (gcm in gcmList) {
#   load(paste0('./Models/Maxent/', gcm, '_geoEvals.Rdata'))
#   write.xlsx(a, file = './Models/Maxent/all_geoEvals.xlsx', sheetName = paste0(gcm, '_auc'), 
#              append = T, row.names = F)
#   write.xlsx(c, file = './Models/Maxent/all_geoEvals.xlsx', sheetName = paste0(gcm, '_cbi'),
#              append = T, row.names = F)
# }
