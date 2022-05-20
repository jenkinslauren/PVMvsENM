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
gcm <- args[1]

# setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
setwd('/mnt/research/TIMBER/PVMvsENM')

pc <- 5
predictors <- c(paste0('pca', 1:pc))

# if (gcm == 'Lorenz_ccsm') speciesList <- paste('Fraxinus', c('cuspidata', 'greggii'))
# if (gcm == 'ecbilt') speciesList <- paste('Fraxinus', c('americana', 'cuspidata',
#                                                         'greggii', 'profunda'))

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

for(s in 1:length(speciesList)) {
  sp <- speciesList[s]
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  print(paste0("GCM = ", gcm))
  print(paste0("Species = ", species))
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  
  # load bg sites, records, and rangeMap
  # load(paste0('./in/bg_sites/Background Sites/Random Background Sites across Study Region - ', 
  #             speciesAb, '.Rdata'))
  load('./in/bg_sites/Background Sites/Random Background Sites across Study Region.Rdata')
  load(paste0('./in/models/maxent/all_model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.Rdata'))
  load(paste0('./in/workspaces/05 - Modeling Workspace - ', speciesAb_,
              ' Model Output - PC', pc, '_GCM_', gcm))
  # load(paste0('./regions/little_range_map/', rangeName, '.Rdata'))
  folderName <- paste0('./in/models/maxent/', speciesAb_,
                       '_Maxent/Model Evaluation - Geographic K-Folds - ', gcm)
  # create output directory for geographically distributed model evaluation
  dir.create(folderName, recursive = TRUE, showWarnings = FALSE)
  
  # for storing model evaluation metrics for each k-fold
  aucGeog <- cbiGeog <- rep(NA, 5)
  
  # create g-folds
  gPres <- geoFold(x = records, k = 5, minIn = 5, minOut = 10, longLat = c('longitude', 'latitude'))
  
  # now, we have our folds, but we want to divide the bg sites into fold based
  # on where they are in relation to records
  
  # initialize vectors to store g-fold assignments
  # gTrainBg <- rep(NA, nrow(randomBg))
  gTestBg <- rep(NA, nrow(bg))
  
  # convert records to sp object for gDistance function
  sp.records <- records
  coordinates(sp.records) <- ~longitude+latitude
  sp.randomBg <- bg
  coordinates(sp.randomBg) <- ~longitude+latitude
  
  # divide bg sites between training & test
  nearest <- apply(gDistance(sp.records, sp.randomBg, byid = T), 1, which.min)
  for (k in 1:nrow(bg)) {
    gTestBg[k] <- gPres[nearest[k]]
  }
  
  for (m in 1:5) {
    # make training data frame with predictors and vector of 1/0
    # for presence/background
    envData <- rbind(
      records[gPres!=m, predictors],
      bg[gTestBg!=m, predictors]
    )
    
    presBg <- c(rep(1, sum(gPres!=m)), rep(0, sum(gTestBg!=m)))
    
    trainData <- cbind(presBg, envData)
    
    # Maxent model
    # model <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg')
    model_tune <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg', 
                                      out = c('models', 'tuning'))
    model <- model_tune$models[[96]]
    
    # predict presences & background sites
    predPres <- raster::predict(model, 
                                newdata = records[gPres == m,],
                                clamp = F,
                                type = 'cloglog')
    predBg <- raster::predict(model, 
                              newdata = bg[gTestBg == m,],
                              clamp = F,
                              type = 'cloglog')
    
    save(model, predPres, predBg, gPres, gTestBg,
         model_tune, file = paste0(folderName, '/Model ', m, '.Rdata'), compress = T)
    
    # evaluate
    thisEval <- evaluate(p = as.vector(predPres), a = as.vector(predBg))
    thisAuc <- thisEval@auc
    thisCbi <- contBoyce(pres = predPres, bg = predBg)
    
    # print(paste('AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2)))
    
    aucGeog[m] <- thisAuc
    cbiGeog[m] <- thisCbi
  }
  save(aucGeog, cbiGeog, file = paste0(folderName, '/auc_cbi_vals.Rdata'))
}

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')
gcmList <- c('Beyer', 'ecbilt', 'Lorenz_ccsm')
for(gcm in gcmList) {
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
}


# for (gcm in gcmList) {
#   load(paste0('./Models/Maxent/', gcm, '_geoEvals.Rdata'))
#   write.xlsx(a, file = './Models/Maxent/geo_evals.xlsx', sheetName = paste0(gcm, '_auc'),
#              append = T, row.names = F)
#   write.xlsx(c, file = './Models/Maxent/geo_evals.xlsx', sheetName = paste0(gcm, '_cbi'),
#              append = T, row.names = F)
# }
