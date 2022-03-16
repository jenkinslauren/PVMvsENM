# Model Evaluation
# Date: 18 Feb 2022
# Author: Lauren Jenkins

rm(list = ls())
library(dismo)
library(sp)
library(enmSdm)
library(xlsx)
library(geosphere)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
# setwd('/mnt/research/TIMBER/PVMvsENM')

ll <- c('longitude', 'latitude')

gcmList <- c('Beyer','Lorenz_ccsm', 'ecbilt')
pc <- 5
predictors <- c(paste0('pca', 1:pc))
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

for (gcm in gcmList) {
  gcm <- gcm
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
    load(paste0('./Background Sites/Random Background Sites across Study Region - ', 
                speciesAb, '_', gcm, '.Rdata'))
    load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
                '_PC', pc, '.rData'))
    load(paste0('./regions/little_range_map/', rangeName, '.Rdata'))
    
    # calculate k-folds for presences and background sites
    kPres <- kfold(records, k = 5)
    kBg <- kfold(randomBg, k = 5)
    
    # visualize fold #1
    plot(rangeMap, main = paste0(sp, ', k-fold #1'))
    points(records$longitude, records$latitude)
    points(records$longitude[kPres==1],
           records$latitude[kPres==1],
           bg='red',
           pch=21
    )
    
    legend('topright',
           legend=c('Training presence', 'Test presence'),
           pch=c(1, 16),
           col=c('black', 'red'),
           bg='white',
           cex=0.8
    )
    
    folderName <- paste0('./Models/Maxent/', speciesAb_,
                         '_Maxent/Model Evaluation - Random K-Folds - ', gcm)
    
    # folderName <- paste0('./in/models/maxent/', speciesAb_, 
    #                      '_Maxent/Model Evaluation - Random K-Folds - ', gcm)
    # create an output directory
    dir.create(folderName, recursive = TRUE, showWarnings = FALSE)
    
    # place to store auc & cbi output
    aucRandom <- cbiRandom <- rep(NA, 5)
    
    # for each k-fold
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
      
      # evaluate
      thisEval <- evaluate(p = as.vector(predPres), a = as.vector(predBg))
      thisAuc <- thisEval@auc
      thisCbi <- contBoyce(pres = predPres, bg = predBg)
      
      # print(paste('AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2)))
      
      aucRandom[i] <- thisAuc
      cbiRandom[i] <- thisCbi
    }
    
    save(aucRandom, cbiRandom, file = paste0(folderName, '/auc_cbi_vals.Rdata'))
    a <- a + 1
    
  }
}

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')
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
                         '_Maxent/Model Evaluation - Random K-Folds - ', gcm)
    load(paste0(folderName, '/auc_cbi_vals.Rdata'))
    
    a <- cbind(a, aucRandom)
    c <- cbind(c, cbiRandom)
    n <- ncol(a)
    colnames(a)[n] <- colnames(c)[n] <- sp
  }
  # save(a, c, file = paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
  save(a, c, file = paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
}

for (gcm in gcmList) {
  load(paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
  write.xlsx(a, file = './Models/Maxent/all_evals.xlsx', sheetName = paste0(gcm, '_auc'), 
             append = T, row.names = F)
  write.xlsx(c, file = './Models/Maxent/all_evals.xlsx', sheetName = paste0(gcm, '_cbi'),
             append = T, row.names = F)
}

# also predict to presences & background sites of non-k-fold sites
# predPres <- predict(tunedModel, records, type='cloglog')
# predBg <- predict(tunedModel, randomBg, type='cloglog')

# geographic
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

gcm <- 'Beyer'
sp <- speciesList[1]
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

gPres <- geoFold(x = records, k = 5, minIn = 5, minOut = 10, longLat = ll)

# maps
par(mfrow=c(1, 1), pty='s')
for (i in 1:5) {
  plot(rangeMap, main=paste0('g-fold #', i))
  # sp::plot(countries, add=TRUE)
  points(records$longitude, records$latitude)
  points(records$longitude[gPres==i],
         records$latitude[gPres==i],
         bg='red',
         pch=21
  )
  legend('topright',
         legend=c('Training presence', 'Test presence'),
         pch=c(1, 16),
         col=c('black', 'red'),
         bg='white',
         cex=0.8
  )
}

# initialize vectors to store g-fold assignments
gTrainBg <- rep(NA, nrow(records))
gTestBg <- rep(NA, nrow(randomBg))

# divide bg sites between training & test

