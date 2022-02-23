# Model Evaluation
# Date: 18 Feb 2022
# Author: Lauren Jenkins

rm(list = ls())
library(dismo)
library(sp)
library(enmSdm)

gcm <- 'Beyer'
pc <- 5
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')
sp <- 'Fraxinus pennsylvanica'
species <- gsub(tolower(sp), pattern=' ', replacement='_')
speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)

# load bg sites 
load(paste0('./Background Sites/Random Background Sites across Study Region - ', 
       speciesAb, '.Rdata'))
load(paste0('./workspaces/04 - Modeling Workspace - Clipping ', sp, '_PC_', pc, '_GCM_', gcm))

# calculate k-folds for presences and background sites
kPres <- kfold(records, k = 5)
kBg <- kfold(randomBg, k = 5)

# visualize fold #1
# rangeMap <- get(rangeName)
# plot(rangeMap, main='k-fold #1')
# points(records$longitude, records$latitude)
# points(records$longitude[kPres==1],
#        records$latitude[kPres==1],
#        bg='red',
#        pch=21
# )
# 
# legend('topright',
#        legend=c('Training presence', 'Test presence'),
#        pch=c(1, 16),
#        col=c('black', 'red'),
#        bg='white',
#        cex=0.8
# )

# create an output directory
dir.create(paste0('./Models/Maxent/', speciesAb_, '_Maxent/Model Evaluation - Random K-Folds'), 
           recursive = TRUE, showWarnings = FALSE)

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
  
  model <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg')
  save(model, file = paste0('./Models/Maxent/', speciesAb_, 
                                '_Maxent/Model Evaluation - Random K-Folds/Model ', 
                                i, '.Rdata'))
  
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
  
  print(paste('AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2)))
  
  aucRandom[i] <- thisAuc
  cbiRandom[i] <- thisCbi
}
