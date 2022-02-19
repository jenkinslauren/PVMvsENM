# Model Evaluation
# Date: 18 Feb 2022
# Author: Lauren Jenkins

rm(list = ls())
library(dismo)

gcm <- 'Beyer'

# load image for evaluation for a particular species
load(paste0('./workspaces/06 - ', gcm, ' Projections'))

# calculate k-folds for presences and background sites
kPres <- kfold(records, k = 5)
kBg <- kfold(targetSites, k=5)

aucRandom <- cbiRandom <- rep(NA, 5)