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

# load(paste0('./workspaces/06 - ', gcm, ' Projections'))
# load(paste0('./workspaces/PCA_', gcm, '_PC', pc))

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# set constants
climYears <- seq(21000, 0, by=-1000)

studyRegionFileName <- './regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

getClimRasts <- function(pc, climYear) {
  
  if (gcm == 'Beyer') { # Beyer
    load('./data_and_analyses/env_data/Beyer/PCA_clim.Rdata')
    load(paste0('./data_and_analyses/env_data/Beyer/pca_pc', pc, '.Rdata'))
    fileName <- paste0('./data_and_analyses/env_data/Beyer/envDataClipped_',
                       climYear, 'KYBP_pc', pc, '.tif')
    vars <- c("BIO1", paste0('BIO', 4:19), "cloudiness", "relative_humidity")
  } else if (gcm == 'Lorenz_ccsm') { # CCSM
    load(paste0('./data_and_analyses/env_data/Lorenz/PCA_', gcm, '_clim.Rdata')) 
    load(paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/pca_pc', pc, '.Rdata')) 
    fileName <- paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/envDataClipped_',
                       climYear, 'KYBP_pc', pc, '.tif')
    clim <- lorenz
    workingFolder <- './data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/vars'
    vars <- sub('\\.tif.*', '', list.files(path = workingFolder, 
                                           pattern='*.tif', all.files = TRUE, full.names = FALSE))
    vars <- vars[lapply(vars, function(x) length(grep("pca_", x, value = F))) == 0]
  } else { # ECBilt
    load(paste0('./data_and_analyses/env_data/Lorenz/PCA_', gcm, '_clim.Rdata')) 
    load(paste0('./data_and_analyses/env_data/Lorenz/V2/ecbilt_21-0k_all_tifs_LJ/pca_pc', pc, '.Rdata'))
    fileName <- paste0('./data_and_analyses/env_data/Lorenz/V2/ecbilt_21-0k_all_tifs_LJ/envDataClipped_',
                       climYear, 'KYBP_pc', pc, '.tif')
    clim <- lorenz
    workingFolder <- paste0('./data_and_analyses/env_data/Lorenz/V2/',
                            gcm, '_21-0k_all_tifs_LJ/vars')
    vars <- sub('\\.tif.*', '', list.files(path = workingFolder, 
                                           pattern='*.tif', all.files = TRUE, full.names = FALSE))
    vars <- vars[lapply(vars, function(x) length(grep("pca_", x, value = F))) == 0]
  }
  
  if (file.exists(fileName)) {
    envData <- brick(fileName)
    # rename raster layers to pc's
    names(envData) <- paste0('pca', 1:pc)
  } else {
    pcPrediction <- list()
    for (i in 1:length(lorenz)) {
      names(lorenz[[i]]) <- vars
      pcPrediction[i] <- raster::predict(clim[[i]], pca, index = 1:pc)
      names(pcPrediction[[i]]) <- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
    }
    envDataPca <- stack(pcPrediction)
    
    # keep only rasters (first five rasters) for climate year
    # climYear <- as.numeric(sub("KYBP.*", "", sub(".*_", "", names(envDataPca)[k*5])))
    #print(paste0("clim year = ", climYear))
    envYr <- pcPrediction[[(climYear/1000) + 1]]
    # plot(envYr) 
    names(envYr) <- paste0('pca', 1:pc)
    
    # check projections = wgs84
    # print("Ensure that the projection of these rasters is WGS84:")
    # print(paste0("Projection of envYr = ", projection(envYr)))
    # 
    # define study region & extent, load if already defined
    if (file.exists('./data_and_analyses/study_region/Study Region.Rdata') & 
        file.exists('./data_and_analyses/study_region/Study Region Extent.Rdata')) {
      load('./data_and_analyses/study_region/Study Region.Rdata')
      load('./data_and_analyses/study_region/Study Region Extent.Rdata')
    } else {
      studyRegion <- rgdal::readOGR('/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region',
                                    'study_region')
      studyRegion <- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
      projection(studyRegion) <- getCRS("WGS84")
      studyExtent <- extent(studyRegion)
      save(studyRegion, 
           file='/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region/Study Region.Rdata',
           compress=TRUE)
      save(studyExtent, 
           file='/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region/Study Region Extent.Rdata', compress=TRUE)
    }
    
    # clip environmental PCAs to study extent for given species, visualize, and save:
    envDataClipped <- list()
    for (i in 1:nlayers(envYr)) {
      x <- envYr[[i]]
      x <- crop(x, studyExtent)
      # x <<- mask(x, studyRegion)
      projection(x) <- getCRS("WGS84")
      envDataClipped[[i]] <- x
      envData <- stack(envDataClipped)
      # plot(envData)
      writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
    }
  }
  
  return(envData)
}

predictions <- function(speciesAb_, pc) {
  
  predictors <<- c(paste0('pca', 1:pc))
  
  for (climYear in climYears) {
    climate <<- getClimRasts(pc = 5, climYear = climYear)
    # plot(climate[[1]], main = paste0(climYear, ' ybp'))
    climate <<- climate[[predictors]]
    thisPred <<- predict(climate, envModel, clamp = F, type='cloglog')
    # plot(thisPred, main = paste0("prediction at ", climYear, " ybp"))
    names(thisPred) <<- (paste0(climYear, ' ybp'))
    preds <<- if (exists('preds')) {
      stack(preds, thisPred)
    } else {
      thisPred
    }
  }
  return(preds)
}

for(sp in speciesList) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  
  fileName <- list.files(path = paste0('./Models/Maxent/', speciesAb_, '_Maxent'),
                         pattern = paste0(gcm, '.Rdata'),
                         full.names = T)
  load(fileName)
  
  if(exists('preds')) rm(preds)
  preds <- predictions(speciesAb_, pc)
  
  # interpFrom <- -1 * climYears
  # interpTo <- seq(-21000, 0, by=1000)
  # predsInterp <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')
  
  # project
  preds <- projectRaster(preds, studyRegionRasts)
  
  # mask by study region and force values to be within [0, 1] 
  # (can get pushed outside this during re-projection)
  preds <- raster::calc(preds, fun = function(x) ifelse(x < 0, 0, x))
  preds <- raster::calc(preds, fun = function(x) ifelse(x > 1, 1, x))
  for (i in 1:nlayers(preds)) {
    landMask <- (1 - studyRegionRasts[[i]])
    preds[[i]] <- preds[[i]] * landMask
  }
  
  # names(preds) <- paste0('ybp', seq(21000, 0, by=-1000))
  
  writeRaster(stack(preds), paste0('./predictions/', gcm, '/', speciesAb_, '_GCM_', gcm, '_PC', pc),
              format = 'GTiff', overwrite = T)
  
  # save.image(paste0('./workspaces/06 - Projections ENM ', speciesAb_))
}

gcm <- 'Beyer'
library(RColorBrewer)
colors <- c('#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850')
mergedRange <- readRDS('./littleMergedRangeFraxinus.rds')

fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)

for(f in fileName) {
  s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
  speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
  load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  pdf(file = paste0('./predictions/pdf/', s, '.pdf'), width = 11, height = 8.5)
  b <- brick(f)
  names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
  title <- gsub('.*/', '', s)
  for (i in 1:22) {
    par(mfrow=c(1,2))
    plot(b[[i]], main = paste0(names(b[[i]]),' ', title),  col = colors, axes = F)
    plot(b[[22]], main = paste0(names(b[[22]]),' ', title),  col = colors, axes = F)
  }
  dev.off()
}

tmp <- list(brick(fileName[[1]]), brick(fileName[[2]]), brick(fileName[[3]]), 
            brick(fileName[[4]]), brick(fileName[[5]]), brick(fileName[[6]]), 
            brick(fileName[[7]]), brick(fileName[[8]]))

meansList <- list()
maxList <- list()
# sumList <- list()

for (i in 1:length(climYears)) {
  climYear <- climYears[i]
  n <- stack(tmp[[1]][[i]], tmp[[2]][[i]], tmp[[3]][[i]], tmp[[4]][[i]],
             tmp[[5]][[i]], tmp[[6]][[i]], tmp[[7]][[i]], tmp[[8]][[i]])
  mn <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]],
                       n[[7]], n[[8]], fun = mean)
  names(mn) <- paste0(climYear, ' ybp')
  meansList <- append(meansList, mn)
  
  mx <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]],
                       n[[7]], n[[8]], fun = max)
  names(mx) <- paste0(climYear, ' ybp')
  maxList <- append(maxList, mx)
  
  # sm <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]], 
  #                      n[[7]], n[[8]], fun = sum)
  # mask <- n[[1]] * 0 + 1
  # sumRasterCorrected <- sm * mask
  # names(sumRasterCorrected) <- paste0(climYear, ' ybp')
  # sumList <- append(sumList, sumRasterCorrected)
}

fileName <- './predictions/pollen/predictions-FRAXINUS_meanpred.tif'
pollenRast <- brick(fileName)
lists <- list(meansList, maxList)
temp <- list()
for(l in lists) {
  resampled <- list()
  for(i in 1:length(l)) {
    thisRast <- raster::resample(l[[i]], pollenRast[[i]], method = 'bilinear')
    resampled <- append(resampled, thisRast)
  }
  temp <- append(stack(resampled), temp)
}

meansList <- stack(temp[1])
maxList <- stack(temp[2])

title <- paste0('Fraxinus, GCM = ', gcm)
pdf(file = paste0('./predictions/pdf/', gcm, '/', gcm, '_meansList.pdf'), width = 11, height = 8.5)
for (i in 1:nlayers(meansList)) {
  par(mfrow=c(1,2))
  plot(meansList[[i]], main = paste0('MEANS, ', title, ', ', names(meansList[[i]])), 
       col = colors, axes = F)
  plot(meansList[[22]], main = paste0('MEANS, ', title, ', ', names(meansList[[22]])), 
       col = colors, axes = F)
}
# plot(stack(meansList))
dev.off()

pdf(file = paste0('./predictions/pdf/', gcm, '/', gcm, '_maxList.pdf'), width = 11, height = 8.5)
for (i in 1:nlayers(maxList)) {
  par(mfrow=c(1,2))
  plot(maxList[[i]], main = paste0('MAX, ', title, ', ', names(maxList[[i]])),  
       col = colors, axes = F)
  plot(maxList[[22]], main = paste0('MAX, ', title, ', ', names(maxList[[22]])),  
       col = colors, axes = F)
}
# plot(stack(maxList))
dev.off()

# pdf(file = paste0('./predictions/pdf/', gcm, '/', gcm, '_sumList.pdf'), width = 11, height = 8.5)
# for (i in 1:length(sumList)) {
#   par(mfrow=c(1,2))
#   plot(sumList[[i]], main = paste0('SUM, ', title, ', ', names(sumList[[i]])),  
#        col = colors, axes = F)
#   plot(sumList[[22]], main = paste0('SUM, ', title, ', ', names(sumList[[22]])),  
#        col = colors, axes = F)
# }
# dev.off()

speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)

tmp <- list(brick(fileName[[1]]), brick(fileName[[2]]), brick(fileName[[3]]), 
            brick(fileName[[4]]), brick(fileName[[5]]), brick(fileName[[6]]), 
            brick(fileName[[7]]), brick(fileName[[8]]))
skip <- 1
for (j in 1:length(tmp)){
  meansSkipList <- list()
  maxSkipList <- list()
  # sumSkipList <- list()
  for(k in 1:length(climYears)) {
    nList <- list()
    for(m in 1:length(tmp)) {
      if(skip != m) { # species to include
        # print(names(tmp[[m]][[k]]))
        thisRast <- raster::resample(brick(tmp[[m]][[k]]), pollenRast[[i]], method = 'bilinear')
        nList <- append(nList, brick(thisRast))
      } else { # species to remove
        skipped <- speciesList[m]
      }
    }
    
    climYear <- climYears[k]
    n <- stack(nList)
    mnSkip <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]],
                             n[[7]], fun = mean)
    names(mnSkip) <- paste0(climYear, ' ybp')
    meansSkipList <- append(meansSkipList, mnSkip)
    
    mxSkip <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]],
                             n[[7]], fun = max)
    names(mxSkip) <- paste0(climYear, ' ybp')
    maxSkipList <- append(maxSkipList, mxSkip)
    
    # smSkip <- raster::mosaic(n[[1]], n[[2]], n[[3]], n[[4]], n[[5]], n[[6]], 
    #                          n[[7]], fun = sum)
    # mask <- n[[1]] * 0 + 1
    # sumSkipRasterCorrected <- smSkip * mask
    # names(sumSkipRasterCorrected) <- paste0(climYear, ' ybp')
    # sumSkipList <- append(sumSkipList, sumSkipRasterCorrected)
  }
  
  title <- paste0('Species skipped = ', skipped, ', GCM = ', gcm)
  sp <- paste0(substr(skipped,1,4), toupper(substr(skipped,10,10)), 
               substr(skipped,11,13))
  pdf(file = paste0('./predictions/pdf/', gcm, '/No_', sp, '_meansList.pdf'), 
      width = 11, height = 8.5)
  for (i in 1:length(meansSkipList)) {
    par(mfrow=c(1,2))
    plot(meansSkipList[[i]], main = paste0('MEANS, ', names(meansSkipList[[i]])), 
         sub = title, col = colors, axes = F)
    plot(meansList[[i]], main = paste0('MEANS without sp removed, ', names(meansList[[i]])), 
         col = colors, axes = F)
  }
  # plot(stack(meansList))
  dev.off()
  
  pdf(file = paste0('./predictions/pdf/', gcm, '/No_', sp, '_maxList.pdf'), width = 11, height = 8.5)
  for (i in 1:length(maxSkipList)) {
    par(mfrow=c(1,2))
    plot(maxSkipList[[i]], main = paste0('MAX, ', names(maxSkipList[[i]])), 
         sub = title, col = colors, axes = F)
    plot(maxList[[i]], main = paste0('MAX without sp removed, ', names(maxList[[i]])), 
         col = colors, axes = F)
  }
  # plot(stack(maxList))
  dev.off()
  
  # pdf(file = paste0('./predictions/pdf/', gcm, '/No_', sp, '_sumList.pdf'), width = 11, height = 8.5)
  # for (i in 1:length(sumSkipList)) {
  #   par(mfrow=c(1,2))
  #   plot(sumSkipList[[i]], main = paste0('SUM, ', names(sumSkipList[[i]])), 
  #        sub = title, col = colors, axes = F)
  #   plot(sumList[[i]], main = paste0('SUM, ', names(sumList[[i]])), 
  #       col = colors, axes = F)
  # }
  # dev.off()
  
  skip <- skip + 1
}

gcmList <- c('Lorenz_ccsm')
library(enmSdm)
library(gtools)

for(gcm in gcmList) {
  # load(paste0('./workspaces/06 - ', gcm, ' Projections'))
  
  stackMeansList <- stack(meansList)
  projection(stackMeansList) <- getCRS('albersNA')
  bvMeans <- bioticVelocity(stackMeansList, times = seq(-21,0, by=1), onlyInSharedCells = T)
  bvMeansF <- bioticVelocity(stackMeansList, times = seq(-21,0, by=1), onlyInSharedCells = F)
  
  stackMaxList <- stack(maxList)
  projection(stackMaxList) <- getCRS('albersNA')
  bvMax <- bioticVelocity(stackMaxList, times = seq(-21,0, by=1), onlyInSharedCells = T)
  bvMaxF <- bioticVelocity(stackMaxList, times = seq(-21,0, by=1), onlyInSharedCells = F)
  
  # stackSumList <- stack(sumList)
  # projection(stackSumList) <- getCRS('albersNA')
  # bvSum <- bioticVelocity(stackSumList, times = seq(-21,0, by=1), onlyInSharedCells = T)
  # bvSumF <- bioticVelocity(stackSumList, times = seq(-21,0, by=1), onlyInSharedCells = F)
  # 
  pdf(file = paste0('./predictions/pdf/', gcm, '_bioticVelocity.pdf'), width = 11, height = 8.5)
  # no y axes limits
  bvMeans$time <- paste0(abs(bvMeans$timeFrom), '-', abs(bvMeans$timeTo), ' kybp')
  bvMeans$time <- factor(bvMeans$time, levels = rev(mixedsort(bvMeans$time)))
  
  bvMeansF$time <- paste0(abs(bvMeansF$timeFrom), '-', abs(bvMeansF$timeTo), ' kybp')
  bvMeansF$time <- factor(bvMeansF$time, levels = rev(mixedsort(bvMeansF$time)))
  
  ggplot(bvMeans, aes(time, centroidVelocity)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("centroid velocity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMeansF, aes(time, centroidVelocity)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means (shared cells = F)") + xlab("time period") + ylab("centroid velocity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  bvMax$time <- paste0(abs(bvMax$timeFrom), '-', abs(bvMax$timeTo), ' kybp')
  bvMax$time <- factor(bvMax$time, levels = rev(mixedsort(bvMax$time)))
  
  bvMaxF$time <- paste0(abs(bvMaxF$timeFrom), '-', abs(bvMaxF$timeTo), ' kybp')
  bvMaxF$time <- factor(bvMaxF$time, levels = rev(mixedsort(bvMaxF$time)))
  
  ggplot(bvMax, aes(time, centroidVelocity)) + geom_bar(stat = 'identity') + 
    ggtitle("Max") + xlab("time period") + ylab("centroid velocity") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMaxF, aes(time, centroidVelocity)) + geom_bar(stat = 'identity') + 
    ggtitle("Max (shared cells = F)") + xlab("time period") + ylab("centroid velocity") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # bvSum$time <- paste0(abs(bvSum$timeFrom), '-', abs(bvSum$timeTo), ' kybp')
  # bvSum$time <- factor(bvSum$time, levels = rev(mixedsort(bvSum$time)))
  # 
  # bvSumF$time <- paste0(abs(bvSumF$timeFrom), '-', abs(bvSumF$timeTo), ' kybp')
  # bvSumF$time <- factor(bvSumF$time, levels = rev(mixedsort(bvSumF$time)))
  # 
  # ggplot(bvSum, aes(time, centroidVelocity)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Sums") + xlab("time period") + ylab("centroid velocity") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # ggplot(bvSumF, aes(time, centroidVelocity)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Sums (shared cells = F)") + xlab("time period") + ylab("centroid velocity") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # all equal y axes limits
  # set limits
  minBV <- min(c(min(bvMeans$centroidVelocity), 
                 min(bvMax$centroidVelocity))) - 1000
  maxBV <- max(c(max(bvMeans$centroidVelocity), 
                 max(bvMax$centroidVelocity))) + 1000
  
  ggplot(bvMeans, aes(time, centroidVelocity)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("centroid velocity") +
    ylim(0, maxBV) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMax, aes(time, centroidVelocity)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Max") + xlab("time period") + ylab("centroid velocity") +
    ylim(0, maxBV) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ggplot(bvSum, aes(time, centroidVelocity)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Sums") + xlab("time period") + ylab("centroid velocity") +
  #   ylim(0, maxBV) + theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
  
  # NS Quant Velocity
  min05 <- min(c(min(bvMeans$nsQuantVelocity_quant0p05), 
               min(bvMax$nsQuantVelocity_quant0p05))) - 1000
  min95 <- min(c(min(bvMeans$nsQuantVelocity_quant0p95), 
                 min(bvMax$nsQuantVelocity_quant0p95))) - 1000
  max05 <- max(c(max(bvMeans$nsQuantVelocity_quant0p05), 
                max(bvMax$nsQuantVelocity_quant0p05))) + 1000
  max95 <- max(c(max(bvMeans$nsQuantVelocity_quant0p95), 
                max(bvMax$nsQuantVelocity_quant0p95))) + 1000
    
  pdf(file = paste0('./predictions/pdf/', gcm, '_nsQuantVelocity.pdf'), width = 11, height = 8.5)
  # no y axes limits
  ggplot(bvMeans, aes(time, nsQuantVelocity_quant0p05)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("NS_quant_05") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMeans, aes(time, nsQuantVelocity_quant0p95)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("NS_quant_95") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMax, aes(time, nsQuantVelocity_quant0p05)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("NS_quant_05") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMax, aes(time, nsQuantVelocity_quant0p95)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("NS_quant_95") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ggplot(bvSum, aes(time, nsQuantVelocity_quant0p05)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Means") + xlab("time period") + ylab("NS_quant_05") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # ggplot(bvSum, aes(time, nsQuantVelocity_quant0p95)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Means") + xlab("time period") + ylab("NS_quant_95") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # all equal y axes limits
  ggplot(bvMeans, aes(time, nsQuantVelocity_quant0p05)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("NS_quant_05") +
    ylim(min05, max05) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMeans, aes(time, nsQuantVelocity_quant0p95)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Means") + xlab("time period") + ylab("NS_quant_95") +
    ylim(min95, max95) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMax, aes(time, nsQuantVelocity_quant0p05)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Max") + xlab("time period") + ylab("NS_quant_05") +
    ylim(min05, max05) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(bvMax, aes(time, nsQuantVelocity_quant0p95)) + 
    geom_bar(stat = 'identity') + 
    ggtitle("Max") + xlab("time period") + ylab("NS_quant_95") +
    ylim(min95, max95) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ggplot(bvSum, aes(time, nsQuantVelocity_quant0p05)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Sums") + xlab("time period") + ylab("NS_quant_05") +
  #   ylim(min05, max05) + theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # ggplot(bvSum, aes(time, nsQuantVelocity_quant0p95)) + 
  #   geom_bar(stat = 'identity') + 
  #   ggtitle("Sums") + xlab("time period") + ylab("NS_quant_95") +
  #   ylim(min95, max95) + theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dev.off()
  
  save.image(file = paste0('./workspaces/06 - ', gcm, ' Projections'))
}

