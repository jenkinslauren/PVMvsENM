rm(list = ls())

library(dismo)
library(sp)
library(enmSdm)
library(geosphere)
library(raster)

library(maps)
library(sf)
library(spatialEco)

library(xlsx)
library(data.table)
library(dplyr)
library(viridis)

genus <- 'fraxinus'
speciesList <- paste0('Fraxinus ', 
                      c('americana', 'caroliniana', 'cuspidata',
                        'greggii', 'nigra', 'pennsylvanica', 
                        'profunda', 'quadrangulata'))

setwd(paste0('/Volumes/lj_mac_22/MOBOT/by_genus/', genus))
# setwd('/mnt/research/TIMBER/PVMvsENM')

ll <- c('longitude', 'latitude')
gcmList <- c('hadley','ccsm', 'ecbilt')

pc <- 5
predictors <- c(paste0('pca', 1:pc))

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- as(world, "Spatial")

colors <- c('gray83', '#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45', '#006d2c', '#00441b')

# function for masking pollen model by land & ice # 
getPollen <- function(times) {
  if(!file.exists(paste0('/Volumes/lj_mac_22/pollen/predictions-', 
                         toupper(genus), '_meanpred_iceMask.tif'))) {
    
    maps <- stack(paste0('/Volumes/lj_mac_22/pollen/predictions-', 
                         toupper(genus), '_meanpred.tif'))
    
    ### mask by glaciers and available land ###
    daltonAges <- read.csv('/Volumes/lj_mac_22/Dalton et al 2020 QSR Ice Layers/Dalton et al 2020 QSR Dates from Shapefile Names.csv')
    
    # mask by land (for visualization) #
    for (countTime in seq_along(times)) {
      time <- times[countTime]
      
      # land mask
      land <- raster(paste0('/Volumes/lj_mac_22/MOBOT/by_genus/env_data/ccsm/tifs/', 
                            -1 * time, 'BP/an_avg_TMAX.tif'))
      # land <- land * 0 + 1
      land <- projectRaster(land, maps)
      land <- land * 0 + 1
      maps[[countTime]] <- maps[[countTime]] * land
      
    }
    
    ### mask by ice (for calculating BV) ### 
    mapsMasked <- maps
    
    for (countTime in seq_along(times)) {
      time <- times[countTime]
      
      # ice mask
      closestDalton <- which.min(abs(-1000 * daltonAges$calKiloYear - time))
      
      load(paste0('/Volumes/lj_mac_22/Dalton et al 2020 QSR Ice Layers/RDA Files/daltonEtAl2020_', 
                  sprintf('%.2f', daltonAges$calKiloYear[closestDalton]), '_kiloCalYBP.rda'))
      daltonIce <- sp::spTransform(daltonIce, getCRS('albersNA', TRUE))
      
      daltonIce <- rasterize(daltonIce, maps)
      daltonIceMask <- calc(daltonIce, fun=function(x) ifelse(is.na(x), 1, NA))
      mapsMasked[[countTime]] <- mapsMasked[[countTime]] * daltonIceMask
      
    }
    
    writeRaster(stack(maps), 
                paste0('/Volumes/lj_mac_22/pollen/predictions-', 
                       toupper(genus), '_meanpred_landMask.tif'), 
                format = 'GTiff', overwrite = T)
    
    writeRaster(stack(mapsMasked), 
                paste0('/Volumes/lj_mac_22/pollen/predictions-', 
                       toupper(genus), '_meanpred_iceMask.tif'), 
                format = 'GTiff', overwrite = T)
  } 
  
  return(stack(paste0('/Volumes/lj_mac_22/pollen/predictions-', 
                      toupper(genus), 
                      '_meanpred_iceMask.tif')))
}

for(gcm in gcmList) {
  thresholds <- list()
  
  for(a in 1:length(speciesList)) {
    sp <- speciesList[a]
    speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
    rangeName <- paste0('littleRange_', speciesAb)
    
    print(paste0("Species = ", species))
    
    folderName <- paste0('./models/predictions/', speciesAb_,
                         '/Model Evaluation - Random K-Folds - ', gcm)
    if(!dir.exists(folderName)) dir.create(folderName)
    
    # load bg sites # 
    load('/Volumes/lj_mac_22/MOBOT/by_genus/background_sites/Random Background Sites across Study Region.Rdata')
    
    # load model #
    load(paste0('./models/', speciesAb_, '_Maxent_PC', pc, '_GCM_', gcm, '.rData'))
   
     # load range map #
    load(paste0('./range_maps/', rangeName, '.Rdata'))
    
    # load(paste0('./Models/Maxent/', gcm, '_evals.Rdata'))
    
    # load k-folds for presences and background sites from model evaluations
    t <- list()
    for(i in 1:5) {
      print(paste0('K-fold ', i, ':'))
      
      load(paste0(folderName, '/Model ', i, '.Rdata'))
      
      temp <- enmSdm::thresholdWeighted(predPres, predBg, na.rm = T)
      t <- append(t, temp['msss'])
    }
    
    thresholds[[a]] <- t
    
  }
  
  # view thresholds as data frame & organize columns #
  thresholds <- data.frame(t(rbindlist(thresholds, fill = T)))
  colnames(thresholds) <- speciesList
  rownames(thresholds) <- paste("K-fold", 1:5) 
  thresholds[thresholds == 0] <- NA
  thresholds <- rbind(thresholds, mean = summarize_all(thresholds, mean, na.rm = T))
  
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
    
    # pdf(file = paste0('./pdf/', gcm, '_refugia', '.pdf'), width = 11, height = 8.5)
    
    # par(mfrow=c(1,2))
    # plot(out$simulationScale[[1]], main = paste0('Refugia\n', speciesAb_, '\n', gcm), 
    #      col = colors, axes = F)
    # map("world", add = T)
    plot(out$simulationScale[[2]], main = paste0('Refugia abundance\n', speciesAb_, '\n', gcm), 
         col = colors, axes = F, box = F)
    maps::map("world", add = T)
  }
  # dev.off()
}

fileName <- '/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred_iceMask.tif'
pollenRast <- stack(fileName)
names(pollenRast) <- c(paste0("Fraxinus_pollen_predictions_", 0:21, "kybp"))
pollenRast <- unstack(pollenRast)
pollenRast <- stack(rev(pollenRast))

for(gcm in gcmList_) {
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
    load('./Background Sites/Random Background Sites across Study Region.Rdata')
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
      t <- append(t, temp['msss'])
    }
    thresholds[[a]] <- t
    
  }
  
  thresholds <- data.frame(t(rbindlist(thresholds, fill = T)))
  colnames(thresholds) <- speciesList
  rownames(thresholds) <- paste("K-fold", 1:5) 
  thresholds[thresholds == 0] <- NA
  
  thresholds <- rbind(thresholds, mean = summarize_all(thresholds, mean, na.rm = T))
  
  threshold <- as.numeric(rowMeans(thresholds)['mean'])
  
  load(paste0('./workspaces/06 - ', gcm, ' Projections'))
  b <- stack(meansList)
  
  b <- b[[1]]
  names(b) <- paste0(21, ' Kybp')
  
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
  # par(mfrow=c(1,2))
  # plot(out$simulationScale[[1]], main = paste0('Refugia\n', speciesAb_, '\n', gcm), 
  #      col = colors, axes = F)
  # map("world", add = T)
  plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = paste0('Refugia abundance\n', gcm), 
       col = colors, axes = F, box = F)
  # plot(raster::crop(sp::spTransform(world, CRS(projection(out$simulationScale[[2]]))), 
  #                   extent(out$simulationScale[[2]])),
  #      border = 'black', add = T)
  save.image(paste0('./workspaces/07 - Analyses, ', gcm, ' Refugia'))
}

library(psych)
library(vegan)

jaccard_fun <- function(x, y) {
  intersection = length(intersect(x, y))
  union = length(x) + length(y) - intersection
  return (intersection/union)
}

pollen_threshold <- 0.03
# par(mfrow=c(1,2))
fileName <- '/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred_iceMask.tif'
pollenRast <- brick(fileName)
pollenRast <- pollenRast$predictions.FRAXINUS_meanpred_iceMask.22
gcmRast <- out$simulationScale[[2]]
thresholds <- seq(0.001, 1, by = 0.001)
kappa <- data.frame('threshold' = thresholds, 'kappa' = NA)
jaccard <- data.frame('threshold' = thresholds, 'j' = NA)

pdf(file = '/Users/laurenjenkins/Downloads/pollen_refugia.pdf', width = 11.5, height = 8.5)
for(t in thresholds) {
  refugia <- pollenRast >= t
  refugiaId <- raster::clump(refugia, directions = 4, gaps = F)
  names(refugiaId) <- 'refugiaId'
  plot(refugiaId, main = paste0(names(pollenRast),' ', title), axes = F)
  abund <- pollenRast * refugia
  names(abund) <- 'refugiaAbund'
  
  # nrows <- nrow(b)
  # ncols <- ncol(b)
  # ncells <- raster::ncell(b)
  # 
  # v <- rep(seq(nrows * (ncols - 1) - 1, 1, by=-ncols), each=ncols) + 0:(ncols - 1)
  # cellNum <- matrix(v, nrow=nrows, ncol=ncols, byrow=TRUE)
  # cellNum <- raster::raster(cellNum, template=b)
  # cellNum <- as.vector(cellNum)
  # simRefugiaBinary <- as.vector(refugia)
  # refugeCellNum <- cellNum[simRefugiaBinary]
  # if (any(is.na(refugeCellNum))) refugeCellNum <- refugeCellNum[!is.na(refugeCellNum)]
  # 
  # # mean refuge abundance
  # meanRefugeAbund <- raster::cellStats(abund, 'sum') / length(refugeCellNum)

  # k <- raster.change(abund, b, stat = c("kappa"), mask = T)
  
  # generate a cohen's kappa value for pollen refugia vs GCM refugia
  k <- cohen.kappa(cbind(as.vector(as.matrix(gcmRast)), as.vector(as.matrix(abund))))
  kappa$kappa[which(kappa$threshold == t)] <- k$kappa
  
  j <- jaccard_fun(as.vector(as.matrix(gcmRast)), as.vector(as.matrix(abund)))
  jaccard$j[which(jaccard$threshold == t)] <- j
  
  # out <- list(
  #   simulationScale = raster::stack(refugiaId, abund),
  #   refugeCellNum = refugeCellNum,
  #   meanRefugeAbund = meanRefugeAbund
  # )
  # write.xlsx(kappa, file = './pollen_refugia_thresholds.xlsx', sheetName = paste0(gcm, '_kappa'),
  #            append = T, row.names = F)
  # write.xlsx(jaccard, file = './pollen_refugia_thresholds.xlsx', sheetName = paste0(gcm, '_jaccard'),
  #            append = T, row.names = F)
}
dev.off()

par(mfrow=c(2,2))
plot(out$simulationScale[[2]], main = paste0('Refugia abundance\nPollen'),
     col = colors, axes = F, box = F)
load('./workspaces/07 - Analyses, ecbilt Refugia')
plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = 'Refugia abundance\nECBilt', 
     col = colors, axes = F, box = F)
load('./workspaces/07 - Analyses, Lorenz_ccsm Refugia')
plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = 'Refugia abundance\nCCSM', 
     col = colors, axes = F, box = F)
load('./workspaces/07 - Analyses, Beyer Refugia')
plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = 'Refugia abundance\nHadAM3H', 
     col = colors, axes = F, box = F)

# plot results of pollen thresholds experiment
# X = threshold, Y = Jaccard value
library(readxl)
gcm <- 'Lorenz_ccsm'
t <- read_excel('./pollen_refugia_thresholds.xlsx', sheet = paste0(gcm, "_jaccard"))

jaccard_ccsm <- ggplot(t, aes(threshold, j)) + 
  geom_point(size = 0.3) + theme_classic() + 
  labs(title = "Pollen threshold\ncompared to CCSM", x = "Threshold", 
       y = "Jaccard Similarity")

gcm <- 'ecbilt'
t <- read_excel('./pollen_refugia_thresholds.xlsx', sheet = paste0(gcm, "_jaccard"))

jaccard_ecbilt <- ggplot(t, aes(threshold, j)) + 
  geom_point(size = 0.3) + theme_classic() + 
  labs(title = "Pollen thresholds\ncompared to ECBilt", x = "Threshold", 
       y = "Jaccard Similarity")

gcm <- 'Beyer'
t <- read_excel('./pollen_refugia_thresholds.xlsx', sheet = paste0(gcm, "_jaccard"))

jaccard_hadley <- ggplot(t, aes(threshold, j)) + 
  geom_point(size = 0.3) + theme_classic() + 
  labs(title = "Pollen threshold\ncompared to HadAM3H", x = "Threshold", 
       y = "Jaccard Similarity")

library(cowplot)
pdf('./pollen_Jaccard.pdf', width = 11, height = 8.5)
plot_grid(jaccard_ccsm, jaccard_ecbilt, jaccard_hadley)
dev.off()

