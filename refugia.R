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
library(maps)
library(sf)
library(spatialEco)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
# setwd('/mnt/research/TIMBER/PVMvsENM')

ll <- c('longitude', 'latitude')

gcmList_ <- c('Beyer','Lorenz_ccsm', 'ecbilt')
pc <- 5
predictors <- c(paste0('pca', 1:pc))
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata',
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

# pollen <- readRDS('/Volumes/lj_mac_22/pollen/bv/FRAXINUS_bvs_n200_v4.1.RDS')
# pollen <- readRDS('/Volumes/lj_mac_22/pollen/polya-gamma-predictions_4.1_overdispersed-001.RDS')
# pollen_locs <- readRDS('/Volumes/lj_mac_22/pollen/pollen_locs_4.1.RDS')
# pollen_dat <- readRDS('/Volumes/lj_mac_22/pollen/pollen_dat_4.1.RDS')
# pollen_taxa <- readRDS('/Volumes/lj_mac_22/pollen/taxa_4.1.RDS')
# pollen_grid <- readRDS('/Volumes/lj_mac_22/pollen/grid_4.1.RDS')
# 
# # keep only fraxinus
# pollen_dat <- pollen_dat[1:465, 8, 1:22]
# pollen_dat <- as.data.frame(pollen_dat)
# View(pollen_dat)
# colnames(pollen_dat) <- c(1:22)
# 
# pollen <- pollen$pi[1:200, 1:3951, 8, 1:22]
# pollen <- as.data.frame(matrix(pollen, prod(dim(pollen)[1:2]), dim(pollen)[3]))
# pollen <- cbind(pollen_grid, pollen)
# 
# alb_proj <-  '+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
# wgs_proj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# 
# # assigning a crs to the pollen coordinates
# spdf <- SpatialPointsDataFrame(coords = pollen_grid[,c('x', 'y')], data = pollen_grid,
#                                proj4string = CRS(alb_proj))
# 
# # transforming to the lat long crs
# dat_transform <- spTransform(spdf, wgs_proj)
# dat <- coordinates(dat_transform)
# dat <- data.frame(dat_transform)
# dat <- dat[,1:4]
# colnames(dat)[3:4] <- ll
# 
# pollen <- merge(dat, pollen, by = c('x', 'y'))
# pollen_1 <- data.frame('longitude' = pollen$longitude, 'latitude' = pollen$latitude, 
#                        'abund' = pollen$V1)
# pollen_sp <- st_as_sf(pollen_1, coords = ll, crs = getCRS('wgs84'))
# plot(pollen_sp)
# 
# pollen_refuge <- pollen_1[pollen_1$abund >= 0.03,]
# pollen_refuge <- st_as_sf(pollen_refuge, coords = ll, crs = getCRS('wgs84'))
# plot(pollen_refuge)
# 
# pollen_sp <- SpatialPointsDataFrame(coords = pollen_1[,ll], data = pollen_1,
#                                proj4string = CRS(alb_proj))
# plot(as(pollen_sp, 'Spatial'))

cols <- c('gray83', '#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45', '#006d2c', '#00441b')
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- as(world, "Spatial")

for(gcm in gcmList_) {
  thresholds <- list()
  # pdf(file = paste0('./PDF_output/', gcm, '_refugia', '.pdf'), 
  #     width = 11, height = 8.5)
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
    # par(mfrow=c(1,2))
    # plot(out$simulationScale[[1]], main = paste0('Refugia\n', speciesAb_, '\n', gcm), 
    #      col = cols, axes = F)
    # map("world", add = T)
    plot(out$simulationScale[[2]], main = paste0('Refugia abundance\n', speciesAb_, '\n', gcm), 
         col = cols, axes = F, box = F)
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
  #      col = cols, axes = F)
  # map("world", add = T)
  plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = paste0('Refugia abundance\n', gcm), 
       col = cols, axes = F, box = F)
  # plot(raster::crop(sp::spTransform(world, CRS(projection(out$simulationScale[[2]]))), 
  #                   extent(out$simulationScale[[2]])),
  #      border = 'black', add = T)
  save.image(paste0('./workspaces/07 - Analyses, ', gcm, ' Refugia'))
}


pollen_threshold <- 0.03
# par(mfrow=c(1,2))
# fileName <- '/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred_iceMask.tif'
# pollenRast <- brick(fileName)
b <- pollenRast$Fraxinus_pollen_predictions_21kybp
thresholds <- seq(0.001, 1, by = 0.001)
kappa <- data.frame('threshold' = thresholds, 'kappa' = NA)

for(t in thresholds) {
  # refugia <- b >= t
  refugia <- b >= pollen_threshold
  refugiaId <- raster::clump(refugia, directions = 8, gaps = F)
  names(refugiaId) <- 'refugiaId'
  # plot(refugiaId, main = paste0(names(b),' ', title), axes = F)
  abund <- b * refugia
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
  k <- cohen.kappa(cbind(as.vector(as.matrix(abund)), as.vector(as.matrix(b))))
  kappa$kappa[which(kappa$threshold == t)] <- k$kappa
  
  # out <- list(
  #   simulationScale = raster::stack(refugiaId, abund),
  #   refugeCellNum = refugeCellNum,
  #   meanRefugeAbund = meanRefugeAbund
  # )

  # par(mfrow=c(1,2))
  # plot(out$simulationScale[[1]], main = paste0('Refugia\n', speciesAb_, '\n', gcm), 
  #      col = cols, axes = F)
  # map("world", add = T)
  plot(out$simulationScale[[2]], main = paste0('Refugia abundance\nPollen'),
       col = cols, axes = F, box = F)
  # map("world", add = T)
}

par(mfrow=c(2,2))
plot(out$simulationScale[[2]], main = paste0('Refugia abundance\nPollen'),
     col = cols, axes = F, box = F)
load('./workspaces/07 - Analyses, ecbilt Refugia')
plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = 'Refugia abundance\nECBilt', 
     col = cols, axes = F, box = F)
load('./workspaces/07 - Analyses, Lorenz_ccsm Refugia')
plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = 'Refugia abundance\nCCSM', 
     col = cols, axes = F, box = F)
load('./workspaces/07 - Analyses, Beyer Refugia')
plot(mask(out$simulationScale[[2]], pollenRast[[22]]), main = 'Refugia abundance\nHadAM3H', 
     col = cols, axes = F, box = F)
