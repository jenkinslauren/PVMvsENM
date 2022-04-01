rm(list = ls())
library(spatialEco)
library(enmSdm)
library(stats)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

gcmList_ <- c('Lorenz_ccsm', 'ecbilt', 'Beyer')
pc <- 5

# set constants
climYears <- seq(21000, 0, by=-1000)

studyRegionFileName <- './regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- as(world, "Spatial")

pollenRast <- brick('/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred.tif')

pdf(file = './PDF/rank_correlation.pdf', width = 11, height = 8.5)
par(mfrow = c(1,3), mar = c(2,1,5,1)+0.1)
for(a in 1:22) {
  for(gcm in gcmList_) {
    load(paste0('./workspaces/06 - ', gcm, ' Projections'))
    if(gcm == 'Lorenz_ccsm') g <- 'CCSM'
    if(gcm == 'ecbilt') g <- 'ECBilt'
    if(gcm == 'Beyer') g <- 'HadAM3H'
    t <- paste0('Pearson Correlation,\nPollen & ', g, '\n', climYears[a], ' ybp')
    thisRast <- raster::resample(meansList[[a]], pollenRast[[a]], method = 'bilinear')
    cor <- rasterCorrelation(thisRast, pollenRast[[a]], type = 'pearson')
    plot(cor, main = t, axes = F, box = F, legend.mar = 10)
    plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)), 
                      border = 'black', add = T)
  }
}
dev.off()

load(paste0('./workspaces/06 - ', 'Beyer', ' Projections'))
beyer <- meansList
load(paste0('./workspaces/06 - ', 'Lorenz_ccsm', ' Projections'))
ccsm <- meansList
load(paste0('./workspaces/06 - ', 'ecbilt', ' Projections'))
ecbilt <- meansList

corr_list <- c()
df <- data.frame(row.names = c('Pollen', 'HadAM3H', 'CCSM', 'ECBilt'))
for(g in 1:22) {
  st <- stack(pollenRast[[g]], beyer[[g]], ccsm[[g]], ecbilt[[g]])
  x <- layerStats(st, 'pearson', na.rm = T)
  corr_matrix <- as.data.frame(x$`pearson correlation coefficient`, 
                               row.names = c('Pollen', 'HadAM3H', 'CCSM', 'ECBilt'))
  df <- cbind(df, corr_matrix[1])
}

colnames(df) <- paste0(climYears, ' YBP')

