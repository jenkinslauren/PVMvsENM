rm(list = ls())
library(spatialEco)
library(enmSdm)
library(stats)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

gcmList_ <- c('Lorenz_ccsm', 'ecbilt', 'Beyer')
pc <- 5

# set constants
climYears <- seq(21000, 0, by=-1000)

studyRegionFileName <- './data_and_analyses/study_region/regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
studyRegionRasts <- brick(studyRegionFileName)

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- as(world, "Spatial")

pollenRast <- brick('/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred_iceMask.tif')
projection(pollenRast) <- getCRS('wgs84')
colors <- c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf',
            '#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')

pdf(file = './PDF_output/rank_correlation.pdf', width = 11, height = 8.5)
par(mfrow = c(1,3), mar = c(2,1,5,1)+0.1)
for(a in 1:22) {
  for(gcm in gcmList_) {
    load(paste0('./workspaces/06 - ', gcm, ' Projections'))
    if(gcm == 'Lorenz_ccsm') g <- 'CCSM'
    if(gcm == 'ecbilt') g <- 'ECBilt'
    if(gcm == 'Beyer') g <- 'HadAM3H'
    t <- paste0('Pearson Correlation,\nPollen & ', g, '\n', climYears[a], ' ybp')
    projection(meansList) <- getCRS('wgs84')
    thisRast <- raster::resample(meansList[[a]], pollenRast[[a]], method = 'bilinear')
    
    cor <- rasterCorrelation(thisRast, pollenRast[[a]], type = 'pearson', s = 5)
    plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
    plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)), 
         border = 'black', add = T)
  }
  
  load('./workspaces/06 - Beyer Projections')
  beyer <- 'HadAM3H'
  bMeans <- meansList
  load('./workspaces/06 - ecbilt Projections')
  ecbilt <- 'ECBilt'
  eMeans <- meansList
  load('./workspaces/06 - Lorenz_ccsm Projections')
  ccsm <- 'CCSM'
  cMeans <- meansList

  t <- paste0('Pearson Correlation,\n', beyer, '& ', ecbilt, '\n', climYears[a], ' ybp')
  cor <- rasterCorrelation(bMeans[[a]], eMeans[[a]], type = 'pearson')
  plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
  plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)),
       border = 'black', add = T)

  t <- paste0('Pearson Correlation,\n', beyer, '& ', ccsm, '\n', climYears[a], ' ybp')
  cor <- rasterCorrelation(bMeans[[a]], cMeans[[a]], type = 'pearson')
  plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
  plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)),
       border = 'black', add = T)

  t <- paste0('Pearson Correlation,\n', ccsm, '& ', ecbilt, '\n', climYears[a], ' ybp')
  cor <- rasterCorrelation(cMeans[[a]], eMeans[[a]], type = 'pearson')
  plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
  plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)),
       border = 'black', add = T)
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

# do this for all GCM + yrs
corDf <- data.frame(climYears)
for (gcm in gcmList_) {
  correlations <- list()
  for(a in 1:22) {
    load(paste0('./workspaces/06 - ', gcm, ' Projections'))
    if(gcm == 'Lorenz_ccsm') g <- 'CCSM'
    if(gcm == 'ecbilt') g <- 'ECBilt'
    if(gcm == 'Beyer') g <- 'HadAM3H'
    thisRast <- raster::resample(meansList[[a]], pollenRast[[a]], method = 'bilinear')
    
    # c <- cor(as.vector(as.matrix(pollenRast[[a]])), as.vector(as.matrix(thisRast)), 
    #          method = c('pearson'), use = 'na.or.complete')
    cNiche <- enmSdm::compareNiches(as.vector(as.matrix(pollenRast[[a]])), 
                                    as.vector(as.matrix(thisRast)), na.rm = T)
    correlations <- append(correlations, cNiche['cor'])
  }
  names(correlations) <- climYears
  corDf <- cbind(corDf, t(data.frame(correlations)))
  colnames(corDf)[ncol(corDf)] <- paste0('Pollen & ', g)
}

load('./workspaces/06 - Beyer Projections')
beyer <- 'HadAM3H'
bMeans <- meansList
load('./workspaces/06 - ecbilt Projections')
ecbilt <- 'ECBilt'
eMeans <- meansList
load('./workspaces/06 - Lorenz_ccsm Projections')
ccsm <- 'CCSM'
cMeans <- meansList

correlations <- list()
for(a in 1:22) {
  # c <- cor(as.vector(as.matrix(pollenRast[[a]])), as.vector(as.matrix(thisRast)), 
  #          method = c('pearson'), use = 'na.or.complete')
  cNiche <- enmSdm::compareNiches(as.vector(as.matrix(bMeans[[a]])), 
                                  as.vector(as.matrix(cMeans[[a]])), na.rm = T)
  correlations <- append(correlations, cNiche['cor'])
}

names(correlations) <- climYears
corDf <- cbind(corDf, t(data.frame(correlations)))
colnames(corDf)[ncol(corDf)] <- paste0(beyer, ' & ', ccsm)

correlations <- list()
for(a in 1:22) {
  cNiche <- enmSdm::compareNiches(as.vector(as.matrix(bMeans[[a]])), 
                                  as.vector(as.matrix(eMeans[[a]])), na.rm = T)
  correlations <- append(correlations, cNiche['cor'])
}

names(correlations) <- climYears
corDf <- cbind(corDf, t(data.frame(correlations)))
colnames(corDf)[ncol(corDf)] <- paste0(beyer, ' & ', ecbilt)

correlations <- list()
for(a in 1:22) {
  cNiche <- enmSdm::compareNiches(as.vector(as.matrix(eMeans[[a]])), 
                                  as.vector(as.matrix(cMeans[[a]])), na.rm = T)
  correlations <- append(correlations, cNiche['cor'])
}

names(correlations) <- climYears
corDf <- cbind(corDf, t(data.frame(correlations)))
colnames(corDf)[ncol(corDf)] <- paste0(ecbilt, ' & ', ccsm)

corDf$climYears <- paste0(corDf$climYears/1000, ' kybp')
colnames(corDf)[1] <- 'Time'
corDf$Time <- factor(corDf$Time, 
                          levels = rev(mixedsort(corDf$Time)))

corDf_plot <- cbind(corDf[1], utils::stack(corDf[2:7]))
corDf_plot_pollen <- cbind(corDf[1], utils::stack(corDf[2:4]))
corDf_plot_gcm <- cbind(corDf[1], utils::stack(corDf[5:ncol(corDf)]))

colnames(corDf_plot) <- c("Time", "cor", "Models")
colnames(corDf_plot_pollen) <- c("Time", "cor", "Models")
colnames(corDf_plot_gcm) <- c("Time", "cor", "Models")

pdf(file = './PDF_output/rank_correlation_line.pdf', width = 11, height = 8.5)
ggplot(corDf_plot, aes(Time, cor, color = Models, group = Models)) + 
  geom_point() + geom_line() + 
  theme_classic() + 
  labs(title = "Rank Correlation", x = "Time", 
       y = "cor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(corDf_plot_pollen, aes(Time, cor, color = Models, group = Models)) + 
  geom_point() + geom_line() + 
  theme_classic() + 
  labs(title = "Rank Correlation", x = "Time", 
       y = "cor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(corDf_plot_gcm, aes(Time, cor, color = Models, group = Models)) + 
  geom_point() + geom_line() + 
  theme_classic() + 
  labs(title = "Rank Correlation", x = "Time", 
       y = "cor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

###################################################################
###### Individual Species
###################################################################
gcm <- 'Lorenz_ccsm'
fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)
load(paste0('./workspaces/06 - ', gcm, ' Projections'))
if(gcm == 'Lorenz_ccsm') g <- 'CCSM'
if(gcm == 'ecbilt') g <- 'ECBilt'
if(gcm == 'Beyer') g <- 'HadAM3H'

for(f in fileName) {
  s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
  speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
  load(paste0('./Models/Maxent/all_model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  b <- stack(f)
  names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
  projection(b) <- getCRS('wgs84')
  
  for(a in 1:22) {
    t <- paste0('Pearson Correlation,\nPollen & ', g, '\n(', speciesAb_,
                ')\n', climYears[a], ' ybp')
    thisRast <- raster::resample(b[[a]], pollenRast[[a]], method = 'bilinear')
    
    cor <- rasterCorrelation(thisRast, pollenRast[[a]], type = 'pearson', s = 5)
    plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
    plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)), 
         border = 'black', add = T)
  }
}

for(a in 1:22) {
  for(gcm in gcmList_) {
    load(paste0('./workspaces/06 - ', gcm, ' Projections'))
    if(gcm == 'Lorenz_ccsm') g <- 'CCSM'
    if(gcm == 'ecbilt') g <- 'ECBilt'
    if(gcm == 'Beyer') g <- 'HadAM3H'
    t <- paste0('Pearson Correlation,\nPollen & ', g, '\n', climYears[a], ' ybp')
    projection(meansList) <- getCRS('wgs84')
    thisRast <- raster::resample(meansList[[a]], pollenRast[[a]], method = 'bilinear')
    
    cor <- rasterCorrelation(thisRast, pollenRast[[a]], type = 'pearson', s = 5)
    plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
    plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)), 
         border = 'black', add = T)
  }
  
  load('./workspaces/06 - Beyer Projections')
  beyer <- 'HadAM3H'
  bMeans <- meansList
  load('./workspaces/06 - ecbilt Projections')
  ecbilt <- 'ECBilt'
  eMeans <- meansList
  load('./workspaces/06 - Lorenz_ccsm Projections')
  ccsm <- 'CCSM'
  cMeans <- meansList
  
  t <- paste0('Pearson Correlation,\n', beyer, '& ', ecbilt, '\n', climYears[a], ' ybp')
  cor <- rasterCorrelation(bMeans[[a]], eMeans[[a]], type = 'pearson')
  plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
  plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)),
       border = 'black', add = T)
  
  t <- paste0('Pearson Correlation,\n', beyer, '& ', ccsm, '\n', climYears[a], ' ybp')
  cor <- rasterCorrelation(bMeans[[a]], cMeans[[a]], type = 'pearson')
  plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
  plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)),
       border = 'black', add = T)
  
  t <- paste0('Pearson Correlation,\n', ccsm, '& ', ecbilt, '\n', climYears[a], ' ybp')
  cor <- rasterCorrelation(cMeans[[a]], eMeans[[a]], type = 'pearson')
  plot(cor, main = t, axes = F, box = F, legend.mar = 10, col = colors)
  plot(raster::crop(sp::spTransform(world, CRS(projection(cor))), extent(cor)),
       border = 'black', add = T)
}



###################################################################
###### JPEGS
###################################################################
jpeg(file = '/Users/laurenjenkins/Downloads/rank_corr_all.jpeg',
     width = 21, height = 21, units = 'cm', res = 300)
ggplot(corDf_plot, aes(Time, cor, color = Models, group = Models)) + 
  geom_point() + geom_line() + 
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) + 
  geom_hline(yintercept = 0.5 , linetype="dashed", color = "red", size=0.5) +
  theme_classic() + 
  labs(title = "Rank Correlation", x = "Time", 
       y = "cor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

jpeg(file = '/Users/laurenjenkins/Downloads/rank_corr_pollen.jpeg',
     width = 21, height = 21, units = 'cm', res = 300)
ggplot(corDf_plot_pollen, aes(Time, cor, color = Models, group = Models)) + 
  geom_point() + geom_line() + 
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) + 
  geom_hline(yintercept = 0.5 , linetype="dashed", color = "red", size=0.5) +
  theme_minimal() + 
  labs(title = "Rank Correlation", x = "Time", 
       y = "cor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

jpeg(file = '/Users/laurenjenkins/Downloads/rank_corr_gcm.jpeg',
     width = 21, height = 21, units = 'cm', res = 300)
ggplot(corDf_plot_gcm, aes(Time, cor, color = Models, group = Models)) + 
  geom_point() + geom_line() + 
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) + 
  geom_hline(yintercept = 0.5 , linetype="dashed", color = "red", size=0.5) +
  theme_minimal() + 
  labs(title = "Rank Correlation", x = "Time", 
       y = "cor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
