rm(list=ls())
library(enmSdm)
library(gtools)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tidyverse)
library(raster)
library(utils)
library(ggrepel)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

gcm <- 'pollen'
fileName <- '/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred.tif'
pollenRast <- brick(fileName)

pollenStack <- stack(fileName)
projection(pollenStack) <- getCRS('albersNA')
# note: there is no difference in output between onlyInSharedCells = T & F
bvPollen <- bioticVelocity(pollenStack, times = seq(-21000,0, by=1000), onlyInSharedCells = T)

bvPollen$time <- paste0(abs(bvPollen$timeFrom)/1000, '-', abs(bvPollen$timeTo)/1000, ' kybp')
bvPollen$time <- factor(bvPollen$time, levels = rev(mixedsort(bvPollen$time)))
# bvPollen$time <- as.factor(bvPollen$time)

p <- ggplot(bvPollen, aes(x = time, y = centroidVelocity, group = 1)) +
  geom_line() + geom_point() +
  theme_classic() +
  ggtitle("Pollen") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

load(paste0('./workspaces/06 - Beyer Projections'))
# Beyer Means
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvBeyerMeans <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = T)

bvBeyerMeans$time <- paste0(abs(bvBeyerMeans$timeFrom)/1000, '-', abs(bvBeyerMeans$timeTo)/1000, ' kybp')
bvBeyerMeans$time <- factor(bvBeyerMeans$time, levels = rev(mixedsort(bvBeyerMeans$time)))

bMean <- ggplot(bvBeyerMeans, aes(time, centroidVelocity, group = 1)) + 
  geom_point() + geom_line() +
  ggtitle("Beyer Means") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Beyer Max
stackMaxList <- stack(maxList)
projection(stackMaxList) <- getCRS('albersNA')
bvBeyerMax <- bioticVelocity(stackMaxList, times = seq(-21000,0, by = 1000), onlyInSharedCells = T) 

bvBeyerMax$time <- paste0(abs(bvBeyerMax$timeFrom)/1000, '-', abs(bvBeyerMax$timeTo)/1000, ' kybp')
bvBeyerMax$time <- factor(bvBeyerMax$time, levels = rev(mixedsort(bvBeyerMax$time)))

bMax <- ggplot(bvBeyerMax, aes(time, centroidVelocity, group = 1)) + 
  geom_point() + geom_line() +
  ggtitle("Beyer Max") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################################################################
######### CCSM
###############################################################################

load(paste0('./workspaces/06 - Lorenz_ccsm Projections'))
# Lorenz_ccsm Means
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvCCSMMean <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = T)

bvCCSMMean$time <- paste0(abs(bvCCSMMean$timeFrom)/1000, '-', abs(bvCCSMMean$timeTo)/1000, ' kybp')
bvCCSMMean$time <- factor(bvCCSMMean$time, levels = rev(mixedsort(bvCCSMMean$time)))

cMean <- ggplot(bvCCSMMean, aes(time, centroidVelocity, group = 1)) + 
  geom_point() + geom_line() +
  ggtitle("CCSM Means") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Lorenz_ccsm Max
stackMaxList <- stack(maxList)
projection(stackMaxList) <- getCRS('albersNA')
bvCCSMMax <- bioticVelocity(stackMaxList, times = seq(-21000,0, by = 1000), onlyInSharedCells = T) 

bvCCSMMax$time <- paste0(abs(bvCCSMMax$timeFrom)/1000, '-', abs(bvCCSMMax$timeTo)/1000, ' kybp')
bvCCSMMax$time <- factor(bvCCSMMax$time, levels = rev(mixedsort(bvCCSMMax$time)))

cMax <- ggplot(bvMax, aes(time, centroidVelocity, group = 1)) + 
  geom_point() + geom_line() +
  ggtitle("CCSM Max") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################################################################
######### ECBILT
###############################################################################

load(paste0('./workspaces/06 - ecbilt Projections'))
# ECBilt Means
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvECBiltMean <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = T)

bvECBiltMean$time <- paste0(abs(bvECBiltMean$timeFrom)/1000, '-', abs(bvECBiltMean$timeTo)/1000, ' kybp')
bvECBiltMean$time <- factor(bvECBiltMean$time, levels = rev(mixedsort(bvECBiltMean$time)))

eMean <- ggplot(bvECBiltMean, aes(time, centroidVelocity, group = 1)) + 
  geom_point() + geom_line() +
  ggtitle("ECBilt Means") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ECBilt Max
stackMaxList <- stack(maxList)
projection(stackMaxList) <- getCRS('albersNA')
bvECBiltMax <- bioticVelocity(stackMaxList, times = seq(-21000,0, by = 1000), onlyInSharedCells = T) 

bvECBiltMax$time <- paste0(abs(bvECBiltMax$timeFrom)/1000, '-', abs(bvECBiltMax$timeTo)/1000, ' kybp')
bvECBiltMax$time <- factor(bvECBiltMax$time, levels = rev(mixedsort(bvECBiltMax$time)))

eMax <- ggplot(bvECBiltMax, aes(time, centroidVelocity, group = 1)) + 
  geom_point() + geom_line() +
  ggtitle("ECBilt Max") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# bv <- data.frame(bvPollen$time, 
#                  bvPollen$centroidVelocity, 
#                  bvBeyerMax$centroidVelocity, 
#                  bvBeyerMeans$centroidVelocity,
#                  bvCCSMMax$centroidVelocity,
#                  bvCCSMMean$centroidVelocity,
#                  bvECBiltMax$centroidVelocity,
#                  bvECBiltMean$centroidVelocity)

bv <- data.frame(bvPollen$time, 
                 bvPollen$centroidVelocity,
                 bvBeyerMeans$centroidVelocity,
                 bvCCSMMean$centroidVelocity,
                 bvECBiltMean$centroidVelocity)
colnames(bv) <- c("Time", "Pollen","Beyer", "CCSM", "ECBilt")
bv <- cbind(bv[1], utils::stack(bv[2:5]))
colnames(bv) <- c("Time", "centroidVelocity", "climateSource")

###############################################################################
###############################################################################
gcm <- 'ecbilt'
fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)
bv_all <- data.frame(bvECBiltMean$time, bvECBiltMean$centroidVelocity)
colnames(bv_all) <- c('Time', 'Entire Genus')

for(f in fileName) {
  s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
  speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
  load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  b <- stack(f)
  names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
  
  projection(b) <- getCRS('albersNA')
  sp_bv <- bioticVelocity(b, times = seq(-21000,0, by = 1000), onlyInSharedCells = T) 
  bv_all <- cbind(bv_all, sp_bv$centroidVelocity)
  colnames(bv_all)[ncol(bv_all)] <- speciesAb_
}

bv_all_sp <- cbind(bv_all[1], utils::stack(bv_all[3:ncol(bv_all)]))
colnames(bv_all_sp) <- c("Time", "centroidVelocity", "Species")

# all species
e_all_sp <- ggplot(bv_all_sp, aes(x = Time, y = centroidVelocity, group = Species, 
                     color = Species)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(vjust = -1)) +
  labs(title = "ECBilt Biotic Velocity\n", x = "Time Period", 
       y = "Centroid Velocity (m/yr)", color = "Species")

bv_select_sp <- data.frame(bv_all[3:4], bv_all[7:9])
bv_select_sp <- cbind(bv_all[1], utils::stack(bv_select_sp))
colnames(bv_select_sp) <- c("Time", "centroidVelocity", "Species")

# select species (best models)
e_select_sp <- ggplot(bv_select_sp, aes(x = Time, y = centroidVelocity, group = Species, 
                      color = Species)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(vjust = -1)) +
  labs(title = "ECBilt Biotic Velocity\n", subtitle = "Select Species Only", 
       x = "Time Period", y = "Centroid Velocity (m/yr)", color = "Species")

gcm <- 'Lorenz_ccsm'
fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)
bv_all <- data.frame(bvECBiltMean$time, bvECBiltMean$centroidVelocity)
colnames(bv_all) <- c('Time', 'Entire Genus')

for(f in fileName) {
  s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
  speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
  load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  b <- stack(f)
  names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
  
  projection(b) <- getCRS('albersNA')
  sp_bv <- bioticVelocity(b, times = seq(-21000,0, by = 1000), onlyInSharedCells = T) 
  bv_all <- cbind(bv_all, sp_bv$centroidVelocity)
  colnames(bv_all)[ncol(bv_all)] <- speciesAb_
}

bv_all_sp <- cbind(bv_all[1], utils::stack(bv_all[3:ncol(bv_all)]))
colnames(bv_all_sp) <- c("Time", "centroidVelocity", "Species")

# all species
c_all_sp <- ggplot(bv_all_sp, aes(x = Time, y = centroidVelocity, group = Species, 
                                  color = Species)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(vjust = -1)) +
  labs(title = "CCSM Biotic Velocity\n", x = "Time Period", 
       y = "Centroid Velocity (m/yr)", color = "Species")

bv_select_sp <- data.frame(bv_all[3:4], bv_all[7:9])
bv_select_sp <- cbind(bv_all[1], utils::stack(bv_select_sp))
colnames(bv_select_sp) <- c("Time", "centroidVelocity", "Species")

# select species (best models)
c_select_sp <- ggplot(bv_select_sp, aes(x = Time, y = centroidVelocity, group = Species, 
                                        color = Species)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(vjust = -1)) +
  labs(title = "CCSM Biotic Velocity\n", subtitle = "Select Species Only", 
       x = "Time Period", y = "Centroid Velocity (m/yr)", color = "Species")

gcm <- 'Beyer'
fileName <- list.files(path = paste0('./predictions/', gcm),
                       pattern = paste0('PC', pc,'.tif'),
                       full.names = T)
bv_all <- data.frame(bvECBiltMean$time, bvECBiltMean$centroidVelocity)
colnames(bv_all) <- c('Time', 'Entire Genus')

for(f in fileName) {
  s <- gsub('\\..*', '', gsub('\\./predictions/*', '', f))
  speciesAb_ <-  gsub('\\_GCM.*', '', gsub(paste0('\\./predictions/', gcm, '/*'), '', f))
  load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
              '_PC', pc, '.rData'))
  b <- stack(f)
  names(b) <- c(paste0(seq(21000, 0, by = -1000), ' ybp'))
  
  projection(b) <- getCRS('albersNA')
  sp_bv <- bioticVelocity(b, times = seq(-21000,0, by = 1000), onlyInSharedCells = T) 
  bv_all <- cbind(bv_all, sp_bv$centroidVelocity)
  colnames(bv_all)[ncol(bv_all)] <- speciesAb_
}

bv_all_sp <- cbind(bv_all[1], utils::stack(bv_all[3:ncol(bv_all)]))
colnames(bv_all_sp) <- c("Time", "centroidVelocity", "Species")

# all species
b_all_sp <- ggplot(bv_all_sp, aes(x = Time, y = centroidVelocity, group = Species, 
                                  color = Species)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(vjust = -1)) +
  labs(title = "HadAM3H Biotic Velocity\n", x = "Time Period", 
       y = "Centroid Velocity (m/yr)", color = "Species")

bv_select_sp <- data.frame(bv_all[3:4], bv_all[7:9])
bv_select_sp <- cbind(bv_all[1], utils::stack(bv_select_sp))
colnames(bv_select_sp) <- c("Time", "centroidVelocity", "Species")

# select species (best models)
b_select_sp <- ggplot(bv_select_sp, aes(x = Time, y = centroidVelocity, group = Species, 
                                        color = Species)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(vjust = -1)) +
  labs(title = "HadAM3H Biotic Velocity\n", subtitle = "Select Species Only", 
       x = "Time Period", y = "Centroid Velocity (m/yr)", color = "Species")

pdf(file = './PDF/bv_plots.pdf', height = 8.5, width = 11)

plot_grid(p, bMean, bMax, cMean, cMax, eMean, eMax)

ggplot(bv, aes(x = Time, y = centroidVelocity, group = climateSource, 
               color = climateSource)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Biotic Velocity\n", x = "Time Period", 
       y = "Centroid Velocity (m/yr)", color = "GCM") +
  scale_color_manual(values = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3'))

plot_grid(e_all_sp, e_select_sp, c_all_sp, c_select_sp, b_all_sp, b_select_sp, ncol = 2)

###############################################################################
###############################################################################
# plot Pollen on X Axis & GCM on Y
x <- data.frame(bvPollen$time, bvPollen$centroidVelocity, 
                bvBeyerMeans$centroidVelocity, bvECBiltMean$centroidVelocity,
                bvCCSMMean$centroidVelocity)

# Beyer vs Pollen
b <- ggplot(x, aes(bvPollen.centroidVelocity, bvBeyerMeans.centroidVelocity, 
              color = bvPollen.time, label = bvPollen.time)) + 
  geom_point() + 
  geom_text_repel(min.segment.length = 0) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(title = "Pollen vs Beyer Biotic Velocities", x = "Pollen Centroid Velocity (m/yr)", 
       y = "Beyer Centroid Velocity (m/yr)")

# ECBilt vs Pollen
e <- ggplot(x, aes(bvPollen.centroidVelocity, bvECBiltMean.centroidVelocity, 
              color = bvPollen.time, label = bvPollen.time)) + 
  geom_point() + 
  geom_text_repel(min.segment.length = 0) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(title = "Pollen vs ECBilt Biotic Velocities", x = "Pollen Centroid Velocity (m/yr)", 
       y = "ECBilt Centroid Velocity (m/yr)")

# CCSM vs Pollen
c <- ggplot(x, aes(bvPollen.centroidVelocity, bvCCSMMean.centroidVelocity, 
              color = bvPollen.time, label = bvPollen.time)) + 
  geom_point() + 
  geom_text_repel(min.segment.length = 0) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(title = "Pollen vs CCSM Biotic Velocities", x = "Pollen Centroid Velocity (m/yr)", 
       y = "CCSM Centroid Velocity (m/yr)")
plot(b)
plot(e)
plot(c)
plot_grid(b, e, c, ncol = 2)

dev.off()
