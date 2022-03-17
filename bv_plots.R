rm(list=ls())
library(enmSdm)
library(gtools)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tidyverse)

gcm <- 'pollen'
fileName <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/0ByjNJEf91IW5SUlEOUJFVGN0Y28/NSF_ABI_2018_2021/data_and_analyses/pg_pollen/matern_overdispersed/predictions-FRAXINUS_meanpred.tif'
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

plot_grid(p, bMean, bMax, cMean, cMax, eMean, eMax)

# bv <- data.frame(bvPollen$time, 
#                  bvPollen$centroidVelocity, 
#                  bvBeyerMax$centroidVelocity, 
#                  bvBeyerMeans$centroidVelocity,
#                  bvCCSMMax$centroidVelocity,
#                  bvCCSMMean$centroidVelocity,
#                  bvECBiltMax$centroidVelocity,
#                  bvECBiltMean$centroidVelocity)

bv <- data.frame(bvPollen$time, 
                 bvBeyerMeans$centroidVelocity,
                 bvCCSMMean$centroidVelocity,
                 bvECBiltMean$centroidVelocity)

colnames(bv)[1] <- 'time'
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bv2 <- tidyr::pivot_longer(bv, -time, names_to = "label", values_to = "centroidVelocity") # New way

ggplot(bv2, aes(x = time, y = centroidVelocity, color = label, group = 1)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bvPollen, aes(x = time, y = centroidVelocity, group = 1)) +
  geom_point(aes(color = cbbPalette[6])) + geom_line(aes(color = cbbPalette[6])) +
  geom_point(data = bvBeyerMeans, aes(y = centroidVelocity, color = cbbPalette[2])) +
  geom_line(data = bvBeyerMeans, aes(y = centroidVelocity, color = cbbPalette[2])) +
  geom_point(data = bvCCSMMean, aes(y = centroidVelocity, color = cbbPalette[3])) +
  geom_line(data = bvCCSMMean, aes(y = centroidVelocity, color = cbbPalette[3])) +
  geom_point(data = bvECBiltMean, aes(y = centroidVelocity, color = cbbPalette[4])) +
  geom_line(data = bvECBiltMean, aes(y = centroidVelocity, color = cbbPalette[4])) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Biotic Velocity for each GCM\n", x = "Time Period", y = "Centroid Velocity (m/yr)", 
       color = "GCM\n") +
  scale_color_manual(labels = c("Pollen", "Beyer", "CCSM", "ECBilt"))

ggplot(bvPollen, aes(x = time, y = centroidVelocity, group = 1)) +
  geom_point() + geom_line() +
  geom_point(data = bvBeyerMeans, aes(y = centroidVelocity)) +
  geom_line(data = bvBeyerMeans, aes(y = centroidVelocity)) +
  geom_point(data = bvCCSMMean, aes(y = centroidVelocity)) +
  geom_line(data = bvCCSMMean, aes(y = centroidVelocity)) +
  geom_point(data = bvECBiltMean, aes(y = centroidVelocity)) +
  geom_line(data = bvECBiltMean, aes(y = centroidVelocity)) +
  labs(title = "Biotic Velocity for each GCM\n", x = "Time Period", y = "Centroid Velocity (m/yr)", 
       color = "GCM\n") +
  scale_color_manual(labels = c("Pollen", "Beyer", "CCSM", "ECBilt"), 
                     values = c(cbbPalette[6], cbbPalette[2], cbbPalette[3], cbbPalette[4])) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = bv, aes(x = time)) +
  geom_point(aes(y = bvPollen.centroidVelocity)) +
  geom_line(aes(y = bvPollen.centroidVelocity)) +
  geom_point(aes(y = bvBeyerMax.centroidVelocity)) +
  geom_line(aes(y = bvBeyerMax.centroidVelocity)) +
  geom_point(aes(y = bvBeyerMeans.centroidVelocity)) +
  geom_line(aes(y = bvBeyerMeans.centroidVelocity)) +
  geom_point(aes(y = bvCCSMMax.centroidVelocity)) +
  geom_line(aes(y = bvCCSMMax.centroidVelocity)) +
  geom_point(aes(y = bvCCSMMean.centroidVelocity)) +
  geom_line(aes(y = bvCCSMMean.centroidVelocity)) +
  geom_point(aes(y = bvECBiltMax.centroidVelocity)) +
  geom_line(aes(y = bvECBiltMax.centroidVelocity)) +
  geom_point(aes(y = bvECBiltMean.centroidVelocity)) +
  geom_line(aes(y = bvECBiltMean.centroidVelocity))

x <- setNames(list(bvPollen, bvBeyerMeans, bvCCSMMean, bvECBiltMean), 
         c("Pollen","Beyer","CCSM", "ECBilt"))
xf <- data.frame(x[[1]]$time, x[[1]]$centroidVelocity, x[[2]]$centroidVelocity, 
                 x[[3]]$centroidVelocity, x[[4]]$centroidVelocity) %>% 
  setNames(c("time", "Pollen", "Beyer", "CCSM", "ECBilt"))

ggplot(aes(x = time, y = centroidVelocity, color = gcmKey, linetype = source)) +
  geom_line() +
  scale_colour_manual(values=c('red','blue','green','pink')) +
  theme_classic()
