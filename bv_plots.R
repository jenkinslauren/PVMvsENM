rm(list=ls())
library(enmSdm)
library(gtools)
library(ggplot2)

gcm <- 'pollen'
fileName <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/0ByjNJEf91IW5SUlEOUJFVGN0Y28/NSF_ABI_2018_2021/data_and_analyses/pg_pollen/matern_overdispersed/predictions-FRAXINUS_meanpred.tif'
pollenRast <- brick(fileName)

pollenStack <- stack(fileName)
projection(pollenStack) <- getCRS('albersNA')
# note: there is no difference in output between onlyInSharedCells = T & F
bvPollen <- bioticVelocity(pollenStack, times = seq(-21000,0, by=1000), onlyInSharedCells = T)

bvPollen$time <- paste0(abs(bvPollen$timeFrom), '-', abs(bvPollen$timeTo), ' kybp')
bvPollen$time <- factor(bvPollen$time, levels = rev(mixedsort(bvPollen$time)))
bvPollen$time <- as.factor(bvPollen$time)

# biotic velocity
p <- ggplot(bvPollen, aes(time, centroidVelocity)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Pollen") + xlab("time period") + ylab("centroid velocity (m/yr)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

load(paste0('./workspaces/06 - Beyer Projections'))
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvMeans <- bioticVelocity(stackMeansList, times = seq(-21,0, by=1), onlyInSharedCells = T)

stackMaxList <- stack(maxList)
projection(stackMaxList) <- getCRS('albersNA')
bvMax <- bioticVelocity(stackMaxList, times = seq(-21,0, by=1), onlyInSharedCells = T)

minBV <- min(c(min(bvMeans$centroidVelocity), 
               min(bvMax$centroidVelocity))) - 1000
maxBV <- max(c(max(bvMeans$centroidVelocity), 
               max(bvMax$centroidVelocity))) + 1000

b <- ggplot(bvMeans, aes(time, centroidVelocity)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Means") + xlab("time period") + ylab("centroid velocity") +
  ylim(0, maxBV) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

load(paste0('./workspaces/06 - Lorenz_ccsm Projections'))
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvMeans <- bioticVelocity(stackMeansList, times = seq(-21,0, by=1), onlyInSharedCells = T)

stackMaxList <- stack(maxList)
projection(stackMaxList) <- getCRS('albersNA')
bvMax <- bioticVelocity(stackMaxList, times = seq(-21,0, by=1), onlyInSharedCells = T)

minBV <- min(c(min(bvMeans$centroidVelocity), 
               min(bvMax$centroidVelocity))) - 1000
maxBV <- max(c(max(bvMeans$centroidVelocity), 
               max(bvMax$centroidVelocity))) + 1000

c <- ggplot(bvMeans, aes(time, centroidVelocity)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Means") + xlab("time period") + ylab("centroid velocity") +
  ylim(0, maxBV) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

load(paste0('./workspaces/06 - ecbilt Projections'))
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvMeans <- bioticVelocity(stackMeansList, times = seq(-21,0, by=1), onlyInSharedCells = T)

stackMaxList <- stack(maxList)
projection(stackMaxList) <- getCRS('albersNA')
bvMax <- bioticVelocity(stackMaxList, times = seq(-21,0, by=1), onlyInSharedCells = T)

minBV <- min(c(min(bvMeans$centroidVelocity), 
               min(bvMax$centroidVelocity))) - 1000
maxBV <- max(c(max(bvMeans$centroidVelocity), 
               max(bvMax$centroidVelocity))) + 1000

e <- ggplot(bvMeans, aes(time, centroidVelocity)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Means") + xlab("time period") + ylab("centroid velocity") +
  ylim(0, maxBV) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))