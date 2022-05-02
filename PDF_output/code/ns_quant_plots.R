rm(list=ls())
library(enmSdm)
library(gtools)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tidyverse)
library(raster)
library(utils)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

pdf <- pdf(file = './PDF_output/ns_quant.pdf', width = 11, height = 8.5)
  
gcm <- 'pollen'
fileName <- '/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred.tif'
pollenRast <- brick(fileName)

pollenStack <- stack(fileName)
projection(pollenStack) <- getCRS('albersNA')
# note: there is no difference in output between onlyInSharedCells = T & F
bvPollen <- bioticVelocity(pollenStack, times = seq(-21000,0, by=1000), onlyInSharedCells = T)
bvPollenF <- bioticVelocity(pollenStack, times = seq(-21000,0, by=1000), onlyInSharedCells = F)

bvPollen$time <- paste0(abs(bvPollen$timeFrom)/1000, '-', abs(bvPollen$timeTo)/1000, ' kybp')
bvPollen$time <- factor(bvPollen$time, levels = rev(mixedsort(bvPollen$time)))

bvPollenF$time <- paste0(abs(bvPollenF$timeFrom)/1000, '-', abs(bvPollenF$timeTo)/1000, ' kybp')
bvPollenF$time <- factor(bvPollenF$time, levels = rev(mixedsort(bvPollenF$time)))

pt05 <- ggplot(bvPollen, aes(time, nsQuantVelocity_quant0p05, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

pf05 <- ggplot(bvPollenF, aes(time, nsQuantVelocity_quant0p05, group = 1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

pt95 <- ggplot(bvPollen, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

pf95 <- ggplot(bvPollenF, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +  
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

title05 <- ggdraw() + 
  draw_label(label = "Biotic Velocity of Southern Edge (5th quantile)\n", 
             fontface = 'bold', hjust = 0.5) +
  draw_label(label = "\nPollen", hjust = 0.5)

title95 <- ggdraw() + 
  draw_label("Biotic Velocity of Northern Edge (95th quantile)",
             fontface = 'bold', hjust = 0.5)

plot_grid(title05, plot_grid(pt05, pf05), title95,
          plot_grid(pt95, pf95), ncol = 1, rel_heights = c(0.2, 1))

load(paste0('./workspaces/06 - Beyer Projections'))
# Beyer Means
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvBeyerMeans <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = T)
bvBeyerMeansF <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = F)

bvBeyerMeans$time <- paste0(abs(bvBeyerMeans$timeFrom)/1000, '-', abs(bvBeyerMeans$timeTo)/1000, ' kybp')
bvBeyerMeans$time <- factor(bvBeyerMeans$time, levels = rev(mixedsort(bvBeyerMeans$time)))

bvBeyerMeansF$time <- paste0(abs(bvBeyerMeansF$timeFrom)/1000, '-', abs(bvBeyerMeansF$timeTo)/1000, ' kybp')
bvBeyerMeansF$time <- factor(bvBeyerMeansF$time, levels = rev(mixedsort(bvBeyerMeansF$time)))

bt05 <- ggplot(bvBeyerMeans, aes(time, nsQuantVelocity_quant0p05, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

bf05 <- ggplot(bvBeyerMeansF, aes(time, nsQuantVelocity_quant0p05, group = 1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

bt95 <- ggplot(bvBeyerMeans, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

bf95 <- ggplot(bvBeyerMeansF, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +  
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

title05 <- ggdraw() + 
  draw_label(label = "Biotic Velocity of Southern Edge (5th quantile)\n", 
             fontface = 'bold', hjust = 0.5) +
  draw_label(label = "\nGCM = HadAM3H", hjust = 0.5)

title95 <- ggdraw() + 
  draw_label("Biotic Velocity of Northern Edge (95th quantile)",
             fontface = 'bold', hjust = 0.5)

plot_grid(title05, plot_grid(bt05, bf05), title95,
          plot_grid(bt95, bf95), ncol = 1, rel_heights = c(0.2, 1))

###############################################################################
######### CCSM
###############################################################################

load(paste0('./workspaces/06 - Lorenz_ccsm Projections'))
# Lorenz_ccsm Means
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvCCSMMean <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = T)
bvCCSMMeanF <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = F)

bvCCSMMean$time <- paste0(abs(bvCCSMMean$timeFrom)/1000, '-', abs(bvCCSMMean$timeTo)/1000, ' kybp')
bvCCSMMean$time <- factor(bvCCSMMean$time, levels = rev(mixedsort(bvCCSMMean$time)))

bvCCSMMeanF$time <- paste0(abs(bvCCSMMeanF$timeFrom)/1000, '-', abs(bvCCSMMeanF$timeTo)/1000, ' kybp')
bvCCSMMeanF$time <- factor(bvCCSMMeanF$time, levels = rev(mixedsort(bvCCSMMeanF$time)))

ct05 <- ggplot(bvCCSMMean, aes(time, nsQuantVelocity_quant0p05, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

cf05 <- ggplot(bvCCSMMeanF, aes(time, nsQuantVelocity_quant0p05, group = 1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

ct95 <- ggplot(bvCCSMMean, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

cf95 <- ggplot(bvCCSMMeanF, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +  
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

title05 <- ggdraw() + 
  draw_label(label = "Biotic Velocity of Southern Edge (5th quantile)\n", 
             fontface = 'bold', hjust = 0.5) +
  draw_label(label = "\nGCM = CCSM", hjust = 0.5)

title95 <- ggdraw() + 
  draw_label("Biotic Velocity of Northern Edge (95th quantile)",
             fontface = 'bold', hjust = 0.5)

plot_grid(title05, plot_grid(ct05, cf05), title95,
          plot_grid(ct95, cf95), ncol = 1, rel_heights = c(0.2, 1))

###############################################################################
######### ECBILT
###############################################################################

load(paste0('./workspaces/06 - ecbilt Projections'))
# ECBilt Means
stackMeansList <- stack(meansList)
projection(stackMeansList) <- getCRS('albersNA')
bvECBiltMean <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = T)
bvECBiltMeanF <- bioticVelocity(stackMeansList, times = seq(-21000, 0, by = 1000), onlyInSharedCells = F)

bvECBiltMean$time <- paste0(abs(bvECBiltMean$timeFrom)/1000, '-', abs(bvECBiltMean$timeTo)/1000, ' kybp')
bvECBiltMean$time <- factor(bvECBiltMean$time, levels = rev(mixedsort(bvECBiltMean$time)))

bvECBiltMeanF$time <- paste0(abs(bvECBiltMeanF$timeFrom)/1000, '-', abs(bvECBiltMeanF$timeTo)/1000, ' kybp')
bvECBiltMeanF$time <- factor(bvECBiltMeanF$time, levels = rev(mixedsort(bvECBiltMeanF$time)))

et05 <- ggplot(bvECBiltMean, aes(time, nsQuantVelocity_quant0p05, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

ef05 <- ggplot(bvECBiltMeanF, aes(time, nsQuantVelocity_quant0p05, group = 1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

et95 <- ggplot(bvECBiltMean, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
  geom_point() + geom_line() +
  ggtitle("Shared Cells Only") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

ef95 <- ggplot(bvECBiltMeanF, aes(time, nsQuantVelocity_quant0p95, group = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +  
  geom_point() + geom_line() +
  ggtitle("All Cells") +
  xlab("time period") + ylab("velocity (m/yr)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5))

title05 <- ggdraw() + 
  draw_label(label = "Biotic Velocity of Southern Edge (5th quantile)\n", 
             fontface = 'bold', hjust = 0.5) +
  draw_label(label = "\nGCM = ECBilt", hjust = 0.5)

title95 <- ggdraw() + 
  draw_label("Biotic Velocity of Northern Edge (95th quantile)",
             fontface = 'bold', hjust = 0.5)

plot_grid(title05, plot_grid(et05, ef05), title95,
          plot_grid(et95, ef95), ncol = 1, rel_heights = c(0.2, 1))

dev.off()

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
colnames(bv) <- c("Time", "Beyer", "CCSM", "ECBilt")
bv2 <- cbind(bv[1], utils::stack(bv[2:4]))
colnames(bv2) <- c("Time", "centroidVelocity", "climateSource")

pdf(file = './bv_plots.pdf', height = 8.5, width = 11)

plot_grid(p, bMean, bMax, cMean, cMax, eMean, eMax)

ggplot(bv2, aes(x = Time, y = centroidVelocity, group = climateSource, color = climateSource)) +
  geom_point() + geom_line() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Biotic Velocity\n", x = "Time Period", 
       y = "Centroid Velocity (m/yr)", color = "GCM")

dev.off()
