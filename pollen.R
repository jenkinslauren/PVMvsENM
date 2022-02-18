# pollen rasters
library(enmSdm)
library(gtools)

gcm <- 'pollen'
fileName <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/0ByjNJEf91IW5SUlEOUJFVGN0Y28/NSF_ABI_2018_2021/data_and_analyses/pg_pollen/matern_overdispersed/predictions-FRAXINUS_meanpred.tif'
pollenRast <- brick(fileName)

pollenStack <- stack(fileName)
projection(pollenStack) <- getCRS('albersNA')
# note: there is no difference in output between onlyInSharedCells = T & F
bvPollen <- bioticVelocity(pollenStack, times = seq(-21,0, by=1), onlyInSharedCells = T)

bvPollen$time <- paste0(abs(bvPollen$timeFrom), '-', abs(bvPollen$timeTo), ' kybp')
bvPollen$time <- factor(bvPollen$time, levels = rev(mixedsort(bvPollen$time)))

# biotic velocity
pdf(file = paste0('./predictions/pdf/', gcm, '_bioticVelocity.pdf'), width = 11, height = 8.5)
ggplot(bvPollen, aes(time, centroidVelocity)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Pollen") + xlab("time period") + ylab("centroid velocity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# NS Quant Velocity
pdf(file = paste0('./predictions/pdf/', gcm, '_nsQuantVelocity.pdf'), width = 11, height = 8.5)
ggplot(bvPollen, aes(time, nsQuantVelocity_quant0p05)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Pollen") + xlab("time period") + ylab("NS_quant_05") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bvPollen, aes(time, nsQuantVelocity_quant0p95)) + 
  geom_bar(stat = 'identity') + 
  ggtitle("Pollen") + xlab("time period") + ylab("NS_quant_95") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
