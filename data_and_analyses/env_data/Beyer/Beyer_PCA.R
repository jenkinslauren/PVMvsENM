rm(list = ls())
library(raster)
library(enmSdm)
library(dplyr)
library(terra)
library(gtools)
setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

# set constants
pc <- 5

# load(paste0('./PCA_Beyer_PC', pc))

# retrieve list of files
fileList <- list.files(path = './data_and_analyses/env_data/Beyer/tifs', 
                   pattern='*.tif', all.files=TRUE, full.names=TRUE)
fileList <- mixedsort(fileList)
clim <- lapply(fileList, brick)
clim <- stack(clim)
climDf <- as.data.frame(clim)

climList <- list()
vars <- c('BIO1', paste0('BIO', 4:18), 'cloudiness', 'relative_humidity')
yrs <- paste0(1:22)
for(i in 1:22) {
  climTimes <- list()
  for(var in vars) {
    index <- paste0(var, ".", i)
    climTimes <- append(climTimes, clim[[index]])
  }
  climList <- append(climList,brick(climTimes))
}

clim <- lapply(climList, stack)

# build climate over all time periods, not just the latest
# using 10,000 random background points, extract climate at these points
# for each time period
buildClim <- function(brick) {
  # df <- data.frame()
  for (i in 1:nlayers(brick)) {
    # extract environment at each site
    randomBgEnv <- raster::extract(brick, randomBgSites)
    randomBgEnv <- as.data.frame(randomBgEnv)
    
    # remove any NAs for at least one variable
    # isNa <- is.na(rowSums(randomBgEnv))
    # if (any(isNa)) {
    #   randomBgSites <- randomBgSites[-which(isNa), ]
    #   randomBgEnv <- randomBgEnv[-which(isNa), ]
    # }
    randomBg <- cbind(randomBgSites, randomBgEnv)
    names(randomBg)[1:3] <- c('longitude', 'latitude', 
                              sub("\\..*", "", names(brick[[1]])))
    # head(randomBg)
    df <- rbind(df, randomBg)
  }
  return(df)
}

listofDf <- list()
bricks <- lapply(fileList, brick)
df <- data.frame()
# list of dataframes containing the data for each variable
# so, listofDf[1] = bioclim1, listofDf[2] = bioclim4
# recall the last 2 variables (18 & 19) = cloudiness & relative_humidity
randomBgSites <- randomPoints(clim[[1]], 20000)
df <- lapply(bricks, buildClim)

stackDfs <- list()
for (n in 1:length(df)) {
  dfn <- df[[n]]
  dfn <- data.frame(dfn[1:2], stack(dfn[3:ncol(dfn)]))
  names(dfn)[3] <- as.character(dfn$ind[1])
  dfn$ind <- NULL
  stackDfs[[n]] <- dfn
}

# dfn <- data.frame(dfn[1:2], stack(dfn[3:ncol(dfn)]))
# names(dfn)[3] <- as.character(dfn$ind[1])
# dfn$ind <- NULL
# stackDfs[[7]] <- dfn
# 
# for(df in stackDfs){
#   dfFinal <- cbind(dfFinal, )
# }

df1 <- stackDfs[[1]]
df2 <- stackDfs[[2]]
dfFinal <- cbind(df1, df2[,3])
names(dfFinal)[4] <- colnames(df2)[3]

for (n in 3:length(stackDfs)) {
  temp <- stackDfs[[n]]
  dfFinal <- cbind(dfFinal, temp[,3])
  names(dfFinal)[n + 2] <- colnames(temp)[3]
}


climxDf <- dfFinal
climRast <- rasterFromXYZ(climxDf)
climx <- stack(climRast)
climxDf <- as.data.frame(climx)

# apply labels to each variable
# recall that bioclim 2 & 3 aren't included
# names(climxDf) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
#                     "bio17", "bio18", "bio19", "bio4", "bio5", "bio6", "bio7",
#                     "bio8", "bio9", "cloudiness", "relative_humidity")
names(climDf) <- c(rep("bio1",22), rep("bio4",22), rep("bio5",22), rep("bio6",22), 
                   rep("bio7",22),rep("bio8",22), rep("bio9",22), rep("bio10",22), 
                   rep("bio11",22), rep("bio12",22), rep("bio13",22), rep("bio14",22), 
                   rep("bio15",22), rep("bio16",22), rep("bio17",22), rep("bio18",22), 
                  rep("cloudiness",22), rep("relative_humidity",22))

# crops to NA, don't need to
# 		clim <- crop(clim, nam0Sp)
# 		clim <- clim * mask
		
# do this after making the PCAs, otherwise it will take a while
# fill NA cells near coasts to account for fact that some records 
#   may not fall in a cell near a coast (don't want trees in the ocean)
# for (i in 1:nlayers(envPca)) {
#   envPca[[i]] <- focal(envPca[[i]], w=matrix(1, nrow=3, ncol=3), fun=mean, na.rm=TRUE, NAonly=TRUE)
# }
# names(envPca) <- names(envPca) <- c("bio1", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10",
#                    "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
#                    "bio18", "bio19", "cloudiness", "relative_humidity")

# PCA on climate (only)
nonNas <- which(complete.cases(climxDf))
climxDf <- climxDf[nonNas, ]
climxDf <- climxDf[1:10000,]

nonNas <- which(complete.cases(climDf))
climDf <- climDf[nonNas, ]
climDf <- climDf[1:10000,]

pca <- prcomp(climxDf, center = TRUE, scale = TRUE)
# prcomp calculates the eigenvectors on correlation matrix
# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab10_PCA.html#Principal_Component_Analysis 
# in pca: the coefficients reflect the contribution of each variable to each principal component
# "proportion of variance" tells us the percentage of variance explained 
#   by that respective principal component

# summary(pca)
# pca_scores <- data.frame(pca$x)
# head(pca_scores)
# df <- cbind(climDf, pca_scores)
# biplot(pca, xlabs=rep('.', nrow(pca$x)))

# want bricks by timestep, not by variable
# brick[["bio1.22"]]
# names(clim) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
#                  "bio17", "bio18", "bio19", "bio4", "bio5", "bio6", "bio7",
#                  "bio8", "bio9", "cloudiness", "relative_humidity")

pcPrediction <- list()
climYears <- seq(0, 21000, by = 1000)
for (i in 1:length(clim)) {
  names(clim[[i]]) <- c("BIO1", paste0('BIO', 4:18), "cloudiness", "relative_humidity")
  pcPrediction[i] <- raster::predict(clim[[i]], pca, index = 1:pc)
  names(pcPrediction[[i]]) <- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
}

for(i in 1:22) {
  writeRaster(pcPrediction[[i]], paste0('./stackPca/Beyer_yr', i, '.tif'), format = 'GTiff')
}
f <- list.files('./stackPca', pattern = '.tif', full.names = T)
f <- mixedsort(f)
redblue <- colorRampPalette(c("red","orange","blue"))
pdf(file = "pcaStackPlots.pdf", width = 11, height = 8.5)
for(i in 1:22) {
  n <- f[i]
  plot(raster(n), main = names(raster(n)), col = redblue(7))
}
dev.off()

stackPca <- stack(pcPrediction)
writeRaster(stackPca, paste0('./stackPca/Beyer_PC', pc, '.tif'), format = 'GTiff',
            bylayer = T, overwrite = T)
f <- list.files('./stackPca', pattern = '.tif', full.names = T)
f <- mixedsort(f)
redblue <- colorRampPalette(c("red","orange","blue"))
pdf(file = "pcaStackPlots.pdf", width = 11, height = 8.5)
for(i in 1:22) {
  n <- f[i*5]
  plot(raster(n), main = names(raster(n)), col = redblue(7))
}
dev.off()

plot(stackPca)

outfile <- paste0('./data_and_analyses/env_data/Beyer/tifs/pca_output_PC', pc, '.tif')
writeRaster(stackPca, outfile, format = 'GTiff', overwrite = T)
save.image(paste0('./PCA_Beyer_PC', pc))
