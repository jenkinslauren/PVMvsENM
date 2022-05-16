rm(list = ls())
library(raster)
library(enmSdm)
library(dplyr)
library(terra)
library(gtools)
setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

# set constants, define gcm & number of pc's to use in predictions
gcm <- 'ecbilt'
pc <- 5
workingFolder <- paste0('./data_and_analyses/env_data/Lorenz/V2/',
                        gcm, '_21-0k_all_tifs_LJ')
# load(paste0('./PCA_Lorenz_', gcm, '_PC', pc))

# retrieve list of files
dirList <- list.dirs(path = workingFolder,
                     recursive = FALSE)
dirList <- mixedsort(dirList)

if(exists('clim')) rm(clim)

for (dir in dirList) {
  fileList <- list.files(path = dir, pattern = '*.tif', 
                         all.files = TRUE, full.names = TRUE)
  
  # remove wdi & etr variables, leaving only 38 variables 
  fileList <- Filter(function(x) !any(grepl("ETR", x)), fileList)
  fileList <- Filter(function(x) !any(grepl("WDI", x)), fileList)
  # print(paste0(sub('.+LJ/(.+)', '\\1', dir), 'YBP ', length(fileList)))
  thisClim <- lapply(fileList, raster)
  thisClim <- brick(thisClim)
  names(thisClim) <- paste0(names(thisClim), '_', 
                            sub('.+LJ/(.+)', '\\1', dir))
  clim <- if(exists('clim')) {
    stack(clim, thisClim)
  } else {
    thisClim
  }
}

bricks <- list()
for (j in 1:38) {
  stack <- list()
  x <- j
  while(x <= nlayers(clim)) {
    stack <- append(stack, clim[[x]])
    x <- x + 38
  }
  bricks <- append(bricks, brick(stack))
}

i <- 1
vars <- sub('\\_[0-9].*', '', names(thisClim))
for (var in vars) {
  outfile <- paste0(workingFolder, '/', var, '.tif')
  writeRaster(bricks[[i]], outfile, format = 'GTiff', overwrite = T)
  i <- i + 1
}

fileList <- list.files(path = workingFolder, 
                       pattern='*.tif', all.files = TRUE, full.names = TRUE)
lorenz <- lapply(fileList, brick)
lorenz <- stack(lorenz)
climDf <- as.data.frame(lorenz)

vars <- sub('\\.tif.*', '', list.files(path = workingFolder, 
                                       pattern='*.tif', all.files = TRUE, full.names = FALSE))

climList <- list()
for(i in 1:22) {
  climTimes <- list()
  for(var in vars) {
    index <- paste0(var, ".", i)
    climTimes <- append(climTimes, lorenz[[i]])
  }
  climList <- append(climList,brick(climTimes))
}

lorenz <- lapply(climList, stack)

for (i in 1:length(lorenz)) {
  names(lorenz[[i]]) <- vars
}

# build climate over all time periods, not just the latest
# using 10,000 random background points, extract climate at these points
# for each time period
randomBgSites <- randomPoints(lorenz[[1]], 10000)
buildClim <- function(brick) {
  # df <- data.frame()
  for (i in 1:nlayers(brick)) {
    # extract environment at each site
    randomBgEnv <- extract(brick, randomBgSites)
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

# list of dataframes containing the data for each variable
# so, listofDf[1] = first variable, listofDf[2] = second variable
df <- data.frame()
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

nonNas <- which(complete.cases(climDf))
climDf <- climDf[nonNas, ]

pca <- prcomp(climxDf, center = TRUE, scale = TRUE)
fileName <- if(gcm == 'Lorenz_ccsm') { 
  paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/pca_pc', pc, '.Rdata') 
  } else {
  paste0('./data_and_analyses/env_data/Lorenz/V2/ecbilt_21-0k_all_tifs_LJ/pca_pc', pc, '.Rdata') 
    }

save(pca, file = fileName)

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
for (i in 1:length(lorenz)) {
  names(lorenz[[i]]) <- vars
  pcPrediction[i] <- raster::predict(lorenz[[i]], pca, index = 1:pc)
  names(pcPrediction[[i]]) <- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
}

stackPca <- stack(pcPrediction)
plot(stackPca)

outfile <- paste0(workingFolder, '/pca_output_pc', pc,'.tif')
writeRaster(stackPca, outfile, format = 'GTiff', overwrite = T)
save.image(paste0('./workspaces/PCA_', gcm, '_PC', pc))
