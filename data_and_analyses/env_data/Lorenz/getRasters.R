rm(list = ls())
library(raster)
library(enmSdm)
library(dplyr)
setwd('/Volumes/LJ MacBook Backup/MOBOT/PVMvsENM')
lorenzPath <- './data_and_analyses/env_data/Lorenz/V2'

rasts <- list()
getClimRasts <- function(gcm, years, variables) {
  
  # gcm		'ccsm' or 'ecbilt'
  # year		year BP (from 0 to 21000 for ccsm or 22000 for ecbilt)
  # variables names of variables
  # rescale	TRUE ==> rescale rasters to [0, 1] using present-day values for min/max
  # fillCoasts FALSE ==> use rasters as-is; TRUE ==> extrapolate to NA cells immediately adjacent to non-NA cells (typically coastal cells)
  
  gcmFolder <- if (gcm == 'ccsm') {
    'ccsm_21-0k_all_tifs_LJ'
  } else if (gcm == 'ecbilt') {
    'ecbilt_21-0k_all_tifs'
  }
  # get current version of each variable
  
  for (variable in variables) {
    for (year in years) {
      rast <- raster(paste0(lorenzPath, '/', gcmFolder, '/', year, 'BP/used/', variable))
      names(rast) <- paste0(variable, '_', year, 'BP')
      # summary(rast)
      rasts <- append(rasts, rast)
    }
  }
  return(rasts)
}

variables <- list.files(path = './data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/0BP/used/', 
                        pattern='*.tif', all.files=TRUE)
years <- c('0', '1000', '2000', '3000', '4000', '5000', '6000', '7000', '8000',
           '9000', '10000', '11000', '12000', '13000', '14000', '15000', '16000',
           '17000', '18000', '19000', '20000', '21000')

rasts <- getClimRasts('ccsm', years, variables)
rastBricks <- list(brick(rasts[1:22]), brick(rasts[23:44]), brick(rasts[45:66]), 
                   brick(rasts[67:88]), brick(rasts[89:110]), brick(rasts[111:132]), 
                   brick(rasts[133:154]), brick(rasts[155:176]), brick(rasts[177:198]))
clim <- stack(rastBricks)
climDf <- as.data.frame(clim)

# build climate over all time periods, not just the latest
# using 10,000 random backgound points, extracted climate at these points
# for each time period
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
df <- data.frame()
# list of dataframes containing the data for each variable
# so, listofDf[1] = bioclim1, listofDf[2] = bioclim4
# recall the last 2 variables (18 & 19) = cloudiness & relative_humidity
randomBgSites <- randomPoints(clim, 10000)
df <- lapply(rastBricks, buildClim)

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
names(climxDf) <- c("etr", "tmax", "tmin", "aet", "gdd0", "gdd5", "pet", "prcp",
                    "wdi")
names(climDf) <- c(rep("etr",22), rep("tmax",22), rep("tmin",22), rep("aet",22), 
                   rep("gdd0",22), rep("gdd5",22), rep("pet",22), rep("prcp",22),
                   rep("wdi",22))
names(clim) <- c(rep("etr",22), rep("tmax",22), rep("tmin",22), rep("aet",22), 
                 rep("gdd0",22), rep("gdd5",22), rep("pet",22), rep("prcp",22),
                 rep("wdi",22))

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

pcPredictionNoNas <- raster::predict(pca, climDf)
colnames(pcPredictionNoNas) <- paste0('pc', 1:9)

# predict PCA back to rasters
pcPrediction <- as.data.frame(clim)
pcPrediction[nonNas, ] <- pcPredictionNoNas

pcaRasts <- clim * NA
for (pc in 1:9) pcaRasts <- setValues(pcaRasts, values=pcPrediction[ , pc], layer=pc)
# names(pcaRasts) <- paste0('pc', 1:9)

envPca <- stack(clim, pcaRasts)
plot(envPca)
