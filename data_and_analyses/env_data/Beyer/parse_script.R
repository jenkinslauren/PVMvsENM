# Beyer et al. 2020 Environmental Data
# Lauren Jenkins, Fall 2021
# # # # # # # # # # # # # # # #

rm(list=ls())
setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')
library(ncdf4)
library(raster)

# File path on Lauren's computer:
file <- "./data_and_analyses/env_data/Beyer/LateQuaternary_Environment.nc"

createRasters <- function(varName) {
  if (grepl('bio', varName)) {
    upperVar <- toupper(varName)
  } else {
    upperVar <- varName
  }
  
  outfile <- paste0('./data_and_analyses/env_data/Beyer/tifs/', upperVar, '/', upperVar, '.tif')
  
  var <- brick(file, varname = upperVar, level = 1)
  
  # can check this by overlaying the US map and checking edges
  projection(var) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(stack(var), outfile, bylayer=TRUE, format='GTiff', overwrite = T)
  
  rasterList <- list.files(path = paste0('./data_and_analyses/env_data/Beyer/tifs/', upperVar), 
                           pattern='*.tif', all.files=TRUE, full.names=TRUE)
  allrasters <- lapply(rasterList, raster)
  allrasters <- stack(allrasters)
}

rasters <- c('bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 
             'bio18', 'bio19')
lapply(rasters, createRasters)

stack <- list()
means <- list()


for(i in 1:12){
  
  f <- list(brick(file, varname = 'relative_humidity', level = i)) 
  stack <- append(stack, f)
}

meansStack <- list()

s <- stack(stack[[1]]$X.21000, stack[[2]]$X.21000, stack[[3]]$X.21000, stack[[4]]$X.21000,
           stack[[5]]$X.21000, stack[[6]]$X.21000, stack[[7]]$X.21000, stack[[8]]$X.21000,
           stack[[9]]$X.21000, stack[[10]]$X.21000, stack[[11]]$X.21000, stack[[12]]$X.21000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.20000, stack[[2]]$X.20000, stack[[3]]$X.20000, stack[[4]]$X.20000,
           stack[[5]]$X.20000, stack[[6]]$X.20000, stack[[7]]$X.20000, stack[[8]]$X.20000,
           stack[[9]]$X.20000, stack[[10]]$X.20000, stack[[11]]$X.20000, stack[[12]]$X.20000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.19000, stack[[2]]$X.19000, stack[[3]]$X.19000, stack[[4]]$X.19000,
           stack[[5]]$X.19000, stack[[6]]$X.19000, stack[[7]]$X.19000, stack[[8]]$X.19000,
           stack[[9]]$X.19000, stack[[10]]$X.19000, stack[[11]]$X.19000, stack[[12]]$X.19000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.18000, stack[[2]]$X.18000, stack[[3]]$X.18000, stack[[4]]$X.18000,
           stack[[5]]$X.18000, stack[[6]]$X.18000, stack[[7]]$X.18000, stack[[8]]$X.18000,
           stack[[9]]$X.18000, stack[[10]]$X.18000, stack[[11]]$X.18000, stack[[12]]$X.18000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.17000, stack[[2]]$X.17000, stack[[3]]$X.17000, stack[[4]]$X.17000,
           stack[[5]]$X.17000, stack[[6]]$X.17000, stack[[7]]$X.17000, stack[[8]]$X.17000,
           stack[[9]]$X.17000, stack[[10]]$X.17000, stack[[11]]$X.17000, stack[[12]]$X.17000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.16000, stack[[2]]$X.16000, stack[[3]]$X.16000, stack[[4]]$X.16000,
           stack[[5]]$X.16000, stack[[6]]$X.16000, stack[[7]]$X.16000, stack[[8]]$X.16000,
           stack[[9]]$X.16000, stack[[10]]$X.16000, stack[[11]]$X.16000, stack[[12]]$X.16000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.15000, stack[[2]]$X.15000, stack[[3]]$X.15000, stack[[4]]$X.15000,
           stack[[5]]$X.15000, stack[[6]]$X.15000, stack[[7]]$X.15000, stack[[8]]$X.15000,
           stack[[9]]$X.15000, stack[[10]]$X.15000, stack[[11]]$X.15000, stack[[12]]$X.15000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.14000, stack[[2]]$X.14000, stack[[3]]$X.14000, stack[[4]]$X.14000,
           stack[[5]]$X.14000, stack[[6]]$X.14000, stack[[7]]$X.14000, stack[[8]]$X.14000,
           stack[[9]]$X.14000, stack[[10]]$X.14000, stack[[11]]$X.14000, stack[[12]]$X.14000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.13000, stack[[2]]$X.13000, stack[[3]]$X.13000, stack[[4]]$X.13000,
           stack[[5]]$X.13000, stack[[6]]$X.13000, stack[[7]]$X.13000, stack[[8]]$X.13000,
           stack[[9]]$X.13000, stack[[10]]$X.13000, stack[[11]]$X.13000, stack[[12]]$X.13000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.12000, stack[[2]]$X.12000, stack[[3]]$X.12000, stack[[4]]$X.12000,
           stack[[5]]$X.12000, stack[[6]]$X.12000, stack[[7]]$X.12000, stack[[8]]$X.12000,
           stack[[9]]$X.12000, stack[[10]]$X.12000, stack[[11]]$X.12000, stack[[12]]$X.12000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.11000, stack[[2]]$X.11000, stack[[3]]$X.11000, stack[[4]]$X.11000,
           stack[[5]]$X.11000, stack[[6]]$X.11000, stack[[7]]$X.11000, stack[[8]]$X.11000,
           stack[[9]]$X.11000, stack[[10]]$X.11000, stack[[11]]$X.11000, stack[[12]]$X.11000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.10000, stack[[2]]$X.10000, stack[[3]]$X.10000, stack[[4]]$X.10000,
           stack[[5]]$X.10000, stack[[6]]$X.10000, stack[[7]]$X.10000, stack[[8]]$X.10000,
           stack[[9]]$X.10000, stack[[10]]$X.10000, stack[[11]]$X.10000, stack[[12]]$X.10000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.9000, stack[[2]]$X.9000, stack[[3]]$X.9000, stack[[4]]$X.9000,
           stack[[5]]$X.9000, stack[[6]]$X.9000, stack[[7]]$X.9000, stack[[8]]$X.9000,
           stack[[9]]$X.9000, stack[[10]]$X.9000, stack[[11]]$X.9000, stack[[12]]$X.9000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.8000, stack[[2]]$X.8000, stack[[3]]$X.8000, stack[[4]]$X.8000,
           stack[[5]]$X.8000, stack[[6]]$X.8000, stack[[7]]$X.8000, stack[[8]]$X.8000,
           stack[[9]]$X.8000, stack[[10]]$X.8000, stack[[11]]$X.8000, stack[[12]]$X.8000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.7000, stack[[2]]$X.7000, stack[[3]]$X.7000, stack[[4]]$X.7000,
           stack[[5]]$X.7000, stack[[6]]$X.7000, stack[[7]]$X.7000, stack[[8]]$X.7000,
           stack[[9]]$X.7000, stack[[10]]$X.7000, stack[[11]]$X.7000, stack[[12]]$X.7000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.6000, stack[[2]]$X.6000, stack[[3]]$X.6000, stack[[4]]$X.6000,
           stack[[5]]$X.6000, stack[[6]]$X.6000, stack[[7]]$X.6000, stack[[8]]$X.6000,
           stack[[9]]$X.6000, stack[[10]]$X.6000, stack[[11]]$X.6000, stack[[12]]$X.6000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.5000, stack[[2]]$X.5000, stack[[3]]$X.5000, stack[[4]]$X.5000,
           stack[[5]]$X.5000, stack[[6]]$X.5000, stack[[7]]$X.5000, stack[[8]]$X.5000,
           stack[[9]]$X.5000, stack[[10]]$X.5000, stack[[11]]$X.5000, stack[[12]]$X.5000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.4000, stack[[2]]$X.4000, stack[[3]]$X.4000, stack[[4]]$X.4000,
           stack[[5]]$X.4000, stack[[6]]$X.4000, stack[[7]]$X.4000, stack[[8]]$X.4000,
           stack[[9]]$X.4000, stack[[10]]$X.4000, stack[[11]]$X.4000, stack[[12]]$X.4000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.3000, stack[[2]]$X.3000, stack[[3]]$X.3000, stack[[4]]$X.3000,
           stack[[5]]$X.3000, stack[[6]]$X.3000, stack[[7]]$X.3000, stack[[8]]$X.3000,
           stack[[9]]$X.3000, stack[[10]]$X.3000, stack[[11]]$X.3000, stack[[12]]$X.3000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.2000, stack[[2]]$X.2000, stack[[3]]$X.2000, stack[[4]]$X.2000,
           stack[[5]]$X.2000, stack[[6]]$X.2000, stack[[7]]$X.2000, stack[[8]]$X.2000,
           stack[[9]]$X.2000, stack[[10]]$X.2000, stack[[11]]$X.2000, stack[[12]]$X.2000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X.1000, stack[[2]]$X.1000, stack[[3]]$X.1000, stack[[4]]$X.1000,
           stack[[5]]$X.1000, stack[[6]]$X.1000, stack[[7]]$X.1000, stack[[8]]$X.1000,
           stack[[9]]$X.1000, stack[[10]]$X.1000, stack[[11]]$X.1000, stack[[12]]$X.1000)
meansStack <- append(meansStack, s)

s <- stack(stack[[1]]$X0, stack[[2]]$X0, stack[[3]]$X0, stack[[4]]$X0,
           stack[[5]]$X0, stack[[6]]$X0, stack[[7]]$X0, stack[[8]]$X0,
           stack[[9]]$X0, stack[[10]]$X0, stack[[11]]$X0, stack[[12]]$X0)
meansStack <- append(meansStack, s)

# this will fill the means list from 21,000 ybp to 0 ybp, so
# means[[1]] is 21,000 ybp, but means [[72]] is 0 ybp
for (i in 1:length(meansStack)) {
  means[[i]] <- calc(meansStack[[i]], fun = mean, na.rm = TRUE)
}

outfile <- './data_and_analyses/env_data/Beyer/tifs/relative_humidity/relative_humidity.tif'
var <- brick(means)

# can check this by overlaying the US map and checking edges
projection(var) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

writeRaster(stack(var), outfile, bylayer = TRUE, format='GTiff', overwrite = T)
run('relative_humidity', 21000)

changeFileName <- function(i, year, var) {
  oldFile <- paste0('./data_and_analyses/env_data/Beyer/tifs/', var, '/', var, '_' , i, '.tif')
  # print(paste0("Old file name: ", oldFile))
  newFile <- paste0('./data_and_analyses/env_data/Beyer/tifs/', var, '/', var, '_' , year, 'ybp.tif')
  # print(paste0("New file name: ", newFile))
  file.rename(oldFile, newFile)
}

run <- function(var, yr) {
  i <- 1
  v <- var
  while(i <= 72) {
    changeFileName(i, yr, v)
    if (yr > 22000) {
      # decrease year by 2000 increments
      yr <- yr - 2000
    } else {
      # decrease year by 1000 increments
      yr <- yr - 1000
    }
    i <- i + 1
  }
}


variables <- c('BIO1', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9',
               'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 
               'BIO17', 'BIO18', 'BIO19')
lapply(variables, run, yr = 120000)

# check multiple timesteps
rasterList <- list.files(path = './data_and_analyses/env_data/Beyer/tifs/BIO6', 
                         pattern='*.tif', all.files=TRUE, full.names=TRUE)
allrasters <- lapply(rasterList, raster)
allrasters <- stack(allrasters)
plot(allrasters)

# make sure they saved correctly
test <- raster('./data_and_analyses/env_data/Beyer/tifs/BIO1/BIO1_1.tif')
world <- ne_countries(scale = "medium", returnclass = "sp")
plot(test)
plot(world, add = TRUE)

# write as stacks for PCA 
writeStacks <- function(var) {
  outfile <- paste0('./data_and_analyses/env_data/Beyer/tifs/', var, '.tif')
  
  rasterList <- list.files(path = paste0('./data_and_analyses/env_data/Beyer/tifs/', var), 
                           pattern='*.tif', all.files=TRUE, full.names=TRUE)
  allrasters <- lapply(rasterList, raster)
  allrasters <- stack(allrasters)
  
  writeRaster(allrasters, outfile, format = 'GTiff', overwrite = T)
}

variables <- c('BIO1', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9',
               'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 
               'BIO17', 'BIO18', 'BIO19', 'cloudiness', 'relative_humidity')
lapply(variables, writeStacks)

