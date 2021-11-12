# Carolina Ash
# Lauren Jenkins, Fall 2021

### CONTENTS ###
### setup ###

### setup ###
rm(list = ls())
setwd('/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/')

# load packages
library(BIEN)
library(dismo)
library(brglm2)
library(cluster)
library(maxnet)
library(raster)
library(rgbif)
library(rgeos)
library(geosphere)
library(rgdal)
library(scales)
library(sp)
library(omnibus)
library(enmSdm)
library(holoSimCell)
library(RColorBrewer)
library(rJava)
library(maps)

library(statisfactory)
library(legendary)

dirCreate('./species_records')
dirCreate('./figures_and_tables')
dirCreate('./regions/bien_range_map')

### obtain & clean data ### 
# use BIEN and BONAP
# visually check BIEN vs Little
# manually inspect/clean
## remove non-natives
## remove naturalized (may be more complex)
## look for outliers

ll <- c('longitude', 'latitude')

species <- 'Fraxinus caroliniana'

BIEN_ranges_species(species, directory = './regions/bien_range_map')
bienRange <- BIEN_ranges_load_species(species)
bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA'), TRUE)

# Little range maps
littleFileName <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_caroliniana/fraxcaro.shp'
littleRange <- shapefile(littleFileName)
projection(littleRange) <- getCRS('wgs84')
# littleRange <- littleRange[littleRange$CODE == 1,]
dirCreate('./regions/little_range_map')

# name file by when it was obtained
speciesFileName <- paste0('./species_records/00_', 
                          gsub(tolower(species), pattern = ' ', 
                               replacement = '_'), '_bien_all_occurrences.rda')

if (file.exists(speciesFileName)) {
  load(speciesFileName)
} else {
  occsRaw <- BIEN_occurrence_species(species = species, cultivated = FALSE,
                                     only.new.world = TRUE, all.taxonomy = FALSE,
                                     native.status = FALSE, natives.only = TRUE,
                                     observation.type = TRUE, political.boundaries = TRUE,
                                     collection.info = TRUE)
  sink('./species_records/bien_download.txt')
  say('Data representing species records were downloaded from BIEN on ', date(), '.')
  say('Bien Version 4.1')
  sink()
  save(occsRaw, file = speciesFileName)
}

# remove records with missing coordinates and dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]

# visualize
# convert to spatial object
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('wgs84', TRUE))
# occsSp <- sp::spTransform(occsSp, getCRS('albersNA', TRUE))

par(mfrow=c(1,2))
# plot BIEN
# plot(occsSp, main = 'BIEN range map')
plot(occsSp, pch = 16, cex = 0.5, col = 'black', main = 'BIEN range map')
plot(bienRange, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# plot Little
# plot(occsSp, main = 'Little range map')
plot(occsSp, pch = 16, cex = 0.5, col = 'black', main = 'Little range map')
plot(littleRange, col = alpha('green', 0.4), border = 'green', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# find which entries are outside the Little range for data cleaning/inspection
extractNA <- extract(littleRange, occsSp)
# determine which indices are outside range
indexNA <- which(!complete.cases(x))
outsideLittleRange <- occsSp$locality[indexNA]

title(sub=date(), outer=TRUE, cex.sub=0.7, line=-2)


## visualize environmental data
climate <- stack(listFiles('./data_and_analyses/env_data/Lorenz et al 2016 North America 21Kybp to 2100 CE/V2/ccsm3_22-0k_all_tifs/0BP'))
names(climate) <- paste0(prefix(1:50, 2))
climate
plot(climate)

