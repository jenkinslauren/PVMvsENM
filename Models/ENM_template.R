# SDM Template
# Author: Lauren Jenkins
# 25 January 2022
# Last updated: 18 February 2022

rm(list = ls())
# Load required packages.
# Load packages
library(enmSdm)

# for handling rasters
library(raster)
library(rgdal)
library(dismo)
# library(rgeos) # not used, but may be helpful when debugging

# library(terra)
# library(maxnet)

# visualization tools
library(sp)
library(sf)
library(maps)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

# additional tools
library(tools)
library(units)
library(dplyr)
library(tidyr)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

# nad27
default_crs = sf::st_crs(4267)

# setwd('/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/')
# load('./workspaces/01 - Modeling Workspace - Fraxinus Range Maps (BIEN + Little)')

# set constants
ll <- c('longitude', 'latitude')
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

pc <- 5
gcmList <- c('Beyer','Lorenz_ccsm', 'ecbilt')
# gcmList <- c('Beyer')
climYear <- 0

lorenzRast <- raster::raster('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/0BP/an_avg_ETR.tif')

for(sp in speciesList) {
  sp <- sp
  print(paste0("SPECIES = ", sp))
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  rangeName <- paste0('littleRange_', speciesAb)
  load(paste0('./data_and_analyses/study_region/regions/little_range_map/', 
              rangeName, '.Rdata'))
  if (file.exists(paste0("./workspaces/03 - Modeling Workspace - ", speciesAb, " Cleaning"))) {
    print('Already cleaned!')
  } else {
    range <- rangeMap
    
    # Standard cleaning procedures: 
    # * remove records without coordinates or dates
    # * remove records before 1900
    # * date_collected column is recorded twice, remove one of them
    
    speciesFileName <- paste0('./species_records/01_', 
                              gsub(tolower(species), pattern = ' ', 
                                   replacement = '_'), '_retained_records.rds')
    
    if (file.exists(speciesFileName)) {
      occs <- readRDS(speciesFileName)
    } else {
      load(paste0('./species_records/00_', 
                  gsub(tolower(species), pattern = ' ', 
                       replacement = '_'), '_bien_all_occurrences.rda'))
      
      # X occurrences
      
      # remove records without coordinates or dates to ensure locational reliability
      occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
      
      # remove records before 1900
      occs$Year <- substr(occs$date_collected, 1, 4)
      occs <- occs[occs$Year >= 1900, ]
      # X occurrences
      
      # date_collected column is there twice, remove one of them
      occs[18] <- NULL
      
      saveRDS(occs, speciesFileName)
    }
    
    
    # Dealing with imprecise coordinates DMS = FALSE.
    
    llOcc <- cbind(occs$longitude, occs$latitude)
    coordPrecis <- coordPrecision(llOcc)
    precision <- data.frame(llOcc, coordPrecis)
    
    hist(precision$coordPrecis,
         breaks=100,
         xlab='Coordinate uncertainty (m)',
         ylab='Number of records',
         main='',
         col='red'
    )
    
    
    # Dealing with imprecise coordinates DMS = TRUE.
    llOcc <- cbind(occs$longitude, occs$latitude)
    coordPrecisDMS <- coordPrecision(llOcc, dms = TRUE)
    precision <- data.frame(llOcc, coordPrecisDMS)
    
    hist(precision$coordPrecisDMS,
         breaks=100,
         xlab='Coord Precision (m)',
         ylab='Number of records',
         main='',
         col='red'
    )
    
    
    # Find the max between two precision measurements and remove occurrences where precision > 1,000m.
    pmaxPrecision <- pmax(coordPrecis, coordPrecisDMS)
    occs <- data.frame(occs, pmaxPrecision)
    print(paste0("Number of occurrences with precision > 1,000 m = ", 
                 length(which(occs$pmaxPrecision > 1000))))
    occs <- occs[which(occs$pmaxPrecision <= 1000), ]
    
    
    # Create spatial object for plotting.
    # Map of all observations:
    # Convert to spatial object for plotting
    par(mfrow=c(1,1))
    occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
    plot(occsSp, pch = 16, cex = 0.3, col = 'red', 
         main = "Little range map, all observations not cleaned")
    plot(range, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
    map("state", add = TRUE)
    map("world", add = TRUE)
    
    
    # Thin data.
    thinnedFileName <- paste0('./species_records/02_', 
                              gsub(tolower(species), pattern = ' ', 
                                   replacement = '_'), '_thinned_records.rData')
    
    if (file.exists(thinnedFileName)) {
      load(thinnedFileName)
    } else {
      # memory.limit() - will give you something with memory (how much is remaining or it's capacity)
      # lsos might also help with identifying how much memory
      # k <- nrow(occs) / 1000
      # # cut_number divides the dataframe into k groups by longitude value
      # occs$folds <- as.numeric(cut_number(occs$longitude, k))
      # folds <- split(occs, occs$folds)
      # 
      # # sanity check
      # print(paste0("Sanity check: does (# of occurrences) ", nrow(occs), 
      #              " = (# of occurrences in folds) ", sum(sapply(folds, nrow)), "?"))
      # 
      # i <- 1 
      # thinnedFolds <- list()
      # 
      # thin <- function(df) {
      #   p <- paste0("*X", i, ".")
      #   names(df) <- gsub(pattern = p, replacement = "", x = names(df))
      #   return(geoThinApprox(df, 1000, c('longitude', 'latitude')))
      # }
      # 
      # while(TRUE) {
      #   currentFold <- thin(data.frame(folds[i]))
      #   thinnedFolds[[i]] <- currentFold
      #   i <- i + 1
      #   if (i > (k + 1)) {
      #     break
      #   }
      # }
      # 
      # thinned <- as.data.frame(do.call(rbind, thinnedFolds))
      # thinned <- geoThinApprox(thinned, 1000, c('longitude', 'latitude'))
      # 
      # # finalThinned <- elimCellDups(thinned, raster::raster(), longLat = c('longitude', 'latitude'))
      # occsLongLat <- data.frame(longitude = occs$longitude, latitude = occs$latitude)
      finalThinned <- elimCellDups(occs, lorenzRast, longLat = ll)
      
      save(finalThinned, file = thinnedFileName)
    }
    
    # Before eliminating cell duplicates (X observations):
    # thinnedSp <- SpatialPointsDataFrame(occs[, ll], data = occs, 
    #                                     proj4 = getCRS('nad27', TRUE))
    # plot(thinnedSp, pch = 16, cex = 0.3, col = "red", 
    #      main = 'BIEN occurrences, before eliminating cell duplicates')
    # plot(range, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
    # map("state", add = TRUE)
    # map("world", add = TRUE)
    
    # Using sf for spatial visualization (easily customizable):
    speciesSf <- st_as_sf(x = finalThinned,
                          coords = c(x = 'longitude',
                                     y = 'latitude'),
                          crs = default_crs)
    speciesSfThin <- st_as_sf(x = occs,
                              coords = c(x = 'longitude',
                                         y = 'latitude'),
                              crs = default_crs)
    
    # load country/world basemap data
    world <- ne_countries(scale = "medium", returnclass = "sf")
    states <- st_as_sf(ne_states(country = 'united states of america'))
    canada <- st_as_sf(ne_states(country = 'canada'))
    
    ggplot(data = world) +
      theme_bw() +
      geom_sf(fill = "white") +
      geom_sf(data = states, fill = NA) + 
      geom_sf(data = canada, fill = NA) + 
      geom_sf(data = speciesSf, col = alpha('red'), cex = 0.5) +
      coord_sf(xlim = c(min(occs$longitude), max(occs$longitude)), 
               ylim = c(min(occs$latitude), max(occs$latitude)), expand = TRUE) +
      xlab("Longitude") + 
      ylab("Latitude") +
      ggtitle(paste0(sp, ' occurences'), subtitle = "(Cleaned and cell duplicated eliminated)") +
      theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                            size = 0.5), panel.background = element_rect(fill = "lavender"))
    
    
    bufferFileName <- paste0('./cleaning_records/', 
                             gsub(tolower(species), pattern = ' ', 
                                  replacement = '_'), '_buffer.rData')
    
    # if (file.exists(bufferFileName)) {
    #   load(bufferFileName)
    # } else {
    # choose a buffer that will include 97% (0.03 * # of occurrences) of the observations
    # you can change this value based on the threshold you want for your buffer
    threshold <- ceiling(0.03 * nrow(speciesSf))
    
    rangeMap <- st_transform(st_as_sf(x = range), getCRS('albersNA'))
    speciesSfThinAlb <- st_transform(speciesSfThin, getCRS('albersNA'))
    speciesSfAlb <- st_transform(speciesSf, getCRS('albersNA'))
    
    calculate_buffer <- function(x) {
      while (TRUE) {
        x <- x + 10
        buffer_distance_temp <- as_units(x, "km")
        range_buffer <- st_buffer(rangeMap, dist = buffer_distance_temp)
        
        speciesSf_filtered <- speciesSfAlb %>%
          mutate(within_range = lengths(st_within(x = speciesSfAlb,
                                                  y = rangeMap)),
                 within_buffer = lengths(st_within(x = speciesSfAlb,
                                                   y = range_buffer)))
        
        # These are for a sanity check when first running
        # Uncomment them if you want to see how the # of occurrences outside the buffer
        # decreases as the buffer distance increases
        
        # print(paste0("Buffer distance = ", x, " km"))
        # print(paste0("Number of occurrences outside buffer = ",
        #              length(which(speciesSf_filtered$within_buffer == 0))))
        
        if(length(which(speciesSf_filtered$within_buffer == 0)) < threshold) {
          if (x < 200) {
            print(paste0('Yes, x < 200, x = ', x))
            buffer_distance <<- as_units((x*2), "km")
            print(paste0('buffer distance = ', buffer_distance))
            range_buffer <<- st_buffer(rangeMap, dist = buffer_distance)
            
            speciesSf_filtered <<- speciesSfAlb %>%
              mutate(within_range = lengths(st_within(x = speciesSfAlb,
                                                      y = rangeMap)),
                     within_buffer = lengths(st_within(x = speciesSfAlb,
                                                       y = range_buffer)))
            thinnedSf_filtered <<- speciesSfThinAlb %>% 
              mutate(within_range = lengths(st_within(x = speciesSfThinAlb,
                                                      y = rangeMap)),
                     within_buffer = lengths(st_within(x = speciesSfThinAlb,
                                                       y = range_buffer)))
            
            save(buffer_distance, range_buffer, speciesSf_filtered, 
                 thinnedSf_filtered, file = bufferFileName)
          } else {
            buffer_distance <<- buffer_distance_temp
            thinnedSf_filtered <<- speciesSfThinAlb %>% 
              mutate(within_range = lengths(st_within(x = speciesSfThinAlb,
                                                      y = rangeMap)),
                     within_buffer = lengths(st_within(x = speciesSfThinAlb,
                                                       y = range_buffer)))
            
            save(buffer_distance, range_buffer, speciesSf_filtered, 
                 thinnedSf_filtered, file = bufferFileName)
          }
          break
        }
      }
      
      # return(buffer_distance)
    }
    
    x <- 0
    if(sp == 'Fraxinus pennsylvanica') calculate_buffer(250) else calculate_buffer(x)
    
    load(bufferFileName)
    
    # Thinned, before elimCellDups with calculated buffer.
    # Points that fall outside a buffer are colored red.
    ggplot(data = world) +
      theme_bw() +
      geom_sf(fill = "white") +
      geom_sf(data = states, fill = NA) + 
      geom_sf(data = canada, fill = NA) +
      geom_sf(data = rangeMap, 
              color = 'darkgreen',
              fill = 'green', 
              alpha = 0.4, 
              inherit.aes = FALSE) +
      geom_sf(data = range_buffer,
              color = 'yellow',
              fill = NA,
              inherit.aes = FALSE) +
      geom_sf(data = thinnedSf_filtered,
              size = 0.75, 
              aes(color = within_buffer > 0),
              inherit.aes = FALSE) +
      scale_colour_manual(values = setNames(c('black','red'),c(T, F)), guide = "none") +
      guides(fill = "none") +
      coord_sf(xlim = c(min(finalThinned$longitude) - 5, max(finalThinned$longitude) + 5), 
               ylim = c(min(finalThinned$latitude), max(finalThinned$latitude) + 5), expand = TRUE) +
      xlab("Longitude") + 
      ylab("Latitude") +
      ggtitle(paste0(sp, ' occurences'), 
              subtitle = paste0("(Cleaned, thinned, before cell duplicates eliminated, buffer = ",
                                buffer_distance, 
                                " km)")) +
      theme(panel.grid.major = element_line(color = gray(0.5), 
                                            linetype = "dashed", 
                                            size = 0.5), 
            panel.background = element_rect(fill = "lavender"),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    # Full range with calculated buffer and cell duplicates eliminated.
    ggplot(data = world) +
      theme_bw() +
      geom_sf(fill = "white") +
      geom_sf(data = states, fill = NA) + 
      geom_sf(data = canada, fill = NA) +
      geom_sf(data = rangeMap, 
              color = 'darkgreen',
              fill = 'green', 
              alpha = 0.4, 
              inherit.aes = FALSE) +
      geom_sf(data = range_buffer,
              color = 'yellow',
              fill = NA,
              inherit.aes = FALSE) +
      geom_sf(data = speciesSf_filtered,
              size = 0.75,
              aes(color = within_buffer > 0),
              inherit.aes = FALSE) +
      scale_colour_manual(values = setNames(c('black','red'),c(T, F)), guide = "none") +
      guides(fill = "none") +
      coord_sf(xlim = c(min(finalThinned$longitude) - 5, max(finalThinned$longitude) + 5), 
               ylim = c(min(finalThinned$latitude), max(finalThinned$latitude) + 5), expand = TRUE) +
      xlab("Longitude") +
      ylab("Latitude") +
      ggtitle(paste0(sp, ' occurences'), subtitle = paste0("(Cleaned, thinned, and duplicates eliminated, with buffer = ",
                                                           buffer_distance, " km)")) +
      theme(panel.grid.major = element_line(color = gray(0.5), 
                                            linetype = "dashed", 
                                            size = 0.5), 
            panel.background = element_rect(fill = "lavender"), 
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    speciesSf_filtered_final <- speciesSf %>%
      mutate(within_range = lengths(st_within(x = speciesSfAlb,
                                              y = rangeMap)),
             within_buffer = lengths(st_within(x = speciesSfAlb,
                                               y = range_buffer))) %>%
      filter(within_buffer > 0)
    
    # range_buffer_final <- range_buffer %>% filter(AREA == max(range_buffer$AREA))
    range_buffer_final <- st_union(range_buffer)
    
    cleanedRecordsFileName <- paste0('./cleaning_records/', 
                                     gsub(tolower(species), pattern = ' ', 
                                          replacement = '_'), '_finalRecords.rData')
    
    save(speciesSf_filtered_final, range_buffer_final, file = cleanedRecordsFileName)
    
    ggplot(data = world) +
      theme_bw() +
      geom_sf(fill = "white") +
      geom_sf(data = states, fill = NA) + 
      geom_sf(data = canada, fill = NA) +
      geom_sf(data = rangeMap, 
              color = 'darkgreen',
              fill = 'green', 
              alpha = 0.4, 
              inherit.aes = FALSE) +
      geom_sf(data = range_buffer_final,
              color = 'yellow',
              fill = NA,
              inherit.aes = FALSE) +
      geom_sf(data = speciesSf_filtered_final,
              size = 0.5,
              aes(color = within_range > 0),
              inherit.aes = FALSE) +
      scale_colour_manual(values = setNames(c('black','red'),c(T, F)), guide = "none") +
      guides(fill = "none") +
      coord_sf(xlim = c(min(finalThinned$longitude) - 5, 
                        max(finalThinned$longitude) + 5), 
               ylim = c(min(finalThinned$latitude), 
                        max(finalThinned$latitude) + 5), 
               expand = TRUE) +
      xlab("Longitude") +
      ylab("Latitude") +
      ggtitle(paste0(sp, ' occurences'), subtitle = paste0("(Cleaned, thinned, and duplicates eliminated, with buffer = ",
                                                           buffer_distance, " km)")) +
      theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                            size = 0.5), panel.background = element_rect(fill = "lavender"))
    
    # How many occurrences are in the final dataset?
    #   How many occurrences are outside of the buffer?
    print(paste0('Final number of occurrences = ', nrow(speciesSf_filtered_final)))
    print(paste0("Number of occurrences outside buffer = ", 
                 length(which(speciesSf_filtered$within_buffer == 0))))
    
    imageName <<- paste0("./workspaces/03 - Modeling Workspace - ", speciesAb, " Cleaning")
    save.image(imageName)
    # load(imageName)
  }
  
  for (gcm in gcmList) {
    gcm <- gcm
    print(paste0('GCM = ', gcm))
    
    # load PCA variables
    # load(paste0('./workspaces/PCA_', gcm, '_PC', pc))
    
    # Load environmental PCA rasters. 
    # If already clipped, load that. 
    # Otherwise, clip the data to present.
    studyRegionFileName <- './data_and_analyses/study_region/regions/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
    studyRegionRasts <- brick(studyRegionFileName)
    
    if (gcm == 'Beyer') { # Beyer
      load('./data_and_analyses/env_data/Beyer/PCA_clim.Rdata')
      load(paste0('./data_and_analyses/env_data/Beyer/pca_pc', pc, '.Rdata'))
      fileName <- paste0('./data_and_analyses/env_data/Beyer/envDataClipped_',
                         climYear, 'KYBP_pc', pc, '.tif')
      vars <- c("BIO1", paste0('BIO', 4:19), "cloudiness", "relative_humidity")
    } else if (gcm == 'Lorenz_ccsm') { # CCSM
      load(paste0('./data_and_analyses/env_data/Lorenz/PCA_', gcm, '_clim.Rdata')) 
      load(paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/pca_pc', pc, '.Rdata')) 
      fileName <- paste0('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/envDataClipped_',
                         climYear, 'KYBP_pc', pc, '.tif')
      clim <- lorenz
      workingFolder <- './data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ'
      vars <- sub('\\.tif.*', '', list.files(path = workingFolder, 
                                             pattern='*.tif', all.files = TRUE, full.names = FALSE))
      vars <- vars[lapply(vars, function(x) length(grep("pca_", x, value = F))) == 0]
    } else { # ECBilt
      load(paste0('./data_and_analyses/env_data/Lorenz/PCA_', gcm, '_clim.Rdata')) 
      load(paste0('./data_and_analyses/env_data/Lorenz/V2/ecbilt_21-0k_all_tifs_LJ/pca_pc', pc, '.Rdata'))
      fileName <- paste0('./data_and_analyses/env_data/Lorenz/V2/ecbilt_21-0k_all_tifs_LJ/envDataClipped_',
                         climYear, 'KYBP_pc', pc, '.tif')
      clim <- lorenz
      workingFolder <- paste0('./data_and_analyses/env_data/Lorenz/V2/',
                              gcm, '_21-0k_all_tifs_LJ')
      vars <- sub('\\.tif.*', '', list.files(path = workingFolder, 
                                             pattern='*.tif', all.files = TRUE, full.names = FALSE))
      vars <- vars[lapply(vars, function(x) length(grep("pca_", x, value = F))) == 0]
    }
    
    if (file.exists(fileName)) {
      envData <- brick(fileName)
      # rename raster layers to pc's
      names(envData) <- paste0('pca', 1:pc)
    } else {
      pcPrediction <<- list()
      # label each env layer by variable and year
      for (i in 1:length(clim)) {
        names(clim[[i]]) <- vars
        pcPrediction[i] <- raster::predict(clim[[i]], pca, index = 1:pc)
        names(pcPrediction[[i]]) <- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
      }
      
      envDataPca <- stack(pcPrediction)
      
      # keep only rasters (first five rasters) for climate year
      
      envYr <- pcPrediction[[(climYear/1000) + 1]]
      # plot(envYr) 
      names(envYr) <- paste0('pca', 1:pc)
      
      # check projections = wgs84
      print("Ensure that the projection of these rasters is WGS84:")
      print(paste0("Projection of envYr = ", projection(envYr)))
      
      # define study region & extent, load if already defined
      if (file.exists('./data_and_analyses/study_region/Study Region.Rdata') & 
          file.exists('./data_and_analyses/study_region/Study Region Extent.Rdata')) {
        load('./data_and_analyses/study_region/Study Region.Rdata')
        load('./data_and_analyses/study_region/Study Region Extent.Rdata')
      } else {
        studyRegion <- rgdal::readOGR('/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region',
                                      'study_region')
        studyRegion <- crop(studyRegion, extent(-178.2166, 83.7759, -55.90223, 83.6236))
        projection(studyRegion) <- getCRS("WGS84")
        studyExtent <- extent(studyRegion)
        save(studyRegion, 
             file='/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region/Study Region.Rdata',
             compress=TRUE)
        save(studyExtent, 
             file='/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/study_region/Study Region Extent.Rdata', compress=TRUE)
      }
      
      # clip environmental PCAs to study extent for given species, visualize, and save:
      envDataClipped <- list()
      for (i in 1:nlayers(envYr)) {
        x <- envYr[[i]]
        x <- crop(x, studyExtent)
        # x <<- mask(x, studyRegion)
        projection(x) <- getCRS("WGS84")
        envDataClipped[[i]] <- x
        envData <- stack(envDataClipped)
        # plot(envData)
        writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
      }
    } 
    
    fileName <- paste0('./data_and_analyses/cleaned_records/', 
                       species, '_finalRecords.rData')
    load(fileName)
    records <- data.frame(speciesSf_filtered_final)
    records$geometry <- gsub("[c()]", "", records$geometry)
    records <- separate(data = records, 
                        col = 'geometry', 
                        into = ll, 
                        sep = "\\,")
    records$longitude <- as.double(records$longitude)
    records$latitude <- as.double(records$latitude)
    
    occsEnv <- raster::extract(envData, 
                               cbind(records$longitude, 
                                     records$latitude))
    occsEnvDf <- as.data.frame(occsEnv)
    records <- cbind(records, occsEnvDf)
    
    # records in the water:
    if (any(is.na(rowSums(occsEnvDf)))) {
      water <- records[which(is.na(rowSums(occsEnvDf))), ] 
      water <- SpatialPointsDataFrame(water[,ll], data = water,
                                      proj4 = getCRS('wgs84', TRUE))
    }
    
    # remove records in water from dataset:
    if (any(is.na(rowSums(occsEnvDf)))) records <-
      records[-which(is.na(rowSums(occsEnvDf))), ]
    
    # visualize the points that fall in the water
    # convert to sp object
    recordsSp <- SpatialPointsDataFrame(records[, ll], data = records,
                                        proj4 = getCRS('wgs84', TRUE))
    
    # visualize
    plot(recordsSp, pch = 16, cex = 0.5, col = "red", 
         main = paste0(sp, ' occurrences (BIEN) thinned'))
    if (exists("water")) {
      plot(water, col = 'blue', add = TRUE)
    }
    map("state", add = TRUE)
    map("world", add = TRUE)
    
    save.image(paste0('./workspaces/04 - Modeling Workspace - Clipping ', sp, '_PC_', pc, '_GCM_', gcm))
    # load(paste0('./workspaces/04 - Modeling Workspace - Clipping ', sp, '_PC_', pc, '_GCM_', gcm))
    
    # load(paste0('./workspaces/03 - Modeling Workspace - ', speciesAb, ' Cleaning'))
    
    bufferFileName <- paste0('./data_and_analyses/cleaned_records/', species, '_buffer.rData')
    load(bufferFileName)
    # load(paste0("./workspaces/03 - Modeling Workspace - ", speciesAb, " Cleaning"))
    
    # calculate calibration region buffer at 160-km to extract bg sites
    calibBuffer <- st_buffer(st_transform(st_as_sf(x = recordsSp), getCRS('albersNA')), 
                             dist = as_units(160, 'km'))
    calibBuffer <- st_union(calibBuffer)
    
    calibRegionSpAlb <- as(calibBuffer, 'Spatial')
    calibRegionSpAlb <- sp::spTransform(calibRegionSpAlb, getCRS('albersNA', TRUE))
    calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
    
    bgFileName <- paste0('./Background Sites/Random Background Sites across Study Region - ', 
                         speciesAb, '_', gcm, '.Rdata')
    bgRawFileName <- paste0('./Background Sites/raw/Random Bg Sites - ', speciesAb_, '.Rdata')
    # get 20,000 random background sites from calibration region
    # will only keep 10,000 points...this accounts for points that may fall in water
    if (file.exists(bgRawFileName)) load(bgRawFileName) else {
      bgTestSpAlb <- suppressWarnings(sp::spsample(calibRegionSpAlb, n=20000, type='random'))
      save(bgTestSpAlb, file = bgRawFileName)
    }
    bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
    randomBgSites <- dismo::randomPoints(envData, 20000)
    # randomBgSites <- as.data.frame(coordinates(bgTestSp))
    names(randomBgSites) <- ll
    
    # sanity check: plot the background sites
    plot(bgTestSp, pch = 16, cex = 0.5, col = "red", 
         main = paste0(sp, ' background sites'))
    plot(calibRegionSpWgs, add = TRUE, border = 'blue')
    map("state", add = TRUE)
    map("world", add = TRUE)
    
    # extract environment at background sites
    climate <- envData
    randomBgEnv <- raster::extract(climate, randomBgSites)
    randomBgEnv <- as.data.frame(randomBgEnv)
    
    # remove any sites with NA for at least one variable
    isNa <- is.na(rowSums(randomBgEnv))
    if (any(isNa)) {
      randomBgSites <- randomBgSites[-which(isNa), ]
      randomBgEnv <- randomBgEnv[-which(isNa), ]
    }
    
    # only keep 10,000 random background sites
    randomBgSites <- randomBgSites[1:10000,]
    randomBgEnv <- randomBgEnv[1:10000,]
    
    # combine with coordinates and rename coordinate fields
    randomBg <- cbind(randomBgSites, randomBgEnv)
    names(randomBg)[1:2] <- ll
    
    # dir.create('./Background Sites', recursive=TRUE, showWarnings=FALSE)
    save(randomBg, randomBgEnv,
         file = bgFileName, compress=TRUE)
    
    presBg <- c(rep(1, nrow(records)), rep(0, nrow(randomBgEnv)))
    occsEnv <- occsEnv[complete.cases(occsEnv), ]
    
    env <- rbind(occsEnv, randomBgEnv)
    env <- cbind(presBg, env)
    env <- as.data.frame(env)
    
    env <- env[complete.cases(env), ]
    # View(env)
    
    # model species
    # try: trainMaxNet(data, resp='presBg', regMult=c(0.5, 1, 2, 3, 4, 5, 7.5, 10))
    envModel <- enmSdm::trainMaxNet(data = env, resp = 'presBg', out = c('models', 'tuning'))
    # envModel <- maxnet(p = presBg, data = trainData)
    modelFileName <- paste0('./Models/Maxent/', speciesAb_, '_Maxent/Model_PC', pc, '_GCM_', gcm, '.Rdata')
    save(envModel, file = modelFileName, compress=TRUE)
    
    predictors <- c(paste0('pca', 1:pc))
    # prediction
    # envMap <- predict(climate[[predictors]], envModel, clamp = F, type = 'cloglog')
    envMap <- predict(
      climate[[predictors]],
      envModel,
      filename = paste0('./Models/Maxent/', speciesAb_, '_Maxent/prediction_PC',
                        pc, '_GCM', gcm, '_', climYear, 'ybp'),
      clamp = F,
      format='GTiff',
      overwrite = TRUE,
      type='cloglog')
    
    envMapSp <- rasterToPolygons(envMap)
    
    plot(rangeMap, border = 'blue', main = paste0('Maxent output, ', 
                                                            sp,
                                                            ' occurrences'))
    plot(envMap, add = TRUE)
    plot(rangeMap, border = 'blue', add = TRUE)
    map("state", add = TRUE)
    map("world", add = TRUE)
    points(records$longitude, records$latitude, pch = 16, cex = 0.6, col = 'red')
    
    plot(envMap, main = paste0('Maxent output, ', 
                                         sp,
                                         ' occurrences'))
    plot(rangeMap, border = 'blue', add = TRUE)
    
    outputFileName <<- paste0('./Models/Maxent/', speciesAb_, 
                              '_Maxent/Model_PC', pc, '_GCM_', gcm, '.rData')
    save(rangeMap, envMap, envModel, records, file = outputFileName, overwrite = T)
    
    outputFileName <<- paste0('./Models/Maxent/model_outputs/', speciesAb_,
                              '_GCM', gcm, '_PC', pc, '.Rdata')
    save(rangeMap, envMap, envModel, records, file = outputFileName, overwrite = T)
    
    save.image(paste0('./workspaces/05 - Modeling Workspace - ', speciesAb_,
                      ' Model Output - PC', pc, '_GCM_', gcm))
  }
}

pdf(file = paste0('./cleaning/cleaned_records.pdf'), width = 11, height = 8.5)
for (sp in speciesList) {
  sp <- sp
  # print(paste0("SPECIES = ", sp))
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  load(paste0("./workspaces/03 - Modeling Workspace - ", speciesAb, " Cleaning"))
  
  p <- ggplot(data = world) +
    theme_bw() +
    geom_sf(fill = "white") +
    geom_sf(data = states, fill = NA) + 
    geom_sf(data = canada, fill = NA) +
    geom_sf(data = rangeMap, 
            color = 'darkgreen',
            fill = 'green', 
            alpha = 0.4, 
            inherit.aes = FALSE) +
    geom_sf(data = range_buffer_final,
            color = 'yellow',
            fill = NA,
            inherit.aes = FALSE) +
    geom_sf(data = speciesSf_filtered_final,
            size = 0.5,
            aes(color = within_range > 0),
            inherit.aes = FALSE) +
    scale_colour_manual(values = setNames(c('black','red'),c(T, F)), guide = "none") +
    guides(fill = "none") +
    coord_sf(xlim = c(min(finalThinned$longitude) - 5, 
                      max(finalThinned$longitude) + 5), 
             ylim = c(min(finalThinned$latitude), 
                      max(finalThinned$latitude) + 5), 
             expand = TRUE) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle(paste0(sp, ' occurences'), subtitle = paste0("(Cleaned, thinned, and duplicates eliminated, with buffer = ",
                                                         buffer_distance, " km)")) +
    theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                          size = 0.5), panel.background = element_rect(fill = "lavender"))
  print(p)
}
dev.off()

# zoomed into range
speciesList <- c('Fraxinus americana','Fraxinus caroliniana', 'Fraxinus cuspidata', 
                 'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica', 
                 'Fraxinus profunda', 'Fraxinus quadrangulata')
gcmList <- c('Beyer', 'Lorenz_ccsm', 'ecbilt')

pdf(file = './Models/Maxent/all_models.pdf', height = 8.5, width = 11)
for(sp in speciesList) {
  sp <- sp
  # print(paste0("SPECIES = ", sp))
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  for(gcm in gcmList) {
    load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
                '_PC', pc, '.rData'))
    plot(rangeMap, border = 'blue', 
         main = paste0('Maxent output, ', sp, ' occurrences,', '\nGCM = ', gcm))
    plot(envMap, add = TRUE)
    plot(rangeMap, border = 'blue', add = TRUE)
    map("state", add = TRUE)
    map("world", add = TRUE)
    points(records$longitude, records$latitude, pch = 16, cex = 0.6, col = 'red')
  }
}
dev.off()

# zoomed out
pdf(file = './Models/Maxent/all_models_zoomedOut.pdf', height = 8.5, width = 11)
for(sp in speciesList) {
  sp <- sp
  # print(paste0("SPECIES = ", sp))
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  for(gcm in gcmList) {
    load(paste0('./Models/Maxent/model_outputs/', speciesAb_, '_GCM', gcm, 
                '_PC', pc, '.rData'))
    plot(envMap, main = paste0('Maxent output, ', sp, ' occurrences,', '\nGCM = ', gcm))
    plot(rangeMap, border = 'blue', add = TRUE)
    points(records$longitude, records$latitude, pch = 16, cex = 0.3, col = 'red')
    
  }
}
dev.off()
