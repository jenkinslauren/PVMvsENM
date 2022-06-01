# Cleaning template
# Author: Lauren Jenkins
# 25 January 2022

rm(list = ls())
# Load required packages.
# Load packages
library(enmSdm)

# for handling rasters
library(raster)
# library(rgeos) # not used, but may be helpful when debugging

# visualization tools
library(sp)
library(sf)
library(maps)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# additional tools
library(tools)
library(units)
library(dplyr)

setwd('/Volumes/lj_mac_22/MOBOT/PVMvsENM')

# nad27
default_crs = sf::st_crs(4267)

# setwd('/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/')
load('./workspaces/01 - Modeling Workspace - Fraxinus Range Maps (BIEN + Little)')

# set constants
ll <- c('longitude', 'latitude')
speciesList <- c('Fraxinus americana', 'Fraxinus caroliniana', 'Fraxinus cuspidata'
                 ,'Fraxinus greggii', 'Fraxinus nigra', 'Fraxinus pennsylvanica',
                 'Fraxinus profunda', 'Fraxinus quadrangulata')

lorenzRast <- raster::raster('./data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/0BP/an_avg_ETR.tif')

occs <- data.frame()
for(sp in speciesList) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesFileName <- paste0('./data_and_analyses/species_records/01_', 
                            gsub(tolower(species), pattern = ' ', 
                                 replacement = '_'), '_retained_records.rds')
  
  occs_current <- readRDS(speciesFileName)
  occs <- rbind(occs, occs_current)
}

clean <- function(sp) {
  species <- gsub(tolower(sp), pattern=' ', replacement='_')
  speciesAb <- paste0(substr(sp,1,4), toupper(substr(sp,10,10)), substr(sp,11,13))
  speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", speciesAb)
  
  if (sp == 'Fraxinus all') {
    range <- rgeos::gUnion(rgeos::gUnion(rgeos::gUnion(littleRange_FraxAmer, littleRange_FraxCaro), 
                           rgeos::gUnion(littleRange_FraxCusp, littleRange_FraxGreg)),
                           rgeos::gUnion(rgeos::gUnion(littleRange_FraxNigr, littleRange_FraxPenn), 
                           rgeos::gUnion(littleRange_FraxProf, littleRange_FraxQuad)))
  } else {
    rangeName <- paste0('littleRange_', speciesAb)
    range <- get(rangeName)
  }

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
  thinnedFileName <- paste0('./data_and_analyses/species_records/02_', 
                            gsub(tolower(species), pattern = ' ', 
                                 replacement = '_'), '_thinned_records.rData')
  
  if (file.exists(thinnedFileName)) {
    load(thinnedFileName)
  } else {
    # # memory.limit() - will give you something with memory (how much is remaining or it's capacity)
    # # lsos might also help with identifying how much memory
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
    # finalThinned <- elimCellDups(thinned, raster::raster(), longLat = c('longitude', 'latitude'))
    
    finalThinned <- elimCellDups(occs, lorenzRast, longLat = ll)
    
    save(finalThinned, file = thinnedFileName)
  }
  
  # Before eliminating cell duplicates (X observations):
  thinnedSp <- SpatialPointsDataFrame(thinned[, ll], data = thinned, 
                                      proj4 = getCRS('nad27', TRUE))
  plot(thinnedSp, pch = 16, cex = 0.3, col = "red", 
       main = 'BIEN occurrences thinned, before eliminating cell duplicates')
  plot(range, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
  map("state", add = TRUE)
  map("world", add = TRUE)
  
  thinned <- finalThinned
  # Using sf for spatial visualization (easily customizable):
  speciesSf <- st_as_sf(x = finalThinned,
                        coords = c(x = 'longitude',
                                   y = 'latitude'),
                        crs = default_crs)
  speciesSfThin <- st_as_sf(x = thinned,
                            coords = c(x = 'longitude',
                                       y = 'latitude'),
                            crs = default_crs)
  
  # load country/world basemap data
  world <- ne_countries(scale = "medium", returnclass = "sf")
  states <- ne_states(returnclass = "sf")
  canada <- ne_states(country = 'canada', returnclass = "sf")
  
  ggplot(data = world) +
    theme_bw() +
    geom_sf(fill = "white") +
    geom_sf(data = states, fill = NA) + 
    geom_sf(data = canada, fill = NA) + 
    geom_sf(data = speciesSf, col = alpha('red'), cex = 0.5) +
    coord_sf(xlim = c(min(thinned$longitude), max(thinned$longitude)), 
             ylim = c(min(thinned$latitude), max(thinned$latitude)), expand = TRUE) +
    xlab("Longitude") + 
    ylab("Latitude") +
    ggtitle(paste0(sp, ' occurences'), subtitle = "(Cleaned and thinned)") +
    theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                          size = 0.5), panel.background = element_rect(fill = "lavender"))
  
  
  # Set buffer.
  bufferFileName <- paste0('./data_and_analyses/cleaned_records/', 
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
      buffer_distance <- as_units(x, "km")
      range_buffer <- st_buffer(rangeMap, dist = buffer_distance)
      
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
        if (x < 300) {
          x <<- 300
          rm(buffer_distance)
          rm(range_buffer)
          buffer_distance <<- as_units(x, "km")
          range_buffer <<- st_buffer(rangeMap, dist = buffer_distance)
          
          speciesSf_filtered <<- speciesSfAlb %>%
            mutate(within_range = lengths(st_within(x = speciesSfAlb,
                                                    y = rangeMap)),
                   within_buffer = lengths(st_within(x = speciesSfAlb,
                                                     y = range_buffer)))
        }
        thinnedSf_filtered <<- speciesSfThinAlb %>% 
          mutate(within_range = lengths(st_within(x = speciesSfThinAlb,
                                                  y = rangeMap)),
                 within_buffer = lengths(st_within(x = speciesSfThinAlb,
                                                   y = range_buffer)))
        
        save(buffer_distance, range_buffer, speciesSf_filtered, 
             thinnedSf_filtered, file = bufferFileName)
        break
      }
    }
    
    # return(buffer_distance)
  }
  
  x <- 0
  calculate_buffer(x)
  # }
  # load(bufferFileName)
  
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
    coord_sf(xlim = c(min(thinned$longitude) - 5, max(thinned$longitude) + 5), 
             ylim = c(min(thinned$latitude), max(thinned$latitude) + 5), expand = TRUE) +
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
    coord_sf(xlim = c(min(thinned$longitude) - 5, max(thinned$longitude) + 5), 
             ylim = c(min(thinned$latitude), max(thinned$latitude) + 5), expand = TRUE) +
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
  
  cleanedRecordsFileName <- paste0('./data_and_analyses/cleaned_records/', 
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

lapply(speciesList, clean)
