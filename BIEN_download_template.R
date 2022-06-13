# Downloading species occurrence data for any genus.
# Author, Lauren Jenkins
# 13 June 2022
# Last updated: 13 June 2022

rm(list=ls())

library(BIEN)

## constants for a particular genus ##
genus <- 'fagus'
speciesList <- c('Fagus grandifolia')

setwd(paste0('/Volumes/lj_mac_22/MOBOT/by_genus/', genus))

dir.create('./species_records')

for(sp in speciesList) {
  
  speciesFileName <- paste0('./species_records/00_', gsub(' ', '_', tolower(sp)), 
                            '_bien_all_occurrences.rda')
  
  occsRaw <- BIEN_occurrence_species(
    species = sp,
    cultivated = F,
    all.taxonomy = F,
    native.status = F,
    natives.only = T,
    observation.type = T,
    political.boundaries = T,
    collection.info = T
  )
  
  save(occsRaw, file = speciesFileName)
  
}

sink('./species_records/bien_download.txt')
print(paste0('Data representing species records were downloaded from BIEN on ', date(), '.'))
print(BIEN_metadata_database_version())
sink()
