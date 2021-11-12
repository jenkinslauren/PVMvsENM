# Load all of BIEN Fraxinus data
# Lauren Jenkins, Fall 2021

### setup ###
## cleaning ##

### setup ###
rm(list = ls())
setwd('/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/')

load('./00 - Modeling Workspace - Fraxinus Range Maps (Little)')

ll <- c('longitude', 'latitude')

species <- c('Fraxinus americana' ,'Fraxinus anomala' ,'Fraxinus berlandieriana' ,
             'Fraxinus caroliniana' ,'Fraxinus cuspidata' ,'Fraxinus dipetala', 
             'Fraxinus gooddingii' ,'Fraxinus greggii', 'Fraxinus latifolia', 
             'Fraxinus nigra', 'Fraxinus papillosa', 'Fraxinus parryi', 
             'Fraxinus pennsylvanica', 'Fraxinus profunda', 'Fraxinus quadrangulata', 
             'Fraxinus velutina')

load_BIEN <- function(x) {
    BIEN_ranges_species(x, directory = './regions/bien_range_map')
    bienRange <- BIEN_ranges_load_species(x)
    bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA'), TRUE)
    
    # name file by when it was obtained
    speciesFileName <- paste0('./species_records/00_', 
                              gsub(tolower(x), pattern = ' ', 
                                   replacement = '_'), '_bien_all_occurrences.rda')
    
    if (file.exists(speciesFileName)) {
      load(speciesFileName)
    } else {
      occsRaw <- BIEN_occurrence_species(species = x, cultivated = FALSE,
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
}

lapply(species, load_BIEN)

load_Occs <- function(x) {
  BIEN_ranges_species(x, directory = './regions/bien_range_map')
  bienRange <- BIEN_ranges_load_species(x)
  bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA'), TRUE)
  
  # name file by when it was obtained
  speciesFileName <- paste0('./species_records/00_', 
                            gsub(tolower(x), pattern = ' ', 
                                 replacement = '_'), '_bien_all_occurrences.rda')
  load(speciesFileName)
  print(paste0("Species: ", x))
  print(paste0("# of occurrences = ", nrow(occsRaw)))
}

lapply(species, load_Occs)

## [1] "Species: Fraxinus americana"
# [1] "# of occurrences = 97302"
# [1] "Species: Fraxinus anomala"
# [1] "# of occurrences = 512"
# [1] "Species: Fraxinus berlandieriana"
# [1] "# of occurrences = 9"
# [1] "Species: Fraxinus caroliniana"
# [1] "# of occurrences = 2287"
# [1] "Species: Fraxinus cuspidata"
# [1] "# of occurrences = 193"
# [1] "Species: Fraxinus dipetala"
# [1] "# of occurrences = 283"
# [1] "Species: Fraxinus gooddingii"
# [1] "# of occurrences = 72"
# [1] "Species: Fraxinus greggii"
# [1] "# of occurrences = 154"
# [1] "Species: Fraxinus latifolia"
# [1] "# of occurrences = 749"
# [1] "Species: Fraxinus nigra"
# [1] "# of occurrences = 48624"
# [1] "Species: Fraxinus papillosa"
# [1] "# of occurrences = 30"
# [1] "Species: Fraxinus parryi"
# [1] "# of occurrences = 0"
# [1] "Species: Fraxinus pennsylvanica"
# [1] "# of occurrences = 93206"
# [1] "Species: Fraxinus profunda"
# [1] "# of occurrences = 1732"
# [1] "Species: Fraxinus quadrangulata"
# [1] "# of occurrences = 1336"
# [1] "Species: Fraxinus velutina"
# [1] "# of occurrences = 955"

# BIEN range maps
# Fraxinus americana
bienRange_FraxAmer <- shapefile('./regions/bien_range_map/Fraxinus_americana.prj')

# Fraxinus anomala
bienRange_FraxAnom <- shapefile('./regions/bien_range_map/Fraxinus_anomala.prj')

# Fraxinus berlandieriana
bienRange_FraxBerl <- shapefile('./regions/bien_range_map/Fraxinus_berlandieriana.prj')

# Fraxinus caroliniana
bienRange_FraxCaro <- shapefile('./regions/bien_range_map/Fraxinus_caroliniana.prj')

# Fraxinus cuspidata
bienRange_FraxCusp <- shapefile('./regions/bien_range_map/Fraxinus_cuspidata.prj')

# Fraxinus dipetala
bienRange_FraxDipe <- shapefile('./regions/bien_range_map/Fraxinus_dipetala.prj')

# Fraxinus gooddingii
bienRange_FraxGood <- shapefile('./regions/bien_range_map/Fraxinus_gooddingii.prj')

# Fraxinus greggii
bienRange_FraxGreg <- shapefile('./regions/bien_range_map/Fraxinus_greggii.prj')

# Fraxinus latifolia
bienRange_FraxLati <- shapefile('./regions/bien_range_map/Fraxinus_latifolia.prj')

# Fraxinus nigra
bienRange_FraxNigr <- shapefile('./regions/bien_range_map/Fraxinus_nigra.prj')

# Fraxinus papillosa
bienRange_FraxPapi <- shapefile('./regions/bien_range_map/Fraxinus_papillosa.prj')

# Fraxinus pennsylvanica
bienRange_FraxPenn <- shapefile('./regions/bien_range_map/Fraxinus_pennsylvanica.prj')

# Fraxinus profunda
bienRange_FraxProf <- shapefile('./regions/bien_range_map/Fraxinus_profunda.prj')

# Fraxinus quadrangulata
bienRange_FraxQuad <- shapefile('./regions/bien_range_map/Fraxinus_quadrangulata.prj')

# Fraxinus velutina
bienRange_FraxVelu <- shapefile('./regions/bien_range_map/Fraxinus_velutina.prj')

bienRange <- rbind(bienRange_FraxAmer, bienRange_FraxAnom)
bienRange <- rbind(bienRange, bienRange_FraxBerl)
bienRange <- rbind(bienRange, bienRange_FraxCaro)
bienRange <- rbind(bienRange, bienRange_FraxCusp)
bienRange <- rbind(bienRange, bienRange_FraxDipe)
bienRange <- rbind(bienRange, bienRange_FraxGood)
bienRange <- rbind(bienRange, bienRange_FraxGreg)
bienRange <- rbind(bienRange, bienRange_FraxLati)
bienRange <- rbind(bienRange, bienRange_FraxNigr)
bienRange <- rbind(bienRange, bienRange_FraxPapi)
bienRange <- rbind(bienRange, bienRange_FraxPenn)
bienRange <- rbind(bienRange, bienRange_FraxProf)
bienRange <- rbind(bienRange, bienRange_FraxQuad)
bienRange <- rbind(bienRange, bienRange_FraxVelu)

# Load entire genus
# This takes a long time, so be prepared
# occsRaw <- BIEN_occurrence_genus('Fraxinus', cultivated = FALSE,
                      # only.new.world = TRUE, all.taxonomy = FALSE,
                      # native.status = FALSE, natives.only = TRUE,
                      # observation.type = TRUE, political.boundaries = TRUE,
                      # collection.info = TRUE)
# sink('./species_records/bien_download.txt')
# say('Data representing species records were downloaded from BIEN on ', date(), '.')
# say('Bien Version 4.2')
# sink()

# remove anything that is outside the range we are looking at (US, Canada, Mexico)
# Fraxinus <- Fraxinus_load_genus[Fraxinus_load_genus$country %in% countries, ]

# save(Fraxinus_entire_genus, file = speciesFileName)
# load("./species_records/00_Fraxinus_bien_all_occurrences.rda")

# save workspace
save.image('./01 - Modeling Workspace - Fraxinus Range Maps (BIEN)')
