# Fraxinus
# Lauren Jenkins, Fall 2021

# List of working species: 
# Fraxinus americana, Fraxinus anomala, Fraxinus berlandieriana,
# Fraxinus caroliniana, Fraxinus cuspidata, Fraxinus dipetala, 
# Fraxinus gooddingii, Fraxinus greggii, Fraxinus latifolia, 
# Fraxinus nigra, Fraxinus papillosa, Fraxinus parryi, 
# Fraxinus pennsylvanica, Fraxinus profunda, Fraxinus quadrangulata, 
# Fraxinus velutina

## setup
rm(list = ls())
setwd('/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/')
ll <- c('longitude', 'latitude')

## obtaining data

# Load Little ranges
# Fraxinus americana
littleFraxAmer <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_americana/fraxamer.shp'
littleRange_FraxAmer <- shapefile(littleFraxAmer)
projection(littleRange_FraxAmer) <- enmSdm::getCRS('nad27')
littleRange_FraxAmer <- littleRang_eFraxAmer[littleRange_FraxAmer$CODE == 1, ] # remove holes

# Fraxinus anomala
littleFraxAnom <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_anomala/fraxanom.shp'
littleRange_FraxAnom <- shapefile(littleFraxAnom)
projection(littleRange_FraxAnom) <- enmSdm::getCRS('nad27')
littleRange_FraxAnom <- littleRange_FraxAnom[littleRange_FraxAnom$CODE == 1, ] # remove holes

# Fraxinus berlandieriana
littleFraxBerl <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_berlandieriana/fraxberl.shp'
littleRange_FraxBerl <- shapefile(littleFraxBerl)
projection(littleRangeFraxBerl) <- enmSdm::getCRS('nad27')
littleRange_FraxBerl <- littleRange_FraxBerl[littleRange_FraxBerl$CODE == 1, ] # remove holes

# Fraxinus caroliniana
littleFraxCaro <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_caroliniana/fraxcaro.shp'
littleRange_FraxCaro <- shapefile(littleFraxCaro)
projection(littleRange_FraxCaro) <- enmSdm::getCRS('nad27')
littleRange_FraxCaro <- littleRange_FraxCaro[littleRange_FraxCaro$CODE == 1, ] # remove holes

# Fraxinus cuspidata
littleFraxCusp <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_cuspidata/fraxcusp.shp'
littleRange_FraxCusp <- shapefile(littleFraxCusp)
projection(littleRange_FraxCusp) <- enmSdm::getCRS('nad27')
littleRange_FraxCusp <- littleRange_FraxCusp[littleRange_FraxCusp$CODE == 1, ] # remove holes

# Fraxinus dipetala
littleFraxDipe <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_dipetala/fraxdipe.shp'
littleRange_FraxDipe <- shapefile(littleFraxDipe)
projection(littleRange_FraxDipe) <- enmSdm::getCRS('nad27')
littleRange_FraxDipe <- littleRange_FraxDipe[littleRange_FraxDipe$CODE == 1, ] # remove holes

# Fraxinus gooddingii
littleFraxGood <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_gooddingii/fraxgood.shp'
littleRange_FraxGood <- shapefile(littleFraxGood)
projection(littleRange_FraxGood) <- enmSdm::getCRS('nad27')
littleRange_FraxGood <- littleRange_FraxGood[littleRange_FraxGood$CODE == 1, ] # remove holes

# Fraxinus greggii
littleFraxGreg <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_greggii/fraxgreg.shp'
littleRange_FraxGreg <- shapefile(littleFraxGreg)
projection(littleRange_FraxGreg) <- enmSdm::getCRS('nad27')
littleRange_FraxGreg <- littleRange_FraxGreg[littleRange_FraxGreg$CODE == 1, ] # remove holes

# Fraxinus latifolia
littleFraxLati <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_latifolia/fraxlati.shp'
littleRange_FraxLati <- shapefile(littleFraxLati)
projection(littleRange_FraxLati) <- enmSdm::getCRS('nad27')
littleRange_FraxLati <- littleRange_FraxLati[littleRange_FraxLati$CODE == 1, ] # remove holes

# Fraxinus nigra
littleFraxNigr <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_nigra/fraxnigr.shp'
littleRange_FraxNigr <- shapefile(littleFraxNigr)
projection(littleRange_FraxNigr) <- enmSdm::getCRS('nad27')
littleRange_FraxNigr <- littleRange_FraxNigr[littleRange_FraxNigr$CODE == 1, ] # remove holes

# Fraxinus papillosa
littleFraxPapi <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_papillosa/fraxpapi.shp'
littleRange_FraxPapi <- shapefile(littleFraxPapi)
projection(littleRange_FraxPapi) <- enmSdm::getCRS('nad27')
littleRange_FraxPapi <- littleRange_FraxPapi[littleRange_FraxPapi$CODE == 1, ] # remove holes

# Fraxinus pennsylvanica
littleFraxPenn <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_pennsylvanica/fraxpenn.shp'
littleRange_FraxPenn <- shapefile(littleFraxPenn)
projection(littleRange_FraxPenn) <- enmSdm::getCRS('nad27')
littleRange_FraxPenn <- littleRange_FraxPenn[littleRange_FraxPenn$CODE == 1, ] # remove holes

# Fraxinus profunda
littleFraxProf <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_profunda/fraxprof.shp'
littleRange_FraxProf <- shapefile(littleFraxProf)
projection(littleRange_FraxProf) <- enmSdm::getCRS('nad27')
littleRange_FraxProf <- littleRange_FraxProf[littleRange_FraxProf$CODE == 1, ] # remove holes

# Fraxinus quadrangulata
littleFraxQuad <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_quadrangulata/fraxquad.shp'
littleRange_FraxQuad <- shapefile(littleFraxQuad)
projection(littleRange_FraxQuad) <- enmSdm::getCRS('nad27')
littleRange_FraxQuad <- littleRange_FraxQuad[littleRange_FraxQuad$CODE == 1, ] # remove holes

# Fraxinus velutina
littleFraxVelu <- '/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/range_maps/little/fraxinus_velutina/fraxvelu.shp'
littleRange_FraxVelu <- shapefile(littleFraxVelu)
projection(littleRange_FraxVelu) <- enmSdm::getCRS('nad27')
littleRange_FraxVelu <- littleRange_FraxVelu[littleRange_FraxVelu$CODE == 1, ] # remove holes

littleRange <- rgeos::gUnion(littleRange_FraxAmer, littleRange_FraxAnom)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxBerl)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxCaro)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxCusp)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxDipe)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxGood)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxGreg)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxLati)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxNigr)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxPapi)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxPenn)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxProf)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxQuad)
littleRange <- rgeos::gUnion(littleRange, littleRange_FraxVelu)

plot(littleRange, col = alpha('green', 0.4), border = 'green', main = 'Little range map')
map("state", add = TRUE)
map("world", add = TRUE)

# save workspace
save.image('./00 - Modeling Workspace - Fraxinus Range Maps (Little)')
saveRDS(littleRange, './littleMergedRangeFraxinus.rds')

# load BIEN ranges
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

plot(bienRange, col = alpha('green', 0.4), border = 'green', main = 'BIEN range map')
map("state", add = TRUE)
map("world", add = TRUE)

# plot side by side
par(mfrow=c(1,2))
plot(bienRange, col = alpha('blue', 0.4), border = 'blue', main = 'BIEN range map')
map("state", add = TRUE)
map("world", add = TRUE)

plot(littleRange, col = alpha('green', 0.4), border = 'green', main = 'Little range map')
map("state", add = TRUE)
map("world", add = TRUE)

# save workspace
save.image('./01 - Modeling Workspace - Fraxinus Range Maps (BIEN + Little)')

## cleaning ##
# Fraxinus americana
load('./species_records/00_fraxinus_americana_bien_all_occurrences.rda')

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 97,136 observations

# remove records with '[no precise location]' locality
occs <- occs[!grepl("no precise location", occs$locality), ]
# 97,113 observations

# remove duplicate records
# can't have 2 columns with the same label
names(occs)[18] <- 'date_collected_2'
occs <- occs %>% distinct()
# 82,888 observations

# looking at year distributions
## occs$Year <- substr(occs$date_collected, 1, 4)
## occsBefore2000 <- occs[occs$Year <= 2000, ]
## occsBefore1980 <- occs[occs$Year <= 1980, ]
## occsBefore1970 <- occs[occs$Year <= 1970, ]
## hist(as.numeric(occsBefore1970$Year), xlab = 'Collection year', 
   ##  ylab = 'Number of records', main = 'Collection year')


occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'blue', main = 'BIEN range map')
map("state", add = TRUE)
map("world", add = TRUE)

plot(littleRange_FraxAmer, col = alpha('green', 0.4), border = 'green', add = TRUE)

# look at Kansas
## kansas <- occs[which(occs$state_province == 'Kansas'), ]

# look at Mexico
## mexico <- occs[which(occs$country == 'Mexico'), ]

# Fraxinus anomala
load('./species_records/00_fraxinus_anomala_bien_all_occurrences.rda')
# 512 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 490 occurrences

# 2 weird occurrences in Mexico
# 1 weird occurrence in central Colorado
plot(occsSp, pch = 16, cex = 0.5, col = 'green', main = 'BIEN range map')
map("state", add = TRUE)
# lineInSand <- click(n=1)
# lineInSand
# occs[which(occs$longitude > lineInSand[1]), ]
## record 31
## "hill above junction before Delores river"...not where the lat/long are

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'green', main = 'BIEN range map')
plot(bienRange_FraxAnom, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus berlandieriana
load('./species_records/00_fraxinus_berlandieriana_bien_all_occurrences.rda')
# 9 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 6 occurrences

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxBerl, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus caroliniana
load('./species_records/00_fraxinus_caroliniana_bien_all_occurrences.rda')
# 2287 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 2279 occurrences

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxCaro, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus cuspidata
load('./species_records/00_fraxinus_cuspidata_bien_all_occurrences.rda')
# 193 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 188 occurrences

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxCusp, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus dipetala
load('./species_records/00_fraxinus_dipetala_bien_all_occurrences.rda')
# 283 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 277 occurrences

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxDipe, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus gooddingii
load('./species_records/00_fraxinus_gooddingii_bien_all_occurrences.rda')
# 72 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 69 occurrences

# lineInSand <- click(n=1)
# lineInSand
# occs[which(occs$longitude < lineInSand[1]), ]

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxGood, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus greggii
load('./species_records/00_fraxinus_greggii_bien_all_occurrences.rda')
# 154 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 149 occurrences

# lineInSand <- click(n=1)
# lineInSand
# occs[which(occs$longitude > lineInSand[1]), ]

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxGreg, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus latifolia
load('./species_records/00_fraxinus_latifolia_bien_all_occurrences.rda')
# 749 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 744 occurrences

lineInSand <- click(n=1)
lineInSand
occs[which(occs$latitude > lineInSand[2]), ]

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(bienRange_FraxLati, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus nigra
load('./species_records/00_fraxinus_nigra_bien_all_occurrences.rda')
# 48,624 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 48,564 occurrences

# southern Illinois occurrences are all FIA
# N. Dakota all FIA, except 1 verified on gbif
# W Iowa all FIA, except 1 verified on gbif #1462
# Virginia most from VegBank

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'Little range map')
plot(littleRange_FraxNigr, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus papillosa
load('./species_records/00_fraxinus_papillosa_bien_all_occurrences.rda')
# 30 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 30 occurrences

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(bienRange_FraxPapi, col = alpha('blue', 0.4), border = 'blue')
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# Fraxinus pennsylvanica
load('./species_records/00_fraxinus_pennsylvanica_bien_all_occurrences.rda')
# 93,206 occurrences

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 93,084 occurrences

# remove records before 1900
occs$Year <- substr(occs$date_collected, 1, 4)
occs <- occs[occs$Year >= 1900, ]
# 93,000 occurrences

# date_collected column is there twice, remove one of them
occs[18] <- NULL

# look at duplicates by date
duplicates <- occs[duplicated(occs$date_collected) | duplicated(occs$date_collected, fromLast = TRUE), ]
# this should remove rows where date & data source are the same
# we don't actually want to remove these records, 
# just trying to isolate where date is same, but data source different
duplicates <- duplicates[!duplicated(duplicates[, c("date_collected", "datasource")]), ] #4,549 obs
# now only keep duplicated dates to isolate these
duplicates <- duplicates[duplicated(duplicates$date_collected) | duplicated(duplicates$date_collected, 
                                                                            fromLast = TRUE), ]
# only want duplicate dates where lat/long are same in ones digit
duplicates$lat_simp <- substr(duplicates$latitude, 1, 2)
duplicates$long_simp <- substr(duplicates$longitude, 1, 3)
duplicates$date_lat_long <- paste0(duplicates$lat_simp, ", ", 
                                   duplicates$long_simp, ", ", 
                                   duplicates$ date_collected) # 976 obs
duplicates <- duplicates %>% group_by(date_lat_long) %>% filter(n() > 1) # 148 obs
# these appear to NOT be duplicates : 17 & 18, 131 & 132, 133 & 134, 135 & 136
false_duplicates <- c(17, 18, 131, 132, 133, 134, 135, 136)
# remove false duplicates
duplicates <- duplicates[-false_duplicates, ]

lineInSand <- click(n=1)
lineInSand
occs[which(occs$latitude > lineInSand[2]), ]

# visualize
par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map')
plot(littleRange_FraxPenn, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# determine which observations are outside range
extractOcc <- extract(bienRange_FraxPenn, occsSp)
# determine which indices are outside range
indexOutsideRange <- which(!complete.cases(extractOcc))
outsideRange <- occs[indexOutsideRange, ]

outsideRangeSp <- SpatialPointsDataFrame(outsideRange[, ll], data = outsideRange, 
                                         proj4 = getCRS('nad27', TRUE))
plot(outsideRangeSp, pch = 16, cex = 0.5, main = 'BIEN occurrences outside range')
plot(bienRange_FraxPenn, col = alpha('blue', 0.4), border = 'blue', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

# records removed: 
# record 2472, very far north, "taxon match fuzzy" GBIF


title(sub=date(), outer=TRUE, cex.sub=0.7, line=-2)

# Entire Fraxinus genus
load('./species_records/00_Fraxinus_bien_all_occurrences.rda')
# 247.681 occurrences
occsRaw <- Fraxinus

# remove records without coordinates or dates
occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
# 247,242 occurrences

# remove records before 1900
occs$Year <- substr(occs$date_collected, 1, 4)
occs <- occs[occs$Year >= 1900, ]
# 247,001 occurrences

par(mfrow=c(1,1))
occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4 = getCRS('nad27', TRUE))
plot(bienRange, col = alpha('blue', 0.4), border = 'blue')
plot(occsSp, pch = 16, cex = 0.5, col = 'red', main = 'BIEN range map', add = TRUE)
map("state", add = TRUE)
map("world", add = TRUE)

precision <- enmSdm::coordPrecision(occsSp, dms = FALSE) # 1st time
precision_dms <- enmSdm::coordPrecision(occsSp, dms = TRUE) # 2nd time

# save workspace
save.image('./02 - Modeling Workspace - Fraxinus Cleaning Data')
