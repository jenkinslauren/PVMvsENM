# Load Little Range Maps
# Lauren Jenkins, Fall 2021

### setup ###
rm(list = ls())
setwd('/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/')

# load packages
library(sp)
library(raster)
library(enmSdm)

library(BIEN)
library(dismo)
library(brglm2)
library(cluster)
library(maxnet)
library(rgbif)
library(rgeos)
library(geosphere)
library(rgdal)
library(scales)
library(omnibus)
library(RColorBrewer)
library(rJava)
library(maps)

library(statisfactory)
library(legendary)

ll <- c('longitude', 'latitude')

# Little range maps
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

# To avoid thinking it's an attack ... Sys.sleep(1)
# utils::download.file(thisUrl, destfile=filePath, mode='wb', quiet=TRUE)
# could also try in dismo package shapefile function just enter the url in there
# or try going to the original repo rather than fork

# sp::spTransform()
# rbind or c to merge raster polygons
# rgeos::gUnion()

