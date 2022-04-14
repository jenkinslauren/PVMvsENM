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