
rm(list=ls())

require("ncdf4")
require("lattice")
require("ggplot2")

file <- "/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/env_data/Beyer/LateQuaternary_Environment.nc";

env_nc      <- ncdf4::nc_open(file)
longitude   <- ncdf4::ncvar_get(env_nc, "longitude")
latitude    <- ncdf4::ncvar_get(env_nc, "latitude")
years       <- ncdf4::ncvar_get(env_nc, "time")
months      <- ncdf4::ncvar_get(env_nc, "month")
temperature <- ncdf4::ncvar_get(env_nc, "temperature")
biome       <- ncdf4::ncvar_get(env_nc, "biome")
ncdf4::nc_close(env_nc)

my_year      <- -10000;   # 10,0000 BP
my_month     <- 6;        # June
my_longitude <- 0.1218;
my_latitude  <- 52.2053;  # Cambridge (UK)

p1 <- print(lattice::levelplot(biome[,,years == my_year], main = "Biome distribution, 10000 BP"))

p2 <- print(lattice::levelplot(temperature[,,months == my_month,years == my_year], main = "Mean June temperature, 10000 BP"))

lonID <- which.min(abs(longitude - my_longitude));
latID <- which.min(abs(latitude - my_latitude));
yearID <- which.min(abs(years - my_year));

p3 <- ggplot2::qplot(months, temperature[lonID,latID,,yearID], xlab = "Month", ylab = "Mean temperature",  geom=c("point", "line"))

mean_annual_temperature <- apply(temperature, c(1,2,4), mean)

p4 <- ggplot2::qplot(years, mean_annual_temperature[lonID,latID,], xlab = "Year", ylab = "Temperature", main = "Mean annual temperature time series, Cambridge (UK)", geom=c("point", "line"))

gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)

rm(list=ls())
library(ncdf4)
library(raster)

# File path on Lauren's computer:
file <- "/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/env_data/Beyer/LateQuaternary_Environment.nc";

ybp <- 21000
temperature <- brick(file, level = 1, varname = "temperature")

while (ybp >= 0) {
  rName <- paste0('temperature$X.', ybp)
  r <- raster(brick(rName))
  ybp <- ybp - 1000
}

temp <- brick(file, level = 1, varname = "temperature")
temp21 <- raster(temp$X.21000)
writeRaster(temp21, filename = "/Users/laurenjenkins/Documents/MOBOT/PVMvsENM/data_and_analyses/env_data/Beyer/tifs/temp.tif", 
            format = "GTiff", overwrite = TRUE)

# Load raster to examine if loaded properly (may be rotated 90 deg cw)
test <- raster('./data_and_analyses/env_data/Beyer/tifs/test_raster.tif')
plot(test)


timesteps <- c()
ybp <- 21000

while(ybp >= 0) {
  rName <- paste0(var, '$X.', ybp)
  timesteps <- append(timesteps, rName)
  ybp <- ybp - 1000
}
