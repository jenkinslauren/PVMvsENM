import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt

climate_file = 'LateQuaternary_Environment.nc';

## Read variables from NetCDF file

nc = nc4.Dataset(climate_file,'r')
longitude   = nc.variables['longitude'][:]
latitude    = nc.variables['latitude'][:]
years       = nc.variables['time'][:]
months      = nc.variables['month'][:]
temperature = nc.variables['temperature'][:]
biome       = nc.variables['biome'][:]

## Specify relevant point in time and space

my_year      = -10000    # 10,000 BP
my_month     = 5         # June
my_longitude = 0.1218
my_latitude  = 52.2053   # Cambridge (UK)

LON,LAT = np.meshgrid(longitude,latitude)

# use nearest grid points (alternatively, use 3D interpolation)
lonID   = np.argmin(np.abs(longitude - my_longitude));
latID   = np.argmin(np.abs(latitude - my_latitude));
yearID  = np.argmin(np.abs(years - my_year));
monthID = np.argmin(np.abs(months - my_month));

## Global biome map
ax = plt.subplot(2,2,1)
plt.imshow(biome[years == my_year,...].squeeze());
plt.gca().invert_yaxis()
ax.set_title('Biome distribution, 10000 BP')

## Global temperature map
ax = plt.subplot(2,2,2)
plt.imshow(temperature[yearID,monthID,...].squeeze());
plt.colorbar();
plt.gca().invert_yaxis()
ax.set_title('Mean June temperature, 10000 BP')

#Monthly temperature distribution
ax = plt.subplot(2,2,3)
ax.plot(months, temperature[yearID,:,latID,lonID], ls='-', marker='.')
ax.set_title('Monthly temperature distribution, Cambridge (UK), 10000 BP')
ax.set_xticks(range(1,13))
ax.set_xlabel("Month")
ax.set_ylabel("Mean temperature")

# Mean annual temperature time series
mean_annual_temperature = np.mean(temperature,axis=1);
ax = plt.subplot(2,2,4)
ax.plot(years, mean_annual_temperature[:, latID, lonID], ls='-', marker='.')
ax.set_title('Mean annual temperature time series, Cambridge (UK)')
ax.set_xlabel('Year');

plt.show();

