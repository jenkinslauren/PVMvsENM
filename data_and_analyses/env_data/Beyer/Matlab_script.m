clear; clc; close all;

file = 'LateQuaternary_Environment.nc';

%% Read variables from NetCDF file

longitude   = ncread(file, 'longitude');
latitude    = ncread(file, 'latitude');
years       = ncread(file, 'time');
months      = ncread(file, 'month');
temperature = ncread(file, 'temperature');
biome       = ncread(file, 'biome');

%% Specify relevant point in time and space

my_year      = -10000;    % 10,000 BP
my_month     = 6;         % June
my_longitude = 0.1218;
my_latitude  = 52.2053;   % Cambridge, UK

figure(1);

%% Global biome map
subplot(2,2,1);
imagesc(longitude, latitude, biome(:,:, years == my_year)'); axis xy;
title(['Biome distribution, 10000 BP']);

%% Global temperature map
subplot(2,2,2);
imagesc(longitude, latitude, temperature(:,:, months == my_month, years == my_year)'); axis xy;
colorbar;
title(['Mean June temperature, 10000 BP']);

% use nearest grid points (alternatively, use 3D interpolation)
[~,lonID]   = min(abs(longitude - my_longitude));
[~,latID]   = min(abs(latitude - my_latitude));
[~,yearID]  = min(abs(years - my_year));

%% Monthly temperature distribution
subplot(2,2,3);
plot(months, squeeze(temperature(lonID, latID, :, yearID)), '-o');
title(['Monthly temperature distribution, Cambridge (UK), 10000 BP']);
xticks(1:12); xlabel("Month");
ylabel('Mean temperature');

%% Mean annual temperature time series
mean_annual_temperature = mean(temperature,3);

subplot(2,2,4);
plot(years, squeeze(mean_annual_temperature(lonID, latID, :)), '-o');
title(['Mean annual temperature time series, Cambridge (UK)']);
xlabel('Year');

