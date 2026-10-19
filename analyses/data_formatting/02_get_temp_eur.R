## Started 5 July 2025
## By Ken

## Get Europe climate data

rm(list = ls())
options(stringsAsFactors = FALSE)

library(sf)
library(terra)

startyr <- 1950
endyr <- 2024
temp <- rast('/shared/data/climate/era5land/tmp_1950.nc')

world_map <- read_sf(dsn = "data/world_continents",
                     layer = "world_continents")
eur_map <- world_map[which(st_drop_geometry(world_map[, "CONTINENT"]) == "Europe"), ]
eur_map <- st_transform(eur_map, crs = crs(temp))

for(year in startyr:endyr){
  temp <- rast(paste0('/shared/data/climate/era5land/tmp_', year ,'.nc'))
  east <- crop(temp, ext(c(-0.05, 180, -91, 90)))
  west <- crop(temp, ext(c(180, 360, -91, 90)))
  west <- shift(west, dx = -360)
  
  temp <- merge(west, east)
  eur_temp <- crop(temp, vect(eur_map))
  eur_temp <- mask(eur_temp, vect(eur_map$geometry))
  writeCDF(eur_temp, paste0('output/formatted_data/eur_temp/', year, '.nc'))
}
