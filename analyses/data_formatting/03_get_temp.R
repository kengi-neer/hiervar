## Started 1 July 2025
## By Ken

## Extract ERA5-Land data for data points

rm(list = ls())
options(stringsAsFactors = FALSE)

library(sf)
library(terra)

startyr <- 1950
endyr <- 2024
pheno_data <- read.csv("output/formatted_data/pheno_data_latlong.csv", header = TRUE)
temp <- rast('output/formatted_data/eur_temp/1950.nc')

coords <- unique(pheno_data[, c('site', 'location', 'longitude', 'latitude')], na.rm = TRUE)
coords <- aggregate(coords[c('longitude', 'latitude')],
                    by = coords['location'],
                    FUN = function(x){return(mean(x, na.rm = TRUE))})
site_coords <- st_as_sf(coords, coords = c('longitude', 'latitude'), crs = 4326)
site_coords <- st_transform(site_coords, crs = crs(temp))

temp_coords <- extract(temp[[1]], vect(site_coords), xy = TRUE, search_radius = 20000)
coords[, c('longitude', 'latitude')] <- temp_coords[, c('x', 'y')]

coords <- st_as_sf(coords, coords = c('longitude', 'latitude'), crs = 4326)
coords <- st_transform(coords, crs = crs(temp))

temps <- list()
for(year in startyr:endyr){
  temp <- rast(paste0('output/formatted_data/eur_temp/', year, '.nc'))
  temp <- extract(temp, vect(coords$geometry), ID = FALSE)
  
  site_temps <- list()
  for(site in 1:nrow(temp)){
    site_temps[[site]] <- data.frame('site' = rep(st_drop_geometry(coords)[site, 'location'], ncol(temp)),
                                     'year' = rep(year, ncol(temp)),
                                     'doy' = 1:ncol(temp),
                                     't_mean' = unlist(temp[site, 1:ncol(temp)]))
  }
  temps[[year - startyr + 1]] <- do.call(rbind, site_temps)
}
temps <- do.call(rbind, temps)
temps[, 't_mean'] <- temps[, 't_mean'] - 273.15

write.csv(temps, 'output/formatted_data/temp_data.csv', row.names = FALSE)
