## Started 18 June 2026
## By Ken

## Project phenology with variation in threshold parameters

# Set working directory to analyses folder

library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

library(geodata)
library(sf)
library(terra)

util <- new.env()
source("source/mcmc_analysis_tools_rstan.R", local=util)

# Model parameters from Morales-Castilla et al. (2020)
phen_params <- as.data.frame(array(NA, dim = c(12, 18)))
colnames(phen_params) <- c('variety', 'su_t1', 'su_topt', 'su_t2', 'su_c_min', 'c_bb',
                           'we_bb_tmin', 'we_bb_topt', 'we_bb_tmax', 'f_bb',
                           'we_flo_tmin', 'we_flo_topt', 'we_flo_tmax', 'f_flo',
                           'we_ver_tmin', 'we_ver_topt', 'we_ver_tmax', 'f_ver')

phen_params['variety'] <- c('Cabernet-Sauvignon', 'Chardonnay', 'Chasselas', 'Grenache',
                            'Merlot', 'Monastrell', 'Pinot noir', 'Riesling', 'Sauvignon',
                            'Syrah', 'Ugni blanc', 'Late-ripening composite')
phen_params['su_t1'] <- rep(-6.7, 12)
phen_params['su_topt'] <- rep(7.94, 12)
phen_params['su_t2'] <- rep(40.58, 12)
phen_params['su_c_min'] <- rep(-0.17, 12)
phen_params['c_bb'] <- c(166.02, 199.01, 196.68, 197.68, 152.14, 216.56, 198.53, 185.38, 203.48,
                         198.63, 199.74, 232.67)

phen_params['we_bb_tmin'] <- rep(0, 12)
phen_params['we_bb_topt'] <- rep(26.1, 12)
phen_params['we_bb_tmax'] <- rep(40, 12)
phen_params['f_bb'] <- c(20.24, 12.50, 12.40, 14.02, 19.91, 7.64, 12.98, 14.08, 13.06, 13.55,
                         17.80, 23.39)

phen_params['we_flo_tmin'] <- rep(0, 12)
phen_params['we_flo_topt'] <- rep(29.34, 12)
phen_params['we_flo_tmax'] <- rep(40, 12)
phen_params['f_flo'] <- c(23.61, 22.26, 24.46, 24.71, 24.96, 28.04, 21.85, 22.13, 25.34, 23.91,
                          24.41, 29.59)

phen_params['we_ver_tmin'] <- rep(0, 12)
phen_params['we_ver_topt'] <- rep(21.84, 12)
phen_params['we_ver_tmax'] <- rep(40, 12)
phen_params['f_ver'] <- c(65.30, 63.28, 53.40, 66.90, 64.45, 62.07, 59.50, 64.63, 59.27, 59.99,
                          67.36, 70.85)

# Calculate daily forcing according to Smoothed Wang & Engel model
swe_model <- function(temp, tmin, topt, tmax, delta){
  forcing <- rep(NA, length(temp))
  
  for(i in 1:length(temp)){
    if(temp[i] < tmin){
      forcing[i] <- 0.001
    } else if(temp[i] > tmax){
      forcing[i] <- 0.001
    } else if(topt > (tmin + tmax) / 2){
      phi <- (tmax - topt) / (topt - tmin)
      gamma <- (delta * tmax + topt - (delta + 1) * tmin) / (tmax - topt)
      
      a <- (temp[i] - tmin) / (topt - tmin)
      b <- (tmax - temp[i]) / (tmax - topt)
      forcing[i] <- (a * b ^ phi) ^ gamma + 0.001
    } else if(topt <= (tmin + tmax) / 2){
      phi <- (topt - tmin) / (tmax - topt)
      gamma <- ((delta + 1) * tmax - topt - delta * tmin) / (topt - tmin)
      
      a <- (tmax - temp[i]) / (tmax - topt)
      b <- (temp[i] - tmin) / (topt - tmin)
      forcing[i] <- (a * b ^ phi) ^ gamma + 0.001
    }
  }
  
  return(forcing)
}

# varieties and years to project
vars <- c('Aramon', 'Cabernet-Sauvignon', 'Carignan', 'Chardonnay', 'Chasselas',
          'Grenache', 'Merlot', 'Pinot noir', 'Sauvignon', 'Syrah', 'Ugni Blanc')
proj_var <- which(tolower(vars) %in% tolower(phen_params$variety))
phen_params <- subset(phen_params, tolower(variety) %in% tolower(vars))
proj_years <- c(1970:1979, 2015:2024)

temp <- rast('output/formatted_data/eur_temp/1950.nc')
coords <- crds(temp, df = TRUE, na.rm = TRUE) # check na.all = TRUE
phen <- read.csv('output/prepared_data/proj/1970.csv')

# Subset sample thermal sum threshold and temporal synchrony estimates, and parameter means for others
fit <- readRDS("fits/fit_main_wide_prior.rds")
samples <- util$extract_expectand_vals(fit)

tmax <- mean(c(samples[['T_max']]))
tmin <- mean(c(samples[['T_min']]))
delta <- mean(c(samples[['delta']]))
topt <- rep(NA, length(proj_var))
psi0 <- list()
sigma <- list()

for(var in 1:length(proj_var)){
  topt[var] <- mean(c(samples[[paste0('T_opt[', proj_var[var], ']')]]))
  psi0[[var]] <- c(samples[[paste0('Psi0[', proj_var[var], ']')]])
  sigma[[var]] <- c(samples[[paste0('sigma[', proj_var[var], ']')]])
}

# Subset Denmark and areas of Spain for projections
country <- world(resolution = 3, path = 'input/worldclim')
reg <- which(country$GID_0 == 'DNK')
reg <- mask(crop(temp, country[reg]), country[reg])
reg_cells <- crds(reg, df = TRUE, na.rm = TRUE)

coords_idx <- sapply(1:nrow(reg_cells), function(x){
  which(round(coords[, 'x'], 1) == round(reg_cells[x, 'x'], 1) & round(coords[, 'y'], 1) == round(reg_cells[x, 'y'], 1))
})

world_map <- world(resolution = 3, path = 'input/worldclim')
country <- 'Spain'
region <- which(tolower(world_map$NAME_0) == tolower(country))
region <- mask(crop(temp, world_map[region]), world_map[region])
wine_map <- read_sf(dsn = 'input/winegrowing_regions',
                    layer = 'WinegrowingRegionsoftheWorld')
wine_map <- st_transform(wine_map, crs = crs(world_map))

spainwinemap <-  crop(vect(wine_map), region)
spainwinemap$uniqid <- 1:dim(spainwinemap)[1]
castilla_extremadura <- subset(spainwinemap, spainwinemap$uniqid %in% c(7,8,27))
wg <- mask(region, castilla_extremadura)
region <- crds(wg, df = TRUE, na.rm = TRUE)

idx <- sapply(1:nrow(region), function(x){
  which(round(coords[, 'x'], 1) == round(region[x, 'x'], 1) & round(coords[, 'y'], 1) == round(region[x, 'y'], 1))
})

coords_idx <- c(coords_idx, idx)

# Project veraison using sample thermal sum threshold and temporal synchrony estimates otherwise parameter
phen_idx <- list()
for(i in 1:length(coords_idx)){
  phen_idx[[i]] <- which(round(phen[, 'x'], 1) == round(coords[coords_idx[i], 'x'], 1) & round(phen[, 'y'], 1) == round(coords[coords_idx[i], 'y'], 1))
}
phen_idx <- do.call(c, phen_idx)

for(year in proj_years){
  print(year)
  temp <- rast(paste0('output/formatted_data/eur_temp/', year, '.nc'))
  temp <- as.data.frame(temp)
  
  phen <- read.csv(paste0('output/projections/', year, '.csv'))
  phen <- phen[phen_idx,]
  
  for(site in 1:length(coords_idx)){
    if(site %% 100 == 0){print(site)}
    site_temp <- unname(unlist(temp[coords_idx[site], ])) - 273.15
    idx <- which(round(phen[, 'x'], 1) == round(coords[coords_idx[site], 'x'], 1) & round(phen[, 'y'], 1) == round(coords[coords_idx[site], 'y'], 1))
    
    for(var in 1:length(proj_var)){
      flo <- phen[idx[var], 'flo']
      ver <- rep(NA, length(psi0[[1]]))
      
      if(!is.na(flo) && flo <= length(site_temp)){
        forcing <- swe_model(site_temp[flo:length(site_temp)], tmin, topt[var], tmax, delta)
      } else {
        next
      }
      
      for(iter in 1:length(psi0[[1]])){
        thresh <- rlogis(1, psi0[[var]][iter], sigma[[var]][iter])
        ver[iter] <- which(cumsum(forcing) >= thresh)[1] + flo - 1
      }
      phen[idx[var], 8:4007] <- ver
    }
  }
  write.csv(phen, paste0('output/projections/', year, '_swe_ver_thresh_var.csv'), row.names = FALSE)
}
