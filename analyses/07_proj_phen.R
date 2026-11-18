## Started 25 June 2025
## By Ken

## Project phenology

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working directory to analyses folder

#stan
library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

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

# Calculate daily chilling according to Smoothed Utah model
su_model <- function(temp, t1, topt, t2, c_min){
  chilling <- rep(NA, length(temp))
  
  for(i in 1:length(temp)){
    if(temp[i] < t1){
      chilling[i] <- 1 / (1 + exp(-4 * (temp[i] - t1) / (topt - t1)))
    } else if(temp[i] < topt){
      chilling[i] <- 1 - (0.5 / ((t1 - topt) ^ 2)) * (temp[i] - topt) ^ 2
    } else if(temp[i] < t2){
      chilling[i] <- 1 - (1 - c_min) * ((temp[i] - topt) ^ 2) / (2 * ((t2 - topt) ^ 2))
    } else {
      chilling[i] <- c_min + (1 - c_min) / (1 + exp(-4 * (t2 - temp[i]) / (t2 - topt)))
    }
  }
  
  return(chilling)
}

# Calculate daily forcing according to Wang & Engel model
we_model <- function(temp, tmin, topt, tmax){
  a <- log(2) / log((tmax - tmin) / (topt - tmin))
  forcing <- rep(NA, length(temp))
  
  for(i in 1:length(temp)){
    if(temp[i] > tmin & temp[i] < tmax){
      forcing[i] <- (2 * (temp[i] - tmin) ^ a * (topt - tmin) ^ a - (temp[i] - tmin) ^ (2 * a)) / ((topt - tmin) ^ (2 * a))
    } else {
      forcing[i] <- 0
    }
  }
  
  return(forcing)
}

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

# Project phenology according to model parameters from Morales-Castilla et al. (2020)
prev_temp <- rast(paste0('output/formatted_data/eur_temp/', proj_years[1] - 1 , '.nc'))
prev_temp <- as.data.frame(prev_temp)
if((proj_years[1] - 1) %% 4 == 0){
  prev_leap <- TRUE
} else {
  prev_leap <- FALSE
}

for(year in proj_years){
  proj <- list()
  counter <- 0
  
  cur_temp <- rast(paste0('output/formatted_data/eur_temp/', year, '.nc'))
  cur_temp <- as.data.frame(cur_temp)
  if(year %% 4 == 0){
    cur_leap <- TRUE
    days <- 366
  } else {
    cur_leap <- FALSE
    days <- 365
  }
  
  for(site in 1:nrow(coords)){
    temp <- rep(NA, 153 + days)
    if(prev_leap){
      temp[1:153] <- unlist(prev_temp[site, 214:366])
    } else {
      temp[1:153] <- unlist(prev_temp[site, 213:365])
    }
    temp[154:length(temp)] <- unlist(cur_temp[site, 1:days])
    temp <- temp - 273.15
    
    chilling <- su_model(temp = temp, t1 = -6.7, topt = 7.94, t2 = 40.58, c_min = -0.17)
    bb_forcing <- we_model(temp = temp, tmin = 0, topt = 26.1, tmax = 40)
    flo_forcing <- we_model(temp = temp, tmin = 0, topt = 29.34, tmax = 40)
    ver_forcing <- we_model(temp = temp, tmin = 0, topt = 21.84, tmax = 40)
    
    for(var in 1:length(proj_var)){
      db <- NA
      bb <- NA
      flo <- NA
      ver <- NA
      
      db <- round(approx(x = cumsum(chilling), y = 1:length(chilling), xout = phen_params$c_bb[var])$y)
      if(!is.na(db)){
        db:length(temp)
        bb <- which(cumsum(bb_forcing[db:length(temp)]) >= phen_params$f_bb[var])[1] + db
      }
      if(!is.na(bb)){
        flo <- which(cumsum(flo_forcing[bb:length(temp)]) >= phen_params$f_flo[var])[1] + bb
      }
      if(!is.na(flo)){
        ver <- which(cumsum(ver_forcing[flo:length(temp)]) >= phen_params$f_ver[var])[1] + flo
      }
      
      counter <- counter + 1
      proj[[counter]] <- data.frame('x' = coords$x[site],
                                    'y' = coords$y[site],
                                    'variety' = proj_var[var],
                                    'db' = db - 153,
                                    'bb' = bb - 153,
                                    'flo' = flo - 153,
                                    'ver' = ver - 153)
    }
  }
  proj <- do.call(rbind, proj)
  write.csv(proj, paste0('output/projections/', year, '.csv'), row.names = FALSE)
  
  prev_temp <- cur_temp
  prev_leap <- cur_leap
}

# Project veraison according to model parameter means
fit <- readRDS("fits/fit_main_wide_prior.rds")
samples <- util$extract_expectand_vals(fit)

tmax <- mean(c(samples[['T_max']]))
tmin <- mean(c(samples[['T_min']]))
delta <- mean(c(samples[['delta']]))
topt <- rep(NA, length(proj_var))
psi0 <- rep(NA, length(proj_var))
sigma <- rep(NA, length(proj_var))

for(var in 1:length(proj_var)){
  topt[var] <- mean(c(samples[[paste0('T_opt[', proj_var[var], ']')]]))
  psi0[var] <- mean(c(samples[[paste0('Psi0[', proj_var[var], ']')]]))
  sigma[var] <- mean(c(samples[[paste0('sigma[', proj_var[var], ']')]]))
}

for(year in proj_years){
  temp <- rast(paste0('output/formatted_data/eur_temp/', year, '.nc'))
  temp <- as.data.frame(temp)
  
  phen <- read.csv(paste0('output/projections/', year, '.csv'))
  
  for(site in 1:nrow(coords)){
    site_temp <- unname(unlist(temp[site, ])) - 273.15
    idx <- which(round(phen[, 'x'], 1) == round(coords[site, 'x'], 1) & round(phen[, 'y'], 1) == round(coords[site, 'y'], 1))
    
    for(var in 1:length(proj_var)){
      flo <- phen[idx[var], 'flo']
      ver <- NA
      
      if(!is.na(flo) && flo <= length(site_temp)){
        forcing <- swe_model(site_temp[flo:length(site_temp)], tmin, topt[var], tmax, delta)
        thresh <- rlogis(1, psi0[var], sigma[var])
        
        ver <- which(cumsum(forcing) >= thresh)[1] + flo - 1
      }
      
      phen[idx[var], 'swe_ver'] <- ver
    }
  }
  write.csv(phen, paste0('output/projections/', year, '_swe_ver.csv'), row.names = FALSE)
}
