## Started 24 January 2025
## By Ken

## Format data and fit model

rm(list = ls())
options(stringsAsFactors = FALSE)

#stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

# Format temperature and phenological data
format_data <- function(temp_data, pheno_data, sites, vars){
  model_pheno_data <- subset(pheno_data, site %in% sites &
                               tolower(variety) %in% tolower(vars) &
                               year >= min(temp_data$year) &
                               year <= max(temp_data$year))
  model_pheno_data$site_var <- paste0(model_pheno_data$site, "_", model_pheno_data$variety)
  for(i in 1:length(vars)){
    model_pheno_data$variety[tolower(model_pheno_data$variety) == tolower(vars[i])] <- vars[i]
  }
  
  #assign IDs by site, variety, and site-variety
  site_ID <- data.frame(site = sites,
                        site_ID = 1:length(sites))
  var_ID <- data.frame(variety = vars,
                       var_ID = 1:length(vars))
  site_var_ID <- data.frame(site_var = unique(model_pheno_data$site_var),
                            site_var_ID = 1:length(unique(model_pheno_data$site_var)))
  
  model_pheno_data <- merge(model_pheno_data, site_ID, by = "site")
  model_pheno_data <- merge(model_pheno_data, var_ID, by = "variety")
  model_pheno_data <- merge(model_pheno_data, site_var_ID, by = "site_var")
  
  #get site years with corresponding sites for subsetting data
  site_years <- list()
  site <- list()
  for(i in 1:length(sites)){
    site_years[[i]] <- unique(model_pheno_data$year[model_pheno_data$site == sites[i]])
    site[[i]] <- rep(sites[i], length(site_years[[i]]))
  }
  site_years <- do.call("c", site_years)
  site <- do.call("c", site)
  N_site_years <- length(site_years)
  
  #arrange temperature data with corresponding indexing object
  N_days <- as.array(sapply(1:length(site_years), function(idx){
    length(temp_data$doy[temp_data$site == site[idx] & temp_data$year == site_years[idx]])
  })
  )
  temp_start_idxs <- as.array(rep(NA, N_site_years))
  
  N_obs_temp <- sum(N_days)
  obs_temp <- rep(NA, N_obs_temp)
  
  start_idx <- 1
  end_idx <- start_idx + N_days[1] - 1
  
  temp_start_idxs[1] <- start_idx
  obs_temp[start_idx:end_idx] <- temp_data$t_mean[temp_data$site == site[1] & temp_data$year == site_years[1]]
  
  for(idx in 2:N_site_years){
    if(N_site_years < 2) {next}
    start_idx <- end_idx + 1
    end_idx <- start_idx + N_days[idx] - 1
    
    temp_start_idxs[idx] <- start_idx
    obs_temp[start_idx:end_idx] <- temp_data$t_mean[temp_data$site == site[idx] & temp_data$year == site_years[idx]]
  }
  
  #arrange observations by site-year (temperature set) with corresponding indexing object
  N_obs <- as.array(rep(NA, N_site_years))
  obs_start_idxs <- as.array(rep(NA, N_site_years))
  
  N <- nrow(model_pheno_data)
  flos <- rep(NA, N)
  vers <- rep(NA, N)
  site_var_IDs <- rep(NA, N)
  site_IDs <- rep(NA, N)
  var_IDs <- rep(NA, N)
  
  N_obs[1] <- nrow(model_pheno_data[model_pheno_data$site == site[1] & model_pheno_data$year == site_years[1],])
  start_idx <- 1
  end_idx <- start_idx + N_obs[1] - 1
  
  obs_start_idxs[1] <- start_idx
  flos[start_idx:end_idx] <- model_pheno_data$flo[model_pheno_data$site == site[1] &
                                                    model_pheno_data$year == site_years[1]]
  vers[start_idx:end_idx] <- model_pheno_data$ver[model_pheno_data$site == site[1] &
                                                    model_pheno_data$year == site_years[1]]
  site_var_IDs[start_idx:end_idx] <- model_pheno_data$site_var_ID[model_pheno_data$site == site[1] &
                                                                    model_pheno_data$year == site_years[1]]
  site_IDs[start_idx:end_idx] <- model_pheno_data$site_ID[model_pheno_data$site == site[1] &
                                                            model_pheno_data$year == site_years[1]]
  var_IDs[start_idx:end_idx] <- model_pheno_data$var_ID[model_pheno_data$site == site[1] &
                                                          model_pheno_data$year == site_years[1]]
  
  for(idx in 2:N_site_years){
    if(N_site_years < 2) {next}
    N_obs[idx] <- nrow(model_pheno_data[model_pheno_data$site == site[idx] &
                                          model_pheno_data$year == site_years[idx],])
    start_idx <- end_idx + 1
    end_idx <- start_idx + N_obs[idx] - 1
    
    obs_start_idxs[idx] <- start_idx
    flos[start_idx:end_idx] <- model_pheno_data$flo[model_pheno_data$site == site[idx] &
                                                      model_pheno_data$year == site_years[idx]]
    vers[start_idx:end_idx] <- model_pheno_data$ver[model_pheno_data$site == site[idx] &
                                                      model_pheno_data$year == site_years[idx]]
    site_var_IDs[start_idx:end_idx] <- model_pheno_data$site_var_ID[model_pheno_data$site == site[idx] &
                                                                      model_pheno_data$year == site_years[idx]]
    site_IDs[start_idx:end_idx] <- model_pheno_data$site_ID[model_pheno_data$site == site[idx] &
                                                              model_pheno_data$year == site_years[idx]]
    var_IDs[start_idx:end_idx] <- model_pheno_data$var_ID[model_pheno_data$site == site[idx] &
                                                            model_pheno_data$year == site_years[idx]]
  }
  
  list("N_temp_sets" = N_site_years,
       "N_temps" = N_obs_temp,
       "temps" = obs_temp,
       "temp_start_idxs" = temp_start_idxs,
       "N_days" = N_days,
       "N" = N,
       "flos" = flos,
       "vers" = vers,
       "site_var_IDs" = site_var_IDs,
       "site_IDs" = site_IDs,
       "var_IDs" = var_IDs,
       "obs_start_idxs" = obs_start_idxs,
       "N_obs" = N_obs,
       "site_var_ID" = site_var_ID,
       "site_ID" = site_ID,
       "var_ID" = var_ID
  )
}

# Create list of prior parameters for specific set
make_prior_set <- function(prior_sets, set_num){
  prior_params <- list(
    "mu_T_min" = prior_sets[['mtmin']][set_num],
    "sigma_T_min" = prior_sets[['stmin']][set_num],
    
    "mu_T_max" = prior_sets[['mtmax']][set_num],
    "sigma_T_max" = prior_sets[['stmax']][set_num],
    
    "mu_log_delta" = prior_sets[['mld']][set_num],
    "sigma_log_delta" = prior_sets[['sld']][set_num],
    
    # logistic hierarchical model-specific parameters
    "mu_nu_logit_x_opt" = prior_sets[['mnlx']][set_num],
    "sigma_nu_logit_x_opt" = prior_sets[['snlx']][set_num],
    "mu_tau_logit_x_opt" = prior_sets[['mtlx']][set_num],
    "sigma_tau_logit_x_opt" = prior_sets[['stlx']][set_num],
    
    "mu_nu_log_Psi0" = prior_sets[['mnlp']][set_num],
    "sigma_nu_log_Psi0" = prior_sets[['snlp']][set_num],
    "mu_tau_log_Psi0" = prior_sets[['mtlp']][set_num],
    "sigma_tau_log_Psi0" = prior_sets[['stlp']][set_num],
    
    "mu_nu_log_sigma" = prior_sets[['mnls']][set_num],
    "sigma_nu_log_sigma" = prior_sets[['snls']][set_num],
    "mu_tau_log_sigma" = prior_sets[['mtls']][set_num],
    "sigma_tau_log_sigma" = prior_sets[['stls']][set_num]
  )
  
  return(prior_params)
}

# Generate sample parameter values from prior parameters for logistic_hierarchical model
inits_lh <- function(N_groups, prior_params){
  T_min <- rnorm(1, prior_params[['mu_T_min']], prior_params[['sigma_T_min']])
  T_max <- rnorm(1, prior_params[['mu_T_max']], prior_params[['sigma_T_max']])
  delta <- rlnorm(1, prior_params[['mu_log_delta']], prior_params[['sigma_log_delta']])
  
  nu_logit_x_opt <- rnorm(1, prior_params[['mu_nu_logit_x_opt']], prior_params[['sigma_nu_logit_x_opt']])
  tau_logit_x_opt <- abs(rnorm(1, prior_params[['mu_tau_logit_x_opt']], prior_params[['sigma_tau_logit_x_opt']]))
  logit_x_opt <- as.array(rnorm(N_groups, nu_logit_x_opt, tau_logit_x_opt))
  
  nu_log_Psi0 <- rnorm(1, prior_params[['mu_nu_log_Psi0']], prior_params[['sigma_nu_log_Psi0']])
  tau_log_Psi0 <- abs(rnorm(1, prior_params[['mu_tau_log_Psi0']], prior_params[['sigma_tau_log_Psi0']]))
  log_Psi0 <- as.array(rnorm(N_groups, nu_log_Psi0, tau_log_Psi0))
  
  nu_log_sigma <- rnorm(1, prior_params[['mu_nu_log_sigma']], prior_params[['sigma_nu_log_sigma']])
  tau_log_sigma <- abs(rnorm(1, prior_params[['mu_tau_log_sigma']], prior_params[['sigma_tau_log_sigma']]))
  log_sigma <- as.array(abs(rnorm(N_groups, nu_log_sigma, tau_log_sigma)))
  
  return(list("T_min" = T_min, "T_max" = T_max, "delta" = delta,
              "logit_x_opt" = logit_x_opt,
              "nu_logit_x_opt" = nu_logit_x_opt, "tau_logit_x_opt" = tau_logit_x_opt,
              "log_Psi0" = log_Psi0,
              "nu_log_Psi0" = nu_log_Psi0, "tau_log_Psi0" = tau_log_Psi0,
              "log_sigma" = log_sigma,
              "nu_log_sigma" = nu_log_sigma, "tau_log_sigma" = tau_log_sigma))
}

# List of prior parameters
priors <- list(
  'mtmin' = rep(0, 2),
  'stmin' = c(7.5, 2.5) / 2.32,
  'mtmax' = rep(40, 2),
  'stmax' = c(7.5, 2.5) / 2.32,
  'mld' = rep(0.8, 2),
  'sld' = rep(0.8, 2),
  
  'mnlp' = rep(3.98, 2),
  'snlp' = rep(0.23, 2) / 2.32,
  'mtlp' = rep(0.49, 2) / 2.32,
  'stlp' = rep(0.23, 2) / 2.32^2,
  
  'mnls' = rep(1.16, 2),
  'snls' = rep(0.08, 2) / 2.32,
  'mtls' = rep(1.08, 2) / 2.32,
  'stls' = rep(0.08, 2) / 2.32^2,
  
  'mnlx' = c(1.11, 0.79),
  'snlx' = c(1.41, 0.59) / 2.32,
  'mtlx' = c(1.75, 0.94) / 2.32,
  'stlx' = c(1.41, 0.59) / 2.32^2
)

prior_set_name <- c("main_wide", "alternative_narrow")

# source('source/model_fxs.R')
# source('source/proj_fxs.R')

#sites and varieties to work on
sites <- c("Bellay", "Bergheim", "Lafite", "Vassal")
vars <- c('Aramon', 'Cabernet-Sauvignon', 'Carignan', 'Chardonnay', 'Chasselas',
          'Grenache', 'Merlot', 'Pinot noir', 'Sauvignon', 'Syrah', 'Ugni Blanc')

#retrieve and format data
temp_data <- read.csv("output/formatted_data/simu_temp_data.csv", header = TRUE)
pheno_data <- read.csv("output/formatted_data/simu_pheno_data.csv", header = TRUE)
pheno_data <- subset(pheno_data, year >= 1981)

data <- format_data(temp_data, pheno_data, sites, vars)

temp_grid <- seq(-10, 50, 0.25)
temp_grid <- list("N_temp_grid" = length(temp_grid), "temp_grid" = temp_grid)

data <- append(data, temp_grid, after = 0)
names(data)[13] <- "groups"
grouping <- list("N_groups" = nrow(data$var_ID))
data <- append(data, grouping)
data_backup <- data

set.seed(49384838)

for(i in 1:2){
  prior_params <- make_prior_set(priors, i)
  data <- data_backup
  data <- append(data, prior_params)
  
  inits <- function(chain_id){
    params <- inits_lh(data$N_groups, prior_params)
    return(params)
  }
  
  # Quantify posterior distribution
  fit <- stan(file = "stan/hierarchical_model.stan", 
              data = data, init = inits, seed = 8438338,
              warmup = 1000, iter = 2000, refresh = 500, chains = 4)
  
  saveRDS(fit, paste0("fits/fit_simu_", prior_set_name[i], "_prior.rds"))
}
