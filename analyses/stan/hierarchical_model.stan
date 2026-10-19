// Started 30 July 2024
// By Ken, based on logistic_heterogeneous.stan

// Model fitting code for logistic model with grouping of choice

functions {
  // Forcing Function
  real forcing(real T, real T_min, real T_opt, real T_max, real delta) {
    if (T < T_min)
      return 0.001;
    if (T > T_max)
      return 0.001;
    
    if (T_opt > 0.5 * (T_min + T_max)) {
      real phi = (T_max - T_opt) / (T_opt - T_min);
      real gamma =   (delta * T_max + T_opt - (1 + delta) * T_min) 
      / (T_max - T_opt);
      real a = (T - T_min) / (T_opt - T_min);
      
      if (T > T_opt) {
        real b = (T_max - T) / (T_max - T_opt);
        real c = pow(b, phi);
        return (pow(a * c, gamma) + 0.001);
      } else {
        real b = (T_max - T_opt) / (T_max - T);
        real c = pow(b, -phi);
        return (pow(a * c, gamma) + 0.001);
      }
    } else {
      real phi = (T_opt - T_min) / (T_max - T_opt);
      real gamma =   ((1 + delta) * T_max - T_opt - delta * T_min) 
      / (T_opt - T_min);
      real b = (T_max - T) / (T_max - T_opt);
      
      if (T < T_opt) {
        real a = (T - T_min) / (T_opt - T_min);
        real c = pow(a, phi);
        return (pow(c * b, gamma) + 0.001);
      } else {
        real a = (T_opt - T_min) / (T - T_min);
        real c = pow(a, -phi);
        return (pow(c * b, gamma) + 0.001);
      }
    }
  }
  
  // Log Forcing Function
  real log_forcing(real T, real T_min, real T_opt, real T_max, real delta) {
    if (T < T_min)
      return log(0.001);
    if (T > T_max)
      return log(0.001);
    
    if (T_opt > 0.5 * (T_min + T_max)) {
      real gamma =   (delta * T_max + T_opt - (1 + delta) * T_min) 
      / (T_max - T_opt);
      real a = (T - T_min) / (T_opt - T_min);
      real b = (T_max - T_opt) / (T_opt - T_min);
      real c = (T_max - T) / (T_max - T_opt);
      return log(pow(a * pow(c, b), gamma) + 0.001);
    } else {
      real gamma =   ((1 + delta) * T_max - T_opt - delta * T_min) 
      / (T_opt - T_min);
      real a = (T - T_min) / (T_opt - T_min);
      real b = (T_opt - T_min) / (T_max - T_opt);
      real c = (T_max - T) / (T_max - T_opt);
      return log(pow(pow(a, b) * c, gamma) + 0.001);
    }
  }
}

data {
  int<lower=1> N_temp_grid;
  array[N_temp_grid] real temp_grid;
  
  //Temperature and divisions
  int<lower=1> N_temp_sets;
  int<lower=1> N_temps;
  
  array[N_temps] real temps;
  array[N_temp_sets] int<lower=1, upper=N_temps> temp_start_idxs;
  array[N_temp_sets] int<lower=1> N_days;
  
  //Phenology and divisions, and grouping
  int<lower=1> N;
  int<lower=1> N_groups;
  
  array[N] int<lower=1> flos;
  array[N] int<lower=1> vers;
  array[N] int<lower=1> groups;
  array[N_temp_sets] int<lower=1, upper=N> obs_start_idxs;
  array[N_temp_sets] int<lower=1> N_obs;
  
  // Prior parameters
  real mu_T_min;
  real sigma_T_min;
  
  real mu_T_max;
  real sigma_T_max;
  
  real mu_log_delta;
  real sigma_log_delta;
  
  // logit_x_opt
  real mu_nu_logit_x_opt;
  real sigma_nu_logit_x_opt;
  
  real mu_tau_logit_x_opt;
  real sigma_tau_logit_x_opt;
  
  // log_Psi0
  real mu_nu_log_Psi0;
  real sigma_nu_log_Psi0;
  
  real mu_tau_log_Psi0;
  real sigma_tau_log_Psi0;
  
  // log_sigma
  real mu_nu_log_sigma;
  real sigma_nu_log_sigma;
  
  real mu_tau_log_sigma;
  real sigma_tau_log_sigma;
}

parameters {
  // Forcing Function Parameters
  real T_min;                  // Minimum temperature for non-zero forcing (C)
  real<lower=T_min> T_max;     // Maximum temperature for non-zero forcing (C)
  real<lower=0> delta;         // Forcing function shape parameter (unitless)
  real nu_logit_x_opt;         // Logit relative temperature with maximal forcing mean (unitless)
  real <lower=0> tau_logit_x_opt;        // Logit relative temperature with maximal forcing sd (unitless)
  array[N_groups] real off;
  
  // Phenological Transition Parameters
  array[N_groups] real log_Psi0;     // Log phenology threshold (Log Forcing Units)
  real nu_log_Psi0;            // Log phenology threshold mean (Log Forcing Units)
  real<lower=0> tau_log_Psi0;  // Log phenology threshold sd (Log Forcing Units)
  
  array[N_groups] real log_sigma;    // Log phenology threshold scale (Log Forcing Units)
  real nu_log_sigma;           // Log phenology threshold scale mean (Log Forcing Units)
  real<lower=0> tau_log_sigma; // Log phenology threshold scale sd (Log Forcing Units)
}

transformed parameters {
  array[N_groups] real logit_x_opt;
  array[N_groups] real T_opt;
  array[N_groups] real Psi0;
  array[N_groups] real sigma;
  
  for(group in 1:N_groups){
    logit_x_opt[group] = nu_logit_x_opt + tau_logit_x_opt * off[group];
    T_opt[group] = (T_max - T_min) * inv_logit(logit_x_opt[group]) + T_min;
    Psi0[group] = exp(log_Psi0[group]);
    sigma[group] = exp(log_sigma[group]);
  }
}

model {
  T_min ~ normal(mu_T_min, sigma_T_min);
  T_max ~ normal(mu_T_max, sigma_T_max);
  delta ~ lognormal(mu_log_delta, sigma_log_delta);
  
  nu_logit_x_opt ~ normal(mu_nu_logit_x_opt, sigma_nu_logit_x_opt);
  tau_logit_x_opt ~ normal(mu_tau_logit_x_opt, sigma_tau_logit_x_opt);
  off ~ std_normal();
  
  log_Psi0 ~ normal(nu_log_Psi0, tau_log_Psi0);
  nu_log_Psi0 ~ normal(mu_nu_log_Psi0, sigma_nu_log_Psi0);
  tau_log_Psi0 ~ normal(mu_tau_log_Psi0, sigma_tau_log_Psi0);
  
  log_sigma ~ normal(nu_log_sigma, tau_log_sigma);
  nu_log_sigma ~ normal(mu_nu_log_sigma, sigma_nu_log_sigma);
  tau_log_sigma ~ normal(mu_tau_log_sigma, sigma_tau_log_sigma);
  
  for(temp_set in 1:N_temp_sets){
    // Subset local temperature data
    int temp_start_idx = temp_start_idxs[temp_set];
    int temp_end_idx = temp_start_idx + N_days[temp_set] - 1;
    
    array[N_days[temp_set]] real local_temps = temps[temp_start_idx:temp_end_idx];
    
    // Subset local phenological data
    int obs_start_idx = obs_start_idxs[temp_set];
    int obs_end_idx = obs_start_idx + N_obs[temp_set] - 1;
    
    array[N_obs[temp_set]] int local_flos = flos[obs_start_idx:obs_end_idx];
    array[N_obs[temp_set]] int local_vers = vers[obs_start_idx:obs_end_idx];
    array[N_obs[temp_set]] int local_groups = groups[obs_start_idx:obs_end_idx];
    
    // Compute daily forcings
    int earliest_forcing = min(local_flos);
    int latest_forcing = max(local_vers);
    array[N_groups, latest_forcing] real daily_forcings;
    
    for(grp in 1:N_groups){
      for(day in earliest_forcing:latest_forcing){
        daily_forcings[grp, day] = forcing(local_temps[day], T_min, T_opt[grp], T_max, delta);
      }
    }
    
    for(obs in 1:N_obs[temp_set]){
      real Psi = sum(daily_forcings[local_groups[obs], local_flos[obs]:local_vers[obs]]);
      real log_dPsidt = log_forcing(local_temps[local_vers[obs]], T_min, T_opt[local_groups[obs]], T_max, delta);
      target += logistic_lpdf(Psi | Psi0[local_groups[obs]], sigma[local_groups[obs]]) + log_dPsidt;
    }
  }
}

generated quantities {
  array[N_groups, N_temp_grid] real forcings;
  array[N] int<lower=1, upper=366> pred_events;
  
  array[N] real sum_forcings;
  array[N] int active_days;
  
  for(grp in 1:N_groups){
    for(tg in 1:N_temp_grid){
      forcings[grp, tg] = forcing(temp_grid[tg], T_min, T_opt[grp], T_max, delta);
    }
  }
  
  for(temp_set in 1:N_temp_sets){
    // Subset local temperature data
    int temp_start_idx = temp_start_idxs[temp_set];
    int temp_end_idx = temp_start_idx + N_days[temp_set] - 1;
    
    array[N_days[temp_set]] real local_temps = temps[temp_start_idx:temp_end_idx];
    
    // Subset local phenological data
    int obs_start_idx = obs_start_idxs[temp_set];
    int obs_end_idx = obs_start_idx + N_obs[temp_set] - 1;
    
    array[N_obs[temp_set]] int local_flos = flos[obs_start_idx:obs_end_idx];
    array[N_obs[temp_set]] int local_vers = vers[obs_start_idx:obs_end_idx];
    array[N_obs[temp_set]] int local_groups = groups[obs_start_idx:obs_end_idx];
    
    // Compute daily forcings
    int earliest_forcing = min(local_flos);
    array[N_groups, N_days[temp_set]] real daily_forcings;
    array[N_groups, N_days[temp_set]] real accum_forcings;
    array[N_groups, N_days[temp_set]] int forcing_days;
    array[N_groups, N_days[temp_set]] int accum_forcing_days;
    
    for(grp in 1:N_groups){
      for(day in 1:earliest_forcing){
        forcing_days[grp, day] = 0;
      }
      daily_forcings[grp] = rep_array(0, N_days[temp_set]);
      for(day in earliest_forcing:N_days[temp_set]){
        daily_forcings[grp, day] = forcing(local_temps[day], T_min, T_opt[grp], T_max, delta);
        if(daily_forcings[grp, day] > 0.001){
          forcing_days[grp, day] = 1;
        }
      }
      accum_forcings[grp] = cumulative_sum(daily_forcings[grp]);
      accum_forcing_days[grp] = cumulative_sum(forcing_days[grp]);
    }
    
    for(n in 1:N_obs[temp_set]){
      int obs_idx = obs_start_idxs[temp_set] + n - 1;
      int flo = local_flos[n];
      int ver = local_vers[n];
      int group = local_groups[n];
      real init_forcing = accum_forcings[group, flo] - daily_forcings[group, flo];
      
      sum_forcings[obs_idx] = accum_forcings[group, ver] - init_forcing;
      active_days[obs_idx] = accum_forcing_days[group, ver] - accum_forcing_days[group, flo] + forcing_days[group, flo];
      
      real l = logistic_rng(Psi0[group], sigma[group]);
      
      pred_events[obs_idx] = N_days[temp_set];
      for(day in flo:N_days[temp_set]){
        real Psi = accum_forcings[group, day] - init_forcing;
        if(Psi >= l){
          pred_events[obs_idx] = day;
          break;
        }
      }
    }
  }
}
