## Started 5 May 2025
## By Ken

## Simulate data

rm(list = ls())
options(stringsAsFactors = FALSE)

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


#sites and varieties to work on from available PNAS 2020 paper
sites <- c("Bellay", "Bergheim", "Lafite", "Vassal")
vars <- c('Aramon', 'Cabernet-Sauvignon', 'Carignan', 'Chardonnay', 'Chasselas',
          'Grenache', 'Merlot', 'Pinot noir', 'Sauvignon', 'Syrah', 'Ugni Blanc')
vars <- vars[which(tolower(vars) %in% tolower(phen_params$variety))]
phen_params <- subset(phen_params, tolower(variety) %in% tolower(vars))

#retrieve and format data
temp_data <- read.csv("output/formatted_data/temp_data.csv", header = TRUE)
pheno_data <- read.csv("output/formatted_data/pheno_data.csv", header = TRUE)
pheno_data <- subset(pheno_data, tolower(variety) %in% tolower(vars) &
                       year >= 1981)

# Fixed parameters for simulation
simu_params <- phen_params[, c('variety', 'we_ver_tmin', 'we_ver_topt', 'we_ver_tmax', 'f_ver')]
delta <- rep(exp(0.8), nrow(simu_params))
sigma <- c(2.71, 2.43, 3.05, 4.46, 2.90, 3.56, 3.01, 0.88, 5.17)
lx <- rep(NA, nrow(simu_params))
nu_lx <- rep(NA, nrow(simu_params))
tau_lx <- rep(NA, nrow(simu_params))
nu_psi0 <- rep(NA, nrow(simu_params))
tau_psi0 <- rep(NA, nrow(simu_params))
nu_sigma <- rep(1.05, nrow(simu_params))
tau_sigma <- rep(0.5, nrow(simu_params))
add <- list(delta, sigma, lx, nu_lx, tau_lx, nu_psi0, tau_psi0, nu_sigma, tau_sigma)
add <- do.call(cbind, add)
simu_params <- cbind(simu_params, add)
colnames(simu_params) <- c('variety', 'tmin', 'topt', 'tmax', 'psi0', 'delta', 'sigma',
                           'lx', 'nu_lx', 'tau_lx', 'nu_psi0', 'tau_psi0',
                           'nu_sigma', 'tau_sigma')
simu_params$topt <- c(24.43, 22.32, 20.76, 20.08, 19.45, 21.19, 23.86, 21.22, 23.25)
simu_params$lx <- logit(simu_params$topt / 40)
simu_params$nu_lx <- mean(simu_params$lx)
simu_params$tau_lx <- sd(simu_params$lx)
simu_params$nu_psi0 <- mean(log(simu_params$psi0))
simu_params$tau_psi0 <- sd(log(simu_params$psi0))
write.csv(simu_params, 'output/formatted_data/simu_params.csv', row.names = FALSE)

#simulate variety-level xopt values
#match <- FALSE
#while(!match){
#  sample <- round(rnorm(9, 21.84, 4 / 2.32), 2)
#  mean <- mean(sample)
#  sd <- round(sd(sample), 3)
#  
#  if(mean == 21.84 & sd == 1.724){
#    match <- TRUE
#  }
#}

#simulate variety-level sigma values
# match <- FALSE
# while(!match){
#   sample <- round(rlnorm(9, 1.05, 0.35), 2)
#   mean <- round(mean(log(sample)), 2)
#   sd <- round(sd(log(sample)), 2)
#   
#   if(mean == 1.05 & sd == 0.50){
#     match = TRUE
#   }
# }
#

#check day of flowering distributions
#pdf('flo.pdf')
#flo <- sapply(sort(unique(pheno_data$variety)), function(x){
#  var_data <- subset(pheno_data, variety == x)
#  hist(var_data$flo, xlab = '', main = paste0('Flowering (', x, ')'))
#  return(nrow(var_data))
#})
#dev.off()

# Simulate all data
startyr <- 1981
endyr <- 2024

temp <- list()
counter <- 0
for(i in 1:1){
  for(j in 1){
    if(j %% 4 == 0){
      days <- 366
    } else {
      days <- 365
    }
    
    for(k in 1:days){
      counter <- counter + 1
      temp[[counter]] <- data.frame('site' = sites[i],
                                    'year' = j,
                                    'doy' = k,
                                    't_mean' = runif(1, 0, 40))
    }
  }
}
temp <- do.call(rbind, temp)
write.csv(temp, 'output/formatted_data/simu_temp_data.csv', row.names = FALSE)

# Simulate flowering
mean <- c(145, 150, 150, 150, 160, 150, 165, 145, 150, 150, 160)
sd <- c(15, 15, 15, 20, 30, 20, 35, 15, 15, 15, 20) / 2.32

phen <- list()
counter <- 0
for(i in 1:1){
  for(j in 1:1){
    for(k in 1){
      counter <- counter + 1
      phen[[counter]] <- data.frame('variety' = vars[j],
                                    'clone' = '',
                                    'site' = sites[i],
                                    'year' = k,
                                    'flo' = round(rnorm(1, mean[i], sd[i])),
                                    'ver' = 365)
    }
  }
}
phen <- do.call(rbind, phen)
phen$flo <- phen$flo - min(phen$flo) + 1

# Simulate veraison
for(i in 1:1){
  for(j in 1:1){
    for(k in 1){
      calc <- subset(temp, site == sites[i] & year == k)
      idx <- which(tolower(phen$variety) == tolower(vars[j]) & phen$site == sites[i] & phen$year == k)
      flo <- phen[idx, 'flo']
      
      calc <- unlist(calc[flo:nrow(calc), 't_mean'])
      
      forcing <- swe_model(calc, simu_params$tmin[j], simu_params$topt[j], simu_params$tmax[j], simu_params$delta[j])
      thresh <- rlogis(1, simu_params$psi0[j], simu_params$sigma[j])
      
      phen[which(tolower(phen$variety) == tolower(vars[j]) &
                   phen$site == sites[i] &
                   phen$year == k), 'ver'] <- which(cumsum(forcing) >= thresh)[1] + flo - 1
    }
  }
}
phen <- subset(phen, !is.na(ver))
write.csv(phen, 'output/formatted_data/simu_pheno_data.csv', row.names = FALSE)
