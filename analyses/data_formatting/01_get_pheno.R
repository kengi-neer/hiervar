## Started 19 July 2024
## By Ken, some work by Lizzie, based on grabdataTEMPO.R

## Format TEMPO phenological data and merge with processed Vassal data from Vassal server

rm(list = ls())
options(stringsAsFactors = FALSE)

# Setting working directory... 
if(length(grep("lizzie", getwd()) > 0)) { 
  setwd("~/Documents/git/projects/vinmisc/survin/analyses")
  par(family="serif")
} else if(length(grep("Ken", getwd()) > 0)){
  setwd("/Users/Ken Michiko Samson/Documents/Temporal Ecology Lab/hierver/analyses")
}

# read data
tempo_data <- read.csv2("../data/common_variables.csv", header = TRUE)
vassal_data <- read.csv("../data/fromvassal_phenology.csv", header = TRUE, fileEncoding = "latin1")
vassal_data <- subset(vassal_data, !is.na(doy_flo) & !is.na(doy_ver) & doy_flo < doy_ver,
               select = c("variety", "accession_ID", "year", "doy_flo", "doy_ver"))
names(vassal_data) <- c("variety", "clone", "year", "flo", "ver")
vassal_data$flo <- as.integer(vassal_data$flo)
vassal_data$ver <- as.integer(vassal_data$ver)

# assign tempo_dataological events
tempo_data$stagesimple <- NA
tempo_data$stagesimple[tempo_data$stage_description == "véraison : les baies deviennent molles au toucher"] <- "veraison"
tempo_data$stagesimple[tempo_data$stage_description == "mi-floraison : 50% des capuchons floraux sont tombés"] <- "flowering"
tempo_data$stagesimple[tempo_data$stage_description == "Veraison date (50%)"] <- "veraison"
tempo_data$stagesimple[tempo_data$stage_description == "stade avancé de la maturation ou de la coloration du fruit - g stade pâteux du caryopse - 50% ds fruits ont atteint la maturité"] <- "veraison"
tempo_data$stagesimple[tempo_data$stage_description == "Environ 50% des fleurs épanouies"] <- "flowering"
tempo_data$stagesimple[tempo_data$stage_description == "Au moins 50% des fruits sont matures"] <- "veraison"
tempo_data$stagesimple[tempo_data$stage_description == "Veraison date (50%)"] <- "veraison"
tempo_data$stagesimple[tempo_data$stage_description == "Flowering date (50%)"] <- "flowering"

# remove data with missing doy and get mean for each unique combination of
# location, time, variety, and clone
tempo_data <- subset(tempo_data, !is.na(day_of_year))
tempo_data_mean <- aggregate(tempo_data["day_of_year"],
                       tempo_data[c("year", "site_name", "site_latitude", "site_longitude",
                              "variety", "clone_or_accession", "stagesimple")],
                       FUN = mean)

# rearrange data to get doy for flowering and veraison for each unique combination of
# location, time, variety, and clone
flo_ver <- reshape(tempo_data_mean,
                   idvar = c("year", "site_name", "site_latitude", "site_longitude",
                             "variety", "clone_or_accession"),
                   timevar = "stagesimple",
                   direction = "wide")

# remove data with missing doy for either flowering or veraison
# remove flowering doy later than veraison doy e.g. INRA Bergheim, Pinot noir, 2002
flo_ver <- subset(flo_ver, !is.na(day_of_year.flowering) & !is.na(day_of_year.veraison))
flo_ver <- subset(flo_ver, day_of_year.flowering <= day_of_year.veraison)
flo_ver$day_of_year.flowering <- as.integer(flo_ver$day_of_year.flowering)
flo_ver$day_of_year.veraison <- as.integer(flo_ver$day_of_year.veraison)

# retrieve data from sites with more than 50 years of total data (all varieties, all clones)
flo_ver_site <- aggregate(flo_ver["year"],
                          flo_ver[c("site_name", "site_latitude", "site_longitude")],
                          FUN = length)
flo_ver_site_50 <- subset(flo_ver_site, year > 50)
# manually adding shorter names of the sites with >50 (old code)
flo_ver_site_50$location <- c("Bellay", "Lafite", "Vassal", "Vassal", "Bergheim", "Bergheim")

# write varietal data from sites with more than 50 years of total data
flo_ver <- merge(flo_ver, flo_ver_site_50, by = c("site_name", "site_latitude", "site_longitude"),
                 all.x = TRUE)
flo_ver <- subset(flo_ver, variety != "" & site_name %in% flo_ver_site_50$site_name,
                  select = c("variety", "clone_or_accession", "site_name", "location",
                             "site_longitude", "site_latitude", "year.x",
                             "day_of_year.flowering", "day_of_year.veraison"))
colnames(flo_ver) <- c("variety", "clone", "site", "location", "longitude", "latitude",
                       "year", "flo", "ver")
flo_ver <- flo_ver[order(flo_ver$variety, flo_ver$site, flo_ver$year, flo_ver$clone),]

# site data summary
sites <- unique(flo_ver$site)
site_info <- list()
site_info_by_year <- list()

for(i in 1:length(sites)){
  site_data <- flo_ver[which(flo_ver$site == sites[i]),]
  site_data_by_year <- as.data.frame(table(site_data$year))
  names(site_data_by_year) <- c("year", "rowsdata")
  site_data_by_year$site <- sites[i]
  site_info_by_year[[i]] <- site_data_by_year
  site_info[[i]] <- data.frame(site = sites[i], start_year = min(site_data$year),
                               end_year = max(site_data$year),
                               total_years = nrow(site_data_by_year),
                               range = max(site_data$year) - min(site_data$year) + 1)
}
site_info <- do.call(rbind, site_info)
site_info_by_year <- do.call(rbind, site_info_by_year)

# merge INRA data to Bergheim data
berg_data <- flo_ver[which(flo_ver$site == "Bergheim Burlenberg"),]
inra_data <- flo_ver[which(flo_ver$site == "INRA Colmar Domaine de Bergheim"),] #253

inra_include <- rep(TRUE, nrow(inra_data))
for(i in 1:nrow(inra_data)){
  for(j in 1:nrow(berg_data)){
    if(inra_data$variety[i] == berg_data$variety[j] &
       inra_data$year[i] == berg_data$year[j]){
      inra_include[i] <- FALSE
      break
    }
  }
}
inra_include <- inra_data[inra_include,] #y = 111, vy = 171, vcy = 253

flo_ver <- subset(flo_ver, site != "INRA Colmar Domaine de Bergheim")
flo_ver <- rbind(flo_ver, inra_include)

# merge INRA data to Vassal data from TEMPO
vue_data <- flo_ver[which(flo_ver$site == "Vassal-UE"),]
inra_data <- flo_ver[which(flo_ver$site == "INRA Domaine de Vassal"),] #160

inra_include <- rep(TRUE, nrow(inra_data))
for(i in 1:nrow(inra_data)){
  for(j in 1:nrow(vue_data)){
    if(inra_data$variety[i] == vue_data$variety[j] &
       inra_data$year[i] == vue_data$year[j]){
      inra_include[i] <- FALSE
      break
    }
  }
}
sum(inra_include)
inra_include <- inra_data[inra_include,] #y = 11, vy = 16, vcy = 160

flo_ver <- subset(flo_ver, site != "INRA Domaine de Vassal")
flo_ver <- rbind(flo_ver, inra_include)

# merge INRA data to Vassal data from Vassal server
vassal_data <- data.frame(variety = vassal_data$variety, clone = vassal_data$clone, site = "Vassal",
                          location = "Vassal", longitude = NA, latitude = NA, year = vassal_data$year,
                          flo = vassal_data$flo, ver = vassal_data$ver)
tempo_data <- flo_ver[which(flo_ver$location == "Vassal"),] #1051

tempo_include <- rep(TRUE, nrow(tempo_data))
for(i in 1:nrow(tempo_data)){
  for(j in 1:nrow(vassal_data)){
    if(tempo_data$variety[i] == vassal_data$variety[j] &
       tempo_data$clone[i] == vassal_data$clone[j] &
       tempo_data$year[i] == vassal_data$year[j]){
      tempo_include[i] <- FALSE
      break
    }
  }
}
sum(tempo_include)
tempo_include <- tempo_data[tempo_include,] #y = 11, vy = 154, vcy = 167

flo_ver <- subset(flo_ver, location != "Vassal")
flo_ver <- rbind(flo_ver, vassal_data)
flo_ver <- rbind(flo_ver, tempo_include)

write.csv(flo_ver, "output/formatted_data/pheno_data_latlong.csv", row.names = FALSE)
flo_ver <- flo_ver[,c("variety", "clone", "location", "year", "flo", "ver")]
colnames(flo_ver)[3] <- "site"
write.csv(flo_ver, "output/formatted_data/pheno_data.csv", row.names = FALSE)
