################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Eliminating the occurrences from the Caribbean Islands
################

setwd("WORKING DIRECTORY")

# loading packages
library(terra)
library(sp)

papi <- list.files(path = "Thinned/Papilionidae",
                   pattern = ".csv$", full.names = TRUE, recursive = T)
papi_names<- list.files(path = "Thinned/Papilionidae",
                        pattern = ".csv$", full.names = F, recursive = T)
papi_names <- paste0("Final/", papi_names)

carib <- rast("Shapefile/country.shp")

for (i in 1:length(papi_names)){
  occ <- read.csv(papi[i])
  occ1 <- occ[, c("acceptedScientificName", "decimalLongitude", "decimalLatitude")]
  
  occ1 <- SpatialPointsDataFrame(coords = occ1[, 2:3], data = occ,
                                 proj4string = carib@proj4string, match.ID = F)
  
  occ1 <- occ1[carib, ]
  
  if (nrow(occ1@data) > 0){
    write.csv(occ1@data, papi_names[i], row.names = F)
  }
}

###-------------------Pieridae-------------------
pier <- list.files(path = "Thinned/Pieridae",
                   pattern = ".csv$", full.names = TRUE, recursive = T)
pier_names<- list.files(path = "Thinned/Pieridae",
                        pattern = ".csv$", full.names = F, recursive = T)
pier_names <- paste0("Final/", pier_names)

for (i in 1:length(pier_names)){
  occ <- read.csv(pier[i])
  occ1 <- occ[, c("scientificName", "decimalLongitude", "decimalLatitude")]
  
  occ1 <- SpatialPointsDataFrame(coords = occ1[, 2:3], data = occ,
                                 proj4string = carib@proj4string, match.ID = F)
  
  occ1 <- occ1[carib, ]
  
  if (nrow(occ1@data) > 0){
    write.csv(occ1@data, pier_names[i], row.names = F)
  }
}
