################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Initial cleaning of occurrence data
         ## Excluding: records with no coordinates, duplicates, and records  
         ##            with coordinates 0°,0°.

# Data needed:
# Species occurrences 
# The species occurrence file must contain the following columns 
# (in that order): ID, Species_name, Longitude, Latitude.
################

# Defining working directory
setwd("WORKING DIRECTORY") # change this to your working directory

# Pieridae
D <- list.files(path = "Data/GBIF_Pieridae", pattern = ".csv$", full.names = T)
nam <- list.files(path = "Data/GBIF_Pieridae", pattern = ".csv$", 
                  full.names = F)
Finalnam <- paste0("Data/Clean/Pieridae/", nam)

for (i in 1:length(D)){
  occurrences <- read.csv(D[i]) # occurrences 
  
  # Keeping only columns of interest
  occurrences <- occurrences[, c("scientificName", "decimalLongitude", 
                                 "decimalLatitude")]
  colnames(occurrences) <- c("Species", "Longitude", "Latitude")
  
  # Excluding duplicates
  occurrences$code <-  paste(occurrences$Species, occurrences$Longitude, # concatenating columns of interest
                             occurrences$Latitude, sep = "_")
  
  occurrences <- occurrences[!duplicated(occurrences$code), 1:4] # erasing duplicates
  occurrences <- na.omit(occurrences[, 1:3])
  
  # Excluding records with 0°,0° coordinates
  occurrences <- occurrences[occurrences$Longitude != 0 & occurrences$Latitude != 0, ]
  
  # Saving the new set of occurrences
  write.csv(occurrences, Finalnam[i], row.names = FALSE)
}

#-------------------------------------------------------------------------------
# Papilionidae
D1 <- list.files(path = "Data/GBIF_Papilionidae", pattern = ".csv$", 
                 full.names = T)
nam1 <- list.files(path = "Data/GBIF_Papilionidae", pattern = ".csv$", 
                   full.names = F)
Finalnam1 <- paste0("Data/Clean/Papilionidae/", nam1)

for (i in 1:length(D1)){
  occurrences <- read.csv(D1[i]) # occurrences 
  
  # Keeping only columns of interest
  occurrences <- occurrences[, c("scientificName", "decimalLongitude", "decimalLatitude")]
  colnames(occurrences) <- c("Species", "Longitude", "Latitude")
  
  # Excluding duplicates
  occurrences$code <-  paste(occurrences$Species, occurrences$Longitude, # concatenating columns of interest
                             occurrences$Latitude, sep = "_")
  
  occurrences <- occurrences[!duplicated(occurrences$code), 1:4] # erasing duplicates
  occurrences <- na.omit(occurrences[, 1:3])
  
  # Excluding records with 0°,0° coordinates
  occurrences <- occurrences[occurrences$Longitude != 0 & occurrences$Latitude != 0, ]
  
  # saving the new set of occurrences inside continents and area of interest
  write.csv(occurrences, Finalnam1[i], row.names = FALSE)
}
