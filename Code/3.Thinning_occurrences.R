################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Spatial thinning of species occurrences
################

library("ellipsenm")

setwd("WORKING DIRECTORY")

sp <- list.files(path = "Clean/Pieridae/", pattern = ".csv$", full.names = T)
nam <- list.files(path = "Clean/Pieridae/", pattern = ".csv$", full.names = F)
Finalnam <- paste0("Thinned/", nam)

for (i in 1:length(sp)){
  occurrences <- read.csv(sp[i])
  
  thinned <-thin_data(occurrences, longitude = "Longitude", 
                      latitude = "Latitude", thin_distance = 15, 
                      save = T, name = Finalnam[i])
  }

