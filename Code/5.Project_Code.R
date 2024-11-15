################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Ecological Niche Modeling with Ellipsoids
################

# Installing and loading packages
remotes::install_github("marlonecobos/ellipsenm")
install.packages("geodata")
install.packages("terra")
library(terra)
library(ellipsenm)
library(raster)

# Establishing working directory
setwd("WORKING DIRECTORY")

# Calling functions (it have to be on the same working directory)
source("Functions_ellip.R")

# Downloading current variables from worldclim version 2.1
cvars <- get_NWC_bio(period = "historical", res = "2.5m", time = NULL, SSP = NULL, 
                     GCM = NULL, output_dir = "variable") 

vars <- list.files("variables/wc2.1_5m_bio/", pattern = ".tif", full.names = TRUE)
var_stack <- terra::rast(vars)

# Reading shapefile to crop the variables to area of interest
shp <- vect("country.shp")

# Cropping the variables
var <- terra::crop(var_stack, shp, mask = T)

# Getting variable names
name <- names(var)
name1 <- gsub("wc2.1_2.5m_", "", name) 

names(var) <- name1
var_name <- paste0("variable/", name1, ".tif")

# Saving cropped variables
for (i in 1:nlayers(var)){
  writeRaster(var[[i]], var_name[i], format = "GTiff", overwrite = TRUE)
}

# Selecting variables of interest
var_use <- var[[c(15, 16, 5, 6)]]

# Reading variables of interest
var_use1 <- list.files(path = "Shapefiles/variable", pattern = ".tif", 
                       full.names = TRUE)
var_use <- stack(var_use1[c(7, 8, 15, 16)])

######
# Reading occurrences
sp <- list.files(path = "Final", pattern = ".csv$", full.names = TRUE)
sp1 <- list.files(path = "Final", pattern = ".csv$", full.names = FALSE)
sp1 <- paste0("Models/", sp1)

errors <- vector()

for (i in 1:length(sp1)) {
  sp2 <- read.csv(file = sp[i])
  ext <- extract(var_use, sp2[, c(2,3)])
  
  non_na <- which(is.na(ext), arr.ind = T)[, 1]
  if (length(non_na) > 0) {
    sp2 <- sp2[-non_na, ]
    ext <- ext[-non_na, ]
  }
  
  if (nrow(sp2) > 5) {  # only speices that have more than 5 occurrences
    # initial test
    covm <- cov(ext) # ext is data with which ellipsoids will be created
    
    # testing if positive definite
    test <- is_pos_def(covm)
    
    # testing for non_singularity
    test1 <- non_sing(covm)
    
    # conditional running
    if (test == TRUE & test1 == TRUE) {
      # creating the model with no replicates
      err <- try(ellipsoid_model(data = sp2, species = "Species",
                                 longitude = "Longitude", 
                                 latitude = "Latitude",
                                 raster_layers = var_use, method = "covmat", 
                                 level = 95, replicates = 10, 
                                 bootstrap_percentage = 75,
                                 prediction = "suitability",
                                 return_numeric = TRUE, format = "GTiff",
                                 overwrite = TRUE, 
                                 output_directory = sp1[i]), silent = TRUE)
      if (class(err) == "try-error") {
        errors[i] <- i 
      }
    } else {
      message("\nNon positive definite or singular cov. matrix for species ", 
              as.character(sp2[1, 1]))
    }
  } else {
    message("\nSpecies ", sp1[i], "less than 5 records.............")
  }
}
