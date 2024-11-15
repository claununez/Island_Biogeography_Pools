################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Creating PAM from ellipses
################

# Installing and loading packages
library(biosurvey)
library(raster)
library(rgdal)

# Establishing working directory
setwd("WORKING DIRECTORY")

##### PAM by 50km grid
# Reading the median of the models
models <- list.files(path = "Models", 
                     pattern = "^mean_suitability_calibration_.*tif$",
                     full.names = TRUE, recursive = T)
models1 <- list.files(path = "Models", 
                     pattern = "^mean_suitability_calibration_.*tif$",
                     full.names = F, recursive = T)
#models1 <- gsub("/mean_suitability_calibration_.*tif$", "", models1)
models1 <- gsub(".csv/mean_suitability_calibration_.*tif$", "", models1)

models2 <- stack(models)
models2 <- models2 > 0

for (i in 1:length(models1)){
  writeRaster(models2[[i]], filename = paste0("Median/", models1[i], ".tif"), 
              format = "GTiff")
}

# Reading region of interest
region <- readOGR(dsn = "Shapefiles/Country1.shp")

# Reading environmental variable (2.5 resolution) used to create the master_matrix
variables <- stack(raster("Shapefiles/variable/bio_1.tif"))

region@proj4string <- variables@crs

# Creating master matrix object
m_matrix <- new_master_matrix(data_matrix = data.frame(ID = 1, var = 2), 
                              region = region, raster_base = variables[[1]])

region2 <- spTransform(region, m_matrix$region@proj4string)

# Getting species p-a by country
rasdata <-  files_2data(path = "Median", format = "GTiff")
sppoints <- sp::SpatialPointsDataFrame(rasdata[, 1:2], data = rasdata, 
                                       proj4string = sp::CRS("+init=epsg:4326"))

sppoints@data[, 3] <- gsub(".csv", "", sppoints@data[, 3])
sppoints <- data.frame(ID = sp::over(methods::as(sppoints, "SpatialPoints"), 
                                     region2[, "CNTRY_NAME"]), 
                       Species = sppoints@data[, 3])

sppoints <- PAM_from_table(sppoints, ID_column = "CNTRY_NAME", 
                            species_column = "Species")
sppoints[is.na(sppoints)] <- 0


# Saving results
write.csv(sppoints, "Results/PAM_from_ellipses.csv", row.names = F)
