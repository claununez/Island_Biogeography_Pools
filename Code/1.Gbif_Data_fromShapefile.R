################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Download all occurrences of Papilionidae and Pieridae from GBIF 
#          from the area of interest 
################

# Installing and loading packages
install.packages("rgbif")
library(rgbif)

setwd("WORKING DIRECTORY")

# Get GBIF keys for family
famKey <- name_backbone("Papilionidae")$usageKey  #Add your taxon
famKey1 <- name_backbone("Pieridae")$usageKey  #Add your taxon

# List GADM_GID for each state or country we need data for 
GADM_ids <- c("BRB","BHS", "BLZ", "BRA","COL", "CRI", "CUB", "DMA", "DOM",
              "GUF", "GRD", "GLP", "GTM", "GUY", "HTI", "HND", "JAM", "MTQ",
              "MEX", "SUR", "NLD", "NIC", "PAN", "PRI", "LCA", "TTO", "USA",
              "VCT", "VEN")

# Perform a download for your desired taxon and retrieve a download key. 
# You can set download filter parameters using pred_in and pred functions
gbif_download_key = occ_download(
  pred("taxonKey", famKey), # insert taxon key for the taxon interested in
  pred_in("gadm",GADM_ids),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  # Specify your GBIF data download user details.
  user = "YOUR USER NAME", 
  pwd = "YOUR PASSWORD", 
  email = "YOUR EMAIL"
)

gbif_download_key1 = occ_download(
  pred("taxonKey", famKey1), # insert taxon key for the taxon interested in
  pred_in("gadm",GADM_ids),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  # Specify your GBIF data download user details.
  user = "YOUR USER NAME", 
  pwd = "YOUR PASSWORD", 
  email = "YOUR EMAIL"
)

# This download may take some time to finish. Check its status with this command. 
occ_download_wait(gbif_download_key)
occ_download_wait(gbif_download_key1)

data_download <- occ_download_get(gbif_download_key, overwrite = TRUE) %>%
  occ_download_import()
data_download1 <- occ_download_get(gbif_download_key1, overwrite = TRUE) %>%
  occ_download_import()

# Papilionidae
pap <- unique(data_download$species)
pap_p <- for (i in 1:length(pap)){
  fila <- data_download[data_download$species == pap[i], ]
  write.csv(fila, paste0("GBIF_Papilionidae/", pap[i], ".csv"), 
            row.names = F)
}

# Pieridae
pie <- unique(data_download1$species)
pie_p <- for (i in 1:length(pie)){
  fila1 <- data_download1[data_download1$species == pie[i], ]
  write.csv(fila1, paste0("GBIF_Pieridae/", pie[i], ".csv"), 
            row.names = F)
}
