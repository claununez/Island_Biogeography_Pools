################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: 
################

#Loading packages
library(terra)
library(ellipse)

#Establishing working directory
setwd("WORKING DIRECTORY")


#Reading environmental variables
a5 <- rast("Rasters/bio_5.tif")
a6 <- rast("Rasters/bio_6.tif")
a13 <- rast("Rasters/bio_13.tif")
a14 <- rast("Rasters/bio_14.tif")

l <- list(a5, a6, a13, a14)
st <- rast(l)
d <- dim(st)[[3]]

# Below is for pool = Countries
dt <- read.csv("all_Bats_thinned.csv") 

#keep species with 10 or more dif coordinates
tt <- table(dt[, 1]) #A tally of how many records per species
ii_many <- which(tt > 5)
length(ii_many)

# Getting the names of the species in a vector of characters
nms <- names(tt[ii_many])

s <- length(ii_many) # number of species with more than 9 occurrences
mus <- matrix(0, ncol = d, nrow = s) # as  many centroids of 'd'dimensions as 
                                     # species with enough points (ii_many)
cvrs <- matrix(0, ncol = d^2, nrow = s) # variances-covariances as vectors

# With a little progress bar
pb <- winProgressBar(title = "progress bar", min = 0, max = s, width = 300)
for(j in 1:s){
  iii <- which(dt[, 1] == nms[j]);
  
  vrs3 <- extract(st, dt[iii, 2:3]) # this extracts the raster values from the 
                                    # coordinates in each one of the iii columns
  vrs2 <- vrs3[, -1]
  
  # remove NAs
  ina <- which(is.na(vrs2))
  if(identical(ina, integer(0))) {
    vrs = vrs2
    } 
  else {
    vrs = vrs2[-ina, ]
    } # notice treatment of 'ina' in case it is 'integer(0)'
  mus[j, ] <- colMeans(vrs, na.rm = TRUE) # the i-th mean
  cvrs[j, ] <- as.vector(cov(vrs)) # the i-th covariance matrix (as a vector)
  setWinProgressBar(pb, j, title =  paste(round(j/s * 100, 0),"% done"))
}
close(pb)

# Check for singularities creating a vector "dts"
# of the determinants
epsilon <- 0.00001
dts <- vector(mode = "numeric", length = s)
for(j in 1:s){
  dts[j] <- abs(det(matrix(cvrs[j, ], ncol = d, nrow = d)))
}
nosings <- which(ifelse(dts < epsilon, 1, 0) == 0) #notice we are checking for dts > epsilon, i.e., for NOT-singular!
musM <- mus[nosings,]
cvrsM <- cvrs[nosings, ]
nmsM <- nms[nosings]
s <- length(nosings)

# Creating a list of covariance matrices
myList <- vector("list", s) #Empty list
for(j in 1:s){
  myList[[j]] = matrix(cvrsM[j, ], ncol = d, nrow = d) # Matrices in the list
}

# Reading shapefiles of the islands
ecrgs <- vect("Shapefiles/PeriOnly.shp")
count <- vect("Shapefiles/country.shp")
leew <- vect("Shapefiles/Leeward.shp")
wind <- vect("Shapefiles/Windward.shp")
jam <- vect("Shapefiles/Jamaica.shp")
Bvrg <- vect("Shapefiles/BritVirginIls.shp")
vrg <- vect("Shapefiles/VirginIls.shp")
crs(Bvrg) <- crs(vrg)

plot(ecrgs, asp = 1)
lines(leew, col = "red")
lines(jam, col = "green")
lines(Bvrg, col = "blue")
lines(wind, col = "gray")

# Obtaining clouds of E-points for every region
clEcrgs <- mask(st, ecrgs)
clEcrgs2 <- as.points(clEcrgs, values = TRUE, na.rm = TRUE)
cloudE <- values(clEcrgs2)

clCoun <- mask(st, count)
clCoun2 <- as.points(clCoun, values = TRUE, na.rm = TRUE)
cloudC <- values(clCoun2)

clleew <- mask(st, leew)
clleew2 <- as.points(clleew, values = TRUE, na.rm = TRUE)
cloudL <- values(clleew2)

clwind <- mask(st, wind)
clwind2 <- as.points(clwind, values = TRUE, na.rm = TRUE)
cloudW <- values(clwind2)

clJam <- mask(st, jam)
clJam2 <- as.points(clJam, values = TRUE, na.rm = TRUE)
cloudJ <- values(clJam2)

clBvrg <- mask(st, Bvrg)
clBvrg2 <- as.points(clBvrg, values = TRUE, na.rm = TRUE)
cloudBv <- values(clBvrg2)

# Getting the identity of species of interest
# for Ecoregions as pool
Ps <- which(nmsM == "Phoebis_sennae")
Am <- which(nmsM == "Ascia_monuste")
Ed <- which(nmsM == "Eurema_daira")
El <- which(nmsM == "Pyrisitia_lisa")
Ni <- which(nmsM == "Nathalis_iole")
En <- which(nmsM == "Abaeis_nicippe")

#'ellip' is a function to plot orthogonal projections on 2Dimensional planes
#in other words, I have my d-dimensional niche, and I would like
#to see the projections on 2-dimensional space, i.e, planes.
#h and k are the indices of the columns of the variables to define the plane

#sp is the index of the species name, as in nmsM. axes is myList[[sp]]
#center is musmM[sp,c(h,k)]
#It returns an ellipse object in 2-Dimensions

ellip <- function(sp, h, k, center, axes)
{
  cn <- center[sp, c(h, k)]
  v <- axes[[sp]]
  pr=v[c(h,k),c(h,k)] # this is the submatrix of the two variables of interest
  el <- ellipse(pr, centre = cn, level = 0.95)
  return(el)
}

####Country
jpeg(filename = "County_PreDry_MinTem.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudC[, c(2, 4)], pch = ".", xlab = "Minimum temperature of coldest month", 
     ylab = "Precipitation of driest month", main = "Country", col = "gray") #Entire region
points(cloudJ[, c(2, 4)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(2, 4)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(2, 4)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(2, 4)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps <- ellip(Ps, 2, 4, musM, myList)
lines(esp_Ps, col = "#313695", lwd = 2)#Phoebis sennae
esp_Ed <- ellip(En, 2, 4, musM, myList)
lines(esp_Ed, col = "#a50026", lwd = 2)#Eurema nicippe
esp_El <- ellip(El, 2, 4, musM, myList)
lines(esp_El, col = "#fdae61", lwd = 2)#Eurema lisa
esp_Ni <- ellip(Ni, 2, 4, musM, myList)
lines(esp_Ni, col = "black", lwd = 2)#Nathalis iole
esp_Ed <- ellip(Ed, 2, 4, musM, myList)
lines(esp_Ed, col = "#74add1", lwd = 2)#Eurema daira
legend("topleft", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                             "Phoebis sennae", "Eurema daira","Eurema nicippe", "Eurema lisa", 
                             "Nathalis iole"), text.font = c(1, 1, 1, 1, 3, 3, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026", "#fdae61", "black"),  bty = "n")
dev.off()


jpeg(filename = "County_PreWet_MaxTem.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudC[, c(1, 3)], pch = ".", xlab = "Maximum temperature of warmest month", 
     ylab = "Precipitation of wettest month", main = "Country", col = "gray") #Entire region
points(cloudJ[, c(1, 3)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(1, 3)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(1, 3)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(1, 3)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps1 <- ellip(Ps, 1, 3, musM, myList)
lines(esp_Ps1, col = "#313695", lwd = 2)#Phoebis sennae
esp_En1 <- ellip(En, 1, 3, musM, myList)
lines(esp_En1, col = "#a50026", lwd = 2)#Eurema nicippe
esp_El1 <- ellip(El, 1, 3, musM, myList)
lines(esp_El1, col = "#fdae61", lwd = 2)#Eurema lisa
esp_Ni1 <- ellip(Ni, 1, 3, musM, myList)
lines(esp_Ni1, col = "black", lwd = 2)#Nathalis iole
esp_Ed1 <- ellip(Ed, 1, 3, musM, myList)
lines(esp_Ed1, col = "#74add1", lwd = 2)#Eurema daira
legend("topleft", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                             "Phoebis sennae", "Eurema daira","Eurema nicippe", "Eurema lisa", 
                             "Nathalis iole"), text.font = c(1, 1, 1, 1, 3, 3, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026", "#fdae61", "black"),  bty = "n")
dev.off()


jpeg(filename = "County_PreDry_PreWet.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudC[, c(1, 4)], pch = ".", xlab = "Precipitation of wettest month", 
     ylab = "Precipitation of driest month", main = "Country", col = "gray") #Entire region
points(cloudJ[, c(1, 4)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(1, 4)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(1, 4)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(1, 4)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps2 <- ellip(Ps, 1, 4, musM, myList)
lines(esp_Ps2, col = "#313695", lwd = 2)#Phoebis sennae
esp_Ed2 <- ellip(En, 1, 4, musM, myList)
lines(esp_Ed2, col = "#a50026", lwd = 2)#Eurema nicippe
esp_El2 <- ellip(El, 1, 4, musM, myList)
lines(esp_El2, col = "#fdae61", lwd = 2)#Eurema lisa
esp_Ni2 <- ellip(Ni, 1, 4, musM, myList)
lines(esp_Ni2, col = "black", lwd = 2)#Nathalis iole
esp_Ed2 <- ellip(Ed, 1, 4, musM, myList)
lines(esp_Ed2, col = "#74add1", lwd = 2)#Eurema daira
legend("topleft", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                             "Phoebis sennae", "Eurema daira","Eurema nicippe", "Eurema lisa", 
                             "Nathalis iole"), text.font = c(1, 1, 1, 1, 3, 3, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026", "#fdae61", "black"),  bty = "n")
dev.off()


jpeg(filename = "County_MinTem_MaxTem.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudC[, c(3, 2)], pch = ".", xlab = "Maximum temperature of warmest month", 
     ylab = "Minimum temperature of coldest month", main = "Country", col = "gray") #Entire region
points(cloudJ[, c(3, 2)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(3, 2)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(3, 2)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(3, 2)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps3 <- ellip(Ps, 3, 2, musM, myList)
lines(esp_Ps3, col = "#313695", lwd = 2)#Phoebis sennae
esp_En3 <- ellip(En, 3, 2, musM, myList)
lines(esp_En3, col = "#a50026", lwd = 2)#Eurema nicippe
esp_El3 <- ellip(El, 3, 2, musM, myList)
lines(esp_El3, col = "#fdae61", lwd = 2)#Eurema lisa
esp_Ni3 <- ellip(Ni, 3, 2, musM, myList)
lines(esp_Ni3, col = "black", lwd = 2)#Nathalis iole
esp_Ed3 <- ellip(Ed, 3, 2, musM, myList)
lines(esp_Ed3, col = "#74add1", lwd = 2)#Eurema daira
legend("bottomright", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                                 "Phoebis sennae", "Eurema daira","Eurema nicippe", "Eurema lisa", 
                                 "Nathalis iole"), text.font = c(1, 1, 1, 1, 3, 3, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026", "#fdae61", "black"),  bty = "n")
dev.off()


####Ecoregions
jpeg(filename = "Eco_PreDry_MinTem.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudE[, c(2, 4)], pch = ".", xlab = "Minimum temperature of coldest month", 
     ylab = "Precipitation of driest month", main = "Ecoregion", col = "gray") #Entire region
points(cloudJ[, c(2, 4)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(2, 4)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(2, 4)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(2, 4)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps_e <- ellip(Ps_e, 2, 4, musM_e, myList_e)
lines(esp_Ps_e, col = "#313695", lwd = 2)#Phoebis sennae
esp_Am_e <- ellip(Am_e, 2, 4, musM_e, myList_e)
lines(esp_Am_e, col = "#a50026", lwd = 2)#Ascia monuste
esp_Ed_e <- ellip(Ed_e, 2, 4, musM_e, myList_e)
lines(esp_Ed_e, col = "#74add1", lwd = 2)#Eurema daira
legend("topleft", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                             "Phoebis sennae", "Eurema daira","Ascia monuste"), 
       text.font = c(1, 1, 1, 1, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026"),  bty = "n")
dev.off()


jpeg(filename = "Eco_PreWet_MaxTem.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudE[, c(1, 3)], pch = ".", xlab = "Maximum temperature of warmest month", 
     ylab = "Precipitation of wettest month", main = "Ecoregion", col = "gray") #Entire region
points(cloudJ[, c(1, 3)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(1, 3)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(1, 3)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(1, 3)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps_e1 <- ellip(Ps_e, 1, 3, musM_e, myList_e)
lines(esp_Ps_e1, col = "#313695", lwd = 2)#Phoebis sennae
esp_Am_e1 <- ellip(Am_e, 1, 3, musM_e, myList_e)
lines(esp_Am_e1, col = "#a50026", lwd = 2)#Ascia monuste
esp_Ed_e1 <- ellip(Ed_e, 1, 3, musM_e, myList_e)
lines(esp_Ed_e1, col = "#74add1", lwd = 2)#Eurema daira
legend("topleft", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                             "Phoebis sennae", "Eurema daira","Ascia monuste"), 
       text.font = c(1, 1, 1, 1, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026"),  bty = "n")
dev.off()


jpeg(filename = "Eco_PreDry_PreWet.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudE[, c(1, 4)], pch = ".", xlab = "Precipitation of wettest month", 
     ylab = "Precipitation of driest month", main = "Ecoregion", col = "gray") #Entire region
points(cloudJ[, c(1, 4)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(1, 4)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(1, 4)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(1, 4)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps_e2 <- ellip(Ps_e, 1, 4, musM_e, myList_e)
lines(esp_Ps_e2, col = "#313695", lwd = 2)#Phoebis sennae
esp_Am_e2 <- ellip(Am_e, 1, 4, musM_e, myList_e)
lines(esp_Am_e2, col = "#a50026", lwd = 2)#Ascia monuste
esp_Ed_e2 <- ellip(Ed_e, 1, 4, musM_e, myList_e)
lines(esp_Ed_e2, col = "#74add1", lwd = 2)#Eurema daira
legend("topleft", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                             "Phoebis sennae", "Eurema daira","Ascia monuste"), 
       text.font = c(1, 1, 1, 1, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026"),  bty = "n")
dev.off()


jpeg(filename = "Eco_MinTem_MaxTem.jpg",
     width = 20, height = 20, units = "cm", pointsize = 14,
     quality = 100, res = 600)
plot(cloudE[, c(3, 2)], pch = ".", xlab = "Maximum temperature of warmest month", 
     ylab = "Minimum temperature of coldest month", main = "Ecoregion", col = "gray") #Entire region
points(cloudJ[, c(3, 2)], pch = 19, col = "#fee090") #Jamaica
points(cloudW[, c(3, 2)], pch = 19, col = "#d73027") #Windward
points(cloudL[, c(3, 2)], pch = 19, col = "#4575b4") #Leeward
points(cloudBv[, c(3, 2)], pch = 19, col = "#f46d43") #Virgin Islands

esp_Ps_e3 <- ellip(Ps_e, 3, 2, musM_e, myList_e)
lines(esp_Ps_e3, col = "#313695", lwd = 2)#Phoebis sennae
esp_Am_e3 <- ellip(Am_e, 3, 2, musM_e, myList_e)
lines(esp_Am_e3, col = "#a50026", lwd = 2)#Ascia monuste
esp_Ed_e3 <- ellip(Ed_e, 3, 2, musM_e, myList_e)
lines(esp_Ed_e3, col = "#74add1", lwd = 2)#Eurema daira
legend("bottomright", legend = c("Jamaica", "Virgin Islands", "Leeward", "Windward",
                                 "Phoebis sennae", "Eurema daira","Ascia monuste"), 
       text.font = c(1, 1, 1, 1, 3, 3, 3),
       fill = c("#fee090", "#f46d43", "#4575b4", "#d73027", "#313695", "#74add1", 
                "#a50026"),  bty = "n")
dev.off()