################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Jaccard analysis and omission and commission errors
################

# Installing and loading packages
library(vegan)
library(biosurvey)

# Establishing working directory
setwd("WORKING DIRECTORY")

# Data
## PAM from atlas
pam1 <- read.csv("PAM_Antilles.csv") 
pam1 <- pam1[pam1$CNTRY_NAME != "", ]
endem <- read.csv("Endemics/AllCaribbean.csv", header = F)
endem <- endem[, 1]

pam1 <- pam1[, !colnames(pam1) %in% endem]
  
## PAM from models (ecoregions)
pam2 <- read.csv("PAM_from_ellipses_ecoregions.csv") 

pam2 <- pam2[, !colnames(pam2) %in% endem]

# Preparing data for analysis
## Keeping all countries that coincide (pam1 must have all that matter)
countries <- pam1$CNTRY_NAME[pam1$CNTRY_NAME %in% pam2$CNTRY_NAME] 

## Excluding columns that are not required (if needed)
colnames(pam1)
colnames(pam2)

## Keeping only countries or areas that match
pam1 <- pam1[pam1$CNTRY_NAME %in% countries, ]
pam2 <- pam2[pam1$CNTRY_NAME %in% countries, ]

## Fixing species list in both maps (pam2 must have all that matter)
species1 <- colnames(pam1)[-1]
species2 <- colnames(pam2)[-1]

species_1in2 <- species1[species1 %in% species2]

species_match <- species2[species2 %in% species_1in2]

miss_species <- species2[!species2 %in% species_1in2] # missing in pam1

### Only if there are missing ones
for (i in miss_species) {
  pam1 <- cbind(pam1, 0)
  colnames(pam1)[ncol(pam1)] <- i
}

pam1a <- pam1[, c(colnames(pam1)[1], sort(c(species_match, miss_species)))]
pam2a <- pam2[, c(colnames(pam2)[1], sort(c(species_match, miss_species)))]


# Running analyses in loop (comparisons per country of both pams)

# Jaccard
# Omission and commission errors
# (omission = species in atlas not predicted by models)
# (commission = species predicted by models but not in atlas)
jacc_omi_comi <- lapply(countries, function(x) {
  ## compacted matrix to compare pams
  pam_mat <- rbind(pam1a[pam1a[, 1] == x, -1], pam2a[pam2a[, 1] == x, -1])
  
  ## Number of species 
  nsp_p <- sum(pam_mat[1, ] == 1)
  nsp_np <- sum(pam_mat[1, ] == 0)
  
  ## Jaccard
  jac <- c(vegdist(pam_mat, method = "jaccard"))
  
  ## Making it comparable by summing
  pam_mat[1, ] <- pam_mat[1, ] * 10
  sums <- colSums(pam_mat)
  
  ## Omissions
  omi <- colnames(pam_mat)[sums == 10]
  omir <- length(omi) / nsp_p
  
  ## Commissions
  comi <- colnames(pam_mat)[sums == 1]
  comir <- length(comi) / nsp_np
  
  list(jaccard = jac, omission_rate = omir, omissions = omi, 
       commision_rate = comir, commisions = comi)
})

island <- c("Bahamas", "Cuba", "Hispaniola", "Jamaica", "Puerto Rico",
            "Virgin Islands, Leeward", "Windward")

names(jacc_omi_comi) <- island

# Indices and rates
jacc <- sapply(jacc_omi_comi, function(x) {x$jaccard})

omissionr <- sapply(jacc_omi_comi, function(x) {x$omission_rate})

commissionr <- sapply(jacc_omi_comi, function(x) {x$commision_rate})

summary(omissionr)

# Plotting results
bah <- jacc_omi_comi$`Bahamas`$omission_rate
cub <- jacc_omi_comi$Cuba$omission_rate
his <- jacc_omi_comi$Hispaniola$omission_rate
jam <- jacc_omi_comi$Jamaica$omission_rate
puer <- jacc_omi_comi$`Puerto Rico`$omission_rate
virg <- jacc_omi_comi$`Virgin`$omission_rate
lw <- jacc_omi_comi$`Leeward`$omission_rate
ww <- jacc_omi_comi$`Windward`$omission_rate
omiss <- c(bah, cub, his, jam, puer, lw, ww)

bah1 <- jacc_omi_comi$`Bahamas`$commision_rate
cub1 <- jacc_omi_comi$Cuba$commision_rate
his1 <- jacc_omi_comi$Hispaniola$commision_rate
jam1 <- jacc_omi_comi$Jamaica$commision_rate
puer1 <- jacc_omi_comi$`Puerto Rico`$commision_rate
virg1 <- jacc_omi_comi$`Virgin`$commision_rate
lw1 <- jacc_omi_comi$`Leeward`$commision_rate
ww1 <- jacc_omi_comi$`Windward`$commision_rate
comiss <- c(bah1, cub1, his1, jam1, puer1, lw1, ww1)

bah2 <- length(jacc_omi_comi$`Bahamas`$omissions)
cub2 <- length(jacc_omi_comi$Cuba$omissions)
his2 <- length(jacc_omi_comi$Hispaniola$omissions)
jam2 <- length(jacc_omi_comi$Jamaica$omissions)
puer2 <- length(jacc_omi_comi$`Puerto Rico`$omissions)
virg2 <- length(jacc_omi_comi$`Virgin`$omissions)
lw2 <- length(jacc_omi_comi$`Leeward`$omissions)
ww2 <- length(jacc_omi_comi$`Windward`$omissions)
omiss_n <- c(bah2, cub2, his2, jam2, puer2, lw2, ww2)

bah3 <- length(jacc_omi_comi$`Bahamas`$commisions)
cub3 <- length(jacc_omi_comi$Cuba$commisions)
his3 <- length(jacc_omi_comi$Hispaniola$commisions)
jam3 <- length(jacc_omi_comi$Jamaica$commisions)
puer3 <- length(jacc_omi_comi$`Puerto Rico`$commisions)
virg3 <- length(jacc_omi_comi$`Virgin`$commisions)
lw3 <- length(jacc_omi_comi$`Leeward`$commisions)
ww3 <- length(jacc_omi_comi$`Windward`$commisions)
comiss_n <- c(bah3, cub3, his3, jam3, puer3, lw3, ww3)

omicimi <- rbind(omiss, comiss)
colnames(omicimi) <- island

omi <- rbind(omiss_n)
colnames(omi) <- island

comi <- rbind(comiss_n)
colnames(comi) <- island

jpeg(filename = "JacOmComThin_NoTT.jpg",
     width = 28, height = 28, units = "cm", pointsize = 12,
     quality = 100, res = 600)
par(mfrow = c(4, 1), mar = c(4.5, 4.5, 1, 0.5), cex = 0.9)
barplot(jacc, ylab = "", xlab = "", main = "", ylim = c(0, 1))
title(ylab = "Jaccard index", cex = 2)
box(bty = "l")

barplot(omicimi, ylab = "", xlab = "", main = "", 
        beside = T, col = c("gray50", "lightgray"), ylim = c(0, 1.4))
legend("topright", legend = c("Omission", "Comission"), 
       fill = c("gray50", "lightgray"), bty = "n", horiz = T)
title(ylab = "Rate", cex = 2)
box(bty = "l")

barplot(omi, ylab = "", xlab = "", col = "gray50")
title(main = "Omissions", font.main = 1, ylab = "Number of species", 
      cex.main = 1)
box(bty = "l")

barplot(comi, ylab = "", xlab = "Islands", col = "lightgray")
title(main = "Comissions", font.main = 1, ylab = "Number of species", 
      cex.main = 1)
box(bty = "l")
dev.off()

