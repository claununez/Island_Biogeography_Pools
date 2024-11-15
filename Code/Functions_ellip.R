# Download worldclim variables version 2.1
#
# argument options
# period: "historical", "future"
# res: 10m, 5m, 2.5m, 30s (30s only for historical)
# time: 2021-2040, 2041-2060, 2061-2080, 2081-2100 (only for future) 
# SSP: 126, 245, 370, 585 (only for future) 
# GCM: BCC-CSM2-MR, CNRM-CM6-1, CNRM-ESM2-1, CanESM5, IPSL-CM6A-LR, MIROC-ES2L, 
#      MIROC6, MRI-ESM2-0 (only for future) 
# output_dir: the directory the user wants to put the results in 

get_NWC_bio <- function(period, res, time = NULL, SSP = NULL, GCM = NULL, 
                        output_dir = NULL) {
  if (is.null(output_dir)) {dir <- getwd()} else {dir <- output_dir}
  if (length(res) > 1) {stop("Argument 'res' must be of length 1.")}
  
  if (period[1] == "historical") {
    if (any(!res %in% c("10m", "5m", "2.5m", "30s"))) {
      stop("Argument 'res' must by any of: 10m, 5m, 2.5m, or 30s")
    }
    
    url <- paste0("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_", 
                  res, "_bio.zip")
    
    dir.create(dir)
    dfile <- paste0(dir, "/wc2.1_", res, "_bio.zip")
    
    download.file(url, destfile = dfile, method = "auto", 
                  quiet = FALSE, mode = "wb", cacheOK = TRUE)
    
    dfol <- paste0(dir, "/wc2.1_", res, "_bio")
    dir.create(dfol)
    unzip(zipfile = dfile, exdir = dfol)
    
    files <- list.files(dfol, pattern = ".tif$", full.names = TRUE, recursive = TRUE)
    return(raster::stack(files))
    
  } else {
    if (any(!res %in% c("10m", "5m", "2.5m"))) {
      stop("Argument 'res' must by any of: 10m, 5m, or 2.5m")
    }
    
    lens <- c(length(time), length(SSP), length(GCM))
    if (any(lens > 1)) {
      stop("Arguments 'time', 'SSP', and 'GCM' must be of length 1.")
    }
    
    url <- paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/", 
                  res, "/wc2.1_", res, "_bioc_", GCM, "_ssp", SSP, "_",
                  time, ".zip")
    
    dir.create(dir)
    dfile <- paste0(dir, "/wc2.1_", res, "_bioc_", GCM, "_ssp", SSP, "_",
                    time, ".zip")
    
    download.file(url, destfile = dfile, method = "auto", 
                  quiet = FALSE, mode = "wb", cacheOK = TRUE)
    
    dfol <- paste0(dir, "/wc2.1_", res, "_bioc_", GCM, "_ssp", SSP, "_", time)
    dir.create(dfol)
    unzip(zipfile = dfile, exdir = dfol)
    
    files <- list.files(dfol, pattern = ".tif$", full.names = TRUE, recursive = TRUE)
    return(raster::stack(files))
  }
}


citation_NWC_bio <- function(period) {
  if (period[1] == "historical") {
    cat("Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: New 1km spatial resolution climate surfaces\n\tfor global land areas. International Journal of Climatology 37 (12): 4302-4315.")
  } else {
    cat("For more information about climatic data derived from the CMIP6 project, please visit:\n\thttps://wcrp-cmip.github.io/CMIP6_CVs/docs/CMIP6_institution_id.html\n\thttp://bit.ly/2gBCuqM")
  }
}

# function to check if matrix is positive definite
is_pos_def <- function(x, tol = 1e-8) {
  eigenvalues <- eigen(x, only.values = TRUE)$values
  eigenvalues <- ifelse(eigenvalues < tol, 0, eigenvalues)
  return(all(eigenvalues > 0))
}

# function to check if matrix is non-singular
non_sing <- function(x, tol = 1e-8) {
  return(det(x) > tol)
}
