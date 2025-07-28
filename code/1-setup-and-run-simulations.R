library(uuid)
library(RangeShiftR)
library(readxl)
library(openxlsx)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ips-shift-functions.R")

## create mock basic landscapes ####

# Create a 4x4 raster with 25m resolution
r <- rast(nrows=4, ncols=4, xmin=0, xmax=100, ymin=0, ymax=100, crs="EPSG:27700")

# Assign values: one cell with 100, others NA
values(r) <- rep(100, ncell(r))
#values(r)[1] <- 100

# Set NA flag
NAflag(r) <- -9999

# View raster summary
r
plot(r)
# Write to ASCII
writeRaster(r, "../data/processed/hab_map_qual.asc", filetype="AAIGrid", NAflag=-9999, overwrite=TRUE)

# Write to GeoTIFF
writeRaster(r, "../data/processed/hab_map_qual.tif", filetype="GTiff", NAflag=-9999, overwrite=TRUE)



# Create a 4x4 raster with 25m resolution
r_cost <- rast(nrows=4, ncols=4, xmin=0, xmax=100, ymin=0, ymax=100, crs="EPSG:27700")

# Assign values: one cell with 100, others NA
values(r_cost) <- rep(1, ncell(r_cost))
#values(r)[1] <- 100

# Set NA flag
NAflag(r_cost) <- -9999

# View raster summary
r_cost
plot(r_cost)
# Write to ASCII
writeRaster(r_cost, "../data/processed/cost_map.asc", filetype="AAIGrid", NAflag=-9999, overwrite=TRUE)

# Write to GeoTIFF
writeRaster(r_cost, "../data/processed/cost_map.tif", filetype="GTiff", NAflag=-9999, overwrite=TRUE)


############# set up batch processing ####
  

  
  library(parallel)
  
  cl <- parallel::makeCluster(detectCores())
  
  batches <- read.csv("../code/batches.csv")
  
  phenips = batches$Phenips
  Batch_num = batches$Batch.name
  Inds_file_name =  batches$Inds_file_name
  habitat_rst_name = batches$habitat_rst_name
  cost_map_name = batches$cost_map_name
  reps = batches$reps
  dens_inflection = batches$dens_inflection
  dens_slope = batches$dens_slope
  DP = batches$DP
  PR = batches$PR
  MemSize= batches$MemSize
  GoalType= batches$GoalType
  GoalBias= batches$GoalBias
  AlphaDB= batches$AlphaDB
  BetaDB= batches$BetaDB
  

  
  envir <- environment()
  parallel::clusterExport(cl, varlist = ls(envir), envir = envir)

  
  
  ### * Start parallel processing run ----
  

  # Load required packages on all cluster nodes
  clusterEvalQ(cl, {
    library(terra)
    library(sf)
    library(uuid)
    library(RangeShiftR)
    library(readxl)
    library(openxlsx)
  })
  
  
  time <- system.time(clusterMap(cl,
                                 Run_sims,
                                 Inds_file_name,
                                 habitat_rst_name,
                                 cost_map_name,
                                 phenips,
                                 Batch_num,
                                 reps,
                                 dens_inflection,
                                 dens_slope,
                                 DP,
                                 PR,
                                 MemSize,
                                 GoalType,
                                 GoalBias,
                                 AlphaDB,
                                 BetaDB))
  
  
  parallel::stopCluster(cl)
  
  
  
  
  

  
