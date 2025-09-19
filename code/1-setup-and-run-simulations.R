library(parallel)
library(RangeShiftR)
library(yaml)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ips-shift-functions.R")


# Load the base config
# This config allows you to specify default values, so that for the batch simulations
# you only need to specify parameters that change from those defaults.
base_config <- yaml.load_file("../config-SMS.yml")

# Load the batches data, and convert to an array so
# processing can be parallelised
batches <- read.csv("../data/batches.csv")
batches_array <- split(batches, seq(nrow(batches)))

output_path <- "../outputs"

# USED FOR DEBUGGING:
# IF FALSE, ALL BATCHES WILL BE RUN USING A CLUSTER;
# IF TRUE, ONLY 1 SIMULATION FROM THE BATCH WILL BE RUN IN THE CURRENT R SESSION.
DEBUG <- FALSE


if (DEBUG) {

  # Run only the first simulation for debugging
  run_sims(
    batches_array$`1`,
    base_config,
    output_path
  )

} else {

  # Set the number of nodes as half of available cores.
  n_nodes <- as.integer(detectCores()/2)
  
  cl <- parallel::makeCluster(n_nodes)
  
  envir <- environment()
  parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
  
  clusterEvalQ(cl, {
    library(sf)
    library(terra)
    library(RangeShiftR)
  })
  
  time <- system.time(clusterMap(
    cl,
    fun=run_sims,
    x=batches_array,
    MoreArgs=list(
      base_config=base_config,
      output_path=output_path
    )
  ))
    
  parallel::stopCluster(cl)

}





##############################
##############################
# create mock basic landscapes
if (FALSE) {
  
  for (qual in c(100, 80, 60, 40, 20)) {
    # Create a 4x4 raster with 25m resolution
    r <- rast(nrows=4, ncols=4, xmin=0, xmax=100, ymin=0, ymax=100, crs="EPSG:27700")
    
    # Assign values: one cell with qual, others NA
    values(r) <- rep(qual, ncell(r))
    NAflag(r) <- -9999
    
    writeRaster(r, sprintf("../data/processed/hab_map_qual_%s.asc", qual), filetype="AAIGrid", NAflag=-9999, overwrite=TRUE)
    writeRaster(r, sprintf("../data/processed/hab_map_qual_%s.tif", qual), filetype="GTiff", NAflag=-9999, overwrite=TRUE)
  }
  
  for (cost in c(1, 2, 3, 4, 5)) {
    # Create a 4x4 raster with 25m resolution
    r <- rast(nrows=4, ncols=4, xmin=0, xmax=100, ymin=0, ymax=100, crs="EPSG:27700")
    
    # Assign values: one cell with qual, others NA
    values(r) <- rep(cost, ncell(r))
    NAflag(r) <- -9999
    
    writeRaster(r, sprintf("../data/processed/cost_map_%s.asc", cost), filetype="AAIGrid", NAflag=-9999, overwrite=TRUE)
    writeRaster(r, sprintf("../data/processed/cost_map_%s.tif", cost), filetype="GTiff", NAflag=-9999, overwrite=TRUE)
  }
}
