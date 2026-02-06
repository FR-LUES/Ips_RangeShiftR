


library(Rcpp)
library(devtools)
library(Rdpack)
devtools::install_github("https://github.com/RangeShifter/RangeShiftR-pkg", ref = "main", subdir="RangeShiftR")
library(parallel)
library(RangeShiftR)
library(yaml)


# Set working directory (Databricks workspace path)
#install.packages("terra")
setwd("/dbfs/mnt/lab/unrestricted/Fiona/DASH_files/FoD_2030/")

# Load custom functions
source("/Workspace/Users/fiona.plenderleith@defra.onmicrosoft.com/RS/code/ips-shift-functions.R")

#############
# PARAMETERS
#############

# Base config file with default simulation parameters
BASE_CONFIG_PATH <- "/Workspace/Users/fiona.plenderleith@defra.onmicrosoft.com/RS/config-SMS.yml"
#
# Batch file with parameter variations
BATCHES_FILE_PATH <- "/dbfs/mnt/lab/unrestricted/Fiona/DASH_files/FoD_2030/Data/Batches_FoD_2030_expanded.csv"

# Debug mode: TRUE = run one simulation locally; FALSE = run all in parallel
DEBUG <- FALSE

#############
# Load Config and Batches
#############

# Load base configuration
base_config <- yaml::read_yaml(BASE_CONFIG_PATH)

# Load batch definitions and split into list for iteration
batches <- read.csv(BATCHES_FILE_PATH)
dynamic <- read.csv("/dbfs/mnt/lab/unrestricted/Fiona/DASH_files/FoD_2030/Data/dynamic_landscapes_2030.csv")
# Define the new base path
new_base <- "/dbfs/mnt/lab/unrestricted/Fiona/DASH_files/FoD_2030/Data"

# Function to replace leading ../data with your absolute path
replace_data_prefix <- function(x) {
  if (is.character(x)) {
    # Replace ONLY when the string starts with ../data
    sub("^\\.\\./data", new_base, x)
  } else {
    x
  }
}

# Apply to all character columns
batches[] <- lapply(batches, replace_data_prefix)
dynamic[] <- lapply(dynamic, replace_data_prefix)
write.csv(dynamic, "/dbfs/mnt/lab/unrestricted/Fiona/DASH_files/FoD_2030/Data/dynamic_landscapes_2030.csv")
str(dynamic)
batches_array <- split(batches, seq(nrow(batches)))

#### create cost landscape
library(terra)

land_path <- dynamic$landscape[1]
asc_out   <- dynamic$cost_map[1]

land_rast <- rast(land_path)

# Create cost raster (all cells = 1)
cost_rast <- setValues(land_rast, 1)

# Set NA flag (RangeShifter requires integer nodata)
NAflag(cost_rast) <- -9999

# ✔ Write clean ASCII in RangeShifter-compatible format
writeRaster(cost_rast, asc_out, filetype = "AAIGrid", overwrite = TRUE, NAflag = -9999)

# Run Simulations
#############

if (DEBUG) {
  # Run only the first simulation for debugging
  run_sims(batches_array[[1]], base_config)
} else {
  # Load parallel packages
  install.packages("furrr")
  library(furrr)
  library(future)

  # Set up parallel plan (adjust workers if needed)
  plan(multisession, workers = parallel::detectCores()/2)

  # Run all simulations in parallel
  future_map(batches_array, ~run_sims(.x, base_config))
}

