#############################################
# LOAD LIBRARIES
#############################################

library(Rcpp)
library(devtools)
library(Rdpack)
devtools::install_github("https://github.com/RangeShifter/RangeShiftR-pkg", 
                         ref = "main", subdir = "RangeShiftR")
library(parallel)
library(RangeShiftR)
library(yaml)

# Load custom functions
source("/Workspace/Users/fiona.plenderleith@defra.onmicrosoft.com/RS/code/ips-shift-functions.R")


#############################################
# USER PARAMETERS (SET HERE)
#############################################

site_val <- "FoD_updated"     # <--- CHANGE SITE HERE
year     <- 2020              # <--- CHANGE YEAR HERE

BASE_CONFIG_PATH  <- "/Workspace/Users/fiona.plenderleith@defra.onmicrosoft.com/RS/config-SMS.yml"
BATCHES_FILE_PATH <- sprintf("/dbfs/mnt/lab/unrestricted/Fiona/Batches_%s_%d_expanded.csv",site_val, year)

DEBUG <- FALSE
#############################################
# LOAD BASE CONFIG & INPUT TABLES
#############################################

base_config <- yaml::read_yaml(BASE_CONFIG_PATH)

# -------------------------------
# Load batches (already year + site dependent)
# -------------------------------
batches <- read.csv(BATCHES_FILE_PATH, stringsAsFactors = FALSE)

# -------------------------------
# Load original dynamic file (year‑dependent)
# -------------------------------
dynamic_in_path <- sprintf(
  "/dbfs/mnt/lab/unrestricted/Fiona/dynamic_landscapes_%d.csv",
  year
)
dynamic <- read.csv(dynamic_in_path, stringsAsFactors = FALSE)


#############################################
# PATH REPLACEMENT FUNCTION
#############################################

replace_data_prefix <- function(x, site) {
  if (!is.character(x)) return(x)

  stripped <- sub("^\\.\\./data(/processed)?/", "", x)
  new_base <- "/dbfs/mnt/lab/unrestricted/Fiona"

  # phenips* and dynamic_landscapes* go in Fiona/
  # everything else goes in Fiona/<site>/
  ifelse(
    grepl("^(phenips|dynamic_landscapes)", stripped, ignore.case = TRUE),
    file.path(new_base, stripped, fsep = "/"),
    file.path(new_base, site, stripped, fsep = "/")
  )
}


#############################################
# APPLY PATH FIXES TO BATCHES
#############################################

batch_path_cols <- intersect(
  c("phenips_file", "init_inds_file", "dynamic_landscapes"),
  names(batches)
)

for (col in batch_path_cols) {
  batches[[col]] <- replace_data_prefix(batches[[col]], site_val)
}

# -------------------------------
# Correct PHENIPS filename explicitly (your required pattern)
# -------------------------------
get_phenips_prefix <- function(site_val) sub("_updated$", "", site_val)

phenips_prefix <- get_phenips_prefix(site_val)

batches$phenips_file <- sprintf(
  "/dbfs/mnt/lab/unrestricted/Fiona/phenips%s%d.csv",
  phenips_prefix, year
)

# -------------------------------
# Point batches to UPDATED dynamic file
# -------------------------------
batches$dynamic_landscapes <- sprintf(
  "/dbfs/mnt/lab/unrestricted/Fiona/dynamic_landscapes_%d_UPDATE.csv",
  year
)


#############################################
# APPLY PATH FIXES TO DYNAMIC FILE
#############################################

char_cols <- vapply(dynamic, is.character, logical(1))
dynamic[char_cols] <- lapply(dynamic[char_cols], function(col)
  replace_data_prefix(col, site_val)
)


#############################################
# SAVE UPDATED DYNAMIC FILE
#############################################

dynamic_out_path <- sprintf(
  "/dbfs/mnt/lab/unrestricted/Fiona/dynamic_landscapes_%d_UPDATE.csv",
  year
)

write.csv(dynamic, dynamic_out_path, row.names = FALSE)


#############################################
# CREATE COST LANDSCAPE (terra)
#############################################

library(terra)

land_path <- dynamic$landscape[1]
asc_out   <- dynamic$cost_map[1]

land_rast <- rast(land_path)

# Cost raster = all 1s
cost_rast <- setValues(land_rast, 1)

# RangeShifter requires integer NA flag
NAflag(cost_rast) <- -9999

writeRaster(
  cost_rast,
  asc_out,
  filetype = "AAIGrid",
  overwrite = TRUE,
  NAflag = -9999
)


batches_array <- split(batches, seq(nrow(batches)))

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