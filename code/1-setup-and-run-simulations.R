library(uuid)
library(RangeShiftR)
library(readxl)
library(openxlsx)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ips-shift-functions.R")

# input_data_path is the folder where different simulation inputs are stored
# these are created by "ips-shift-data-prep.R", and provide habitat quality
# and initial distribution files for different simulation extents
input_data_path <- "../data/processed"

output_path <- "../outputs"


# A set of different habitat extents (small extents used mainly for testing).
# extent(xmin, xmax, ymin, ymax)
# Any extent can be set, but if the files do not exist, it will be created below.
# clipped_extent <- terra::ext(605000, 620000, 135000, 150000) # very small clipped extent
clipped_extent <- terra::ext(425000, 645000, 70000, 240000) # small clipped extent
# clipped_extent <- terra::ext(120000, 650000, 10000, 240000) # medium clipped extent

# Aggregation factor is the factor by which the original raster file pixel height/width is aggregated
# Note the original raster file pixel size is hard coded here as 25 metres (ideally
# this should be sourced from the original raster files if there is any chance
# it will change in the future)
aggregation_factor <- 8
pixel_size <- 25 * aggregation_factor

# A string of the simulation extent. Folders storing the relevant inputs will be named
# using extent_str.
extent_str <- sprintf(
  "xmin%s-xmax%s-ymin%s-ymax%s",
  clipped_extent$xmin,
  clipped_extent$xmax,
  clipped_extent$ymin,
  clipped_extent$ymax
)

# If the habitat files do not exist for this extent, create them.
# NOTE the .asc files created often have a decimal in the first value of the raster
# I'm not sure why this happens. It needs to be manually removed otherwise RangeShiftR errors
habitat_rst_name <- "hab_qual"
if (!file.exists(file.path(input_data_path, extent_str, sprintf("%s.asc", habitat_rst_name)))) {
  
  create_habitat_and_initial_dist(
    species_shp_path="../data/raw/species_distribution/IpsShift_points.shp",
    lcm_rst_path="../data/raw/LCM 2021/data/gblcm25m2021.tif",
    nfi_rst_path="../data/raw/England_NFI_2020/England_NFI_2020_with_woodland_type_code.tif",
    spruce_rst_path="../data/raw/Spruce_Map/Spruce_Map.tif",
    clipped_extent=clipped_extent,
    aggregation_factor=aggregation_factor,
    dirpath=input_data_path
    )
}

# The number of "rangeshifter steps"/years that should be run before it is "paused"
# and changes are made to the landscape.
# This is unlikely to change from 1.
rangeshiftr_steps <- 1

df_phenips <- read_excel(path="../data/Mockdata_ips_climate.xlsx", sheet="PHENIPS")

phenips <- list(df_phenips$Generations)

# Create arrays here of the parameters to be used, and a dataframe (df_params)
# is created with a full factorial design.
# Note the parameters currently listed are the only ones that have been
# configured for control in this batch mode. Any additional parameters will
# need to be set up.
rmax_list <- c(5) # c(5, 10, 20)
D0_list <- c(0.7)
alpha_list <- c(10)
beta_list <- c(0.5)
short_range_dist_list <- c(5000)
long_range_dist_list <- c(50000)
short_dist_prob_list <- c(0.995) # c(0.995, 0.99, 0.95)
carrying_capacity_list <- c(200) # c(100, 200, 300, 400, 1000, 5000)
replicates_list <- c(1)
years_list <- c(20)
# Additional parameters for interaction with the landscape
years_for_destruction_list <- c(20) # c(2, 5, 10, 20)
apply_destruction_annually_list <- c(FALSE) # c(FALSE, TRUE)
phenips_list <- c(phenips) # CANNOT USE MULTIPLE VERSIONS OF PHENIPS PER SIM BATCH
dynamic_hab_quality_landscapes <- c("hab_qual_1", "hab_qual_2", "hab_qual_3", "hab_qual_4", "hab_qual_5")
dynamic_hab_quality_years <- c(0, 4, 8, 12, 16)
delete_output_maps <- TRUE # If TRUE, deletes all output raster maps at the end of each simulation


## CREATE LIST OF LANDSCAPES TO USE FOR EACH YEAR
dynamic_hab_quality_landscape_list <- list()
i<-1
for (idx in 1:length(dynamic_hab_quality_landscapes)) {
  cat(sprintf("idx: %s, ", idx))
  if (length(dynamic_hab_quality_landscapes) < idx+1) {
    max_idx <- length(dynamic_hab_quality_landscapes)
  } else {
    max_idx <- idx+1
  }
  for (year in dynamic_hab_quality_years[idx]:(dynamic_hab_quality_years[max_idx]-1)) {
    cat(sprintf("year: %s, ", year))
    dynamic_hab_quality_landscape_list[[i]] <- dynamic_hab_quality_landscapes[idx]
    i <- i + 1
  }
}

# Create dataframe of all parameter combinations
df_params <- expand.grid(
  rmax = rmax_list,
  short_dist_p = short_dist_prob_list,
  D0 = D0_list,
  a = alpha_list,
  b = beta_list,
  short_range_d = short_range_dist_list,
  long_range_d = long_range_dist_list,
  carry_cap = carrying_capacity_list,
  replicates = replicates_list,
  years = years_list,
  years_for_destruction = years_for_destruction_list,
  apply_destruction_annually = apply_destruction_annually_list,
  phenips = phenips_list
)

# Create a unique ID for each simulation, this is used as the directory
# name for the simulation
df_params["sim_uuid"] = sapply(1:nrow(df_params), UUIDgenerate)
sim_start_time <- format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss")
sim_batch_path <- file.path(output_path, sim_start_time)
dir.create(sim_batch_path)
write.xlsx(df_params, file.path(sim_batch_path, "sim-parameters.xlsx"))


##################
# Run simulations
##################
for(i in 1:nrow(df_params)) {
  row <- df_params[i,]
  
  rmax <- row$rmax
  short_dist_prob <- row$short_dist_p
  D0 <- row$D0
  alpha <- row$a
  beta <- row$b
  short_range_dist <- row$short_range_d
  long_range_dist <- row$long_range_d
  carrying_capacity <- row$carry_cap
  sim_uuid <- row$sim_uuid
  replicates <- row$replicates
  years <- row$years
  years_for_destruction <- row$years_for_destruction
  apply_destruction_annually <- row$apply_destruction_annually
  phenips <- row$phenips
  
  dirpath <- file.path(sim_batch_path, sim_uuid)
  
  if (file.exists(dirpath)) { next }
  
  dir.create(dirpath, showWarnings = TRUE)
  dir.create(file.path(dirpath, "Inputs"), showWarnings = TRUE)
  dir.create(file.path(dirpath, "Outputs"), showWarnings = TRUE)
  dir.create(file.path(dirpath, "Outputs", "file-per-pause"), showWarnings = TRUE)
  dir.create(file.path(dirpath, "Output_Maps"), showWarnings = TRUE)
  
  for (habitat_rst_name in dynamic_hab_quality_landscapes) {
    file.copy(file.path(input_data_path, extent_str, sprintf("%s.asc", habitat_rst_name)), file.path(dirpath, "Inputs"))
    file.copy(file.path(input_data_path, extent_str, sprintf("%s.tif", habitat_rst_name)), file.path(dirpath, "Inputs"))
  }
  file.copy(file.path(input_data_path, extent_str, "inds.txt"), file.path(dirpath, "Inputs"))
  
  
  land <- ImportedLandscape(
    LandscapeFile = dynamic_hab_quality_landscapes[1],
    Resolution = pixel_size,
    HabPercent=TRUE, # Landscape file represents habitat quality as percentage
    K_or_DensDep = carrying_capacity
  )
  
  demo <- Demography(
    ReproductionType=0, # asexual / only female model
    Rmax = rmax
  )
  
  disp <-  Dispersal(
    Emigration = Emigration(DensDep=TRUE, EmigProb=matrix(c(D0, alpha, beta), nrow=1)), 
    Transfer = DispersalKernel(
      DoubleKernel = TRUE,
      Distances = matrix(c(short_range_dist, long_range_dist, short_dist_prob), ncol = 3)
    ),
    Settlement = Settlement()
  )
  
  init <- Initialise(
    InitType=2, # Initialise from individuals list
    InitIndsFile="inds.txt"
  )
  
  sim <- Simulation(
    Simulation = 0,
    Years = rangeshiftr_steps,
    Replicates = 1,
    OutIntPop = 1,
    OutIntRange = 1
  )
  
  sim_params <- RSsim(batchnum = 0, land = land, demog = demo, dispersal = disp, simul = sim, init = init)  # seed to allow replication
  
  # Set up the extra parameters used as part of RangeShiftR pausing
  attributes(sim_params)$pausingextras <- PausingExtras(
    replicates=replicates,
    years=years,
    pop_impacts_quality=T,
    years_for_destruction=years_for_destruction,
    apply_destruction_annually=apply_destruction_annually,
    phenips=phenips,
    dynamic_hab_quality_landscape_list=dynamic_hab_quality_landscape_list
  )

  validateRSparams(sim_params)


  run_pausing_simulations(
    sim_params = sim_params,
    input_data_path = input_data_path,
    dirpath = dirpath,
    delete_output_maps = delete_output_maps
  )


}

