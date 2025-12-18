# Ips_RangeShiftR

This repository contains R code to run RangeShiftR within a "pausing" framework that allows for interactions between the landscape and the modelled species.

# Structure

This repository contains four main directories:
- `code` : contains the code to setup, run and analyse simulations
- `data` : can be used to store the input data to be used for simulations. This should be kept empty in the repository, apart from example files (e.g. the batch csv file used to specify parameter variations across a batch).
- `outputs` : not stored in the repository, but it is the default location for simulation outputs to be saved.
- `tutorial` : contains the earliest files added to this repository, based on [rangeshifter.github.io](rangeshifter.github.io) tutorials

# Dependencies

The `RangeShiftR` R package needs to be installed. It cannot be installed like a standard package, so it needs to be installed from source as explained in the [RangeShiftR GitHub page](https://rangeshifter.github.io/RangeshiftR-tutorials/installing.html).

# R Scripts

There are five R scripts in the `code` directory:
- `1-setup-and-run-simulations.R`
- `2-results-analysis.R`
- `ips-shift-functions.R`
- `Landscape_generator.R`
- `Map_starting.R`

## `1-setup-and-run-simulations.R`

In this file, you should specify the base config and batch file, and it will run the specified batches of RangeShiftR simulations. You can also use the `DEBUG` parameter which runs only a single simulation in the current R session (rather than on a cluster), making it easier to debug errors.

The base config (`BASE_CONFIG_PATH <- "../config-SMS.yml"`) should specify all simulation parameters. The batch file (`BATCHES_FILE_PATH <- "../data/batches.csv"`) should specify a separate simulation in each row (named based on the `batch_name` column), any parameters specified in the batch file will override the base config parameters.

### RangeShiftR parameters that can currently be controlled:
```yaml
# Demography
rmax: 2
# Dispersal
## Emigration
D0: 0.7
alpha: 10
beta: 0.5
## Transfer
directional_persistence: 2
perceptual_range: 2
mem_size: 1
goal_type: 2
goal_bias: 2.5
alpha_density_bias: 0.4
beta_density_bias: 10
step_mortality: 0.1 
## Settlement
max_settlement_prob: 0.4
density_inflection: 10
density_slope: 0.05
min_steps: 1
max_steps: 100
# Initialise
init_inds_file: ../data/processed/inds_init.txt
# Land
carrying_capacity: 1000
```

### Additional "pausing" parameters
```yaml
# Does the habitat destruction impact quality every year, or only after `years_for_destruction`
apply_destruction_annually: TRUE

# path to a csv file specifying the habitat quality rasters and cost map rasters to be used for each year of the simulation
dynamic_landscapes: ../data/dynamic_landscapes_1.csv

# specify whether to keep the individual files generated each time Rangeshifter is run
keep_files_per_pause: FALSE

# specify whether to keep the input files used for the simulation
keep_input_files: FALSE

# the base output path within which the simulations will be saved
output_path: ../outputs

# The interval at which to save the individuals file (0 means never, 1 means every year)
out_int_ind: 0

# The interval at which to save the in SMS path file (0 means never, 1 means every year)
out_int_paths: 0

# The interval at which to save the species range file (0 means never, 1 means every year)
out_int_range: 5

# The phenips file specifying the number of generations for each year of the simulation
phenips_file: ../data/phenips1.csv

# Does the population impact the quality of the landscape?
pop_impacts_quality: TRUE

# How many replicates to run for this simulation
replicates: 5

# Number of years for the simulation
years: 20

# Number of years it would take for a fully occupied cell of 100% quality to be fully degraded
years_for_destruction: 10
```

## `2-results-analysis.R`

Contains the code to perform initial analysis of simulation results.

This file may now be redundant since major changes to how simulations are run.


## `ips-shift-functions.R`

Contains all the functions required to run the pausing RangeShifter simulations - loading dynamic landscapes and applying habitat destruction by the species.
