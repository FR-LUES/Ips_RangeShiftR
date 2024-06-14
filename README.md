# Ips_RangeShiftR

This repository contains R code to run RangeShiftR within a "pausing" framework that allows for interactions between the landscape and the modelled species.

# Structure

This repository contains four main directories:
- `code` : contains the code to setup, run and analyse simulations
- `data` : not stored in the repository, but needed for data storage. It contains the `raw` and `processed` raster and vector data to be used for the simulations (this data is not stored in the repository due to size)
- `outputs` : not stored in the repository, but it is the default location for simulation outputs to be saved.
- `tutorial` : contains the earliest files added to this repository, based on [rangeshifter.github.io](rangeshifter.github.io) tutorials

# Dependencies

The `RangeShiftR` R package needs to be installed. It cannot be installed like a standard package, so it needs to be installed from source as explained in the [RangeShiftR GitHub page](https://rangeshifter.github.io/RangeshiftR-tutorials/installing.html).

# R Scripts

There are three R scripts in the `code` directory:
- `1-setup-and-run-simulations.R`
- `2-results-analysis.R`
- `ips-shift-functions.R`

## `1-setup-and-run-simulations.R`

Contains the code to set up and run batches of RangeShiftR simulations.

This includes the creation of the landscape and species files required to initialise simulations. The user can set the extent and provide a species distribution shapefile, and this script will create the required initialisation files.

This file also allows you to create a list of parameter values to experiment with common RangeShifter parameters that we have used so far (e.g. Rmax, D0, alpha, beta, dispersal distance, ...).

### Existing RangeShiftR parameters that can currently be controlled:
- `rmax`
- `D0`
- `alpha`
- `beta`
- `short_range_dist`
- `long_range_dist`
- `short_dist_prob`
- `carrying_capacity`
- `replicates`
- `years`

### Additional parameters for interaction with the landsacpe
- `years_for_destruction`: how many years of maximum occupancy does it take for the species to destroy the habitat
- `apply_destruction_annually: does habitat degradation happen each year, or only after `year_for_destruction` years
- `phenips`: rmax multiplier - effecively broods per year, as influence by climate

The file will generate a table displaying the parameters being tested, and will iterate through and run simulations for all parameter combinations.

## `2-results-analysis.R`

Contains the code to perform initial analysis of simulation results.

This file is an early work in progress, and we may want to add a number of other analyses and plots to investigate simulation outputs further.


## `ips-shift-functions.R`

Contains a range of functions required to move between different data formats required by RangeShiftR, and to run the simulations.

Functions in this file should not generally require any manipulation until we add extra functionality to our simulation set up. At the moment simulations allow species presence to impact on habitat quality, but temporal changes to habitat quality have not been implemented (i.e. habitat changes due to climate, or wind damage).