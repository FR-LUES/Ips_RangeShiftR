# Ips_RangeShiftR

This repository contains R code to run RangeShiftR in within a "pausing" framework that allows for interactions between the landscape and the modelled species.

# Structure

This repository contains four main directories:
- `code` : contains the code to setup, run and analyse simulations
- `data` : contains the `raw` and `processed` raster and vector data to be used for the simulations (this data is not stored in the repository due to size)
- `outputs` : default location for simulation outputs to be saved - this folder will be remain empty in the online repository
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

This file also allows you to create a list of parameter values to experiment with for common RangeShifter parameters that we have used so far (e.g. Rmax, D0, alpha, beta, dispersal distance, ...).

The file will generate a table displaying the parameters being tested, and will iterate through and run simulations for all parameter combinations.

## `2-results-analysis.R`

Contains the code to perform initial analysis of simulation results.

This file is an early work in progress, and we may want to add a number of other analyses and plots to investigate simulation outputs further.


## `ips-shift-functions.R`

Contains a range of functions required to move between different data formats required by RangeShiftR, and to run the simulations.

Functions in this file should not generally require any manipulation until we add extra functionality to our simulation set up. At the moment simulations allow species presence to impact on habitat quality, but temporal changes to habitat quality have not been implemented (i.e. habitat changes due to climate, or wind damage).