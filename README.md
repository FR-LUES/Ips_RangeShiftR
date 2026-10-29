# Ips_RangeShiftR

This repository contains the R scripts, configuration files, and analysis notebooks used to build habitat inputs, run RangeShiftR models for *Ips typographus*, and analyse simulation outputs. The workflow is organised in three main stages: habitat model and inputs, RangeShiftR modelling including the DASH application, and post-processing of outputs.

## Repository structure

* `code/` — core scripts for habitat preparation, simulation setup, execution, and helper functions.
* `data/` — tracked example input files used by the modelling workflow. Most other data files are ignored by `.gitignore`.
* `IpsShift_DASH/` — DASH-specific run scripts and helper functions for site- and time-slice-specific model runs.
* `Data_analysis/` — notebooks and scripts for analysing and visualising model outputs.
* `tutorials/` — earlier tutorial material and the R project file.
* `outputs/` — default output location for simulations; ignored by Git.
* `config-SMS.yml` and `config-dispersal-kernel.yml` — configuration files for the main RangeShiftR workflows.

## Habitat model and inputs

This stage creates the habitat suitability and other spatial or phenology inputs required by the simulations.

### Main scripts

* `code/Landscape_generator.R` — creates yearly habitat suitability maps from spruce, wetness, wind-damage, felling, and drought inputs; these maps are used as landscape inputs to RangeShiftR.
* `code/Map_starting.R` — prepares starting maps for simulations.
* `code/SPEI.R` — extracts site-level SPEI summaries and converts drought values into habitat-quality predictions.
* `code/code_extractphenips.R` — reads PHENIPS Excel grids, crops them to study sites, and summarises values by year.
* `code/Init_sites_PAM.R` — generates candidate starting locations from habitat and quality rasters and summarises surrounding habitat composition.

### Example tracked inputs in `data/`

* `Batches.csv` — example batch-definition file; each row defines one simulation and can override selected parameters from the base configuration.
* `dynamic_landscapes_1.csv` — example schedule of yearly habitat and cost-map files used by the pausing workflow.
* `phenips1.csv` — example PHENIPS generations file mapping calendar years to model years and generation counts.

### Other main inputs used by the habitat model

* Spruce distribution layers — used to define host-tree availability in the landscape-generation workflow.
* TWI rasters — used to identify wetness or drought-risk conditions in landscape preparation.
* Subcompartment and felling layers — used to represent windblow and felled areas in habitat maps.
* ForestGALES outputs — used to incorporate wind-damage predictions into habitat suitability maps.
* SPEI summaries — used to represent drought stress through time.
* PHENIPS Excel grids — used to derive year-specific phenology inputs for simulations.

Additional input files can be stored under `data/`, but most are excluded from version control.

## RangeShiftR models and DASH application

This stage runs the main RangeShiftR simulations and the DASH-specific applications of the model.

### Standard RangeShiftR workflow

* `code/1-setup-and-run-simulations.R` — main entry point for batch simulation runs; loads the base configuration and batch definitions, and supports a `DEBUG` mode for running a single simulation locally.
* `code/ips-shift-functions.R` — core helper functions for the pausing workflow, including dynamic landscapes and habitat-destruction logic.
* `config-SMS.yml` — base simulation configuration for the standard workflow.
* `config-dispersal-kernel.yml` — alternative configuration file for dispersal-kernel settings.

### DASH application

* `code/Run_sims_DASH.R` — runs a DASH-configured simulation workflow for a selected site and year, updating paths and output locations.
* `IpsShift_DASH/ips-shift-functions_DASH.R` — helper functions for DASH runs, including population extraction and habitat-quality updates between pauses.
* `IpsShift_DASH/RunRS_FoD2020.R`, `IpsShift_DASH/RunRS_FoD2050.R`, `IpsShift_DASH/RunRS_FoD2070.R` — Forest of Dean DASH runs for different time slices.
* `IpsShift_DASH/RunRS_Kielder2020.R`, `IpsShift_DASH/RunRS_Kielder2050.R`, `IpsShift_DASH/RunRS_Kielder2070.R` — Kielder DASH runs for different time slices.
* `IpsShift_DASH/RunRS_NYM2020.R`, `IpsShift_DASH/RunRS_NYM2050.R`, `IpsShift_DASH/RunRS_NYM2070.R` — North York Moors DASH runs for different time slices.

## Analysis of outputs

This stage summarises, visualises, and interprets the outputs produced by the simulations.

### Main analysis files

* `code/2-results-analysis.R` — initial analysis script for simulation outputs.
* `Data_analysis/occupancy_dynamics_habitat_suitability` — notebook for mapping occupancy through time and overlaying occupancy on habitat-suitability rasters.
* `Data_analysis/time_series_plots` — notebook for plotting occupancy, spread, and extinction summaries through time.
* `Data_analysis/plot_starting_map.R` — crops a suitability raster and plots the chosen starting locations over it.
* `Data_analysis/initialise_sites_pam.R` — analysis copy of the site-initialisation workflow for generating candidate starting locations.
* `Data_analysis/extract_phenips_data.R` — analysis copy of the PHENIPS extraction workflow for site and year summaries.

## Dependencies

This project requires R plus the packages used by the simulation scripts, including:

* `RangeShiftR`
* `yaml`
* `parallel`
* `sf`
* `terra`

`RangeShiftR` needs to be installed from source. See the [RangeShiftR installation guide](https://rangeshifter.github.io/RangeshiftR-tutorials/installing.html).

## Running simulations

The main workflow is driven from `code/1-setup-and-run-simulations.R`:

1. Update the working paths in the script if needed for your environment.
2. Set the base configuration file.
3. Set the batch definition file.
4. Use `DEBUG <- TRUE` to run a single simulation in the current R session, or `DEBUG <- FALSE` to run the full batch in parallel.
5. Review outputs in the default `outputs/` location.

## Notes

* The `data/` directory is mostly ignored by Git; only selected example files are tracked.
* `Data_analysis/` and `IpsShift_DASH/` contain additional analysis and alternative run workflows beyond the main `code/` pipeline.
* The `outputs/` directory is not stored in the repository.
* The `tutorials/` directory contains earlier tutorial-based material.
