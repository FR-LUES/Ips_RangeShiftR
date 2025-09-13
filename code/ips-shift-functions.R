library(terra)
library(sf)

setwd(getSrcDirectory(function(){})[1])

# Print numbers to strings without exponential notion
options(scipen = 999)


create_habitat_and_initial_dist <- function(
    species_shp_path,
    lcm_rst_path,
    nfi_rst_path,
    spruce_rst_path,
    clipped_extent,
    aggregation_factor,
    dirpath
) {
  
  extent_str <- sprintf("xmin%s-xmax%s-ymin%s-ymax%s", clipped_extent$xmin, clipped_extent$xmax, clipped_extent$ymin, clipped_extent$ymax)
  dir.create(file.path(dirpath, extent_str), recursive=T)
  
  shp_species <- read_sf(dsn = species_shp_path)
  
  rst_lcm <- terra::rast(lcm_rst_path)
  rst_lcm <- rst_lcm[[1]]
  names(rst_lcm) <- c("lcm")
  rst_nfi <-  terra::rast(nfi_rst_path)
  rst_spruce <-  terra::rast(spruce_rst_path)
  
  # all input rasters should be the same resolution
  stopifnot(res(rst_lcm) == res(rst_nfi), res(rst_lcm) == res(rst_spruce))
  original_res <- res(rst_lcm)[1]
  new_res <- original_res * aggregation_factor
  
  # Merge the habitat maps into single ascii raster for rangeshifter
  rst_lcm <- crop(rst_lcm, clipped_extent)
  rst_nfi <- crop(rst_nfi, clipped_extent)
  rst_spruce <- crop(rst_spruce, clipped_extent)
  
  if (aggregation_factor > 1) {
    rst_lcm <- terra::aggregate(rst_lcm, fact=aggregation_factor, fun="modal", na.rm=T)
    rst_nfi <- terra::aggregate(rst_nfi, fact=aggregation_factor, fun="modal", na.rm=T)
    rst_spruce_hab <- terra::aggregate(rst_spruce, fact=aggregation_factor, fun="max", na.rm=T)
    # take sum here because spruce=1, so the sum is the number of 25m pixels occurring in
    # the aggregated pixel
    rst_hab_quality <- terra::aggregate(rst_spruce, fact=aggregation_factor, fun="sum", na.rm=T)
  } else {
    rst_spruce_hab <- rst_spruce
    rst_hab_quality <- rst_spruce
  }
  
  rst_lcm[rst_nfi == 5]<- 22 # conifer
  rst_lcm[rst_nfi == 14] <- 23 # mixed mainly conifer
  rst_lcm[rst_spruce_hab == 1] <- 24 # spruce
  
  # rst_lcm <- terra::classify(rst_lcm, cbind(NA, -9999)) 
  terra::writeRaster(rst_lcm, file.path(dirpath, extent_str, "hab_code.tif"), overwrite=TRUE)
  
  # convert quantity of spruce pixels to a percentage per aggregated pixel
  rst_hab_quality <- (rst_hab_quality / aggregation_factor**2) * 100
  rst_hab_quality <- round(rst_hab_quality)
  rst_hab_quality[rst_lcm < 24] <- 0
  terra::writeRaster(rst_hab_quality, file.path(dirpath, extent_str, "hab_quality.tif"), datatype="INT2S", overwrite=TRUE)
  
  # Convert species distribution point to ascii raster for rangeshifter
  # rst_species <- terra::classify(rst_lcm, cbind(0, 25, NA))
  rst_species <- terra::classify(rst_lcm, cbind(0, 25, 0))
  rst_species <- terra::rasterize(shp_species, rst_species, field="pop")
  names(rst_species) <- c("pop")
  terra::writeRaster(rst_species, file.path(dirpath, extent_str, "initial_dist.tif"), overwrite=TRUE)
  
  # Write a individuals file for initialisation
  pop_size <- 1000
  df_individuals <- create_inds_dataframe(shp_species, rst_species)
  df_individuals$Ninds <- df_individuals$Ninds * pop_size
  df_individuals <- df_individuals[, c("Year", "Species", "X", "Y", "Ninds")]
  write.table(df_individuals,  file.path(dirpath, extent_str, "inds.txt"), row.names = F, quote = F, sep = '\t')
  
  extent_lcm <- terra::ext(rst_lcm)
  # terra::ext(rst_lcm) <- c(0, terra::ext(rst_lcm)[2]-terra::ext(rst_lcm)[1], 0, terra::ext(rst_lcm)[4]-terra::ext(rst_lcm)[3])
  # terra::ext(rst_species) <- c(0, terra::ext(rst_species)[2]-terra::ext(rst_species)[1], 0, terra::ext(rst_species)[4]-terra::ext(rst_species)[3])
  
  terra::writeRaster(
    rst_lcm,
    file.path(dirpath, extent_str, "hab_code.asc"),
    datatype="INT2S",
    NAflag=-9999
  )
  terra::writeRaster(
    rst_hab_quality,
    file.path(dirpath, extent_str, "hab_quality.asc"),
    datatype="INT2S",
    NAflag=-9999
  )
  terra::writeRaster(
    rst_species,
    file.path(dirpath, extent_str, "initial_dist.asc"),
    datatype="INT2S",
    NAflag=-9999
  )
}


create_inds_dataframe <- function(
    shp_species,
    rst_species
) {
  
  
  df_inds <- terra::extract(rst_species, shp_species, cells=T, df=TRUE)
  
  df_inds$Year <- 0
  df_inds$Species <- 0
  df_inds$X <- df_inds$cell %% ncol(rst_species) - 1
  df_inds$Y <- nrow(rst_species) - (floor(df_inds$cell / ncol(rst_species))) - 1
  df_inds$Ninds <- df_inds$pop
  df_inds <- df_inds[ , !(names(df_inds) %in% c('ID', 'cells', 'pop'))]
  df_inds <- df_inds[!is.na(df_inds$Ninds),]
  
  # drop duplicates
  df_inds <- unique(df_inds)
  # add another catch where Ninds == 0 (remove)
  if (nrow(df_inds[which(df_inds$Ninds==0),])>0){
    df_inds <- df_inds[-which(df_inds$Ninds==0),]
  }
  
  return(df_inds)
  
}



extract_pop_raster_to_shp <- function(
    pop_raster
) {
  
  shp_file <- rasterToPoints(pop_raster, fun=function(x){x > 0}, spatial=TRUE) %>% st_as_sf()
  shp_file <- shp_file %>% st_set_crs(st_crs(pop_raster))
  shp_file <- st_transform(shp_file, crs(pop_raster))
  shp_file$id <- 1:nrow(shp_file)
  shp_file$layer <- ceiling(shp_file$layer)
  
  return(shp_file)
}




extract_pop_file_to_shp <- function(
    sim_params,
    coord_ref_sys,
    geom_extent,
    dirpath
) {
  
  df_pop <- read.table(file.path(dirpath, "Outputs", "Batch0_Sim0_Land1_Pop.txt"), sep="\t", header=T)
  df_pop <- subset(df_pop, Year==sim_params@simul@Years & NInd > 0)
  df_pop$x_coord <- (df_pop$x + 0.5) * sim_params@land@Resolution + geom_extent$xmin
  df_pop$y_coord <- (df_pop$y + 0.5) * sim_params@land@Resolution + geom_extent$ymin
  if (nrow(df_pop) > 0) {
    df_pop <- st_as_sf(df_pop, coords = c("x_coord","y_coord"))
    df_pop <- df_pop %>% st_set_crs(st_crs(coord_ref_sys))
  }
  
  return(df_pop)
  
}




calculate_pop_impact_on_hab_quality <- function(
    pop_rst,
    years_for_destruction,
    apply_destruction_annually,
    sim_params
) {
  
  # Calculate carrying capacity per pixel
  # Rangeshifter stores it as per hectare
  carry_cap_per_pixel <- sim_params@land@K_or_DensDep * (sim_params@land@Resolution^2 / 10000)
  
  if (apply_destruction_annually) {
    quality_impact <- ((pop_rst / carry_cap_per_pixel) * 100) / years_for_destruction
  } else {
    quality_impact <- ((pop_rst / carry_cap_per_pixel) * 100)
  }
  
  quality_impact <- round(quality_impact)
  
  return(quality_impact)
  
}


apply_pop_impact_on_hab_quality <- function(
    quality_impact,
    landscape,
    sim_params,
    coord_ref_sys,
    geom_extent,
    dirpath
    
) {
  
  initial_hab_qual_rst  <- terra::rast(file.path(dirpath, "Inputs", landscape))
  crs(initial_hab_qual_rst) <- coord_ref_sys
  terra::ext(initial_hab_qual_rst) <- geom_extent
  
  hab_qual_rst <- round(initial_hab_qual_rst - quality_impact)
  hab_qual_rst[hab_qual_rst < 0] <- 0
  
  return(hab_qual_rst)
  
}



extract_pop_shp_to_table <- function(
    pop_shp
) {
  
  if ("Nfemales" %in% colnames(pop_shp)) {
    pop_shp <- pop_shp[, c("Year", "Species", "x", "y", "NInd", "Nfemales", "Nmales")]
    # Create copy for males and females, only including rows with a population
    pop_shp_female <- subset(pop_shp, Nfemales > 0)
    pop_shp_male <- subset(pop_shp, Nmales > 0)
    
    # The oututs have "Nfemales" and "Nmales", but the input Inds file needs to have
    # "NInd" and "Sex" columns
    pop_shp_female$NInd <- pop_shp_female$Nfemales
    pop_shp_female <- pop_shp_female[, c("Year", "Species", "x", "y", "NInd")]
    if (nrow(pop_shp_female) > 0) {
      pop_shp_female$Sex <- 0
    }
    
    pop_shp_male$NInd <- pop_shp_male$Nmales
    pop_shp_male <- pop_shp_male[, c("Year", "Species", "x", "y", "NInd")]
    if (nrow(pop_shp_male) > 0) {
      pop_shp_male$Sex <- 1
    }
    
    pop_shp <- rbind(pop_shp_female, pop_shp_male)
    
    if (nrow(pop_shp) > 0) {
      pop_shp <- pop_shp[order(pop_shp$x, pop_shp$y, pop_shp$Sex),]
    }
    
  } else {
    pop_shp <- pop_shp[, c("Year", "Species", "x", "y", "NInd")]
  }
  
  names(pop_shp)[names(pop_shp) == "x"] <- "X"
  names(pop_shp)[names(pop_shp) == "y"] <- "Y"
  names(pop_shp)[names(pop_shp) == "NInd"] <- "Ninds"
  
  if (nrow(pop_shp) > 0) {
    pop_shp$Year <- 0
  }
  
  pop_shp <- pop_shp %>% st_drop_geometry()
  
  return(pop_shp)
}




extract_pop_shp_to_rst <- function(
    pop_shp,
    rst_lcm
) {
  
  pop_rst <- terra::classify(rst_lcm, cbind(minmax(rst_lcm)[1], minmax(rst_lcm)[2], 0))
  if (sum(pop_shp$NInd) > 0) {
    pop_rst <- terra::rasterize(pop_shp, pop_rst, field="NInd")
  }
  pop_rst <- classify(pop_rst, cbind(NA, 0))
  
  return(pop_rst)
}




PausingExtras <- function(
    replicates,
    years,
    pop_impacts_quality,
    years_for_destruction,
    apply_destruction_annually,
    phenips,
    dynamic_landscapes_list,
    dynamic_cost_maps_list
) {
  
  setClass(
    "pausingextras",
    representation(
      replicates="numeric",
      years="numeric",
      pop_impacts_quality="logical",
      years_for_destruction="numeric",
      apply_destruction_annually="logical",
      phenips="list",
      dynamic_landscapes_list="list",
      dynamic_cost_maps_list="list"
    )
  )
  
  pausingextras <- new(
    "pausingextras",
    replicates=replicates,
    years=years,
    pop_impacts_quality=pop_impacts_quality,
    years_for_destruction=years_for_destruction,
    apply_destruction_annually=apply_destruction_annually,
    phenips=phenips,
    dynamic_landscapes_list=dynamic_landscapes_list,
    dynamic_cost_maps_list=dynamic_cost_maps_list
  )
  
  return(pausingextras)
  
}



run_pausing_simulations <- function(
    sim_params,
    input_data_path,
    dirpath
) {
  
  save(sim_params, file=file.path(dirpath, "Inputs", "RS-sim-param-object"))
  
  replicates <- sim_params@pausingextras@replicates
  years <- sim_params@pausingextras@years
  carrying_capacity <- sim_params@land@K_or_DensDep
  phenips_list <- sim_params@pausingextras@phenips
  base_rmax <- sim_params@demog@Rmax
  dynamic_landscapes_list <- sim_params@pausingextras@dynamic_landscapes_list
  dynamic_cost_maps_list <- sim_params@pausingextras@dynamic_cost_maps_list
  
  rst_lcm <- terra::rast(file.path(
    dirpath,
    "Inputs",
    sprintf("%s.tif", tools::file_path_sans_ext(basename(dynamic_landscapes_list[[1]])))
  ))
  
  start_init <- sim_params@init
  start_land <- sim_params@land
  
  for (replicate in 1:replicates) {
    
    # Reset initialisation and land parameters ready to start next replicate
    sim_params@init <- start_init
    sim_params@land <- start_land
    sim_params@land@LandscapeFile <- sprintf("%s.asc", tools::file_path_sans_ext(basename(dynamic_landscapes_list[[1]])))
    sim_params@land@CostsFile <- sprintf("%s.asc", tools::file_path_sans_ext(basename(dynamic_cost_maps_list[[1]])))
    total_impact_hab_quality <- rst_lcm
    total_impact_hab_quality[total_impact_hab_quality>=0] <- 0
    
    # Empty raster list to store a stack of population rasters
    out_rasters <- c()
    
    for (year in 1:years) {
      
      # If the year is less than or equal to the length of the landscape list,
      # then use the correct landscape in sequence.
      # Else use the last landscape.
      if (year+1 <= length(dynamic_landscapes_list)) {
        next_base_landscape <- sprintf("%s.asc", tools::file_path_sans_ext(basename(dynamic_landscapes_list[[year+1]])))
      } else {
        next_base_landscape <- sprintf("%s.asc", tools::file_path_sans_ext(basename(dynamic_landscapes_list[[length(dynamic_landscapes_list)]])))
      }
      if (year+1 <= length(dynamic_cost_maps_list)) {
        next_cost_map <- sprintf("%s.asc", tools::file_path_sans_ext(basename(dynamic_cost_maps_list[[year+1]])))
      } else {
        next_cost_map <- sprintf("%s.asc", tools::file_path_sans_ext(basename(dynamic_cost_maps_list[[length(dynamic_cost_maps_list)]])))
      }
      
      # Load the next phenips value
      phenips <- phenips_list[[year]]
      if (is.na(phenips)) {
        stop(sprintf("No phenips value for year %s", year))
      }
      if (!is.numeric(phenips)) {
        stop(sprintf("phenips value must be numeric"))
      }
      # Apply phenips multiplier to the base Rmax
      sim_params@demog@Rmax <- phenips * base_rmax
      OUTINTRANGE <- 3
      if ((year == 1) |
        (year %% OUTINTRANGE == 0)
      ) {
        sim_params@simul@OutIntRange <- 1
      } else {
        sim_params@simul@OutIntRange <- 0
      }
      
      
      # RangeShiftR needs a backslash at end of path
      result <- RunRS(sim_params, sprintf("%s/", dirpath))
      
      pop_shp <- extract_pop_file_to_shp(
        sim_params,
        crs(rst_lcm),
        terra::ext(rst_lcm),
        dirpath
      )
      st_write(pop_shp, file.path(dirpath, "Output_Maps", sprintf("rep%s_year%s_pop.geojson", replicate,  year)))
      
      pop_rst <- extract_pop_shp_to_rst(
        pop_shp,
        rst_lcm
      )
      terra::writeRaster(pop_rst, file.path(dirpath, "Output_Maps", sprintf("rep%s_year%s_pop.tif", replicate,  year)))
      out_rasters <- c(out_rasters, pop_rst)
      
      # Create new inidividuals initialisation file to reinitialise rangeshiftR
      # with same distribution of individuals as this year
      df_pop_table <- extract_pop_shp_to_table(pop_shp)
      
      new_init_inds_txt <- sprintf('rep%s_year%s_inds.txt', replicate, year)
      write.table(df_pop_table, file.path(dirpath, "Inputs", new_init_inds_txt), row.names = F, quote = F, sep = '\t')
      sim_params@init@InitIndsFile <- new_init_inds_txt
      
      
      # Create new adjusted habitat quality raster
      if ((sim_params@pausingextras@pop_impacts_quality) &
          ((sim_params@pausingextras@apply_destruction_annually) |
           (year %% sim_params@pausingextras@years_for_destruction == 0))
      ) {
        
        new_quality_impact_rst <- calculate_pop_impact_on_hab_quality(
          pop_rst=pop_rst,
          years_for_destruction=sim_params@pausingextras@years_for_destruction,
          apply_destruction_annually=sim_params@pausingextras@apply_destruction_annually,
          sim_params=sim_params
        )
        
        terra::writeRaster(new_quality_impact_rst, file.path(dirpath, "Output_maps", sprintf("rep%s_year%s_new_quality_impact.tif", replicate, year)))
        
        total_impact_hab_quality <- total_impact_hab_quality + new_quality_impact_rst
        total_impact_hab_quality[total_impact_hab_quality>100] <- 100
        
        terra::writeRaster(total_impact_hab_quality, file.path(dirpath, "Output_maps", sprintf("rep%s_year%s_total_quality_impact.tif", replicate, year)))
        
        new_hab_qual_rst <- apply_pop_impact_on_hab_quality(
          quality_impact=total_impact_hab_quality,
          landscape=next_base_landscape,
          sim_params=sim_params,
          coord_ref_sys=crs(rst_lcm),
          geom_extent=terra::ext(rst_lcm),
          dirpath=dirpath
        )
        
        if (year+1 <= length(dynamic_landscapes_list)) {
          next_landscape_name <- sprintf("rep%s_year%s_%s", replicate, year,  tools::file_path_sans_ext(basename(dynamic_landscapes_list[[year+1]])))
        } else {
          next_landscape_name <- sprintf("rep%s_year%s_%s", replicate, year,  tools::file_path_sans_ext(basename(dynamic_landscapes_list[[length(dynamic_landscapes_list)]])))
        }
        
        terra::writeRaster(new_hab_qual_rst, file.path(dirpath, "Output_maps", sprintf("%s.tif", next_landscape_name)))
        
        next_landscape <- sprintf("%s.asc", next_landscape_name)
        
        terra::writeRaster(
          new_hab_qual_rst,
          file.path(dirpath, "Inputs", next_landscape),
          datatype="INT2S",
          NAflag=-9999
        )
        
      } else {
        next_landscape <- next_base_landscape
      }
      
      sim_params@land@LandscapeFile <- next_landscape
      sim_params@land@CostsFile <- next_cost_map
      
      for (out_file in c("Pop", "Range", "Parameters")) {
        existing_file_path <- file.path(dirpath, "Outputs", sprintf("Batch0_Sim0_Land1_%s.txt", out_file))
        if (file.exists(existing_file_path)) {
          file.rename(
            existing_file_path,
            file.path(dirpath, "Outputs", "file-per-pause", sprintf("rep%s_year%s_Batch0_Sim0_Land1_%s.txt", replicate, year, out_file))
          )
        }
      }
      
      # If Population is zero, then break out of this simulation
      if (nrow(pop_shp) == 0) {
        cat("########### POPULATION OF ZERO -> SIMULATION ENDED EARLY ##########\n")
        break
      }
      
    }
    out_rst_stack <- terra::rast(out_rasters)
    terra::writeRaster(out_rst_stack, file.path(dirpath, "Output_maps", sprintf("rep%s_pop_stack.tif", replicate)))
    
  }
  
  # Merge all the pop and range files available for these simulations
  merge_output_files(dirpath, "Pop")
  merge_output_files(dirpath, "Range")
  
  # Quick plot of abundance and occupancy for these simulations
  jpeg(file=sprintf("%s/Outputs/%s", dirpath, "abundance-and-occupancy.jpg"))
  par(mfrow=c(1,2))
  plotAbundance(sim_params, sprintf("%s/", dirpath))
  plotOccupancy(sim_params, sprintf("%s/", dirpath))
  dev.off()
  
}


merge_output_files <- function(
    dirpath,
    out_save_type
) {
  
  # Get lists of all reps, years, batches, sims and land numbers
  list_sim_files <- list.files(
    path = file.path(dirpath, "Outputs", "file-per-pause"),
    pattern = sprintf("rep\\d+_year\\d+_Batch\\d+_Sim\\d+_Land\\d+_%s.txt", out_save_type)
  )
  output_file_regex <- sprintf("rep(\\d+)_year(\\d+)_Batch(\\d+)_Sim(\\d+)_Land(\\d+)_%s.txt", out_save_type)
  rep_list <- unique(as.numeric(gsub(output_file_regex, "\\1", list_sim_files)))
  rep_list= rep_list[order(rep_list)]
  year_list <- unique(as.numeric(gsub(output_file_regex, "\\2", list_sim_files)))
  year_list = year_list[order(year_list)]
  batch_list <- unique(as.numeric(gsub(output_file_regex, "\\3", list_sim_files)))
  batch_list = batch_list[order(batch_list)]
  sim_list <- unique(as.numeric(gsub(output_file_regex, "\\4", list_sim_files)))
  sim_list = sim_list[order(sim_list)]
  land_list <- unique(as.numeric(gsub(output_file_regex, "\\5", list_sim_files)))
  land_list = land_list[order(land_list)]
  
  for (land_n in land_list) {
    for (sim_n in sim_list) {
      for (batch_n in batch_list) {
        
        data_list <- list()
        i <- 1
        
        # For each replicate (i.e. single simulation, single land, single batch)
        # get list of all output files with suffix `out_save_type` (Pop or Range)
        for (rep_n in rep_list) {
          
          current_file_list <- list_sim_files[grepl(
            sprintf("rep%s_year(\\d+)_Batch%s_Sim%s_Land%s_%s.txt", rep_n, batch_n, sim_n, land_n, out_save_type),
            list_sim_files
          )]
          
          for (year in year_list) {
            
            data_txt_path <- file.path(
              dirpath,
              "Outputs",
              "file-per-pause",
              sprintf("rep%s_year%s_Batch%s_Sim%s_Land%s_%s.txt", rep_n, year, batch_n, sim_n, land_n, out_save_type)
            )
            
            # If an output file exists for this year, load it, keep year 0 data
            # (as this is the first year of change), and renumber the year
            # to match the pausing count of years from the file name
            if (file.exists(data_txt_path)) {
              df_data <- read.table(data_txt_path, sep="\t", header=T)
              if (year < max(year_list)) {
                df_data <- subset(df_data, Year == 0)
              }
              if (nrow(df_data) > 0) {
                df_data$Year <- df_data$Year + (year - 1)
                df_data$Rep <- rep_n
              }
              data_list[[i]] <- df_data
              
              i <- i + 1
            }
          }
        }
        
        if (length(data_list) > 0) {
          df_data = do.call(rbind, data_list)
          write.table(
            df_data,
            file.path(dirpath, "Outputs", sprintf("Batch%s_Sim%s_Land%s_%s.txt", batch_n, sim_n, land_n, out_save_type)),
            row.names = F,
            quote = F,
            sep = '\t'
          )
        }
      }
    }
  }
  
}


load_param <- function(param_name, batch_row, base_config) {
  # Load the param, trying first to load it from the batch.
  # If it is not specified in the batch (NULL) then load it
  # from the base config.
  #param <- batch_row[param_name] %||% base_config[param_name]
  param <- tryCatch(batch_row[param_name], error = function(e) return(base_config[param_name]))
  if (is.null(param)) {
    stop(sprintf("Error: value for `%s` not provided in config or batch file", param_name))
  }
  return(param)
}


prepare_sim_dirs <- function(dirpath, dynamic_landscapes, dynamic_cost_maps, params) {
  
  dir.create(dirpath, showWarnings = TRUE)
  dir.create(file.path(dirpath, "Inputs"), showWarnings = TRUE)
  dir.create(file.path(dirpath, "Outputs"), showWarnings = TRUE)
  dir.create(file.path(dirpath, "Outputs", "file-per-pause"),
             showWarnings = TRUE)
  dir.create(file.path(dirpath, "Output_Maps"), showWarnings = TRUE)
  
  for (landscape_path in unique(dynamic_landscapes)) {
    file.copy(file.path(landscape_path),
              file.path(dirpath, "Inputs"))
    file.copy(sprintf("%s.tif", tools::file_path_sans_ext(landscape_path)),
              file.path(dirpath, "Inputs"))
  }
  for (cost_map_path in unique(dynamic_cost_maps)) {
    file.copy(file.path(cost_map_path),
              file.path(dirpath, "Inputs"))
    file.copy(sprintf("%s.tif", tools::file_path_sans_ext(cost_map_path)),
              file.path(dirpath, "Inputs"))
  }
  
  file.copy(file.path(params$init_inds_file),
            file.path(dirpath, "Inputs"))
}


run_sims <- function(
  x,
  base_config,
  input_data_path,
  output_path
) {
  
  params <- modifyList(base_config, x)
  
  dirpath <- file.path(output_path, params$batch_name)
  
  if (!file.exists(dirpath)) {
    
    aggregation_factor <- 1
    pixel_size <- 25 * aggregation_factor
    rangeshiftr_steps <- 1
    
    df_phenips <- read.csv(params$phenips_file)
    phenips <- df_phenips$Generations
    
    df_landscapes <- read.csv(params$dynamic_landscapes)
    dynamic_landscapes <- df_landscapes$landscape
    dynamic_cost_maps <- df_landscapes$cost_map
    
    prepare_sim_dirs(dirpath, dynamic_landscapes, dynamic_cost_maps, params)
      
    land <- ImportedLandscape(
      LandscapeFile = basename(dynamic_landscapes[1]),
      CostsFile = basename(dynamic_cost_maps[1]),
      Resolution = pixel_size,
      HabPercent = TRUE,
      K_or_DensDep = params$carrying_capacity
    )
    
    demo <- Demography(
      ReproductionType = 1,  # sexual
      Rmax = params$rmax,
      bc = 2
    )
    
    disp <- Dispersal(
      Emigration = Emigration(
        DensDep = TRUE,
        EmigProb = matrix(c(params$D0, params$alpha, params$beta), nrow = 1)
      ),
      Transfer = SMS(
        DP = params$directional_persistence,
        PR = params$perceptual_range,
        MemSize = params$mem_size,
        GoalType = params$goal_type,
        GoalBias = params$goal_bias,
        AlphaDB = params$alpha_density_bias,
        BetaDB = params$beta_density_bias,
        Costs = "file",
        StepMort = params$step_mortality
      ),
      Settlement = Settlement(
        StageDep = FALSE,
        SexDep = TRUE,
        DensDep = TRUE,
        FindMate = c(TRUE, FALSE),
        # Female requires mate, male does not
        Settle = matrix(
          c(0, params$max_settlement_prob, params$density_inflection, params$density_slope,  # Female: αS = dens_inflection, βS = dens_slope
            1, params$max_settlement_prob, params$density_inflection, params$density_slope),   # Male: same),
          ncol = 4,
          byrow = TRUE
        ),
        MinSteps = params$min_steps,
        MaxSteps = params$max_steps,
        MaxStepsYear = params$max_steps
      )
    )
    
    
    init <- Initialise(
      InitType = 2, # Initialise from individuals list
      InitIndsFile = basename(params$init_inds_file)
    )
    
    
    sim <- Simulation(
      Simulation = 0,
      Years = rangeshiftr_steps,
      Replicates = 1, # one here because of pausing (controlled with pausingextras)
      OutIntPop = 1, # one here because of pausing (controlled with pausingextras)
      OutIntRange = 1 # one here because of pausing (controlled with pausingextras)
    )
    
    sim_params <-
      RSsim(
        batchnum = 0,
        land = land,
        demog = demo,
        dispersal = disp,
        simul = sim,
        init = init
      )  # seed to allow replication
    
    # Set up the extra parameters used as part of RangeShiftR pausing
    attributes(sim_params)$pausingextras <- PausingExtras(
      replicates = params$replicates,
      years = params$years,
      pop_impacts_quality = params$pop_impacts_quality,
      years_for_destruction = params$years_for_destruction,
      apply_destruction_annually = params$apply_destruction_annually,
      phenips = as.list(phenips),
      dynamic_landscapes_list = as.list(dynamic_landscapes),
      dynamic_cost_maps_list = as.list(dynamic_cost_maps)
    )
    
    validateRSparams(sim_params)
    
    run_pausing_simulations(
      sim_params = sim_params,
      input_data_path = input_data_path,
      dirpath = dirpath
    )
  }
  
}
