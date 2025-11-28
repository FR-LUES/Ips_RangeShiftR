
# ============================================================
# Script Name: Landscape Habitat Suitability 
# Author: Fiona Plenderleith

# Description:
# This script generates a time series of habitat suitability maps
# for Ips typographus using spatial data on spruce distribution,
# topographic wetness index (TWI), wind damage (ForestGALES),
# and drought stress (SPEI). 
#
# Inputs:
#   - Spruce shapefile
#   - TWI raster
#   - Subcompartment shapefile
#   - ForestGALES wind damage predictions
#   - SPEI drought index CSV
#
# Outputs:
#   - Yearly habitat suitability rasters (.asc and .tif)
#   - Initial population distribution files for simulation
#
# Dependencies:
#   - terra
#   - sf
#   - dplyr
#   - splines
#   - Ecdat
#
# Notes:
#   - CRS should be British National Grid (EPSG:27700)


# to update
# we need to only have patches exist for 2 or 3 years in a landscape ///?

# Main wrapper function to generate habitat suitability maps over time#function that creates landscape


# ===============================# Main Function: Landscape Time Series Generator# ===============================
create_time_series_maps <- function(Spruce_path,
                                    #path to spruce data
                                    TWI_path,
                                    #path to TWI
                                    Subcompart_path,
                                    #path to subcompartment data
                                    ForestGales_path,
                                    #path to forestgales data
                                    SPEI_data,
                                    # path to SPEI data summary by year
                                    Output_path,
                                    # path to store outputs
                                    extent,
                                    #extent of landscape
                                    crs,
                                    #crs, recoment BNG "EPSG:27700"
                                    TWI_min,
                                    #min cut off for TWI
                                    TWI_max,
                                    #max cut off for TWI
                                    years,
                                    #which years to run model,
                                    threshold) {
  # this function will take .geojson or .tif file and rasterize and clip to specified extent and then convert to specified crs
  
  # --- Rasterization and Clipping Utility ---
  
  Clip_and_rasterize <- function(Map_input,
                                 extent,
                                 crs,
                                 landuse = NA,
                                 res = 25) {
    library(terra)
    library(sf)
    
    # Create target raster template
    clipped_extent <- terra::ext(extent)
    r <- terra::rast(res = res, ext = clipped_extent, crs = crs)
    
    # Detect file type
    is_geo <- grepl("\\.geo", Map_input, ignore.case = TRUE)
    is_shp <- grepl("\\.shp", Map_input, ignore.case = TRUE)
    
    if (is_geo || is_shp) {
      # Read vector data
      OGR_Data <- st_read(dsn = Map_input, quiet = TRUE)
      Vect_Data <- vect(OGR_Data)
      
      # Check for columns
      has_PRILANDUSE <- "PRILANDUSE" %in% names(Vect_Data)
      has_failed <- any(grepl("^failed_", names(Vect_Data)))
      
      # Rasterization logic
      if (has_PRILANDUSE && !is.na(landuse)) {
        # Filter by landuse if provided
        Vect_Data_subset <- subset(Vect_Data, Vect_Data$PRILANDUSE == landuse)
        Map_rast <- terra::rasterize(Vect_Data_subset, r, background = NA)
        
      } else if (has_failed && !is.na(landuse)) {
        # Use 'failed_' field if relevant
        Map_rast <- terra::rasterize(Vect_Data, r, field = landuse, background = NA)
        
      } else {
        # No filtering, rasterize entire dataset
        Map_rast <- terra::rasterize(Vect_Data, r, background = NA)
      }
      
    } else {
      # Input is raster
      Map_rast <- rast(Map_input)
      Map_rast <- resample(Map_rast, r, method = "near")
    }
    
    # Assign CRS and return
    crs(Map_rast) <- crs
    return(Map_rast)
  }
  
  #create baseline map that combines the healthy spruce, subcompartment details and the TWI
  
  # --- Baseline Habitat Map Creation ---
  Create_baseline_map <-
    function(Spruce_path,
             TWI_path,
             Subcompart_path,
             extent,
             crs,
             TWI_min,
             TWI_max) {
      Spruce <- Clip_and_rasterize(Spruce_path, extent, crs) # 21 or 29
      
      # Spruce[!(Spruce %in% c(21, 29))] <- NA
      # Spruce[is.na(Spruce)] <- 0
      #
      
      TWI <- Clip_and_rasterize(TWI_path, extent, crs)
      
      # Extract subcompartment types: windblow and felled      # get habtiat types
      sub_compart_wind <-
        Clip_and_rasterize(Subcompart_path, extent, crs, landuse = "windblow")
      sub_compart_felled <-
        Clip_and_rasterize(Subcompart_path, extent, crs, landuse = "Felled")
      
      sub_compart_felled[is.na(Spruce)] <- 0 # make cells that are not spruce NA
      sub_compart_wind[is.na(Spruce)] <- 0 # make cells that are not spruce NA
      
      # Process TWI to identify drought-prone areas      # TWI data - these are going to be drought risk areas
      TWI[TWI[] < TWI_min] <- 100 # select only extreme values
      TWI[TWI[] > TWI_max] <- 100
      TWI[TWI[] != 100] <- 0 # everything inbetween is not drouight risk
      TWI[is.na(Spruce)] <- 0 # make cells that are not spruce NA
      
      Spruce[TWI == 100] <- 4 # these are drought risk areas
      Spruce[!is.na(sub_compart_wind) & sub_compart_wind == 1] <- 2
      Spruce[!is.na(sub_compart_felled) &
               sub_compart_felled == 1] <- 3
      
      
      
      return(Spruce)
    }
  
  # Integrate wind damage predictions into baseline map  # now add the predicted yearly wind damage to the baseline map and save each of these
  
  # --- Annual Wind Damage Integration ---
  
  
  
  Add_yearly_wind_damage <- function(years,
                                     ForestGales_path,
                                     Output_path,
                                     Baseline_map,
                                     threshold = 2) {
    map_extent <- terra::ext(Baseline_map)
    crs_map <- crs(Baseline_map)
    
    # Initialize age raster
    windblow_age <- rast(Baseline_map)
    windblow_age[] <- NA
    
    for (y in years) {
      # Increment age for windblown cells
      windblow_age <- ifel(!is.na(windblow_age), windblow_age + 1, NA)
      
      # Expire old windblown patches: if age >= expiry_years, revert to NA
      Baseline_map <- ifel(windblow_age >= threshold, NA, Baseline_map)
      
      # Clip and rasterize ForestGales for this year
      field_choice <- paste("failed", y, sep = "_")
      ForestGales <- Clip_and_rasterize(ForestGales_path, map_extent, crs_map, landuse = field_choice)
      
      # Recode ForestGales values
      ForestGales <- classify(ForestGales, rbind(c(0, NA), c(1, 5)))
      ForestGales[is.na(Baseline_map)] <- NA
      
      # Apply new windblow only to cells that are still trees (Baseline_map > 0) and never windblown before
      new_windblow <- ForestGales == 5 &
        Baseline_map > 0 & is.na(windblow_age)
      Baseline_map <- ifel(new_windblow, 5, Baseline_map)
      
      # Set age = 0 for new windblown cells
      windblow_age <- ifel(new_windblow, 0, windblow_age)
      
      # Save outputs
      filename_tif <- paste0(Output_path, y, "_habcodes.tif")
      filename_asc <- paste0(Output_path, y, "_habcodes.asc")
      
      writeRaster(Baseline_map, filename_tif, overwrite = TRUE)
      writeRaster(
        Baseline_map,
        filename_asc,
        overwrite = TRUE,
        NAflag = -999,
        datatype = "INT2S"
      )
      
      plot(Baseline_map, main = paste("Year:", y))
    }
  }
  
  # Convert habitat codes to quality scores based on wind and drought  #This duntion will convert the habitat codes to quality
  
  # requires average SPEI projection and year for area
  
  # --- Habitat Quality Conversion: Wind + Drought ---
  HabitatCodetoIpsQual_both <- function(years, SPEI_data, Output_path) {
    for (y in years) {
      filename <-
        paste0(Output_path, y, "_habcodes.asc", sep = "")
      
      Map_rast <- rast(filename)
      Map_rast[Map_rast[] == 1] <- 0 # LIVING wont sustain a pop
      Map_rast[Map_rast[] == 2] <- 100 # wind blown
      Map_rast[Map_rast[] == 3] <- 10 # felled
      Map_rast[Map_rast[] == 5] <- 100 # high chance of wind damage
      #
      # # now do the drought damage which depends on the year
      library(splines)
      library(Ecdat)
      library(dplyr)
      
      
      # Training data for spline
      SPEI <- seq(from = -4,
                  to = 0,
                  by = 0.5)
      qual <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)
      
      # Fit spline model
      model <- lm(qual ~ bs(SPEI, degree = 3))
      
      # Read real SPEI data
      SPEI_local <- read.csv(SPEI_data)
      
      # Predict quality
      SPEI_local$pred <- predict(model, newdata = list(SPEI = SPEI_local$SPEI))
      
      # Round predictions
      SPEI_local$pred <- round(SPEI_local$pred, digits = 0)
      
      # Clamp predictions to 0–100
      SPEI_local$pred <- pmax(pmin(SPEI_local$pred, 100), 0)
      
      # Apply rule: if SPEI > 0, set quality = 0
      SPEI_local$pred[SPEI_local$SPEI > 0] <- 0
      
      # Handle NA
      SPEI_local$pred[is.na(SPEI_local$pred)] <- 0
      
      # Extract value for current year
      SPEI_value <- subset(SPEI_local, SPEI_local$year == y)$pred
      
      
      Map_rast[Map_rast[] == 4] <- SPEI_value # high chance of wind damage
      # spruce_raster[spruce_raster[] == 5] <- SPEI_value
      
      # now just save it
      filenames_hab <-
        paste0(Output_path, y, "_quality_both.asc", sep = "")
      filenames_hab_tif <-
        paste0(Output_path, y, "_quality_both.tif", sep = "")
      writeRaster(
        Map_rast,
        filenames_hab,
        NAflag = -999,
        overwrite = TRUE,
        datatype = "INT2S"
      )
      writeRaster(Map_rast, filenames_hab_tif, overwrite = TRUE)
      freq(Map_rast)
      plot(
        Map_rast,
        col = c("gray94", "green3", "darkgreen", "darkgoldenrod3", "red"),
        main = y
      )
    }
  }
  
  # --- Habitat Quality Conversion: Wind Only ---
  HabitatCodetoIpsQual_wind <- function(years, SPEI_data, Output_path) {
    for (y in years) {
      filename <-
        paste0(Output_path, y, "_habcodes.asc", sep = "")
      library(splines)
      library(Ecdat)
      library(dplyr)
      
      yearone <- min(years)
      Map_rast <- rast(filename)
      Map_rast[Map_rast[] == 1] <- 0 # LIVING wont sustain a pop
      Map_rast[Map_rast[] == 2] <- 100 # wind blown
      Map_rast[Map_rast[] == 3] <- 10 # felled
      Map_rast[Map_rast[] == 5] <- 100 # high chance of wind damage
      
      
      
      Map_rast[Map_rast[] == 4] <- 0 # high chance of wind damage
      #
      
      # now just save it
      filenames_hab <-
        paste0(Output_path, y, "_quality_wind.asc", sep = "")
      filenames_hab_tif <-
        paste0(Output_path, y, "_quality_wind.tif", sep = "")
      writeRaster(
        Map_rast,
        filenames_hab,
        NAflag = -999,
        overwrite = TRUE,
        datatype = "INT2S"
      )
      writeRaster(Map_rast, filenames_hab_tif, overwrite = TRUE)
      freq(Map_rast)
      plot(
        Map_rast,
        col = c("gray94", "green3", "darkgreen", "darkgoldenrod3", "red"),
        main = y
      )
    }
  }
  
  # --- Habitat Quality Conversion: Drought Only ---
  HabitatCodetoIpsQual_drought <- function(years, SPEI_data, Output_path) {
    yearone <- min(years)
    filename <-
      paste0(Output_path, yearone, "_habcodes.asc", sep = "")
    
    Map_rast <- rast(filename)
    Map_rast[Map_rast[] == 1] <- 0 # LIVING wont sustain a pop
    Map_rast[Map_rast[] == 2] <- 0 # wind blown
    Map_rast[Map_rast[] == 3] <- 10 # felled
    Map_rast[Map_rast[] == 5] <- 0 # high chance of wind damage
    #
    # # now do the drought damage which depends on the year
    library(splines)
    library(Ecdat)
    library(dplyr)
    
    
    for (y in years) {
      Map_rast_y <- Map_rast
      
      # Training data for spline
      SPEI <- seq(from = -4,
                  to = 0,
                  by = 0.5)
      qual <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)
      
      # Fit spline model
      model <- lm(qual ~ bs(SPEI, degree = 3))
      
      # Read real SPEI data
      SPEI_local <- read.csv(SPEI_data)
      
      # Predict quality
      SPEI_local$pred <- predict(model, newdata = list(SPEI = SPEI_local$SPEI))
      
      # Round predictions
      SPEI_local$pred <- round(SPEI_local$pred, digits = 0)
      
      # Clamp predictions to 0–100
      SPEI_local$pred <- pmax(pmin(SPEI_local$pred, 100), 0)
      
      # Apply rule: if SPEI > 0, set quality = 0
      SPEI_local$pred[SPEI_local$SPEI > 0] <- 0
      
      # Handle NA
      SPEI_local$pred[is.na(SPEI_local$pred)] <- 0
      
      # Extract value for current year
      SPEI_value <- subset(SPEI_local, SPEI_local$year == y)$pred
      
      Map_rast_y[Map_rast_y[] == 4] <- SPEI_value # high chance of wind damage
      # spruce_raster[spruce_raster[] == 5] <- SPEI_value
      
      # now just save it
      filenames_hab <-
        paste0(Output_path, y, "_quality_drought.asc", sep = "")
      filenames_hab_tif <-
        paste0(Output_path, y, "_quality_drought.tif", sep = "")
      writeRaster(
        Map_rast_y,
        filenames_hab,
        NAflag = -999,
        overwrite = TRUE,
        datatype = "INT2S"
      )
      writeRaster(Map_rast_y, filenames_hab_tif, overwrite = TRUE)
      
      plot(
        Map_rast_y,
        col = c("gray94", "green3", "darkgreen", "darkgoldenrod3", "red"),
        main = y
      )
    }
  }
  # code  1 = 0 % # LIVING wont sustain a pop
  # code  2 = 100 % # wind blown
  # code  3 = 10 % # felled
  # code  5 = 100 % # high chance of wind damage
  # code  4 = XX % # drought risk area is dependent on the relationship modelled between SPEI and qual
  # SPEI <- seq(from = -4,to = 0,by =0.5)
  # qual <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)
  # Generate initial species distribution files for simulation input  # finally a function to create Init sp dist files - pick random suitablle pop from first year
  
  Baseline <- Create_baseline_map(Spruce_path,
                                  TWI_path,
                                  Subcompart_path,
                                  extent,
                                  crs,
                                  TWI_min,
                                  TWI_max)
  
  Add_yearly_wind_damage(years, ForestGales_path, Output_path, Baseline, threshold)
  
  HabitatCodetoIpsQual_both(years, SPEI_data, Output_path)
  HabitatCodetoIpsQual_drought(years, SPEI_data, Output_path)
  HabitatCodetoIpsQual_wind(years, SPEI_data, Output_path)
  
}


create_time_series_maps(
  Spruce_path =  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/FR tree spp/spruce_BNG25_near.tif",
  TWI_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/TWI/TWI_focal.tif",
  Subcompart_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/subcomprtments_25/Forestry_England_Subcompartments/Forestry_England_Subcompartments.shp", #update
  ForestGales_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/ForestDean_Failure_40%.geojson",
  SPEI_data = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Forest_of_Dean.csv", #to update'
  Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoDTest1/", 
  extent = c(350000, 375000, 200000, 225000),
  crs =  "EPSG:27700",
  TWI_min = 5,
  TWI_max = 14,
  years = 2020:2050, 
  threshold = 2
)


species_init_map <- function(Output_path, NumberCells, start_years, reps, buffer_dist) {
  all_summary <- data.frame()
  
  for (start_year in start_years) {
    # Load quality and habitat maps for this start year
    quality_file <- paste0(Output_path, start_year, "_quality_both.asc")
    habitat_file <- paste0(Output_path, start_year, "_habcodes.asc")
    
    if (!file.exists(quality_file) || !file.exists(habitat_file)) {
      warning("Missing files for year: ", start_year)
      next
    }
    
    Map_rast <- rast(quality_file)
    habitat_rast <- rast(habitat_file)
    
    # Align habitat raster if needed
    if (!compareGeom(habitat_rast, Map_rast, stopOnError = FALSE)) {
      habitat_rast <- project(habitat_rast, Map_rast)
    }
    
    for (rep in 1:reps) {
      # Select valid cells (>0 and not NA in habitat)
      valid_cells <- which(Map_rast[] > 0 & !is.na(habitat_rast[]))
      if (length(valid_cells) < NumberCells) {
        stop("Not enough cells with value > 0 to select.")
      }
      
      selected_indices <- sample(valid_cells, NumberCells)
      
      # Create init raster
      spruce_raster_init <- Map_rast
      spruce_raster_init[] <- NA
      spruce_raster_init[selected_indices] <- 1
      names(spruce_raster_init) <- "pop"
      
      # Save initial distribution raster
      raster_name <- paste0("initial_dist_", start_year, "_rep", rep, ".tif")
      terra::writeRaster(spruce_raster_init, file.path(Output_path, raster_name), overwrite = TRUE, datatype = "INT2S")
      
      # Buffer summary
      init_points <- as.points(spruce_raster_init, values = TRUE, na.rm = TRUE)
      buffers <- terra::buffer(init_points, width = buffer_dist)
      
      summary_list <- lapply(1:nrow(buffers), function(i) {
        cells <- terra::extract(habitat_rast, buffers[i, ], cells = TRUE)
        codes <- cells[[1]][!is.na(cells[[1]])]
        coords <- terra::geom(init_points)[i, c("x", "y")]
        
        data.frame(
          StartYear = start_year,
          Rep = rep,
          InitSite = i,
          X = coords[1],
          Y = coords[2],
          Code1 = sum(codes == 1),
          Code2 = sum(codes == 2),
          Code3 = sum(codes == 3),
          Code4 = sum(codes == 4),
          Code5 = sum(codes == 5),
          TotalCells = length(codes)
        )
      })
      
      all_summary <- rbind(all_summary, do.call(rbind, summary_list))
      
      # Individuals file
      pop_size <- 2
      spruce_raster_init_pnt <- as.points(spruce_raster_init, values = TRUE, na.rm = TRUE)
      
      create_inds_dataframe <- function(shp_species, rst_species) {
        df_inds <- terra::extract(rst_species, shp_species, cells = TRUE, df = TRUE)
        df_inds$Year <- 0
        df_inds$Species <- 0
        df_inds$X <- df_inds$cell %% ncol(rst_species) - 1
        df_inds$Y <- nrow(rst_species) - (floor(df_inds$cell / ncol(rst_species))) - 1
        df_inds$Ninds <- df_inds$pop
        df_inds <- df_inds[, !(names(df_inds) %in% c('ID', 'cells', 'pop'))]
        df_inds <- df_inds[!is.na(df_inds$Ninds), ]
        df_inds <- unique(df_inds)
        df_inds <- df_inds[df_inds$Ninds != 0, ]
        return(df_inds)
      }
      
      df_individuals_females <- create_inds_dataframe(spruce_raster_init_pnt, spruce_raster_init)
      df_individuals_females$Ninds <- df_individuals_females$Ninds * pop_size / 2
      df_individuals_females$Sex <- 1
      df_individuals_females <- df_individuals_females[, c("Year", "Species", "X", "Y", "Ninds", "Sex")]
      
      df_individuals_males <- create_inds_dataframe(spruce_raster_init_pnt, spruce_raster_init)
      df_individuals_males$Ninds <- df_individuals_males$Ninds * pop_size / 2
      df_individuals_males$Sex <- 0
      df_individuals_males <- df_individuals_males[, c("Year", "Species", "X", "Y", "Ninds", "Sex")]
      
      df_individuals <- rbind(df_individuals_females, df_individuals_males)
      inds_file <- paste0("inds_", start_year, "_rep", rep, ".txt")
      write.table(df_individuals, file.path(Output_path, inds_file), row.names = FALSE, quote = FALSE, sep = '\t')
    }
  }
  
  # Save combined summary
  summary_file <- paste0(Output_path, "init_summary_all.csv")
  write.csv(all_summary, summary_file, row.names = FALSE)
  
  return(all_summary)
}

species_init_map(Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoDTest1/", 
                 start_years = c(2020, 2030, 2040), 
                 NumberCells = 1, 
                 reps = 2, 
                 buffer_dist = 500)

