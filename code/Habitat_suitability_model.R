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
                                    PatchAge_path,
                                    #path to forestgales data - processed data - age raster
                                    Felling_raster_path,
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
                                    years) {
                                    #which years to run model) {
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
                                     PatchAge_path,
                                     Output_path,
                                     Baseline_map,
                                     final_year = 2050) {
    
    # Load patch age raster (age in final year)
    patch_age_final <- rast(PatchAge_path)
    patch_age_final <- crop(patch_age_final, Baseline_map)
    patch_age_final <- resample(patch_age_final, Baseline_map, method = "near")
    
    # Loop through simulation years
    for (y in years) {
      
      # Compute age in this year
      # Example: if age=30 in 2050 → age=0 in 2020
      # Age in this simulation year
      age_this_year <- patch_age_final - (final_year - y)
      
      # Negative values become NA
      age_this_year[age_this_year < 0] <- NA
      
      # Start with NA
      wind_class <- age_this_year
      wind_class[] <- NA
      
      # Fresh windblow (0–2)
      wind_class <- ifel(age_this_year >= 0 & age_this_year <= 2, 5, wind_class)
      
      # Rotten (>=3)
      wind_class <- ifel(age_this_year >= 3, 6, wind_class)
      
      # Start from baseline each year
      out_map <- Baseline_map
      
      # Only living (1) or drought (4) can be overwritten
      eligible <- out_map %in% c(1, 4)
      
      # Apply windblow class
      out_map <- ifel(eligible & !is.na(wind_class), wind_class, out_map)
      
      # Save outputs
      filename_tif <- paste0(Output_path, y, "_habcodes.tif")
      filename_asc <- paste0(Output_path, y, "_habcodes.asc")
      
      writeRaster(out_map, filename_tif, overwrite = TRUE)
      writeRaster(out_map, filename_asc,
                  overwrite = TRUE,
                  NAflag = -999,
                  datatype = "INT2S")
      
      plot(out_map, main = paste("Year:", y))
    }
  }
  Add_yearly_felling <- function(years,
                                 Felling_raster_path,
                                 Baseline,
                                 Output_path) {
    
    # Load felling-year raster
    r_fell_year <- rast(Felling_raster_path)
    r_fell_year <- crop(r_fell_year, Baseline)
    r_fell_year <- resample(r_fell_year, Baseline, method = "near")
    
    # First simulation year
    first_year <- min(years)
    
    # Combine scheduled felling with pre-existing 3s
    fell_year_combined <- r_fell_year
    
    # Pre-existing 3s in baseline are treated as felled in first year
    preexisting_3 <- Baseline == 3
    fell_year_combined <- ifel(preexisting_3, first_year, fell_year_combined)
    
    for (y in years) {
      
      # Load existing habcodes (already includes wind)
      hab <- rast(paste0(Output_path, y, "_habcodes.tif"))
      
      # Identify stands scheduled for felling this year
      felling_now <- fell_year_combined == y
      
      # Eligible habitat codes for felling
      eligible <- hab %in% c(1, 4)
      
      # Apply fresh felling (class 3)
      hab <- ifel(eligible & felling_now, 3, hab)
      
      # --- AGE FELLED STANDS (3 → 6 after 3 years) ---
      
      years_since_fell <- y - fell_year_combined
      
      # Fresh felling (0–2 years old)
      fresh_fell <- years_since_fell >= 0 & years_since_fell <= 2
      
      # Rotten felling (≥3 years old)
      rotten_fell <- years_since_fell >= 3
      
      # Only age cells that are currently 3
      is_currently_3 <- hab == 3
      
      # Fresh felling stays 3
      hab <- ifel(is_currently_3 & fresh_fell, 3, hab)
      
      # Rotten felling becomes 6
      hab <- ifel(is_currently_3 & rotten_fell, 6, hab)
      
      
      # ------------------------------------------------
      
      # Save updated maps
      writeRaster(hab,
                  paste0(Output_path, y, "_habcodes_fell.tif"),
                  overwrite = TRUE)
      
      writeRaster(hab,
                  paste0(Output_path, y, "_habcodes_fell.asc"),
                  overwrite = TRUE,
                  NAflag = -999,
                  datatype = "INT2S")
      
      plot(hab, main = paste("Year:", y, "- Wind + Felling + Ageing"))
    }
  }

  # --- Habitat Quality Conversion: Wind + Drought ---
  HabitatCodetoIpsQual_both <- function(years, SPEI_data, Output_path) {
    for (y in years) {
      filename <-
        paste0(Output_path, y, "_habcodes_fell.asc", sep = "")
      
      Map_rast <- rast(filename)
      Map_rast[Map_rast[] == 1] <- 0 # LIVING wont sustain a pop
      Map_rast[Map_rast[] == 2] <- 100 # wind blown
      Map_rast[Map_rast[] == 3] <- 10 # felled
      Map_rast[Map_rast[] == 5] <- 100 # high chance of wind damage
      Map_rast[Map_rast[] == 6] <- 0 # Rotton wood
      #
      # Load real SPEI data
      SPEI_local <- read.csv(SPEI_data)
      SPEI_value <- subset(SPEI_local, year == y)$pred
      
      
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
        paste0(Output_path, y, "_habcodes_fell.asc", sep = "")
      library(splines)
      library(Ecdat)
      library(dplyr)
      
      yearone <- min(years)
      Map_rast <- rast(filename)
      Map_rast[Map_rast[] == 1] <- 0 # LIVING wont sustain a pop
      Map_rast[Map_rast[] == 2] <- 100 # wind blown
      Map_rast[Map_rast[] == 3] <- 10 # felled
      Map_rast[Map_rast[] == 5] <- 100 # high chance of wind damage
      
      Map_rast[Map_rast[] == 6] <- 0 # Rotten wood
      
      Map_rast[Map_rast[] == 4] <- 0 # drought risk - treated as healpthy
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

  
  HabitatCodetoIpsQual_10yWindFell <- function(
    years,         # target output years, e.g. 2020:2049
    SPEI_data,     # CSV file path with columns: year, pred (already 0–100)
    Output_path,   # output directory
    base_years     # source structure years to replay, e.g. 2020:2029 (10-year window)
  ) {
  
       SPEI_local <- read.csv(SPEI_data)
    if (!all(c("year", "pred") %in% names(SPEI_local))) {
      stop("SPEI_data must contain columns 'year' and 'pred'.")
    }
    
    # Small helper: for any target year y, pick the source baseline year ys
    # using modulo alignment over the base 10-year window
    pick_source_year <- function(y) {
      bmin <- min(base_years)            # e.g., 2020
      wlen <- length(base_years)         # should be 10
      idx  <- (y - bmin) %% wlen + 1     # 1..10 for 2020..2049
      base_years[idx]
    }
    
    for (i in seq_along(years)) {
      y  <- years[i]              # target year (climate)
      ys <- pick_source_year(y)   # source year (structure)
      
      # Read the source structure map (wind + felling + ageing) from the base window
      filename_src <- paste0(Output_path, ys, "_habcodes_fell.asc")
      if (!file.exists(filename_src)) {
        stop("Missing base structure file: ", filename_src,
             "\nRun your wind+felling pipeline for the base window first.")
      }
      
      Map_rast <- rast(filename_src)
      
      # Fixed class mappings (same as HabitatCodetoIpsQual_both)
      Map_rast[Map_rast[] == 1] <- 0   # LIVING won't sustain a pop
      Map_rast[Map_rast[] == 2] <- 100 # wind blown
      Map_rast[Map_rast[] == 3] <- 10  # felled
      Map_rast[Map_rast[] == 5] <- 100 # high chance of wind damage
      Map_rast[Map_rast[] == 6] <- 0   # Rotten wood
      
      # Target-year SPEI (already the final 0–100 value)
      SPEI_value <- subset(SPEI_local, year == y)$pred
  
      # Apply drought-only difference (code 4) via target year SPEI
      Map_rast[Map_rast[] == 4] <- SPEI_value
      
      # Save outputs (distinct suffix to keep separate from *_quality_both.*)
      filenames_hab     <- paste0(Output_path, y, "_quality_10yWindFell.asc")
      filenames_hab_tif <- paste0(Output_path, y, "_quality_10yWindFell.tif")
      
      writeRaster(
        Map_rast,
        filenames_hab,
        NAflag = -999,
        overwrite = TRUE,
        datatype = "INT2S"
      )
      writeRaster(Map_rast, filenames_hab_tif, overwrite = TRUE)
      
      # Optional diagnostics (kept in your style)
      freq(Map_rast)
      plot(
        Map_rast,
        col = c("gray94", "green3", "darkgreen", "darkgoldenrod3", "red"),
        main = paste0(y, " (src=", ys, ")")
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
  
  Add_yearly_wind_damage(years, PatchAge_path, Output_path, Baseline,
                         final_year = 2050)
  
  
  Add_yearly_felling(years, Felling_raster_path,Baseline, Output_path)
  HabitatCodetoIpsQual_both(years, SPEI_data, Output_path)
  HabitatCodetoIpsQual_wind(years, SPEI_data, Output_path)
  HabitatCodetoIpsQual_10yWindFell(years = 2020:2049, SPEI_data, Output_path, base_years = 2020:2029)
}


create_time_series_maps(
  Spruce_path =  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/FR tree spp/spruce_BNG25_near.tif",
  TWI_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/TWI/TWI_focal.tif",
  Subcompart_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/subcomprtments_25/Forestry_England_Subcompartments/Forestry_England_Subcompartments.shp", #update
  PatchAge_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/IpsShift_Kielder.tif",
  Felling_raster_path <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/felled_rast_EoR.tif",
  SPEI_data = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Kielder.csv", #to update'
  Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/Kielder/", 
  extent = c(352643.5, 369793.5, 578384.812999999, 595534.812999999 ),
  crs =  "EPSG:27700",
  TWI_min = 4,
  TWI_max = 13,
  years = 2020:2050
)

create_time_series_maps(
  Spruce_path =  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/FR tree spp/spruce_BNG25_near.tif",
  Felling_raster_path <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/felled_rast_EoR.tif",
  TWI_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/TWI/TWI_focal.tif",
  Subcompart_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/subcomprtments_25/Forestry_England_Subcompartments/Forestry_England_Subcompartments.shp", #update
  PatchAge_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/IpsShift_FoD.tif",
  SPEI_data = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Forest_of_Dean.csv", #to update'
  Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoD1/", 
  extent = c(352888.1876, 370038.1876, 205389.812999999, 222539.812999999 ),
  crs =  "EPSG:27700",
  TWI_min = 4,
  TWI_max = 13,
  years = 2020:2050
)

create_time_series_maps(
  Spruce_path =  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/FR tree spp/spruce_BNG25_near.tif",
  TWI_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/TWI/TWI_focal.tif",
  Felling_raster_path <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/felled_rast_EoR.tif",
  Subcompart_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/subcomprtments_25/Forestry_England_Subcompartments/Forestry_England_Subcompartments.shp", #update
  PatchAge_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/IpsShift_NYM.tif",
  SPEI_data = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_NYM.csv", #to update'
  Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/NYM/", 
  extent = c(476175, 493325, 484860.344000001, 502010.344000001),
  crs =  "EPSG:27700",
  TWI_min = 4,
  TWI_max = 13,
  years = 2020:2050
)

species_init_map <- function(Output_path, NumberCells, start_years, reps, buffer_dist) {
  
  library(terra)
  library(dplyr)
  
  all_summary <- data.frame()
  
  for (start_year in start_years) {
    
    quality_file  <- paste0(Output_path, start_year, "_quality_10yWindFell.tif")
    habitat_file  <- paste0(Output_path, start_year, "_habcodes_fell.asc")
    
    if (!file.exists(quality_file) || !file.exists(habitat_file)) {
      warning("Missing files for year: ", start_year)
      next
    }
    
    Map_rast     <- rast(quality_file)
    habitat_rast <- rast(habitat_file)
    
    if (!compareGeom(habitat_rast, Map_rast, stopOnError = FALSE)) {
      habitat_rast <- project(habitat_rast, Map_rast)
    }
    
    for (rep in 1:reps) {
      
      valid_cells <- which(Map_rast[] > 0 & !is.na(habitat_rast[]))
      if (length(valid_cells) < NumberCells) {
        stop("Not enough valid cells to select.")
      }
      
      selected_indices <- sample(valid_cells, NumberCells)
      
      spruce_raster_init <- Map_rast
      spruce_raster_init[] <- NA
      spruce_raster_init[selected_indices] <- 1
      names(spruce_raster_init) <- "pop"
      
      raster_name <- paste0("initial_dist_", start_year, "_rep", rep, ".tif")
      writeRaster(
        spruce_raster_init,
        file.path(Output_path, raster_name),
        overwrite = TRUE,
        datatype = "INT2S"
      )
      
      # ---- BUFFER SUMMARY (FIXED) ----
      
      init_points <- as.points(spruce_raster_init, values = TRUE, na.rm = TRUE)
      buffers     <- buffer(init_points, width = buffer_dist)
      
      # ✅ Extract ONCE for all buffers
      ext <- extract(habitat_rast, buffers, cells = TRUE)
      
      # ✅ Identify correct habcode column
      hab_col <- names(ext)[!(names(ext) %in% c("ID", "cell"))]
      ext$habcode <- ext[[hab_col]]
      
      # ✅ Remove NA explicitly
      ext <- ext[!is.na(ext$habcode), ]
      
      # ✅ Per-buffer freq (like freq())
      buf_freq <- ext |>
        count(ID, habcode, name = "n_cells")
      
      # ✅ Build summary rows per buffer
      summary_rows <- lapply(seq_len(nrow(buffers)), function(i) {
        
        freq_i <- buf_freq[buf_freq$ID == i, ]
        
        coords <- geom(init_points)[i, c("x", "y")]
        
        data.frame(
          StartYear  = start_year,
          Rep        = rep,
          InitSite   = i,
          X          = coords[1],
          Y          = coords[2],
          
          Code1 = sum(freq_i$n_cells[freq_i$habcode == 1]),
          Code2 = sum(freq_i$n_cells[freq_i$habcode == 2]),
          Code3 = sum(freq_i$n_cells[freq_i$habcode == 3]),
          Code4 = sum(freq_i$n_cells[freq_i$habcode == 4]),
          Code5 = sum(freq_i$n_cells[freq_i$habcode == 5]),
          Code6 = sum(freq_i$n_cells[freq_i$habcode == 6]),
          
          TotalCells = sum(freq_i$n_cells)
        )
      })
      
      all_summary <- rbind(all_summary, do.call(rbind, summary_rows))
      
      # ---- Individuals file (UNCHANGED) ----
      
      pop_size <- 2
      spruce_raster_init_pnt <- as.points(spruce_raster_init, values = TRUE, na.rm = TRUE)
      
      create_inds_dataframe <- function(shp_species, rst_species) {
        df_inds <- extract(rst_species, shp_species, cells = TRUE, df = TRUE)
        df_inds$Year <- 0
        df_inds$Species <- 0
        df_inds$X <- df_inds$cell %% ncol(rst_species) - 1
        df_inds$Y <- nrow(rst_species) - (floor(df_inds$cell / ncol(rst_species))) - 1
        df_inds$Ninds <- df_inds$pop
        df_inds <- df_inds[, !(names(df_inds) %in% c('ID', 'cells', 'pop'))]
        df_inds <- df_inds[!is.na(df_inds$Ninds), ]
        df_inds <- unique(df_inds)
        df_inds <- df_inds[df_inds$Ninds != 0, ]
        df_inds
      }
      
      df_f <- create_inds_dataframe(spruce_raster_init_pnt, spruce_raster_init)
      df_f$Ninds <- df_f$Ninds * pop_size / 2
      df_f$Sex <- 1
      
      df_m <- create_inds_dataframe(spruce_raster_init_pnt, spruce_raster_init)
      df_m$Ninds <- df_m$Ninds * pop_size / 2
      df_m$Sex <- 0
      
      df_individuals <- rbind(df_f, df_m)[, c("Year","Species","X","Y","Ninds","Sex")]
      
      inds_file <- paste0("inds_", start_year, "_rep", rep, ".txt")
      write.table(
        df_individuals,
        file.path(Output_path, inds_file),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
    }
  }
  
  write.csv(all_summary, paste0(Output_path, "init_summary_all.csv"), row.names = FALSE)
  
  return(all_summary)
}

# 
# species_init_map <- function(Output_path, NumberCells, start_years, reps, buffer_dist) {
#   all_summary <- data.frame()
#   
#   for (start_year in start_years) {
#     # Load quality and habitat maps for this start year
#     quality_file <- paste0(Output_path, start_year, "_quality_10yWindFell.tif")
#     habitat_file <- paste0(Output_path, start_year, "_habcodes_fell.asc")
#     
#     if (!file.exists(quality_file) || !file.exists(habitat_file)) {
#       warning("Missing files for year: ", start_year)
#       next
#     }
#     
#     Map_rast <- rast(quality_file)
#     habitat_rast <- rast(habitat_file)
#     
#     # Align habitat raster if needed
#     if (!compareGeom(habitat_rast, Map_rast, stopOnError = FALSE)) {
#       habitat_rast <- project(habitat_rast, Map_rast)
#     }
#     
#     for (rep in 1:reps) {
#       # Select valid cells (>0 and not NA in habitat)
#       valid_cells <- which(Map_rast[] > 0 & !is.na(habitat_rast[]))
#       if (length(valid_cells) < NumberCells) {
#         stop("Not enough cells with value > 0 to select.")
#       }
#       
#       selected_indices <- sample(valid_cells, NumberCells)
#       
#       # Create init raster
#       spruce_raster_init <- Map_rast
#       spruce_raster_init[] <- NA
#       spruce_raster_init[selected_indices] <- 1
#       names(spruce_raster_init) <- "pop"
#       
#       # Save initial distribution raster
#       raster_name <- paste0("initial_dist_", start_year, "_rep", rep, ".tif")
#       terra::writeRaster(spruce_raster_init, file.path(Output_path, raster_name), overwrite = TRUE, datatype = "INT2S")
#       
#       # Buffer summary
#       init_points <- as.points(spruce_raster_init, values = TRUE, na.rm = TRUE)
#       buffers <- terra::buffer(init_points, width = buffer_dist)
#       
#       summary_list <- lapply(1:nrow(buffers), function(i) {
#         cells <- terra::extract(habitat_rast, buffers[i, ], cells = TRUE)
#         codes <- cells[[1]][!is.na(cells[[1]])]
#         coords <- terra::geom(init_points)[i, c("x", "y")]
#         
#         data.frame(
#           StartYear = start_year,
#           Rep = rep,
#           InitSite = i,
#           X = coords[1],
#           Y = coords[2],
#           Code1 = sum(codes == 1),
#           Code2 = sum(codes == 2),
#           Code3 = sum(codes == 3),
#           Code4 = sum(codes == 4),
#           Code5 = sum(codes == 5),
#           Code6 = sum(codes == 6),
#           TotalCells = length(codes)
#         )
#       })
#       
#       all_summary <- rbind(all_summary, do.call(rbind, summary_list))
#       
#       # Individuals file
#       pop_size <- 2
#       spruce_raster_init_pnt <- as.points(spruce_raster_init, values = TRUE, na.rm = TRUE)
#       
#       create_inds_dataframe <- function(shp_species, rst_species) {
#         df_inds <- terra::extract(rst_species, shp_species, cells = TRUE, df = TRUE)
#         df_inds$Year <- 0
#         df_inds$Species <- 0
#         df_inds$X <- df_inds$cell %% ncol(rst_species) - 1
#         df_inds$Y <- nrow(rst_species) - (floor(df_inds$cell / ncol(rst_species))) - 1
#         df_inds$Ninds <- df_inds$pop
#         df_inds <- df_inds[, !(names(df_inds) %in% c('ID', 'cells', 'pop'))]
#         df_inds <- df_inds[!is.na(df_inds$Ninds), ]
#         df_inds <- unique(df_inds)
#         df_inds <- df_inds[df_inds$Ninds != 0, ]
#         return(df_inds)
#       }
#       
#       df_individuals_females <- create_inds_dataframe(spruce_raster_init_pnt, spruce_raster_init)
#       df_individuals_females$Ninds <- df_individuals_females$Ninds * pop_size / 2
#       df_individuals_females$Sex <- 1
#       df_individuals_females <- df_individuals_females[, c("Year", "Species", "X", "Y", "Ninds", "Sex")]
#       
#       df_individuals_males <- create_inds_dataframe(spruce_raster_init_pnt, spruce_raster_init)
#       df_individuals_males$Ninds <- df_individuals_males$Ninds * pop_size / 2
#       df_individuals_males$Sex <- 0
#       df_individuals_males <- df_individuals_males[, c("Year", "Species", "X", "Y", "Ninds", "Sex")]
#       
#       df_individuals <- rbind(df_individuals_females, df_individuals_males)
#       inds_file <- paste0("inds_", start_year, "_rep", rep, ".txt")
#       write.table(df_individuals, file.path(Output_path, inds_file), row.names = FALSE, quote = FALSE, sep = '\t')
#     }
#   }
#   
#   # Save combined summary
#   summary_file <- paste0(Output_path, "init_summary_all.csv")
#   write.csv(all_summary, summary_file, row.names = FALSE)
#   
#   return(all_summary)
# }

species_init_map(Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoD1/", 
                 start_years = c(2020, 2030, 2040), 
                 NumberCells = 1, 
                 reps = 100, 
                 buffer_dist = 500)
species_init_map(Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/NYM/", 
                 start_years = c(2020, 2030, 2040), 
                 NumberCells = 1, 
                 reps = 100, 
                 buffer_dist = 500)
species_init_map(Output_path = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/Kielder/", 
                 start_years = c(2020, 2030, 2040), 
                 NumberCells = 1, 
                 reps = 100, 
                 buffer_dist = 500)






















# felling data priocessing 

library(sf)
library(terra)
# 1. Read vector data
felled <- st_read(
  "Z:/CESB/Land Use and Ecosystem Service/LUES_Sware/PersonalFolders/Joe/Data/SCDB_20250331/EoY2425_FE_LMP.gdb",
  layer = "LMP_MANAGEMENT_COUPES",
  quiet = TRUE
)

# 2. Build template raster (extent, CRS, resolution already known)
template <- rast(
  ext = ext(felled),
  res = 25,
  crs = "EPSG:27700"
)

# 3. Convert sf → SpatVector and project to template CRS
v <- vect(felled)
v <- project(v, template)

# 4. Rasterize chosen field (maxmaifellyr)
r_maxmai <- rasterize(
  v,
  template,
  field = "maxmaifellyr",
  background = NA
)
writeRaster(r_maxmai, "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/felled_rast_EoR.tif")



# specific sample area




side_m <- 17150 # 295 km².


library(terra)
make_sample_square <- function(r, side_m, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  e <- ext(r)
  
  # allowable range for square origin
  xmin_max <- e$xmax - side_m
  ymin_max <- e$ymax - side_m
  
  x0 <- runif(1, e$xmin, xmin_max)
  y0 <- runif(1, e$ymin, ymin_max)
  
  ext_square <- ext(x0, x0 + side_m, y0, y0 + side_m)
  
  # as.polygons() already returns a SpatVector
  sq <- as.polygons(ext_square, crs = crs(r))
  
  return(sq)
}

# use this side_m for all landscapes

# paths
kielder_tif <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/IpsShift_Kielder.tif"
nym_tif     <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/IpsShift_NYM.tif"
fod_tif     <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/IpsShift_FoD.tif"

r_kielder <- rast(kielder_tif)
r_nym     <- rast(nym_tif)
r_fod     <- rast(fod_tif)

sq_kielder <- make_sample_square(r_kielder, side_m, seed = 1)
sq_nym     <- make_sample_square(r_nym,     side_m, seed = 1)
sq_fod     <- make_sample_square(r_fod,     side_m, seed = 1)

# crop rasters to the square
r_kielder_crop <- crop(r_kielder, sq_kielder)
r_nym_crop     <- crop(r_nym,     sq_nym)
r_fod_crop     <- crop(r_fod,     sq_fod)
# > r_kielder
# SpatExtent : 352643.5, 369793.5, 578384.812999999, 595534.812999999 (xmin, xmax, ymin, ymax)
# > r_nym <- ext(r_nym_crop)
# > r_nym
# SpatExtent : 476175, 493325, 484860.344000001, 502010.344000001 (xmin, xmax, ymin, ymax)
# > r_fod <- ext(r_fod_crop)
# > r_fod
# SpatExtent : 352888.1876, 370038.1876, 205389.812999999, 222539.812999999 (xmin, xmax, ymin, ymax)
# 
# 




# plotting code ####


# ------------------------------------------------------------
# Interactive map for habitat codes, incursion points & buffers
# ------------------------------------------------------------
# Packages
library(terra)
library(sf)
library(mapview)
library(htmlwidgets)

# -----------------------------
# Parameters (edit these)
# -----------------------------
Output_path <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoD_droughtonly/"
start_year  <- 2020
buffer_dist <- 500                                   # metres
save_html   <- TRUE
html_path   <- file.path(Output_path, paste0("interactive_map_", start_year, ".html"))

# -----------------------------
# 1) Load raster + init summary
# -----------------------------
habitat_file <- file.path(Output_path, paste0(start_year, "_habcodes_fell.asc"))
summary_file <- file.path(Output_path, "init_summary_all.csv")

stopifnot(file.exists(habitat_file), file.exists(summary_file))

hab_r  <- rast(habitat_file)
hab_crs <- crs(hab_r)

init_df <- read.csv(summary_file)
pts_df  <- subset(init_df, StartYear == start_year)

if (nrow(pts_df) == 0) stop("No points for StartYear = ", start_year)

# -----------------------------
# 2) Points to sf, create buffers (meters)
# -----------------------------
# If raster CRS exists, use it for points; otherwise assume BNG (EPSG:27700).
if (!is.na(hab_crs)) {
  pts_sf <- st_as_sf(pts_df, coords = c("X","Y"), crs = hab_crs)
} else {
  warning("hab_r has no CRS; assuming EPSG:27700 for points.")
  pts_sf <- st_as_sf(pts_df, coords = c("X","Y"), crs = 27700)
}

# If points are geographic (degrees), buffer in a metric CRS then come back
if (st_is_longlat(pts_sf)) {
  pts_metric <- st_transform(pts_sf, 27700)
  buffers_sf_m <- st_buffer(pts_metric, buffer_dist)
  buffers_sf <- if (!is.na(hab_crs)) st_transform(buffers_sf_m, st_crs(hab_crs)) else buffers_sf_m
  pts_sf <- if (!is.na(hab_crs)) st_transform(pts_sf, st_crs(hab_crs)) else pts_sf
} else {
  buffers_sf <- st_buffer(pts_sf, buffer_dist)
}

# Ensure both vectors match raster CRS (if defined)
if (!is.na(hab_crs)) {
  pts_sf     <- st_transform(pts_sf,     st_crs(hab_crs))
  buffers_sf <- st_transform(buffers_sf, st_crs(hab_crs))
}

# -----------------------------
# 3) Treat habcodes as categorical + labels
# -----------------------------
hab_f <- as.factor(hab_r)

# Use your known schema
levels(hab_f) <- data.frame(
  ID    = c(1, 2, 3, 4, 5, 6),
  class = c(
    "Living woodland (0%)",
    "Windblown (100%)",
    "Felled (10%)",
    "Drought risk (modelled)",
    "High wind damage risk (100%)",
    "Rotten"
  )
)

# Palette: adjust as you like
pal <- c(
  "Living woodland (0%)"         = "#2B8C3E",
  "Windblown (100%)"             = "#D7301F",
  "Felled (10%)"                 = "#BDBDBD",
  "Drought risk (modelled)"      = "#FDAE61",
  "High wind damage risk (100%)" = "#7F0000",
  "Rotten"                       = "#8C510A"
)

# -----------------------------
# 4) Build interactive map
# -----------------------------
m <- mapview(
  hab_f,
  col.regions = pal,
  alpha = 0.8,
  legend = TRUE
) +
  mapview(
    buffers_sf,
    color = "black",
    lwd = 1.2,
    alpha.regions = 0
  ) +
  mapview(
    pts_sf,
    col.regions = "cyan",
    cex = 2
  )

m  # show in viewer

# -----------------------------
# 5) Optional: save to HTML
# -----------------------------
if (save_html) {
  leaf <- m@map
  saveWidget(leaf, file = html_path, selfcontained = TRUE)
  message("Saved map to: ", html_path)
}

# -----------------------------
# 6) Quick diagnostics (optional)
# -----------------------------
# Check that we have multiple classes present
print(freq(hab_r))


