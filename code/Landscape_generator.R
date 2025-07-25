
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




# Main wrapper function to generate habitat suitability maps over time#function that creates landscape

# ===============================# Main Function: Landscape Time Series Generator# ===============================
create_time_series_maps <- function(Spruce_path,#path to spruce data
                                    TWI_path,#path to TWI
                                    Subcompart_path,#path to subcompartment data
                                    ForestGales_path,#path to forestgales data
                                    SPEI_data,# path to SPEI data summary by year
                                    Output_path,# path to store outputs
                                    extent,#extent of landscape
                                    crs,#crs, recoment BNG "EPSG:27700"
                                    TWI_min, #min cut off for TWI
                                    TWI_max, #max cut off for TWI
                                    years, #which years to run model
                                    Wind_prob #wind probably cut of for habitat suitability
) {
  # this function will take .geojson or .tif file and rasterize and clip to specified extent and then convert to specified crs 

# --- Rasterization and Clipping Utility ---
  Clip_and_rasterize <- function(Map_input, landuse, extent, crs) {
    library(terra)
    library(sf)
    
    clipped_extent <- terra::ext(extent) # Forest of Dean
    r <-     terra::rast(res = 25, ext = clipped_extent, crs = crs)
    
    if_statement_geo <-
      grepl(".geo", Map_input) # what is the file type
    if_statement_shp <-
      grepl(".shp", Map_input) # what is the file type
    
    if (if_statement_geo == TRUE | if_statement_shp == TRUE) {
      OGR_Data <-
        st_read(dsn = Map_input) # read in the forest gales data
      Vect_Data <- vect(OGR_Data)
      col_exist <- any(grepl('PRILANDUSE', names(Vect_Data)))
      col_exist_FG <- any(grepl('prob_damage.be', names(Vect_Data)))
      if (col_exist == TRUE) {
        Vect_Data_subset <-
          subset(Vect_Data, Vect_Data$PRILANDUSE == landuse)
        Map_rast <-
          terra::rasterize(Vect_Data_subset, r, background = NA)
        
      } else if (col_exist_FG == TRUE) {
        Map_rast <-
          terra::rasterize(Vect_Data, r, field = landuse, background = NA)
        
      } else {
        Map_rast <-
          terra::rasterize(Vect_Data, r, background = NA) # rasterize the forestGALES data
      }
      
    } else{
      Map_rast <- rast(Map_input) # read in raster
      Map_rast <- resample(Map_rast, r) # make same extent
      
      
    }
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
      Spruce <- Clip_and_rasterize(Spruce_path, extent, crs, landuse = NA)
      TWI <- Clip_and_rasterize(TWI_path, extent, crs, landuse = NA)
      
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
      
      Spruce[TWI[] == 100] <- 4 # these are drought risk areas
      Spruce[sub_compart_wind[] == 1] <-
        2 # cells that are windblow are 2
      Spruce[sub_compart_felled[] == 1] <-
        3 # cells that are felled are 3
      
      
      
      return(Spruce)
    }
  
# Integrate wind damage predictions into baseline map  # now add the predicted yearly wind damage to the baseline map and save each of these

# --- Annual Wind Damage Integration ---
  Add_yearly_wind_damage <-
    function(years,
             ForestGales_path,
             Output_path,
             Baseline_map,
             Wind_prob) {
      map_extent <- terra::ext(Baseline_map)
      crs <- crs(Baseline_map)
      
      for (y in years) {
        #create copy of spruce raster to play with
        field_choice = paste("prob_damage.be", y, sep = "_")
        ForestGales <-
          Clip_and_rasterize(ForestGales_path, map_extent, crs, landuse = field_choice)
        
        ForestGales[ForestGales[] < Wind_prob] <- NA
        ForestGales[ForestGales[] > Wind_prob] <-
          5 # high wind damage chance
        ForestGales[is.na(Baseline_map)] <- NA # but NA if its not spruce
        Baseline_map_y <- Baseline_map
        Baseline_map_y[ForestGales[] == 5] <-
          5  # in spruce_raster if its a 5 its windblown
        
        ## now lets get previous years map  - uness this is the first year
        if (y < min(years)) {
          filename_previous <-
            paste0(Output_path, y - 1, "habcodes.asc", sep = "")
          filename_previous_tif <-
            paste0(Output_path, y - 1, "habcodes.tif", sep = "")
          last_year <- rast(filename_previous)
          Baseline_map_y[last_year[] == 5] <- 5 #windblow from the last year is now added
          plot(last_year)
          
          
          writeRaster(last_year, filename_previous_tif, overwrite = TRUE)
          
          writeRaster(
            last_year,
            filename,
            overwrite = TRUE,
            # datatype = 'INT1U',
            NAflag = -999,
            datatype = "INT2S"
          )
          
        } else{
          filename <- paste0(Output_path, y, "habcodes.asc", sep = "")
          filename_tif <- paste0(Output_path, y, "habcodes.tif", sep = "")
          
          
          writeRaster(Baseline_map_y, filename_tif, overwrite = TRUE)
          
          writeRaster(
            Baseline_map_y,
            filename,
            overwrite = TRUE,
            #  datatype = 'INT1U',
            NAflag = -999,
            datatype = "INT2S"
          )
        }
        
        plot(Baseline_map_y)
      } # loops through every year specified and then add storm damage (cumulative)
    }
  
  
  
# Convert habitat codes to quality scores based on wind and drought  #This duntion will convert the habitat codes to quality
  
  # requires average SPEI projection and year for area

# --- Habitat Quality Conversion: Wind + Drought ---
  HabitatCodetoIpsQual_both <- function(years, SPEI_data, Output_path) {
    for (y in years) {
      filename <-
        paste0(Output_path, y, "habcodes.asc", sep = "")
      
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
      
      SPEI <- seq(from = -4, to = 0, by = 0.5)
      qual <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)
      model <- lm(qual ~ bs(SPEI, knots = c(12000, 60000, 150000))) # model to fit line through data
      
      SPEI_local <- read.csv(SPEI_data) #year;y
      SPEI_local$pred  <- predict(model, newdata = list(SPEI = SPEI_local$SPEI)) # now predict for our SPEI
      SPEI_local$pred  <- round (SPEI_local$pred, digits = 0)
      SPEI_local$pred[SPEI_local$pred < 0] <- 0 #if SPEI is less than 0 then no qual associated
      SPEI_local$pred[SPEI_local$pred == "NaN"] <- 0 # NA values are also 0
      SPEI_value <- subset(SPEI_local, SPEI_local$year == y)
      SPEI_value <- SPEI_value$pred
      
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
        paste0(Output_path, y, "habcodes.asc", sep = "")
      library(splines)
      library(Ecdat)
      library(dplyr)
      
      yearone <- min(years)
      Map_rast <- rast(filename)
      Map_rast[Map_rast[] == 1] <- 0 # LIVING wont sustain a pop
      Map_rast[Map_rast[] == 2] <- 100 # wind blown
      Map_rast[Map_rast[] == 3] <- 10 # felled
      Map_rast[Map_rast[] == 5] <- 100 # high chance of wind damage
      
      
      SPEI <- seq(from = -4, to = 0, by = 0.5)
      qual <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)
      model <- lm(qual ~ bs(SPEI, knots = c(12000, 60000, 150000))) # model to fit line through data
      
      SPEI_local <- read.csv(SPEI_data) #year;y
      SPEI_local$pred  <- predict(model, newdata = list(SPEI = SPEI_local$SPEI)) # now predict for our SPEI
      SPEI_local$pred  <- round (SPEI_local$pred, digits = 0)
      SPEI_local$pred[SPEI_local$pred < 0] <- 0 #if SPEI is less than 0 then no qual associated
      SPEI_local$pred[SPEI_local$pred == "NaN"] <- 0 # NA values are also 0
      SPEI_value <- subset(SPEI_local, SPEI_local$year == yearone)
      SPEI_value <- SPEI_value$pred
      
      Map_rast[Map_rast[] == 4] <- SPEI_value # high chance of wind damage
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
        paste0(Output_path, yearone, "habcodes.asc", sep = "")
      
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
      
      
      for (y in years) {
        
        Map_rast_y <- Map_rast
      SPEI <- seq(from = -4, to = 0, by = 0.5)
      qual <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)
      model <- lm(qual ~ bs(SPEI, knots = c(12000, 60000, 150000))) # model to fit line through data
      
      SPEI_local <- read.csv(SPEI_data) #year;y
      SPEI_local$pred  <- predict(model, newdata = list(SPEI = SPEI_local$SPEI)) # now predict for our SPEI
      SPEI_local$pred  <- round (SPEI_local$pred, digits = 0)
      SPEI_local$pred[SPEI_local$pred < 0] <- 0 #if SPEI is less than 0 then no qual associated
      SPEI_local$pred[SPEI_local$pred == "NaN"] <- 0 # NA values are also 0
      SPEI_value <- subset(SPEI_local, SPEI_local$year == y)
      SPEI_value <- SPEI_value$pred
      
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
  
  Add_yearly_wind_damage(years,
                         ForestGales_path,
                         Output_path,
                         Baseline,
                         Wind_prob)
  
  HabitatCodetoIpsQual_both(years, SPEI_data, Output_path)
  HabitatCodetoIpsQual_drought(years, SPEI_data, Output_path)
  HabitatCodetoIpsQual_wind(years, SPEI_data, Output_path)
  
}
 

create_time_series_maps(
  Spruce_path =  "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/Nicks model/data/raw/Spruce_Map/Spruce_Map.shp",
  TWI_path = "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/Ips_dynamic/data/TWI_FionaP_crop.tif",
  Subcompart_path = "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/FOD_Feb24/National_Forest_Estate_Subcompartments_England_2019/National_Forest_Estate_Subcompartments_England_2019.shp",
  ForestGales_path = "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/Ips_dynamic/data/June_FORESTGALES/2020-2035_outputs/dean_fg_outputs_2020-2035_Brown_Edge.geojson",
  SPEI_data = "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/PHENIPS/SPEI2.csv",
  Output_path = "E:/Ips_dynamic/data/landscapes/",
  extent = c(350000, 375000, 200000, 225000),
  crs =  "EPSG:27700",
  TWI_min = 2,
  TWI_max = 15,
  years = 2020:2035,
  Wind_prob = 0.9
)

# Load year-one quality map and select initial population cells#reads in baseline map from year one and selects starting population 


# --- Species Initialization Map Generator ---
species_init_map <- function(Output_path, #path to time series maps
                             NumberCells, #number opf starting populations
                             rep, #replicate number for model
                             years #years of model
                             )  {
  first_year <- min(years)
  filename <-
    paste0(Output_path, first_year, "_quality.asc", sep = "")
  
  Map_rast <- rast(filename)
  
  high_value_cells <- Map_rast[Map_rast[] == 100] # list of cells with value greater than cutoff
  selected_indices <- sample(length(high_value_cells[high_value_cells[] == 100]), NumberCells)
  high_value_cells[selected_indices, ] <- 1000
  
  
  spruce_raster_init <- Map_rast
  spruce_raster_init[spruce_raster_init[] == 100] <- high_value_cells
  spruce_raster_init[spruce_raster_init[] < 1000] <- NA
  spruce_raster_init[spruce_raster_init[] == 1000] <- 1
  
  # Convert species distribution point to ascii raster for rangeshifter
  names(spruce_raster_init) <- c("pop")
  terra::writeRaster(
    spruce_raster_init,
    file.path(Output_path, "initial_dist.tif"),
    overwrite = TRUE,
    datatype = "INT2S"
  )
  
  # Write a individuals file for initialisation
  pop_size <- 1000
  spruce_raster_init_pnt <- as.points(spruce_raster_init, values = TRUE, na.rm =
                                        TRUE)
  
  create_inds_dataframe <- function(shp_species, rst_species) {
    df_inds <- terra::extract(rst_species, shp_species, cells = T, df = TRUE)
    
    df_inds$Year <- 0
    df_inds$Species <- 0
    df_inds$X <- df_inds$cell %% ncol(rst_species) - 1
    df_inds$Y <- nrow(rst_species) - (floor(df_inds$cell / ncol(rst_species))) - 1
    df_inds$Ninds <- df_inds$pop
    df_inds <- df_inds[, !(names(df_inds) %in% c('ID', 'cells', 'pop'))]
    df_inds <- df_inds[!is.na(df_inds$Ninds), ]
    
    # drop duplicates
    df_inds <- unique(df_inds)
    # add another catch where Ninds == 0 (remove)
    if (nrow(df_inds[which(df_inds$Ninds == 0), ]) > 0) {
      df_inds <- df_inds[-which(df_inds$Ninds == 0), ]
    }
    
    return(df_inds)
    
  }
  
  df_individuals_females <- create_inds_dataframe(shp_species = spruce_raster_init_pnt, rst_species = spruce_raster_init)
  df_individuals_females$Ninds <- df_individuals_females$Ninds * pop_size / 2
  df_individuals_females$Sex <- 1
  df_individuals_females <- df_individuals_females[, c("Year", "Species", "X", "Y", "Ninds", "Sex")]
  
  df_individuals_males <- create_inds_dataframe(shp_species = spruce_raster_init_pnt, rst_species = spruce_raster_init)
  df_individuals_males$Ninds <- df_individuals_males$Ninds * pop_size / 2
  df_individuals_males$Sex <- 0
  df_individuals_males <- df_individuals_males[, c("Year", "Species", "X", "Y", "Ninds", "Sex")]
  
  df_individuals <- rbind(df_individuals_females, df_individuals_males)
  file <- paste(NumberCells, rep, "inds.txt", sep = "_")
  write.table(
    df_individuals,
    file.path(Output_path, file),
    row.names = F,
    quote = F,
    sep = '\t'
  )
  
  
  return(spruce_raster_init)
}

species_init_map(Output_path = "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/Ips_dynamic/data/forestsGales/", 
                 years = 2020:2035, 
                 NumberCells = 1, 
                 rep = 1)




