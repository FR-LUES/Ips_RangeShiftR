
####### LCM
# 1	Broadleaved woodland
# 2	Coniferous woodland
# 3	Arable
# 4	Improved grassland
# 5	Neutral grassland
# 6	Calcareous grassland
# 7	Acid grassland
# 8	Fen, Marsh and Swamp
# 9	Heather and shrub
# 10	Heather grassland
# 11	Bog
# 12	Inland rock
# 13	Saltwater
# 14	Freshwater
# 15	Supralittoral rock
# 16	Supralittoral sediment
# 17	Littoral rock
# 18	Littoral sediment
# 19	Saltmarsh
# 20	Urban
# 21	Suburban

####### NFI
# 1	Agriculture land
# 2	Assumed woodland
# 3	Bare area
# 4	Broadleaved
# 5	Conifer <<<<<<<<<<
# 6	Coppice
# 7	Coppice with standards
# 8	Failed
# 9	Felled
# 10	Grassland
# 11	Ground prep
# 12	Low density
# 13	Mixed mainly broadleaved
# 14	Mixed mainly conifer <<<<<<<
# 15	Open water
# 16	Other vegetation
# 17	Powerline
# 18	Quarry
# 19	River
# 20	Road
# 21	Shrub
# 22	Uncertain
# 23	Urban
# 24	Windblow
# 25	Young trees


####### Spruce Map
# 1 Spruce



library(raster)
library(rgdal)
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
  
  extent_str <- sprintf("xmin%s-xmax%s-ymin%s-ymax%s", clipped_extent@xmin, clipped_extent@xmax, clipped_extent@ymin, clipped_extent@ymax)
  dir.create(file.path(dirpath, extent_str), recursive=T)
  
  shp_species <- read_sf(dsn = species_shp_path)
  
  rst_lcm <- raster(lcm_rst_path)
  rst_nfi <- raster(nfi_rst_path)
  rst_spruce <- raster(spruce_rst_path)
  
  # all input rasters should be the same resolution
  stopifnot(res(rst_lcm) == res(rst_nfi), res(rst_lcm) == res(rst_spruce))
  original_res <- res(rst_lcm)[1]
  new_res <- original_res * aggregation_factor
  
  # Merge the habitat maps into single ascii raster for rangeshifter
  rst_lcm <- crop(rst_lcm, clipped_extent)
  rst_nfi <- crop(rst_nfi, clipped_extent)
  rst_spruce <- crop(rst_spruce, clipped_extent)
  
  if (aggregation_factor > 1) {
    rst_lcm <- aggregate(rst_lcm, fact=aggregation_factor, fun=modal)
    rst_nfi <- aggregate(rst_nfi, fact=aggregation_factor, fun=modal)
    rst_spruce_hab <- aggregate(rst_spruce, fact=aggregation_factor, fun=max)
    # take sum here because spruce=1, so the sum is the number of 25m pixels occurring in
    # the aggregated pixel
    rst_hab_quality <- aggregate(rst_spruce, fact=aggregation_factor, fun=sum)
  }
  
  rst_lcm[rst_nfi == 5]<- 22 # conifer
  rst_lcm[rst_nfi == 14] <- 23 # mixed mainly conifer
  rst_lcm[rst_spruce_hab == 1] <- 24 # spruce
  
  # rst_lcm <- reclassify(rst_lcm, cbind(NA, -9999)) 
  writeRaster(rst_lcm, file.path(dirpath, extent_str, "hab_code.tif"), overwrite=TRUE)
  
  # convert quantity of spruce pixels to a percentage per aggregated pixel
  rst_hab_quality <- (rst_hab_quality / aggregation_factor**2) * 100
  rst_hab_quality <- as.integer(rst_hab_quality)
  rst_hab_quality[rst_lcm < 24] <- 0
  writeRaster(rst_hab_quality, file.path(dirpath, extent_str, "hab_quality.tif"), overwrite=TRUE)
  
  # Convert species distribution point to ascii raster for rangeshifter
  # rst_species <- reclassify(rst_lcm, cbind(0, 25, NA))
  rst_species <- reclassify(rst_lcm, cbind(0, 25, 0))
  rst_species <- rasterize(shp_species, rst_species, field="pop", update=TRUE)
  writeRaster(rst_species, file.path(dirpath, extent_str, "initial_dist.tif"), overwrite=TRUE)
  
  # Write a individuals file for initialisation
  pop_size <- 1000
  df_individuals <- create_inds_dataframe(shp_species, rst_species)
  df_individuals$Ninds <- df_individuals$Ninds * pop_size
  write.table(df_individuals,  file.path(dirpath, extent_str, "inds.txt"), row.names = F, quote = F, sep = '\t')
  
  extent_lcm <- extent(rst_lcm)
  # extent(rst_lcm) <- c(0, extent(rst_lcm)[2]-extent(rst_lcm)[1], 0, extent(rst_lcm)[4]-extent(rst_lcm)[3])
  # extent(rst_species) <- c(0, extent(rst_species)[2]-extent(rst_species)[1], 0, extent(rst_species)[4]-extent(rst_species)[3])
  writeRaster(rst_lcm, file.path(dirpath, extent_str, "hab_code.asc"), format="ascii", overwrite=TRUE)
  writeRaster(rst_hab_quality, file.path(dirpath, extent_str, "hab_quality.asc"), format="ascii", overwrite=TRUE)
  writeRaster(rst_species, file.path(dirpath, extent_str, "initial_dist.asc"), format="ascii", overwrite=TRUE)
  
}


create_inds_dataframe <- function(
    shp_species,
    rst_species
) {
  
  
  df_inds <- raster::extract(rasterize(shp_species, rst_species, fun="last", field='pop'), shp_species, cellnumbers=T, df=TRUE)
  
  df_inds$Year <- 0
  df_inds$Species <- 0
  df_inds$X <- df_inds$cells %% ncol(rst_species) - 1
  df_inds$Y <- nrow(rst_species) - (floor(df_inds$cells / ncol(rst_species))) - 1
  df_inds$Ninds <- df_inds$layer
  df_inds <- df_inds[ , !(names(df_inds) %in% c('ID', 'cells', 'layer'))]
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
  df_pop$x_coord <- (df_pop$x + 0.5) * sim_params@land@Resolution + geom_extent@xmin
  df_pop$y_coord <- (df_pop$y + 0.5) * sim_params@land@Resolution + geom_extent@ymin
  if (nrow(df_pop) > 0) {
    df_pop <- st_as_sf(df_pop, coords = c("x_coord","y_coord"))
    df_pop <- df_pop %>% st_set_crs(st_crs(coord_ref_sys))
  }
  
  return(df_pop)
  
}




apply_pop_impact_on_hab_quality <- function(
    pop_rst,
    years_for_destruction,
    apply_destruction_annually,
    sim_params,
    coord_ref_sys,
    geom_extent,
    dirpath
  ) {
  
  hab_qual_rst  <- raster(file.path(dirpath, "Inputs", sim_params@land@LandscapeFile))
  crs(hab_qual_rst) <- coord_ref_sys
  extent(hab_qual_rst) <- geom_extent
  
  if (apply_destruction_annually) {
    quality_impact <- ((pop_rst / sim_params@land@K_or_DensDep) * 100)
  } else {
    quality_impact <- ((pop_rst / sim_params@land@K_or_DensDep) * 100) / years_for_destruction
  }
  
  hab_qual_rst <- as.integer(hab_qual_rst - quality_impact)
  hab_qual_rst[hab_qual_rst < 0] <- 0
  
  return(hab_qual_rst)
  
}




extract_pop_shp_to_table <- function(
    pop_shp
) {
  
  pop_shp <- pop_shp[, c("Year", "Species", "x", "y", "NInd")]
  
  names(pop_shp)[names(pop_shp) == "x"] <- "X"
  names(pop_shp)[names(pop_shp) == "y"] <- "Y"
  names(pop_shp)[names(pop_shp) == "NInd"] <- "Ninds"
  
  pop_shp$Year <- 0
  
  pop_shp <- pop_shp %>% st_drop_geometry()
  
  return(pop_shp)
}




extract_pop_shp_to_rst <- function(
    pop_shp,
    rst_lcm
  ) {
  
  pop_rst <- reclassify(rst_lcm, cbind(minValue(rst_lcm), maxValue(rst_lcm), 0))
  if (sum(pop_shp$NInd) > 0) {
    pop_rst <- rasterize(pop_shp, pop_rst, field="NInd", update=TRUE)
  } 
  
  return(pop_rst)
  
}




PausingExtras <- function(
    replicates,
    years,
    pop_impacts_quality,
    years_for_destruction,
    apply_destruction_annually
) {
  
  setClass(
    "pausingextras",
    representation(
      replicates="numeric",
      years="numeric",
      pop_impacts_quality="logical",
      years_for_destruction="numeric",
      apply_destruction_annually="logical"
    )
  )
  
  pausingextras <- new(
    "pausingextras",
    replicates=replicates,
    years=years,
    pop_impacts_quality=pop_impacts_quality,
    years_for_destruction=years_for_destruction,
    apply_destruction_annually=apply_destruction_annually
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
  
  rst_lcm <- raster(file.path(dirpath, "Inputs", sprintf("%s.tif", habitat_rst_name)))
  
  start_init <- sim_params@init
  start_land <- sim_params@land
  
  for (replicate in 1:replicates) {
    
    sim_params@init <- start_init
    sim_params@land <- start_land
    
    out_rst_stack <- stack()
    
    for (year in 1:years) {
      
      # RangeShiftR needs a backslash at end of path
      result <- RunRS(sim_params, sprintf("%s/", dirpath))
      
      # Sleep for a second (not sure if needed?)
      Sys.sleep(1)
      
      pop_shp <- extract_pop_file_to_shp(
        sim_params,
        crs(rst_lcm),
        extent(rst_lcm),
        dirpath
      )
      st_write(pop_shp, file.path(dirpath, "Output_Maps", sprintf("rep%s_year%s_pop.geojson", replicate,  year)))
      
      pop_rst <- extract_pop_shp_to_rst(
        pop_shp,
        rst_lcm
      )
      writeRaster(pop_rst, file.path(dirpath, "Output_Maps", sprintf("rep%s_year%s_pop.tif", replicate,  year)))
      out_rst_stack <- addLayer(out_rst_stack, pop_rst)
      
      # If Population is zero, then break out of this simulation
      if (nrow(pop_shp) > 0) {
        
        df_pop_table <- extract_pop_shp_to_table(pop_shp)
        new_init_inds_txt <- sprintf('rep%s_year%s_inds.txt', replicate, year)
        write.table(df_pop_table, file.path(dirpath, "Inputs", new_init_inds_txt), row.names = F, quote = F, sep = '\t')
        sim_params@init@InitIndsFile <- new_init_inds_txt
        
        # Create new adjusted habitat quality raster
        if (sim_params@pausingextras@pop_impacts_quality) {
          if (
            (sim_params@pausingextras@apply_destruction_annually) |
            (year %% sim_params@pausingextras@years_for_destruction == 0)
          ) {
            
            new_hab_qual_rst <- apply_pop_impact_on_hab_quality(
              pop_rst=pop_rst,
              years_for_destruction=sim_params@pausingextras@years_for_destruction,
              apply_destruction_annually=sim_params@pausingextras@apply_destruction_annually,
              sim_params=sim_params,
              coord_ref_sys=crs(rst_lcm),
              geom_extent=extent(rst_lcm),
              dirpath=dirpath
            )
            
            new_landscape_asc <- sprintf("rep%s_year%s_%s.asc", replicate, year, habitat_rst_name)
            writeRaster(new_hab_qual_rst, file.path(dirpath, "Inputs", new_landscape_asc))
            sim_params@land@LandscapeFile <- new_landscape_asc
            
          }
        }
        
        for (out_file in c("Pop", "Range", "Parameters")) {
          file.rename(
            file.path(dirpath, "Outputs", sprintf("Batch0_Sim0_Land1_%s.txt", out_file)),
            file.path(dirpath, "Outputs", "file-per-pause", sprintf("rep%s_year%s_Batch0_Sim0_Land1_%s.txt", replicate, year, out_file))
          )
        }
        
      }
      
      # If Population is zero, then break out of this simulation
      if (nrow(pop_shp) == 0) break
      
    }
    
    writeRaster(out_rst_stack, file.path(dirpath, "Output_maps", sprintf("rep%s_pop_stack.tif", replicate)))
    
  }
  
  merge_output_pop_and_range(dirpath, replicates, years)
  
  
  jpeg(file=sprintf("%s/Outputs/%s", dirpath, "abundance-and-occupancy.jpg"))
  par(mfrow=c(1,2))
  plotAbundance(sim_params, sprintf("%s/", dirpath))
  plotOccupancy(sim_params, sprintf("%s/", dirpath))
  dev.off()
  
  jpeg(file=sprintf("%s/Outputs/%s", dirpath, "abundance-and-occupancy-SD.jpg"))
  par(mfrow=c(1,2))
  plotAbundance(sim_params, sprintf("%s/", dirpath), sd=T, replicates=F)
  plotOccupancy(sim_params, sprintf("%s/", dirpath), sd=T, replicates=F)
  dev.off()
  
}


merge_output_pop_and_range <- function(
    dirpath,
    replicates,
    years
) {
  pop_list <- list()
  range_list <- list()
  i <- 1
  for (replicate in 1:replicates) {
    for (year in 1:years) {
      
      pop_txt_path <- file.path(dirpath, "Outputs", "file-per-pause", sprintf("rep%s_year%s_Batch0_Sim0_Land1_Pop.txt", replicate, year))
      range_txt_path <- file.path(dirpath, "Outputs", "file-per-pause", sprintf("rep%s_year%s_Batch0_Sim0_Land1_Range.txt", replicate, year))
      
      if (file.exists(pop_txt_path) & file.exists(range_txt_path)) {
        
        df_pop <- read.table(pop_txt_path, sep="\t", header=T)
        df_pop <- subset(df_pop, Year == 0)
        df_pop$Year <- year-1
        df_pop$Rep <- replicate
        pop_list[[i]] <- df_pop
        
        df_range <- read.table(range_txt_path, sep="\t", header=T)
        df_range <- subset(df_range, Year == 0)
        df_range$Year <- year-1
        df_range$Rep <- replicate
        range_list[[i]] <- df_range
        
        i <- i + 1
        
      }
    }
  }
  df_pop = do.call(rbind, pop_list)
  df_range = do.call(rbind, range_list)
  write.table(df_pop,  file.path(dirpath, "Outputs", "Batch0_Sim0_Land1_Pop.txt"), row.names = F, quote = F, sep = '\t')
  write.table(df_range,  file.path(dirpath, "Outputs", "Batch0_Sim0_Land1_Range.txt"), row.names = F, quote = F, sep = '\t')
}



