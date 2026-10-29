# Databricks notebook source
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------
batch_dir <- "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/2020_batch20260306"
years_to_plot <- c(1:10)

# Use the first Pop.txt in the batch as an example
pop_file <- list.files(
  batch_dir,
  pattern = "Pop\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)[1]

cat("Using Pop file:\n", pop_file, "\n\n")

# ------------------------------------------------------------
# READ POP DATA
# ------------------------------------------------------------
df <- read_tsv(
  pop_file,
  col_types = cols(
    Rep       = col_double(),
    Year      = col_double(),
    RepSeason = col_double(),
    x         = col_double(),
    y         = col_double(),
    Species   = col_double(),
    NInd      = col_double(),
    Nfemales  = col_double(),
    Nmales    = col_double()
  )
)

# Keep occupied cells only
df_occ <- df %>% filter(NInd > 0)

# ------------------------------------------------------------
# INITIALISATION SITE (centroid of Year 0, averaged over reps)
# ------------------------------------------------------------
init_site <- df_occ %>%
  filter(Year == 0) %>%
  group_by(Rep) %>%
  summarise(x0 = mean(x), y0 = mean(y), .groups = "drop") %>%
  summarise(
    init_x = mean(x0),
    init_y = mean(y0)
  )

# ------------------------------------------------------------
# FIXED DENOMINATOR: total reps present in Year 1
# ------------------------------------------------------------
n_reps_total <- df %>%
  filter(Year == 1) %>%
  distinct(Rep) %>%
  nrow()

cat("Total replicates (Year 1):", n_reps_total, "\n\n")

# ------------------------------------------------------------
# CELL-LEVEL OCCUPANCY PROPORTIONS
# ------------------------------------------------------------
cell_props <- df_occ %>%
  filter(Year %in% years_to_plot) %>%
  distinct(Rep, Year, x, y) %>%     # occupancy only once per rep
  count(Year, x, y, name = "n_occ") %>%
  mutate(prop = n_occ / n_reps_total)

# ------------------------------------------------------------
# PLOT: ONE MAP PER YEAR
# ------------------------------------------------------------

ggplot(cell_props, aes(x = x, y = y, fill = prop)) +
  geom_tile() +
  geom_point(
    data = init_site,
    aes(x = init_x, y = init_y),
    inherit.aes = FALSE,
    shape = 4,
    size = 3,
    stroke = 1,
    colour = "red"
  ) +
  facet_wrap(~ Year) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    name = "Proportion\nof replicates"
  ) +
  coord_equal() +
  labs(
    title = "Occupancy proportion through time",
    subtitle = "Proportion of initial replicates occupying each cell",
    x = "x (grid)",
    y = "y (grid)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines")
  )

# COMMAND ----------

batch_dir <- "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/"
list.files(batch_dir)

# COMMAND ----------

install.packages("viridis")

library(tidyverse)
library(terra)
library(viridis)
library(ggplot2)

# ------------------------------------------------------------
# SETTINGS — choose your batch here
# ------------------------------------------------------------
batch_dir <- "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/2020_batch20260413/Batch1/Outputs"
pop_file  <- file.path(batch_dir, "Batch0_Sim0_Land1_Pop.txt")

raster_file <- "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/2020_quality_10yWindFell.asc"

# ------------------------------------------------------------
# READ RASTER
# ------------------------------------------------------------
r <- rast(raster_file)

r_df <- as.data.frame(r, xy = TRUE)
colnames(r_df) <- c("x_coord", "y_coord", "value")

# ------------------------------------------------------------
# READ POP FILE
# ------------------------------------------------------------
df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)

# Keep only occupied cells
df_occ <- df %>% filter(NInd > 0)

# ------------------------------------------------------------
# CONVERT GRID INDICES → RASTER COORDINATES
# ------------------------------------------------------------
xy <- xyFromCell(r, cellFromRowCol(r, df_occ$y, df_occ$x))
df_occ$x_coord <- xy[,1]
df_occ$y_coord <- xy[,2]

# ------------------------------------------------------------
# SUMMARISE PRESENCE PER YEAR
# ------------------------------------------------------------
presence_summary <- df_occ %>%
  group_by(x_coord, y_coord, Year) %>%
  summarise(Presence = n(), .groups = "drop")

# Choose years to plot
years_to_plot <- c(1, 2, 3)

presence_filtered <- presence_summary %>%
  filter(Year %in% years_to_plot)

# ------------------------------------------------------------
# PLOT RASTER + OCCUPANCY
# ------------------------------------------------------------
ggplot() +
  geom_raster(data = r_df, aes(x = x_coord, y = y_coord, fill = value)) +
  scale_fill_viridis(name = "Wind-fell\nquality") +
  
  geom_point(
    data = presence_filtered,
    aes(x = x_coord, y = y_coord, colour = Presence),
    size = 2
  ) +
  scale_colour_viridis(option = "plasma", name = "Presence") +
  
  facet_wrap(~Year) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Wind-fell raster with occupancy overlay",
    x = "X coordinate",
    y = "Y coordinate"
  )


# COMMAND ----------

library(tidyverse)
library(terra)
library(viridis)
library(ggplot2)
install.packages("ggnewscale")
# ------------------------------------------------------------
# SETTINGS — choose your batch here
# ------------------------------------------------------------
batch_dir <- "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/2020_batch20260413/Batch1/Outputs"
pop_file  <- file.path(batch_dir, "Batch0_Sim0_Land1_Pop.txt")

raster_file <- "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/2020_quality_10yWindFell.asc"

# ------------------------------------------------------------
# READ RASTER — FORCE SINGLE LAYER
# ------------------------------------------------------------
r <- rast(raster_file)[[1]]   # <— CRITICAL FIX

# Convert to binary suitability (1 = suitable, 0 = unsuitable)
suit <- classify(r, rcl = matrix(c(-Inf, 0, 0,
                                   0, Inf, 1), ncol = 3, byrow = TRUE))

suit_df <- as.data.frame(suit, xy = TRUE)
colnames(suit_df) <- c("x_coord", "y_coord", "suitable")

# ------------------------------------------------------------
# READ POP FILE
# ------------------------------------------------------------
df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)

df_occ <- df %>% filter(NInd > 0)

# ------------------------------------------------------------
# CONVERT GRID INDICES → RASTER COORDINATES (robust)
# ------------------------------------------------------------
cells <- cellFromRowCol(r, df_occ$y, df_occ$x)
cells <- as.numeric(cells)

xy <- xyFromCell(r, cells)

df_occ$x_coord <- xy[,1]
df_occ$y_coord <- xy[,2]

# ------------------------------------------------------------
# RASTERISE OCCUPANCY (heatmap)
# ------------------------------------------------------------
occ_rast <- r[[1]]      # ensure single layer
values(occ_rast) <- NA

coords <- as.matrix(df_occ[, c("x_coord", "y_coord")])
cell_ids <- cellFromXY(r, coords)

cell_ids <- as.numeric(cell_ids)

occ_counts <- table(cell_ids)



occ_counts <- table(cell_ids)

occ_rast[as.numeric(names(occ_counts))] <- as.numeric(occ_counts)

occ_df <- as.data.frame(occ_rast, xy = TRUE)
colnames(occ_df) <- c("x_coord", "y_coord", "occ")

# ------------------------------------------------------------
# PLOT SUITABILITY + OCCUPANCY HEATMAP
# ------------------------------------------------------------
suit_df <- suit_df %>%
  mutate(suitability_fill = ifelse(suitable == 1, "Suitable", "Unsuitable"))

occ_df <- occ_df %>%
  mutate(occ_fill = occ)


library(ggnewscale)

ggplot() +
  # ------------------------------------------------------------
  # 1. Suitability raster (discrete fill)
  # ------------------------------------------------------------
  geom_raster(
    data = suit_df,
    aes(x = x_coord, y = y_coord, fill = suitability_fill)
  ) +
  scale_fill_manual(
    name = "Suitability",
    values = c("Unsuitable" = "white", "Suitable" = "grey80")
  ) +

  # ------------------------------------------------------------
  # NEW FILL SCALE (this resets the fill aesthetic)
  # ------------------------------------------------------------
  new_scale_fill() +

  # ------------------------------------------------------------
  # 2. Occupancy raster (continuous fill)
  # ------------------------------------------------------------
  geom_raster(
    data = occ_df,
    aes(x = x_coord, y = y_coord, fill = occ_fill),
    alpha = 0.7
  ) +
  scale_fill_viridis(
    name = "Occupancy",
    option = "plasma",
    na.value = NA
  ) +

  coord_equal() +
  theme_minimal() +
  labs(
    title = "Binary Suitability Raster + Rasterised Occupancy Heatmap",
    x = "X coordinate",
    y = "Y coordinate"
  )



# COMMAND ----------

install.packages("ggnewscale")
install.packages("viridis")

  library(tidyverse)
  library(terra)
  library(ggplot2)
  library(viridis)
  library(ggnewscale)
plot_time_series <- function(
  site = "2020", 
  batch_id = "Batch20",
  years = 1:10,
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/"
) {
  
  library(tidyverse)
  library(terra)
  library(ggplot2)
  library(viridis)
  library(ggnewscale)
  
  # ------------------------------------------------------------
  # 1. Identify all date-stamped folders for this site
  #    and select ONLY the most recent one
  # ------------------------------------------------------------
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  
  if (length(site_dirs) == 0) stop("No site folders found.")
  
  site_dir <- sort(site_dirs, decreasing = TRUE)[1]   # <-- CRITICAL
  
  # ------------------------------------------------------------
  # 2. Build the path to the chosen batch inside this site
  # ------------------------------------------------------------
  batch_path <- file.path(site_dir, batch_id, "Outputs")
  if (!dir.exists(batch_path)) stop("Batch folder not found in latest site folder.")
  
  # ------------------------------------------------------------
  # 3. Find Pop file ONLY in the main Outputs folder
  # ------------------------------------------------------------
  pop_files <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(pop_files) == 0) {
    stop("No valid Pop files found (excluding file-per-pause).")
  }
  
  # ------------------------------------------------------------
  # 4. Load raster (force single layer)
  # ------------------------------------------------------------
  raster_file <- file.path(raster_root, paste0(site, "_quality_10yWindFell.asc"))
  r <- rast(raster_file)[[1]]
  
  # Binary suitability
  suit <- classify(r, rcl = matrix(c(-Inf, 0, 0,
                                     0, Inf, 1), ncol = 3, byrow = TRUE))
  
  suit_df <- as.data.frame(suit, xy = TRUE)
  value_col <- names(suit_df)[3]
  
  suit_df <- suit_df %>%
    rename(suitable = !!value_col) %>%
    mutate(
      x_coord = x,
      y_coord = y,
      suitability_fill = ifelse(suitable == 1, "Suitable", "Unsuitable")
    ) %>%
    select(x_coord, y_coord, suitable, suitability_fill)
  
  # ------------------------------------------------------------
  # 5. Load Pop file and normalise numeric-like columns
  # ------------------------------------------------------------
  all_df <- map_dfr(pop_files, ~ {
    df <- read_delim(.x, delim = "\t", show_col_types = FALSE)
    
    # Auto-detect numeric-like columns
    numeric_like <- names(df)[sapply(df, function(col) {
      all(grepl("^[-0-9.]*$", col) | is.na(col))
    })]
    
    # Convert them safely
    for (col in numeric_like) {
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    }
    
    df
  })
  
  df_occ <- all_df %>% filter(NInd > 0)
  
  n_reps_total <- length(unique(df_occ$Rep))
  
  # Convert grid indices → XY
  cells <- cellFromRowCol(r, df_occ$y, df_occ$x)
  xy <- xyFromCell(r, cells)
  
  df_occ$x_coord <- xy[,1]
  df_occ$y_coord <- xy[,2]
  
  # ------------------------------------------------------------
  # 6. Auto-detect initialisation site
  # ------------------------------------------------------------
  first_year <- min(df_occ$Year, na.rm = TRUE)
  init_row <- df_occ %>% filter(Year == first_year) %>% slice(1)
  init_x <- init_row$x_coord
  init_y <- init_row$y_coord
  
  # ------------------------------------------------------------
  # 7. Build occupancy PROPORTION rasters for each year
  # ------------------------------------------------------------
  occ_list <- list()
  
  for (yr in years) {
    df_y <- df_occ %>% filter(Year == yr)
    
    coords <- as.matrix(df_y[, c("x_coord", "y_coord")])
    cell_ids <- cellFromXY(r, coords)
    cell_ids <- as.numeric(cell_ids)
    
    occ_counts <- table(cell_ids)
    occ_prop <- as.numeric(occ_counts) / n_reps_total
    
    occ_r <- r[[1]]
    values(occ_r) <- NA
    occ_r[as.numeric(names(occ_counts))] <- occ_prop
    
    occ_df <- as.data.frame(occ_r, xy = TRUE)
    value_col <- names(occ_df)[3]
    
    occ_df <- occ_df %>%
      rename(occ = !!value_col) %>%
      mutate(
        x_coord = x,
        y_coord = y,
        Year = yr,
        occ_type = case_when(
          is.na(occ) ~ "NA",
          occ == 0   ~ "Zero",
          occ > 0    ~ "Positive"
        )
      ) %>%
      select(x_coord, y_coord, occ, occ_type, Year)
    
    occ_list[[as.character(yr)]] <- occ_df
  }
  
  occ_all <- bind_rows(occ_list)
  
  # ------------------------------------------------------------
  # 8. Plot facetted time series
  # ------------------------------------------------------------
  p <- ggplot() +
    
    geom_tile(
      data = suit_df,
      aes(x = x_coord, y = y_coord, fill = suitability_fill)
    ) +
    scale_fill_manual(
      name = "Suitability",
      values = c("Unsuitable" = "white", "Suitable" = "grey80")
    ) +
    
    new_scale_fill() +
    
    geom_tile(
      data = occ_all %>% filter(occ_type == "Zero"),
      aes(x = x_coord, y = y_coord, fill = occ_type),
      alpha = 0.8
    ) +
    scale_fill_manual(
      name = "Occupancy = 0",
      values = c("Zero" = "grey60")
    ) +
    
    new_scale_fill() +
    
    geom_tile(
      data = occ_all %>% filter(occ_type == "Positive"),
      aes(x = x_coord, y = y_coord, fill = occ),
      alpha = 0.8
    ) +
    scale_fill_viridis(
      name = "Occupancy (proportion)",
      option = "plasma",
      limits = c(0, 1)
    ) +
    
    geom_point(aes(x = init_x, y = init_y), colour = "red", size = 3, shape = 4) +
    
    coord_equal() +
    facet_wrap(~Year) +
    theme_minimal() +
    labs(
      title = paste("Site", site, "|", batch_id, "| Occupancy Proportion Time Series"),
      x = "X coordinate",
      y = "Y coordinate"
    )
  
  return(p)
}



plot_time_series(site = "2020", batch_id = "Batch2", years = 8)




# COMMAND ----------

install.packages("ggnewscale")
install.packages("viridis")

  library(tidyverse)
  library(terra)
  library(ggplot2)
  library(viridis)
  library(ggnewscale)

plot_time_series_sum <- function(
  site = "2020", 
  batch_id = "Batch20",
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/"
) {
  
  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  
  # 1. Latest site folder
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  site_dir <- sort(site_dirs, decreasing = TRUE)[1]
  
  # 2. Batch path
  batch_path <- file.path(site_dir, batch_id, "Outputs")
  
  # 3. Pop file
  pop_file <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)
  
  # Normalise numeric-like columns
  numeric_like <- names(df)[sapply(df, function(col) {
    all(grepl("^[-0-9.]*$", col) | is.na(col))
  })]
  for (col in numeric_like) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  
  # Keep only occupied rows
  df_occ <- df %>% filter(NInd > 0)
  if (nrow(df_occ) == 0) {
    message("No occupancy detected.")
    return(NULL)
  }
  
  # 4. Load raster
  r <- rast(file.path(raster_root, paste0(site, "_quality_10yWindFell.asc")))[[1]]
  
  # 5. Convert row/col indices → cell IDs (CRITICAL FIX)
  df_occ$x <- as.numeric(df_occ$x)
  df_occ$y <- as.numeric(df_occ$y)
  df_occ$cell_id <- cellFromRowCol(r, df_occ$y, df_occ$x)
  
  # 6. Compute occupancy proportion per cell
  total_reps <- length(unique(df$Rep))
  
  occ_prop <- df_occ %>%
    group_by(cell_id) %>%
    summarise(prop = n_distinct(Rep) / total_reps)
  
  # 7. Build raster
  prop_r <- r[[1]]
  values(prop_r) <- NA
  prop_r[occ_prop$cell_id] <- occ_prop$prop
  
  prop_df <- as.data.frame(prop_r, xy = TRUE)
  names(prop_df)[3] <- "prop"
  
  # 8. Suitability background
  suit <- classify(r, rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol=3, byrow=TRUE))
  suit_df <- as.data.frame(suit, xy = TRUE)
  names(suit_df)[3] <- "suitable"
  suit_df <- suit_df %>% mutate(suitability_fill = ifelse(suitable==1,"Suitable","Unsuitable"))
  
  # 9. Detect initialisation site (first occupied year)
  df_any_occ <- df %>% filter(NInd > 0)
  first_year <- min(df_any_occ$Year, na.rm = TRUE)
  init_row <- df_any_occ %>% filter(Year == first_year) %>% slice(1)
  
  init_cell <- cellFromRowCol(r, init_row$y, init_row$x)
  init_xy <- xyFromCell(r, init_cell)
  init_x <- init_xy[1]
  init_y <- init_xy[2]
  
  # 10. Plot
  ggplot() +
    geom_tile(data=suit_df, aes(x=x, y=y, fill=suitability_fill)) +
    scale_fill_manual(values=c("Unsuitable"="white","Suitable"="grey80")) +
    
    new_scale_fill() +
    
    geom_tile(
      data = prop_df %>% filter(!is.na(prop)),
      aes(x=x, y=y, fill=prop),
      alpha=0.9
    ) +
    scale_fill_viridis(
      name="Occupancy proportion",
      option="inferno",
      limits=c(0,1)
    ) +
    
    geom_point(aes(x=init_x, y=init_y), colour="red", size=3, shape=4) +
    
    coord_equal() +
    theme_minimal() +
    labs(
      title = paste("Site", site, "|", batch_id, "| Occupancy Proportion Across Replicates"),
      x = "X coordinate",
      y = "Y coordinate"
    )
}

plot_time_series_sum(site = "2070", batch_id = "Batch10")



# COMMAND ----------

plot_time_series_cummean <- function(
  site = "2020", 
  batch_id = "Batch20",
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/Kielder_updated2/",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/Kielder_updated2/"
) {
  
  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  
  # 1. Latest site folder
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  if (length(site_dirs) == 0) stop("No site folders found.")
  site_dir <- sort(site_dirs, decreasing = TRUE)[1]
  
  # 2. Batch path
  batch_path <- file.path(site_dir, batch_id, "Outputs")
  if (!dir.exists(batch_path)) stop("Batch folder not found.")
  
  # 3. Pop file
  pop_file <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE,
    recursive = FALSE
  )
  if (length(pop_file) == 0) stop("No Pop file found.")
  
  # 4. Load Pop + normalise numeric-like
  df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)
  
  numeric_like <- names(df)[sapply(df, function(col) {
    all(grepl("^[-0-9.]*$", col) | is.na(col))
  })]
  for (col in numeric_like) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  
  # 5. Occupancy: use ALL rows with NInd > 0
  df_occ <- df %>% filter(NInd > 0)
  if (nrow(df_occ) == 0) {
    message("No occupancy detected in any year.")
  }
  
  # 6. Raster
  r <- rast(file.path(raster_root, paste0(site, "_quality_10yWindFell.asc")))[[1]]
  
  # 7. Init site (first occupied year)
  df_any_occ <- df %>% filter(NInd > 0)
  if (nrow(df_any_occ) > 0) {
    first_year <- min(df_any_occ$Year, na.rm = TRUE)
    init_row <- df_any_occ %>% filter(Year == first_year) %>% slice(1)
    init_cell <- cellFromRowCol(r, init_row$y, init_row$x)
    init_xy <- xyFromCell(r, init_cell)
    init_x <- init_xy[1]
    init_y <- init_xy[2]
  } else {
    init_x <- NA
    init_y <- NA
  }
  
  # If still no occupancy at all, just plot habitat + init
  if (nrow(df_occ) == 0) {
    suit <- classify(r, rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol=3, byrow=TRUE))
    suit_df <- as.data.frame(suit, xy = TRUE)
    names(suit_df)[3] <- "suitable"
    suit_df <- suit_df %>% mutate(suitability_fill = ifelse(suitable==1,"Suitable","Unsuitable"))
    
    return(
      ggplot() +
        geom_tile(data = suit_df, aes(x=x, y=y, fill=suitability_fill)) +
        scale_fill_manual(values=c("Unsuitable"="white","Suitable"="grey80")) +
        geom_point(aes(x=init_x, y=init_y), colour="red", size=3, shape=4) +
        coord_equal() +
        theme_minimal() +
        labs(
          title = paste("Site", site, "|", batch_id, "| No occupancy detected"),
          x = "X coordinate",
          y = "Y coordinate"
        )
    )
  }
  
  # 8. XY → cell IDs
df_occ$x <- as.numeric(df_occ$x)
df_occ$y <- as.numeric(df_occ$y)

df_occ$cell_id <- cellFromRowCol(r, df_occ$y, df_occ$x)

  
  # 9. Cumulative mean pop size per cell
  mean_pop <- df_occ %>%
    group_by(cell_id) %>%
    summarise(mean_N = mean(NInd, na.rm = TRUE))
  
  mean_r <- r[[1]]
  values(mean_r) <- NA
  mean_r[mean_pop$cell_id] <- mean_pop$mean_N
  
  mean_df <- as.data.frame(mean_r, xy = TRUE)
  names(mean_df)[3] <- "mean_N"
  
  if (all(is.na(mean_df$mean_N))) {
    message("Mean population is NA everywhere — check NInd and Year in the Pop file.")
  }
  
  # 10. Suitability background
  suit <- classify(r, rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol=3, byrow=TRUE))
  suit_df <- as.data.frame(suit, xy = TRUE)
  names(suit_df)[3] <- "suitable"
  suit_df <- suit_df %>% mutate(suitability_fill = ifelse(suitable==1,"Suitable","Unsuitable"))
  
  # 11. Plot
  ggplot() +
    geom_tile(data=suit_df, aes(x=x, y=y, fill=suitability_fill)) +
    scale_fill_manual(values=c("Unsuitable"="white","Suitable"="grey80")) +
    new_scale_fill() +
    geom_tile(data = mean_df %>% filter(!is.na(mean_N)),
              aes(x=x, y=y, fill=mean_N), alpha=0.9) +
    scale_fill_viridis(name="Mean population size", option="plasma") +
    geom_point(aes(x=init_x, y=init_y), colour="red", size=3, shape=4) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = paste("Site", site, "|", batch_id, "| Cumulative Mean Population Size"),
      x = "X coordinate",
      y = "Y coordinate"
    )
}


plot_time_series_cummean(site = "2070", batch_id = "Batch10")




# COMMAND ----------

install.packages("ggnewscale")
install.packages("viridis")

  library(tidyverse)
  library(terra)
  library(ggplot2)
  library(viridis)
  library(ggnewscale)
install.packages("viridis")
install.packages("patchwork")
install.packages("ggnewscale")
#library(scico)
plot_time_series_sum <- function(
  site = "2020", 
  batch_id = "Batch20",
  buffer_km = "full",   # "full" or numeric (km)
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/"
) {
  
  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
   library(patchwork)
  # ------------------------------------------------------------
  # 1. Identify latest site folder
  # ------------------------------------------------------------
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  site_dir <- sort(site_dirs, decreasing = TRUE)[1]
  
  # ------------------------------------------------------------
  # 2. Batch path
  # ------------------------------------------------------------
  batch_path <- file.path(site_dir, batch_id, "Outputs")
  
  # ------------------------------------------------------------
  # 3. Pop file
  # ------------------------------------------------------------
  pop_file <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)
  
  # Normalise numeric-like columns
  numeric_like <- names(df)[sapply(df, function(col) {
    all(grepl("^[-0-9.]*$", col) | is.na(col))
  })]
  for (col in numeric_like) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  
  # Keep only occupied rows
  df_occ <- df %>% filter(NInd > 0)
  if (nrow(df_occ) == 0) {
    message("No occupancy detected.")
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # 4. Load raster
  # ------------------------------------------------------------
  r <- rast(file.path(raster_root, paste0(site, "_habcodes_fell.asc")))[[1]]
  #2024_habcodes_fell
  # ------------------------------------------------------------
  # 5. Convert row/col indices → cell IDs (CRITICAL FIX)
  # ------------------------------------------------------------
  df_occ$x <- as.numeric(df_occ$x)
  df_occ$y <- as.numeric(df_occ$y)
  df_occ$cell_id <- cellFromRowCol(r, df_occ$y, df_occ$x)
  
  # ------------------------------------------------------------
  # 6. Compute occupancy proportion per cell
  # ------------------------------------------------------------
  total_reps <- length(unique(df$Rep))
  
  occ_prop <- df_occ %>%
    group_by(cell_id) %>%
    summarise(prop = n_distinct(Rep) / total_reps)
  
  prop_r <- r[[1]]
  values(prop_r) <- NA
  prop_r[occ_prop$cell_id] <- occ_prop$prop
  
  # ------------------------------------------------------------
  # 7. Detect initialisation site
  # ------------------------------------------------------------
  df_any_occ <- df %>% filter(NInd > 0)
  first_year <- min(df_any_occ$Year, na.rm = TRUE)
  init_row <- df_any_occ %>% filter(Year == first_year) %>% slice(1)
  
  init_cell <- cellFromRowCol(r, init_row$y, init_row$x)
  init_xy <- xyFromCell(r, init_cell)
  init_x <- init_xy[1]
  init_y <- init_xy[2]
  
  # ------------------------------------------------------------
  # 8. Optional buffer clipping
  # ------------------------------------------------------------
  if (!is.character(buffer_km)) {
    init_vect <- vect(cbind(init_x, init_y), crs = crs(r))
    buf <- buffer(init_vect, width = buffer_km * 1000)
    
    # Mask occupancy raster
    prop_r <- mask(prop_r, buf)
    
    # Mask suitability raster
    suit_r <- classify(r, rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol=3, byrow=TRUE))
    suit_r <- mask(suit_r, buf)
  } else {
    suit_r <- classify(r, rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol=3, byrow=TRUE))
  }
  
  # Convert rasters to dataframes
  prop_df <- as.data.frame(prop_r, xy = TRUE)
  names(prop_df)[3] <- "prop"
  
  suit_df <- as.data.frame(suit_r, xy = TRUE)
  names(suit_df)[3] <- "suitable"
  suit_df <- suit_df %>% mutate(suitability_fill = ifelse(suitable==1,"Unoccupied","Unsuitable"))
  
  # ------------------------------------------------------------
  # 9. Plot
  # ------------------------------------------------------------
  ggplot() +
    geom_tile(data=suit_df, aes(x=x, y=y, fill=suitability_fill)) +
    scale_fill_manual(
  values = c("Unsuitable" = "grey99", "Unoccupied" = "grey75")
)+

    new_scale_fill() +
    
    geom_tile(
      data = prop_df %>% filter(!is.na(prop)),
      aes(x=x, y=y, fill=prop),
      alpha=0.9
    ) +
scale_fill_viridis(
  name = "Probability of establishment",
  option = "viridis",
  begin = 0,   # skip the yellow
  end   = 0.9,   # keep green → blue
  direction = -1,  # reverse so light → dark
  limits = c(0.01, 1)
)+

    
   # geom_point(aes(x=init_x, y=init_y), colour="red", size=3, shape=4) +
    
    coord_equal() +
    theme_minimal() +
    labs(title =  site)
}

p1 <- plot_time_series_sum("2020", "Batch3", buffer_km = 2)
p2 <- plot_time_series_sum("2050", "Batch3", buffer_km = 2)
p3 <- plot_time_series_sum("2070", "Batch3", buffer_km = 2)

(p1 + p3) + plot_layout(guides = "collect")



# COMMAND ----------

install.packages("ggnewscale")
  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  library(patchwork)
  library(ggspatial)

plot_time_series_sum <- function(
  site = "2020", 
  batch_id = "Batch20",
  buffer_km = "full",
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/Kielder_updated2/",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/Kielder_updated2/"
) {

  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  library(patchwork)
  library(ggspatial)

  # ------------------------------------------------------------
  # 1. Identify latest site folder
  # ------------------------------------------------------------
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  site_dir  <- sort(site_dirs, decreasing = TRUE)[1]

  batch_path <- file.path(site_dir, batch_id, "Outputs")

  pop_file <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE
  )

  df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)

  numeric_like <- names(df)[sapply(df, function(col)
    all(grepl("^[-0-9.]*$", col) | is.na(col))
  )]
  for (col in numeric_like) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))

  df_occ <- df %>% filter(NInd > 0)
  if (nrow(df_occ) == 0) return(NULL)

  # ------------------------------------------------------------
  # 2. Raster + cell mapping
  # ------------------------------------------------------------
  r <- rast(file.path(
    raster_root,
    paste0(site, "_habcodes_fell.asc")
  ))[[1]]

  df_occ$cell_id <- cellFromRowCol(r, df_occ$y, df_occ$x)

  total_reps <- length(unique(df$Rep))

  occ_prop <- df_occ %>%
    group_by(cell_id) %>%
    summarise(prop = n_distinct(Rep) / total_reps, .groups = "drop")

  prop_r <- r
  values(prop_r) <- NA
  prop_r[occ_prop$cell_id] <- occ_prop$prop

  # 🔧 CRITICAL FIX
  names(prop_r) <- "prop"

  # ------------------------------------------------------------
  # 3. Identify introduction point
  # ------------------------------------------------------------
  first_year <- min(df$Year[df$NInd > 0], na.rm = TRUE)
  init_row   <- df %>% filter(Year == first_year, NInd > 0) %>% slice(1)
  init_cell  <- cellFromRowCol(r, init_row$y, init_row$x)
  init_xy    <- xyFromCell(r, init_cell)

  # ------------------------------------------------------------
  # 4. Suitability + optional buffer
  # ------------------------------------------------------------
  suit_r <- classify(
    r,
    rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol = 3, byrow = TRUE)
  )

  if (!is.character(buffer_km)) {
    init_vect <- vect(init_xy, crs = crs(r))
    buf <- buffer(init_vect, width = buffer_km * 1000)
    prop_r <- mask(prop_r, buf)
    suit_r <- mask(suit_r, buf)
  }

  # ------------------------------------------------------------
  # 5. Data frames
  # ------------------------------------------------------------
  prop_df <- as.data.frame(prop_r, xy = TRUE, na.rm = TRUE) |>
    mutate(
      prob_class = cut(
        prop,
        breaks = c(0.01, 0.05, 0.2, 0.5, 1),
        labels = c("1–5%", "5–20%", "20–50%", ">50%"),
        include.lowest = TRUE
      )
    )
names(suit_r) <- "suit"
suit_df <- as.data.frame(suit_r, xy = TRUE) |>
  mutate(
    habitat = factor(
      ifelse(suit == 1, "Spruce", "Not spruce"),
      levels = c("Spruce", "Not spruce")
    )
  )

  # ------------------------------------------------------------
  # 6. Plot
  # ------------------------------------------------------------
  ggplot() +

    # Background habitat (very faint, but in legend)
    geom_tile(
      data = suit_df,
      aes(x = x, y = y, fill = habitat)
    ) +
    scale_fill_manual(
      name = "Habitat",
      values = c(
        "Spruce" = "grey90",
        "Not spruce"             = "grey97"
      ),
      guide = guide_legend(order = 2)
    ) +

    new_scale_fill() +

    # Probability of establishment
    geom_tile(
      data = prop_df,
      aes(x = x, y = y, fill = prob_class),
      alpha = 0.95
    ) +
    scale_fill_viridis_d(
      name = "Probability of establishment",
      option = "viridis",
      direction = 1,
      guide = guide_legend(order = 1)
    ) +

    # Introduction point
    geom_point(
      aes(x = init_xy[1], y = init_xy[2]),
      shape = 3,
      size = 1,
      fill = "white",
      colour = "red",
      stroke = 0.8
    ) +

    # Scale bar only (no grids)
    annotation_scale(
      location = "br",
      width_hint = 0.25,
      text_cex = 0.7
    ) +

    coord_equal() +

    labs(
      title = site
    ) +

    theme_minimal() +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.title   = element_blank(),
      axis.ticks   = element_blank(),
      plot.title   = element_text(face = "bold"),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
}

p1 <- plot_time_series_sum("2020", "Batch5", buffer_km = 2)
#p2 <- plot_time_series_sum("2050", "Batch3", buffer_km = 4)
p3 <- plot_time_series_sum("2070", "Batch5", buffer_km = 2)

(p1 + p3) + plot_layout(guides = "collect")


# COMMAND ----------

install.packages("ggnewscale")
  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  library(patchwork)
  library(ggspatial)

plot_time_series_sum <- function(
  site = "2020", 
  batch_id = "Batch20",
  buffer_km = "full",
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/Kielder_updated2/",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/Kielder_updated2/"
) {

  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  library(patchwork)
  library(ggspatial)

  # ------------------------------------------------------------
  # 1. Identify latest site folder
  # ------------------------------------------------------------
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  site_dir  <- sort(site_dirs, decreasing = TRUE)[1]

  batch_path <- file.path(site_dir, batch_id, "Outputs")

  pop_file <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE
  )

  df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)

  numeric_like <- names(df)[sapply(df, function(col)
    all(grepl("^[-0-9.]*$", col) | is.na(col))
  )]
  for (col in numeric_like) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))

  df_occ <- df %>% filter(NInd > 0)
  if (nrow(df_occ) == 0) return(NULL)

  # ------------------------------------------------------------
  # 2. Raster + cell mapping
  # ------------------------------------------------------------
  r <- rast(file.path(
    raster_root,
    paste0(site, "_habcodes_fell.asc")
  ))[[1]]

  df_occ$cell_id <- cellFromRowCol(r, df_occ$y, df_occ$x)

  total_reps <- length(unique(df$Rep))

  occ_prop <- df_occ %>%
    group_by(cell_id) %>%
    summarise(prop = n_distinct(Rep) / total_reps, .groups = "drop")

  prop_r <- r
  values(prop_r) <- NA
  prop_r[occ_prop$cell_id] <- occ_prop$prop

  # 🔧 CRITICAL FIX
  names(prop_r) <- "prop"

  # ------------------------------------------------------------
  # 3. Identify introduction point
  # ------------------------------------------------------------
  first_year <- min(df$Year[df$NInd > 0], na.rm = TRUE)
  init_row   <- df %>% filter(Year == first_year, NInd > 0) %>% slice(1)
  init_cell  <- cellFromRowCol(r, init_row$y, init_row$x)
  init_xy    <- xyFromCell(r, init_cell)

  # ------------------------------------------------------------
  # 4. Suitability + optional buffer
  # ------------------------------------------------------------
  suit_r <- classify(
    r,
    rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol = 3, byrow = TRUE)
  )

  if (!is.character(buffer_km)) {
    init_vect <- vect(init_xy, crs = crs(r))
    buf <- buffer(init_vect, width = buffer_km * 1000)
    prop_r <- mask(prop_r, buf)
    suit_r <- mask(suit_r, buf)
  }

  # ------------------------------------------------------------
  # 5. Data frames
  # ------------------------------------------------------------
  prop_df <- as.data.frame(prop_r, xy = TRUE, na.rm = TRUE) |>
    mutate(
      prob_class = cut(
        prop,
        breaks = c(0.01, 0.05, 0.2, 0.5, 1),
        labels = c("1–5%", "5–20%", "20–50%", ">50%"),
        include.lowest = TRUE
      )
    )
names(suit_r) <- "suit"
suit_df <- as.data.frame(suit_r, xy = TRUE) |>
  mutate(
    habitat = factor(
      ifelse(suit == 1, "Spruce", "Not spruce"),
      levels = c("Spruce", "Not spruce")
    )
  )

  # ------------------------------------------------------------
  # 6. Plot
  # ------------------------------------------------------------
  ggplot() +

    # Background habitat (very faint, but in legend)
    geom_tile(
      data = suit_df,
      aes(x = x, y = y, fill = habitat)
    ) +
    scale_fill_manual(
      name = "Habitat",
      values = c(
        "Spruce" = "grey90",
        "Not spruce"             = "grey97"
      ),
      guide = guide_legend(order = 2)
    ) +

    new_scale_fill() +

    # Probability of establishment
    geom_tile(
      data = prop_df,
      aes(x = x, y = y, fill = prob_class),
      alpha = 0.95
    ) +
    scale_fill_viridis_d(
      name = "Probability of establishment",
      option = "viridis",
      direction = 1,
      guide = guide_legend(order = 1)
    ) +

    # Introduction point
    geom_point(
      aes(x = init_xy[1], y = init_xy[2]),
      shape = 3,
      size = 1,
      fill = "white",
      colour = "red",
      stroke = 0.8
    ) +

    # Scale bar only (no grids)
    annotation_scale(
      location = "br",
      width_hint = 0.25,
      text_cex = 0.7
    ) +

    coord_equal() +

    labs(
      title = site
    ) +

    theme_minimal() +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.title   = element_blank(),
      axis.ticks   = element_blank(),
      plot.title   = element_text(face = "bold"),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
}

p1 <- plot_time_series_sum("2020", "Batch5", buffer_km = 2)
#p2 <- plot_time_series_sum("2050", "Batch3", buffer_km = 4)
p3 <- plot_time_series_sum("2070", "Batch5", buffer_km = 2)

(p1 + p3) + plot_layout(guides = "collect")

# COMMAND ----------

# ---- SAVE TIME SERIES ----

plot_time_series_sum <- function(
  site = "2020", 
  batch_id = "Batch20",
  buffer_km = "full",
  root = "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2//",
  raster_root = "/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/"
) {

  library(tidyverse)
  library(terra)
  library(ggnewscale)
  library(viridis)
  library(ggspatial)

  # ------------------------------------------------------------
  # 1. Identify latest site folder
  # ------------------------------------------------------------
  site_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  site_dirs <- site_dirs[grepl(paste0("^", site, "_batch"), basename(site_dirs))]
  site_dir  <- sort(site_dirs, decreasing = TRUE)[1]

  batch_path <- file.path(site_dir, batch_id, "Outputs")

  pop_file <- list.files(
    batch_path,
    pattern = "^Batch0_Sim0_Land1_Pop\\.txt$",
    full.names = TRUE
  )

  df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)

  numeric_like <- names(df)[sapply(df, function(col)
    all(grepl("^[-0-9.]*$", col) | is.na(col))
  )]
  for (col in numeric_like) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))

  df_occ <- df %>% filter(NInd > 0)
  if (nrow(df_occ) == 0) return(NULL)

  # ------------------------------------------------------------
  # 2. Raster + cell mapping (ROW/COL GRID COORDINATES)
  # ------------------------------------------------------------
  r <- rast(file.path(
    raster_root,
    paste0(site, "_habcodes_fell.asc")
  ))[[1]]

  if (ncell(r) == 0) stop("Raster failed to load — zero cells.")
  if (is.na(crs(r)) || crs(r) == "") stop("Raster has no CRS.")

  # Pop.txt uses row/col grid coordinates
  df_occ$cell_id <- cellFromRowCol(r, df_occ$y, df_occ$x)

  total_reps <- length(unique(df$Rep))

  occ_prop <- df_occ %>%
    group_by(cell_id) %>%
    summarise(prop = n_distinct(Rep) / total_reps, .groups = "drop")

  prop_r <- r
  values(prop_r) <- NA
  prop_r[occ_prop$cell_id] <- occ_prop$prop
  names(prop_r) <- "prop"

  # ------------------------------------------------------------
  # 3. Identify introduction point (ROW/COL → XY)
  # ------------------------------------------------------------
  first_year <- min(df$Year[df$NInd > 0], na.rm = TRUE)
  init_row   <- df %>% filter(Year == first_year, NInd > 0) %>% slice(1)

  init_cell  <- cellFromRowCol(r, init_row$y, init_row$x)
  init_xy    <- xyFromCell(r, init_cell)

  # ------------------------------------------------------------
  # 4. Suitability (Spruce vs Not spruce)
  # ------------------------------------------------------------
  suit_r <- classify(
    r,
    rcl = matrix(c(-Inf,0,0, 0,Inf,1), ncol = 3, byrow = TRUE)
  )
  names(suit_r) <- "suit"

  suit_df <- as.data.frame(suit_r, xy = TRUE) |>
    mutate(
      habitat = factor(
        ifelse(suit == 1, "Spruce", "Not spruce"),
        levels = c("Spruce", "Not spruce")
      )
    )

  # ------------------------------------------------------------
  # 5. Distance‑based risk rings (SAFE + STANDARDISED)
  # ------------------------------------------------------------
  init_vect <- vect(matrix(init_xy, ncol = 2), crs = crs(r))

  # compute distance raster ON THE RASTER GRID
  dist_r <- try(distance(r, init_vect), silent = TRUE)

  if (inherits(dist_r, "try-error") || ncell(dist_r) == 0) {
    stop("Distance raster computation failed — intro point likely outside raster.")
  }

  dist_r <- dist_r / 1000
  names(dist_r) <- "dist_km"

  # distance bands (0–10 km in 0.2 km steps)
  breaks_km <- seq(0, 10, by = 0.05)

  dist_band_r <- classify(
    dist_r,
    rcl = cbind(
      head(breaks_km, -1),
      tail(breaks_km, -1),
      seq_len(length(breaks_km) - 1)
    )
  )
  names(dist_band_r) <- "dist_band"

  # per-band risk = max establishment probability in that band
  df_risk <- as.data.frame(c(prop_r, dist_band_r), xy = FALSE, na.rm = TRUE) %>%
    group_by(dist_band) %>%
    summarise(risk = max(prop, na.rm = TRUE), .groups = "drop")

  risk_r <- dist_band_r
  values(risk_r) <- df_risk$risk[match(values(dist_band_r), df_risk$dist_band)]
  names(risk_r) <- "risk_ring"

  # Standardised categories
  risk_levels <- c("<2%", "2–5%", "5–20%", "20–50%", ">50%")

risk_df <- as.data.frame(risk_r, xy = TRUE, na.rm = TRUE) |>
  mutate(
    # ensure missing risk values become 0
    risk_ring = ifelse(is.na(risk_ring), 0, risk_ring),

    # categorise risk
    risk_cat = cut(
      risk_ring,
      breaks = c(0, 0.02, 0.05, 0.2, 0.5, 1),
      labels = risk_levels,
      include.lowest = TRUE,
      right = TRUE
    ),

    # force all categories to appear
    risk_cat = factor(risk_cat, levels = risk_levels)
  )

  # ------------------------------------------------------------
  # 6. Plot
  # ------------------------------------------------------------
  ggplot() +

    # Habitat
    geom_tile(data = suit_df, aes(x = x, y = y, fill = habitat)) +
    scale_fill_manual(
      name = "Habitat",
      values = c("Spruce" = "grey80", "Not spruce" = "grey97"),
      drop = FALSE,
      guide = guide_legend(order = 2)
    ) +

    new_scale_fill() +

    # Risk rings
    geom_tile(data = risk_df, aes(x = x, y = y, fill = risk_cat), alpha = 0.55) +
    scale_fill_manual(
      name = "Risk of spread",
      values = c(
        "<2%"   = "#fde725",
        "2–5%"   = "#5ec962",
        "5–20%"  = "#21918c",
        "20–50%" = "#3b528b",
        ">50%"   = "#440154"
      ),
      drop = FALSE,   # <<— KEEP ALL LEVELS IN LEGEND
      guide = guide_legend(order = 1)
    ) +

    # Intro point
    geom_point(
      aes(x = init_xy[1], y = init_xy[2]),
      shape = 3,
      size = 1.5,
      colour = "red",
      stroke = 0.9
    ) +

    annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.7) +
    coord_equal() +
    labs(title = site) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
}

p1 <- plot_time_series_sum("2020", "Batch9", buffer_km = 4)
#p2 <- plot_time_series_sum("2050", "Batch3", buffer_km = 4)
p3 <- plot_time_series_sum("2070", "Batch9", buffer_km = 4)

(p1 + p3) + plot_layout(guides = "collect")

# COMMAND ----------

library(tidyverse)
library(terra)

# ---- 1. Load raster ----
r <- rast("/dbfs/mnt/lab/unrestricted/Fiona/Kielder_updated2/2020_habcodes_fell.asc")[[1]]

cat("\nRASTER EXTENT:\n")
print(ext(r))

cat("\nRASTER CRS:\n")
print(crs(r))

# ---- 2. Load a Pop.txt ----
pop_file <- list.files(
  "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/Kielder_updated2/2020_batch20260306/Batch1/Outputs",
  pattern = "Pop\\.txt$",
  full.names = TRUE
)[1]

df <- read_delim(pop_file, delim = "\t", show_col_types = FALSE)

cat("\nFIRST 10 ROWS OF POP.TXT:\n")
print(head(df, 10))

cat("\nRANGE OF X AND Y IN POP.TXT:\n")
print(range(df$x, na.rm = TRUE))
print(range(df$y, na.rm = TRUE))

# ---- 3. Test whether Pop.txt coordinates fall inside raster extent ----
ext_r <- ext(r)

inside_x <- df$x >= ext_r$xmin & df$x <= ext_r$xmax
inside_y <- df$y >= ext_r$ymin & df$y <= ext_r$ymax

cat("\nPERCENT OF X INSIDE RASTER EXTENT:\n")
print(mean(inside_x, na.rm = TRUE))

cat("\nPERCENT OF Y INSIDE RASTER EXTENT:\n")
print(mean(inside_y, na.rm = TRUE))

# ---- 4. Test both interpretations: row/col vs XY ----

# A) Treat x,y as row/col indices
test_cell_rowcol <- try(cellFromRowCol(r, df$y[1], df$x[1]), silent = TRUE)

cat("\nTEST cellFromRowCol(r, y, x):\n")
print(test_cell_rowcol)

# B) Treat x,y as projected XY coordinates
test_cell_xy <- try(cellFromXY(r, cbind(df$x[1], df$y[1])), silent = TRUE)

cat("\nTEST cellFromXY(r, cbind(x, y)):\n")
print(test_cell_xy)

# ---- 5. Compute intro point using both methods ----
first_year <- min(df$Year[df$NInd > 0], na.rm = TRUE)
init_row <- df %>% filter(Year == first_year, NInd > 0) %>% slice(1)

cat("\nINTRO ROW FROM POP.TXT:\n")
print(init_row)

# A) intro cell using row/col
intro_cell_rowcol <- try(cellFromRowCol(r, init_row$y, init_row$x), silent = TRUE)
intro_xy_rowcol   <- try(xyFromCell(r, intro_cell_rowcol), silent = TRUE)

cat("\nINTRO CELL (row/col interpretation):\n")
print(intro_cell_rowcol)
cat("INTRO XY (row/col interpretation):\n")
print(intro_xy_rowcol)

# B) intro cell using XY
intro_cell_xy <- try(cellFromXY(r, cbind(init_row$x, init_row$y)), silent = TRUE)
intro_xy_xy   <- try(xyFromCell(r, intro_cell_xy), silent = TRUE)

cat("\nINTRO CELL (XY interpretation):\n")
print(intro_cell_xy)
cat("INTRO XY (XY interpretation):\n")
print(intro_xy_xy)


# COMMAND ----------

library(dplyr)
library(readr)
library(purrr)

test_batches <- c(
  "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/2020_batch20260306",
  "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/2020_batch20260308",
  "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/FoD_updated2/2020_batch20260310"
)

cat("Running TEST summarisation on", length(test_batches), "batches\n\n")

test_results <- map(test_batches, function(bd) {

  cat("--------------------------------------------------\n")
  cat("Summarising:", bd, "\n")

  res <- summarise_rs_batch(bd)

  cat("  Replicates rows :", nrow(res$replicates), "\n")
  cat("  Batch rows      :", nrow(res$by_batch), "\n")

  ts_file <- file.path(bd, "pop_metrics_timeseries.csv")

  if (file.exists(ts_file)) {
    ts <- read_csv(ts_file, show_col_types = FALSE)

    cat("  Time series rows:", nrow(ts), "\n")
    cat("  Unique reps     :", length(unique(ts$Rep)), "\n")
    cat("  Year range      :", min(ts$Year), "–", max(ts$Year), "\n")

    cat("  Max spread (m)              :", round(max(ts$max_spread_m, na.rm = TRUE), 1), "\n")
    cat("  Max cumulative unique cells :", max(ts$cum_unique_cells, na.rm = TRUE), "\n")
    cat("  Max new cells colonised     :", max(ts$new_cells_colonised, na.rm = TRUE), "\n")

    cat("\n--- Head of pop_metrics_timeseries.csv ---\n")
    print(head(ts, 10))
    cat("\n")

  } else {
    cat("  ⚠ pop_metrics_timeseries.csv NOT FOUND\n")
  }

  res
})

cat("--------------------------------------------------\n")
cat("✅ TEST summarisation complete\n")
