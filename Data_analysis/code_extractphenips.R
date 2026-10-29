
library(readxl)
library(terra)

setwd("C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/PHENIPS")
LatitudeGrid=as.matrix(read_excel('UKCP_PHENIPS_Generations_ProjNum_05.xlsx',sheet = "LatGrid", col_names = FALSE))
LongitudeGrid=as.matrix(read_excel('UKCP_PHENIPS_Generations_ProjNum_05.xlsx',sheet = "LonGrid", col_names = FALSE)) 

Xmin = min(LongitudeGrid)
Xmax = max(LongitudeGrid)
Ymin = min(LatitudeGrid)
Ymax = max(LatitudeGrid)

# 
# image(baseline)
# 
# Generations2021 <- read_excel('UKCP_PHENIPS_test.xlsx',sheet = "2035", col_names = FALSE,  na = "NA")
# Generations2021$...111 <- NA
# Generations2021$...112 <- NA
# for(i in 1:3){Generations2021 <- rbind(NA,Generations2021)}
# 
# Generations2021 <- as.matrix(Generations2021)
# Generations2021 <- t(Generations2021)[ncol(Generations2021):1,]
# image(Generations2021)
# 
# baseline <- rast(ymin = Ymin, xmin = Xmin, ymax = Ymax, xmax = Xmax, nrow = 112, ncol = 82, vals = Generations2020)
# 
#  # Forest of Dean - latlong
# 
# FoD <- crop(baseline, clipped_extent)
# 
# plot(baseline)
# plot(FoD, add=T, col=terrain.colors(3))
# 
# df <- freq(FoD)
# 



# now make a loop
years <- as.character(c(2020:2035))
FoD_full <- data.frame(layer = numeric(0),
            value = numeric(0),
            count = numeric(0),
            Year = numeric(0))

for(year in years){
  y <- as.character(year)
  shape <- read_excel('UKCP_PHENIPS_Generations_ProjNum_01.xlsx', sheet = y, col_names = FALSE)
  shape$...111 <- NA # NAs were removed when data laods 
  shape$...112 <- NA # NAs were removed when data laods 
  for(i in 1:3){shape <- rbind(NA,shape)} # NAs were removed when data laods 
  shape <- as.matrix(shape) # convert to matric
  shape <- t(shape)[ncol(shape):1,]
  shape_rast <- rast(ymin = Ymin, xmin = Xmin, ymax = Ymax, xmax = Xmax, nrow = 112, ncol = 82, vals = shape)
  FoD <- crop(shape_rast, terra::ext(-2.724827, -2.364937, 51.696508, 51.922954)) #clip raster for shape
  ### Calculate the Class statistics in each grid
  FoD_year<- as.data.frame(freq(FoD))
  FoD_year$Year <- y
  
  FoD_full <- rbind(FoD_full, FoD_year)
  
}



library(terra)
library(readxl)

# Years
years <- as.character(2020:2040)

# Site extents (BNG)
site_extents_bng <- list(
  Kielder = ext(352643.5, 369793.5, 578384.813, 595534.813),
  FoD     = ext(352888.188, 370038.188, 205389.813, 222539.813),
  NYM     = ext(476175.000, 493325.000, 484860.344, 502010.344)
)


ext_bng_to_wgs84 <- function(e) {
  # Corner points of the BNG extent
  pts <- vect(
    cbind(
      x = c(e$xmin, e$xmax, e$xmin, e$xmax),
      y = c(e$ymin, e$ymin, e$ymax, e$ymax)
    ),
    crs = "EPSG:27700"
  )
  
  # Reproject to WGS84
  pts_wgs <- project(pts, "EPSG:4326")
  
  # Build WGS84 extent
  ext(min(pts_wgs$x), max(pts_wgs$x),
      min(pts_wgs$y), max(pts_wgs$y))
}

site_extents_wgs <- lapply(site_extents_bng, ext_bng_to_wgs84)




# Output container
FoD_full <- data.frame(
  layer = numeric(0),
  value = numeric(0),
  count = numeric(0),
  Year  = character(0),
  Site  = character(0)
)

# ---- MAIN LOOP ----
for (year in years) {
  
  # Read the sheet
  shape <- read_excel("UKCP_PHENIPS_Generations_ProjNum_11.xlsx",
                      sheet = year, col_names = FALSE)
  
  # Add missing NA columns (removed on load)
  shape$...111 <- NA
  shape$...112 <- NA
  
  # Add 3 NA rows at top
  for (i in 1:3) shape <- rbind(NA, shape)
  
  # Convert to matrix and flip/rotate
  shape <- as.matrix(shape)
  shape <- t(shape)[ncol(shape):1,]
  
  # Create full raster
  shape_rast <- rast(ymin = Ymin, xmin = Xmin,
                     ymax = Ymax, xmax = Xmax,
                     nrow = 112, ncol = 82,
                     vals = shape)
  
  # --- LOOP THROUGH SITES ---
  for (site_name in names(site_extents_bng)) {
    
    site_ext <- site_extents_bng[[site_name]]
    
    # Crop to site
    site_r <- crop(shape_rast, site_ext)
    
    # Frequency table
    site_freq <- as.data.frame(freq(site_r))
    site_freq$Year <- year
    site_freq$Site <- site_name
    
    # Append
    FoD_full <- rbind(FoD_full, site_freq)
  }
}

# Resulting table:
FoD_full




# SPEI analysis for IPs

library(plyr)
library(terra)
library(dplyr)
library(tidyr)
FoD <- crop(shape_rast, terra::ext(-2.724827, -2.364937, 51.696508, 51.922954))


filesp_folder <- "Z:Common/FelixTrotter/Data/UKCP18_data/Europe/rcp85/V01/SPEI/6_mths/interp/y"
filesp_folder_sub <- "/SPEI_6m_V01_"

filetmx_folder <- "Z:Common/FelixTrotter/Data/UKCP18_data/Europe/rcp85/V01/tasmax/y"
filetmx_folder_sub <- "/tasmax_V01_"



result_final <- data.frame(
  med = numeric(0),
  mean = numeric(0),
  min = numeric(0),
  max = numeric(0),
  year = numeric(0)
)
years <- c(2015:2040)
months <- c(6:8)

for (year in years) {
  filesp <-
    paste0(filesp_folder,
           year,
           filesp_folder_sub,
           year,
           "_09_interp.tif")
  r <- rast(filesp)
  r <- project(r, shape_rast)
  clipped <- crop(r, terra::ext(-2.724827, -2.364937, 51.696508, 51.922954))
  grid.data <- as.data.frame(freq(clipped))
  
  full_data <- uncount(grid.data, weights = count)
  result <- data.frame(
    med = median(full_data$value),
    mean = mean(full_data$value),
    min = min(full_data$value),
    max = max(full_data$value)
  )
  
  result_year <- result
  result_year$year <- year
  result_final <- rbind(result_year, result_final)
}
write.csv(result_final, "C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/PHENIPS/SPEI.csv")








# ==========================================================
# PHENIPS summary by site-year (mode class -> Rmax), WGS84
# ==========================================================
library(terra)
library(readxl)

setwd("C:/Users/fiona.plenderleith/OneDrive - Forest Research/Documents/IpsShift/PHENIPS")
LatitudeGrid=as.matrix(read_excel('UKCP_PHENIPS_Generations_ProjNum_05.xlsx',sheet = "LatGrid", col_names = FALSE))
LongitudeGrid=as.matrix(read_excel('UKCP_PHENIPS_Generations_ProjNum_05.xlsx',sheet = "LonGrid", col_names = FALSE)) 

Xmin = min(LongitudeGrid)
Xmax = max(LongitudeGrid)
Ymin = min(LatitudeGrid)
Ymax = max(LatitudeGrid)
# ---------------------------------------------------------
# 0) Helper functions
# ---------------------------------------------------------
# Mode function (most frequent value in a vector; ignores NA)

get_mode_from_freq <- function(r) {
  fr <- terra::freq(r)
  fr[which.max(fr[, "count"]), "value"]
}


# PHENIPS class -> Rmax lookup
rmax_lookup <- c(
  `1` = 5,
  `2` = 7.5,
  `3` = 32.5,
  `4` = 45,
  `5` = 165
)

# ---------------------------------------------------------
# 1) Load site TIFFs (BNG)
#    (paths as per your last message; names corrected to match sites)
# ---------------------------------------------------------
r_fod_bng     <- rast("Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoD/2050_quality_both.tif")
r_nym_bng     <- rast("Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/NYM/2050_quality_both.tif")
r_kielder_bng <- rast("Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/Kielder/2050_quality_both.tif")

# ---------------------------------------------------------
# 2) Reproject site rasters to WGS84 and take their extents
# ---------------------------------------------------------
r_fod_wgs     <- project(r_fod_bng,     "EPSG:4326")
r_nym_wgs     <- project(r_nym_bng,     "EPSG:4326")
r_kielder_wgs <- project(r_kielder_bng, "EPSG:4326")

site_extents_wgs <- list(
  FoD     = ext(r_fod_wgs),
  NYM     = ext(r_nym_wgs),
  Kielder = ext(r_kielder_wgs)
)

# ---------------------------------------------------------
# 3) Years and output containers
# ---------------------------------------------------------
years <- as.character(2020:2040)

FoD_full     <- data.frame(Year = character(), PHENIPS = numeric(), Rmax = numeric(), stringsAsFactors = FALSE)
NYM_full     <- data.frame(Year = character(), PHENIPS = numeric(), Rmax = numeric(), stringsAsFactors = FALSE)
Kielder_full <- data.frame(Year = character(), PHENIPS = numeric(), Rmax = numeric(), stringsAsFactors = FALSE)

# ---------------------------------------------------------
# 4) Loop through PHENIPS sheets and summarize per site
# ---------------------------------------------------------
for (year in years) {
  message("Processing year: ", year)

  # ---- Load PHENIPS sheet ----
  shape <- read_excel("UKCP_PHENIPS_Generations_ProjNum_01.xlsx",
                      sheet = year, col_names = FALSE)

  # Restore NA rows/columns removed on import
  shape$...111 <- NA
  shape$...112 <- NA
  for (i in 1:3) shape <- rbind(NA, shape)

  # Matrix grid transformation
  shape <- as.matrix(shape)
  shape <- t(shape)[ncol(shape):1, ]

  # Build PHENIPS raster (assumes WGS84 extents below are degrees)
  r_year <- rast(
    ymin = Ymin, xmin = Xmin,
    ymax = Ymax, xmax = Xmax,
    nrow = 112, ncol = 82,
    vals = shape
  )

  # Ensure CRS is WGS84
  if (is.na(crs(r_year))) crs(r_year) <- "EPSG:4326"

  # ---- Extract for each site: compute mode PHENIPS and map to Rmax ----
  
  ## FoD
  r_fod_crop <- crop(r_year, site_extents_wgs$FoD)
  mode_fod   <- get_mode_from_freq(r_fod_crop)
  rmax_fod   <- unname(rmax_lookup[as.character(mode_fod)])
  FoD_full   <- rbind(
    FoD_full,
    data.frame(Year = year, PHENIPS = mode_fod, Rmax = rmax_fod)
  )
  
  ## NYM
  r_nym_crop <- crop(r_year, site_extents_wgs$NYM)
  mode_nym   <- get_mode_from_freq(r_nym_crop)
  rmax_nym   <- unname(rmax_lookup[as.character(mode_nym)])
  NYM_full   <- rbind(
    NYM_full,
    data.frame(Year = year, PHENIPS = mode_nym, Rmax = rmax_nym)
  )
  
  ## Kielder
  r_kielder_crop <- crop(r_year, site_extents_wgs$Kielder)
  mode_kie       <- get_mode_from_freq(r_kielder_crop)
  rmax_kie       <- unname(rmax_lookup[as.character(mode_kie)])
  Kielder_full   <- rbind(
    Kielder_full,
    data.frame(Year = year, PHENIPS = mode_kie, Rmax = rmax_kie)
  )
  
}

# ---------------------------------------------------------
# 5) (Optional) Combine to one long table if useful
# ---------------------------------------------------------
all_sites <- rbind(
  transform(FoD_full,     Site = "FoD"),
  transform(NYM_full,     Site = "NYM"),
  transform(Kielder_full, Site = "Kielder")
)

# ---------------------------------------------------------
# 6) (Optional) Write CSVs
# ---------------------------------------------------------
# write.csv(FoD_full,     "FoD_PHENIPS_Rmax_summary.csv", row.names = FALSE)
# write.csv(NYM_full,     "NYM_PHENIPS_Rmax_summary.csv", row.names = FALSE)
# write.csv(Kielder_full, "Kielder_PHENIPS_Rmax_summary.csv", row.names = FALSE)
# write.csv(all_sites,    "AllSites_PHENIPS_Rmax_summary.csv", row.names = FALSE)

# ==========================================================
# 7) FINAL MAP: plot last computed PHENIPS + site extents
# ==========================================================
plot(r_year, main = paste0("PHENIPS ", year, " with Case Study Extents (WGS84)"))

cols <- c("#D55E00", "#CC79A7", "#A52A2A")

i <- 1
for (nm in names(site_extents_wgs)) {
  e <- site_extents_wgs[[nm]]
  rect(e$xmin, e$ymin, e$xmax, e$ymax, border = cols[i], lwd = 3)
  text((e$xmin + e$xmax)/2, (e$ymin + e$ymax)/2, labels = nm, col = cols[i], cex = 1.2)
  i <- i + 1
}

# ----------------------------------------------------------
# Save PHENIPS summaries by Site × Decade
# ----------------------------------------------------------

# Ensure Year is numeric
all_sites$Year <- as.numeric(all_sites$Year)

# Compute decade block (e.g. 2020–2029 → 2020; 2030–2039 → 2030)
all_sites$Decade <- floor(all_sites$Year / 10) * 10

# Base output directory
base_dir <- "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data"

# Unique Sites and Decades
sites   <- unique(all_sites$Site)
decades <- unique(all_sites$Decade)

for (s in sites) {
  
  # Safe file name for spaces
  s_safe <- gsub(" ", "_", s)
  
  for (d in decades) {
    
    # Filter table
    subdf <- subset(all_sites, Site == s & Decade == d)
    
    # Skip empty subsets
    if (nrow(subdf) == 0) next
    
    # Output file name
    outfile <- file.path(base_dir,
                         paste0("phenips", s_safe, "", d, ".csv"))
    
    # Save
    write.csv(subdf, outfile, row.names = FALSE)
    
    message("Saved: ", outfile)
  }
}

# ----------------------------------------------------------
# End of save‑segment
# ----------------------------------------------------------
