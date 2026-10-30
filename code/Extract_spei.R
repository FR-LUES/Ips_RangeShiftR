library(terra)
library(dplyr)
library(sf)  # for KMZ support

# Paths
case_study_path <- "C:/Users/fiona.plenderleith/OneDrive - Forest Research/TW IpsShift - IpsShift2 - Habitat model/Case_study_sites"
forest_shp <- file.path(case_study_path, "Forest_of_Dean_bbox.shp")
nym_shp <- file.path(case_study_path, "NYM_bbox.shp")
kielder_kmz <- file.path(case_study_path, "Kielder_smallerarea.kmz")

# Read polygons
fod_poly <- vect(forest_shp)
nym_poly <- vect(nym_shp)

# kielder
kmz_path <- "C:/Users/fiona.plenderleith/OneDrive - Forest Research/TW IpsShift - IpsShift2 - Habitat model/Case_study_sites/Kielder_smallerarea.kmz"
unzip(kmz_path, exdir = tempdir())  # extract to temp folder

# Find the KML file
kml_file <- list.files(tempdir(), pattern = "\\.kml$", full.names = TRUE)

# Read with terra
library(terra)
kielder_poly <- vect(kml_file)
#plot(kielder_poly)

# SPEI paths
filesp_folder <- "Z:/Common/FelixTrotter/Data/UKCP18_data/Europe/rcp85/V01/SPEI/6_mths/interp/y"
filesp_folder_sub <- "/SPEI_6m_V01_"
years <- 2015:2055

# Function to process polygon
process_polygon <- function(poly, site_name) {
  out <- data.frame()
  for (year in years) {
    filesp <- paste0(filesp_folder, year, filesp_folder_sub, year, "_09_interp.tif")
    r <- rast(filesp)
    r <- project(r, poly)
    clipped <- crop(r, poly)
    vals <- terra::extract(clipped, poly)[,2]
    
    out <- rbind(out, data.frame(
      site = site_name,
      year = year,
      med = median(vals, na.rm = TRUE),
      mean = mean(vals, na.rm = TRUE),
      min = min(vals, na.rm = TRUE),
      max = max(vals, na.rm = TRUE)
    ))
  }
  return(out)
}

# Process each and save separately
fod_result <- process_polygon(fod_poly, "Forest_of_Dean")
nym_result <- process_polygon(nym_poly, "NYM")
kielder_result <- process_polygon(kielder_poly, "Kielder_smallerarea")




write.csv(fod_result, "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Forest_of_Dean.csv", row.names = FALSE)
write.csv(nym_result, "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_NYM.csv", row.names = FALSE)
write.csv(kielder_result, "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Kielder.csv", row.names = FALSE)


### spei - addition of quality predcitiomn



# ---------------------------------------------------------
# FULL SCRIPT TO ADD 'pred' COLUMN TO ALL THREE CSV FILES
# ---------------------------------------------------------


SPEI_train <- seq(-4, 0, 0.5)
qual_train <- c(100, 90, 84, 73, 50, 20, 10, 4, 0)


# ---------------------------------------------------------
# Prediction function
# ---------------------------------------------------------
predict_quality <- function(values, SPEI_train, qual_train) {
  
  pred <- approx(
    x = SPEI_train,
    y = qual_train,
    xout = values,
    rule = 2
  )$y
  
  pred <- round(pred)
  pred <- pmax(pmin(pred, 100), 0)
  pred[values > 0] <- 0
  pred[is.na(pred)] <- 0
  
  return(pred)
}

# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------
csv_files <- c(
  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Kielder.csv",
  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_Forest_of_Dean.csv",
  "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/SPEI_NYM.csv"
)

# ---------------------------------------------------------
# Process each CSV
# ---------------------------------------------------------
for (f in csv_files) {
  
  df <- read.csv(f)
  
  if (!"SPEI" %in% names(df)) {
    stop(paste("Column 'mean' not found in:", f))
  }
  
  df$pred <- predict_quality(df$SPEI, SPEI_train, qual_train)
  
  write.csv(df, f, row.names = FALSE)
  message("Updated CSV saved: ", f)
}

# ---------------------------------------------------------
# Done
# ---------------------------------------------------------
message("All CSV files processed successfully.")


