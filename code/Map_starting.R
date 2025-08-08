library(terra)
library(ggplot2)
library(dplyr)



library(terra)

# Load the full raster
r <- rast("Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/processed/2035_quality.asc")
library(terra)
# Get resolution and extent
res_x <- res(r)[1]
res_y <- res(r)[2]
ext_r <- ext(r)

# Calculate center coordinates
center_x <- (ext_r[1] + ext_r[2]) / 2
center_y <- (ext_r[3] + ext_r[4]) / 2

# Calculate half-width and half-height of 200 cells
half_width <- (200 / 2) * res_x
half_height <- (200 / 2) * res_y

# Define cropping extent
crop_extent <- ext(center_x - half_width, center_x + half_width,
                   center_y - half_height, center_y + half_height)


# Crop raster
r_crop <- crop(r, crop_extent)


r_crop_clean <- ifel(is.na(r_crop), 0, r_crop)
# Save cropped raster (optional)
writeRaster(r_crop_clean, "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/processed/quality_cropped.asc", filetype="AAIGrid", NAflag=-9999, overwrite=TRUE)
writeRaster(r_crop_clean, "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/processed/quality_cropped.tif")

# Plot to check
plot(r_crop_clean, main = "Cropped 200x200 Center of 2035 Suitability Map")



# Load raster
r <- rast("Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/processed/quality_cropped.asc")

# Convert raster to data frame for ggplot
r_df <- as.data.frame(r, xy = TRUE)
colnames(r_df) <- c("x", "y", "value")

# Load or create starting point raster (e.g., from initial_dist.tif)
init <- rast("Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/processed/initial_dist.tif")
init_df <- as.data.frame(init, xy = TRUE)
colnames(init_df) <- c("x", "y", "pop")

# Filter to only show starting points
start_points <- init_df %>% filter(pop == 1)

# Plot
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(name = "Suitability", na.value = "transparent") +
  geom_point(data = start_points, aes(x = x, y = y), color = "red", size = 2) +
  coord_equal() +
  theme_minimal() +
  labs(title = "2035 Suitability Map with Starting Point(s)",
       x = "Easting", y = "Northing")

