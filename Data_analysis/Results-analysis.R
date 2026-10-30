library(ggplot2)
library(sf)

# These are just some ideas for potential visualisations of the simulation outputs
# These may need to be modified and others added

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ips-shift-functions.R")


# Use the timestamp used as the name for the batch of simulations
sim_batch_timestamp <- "2023-09-21-15h32m50s"

output_path <- "../outputs"

sim_batch_path <- file.path(output_path, sim_batch_timestamp)

df_sims <- read.csv(file.path(sim_batch_path, "sim-parameters.csv"))

# Arrays to store mean and SD of population per sim
mean_pop_array <- c()
sd_pop_array <- c()

# Arrays to store mean and SD of distance travelled east per sim
mean_east_dist_array <- c()
sd_east_dist_array <- c()

# List of final populations for each simulation
# Currently used to generate violin plot of population distribution
list_final_pop <- vector(mode = "list", length = 100)

for(sim_n in 1:nrow(df_sims)) {
  row <- df_sims[sim_n,]
  
  # Read the population file and calculate mean and SD of population in final year
  df_pop <- read.csv(file.path(sim_batch_path, row$sim_uuid, "Outputs", "Batch0_Sim0_Land1_Range.txt"), sep="\t", header=T)
  df_summary <- aggregate(data = subset(df_pop, Year == max(df_pop$Year)), NInds ~ Rep, function(x) sum(x))
  mean_pop_array <- c(mean_pop_array, mean(df_summary$NInd))
  sd_pop_array <- c(sd_pop_array, sd(df_summary$NInd))
  
  # Read the range file and calculate mean and SD of total distance travelled east in final year
  df_range <- read.csv(file.path(sim_batch_path, row$sim_uuid, "Outputs", "Batch0_Sim0_Land1_Range.txt"), sep="\t", header=T)
  START_X <- subset(df_range, Year == 0)$min_X[1]
  LAST_YEAR <- max(df_range$Year)
  REPLICATES <- max(df_range$Rep)
  df_summary <- aggregate(data = subset(df_range, Year == LAST_YEAR), min_X ~ Rep, function(x) sum(x))
  df_summary$east_dist <- START_X - df_summary$min_X
  mean_east_dist_array <- c(mean_east_dist_array, mean(df_summary$east_dist))
  sd_east_dist_array <- c(sd_east_dist_array, sd(df_summary$east_dist))
  
  
  # Save the final population for each replicate in each sim, to provide
  # population data for violin plots
  for (rep in 1:REPLICATES) {
    path_geojson <- file.path(sim_batch_path, row$sim_uuid, "Output_Maps", sprintf("rep1_year%s_pop.geojson", LAST_YEAR))
    if (file.exists(path_geojson)) {
      pop_geom <- read_sf(dsn = path_geojson)
      pop_geom$sim_uuid <- row$sim_uuid
      pop_geom$label <- paste("rmax:", row$rmax, "YFD:", row$years_for_destruction, "Annual:", row$apply_destruction_annually)
      pop_geom$rmax <- row$rmax
      pop_geom$years_for_destruction <- row$years_for_destruction
      pop_geom$apply_destruction_annually <- row$apply_destruction_annually
      list_final_pop[[rep * sim_n]] <- pop_geom
    }
    
  }
}

# Save all the summary population and distance data to the simulation dataframe
df_sims$mean_pop <- mean_pop_array
df_sims$sd_pop <- sd_pop_array
df_sims$mean_east_dist <- mean_east_dist_array
df_sims$sd_east_dist <- sd_east_dist_array
write.csv(df_sims, file.path(sim_batch_path, "sim-parameters-with-summary-results.csv"), row.names = FALSE)

# Create a label for the output plot
# THIS MAY NEED CHANGING DEPENDING ON WHAT PARAMETERS ARE BEING TESTED
df_sims$label <- paste("rmax:", df_sims$rmax, "YFD:", df_sims$years_for_destruction, "Annual:", df_sims$apply_destruction_annually)
df_sims <- df_sims[order(df_sims$rmax, df_sims$years_for_destruction, df_sims$apply_destruction_annually),]
df_sims$label <- factor(df_sims$label, levels = df_sims$label)


ggplot(df_sims) +
  geom_bar(aes(x=label, y=mean_pop, fill=rmax), stat="identity") +
  geom_errorbar(aes(x=label, ymin=mean_pop-sd_pop, ymax=mean_pop+sd_pop)) +
  scale_y_continuous(name=sprintf("Population after %s years", unique(df_sims$years)), labels = scales::comma) +
  coord_flip() +
  theme_bw()

ggsave(file.path(sim_batch_path, "population-barplot.jpg"), height=10, width=8, dpi=300)


ggplot(df_sims) +
  geom_bar(aes(x=label, y=mean_east_dist/1000, fill=rmax), stat="identity") +
  geom_errorbar(aes(x=label, ymin=(mean_east_dist-sd_east_dist)/1000, ymax=(mean_east_dist+sd_east_dist)/1000)) +
  scale_y_continuous(name=sprintf("Distance travelled east after %s years (km)", unique(df_sims$years)), labels = scales::comma) +
  coord_flip() +
  theme_bw()

ggsave(file.path(sim_batch_path, "east-distance-barplot.jpg"), height=10, width=8, dpi=300)


df_final_pop <- do.call("rbind", list_final_pop)
df_final_pop = df_final_pop[, c("label", "rmax", "years_for_destruction", "NInd","geometry")]
df_final_pop <- df_final_pop[order(df_final_pop$rmax, df_final_pop$years_for_destruction),]
df_final_pop$label <- factor(df_final_pop$label, levels = unique(df_final_pop$label))

df_final_pop$longitude <- st_coordinates(df_final_pop)[,1]
df_final_pop$latitude <- st_coordinates(df_final_pop)[,2]

df_final_pop_sum <- aggregate(data = df_final_pop, NInd ~ label + rmax + years_for_destruction + longitude + latitude, function(x) sum(x))

df <- data.frame(df_final_pop_sum[rep(seq_len(dim(df_final_pop_sum)[1]), df_final_pop_sum$NInd), c("label", "rmax", "latitude", "longitude"), drop = FALSE], row.names=NULL)

ggplot(df, aes(x=label, y=longitude, fill=rmax)) +
  geom_violin(bw = 500) +
  coord_flip() +
  theme_bw()

ggsave(file.path(sim_batch_path, "population-longitude-violin.jpg"), height=10, width=8, dpi=300)

# Read the population file and calculate mean and SD of population in final year
