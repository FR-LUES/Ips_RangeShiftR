# Databricks notebook source
# MAGIC %r
# MAGIC library(dplyr)
# MAGIC library(ggplot2)
# MAGIC library(scales)
# MAGIC library(gridExtra)
# MAGIC library(tidyr)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 1) Load yearly dataset WITH metadata
# MAGIC # ==================================================
# MAGIC recode_site <- function(x) {
# MAGIC   dplyr::recode_factor(
# MAGIC     x,
# MAGIC     "FoD_updated2"     = "FoD",
# MAGIC     "Kielder_updated2" = "Kielder",
# MAGIC     "NYM_updated2"     = "NYM"
# MAGIC   )
# MAGIC }
# MAGIC
# MAGIC pd <- position_dodge(width = 0.5)
# MAGIC
# MAGIC the_theme <- theme_classic(base_size = 12) +
# MAGIC   theme(
# MAGIC     legend.position = "right",
# MAGIC     axis.title.x = element_text(margin = margin(t = 12)),
# MAGIC     axis.title.y = element_text(margin = margin(r = 12)),
# MAGIC     axis.text     = element_text(size = 12),
# MAGIC     legend.text   = element_text(size = 12),
# MAGIC     legend.title  = element_text(size = 12)
# MAGIC   )
# MAGIC
# MAGIC ts_all <- read.csv(
# MAGIC   "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/ALL_pop_metrics_timeseries_with_meta.csv"
# MAGIC )
# MAGIC
# MAGIC ts_all %>%
# MAGIC   mutate(
# MAGIC     Rep_global = paste(batch, Rep, sep = "_")   # ✅ ADD THIS
# MAGIC   )
# MAGIC
# MAGIC ts_sub <- ts_all %>%
# MAGIC   filter(Year >= 0, Year <= 10) %>%
# MAGIC   mutate(
# MAGIC     site = recode_site(site),
# MAGIC     year = factor(year),
# MAGIC     Rep_global = paste(batch, Rep, sep = "_")   # ✅ ADD THIS
# MAGIC   )
# MAGIC
# MAGIC
# MAGIC baseline_totals <- ts_sub %>%
# MAGIC   filter(Year == min(Year)) %>%
# MAGIC   group_by(site, year, out_subdir) %>%
# MAGIC   summarise(reps_total = n_distinct(Rep_global), .groups = "drop")
# MAGIC
# MAGIC # ==================================================
# MAGIC # 2) Compute ONE summary per REAL batch (Batch1, Batch2, …)
# MAGIC # ==================================================
# MAGIC
# MAGIC batch_summary <- ts_sub %>%
# MAGIC   group_by(site, year, out_subdir, Year) %>%
# MAGIC   summarise(
# MAGIC     occ_median    = median(occupied_cells, na.rm = TRUE),
# MAGIC     spread_median = median(max_spread_m_monotonic, na.rm = TRUE),
# MAGIC     reps_alive = n_distinct(Rep_global),
# MAGIC     .groups = "drop"
# MAGIC   ) %>%
# MAGIC   left_join(baseline_totals,
# MAGIC             by = c("site", "year", "out_subdir")) %>%
# MAGIC   mutate(
# MAGIC     p_extinct = 1 - reps_alive / reps_total
# MAGIC   )
# MAGIC
# MAGIC
# MAGIC
# MAGIC
# MAGIC # ==================================================
# MAGIC # 3) Between‑batch summaries (SURVIVORS ONLY for occupancy/spread)
# MAGIC # ==================================================
# MAGIC
# MAGIC between_batch_survivors <- batch_summary %>%
# MAGIC   group_by(site, year, Year) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = median(occ_median),
# MAGIC     occ_lo  = quantile(occ_median, 0.25),
# MAGIC     occ_hi  = quantile(occ_median, 0.75),
# MAGIC
# MAGIC     spread_medi = median(spread_median),
# MAGIC     spread_lo  = quantile(spread_median, 0.25),
# MAGIC     spread_hi  = quantile(spread_median, 0.75),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 4) Between‑batch extinction (ALL batches)
# MAGIC # ==================================================
# MAGIC
# MAGIC between_batch_extinction <- batch_summary %>%
# MAGIC   group_by(site, year, Year) %>%  
# MAGIC   summarise(
# MAGIC     ext_median = median(p_extinct, na.rm = TRUE),
# MAGIC     ext_lo  = quantile(p_extinct, 0.25, na.rm = TRUE),
# MAGIC     ext_hi  = quantile(p_extinct, 0.75, na.rm = TRUE),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 5) Collapse to site × decade × year
# MAGIC # ==================================================
# MAGIC
# MAGIC site_decade_summary <- between_batch_survivors %>%
# MAGIC   left_join(
# MAGIC     between_batch_extinction,
# MAGIC     by = c("site", "year", "Year")
# MAGIC   )
# MAGIC # ==================================================
# MAGIC # 6) PLOTS
# MAGIC # ==================================================
# MAGIC # ---- Extinction plot ----
# MAGIC ext_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, ext_median, colour = year, group = year)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = year),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site, nrow = 1) +
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
# MAGIC   scale_y_continuous(labels = percent_format(accuracy = 1)) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Climate decade",
# MAGIC     fill   = "Climate decade"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme
# MAGIC
# MAGIC # ---- Occupancy plot ----
# MAGIC occ_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, occ_medi, colour = year, group = year)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = year),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site, nrow = 1) +
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
# MAGIC   #scale_y_log10(labels = label_comma()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Climate decade",
# MAGIC     fill   = "Climate decade"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme
# MAGIC
# MAGIC # ---- Spread plot ----
# MAGIC spread_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, spread_medi, colour = year, group = year)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = year),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site, nrow = 1) +
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Climate decade",
# MAGIC     fill   = "Climate decade"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme
# MAGIC
# MAGIC # ==================================================
# MAGIC # 7) Arrange vertically
# MAGIC # ==================================================
# MAGIC install.packages("patchwork")
# MAGIC library(patchwork)
# MAGIC
# MAGIC (ext_plot / occ_plot / spread_plot) + plot_layout(guides = "collect")
# MAGIC
# MAGIC
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC # ==================================================
# MAGIC # RESULTS SECTION SUMMARY STATS (clean + no warnings)
# MAGIC # ==================================================
# MAGIC
# MAGIC
# MAGIC display(site_decade_summary)
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC batch_summary_full %>%
# MAGIC   summarise(n_missing = sum(is.na(reps_total)))

# COMMAND ----------

# MAGIC %r
# MAGIC install.packages("cowplot")
# MAGIC library(cowplot)
# MAGIC # -----------------------------
# MAGIC # Remove x-axis from top two
# MAGIC # -----------------------------
# MAGIC ext_plot_nox <- ext_plot +
# MAGIC   theme(
# MAGIC     axis.title.x = element_blank(),
# MAGIC     axis.text.x  = element_blank(),
# MAGIC     axis.ticks.x = element_blank()
# MAGIC   )
# MAGIC
# MAGIC occ_plot_nox <- occ_plot +
# MAGIC   theme(
# MAGIC     axis.title.x = element_blank(),
# MAGIC     axis.text.x  = element_blank(),
# MAGIC     axis.ticks.x = element_blank()
# MAGIC   )
# MAGIC
# MAGIC # -----------------------------
# MAGIC # Extract shared legend
# MAGIC # -----------------------------
# MAGIC shared_legend <- get_legend(
# MAGIC   ext_plot +
# MAGIC     theme(legend.position = "bottom",
# MAGIC           legend.box.margin = margin(t = 10))
# MAGIC )
# MAGIC
# MAGIC # -----------------------------
# MAGIC # Remove legends from individual plots
# MAGIC # -----------------------------
# MAGIC ext_clean  <- ext_plot_nox  + theme(legend.position = "none")
# MAGIC occ_clean  <- occ_plot_nox  + theme(legend.position = "none")
# MAGIC spread_clean <- spread_plot + theme(legend.position = "none")
# MAGIC
# MAGIC # -----------------------------
# MAGIC # Arrange vertically with shared legend
# MAGIC # -----------------------------
# MAGIC final_plot <- plot_grid(
# MAGIC   ext_clean,
# MAGIC   occ_clean,
# MAGIC   spread_clean,
# MAGIC   shared_legend,
# MAGIC   ncol = 1,
# MAGIC   rel_heights = c(1, 1, 1, 0.25)
# MAGIC )
# MAGIC
# MAGIC final_plot
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC library(dplyr)
# MAGIC
# MAGIC metrics <- batch_summary %>%
# MAGIC filter(site == "Kielder" )%>%
# MAGIC   group_by(site, year, Year) %>%
# MAGIC   summarise(
# MAGIC     ext_median = median(p_extinct, na.rm = TRUE),
# MAGIC     ext_lo     = quantile(p_extinct, 0.25, na.rm = TRUE),
# MAGIC     ext_hi     = quantile(p_extinct, 0.75, na.rm = TRUE),
# MAGIC
# MAGIC     occ_median = median(occ_median, na.rm = TRUE),
# MAGIC     occ_lo     = quantile(occ_median, 0.25, na.rm = TRUE),
# MAGIC     occ_hi     = quantile(occ_median, 0.75, na.rm = TRUE),
# MAGIC
# MAGIC     spread_median = median(spread_median, na.rm = TRUE),
# MAGIC     spread_lo     = quantile(spread_median, 0.25, na.rm = TRUE),
# MAGIC     spread_hi     = quantile(spread_median, 0.75, na.rm = TRUE),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC return(metrics)

# COMMAND ----------

# MAGIC %r
# MAGIC library(dplyr)
# MAGIC library(readr)
# MAGIC library(stringr)
# MAGIC library(purrr)
# MAGIC pick_col <- function(df, candidates) {
# MAGIC   cols <- intersect(candidates, names(df))
# MAGIC   if (length(cols)) cols[1] else NA_character_
# MAGIC }
# MAGIC
# MAGIC rep_all <- read.csv("/dbfs/mnt/lab/unrestricted/Fiona/Outputs/ALL_summary_replicates.csv")
# MAGIC nrow(rep_all)
# MAGIC str(rep_all)
# MAGIC # ------------------------------------------------------------
# MAGIC # 1) Identify site directories from rep_all
# MAGIC # ------------------------------------------------------------
# MAGIC inds_col <- pick_col(rep_all, c("init_inds_file", "inds_file", "file"))
# MAGIC if (is.na(inds_col)) stop("No inds file column found in rep_all.")
# MAGIC
# MAGIC site_dirs <- rep_all %>%
# MAGIC   mutate(site_dir = dirname(.data[[inds_col]])) %>%
# MAGIC   distinct(site_dir) %>%
# MAGIC   pull(site_dir)
# MAGIC
# MAGIC # ------------------------------------------------------------
# MAGIC # 2) Read RS_index_summary_*.csv from each site_dir
# MAGIC # ------------------------------------------------------------
# MAGIC rsidx_all <- map_dfr(site_dirs, function(sd) {
# MAGIC   f <- list.files(sd, pattern = "^RS_index_summary_.*\\.csv$", full.names = TRUE)
# MAGIC   if (!length(f)) return(NULL)
# MAGIC
# MAGIC   read_csv(f[1], show_col_types = FALSE) %>%
# MAGIC     mutate(site_dir = sd,
# MAGIC            rsidx_file = f[1])
# MAGIC })
# MAGIC
# MAGIC # ------------------------------------------------------------
# MAGIC # 3) Build rsidx_all2 (summarised per inds file)
# MAGIC # ------------------------------------------------------------
# MAGIC file_col <- pick_col(rsidx_all, c("file", "init_inds_file", "inds_file", "path"))
# MAGIC if (is.na(file_col)) stop("RS-index file column not found.")
# MAGIC
# MAGIC rsidx_all2 <- rsidx_all %>%
# MAGIC   mutate(
# MAGIC     file_basename = basename(.data[[file_col]]),
# MAGIC     year          = as.integer(year),
# MAGIC     raster_value  = if ("raster_value" %in% names(.)) raster_value else NA_character_
# MAGIC   ) %>%
# MAGIC   group_by(site_dir, file_basename, year) %>%
# MAGIC   summarise(
# MAGIC     raster_values = paste(unique(na.omit(raster_value)), collapse = ","),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC cat("RS-index rows:", nrow(rsidx_all2), "\n")
# MAGIC str(rsidx_all2)

# COMMAND ----------

# MAGIC %r
# MAGIC library(dplyr)
# MAGIC library(ggplot2)
# MAGIC library(scales)
# MAGIC library(gridExtra)
# MAGIC library(stringr)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 0) Recode site helper
# MAGIC # ==================================================
# MAGIC recode_site <- function(x) {
# MAGIC   dplyr::recode_factor(
# MAGIC     x,
# MAGIC     "FoD_updated2"     = "FoD",
# MAGIC     "Kielder_updated2" = "Kielder",
# MAGIC     "NYM_updated2"     = "NYM"
# MAGIC   )
# MAGIC }
# MAGIC
# MAGIC # ==================================================
# MAGIC # 1) Load yearly dataset WITH metadata
# MAGIC # ==================================================
# MAGIC
# MAGIC ts_all <- read.csv(
# MAGIC   "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/ALL_pop_metrics_timeseries_with_meta.csv"
# MAGIC )
# MAGIC
# MAGIC ts_sub <- ts_all %>%
# MAGIC   filter(Year >= 0, Year <= 10) %>%
# MAGIC   mutate(
# MAGIC     site = recode_site(site),
# MAGIC     year = factor(year)
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 2) Build habitat lookup from init_inds_file + rsidx_all2
# MAGIC # ==================================================
# MAGIC
# MAGIC hab_lookup <- ts_all %>%
# MAGIC   mutate(
# MAGIC     site = recode_site(site),
# MAGIC     site_dir      = dirname(init_inds_file),
# MAGIC     file_basename = basename(init_inds_file),
# MAGIC     file_year     = as.integer(str_extract(file_basename, "\\d{4}"))
# MAGIC   ) %>%
# MAGIC   left_join(
# MAGIC     rsidx_all2 %>%
# MAGIC       mutate(file_year = as.integer(str_extract(file_basename, "\\d{4}"))) %>%
# MAGIC       select(site_dir, file_basename, file_year, raster_values),
# MAGIC     by = c("site_dir", "file_basename", "file_year")
# MAGIC   ) %>%
# MAGIC   mutate(
# MAGIC     init_code_int = suppressWarnings(as.integer(str_extract(raster_values, "\\d+"))),
# MAGIC     initial_site_short = case_when(
# MAGIC       init_code_int == 10  ~ "Clearfell (managed felling)",
# MAGIC       init_code_int == 100 ~ "Windblow (storm-felled)",
# MAGIC       is.na(init_code_int) ~ NA_character_,
# MAGIC       TRUE                 ~ "Other"
# MAGIC     ),
# MAGIC     initial_site_short = factor(
# MAGIC       initial_site_short,
# MAGIC       levels = c("Clearfell (managed felling)",
# MAGIC                  "Windblow (storm-felled)",
# MAGIC                  "Other")
# MAGIC     )
# MAGIC   ) %>%
# MAGIC   distinct(site, out_subdir, batch, initial_site_short)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 3) Attach habitat to ts_sub
# MAGIC # ==================================================
# MAGIC
# MAGIC ts_sub2 <- ts_sub %>%
# MAGIC   left_join(hab_lookup, by = c("site", "out_subdir", "batch"))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 4) Batch-level summaries
# MAGIC # ==================================================
# MAGIC
# MAGIC batch_summary <- ts_sub2 %>%
# MAGIC   group_by(site, year, out_subdir, batch, Year, initial_site_short) %>%
# MAGIC   summarise(
# MAGIC     occ_median    = median(occupied_cells, na.rm = TRUE),
# MAGIC     spread_median = median(max_spread_m_monotonic, na.rm = TRUE),
# MAGIC     reps_alive  = n_distinct(Rep),
# MAGIC     reps_total  = max(replicates),
# MAGIC     p_extinct   = 1 - reps_alive / reps_total,
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 5) Between-batch summaries (survivors only)
# MAGIC # ==================================================
# MAGIC
# MAGIC between_batch_survivors <- ts_sub2 %>%
# MAGIC   filter(occupied_cells > 0) %>%
# MAGIC   group_by(site, year, out_subdir, Year, initial_site_short) %>%
# MAGIC   summarise(
# MAGIC     occ_med = median(occupied_cells),
# MAGIC     occ_lo  = quantile(occupied_cells, 0.25),
# MAGIC     occ_hi  = quantile(occupied_cells, 0.75),
# MAGIC     spread_med = median(max_spread_m_monotonic),
# MAGIC     spread_lo  = quantile(max_spread_m_monotonic, 0.25),
# MAGIC     spread_hi  = quantile(max_spread_m_monotonic, 0.75),
# MAGIC     n_batches_surviving = n(),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 6) Between-batch extinction (all batches)
# MAGIC # ==================================================
# MAGIC
# MAGIC between_batch_extinction <- batch_summary %>%
# MAGIC   group_by(site, year, out_subdir, Year, initial_site_short) %>%
# MAGIC   summarise(
# MAGIC     ext_med = median(p_extinct),
# MAGIC     ext_lo  = quantile(p_extinct, 0.25),
# MAGIC     ext_hi  = quantile(p_extinct, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 7) Collapse to site × decade × year × habitat
# MAGIC # ==================================================
# MAGIC
# MAGIC site_decade_summary <- between_batch_survivors %>%
# MAGIC   group_by(site, year, Year, initial_site_short) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = median(occ_med),
# MAGIC     occ_lo  = quantile(occ_med, 0.25),
# MAGIC     occ_hi  = quantile(occ_med, 0.75),
# MAGIC     spread_medi = median(spread_med),
# MAGIC     spread_lo  = quantile(spread_med, 0.25),
# MAGIC     spread_hi  = quantile(spread_med, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   ) %>%
# MAGIC   left_join(
# MAGIC     between_batch_extinction %>%
# MAGIC       group_by(site, year, Year, initial_site_short) %>%
# MAGIC       summarise(
# MAGIC         ext_medi = median(ext_med),
# MAGIC         ext_lo  = quantile(ext_med, 0.25),
# MAGIC         ext_hi  = quantile(ext_med, 0.75),
# MAGIC         .groups = "drop"
# MAGIC       ),
# MAGIC     by = c("site", "year", "Year", "initial_site_short")
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 8) PLOTS with habitat as shape
# MAGIC # ==================================================
# MAGIC
# MAGIC shape_vals <- c(16, 17, 15)  # adjust if needed
# MAGIC
# MAGIC ext_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, ext_medi, colour = year, shape = initial_site_short)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, initial_site_short)), size = 0.8) +
# MAGIC   geom_point(size = 3) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = year),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site, nrow = 1) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   scale_y_continuous(labels = percent_format(accuracy = 1)) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Climate decade",
# MAGIC     fill   = "Climate decade",
# MAGIC     shape  = "Initial habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC occ_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, occ_medi, colour = year, shape = initial_site_short)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, initial_site_short)), size = 0.8) +
# MAGIC   geom_point(size = 3) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = year),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site, nrow = 1) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   #scale_y_log10(labels = label_comma()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Climate decade",
# MAGIC     fill   = "Climate decade",
# MAGIC     shape  = "Initial habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC spread_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, spread_medi, colour = year, shape = initial_site_short)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, initial_site_short)), size = 0.8) +
# MAGIC   geom_point(size = 3) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = year),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site, nrow = 1) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance",
# MAGIC     colour = "Climate decade",
# MAGIC     fill   = "Climate decade",
# MAGIC     shape  = "Initial habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 9) Arrange vertically
# MAGIC # ==================================================
# MAGIC
# MAGIC
# MAGIC (ext_plot / occ_plot / spread_plot) + plot_layout(guides = "collect")
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC display(site_decade_summary)

# COMMAND ----------

# MAGIC %r
# MAGIC
# MAGIC # ==================================================
# MAGIC # 1) Average across sites
# MAGIC # ==================================================
# MAGIC
# MAGIC avg_across_sites <- site_decade_summary %>%
# MAGIC   group_by(year, Year, initial_site_short) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = mean(occ_medi, na.rm = TRUE),
# MAGIC     occ_lo  = mean(occ_lo,  na.rm = TRUE),
# MAGIC     occ_hi  = mean(occ_hi,  na.rm = TRUE),
# MAGIC
# MAGIC     spread_medi = mean(spread_medi, na.rm = TRUE),
# MAGIC     spread_lo  = mean(spread_lo,  na.rm = TRUE),
# MAGIC     spread_hi  = mean(spread_hi,  na.rm = TRUE),
# MAGIC
# MAGIC     ext_medi = mean(ext_medi, na.rm = TRUE),
# MAGIC     ext_lo  = mean(ext_lo,  na.rm = TRUE),
# MAGIC     ext_hi  = mean(ext_hi,  na.rm = TRUE),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC shape_vals <- c(16, 17, 15)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 2) Extinction plot (averaged)
# MAGIC # ==================================================
# MAGIC
# MAGIC ext_plot_avg <- ggplot(
# MAGIC   avg_across_sites,
# MAGIC   aes(Year, ext_medi, colour = initial_site_short)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, initial_site_short)), size = 0.9) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = initial_site_short),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   scale_y_continuous(labels = percent_format(accuracy = 1)) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Initial habitat",
# MAGIC     fill   = "Initial habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   facet_wrap(~year)+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 3) Occupancy plot (averaged)
# MAGIC # ==================================================
# MAGIC
# MAGIC occ_plot_avg <- ggplot(
# MAGIC   avg_across_sites,
# MAGIC   aes(Year, occ_medi, colour = initial_site_short)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, initial_site_short)), size = 0.9) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = initial_site_short),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   #scale_y_log10(labels = label_comma()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Initial habitat",
# MAGIC     fill   = "Initial habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   facet_wrap(~year)+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 4) Spread plot (averaged)
# MAGIC # ==================================================
# MAGIC
# MAGIC spread_plot_avg <- ggplot(
# MAGIC   avg_across_sites,
# MAGIC   aes(Year, spread_medi, colour = initial_site_short)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, initial_site_short)), size = 0.9) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = initial_site_short),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Initial habitat",
# MAGIC     fill   = "Initial habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme +
# MAGIC   facet_wrap(~year)+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 5) Arrange vertically
# MAGIC # ==================================================
# MAGIC
# MAGIC grid.arrange(
# MAGIC   ext_plot_avg,
# MAGIC   occ_plot_avg,
# MAGIC   spread_plot_avg,
# MAGIC   nrow = 3
# MAGIC )
# MAGIC
# MAGIC
# MAGIC
# MAGIC (ext_plot_avg / occ_plot_avg / spread_plot_avg) + plot_layout(guides = "collect")
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC
# MAGIC ts_all_mgmt <- ts_all %>%
# MAGIC   mutate(
# MAGIC     decade = factor(year),
# MAGIC     site = factor(site),
# MAGIC
# MAGIC     mgmt_class = case_when(
# MAGIC       apply_destruction_annually == FALSE & years_for_destruction == 10 ~ "No management",
# MAGIC       apply_destruction_annually == FALSE & years_for_destruction == 2  ~ "2-year management",
# MAGIC       apply_destruction_annually == FALSE & years_for_destruction == 5  ~ "5-year management",
# MAGIC       TRUE ~ NA_character_
# MAGIC     )
# MAGIC   ) %>%
# MAGIC   filter(!is.na(mgmt_class))
# MAGIC
# MAGIC library(dplyr)
# MAGIC library(ggplot2)
# MAGIC library(scales)
# MAGIC library(gridExtra)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 1) Filter to management dataset
# MAGIC # ==================================================
# MAGIC
# MAGIC ts_mgmt <- ts_all_mgmt %>%
# MAGIC   filter(Year >= 0, Year <= 10) %>%
# MAGIC   mutate(
# MAGIC     site = recode_site(site),
# MAGIC     year = factor(year),
# MAGIC     mgmt_class = factor(mgmt_class)
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 2) Batch-level summaries
# MAGIC # ==================================================
# MAGIC
# MAGIC batch_mgmt <- ts_mgmt %>%
# MAGIC   group_by(site, year, out_subdir, batch, Year, mgmt_class) %>%
# MAGIC   summarise(
# MAGIC     occ_median    = median(occupied_cells, na.rm = TRUE),
# MAGIC     spread_median = median(max_spread_m_monotonic, na.rm = TRUE),
# MAGIC     reps_alive  = n_distinct(Rep),
# MAGIC     reps_total  = max(replicates),
# MAGIC     p_extinct   = 1 - reps_alive / reps_total,
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 3) Between-batch survivors
# MAGIC # ==================================================
# MAGIC
# MAGIC surv_mgmt <- ts_mgmt %>%
# MAGIC   filter(occupied_cells > 0) %>%
# MAGIC   group_by(site, year, out_subdir, Year, mgmt_class) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = median(occupied_cells),
# MAGIC     occ_lo  = quantile(occupied_cells, 0.25),
# MAGIC     occ_hi  = quantile(occupied_cells, 0.75),
# MAGIC     spread_medi = median(max_spread_m_monotonic),
# MAGIC     spread_lo  = quantile(max_spread_m_monotonic, 0.25),
# MAGIC     spread_hi  = quantile(max_spread_m_monotonic, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 4) Between-batch extinction
# MAGIC # ==================================================
# MAGIC
# MAGIC ext_mgmt <- batch_mgmt %>%
# MAGIC   group_by(site, year, out_subdir, Year, mgmt_class) %>%
# MAGIC   summarise(
# MAGIC     ext_medi = median(p_extinct),
# MAGIC     ext_lo  = quantile(p_extinct, 0.25),
# MAGIC     ext_hi  = quantile(p_extinct, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 5) Collapse to site × decade × year × mgmt_class
# MAGIC # ==================================================
# MAGIC
# MAGIC site_decade_mgmt <- surv_mgmt %>%
# MAGIC   group_by(site, year, Year, mgmt_class) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = median(occ_medi),
# MAGIC     occ_lo  = quantile(occ_lo, 0.25),
# MAGIC     occ_hi  = quantile(occ_hi, 0.75),
# MAGIC     spread_medi = median(spread_medi),
# MAGIC     spread_lo  = quantile(spread_lo, 0.25),
# MAGIC     spread_hi  = quantile(spread_hi, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   ) %>%
# MAGIC   left_join(
# MAGIC     ext_mgmt %>%
# MAGIC       group_by(site, year, Year, mgmt_class) %>%
# MAGIC       summarise(
# MAGIC         ext_medi = median(ext_medi),
# MAGIC         ext_lo  = quantile(ext_medi, 0.25),
# MAGIC         ext_hi  = quantile(ext_medi, 0.75),
# MAGIC         .groups = "drop"
# MAGIC       ),
# MAGIC     by = c("site", "year", "Year", "mgmt_class")
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 6) Average across sites
# MAGIC # ==================================================
# MAGIC
# MAGIC avg_mgmt <- site_decade_mgmt %>%
# MAGIC   group_by(year, Year, mgmt_class) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = mean(occ_medi, na.rm=T),
# MAGIC     occ_lo  = mean(occ_lo, na.rm=T),
# MAGIC     occ_hi  = mean(occ_hi, na.rm=T),
# MAGIC     spread_medi = mean(spread_medi, na.rm=T),
# MAGIC     spread_lo  = mean(spread_lo, na.rm=T),
# MAGIC     spread_hi  = mean(spread_hi, na.rm=T),
# MAGIC     ext_medi = mean(ext_medi, na.rm=T),
# MAGIC     ext_lo  = mean(ext_lo, na.rm=T),
# MAGIC     ext_hi  = mean(ext_hi, na.rm=T),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC shape_vals <- c(16, 17, 15)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 7) Plots (no facets)
# MAGIC # ==================================================
# MAGIC
# MAGIC # ---- Extinction ----
# MAGIC ext_plot_mgmt <- ggplot(
# MAGIC   avg_mgmt,
# MAGIC   aes(Year, ext_medi, colour = mgmt_class)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, mgmt_class)), linewidth = 1.1) +
# MAGIC   geom_point(size = 3) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = mgmt_class),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   scale_y_continuous(labels = percent_format(accuracy = 1)) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinct",
# MAGIC     colour = "Management",
# MAGIC     fill   = "Management"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   facet_wrap(~year)+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ---- Occupancy ----
# MAGIC occ_plot_mgmt <- ggplot(
# MAGIC   avg_mgmt,
# MAGIC   aes(Year, occ_medi, colour = mgmt_class)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, mgmt_class)), linewidth = 1.1) +
# MAGIC   geom_point(size = 3) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = mgmt_class),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC  # scale_y_log10(labels = label_comma()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Management",
# MAGIC     fill   = "Management"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   facet_wrap(~year)+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ---- Spread ----
# MAGIC spread_plot_mgmt <- ggplot(
# MAGIC   avg_mgmt,
# MAGIC   aes(Year, spread_medi, colour = mgmt_class)
# MAGIC ) +
# MAGIC   geom_line(aes(group = interaction(year, mgmt_class)), linewidth = 1.1) +
# MAGIC   geom_point(size = 3) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = mgmt_class),
# MAGIC     alpha = 0.15, colour = NA
# MAGIC   ) +
# MAGIC   scale_shape_manual(values = shape_vals) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Management",
# MAGIC     fill   = "Management"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   the_theme+
# MAGIC   facet_wrap(~year)+
# MAGIC   scale_x_continuous(breaks = c(0,2,4,6,8,10))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 8) Arrange vertically
# MAGIC # ==================================================
# MAGIC
# MAGIC
# MAGIC (ext_plot_mgmt / occ_plot_mgmt / spread_plot_mgmt) + plot_layout(guides = "collect")
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC install.packages("tidyverse")
# MAGIC install.packages("patchwork")
# MAGIC install.packages("scales")
# MAGIC library(dplyr)
# MAGIC library(ggplot2)
# MAGIC library(scales)
# MAGIC library(patchwork)
# MAGIC #library(data.table)
# MAGIC
# MAGIC # ==================================================
# MAGIC # 1) LOAD DATA
# MAGIC # ==================================================
# MAGIC
# MAGIC recode_site <- function(x) {
# MAGIC   dplyr::recode_factor(
# MAGIC     x,
# MAGIC     "FoD_updated2"     = "FoD",
# MAGIC     "Kielder_updated2" = "Kielder",
# MAGIC     "NYM_updated2"     = "NYM"
# MAGIC   )
# MAGIC }
# MAGIC
# MAGIC ts_all <- read.csv(
# MAGIC   "/dbfs/mnt/lab/unrestricted/Fiona/Outputs/ALL_pop_metrics_timeseries_with_meta.csv"
# MAGIC )
# MAGIC
# MAGIC ts_sub <- ts_all %>%
# MAGIC   filter(Year >= 0, Year <= 10) %>%
# MAGIC   mutate(
# MAGIC     site = recode_site(site),
# MAGIC     year = factor(year)
# MAGIC   )
# MAGIC ts_sub <- ts_sub %>%
# MAGIC   mutate(
# MAGIC     Rep_global = paste(batch, Rep, sep = "_")
# MAGIC   )
# MAGIC str(ts_sub)
# MAGIC # ==================================================
# MAGIC # 2) READ + ATTACH HABITAT (SITE + REP LEVEL)
# MAGIC # ==================================================
# MAGIC
# MAGIC init_all <- bind_rows(
# MAGIC   read.csv("/dbfs/mnt/lab/unrestricted/Fiona/FoD_updated2/init_summary_all.csv") %>%
# MAGIC     mutate(site = "FoD_updated2"),
# MAGIC   read.csv("/dbfs/mnt/lab/unrestricted/Fiona/Kielder_updated2/init_summary_all.csv") %>%
# MAGIC     mutate(site = "Kielder_updated2"),
# MAGIC   read.csv("/dbfs/mnt/lab/unrestricted/Fiona/NYM_updated2/init_summary_all.csv") %>%
# MAGIC     mutate(site = "NYM_updated2")
# MAGIC ) %>%
# MAGIC   rename(rep = Rep)
# MAGIC
# MAGIC init_all <- init_all %>%
# MAGIC   mutate(site = recode_site(site)) %>%
# MAGIC   group_by(site, rep) %>%
# MAGIC   slice(1) %>%  # habitat fixed per landscape
# MAGIC   ungroup() %>%
# MAGIC   select(-StartYear)
# MAGIC
# MAGIC ts_sub <- ts_sub %>%
# MAGIC   left_join(init_all, by = c("site", "rep"))
# MAGIC
# MAGIC # ==================================================
# MAGIC # 3) CREATE HABITAT METRICS
# MAGIC # ==================================================
# MAGIC
# MAGIC ts_sub <- ts_sub %>%
# MAGIC   mutate(
# MAGIC     disturbed_p = (Code3 + Code4 + Code5) / TotalCells,
# MAGIC     disturbed_bin = ntile(disturbed_p, 4)
# MAGIC   ) %>%
# MAGIC   mutate(
# MAGIC     disturbed_bin = factor(
# MAGIC       disturbed_bin,
# MAGIC       levels = 1:4,
# MAGIC       labels = c("Very low", "Low", "Medium", "High")
# MAGIC     )
# MAGIC   )
# MAGIC # ==================================================
# MAGIC # 4) BATCH SUMMARY (CORRECT EXTINCTION)
# MAGIC # ==================================================
# MAGIC
# MAGIC batch_summary <- ts_sub %>%
# MAGIC   group_by(site, year, out_subdir, disturbed_bin) %>%
# MAGIC   mutate(
# MAGIC     reps_total = n_distinct(Rep_global[Year == min(Year)])
# MAGIC   ) %>%
# MAGIC   ungroup() %>%
# MAGIC   
# MAGIC   group_by(site, year, out_subdir, Year, disturbed_bin) %>%
# MAGIC   summarise(
# MAGIC     occ_median    = median(occupied_cells, na.rm = TRUE),
# MAGIC     spread_median = median(max_spread_m_monotonic, na.rm = TRUE),
# MAGIC
# MAGIC     reps_alive = n_distinct(Rep_global),
# MAGIC     reps_total = first(reps_total),
# MAGIC
# MAGIC     p_extinct  = 1 - reps_alive / reps_total,
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC # ==================================================
# MAGIC # 5) BETWEEN-BATCH SUMMARIES
# MAGIC # ==================================================
# MAGIC
# MAGIC between_batch_survivors <- ts_sub %>%
# MAGIC   group_by(site, year, out_subdir, Year, disturbed_bin) %>%
# MAGIC   summarise(
# MAGIC     occ_med = median(occupied_cells),
# MAGIC     occ_lo  = quantile(occupied_cells, 0.25),
# MAGIC     occ_hi  = quantile(occupied_cells, 0.75),
# MAGIC
# MAGIC     spread_med = median(max_spread_m_monotonic),
# MAGIC     spread_lo  = quantile(max_spread_m_monotonic, 0.25),
# MAGIC     spread_hi  = quantile(max_spread_m_monotonic, 0.75),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC between_batch_extinction <- batch_summary %>%
# MAGIC   group_by(site, year, out_subdir, Year, disturbed_bin) %>%
# MAGIC   summarise(
# MAGIC     ext_med = median(p_extinct),
# MAGIC     ext_lo  = quantile(p_extinct, 0.25),
# MAGIC     ext_hi  = quantile(p_extinct, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # 6) SITE × DECADE SUMMARY
# MAGIC # ==================================================
# MAGIC # ==================================================
# MAGIC # LABEL DISTURBANCE BINS WITH REAL RANGES
# MAGIC # ==================================================
# MAGIC
# MAGIC # compute ranges once (before summaries)
# MAGIC bin_ranges <- ts_sub %>%
# MAGIC   group_by(disturbed_bin) %>%
# MAGIC   summarise(
# MAGIC     min_p = min(disturbed_p),
# MAGIC     max_p = max(disturbed_p),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC bin_labels <- paste0(
# MAGIC   c("Very low", "Low", "Medium", "High"),
# MAGIC   " (",
# MAGIC   round(bin_ranges$min_p, 2),
# MAGIC   "–",
# MAGIC   round(bin_ranges$max_p, 2),
# MAGIC   ")"
# MAGIC )
# MAGIC
# MAGIC ts_sub$disturbed_bin <- factor(
# MAGIC   ts_sub$disturbed_bin,
# MAGIC   levels = 1:4,
# MAGIC   labels = bin_labels
# MAGIC )
# MAGIC
# MAGIC # ==================================================
# MAGIC # SITE × DECADE SUMMARY (FIXED VARIANCE)
# MAGIC # ==================================================
# MAGIC
# MAGIC site_decade_summary <- between_batch_survivors %>%
# MAGIC   group_by(year, Year, disturbed_bin) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = median(occ_med),
# MAGIC     occ_lo  = min(occ_lo),    # ✅ FIX
# MAGIC     occ_hi  = max(occ_hi),    # ✅ FIX
# MAGIC
# MAGIC     spread_medi = median(spread_med),
# MAGIC     spread_lo  = min(spread_lo),   # ✅ FIX
# MAGIC     spread_hi  = max(spread_hi),   # ✅ FIX
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   ) %>%
# MAGIC   left_join(
# MAGIC     between_batch_extinction %>%
# MAGIC       group_by(year, Year, disturbed_bin) %>%
# MAGIC       summarise(
# MAGIC         ext_medi = median(ext_med),
# MAGIC         ext_lo  = min(ext_lo),  # ✅ FIX
# MAGIC         ext_hi  = max(ext_hi),  # ✅ FIX
# MAGIC         .groups = "drop"
# MAGIC       ),
# MAGIC     by = c("year", "Year", "disturbed_bin")
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # COMMON THEME
# MAGIC # ==================================================
# MAGIC common_theme <- 
# MAGIC   theme(
# MAGIC     legend.position = "right",
# MAGIC     legend.direction = "vertical"
# MAGIC   )
# MAGIC
# MAGIC # ---- EXTINCTION ----
# MAGIC ext_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, ext_medi, colour = disturbed_bin, group = disturbed_bin)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = disturbed_bin),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_grid(~ year) +
# MAGIC   scale_y_continuous(labels = percent_format()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Disturbed habitat proportion"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- OCCUPANCY ----
# MAGIC occ_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, occ_medi, colour = disturbed_bin, group = disturbed_bin)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = disturbed_bin),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_grid(~ year) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Disturbed habitat proportion"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- SPREAD ----
# MAGIC spread_plot <- ggplot(
# MAGIC   site_decade_summary,
# MAGIC   aes(Year, spread_medi, colour = disturbed_bin, group = disturbed_bin)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = disturbed_bin),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_grid(~ year) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Disturbed habitat proportion"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Dark2") +
# MAGIC   scale_fill_brewer(palette = "Dark2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ==================================================
# MAGIC # COMBINE
# MAGIC # ==================================================
# MAGIC
# MAGIC (ext_plot / occ_plot / spread_plot) +
# MAGIC   plot_layout(guides = "collect")

# COMMAND ----------

# MAGIC %r
# MAGIC init_all %>%
# MAGIC   summarise(
# MAGIC     mean_Code2 = mean(Code2),
# MAGIC     mean_Code3 = mean(Code3),
# MAGIC     mean_Code4 = mean(Code4),
# MAGIC     mean_5 = mean(Code5)
# MAGIC   )

# COMMAND ----------

# MAGIC %r
# MAGIC # ==================================================
# MAGIC # CLASSIFY DOMINANT HABITAT (FINAL VERSION)
# MAGIC # ==================================================
# MAGIC common_theme <- theme_minimal() +
# MAGIC   theme(
# MAGIC     panel.grid = element_blank(),        # remove all gridlines
# MAGIC     axis.line = element_line(colour = "black"),
# MAGIC     axis.ticks = element_line(colour = "black"),
# MAGIC     axis.text = element_text(colour = "black"),
# MAGIC     axis.title = element_text(colour = "black"),
# MAGIC     panel.border = element_blank()       # no box unless you want one
# MAGIC   )
# MAGIC
# MAGIC
# MAGIC ts_sub <- ts_sub %>%
# MAGIC   mutate(
# MAGIC     wind    = Code2 + Code5,
# MAGIC     felled  = Code3,
# MAGIC     drought = Code4,
# MAGIC
# MAGIC     max_val = pmax(wind, felled, drought),
# MAGIC
# MAGIC     dominant_habitat = case_when(
# MAGIC       wind    == max_val ~ "Windblown",
# MAGIC       felled  == max_val ~ "Felled",
# MAGIC       drought == max_val ~ "Drought risk"
# MAGIC     )
# MAGIC   )
# MAGIC
# MAGIC all_rep_years <- ts_sub %>%
# MAGIC   distinct(site, year, out_subdir, Rep_global, dominant_habitat) %>%
# MAGIC   tidyr::crossing(Year = 0:10)
# MAGIC
# MAGIC ts_full <- all_rep_years %>%
# MAGIC   left_join(ts_sub,
# MAGIC             by = c("site", "year", "out_subdir",
# MAGIC                    "Rep_global", "Year", "dominant_habitat"))
# MAGIC
# MAGIC ts_full <- ts_full %>%
# MAGIC   mutate(
# MAGIC     alive = !is.na(pop_size)
# MAGIC   )
# MAGIC
# MAGIC
# MAGIC # ==================================================
# MAGIC # BATCH SUMMARY
# MAGIC # ==================================================
# MAGIC
# MAGIC batch_summary_hab <- ts_full %>%
# MAGIC   group_by(site, year, out_subdir, Year, dominant_habitat) %>%
# MAGIC   summarise(
# MAGIC     reps_alive = sum(alive),
# MAGIC     reps_total = n_distinct(Rep_global),
# MAGIC     p_extinct  = 1 - reps_alive / reps_total,
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # BETWEEN-BATCH SUMMARIES
# MAGIC # ==================================================
# MAGIC
# MAGIC between_batch_survivors_hab <- ts_sub %>%
# MAGIC   group_by(site, year, out_subdir, Year, dominant_habitat) %>%
# MAGIC   summarise(
# MAGIC     occ_med = median(occupied_cells),
# MAGIC     occ_lo  = quantile(occupied_cells, 0.25),
# MAGIC     occ_hi  = quantile(occupied_cells, 0.75),
# MAGIC
# MAGIC     spread_med = median(max_spread_m_monotonic),
# MAGIC     spread_lo  = quantile(max_spread_m_monotonic, 0.25),
# MAGIC     spread_hi  = quantile(max_spread_m_monotonic, 0.75),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC between_batch_extinction_hab <- batch_summary_hab %>%
# MAGIC   group_by(site, year, out_subdir, Year, dominant_habitat) %>%
# MAGIC   summarise(
# MAGIC     ext_med = median(p_extinct),
# MAGIC     ext_lo  = quantile(p_extinct, 0.25),
# MAGIC     ext_hi  = quantile(p_extinct, 0.75),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # SITE × DECADE SUMMARY
# MAGIC # ==================================================
# MAGIC
# MAGIC site_decade_summary_hab <- between_batch_survivors_hab %>%
# MAGIC   group_by(year, Year, dominant_habitat) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = mean(occ_med),
# MAGIC     occ_lo  = quantile(occ_med, 0.25),
# MAGIC     occ_hi  = quantile(occ_med, 0.75),
# MAGIC
# MAGIC     spread_medi = mean(spread_med),
# MAGIC     spread_lo  = quantile(spread_med, 0.25),
# MAGIC     spread_hi  = quantile(spread_med, 0.75),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   ) %>%
# MAGIC   left_join(
# MAGIC     between_batch_extinction_hab %>%
# MAGIC       group_by(year, Year, dominant_habitat) %>%
# MAGIC       summarise(
# MAGIC         ext_medi = median(ext_med),
# MAGIC         ext_lo  = quantile(ext_med, 0.25),
# MAGIC         ext_hi  = quantile(ext_med, 0.75),
# MAGIC         .groups = "drop"
# MAGIC       ),
# MAGIC     by = c("year", "Year", "dominant_habitat")
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # PLOTS
# MAGIC # ==================================================
# MAGIC
# MAGIC # ---- EXTINCTION ----
# MAGIC ext_plot_hab <- ggplot(
# MAGIC   site_decade_summary_hab,
# MAGIC   aes(Year, ext_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_grid(~ year) +
# MAGIC   scale_y_continuous(labels = percent_format()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- OCCUPANCY ----
# MAGIC occ_plot_hab <- ggplot(
# MAGIC   site_decade_summary_hab,
# MAGIC   aes(Year, occ_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_grid(~ year) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- SPREAD ----
# MAGIC spread_plot_hab <- ggplot(
# MAGIC   site_decade_summary_hab,
# MAGIC   aes(Year, spread_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_grid(~ year) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Dominant surrounding habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ==================================================
# MAGIC # COMBINE
# MAGIC # ==================================================
# MAGIC
# MAGIC (ext_plot_hab / occ_plot_hab / spread_plot_hab) +
# MAGIC   plot_layout(guides = "collect")
# MAGIC #str(ts_sub)

# COMMAND ----------

# MAGIC %r
# MAGIC # ==================================================
# MAGIC # NEW SITE-LEVEL SUMMARY (KEEPING SITE)
# MAGIC # ==================================================
# MAGIC
# MAGIC site_summary_hab <- between_batch_survivors_hab %>%
# MAGIC   group_by(site, year, Year, dominant_habitat) %>%
# MAGIC   summarise(
# MAGIC     occ_medi = mean(occ_med),
# MAGIC     occ_lo  = quantile(occ_med, 0.25),
# MAGIC     occ_hi  = quantile(occ_med, 0.75),
# MAGIC
# MAGIC     spread_medi = mean(spread_med),
# MAGIC     spread_lo  = quantile(spread_med, 0.25),
# MAGIC     spread_hi  = quantile(spread_med, 0.75),
# MAGIC
# MAGIC     .groups = "drop"
# MAGIC   ) %>%
# MAGIC   left_join(
# MAGIC     between_batch_extinction_hab %>%
# MAGIC       group_by(site, year, Year, dominant_habitat) %>%
# MAGIC       summarise(
# MAGIC         ext_medi = median(ext_med),
# MAGIC         ext_lo  = quantile(ext_med, 0.25),
# MAGIC         ext_hi  = quantile(ext_med, 0.75),
# MAGIC         .groups = "drop"
# MAGIC       ),
# MAGIC     by = c("site", "year", "Year", "dominant_habitat")
# MAGIC   )
# MAGIC
# MAGIC # ==================================================
# MAGIC # NEW PLOTS (FACET BY SITE)
# MAGIC # ==================================================
# MAGIC
# MAGIC # ---- EXTINCTION ----
# MAGIC ext_plot_site <- ggplot(
# MAGIC   site_summary_hab,
# MAGIC   aes(Year, ext_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   scale_y_continuous(labels = percent_format()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- OCCUPANCY ----
# MAGIC occ_plot_site <- ggplot(
# MAGIC   site_summary_hab,
# MAGIC   aes(Year, occ_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- SPREAD ----
# MAGIC spread_plot_site <- ggplot(
# MAGIC   site_summary_hab,
# MAGIC   aes(Year, spread_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Dominant surrounding habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ==================================================
# MAGIC # COMBINE (NEW FIGURE)
# MAGIC # ==================================================
# MAGIC
# MAGIC (ext_plot_site / occ_plot_site / spread_plot_site) +
# MAGIC   plot_layout(guides = "collect")

# COMMAND ----------

# MAGIC %r
# MAGIC # ==================================================
# MAGIC # FILTER TO 2050 ONLY
# MAGIC # ==================================================
# MAGIC
# MAGIC site_summary_2050 <- site_summary_hab %>%
# MAGIC   filter(year == "2070")
# MAGIC
# MAGIC # ==================================================
# MAGIC # PLOTS (FACET = SITE, COLOUR = HABITAT)
# MAGIC # ==================================================
# MAGIC
# MAGIC # ---- EXTINCTION ----
# MAGIC ext_plot_2050 <- ggplot(
# MAGIC   site_summary_2050,
# MAGIC   aes(Year, ext_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = ext_lo, ymax = ext_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   scale_y_continuous(labels = percent_format()) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Probability extinction",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- OCCUPANCY ----
# MAGIC occ_plot_2050 <- ggplot(
# MAGIC   site_summary_2050,
# MAGIC   aes(Year, occ_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = occ_lo, ymax = occ_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Occupied cells",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ---- SPREAD ----
# MAGIC spread_plot_2050 <- ggplot(
# MAGIC   site_summary_2050,
# MAGIC   aes(Year, spread_medi, colour = dominant_habitat, group = dominant_habitat)
# MAGIC ) +
# MAGIC   geom_line(linewidth = 1.2) +
# MAGIC   geom_point(size = 2) +
# MAGIC   geom_ribbon(
# MAGIC     aes(ymin = spread_lo, ymax = spread_hi, fill = dominant_habitat),
# MAGIC     alpha = 0.15,
# MAGIC     colour = NA
# MAGIC   ) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   labs(
# MAGIC     x = "Simulation year",
# MAGIC     y = "Spread distance (m)",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   scale_color_brewer(palette = "Set2") +
# MAGIC   scale_fill_brewer(palette = "Set2") +
# MAGIC   guides(fill = "none") +
# MAGIC   common_theme
# MAGIC
# MAGIC # ==================================================
# MAGIC # COMBINE
# MAGIC # ==================================================
# MAGIC
# MAGIC (ext_plot_2050 / occ_plot_2050 / spread_plot_2050) +
# MAGIC   plot_layout(guides = "collect")

# COMMAND ----------

# MAGIC %r
# MAGIC ggplot(batch_summary,
# MAGIC        aes(Year, p_extinct,
# MAGIC            group = out_subdir,
# MAGIC            colour = out_subdir)) +
# MAGIC   geom_line(alpha = 0.4) +
# MAGIC   facet_grid(site ~ year) +
# MAGIC   labs(
# MAGIC     x = "Year",
# MAGIC     y = "Extinction probability",
# MAGIC     colour = "Batch"
# MAGIC   ) +
# MAGIC   theme_minimal()+ 
# MAGIC   theme(legend.position = "none")
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC
# MAGIC library(dplyr)
# MAGIC library(tidyr)
# MAGIC library(ggplot2)
# MAGIC
# MAGIC # ==============================
# MAGIC # 1) FULL REP × YEAR (KEEP EXTINCT)
# MAGIC # ==============================
# MAGIC
# MAGIC ts_full <- ts_sub %>%
# MAGIC   distinct(site, year, out_subdir, Rep_global) %>%
# MAGIC   crossing(Year = 0:10) %>%
# MAGIC   left_join(ts_sub,
# MAGIC             by = c("site", "year", "out_subdir", "Rep_global", "Year")) %>%
# MAGIC   mutate(
# MAGIC     occupied_cells = ifelse(is.na(occupied_cells), 0, occupied_cells),
# MAGIC     max_spread_m_monotonic = ifelse(is.na(max_spread_m_monotonic), 0,
# MAGIC                                     max_spread_m_monotonic)
# MAGIC   )
# MAGIC
# MAGIC # ==============================
# MAGIC # 2) BASELINE DISTURBANCE (FIXED)
# MAGIC # ==============================
# MAGIC
# MAGIC ts_full <- ts_full %>%
# MAGIC   group_by(Rep_global) %>%
# MAGIC   arrange(Year) %>%
# MAGIC   mutate(
# MAGIC     disturbed_p_baseline = first(disturbed_p)
# MAGIC   ) %>%
# MAGIC   ungroup()
# MAGIC
# MAGIC # ==============================
# MAGIC # 3) EXTRACT YEAR 10
# MAGIC # ==============================
# MAGIC
# MAGIC spread_10 <- ts_full %>%
# MAGIC   filter(Year == 10)
# MAGIC
# MAGIC # ==============================
# MAGIC # 4) PLOT (CONTINUOUS → BOXPLOT)
# MAGIC # ==============================
# MAGIC
# MAGIC ggplot(spread_10,
# MAGIC        aes(x = year,
# MAGIC            y = max_spread_m_monotonic,
# MAGIC            fill = year)) +
# MAGIC   geom_boxplot(outlier.alpha = 0.2) +
# MAGIC   facet_wrap(~ site) +
# MAGIC   labs(
# MAGIC     x = "Climate decade",
# MAGIC     y = "Spread distance at Year 10 (m)",
# MAGIC     fill = "Climate decade"
# MAGIC   ) +
# MAGIC   theme_minimal()
# MAGIC
# MAGIC ggplot(spread_10,
# MAGIC        aes(x = disturbed_p_baseline,
# MAGIC            y = max_spread_m_monotonic)) +
# MAGIC   geom_point(alpha = 0.2) +
# MAGIC   geom_smooth() +
# MAGIC   #facet_wrap(~ site + year) +
# MAGIC   labs(
# MAGIC     x = "Disturbed habitat proportion (baseline)",
# MAGIC     y = "Spread distance at Year 10 (m)"
# MAGIC   ) +
# MAGIC   theme_minimal()
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC library(dplyr)
# MAGIC library(tidyr)
# MAGIC library(ggplot2)
# MAGIC
# MAGIC ts_sub %>%
# MAGIC   # get FIRST incursion per incursion point
# MAGIC   group_by(init_inds_file) %>%
# MAGIC   arrange(Year) %>%
# MAGIC   filter(occupied_cells > 0) %>%
# MAGIC   slice(1) %>%
# MAGIC   ungroup() %>%
# MAGIC   
# MAGIC   # keep disturbed codes only
# MAGIC   select(site, year, init_inds_file,
# MAGIC          Code2, Code3, Code4, Code5) %>%
# MAGIC   
# MAGIC   # reshape long
# MAGIC   pivot_longer(
# MAGIC     cols = starts_with("Code"),
# MAGIC     names_to = "code",
# MAGIC     values_to = "value"
# MAGIC   ) %>%
# MAGIC   
# MAGIC   # convert to proportions
# MAGIC   group_by(init_inds_file) %>%
# MAGIC   mutate(
# MAGIC     prop = value / sum(value)
# MAGIC   ) %>%
# MAGIC   ungroup() %>%
# MAGIC   
# MAGIC   # plot stacked bars
# MAGIC   ggplot(aes(x = init_inds_file, y = prop, fill = code)) +
# MAGIC   geom_col() +
# MAGIC   
# MAGIC   facet_grid(site ~ year, scales = "free_x") +
# MAGIC   
# MAGIC   labs(
# MAGIC     x = "Incursion point",
# MAGIC     y = "Proportion of disturbed habitat",
# MAGIC     fill = "Code"
# MAGIC   ) +
# MAGIC   
# MAGIC   theme_minimal() +
# MAGIC   theme(
# MAGIC     axis.text.x = element_blank(),
# MAGIC     axis.ticks.x = element_blank()
# MAGIC   )

# COMMAND ----------

# MAGIC %r
# MAGIC library(dplyr)
# MAGIC library(tidyr)
# MAGIC
# MAGIC # ==============================
# MAGIC # 1) DEFINE HABITAT PER INCURSION POINT
# MAGIC # ==============================
# MAGIC
# MAGIC habitat_incursion <- ts_sub %>%
# MAGIC   group_by(init_inds_file) %>%
# MAGIC   arrange(Year) %>%
# MAGIC   slice(1) %>%   # baseline conditions
# MAGIC   ungroup() %>%
# MAGIC   mutate(
# MAGIC     wind    = Code2 + Code5,
# MAGIC     felled  = Code3,
# MAGIC     drought = Code4,
# MAGIC     dominant_habitat = case_when(
# MAGIC       wind    >= felled & wind    >= drought ~ "Windblown",
# MAGIC       felled  >= wind   & felled >= drought ~ "Felled",
# MAGIC       drought >= wind   & drought >= felled ~ "Drought risk"
# MAGIC     )
# MAGIC   ) %>%
# MAGIC   select(init_inds_file, dominant_habitat)
# MAGIC
# MAGIC # ==============================
# MAGIC # 2) JOIN BACK (FIX .x/.y ISSUE)
# MAGIC # ==============================
# MAGIC
# MAGIC ts_with_hab <- ts_sub %>%
# MAGIC   left_join(habitat_incursion, by = "init_inds_file") %>%
# MAGIC   mutate(
# MAGIC     dominant_habitat = coalesce(dominant_habitat.y,
# MAGIC                                 dominant_habitat.x)
# MAGIC   ) %>%
# MAGIC   select(-dominant_habitat.x, -dominant_habitat.y)
# MAGIC
# MAGIC # ==============================
# MAGIC # 3) FORCE COMPLETE TRAJECTORY + GET TRUE YEAR 10
# MAGIC # ==============================
# MAGIC
# MAGIC spread_10 <- ts_with_hab %>%
# MAGIC   group_by(Rep_global, site, year, dominant_habitat) %>%
# MAGIC   complete(Year = 0:10) %>%
# MAGIC   ungroup() %>%
# MAGIC   filter(Year == 9) %>%
# MAGIC   mutate(
# MAGIC     max_spread_m_monotonic = coalesce(max_spread_m_monotonic, 0)
# MAGIC   )
# MAGIC
# MAGIC # ==============================
# MAGIC # 4) SUMMARISE BY HABITAT
# MAGIC # ==============================
# MAGIC
# MAGIC spread_summary <- spread_10 %>%
# MAGIC   group_by(year, dominant_habitat, site) %>%
# MAGIC   summarise(
# MAGIC     spread_med = median(max_spread_m_monotonic, na.rm = TRUE),
# MAGIC     spread_lo  = quantile(max_spread_m_monotonic, 0.25, na.rm = TRUE),
# MAGIC     spread_hi  = quantile(max_spread_m_monotonic, 0.75, na.rm = TRUE),
# MAGIC     .groups = "drop"
# MAGIC   )
# MAGIC
# MAGIC # ==============================
# MAGIC # 5) PLOT
# MAGIC # ==============================
# MAGIC
# MAGIC pos <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)
# MAGIC
# MAGIC ggplot(spread_summary,
# MAGIC        aes(x = year,
# MAGIC            y = spread_med / 10,
# MAGIC            colour = dominant_habitat,
# MAGIC            group = dominant_habitat)) +
# MAGIC   
# MAGIC   geom_point(size = 2, position = pos) +
# MAGIC   geom_errorbar(aes(ymin = spread_lo / 10,
# MAGIC                     ymax = spread_hi / 10),
# MAGIC                 width = 0.2,
# MAGIC                 position = pos) +
# MAGIC   facet_wrap(~site) +
# MAGIC   
# MAGIC   labs(
# MAGIC     x = "Climate decade",
# MAGIC     y = "Spread rate (m/year)",
# MAGIC     colour = "Dominant habitat"
# MAGIC   ) +
# MAGIC   
# MAGIC   theme_minimal() +
# MAGIC   theme(
# MAGIC     panel.grid = element_blank(),
# MAGIC     axis.line  = element_line(colour = "black"),
# MAGIC     axis.ticks = element_line(colour = "black"),
# MAGIC     axis.text  = element_text(colour = "black"),
# MAGIC     axis.title = element_text(colour = "black")
# MAGIC   )
# MAGIC

# COMMAND ----------

# MAGIC %r
# MAGIC library(dplyr)
# MAGIC library(tidyr)
# MAGIC library(ggplot2)
# MAGIC
# MAGIC # ==============================
# MAGIC # 1) GET BASELINE HABITAT PROPORTIONS PER INCURSION
# MAGIC # ==============================
# MAGIC
# MAGIC habitat_props <- ts_sub %>%
# MAGIC   group_by(init_inds_file) %>%
# MAGIC   arrange(Year) %>%
# MAGIC   slice(1) %>%   # baseline
# MAGIC   ungroup() %>%
# MAGIC   mutate(
# MAGIC     wind    = Code2 + Code5,
# MAGIC     felled  = Code3,
# MAGIC     drought = Code4,
# MAGIC     total   = wind + felled + drought,
# MAGIC     
# MAGIC     wind_p    = wind / total,
# MAGIC     felled_p  = felled / total,
# MAGIC     drought_p = drought / total
# MAGIC   ) %>%
# MAGIC   select(init_inds_file, wind_p, felled_p, drought_p)
# MAGIC
# MAGIC # ==============================
# MAGIC # 2) JOIN TO DATA + GET FINAL SPREAD
# MAGIC # ==============================
# MAGIC
# MAGIC spread_10 <- ts_sub %>%
# MAGIC   left_join(habitat_props, by = "init_inds_file") %>%
# MAGIC   group_by(Rep_global) %>%
# MAGIC   arrange(Year) %>%
# MAGIC   slice_tail(n = 1) %>%   # Year 10
# MAGIC   ungroup() %>%
# MAGIC   mutate(
# MAGIC     spread_rate = ifelse(is.na(max_spread_m_monotonic), 0,
# MAGIC                          max_spread_m_monotonic / 10) # per year
# MAGIC   )
# MAGIC
# MAGIC # ==============================
# MAGIC # 3) LONG FORMAT FOR FACETING BY HABITAT TYPE
# MAGIC # ==============================
# MAGIC
# MAGIC spread_long <- spread_10 %>%
# MAGIC   pivot_longer(
# MAGIC     cols = c(wind_p, felled_p, drought_p),
# MAGIC     names_to = "habitat",
# MAGIC     values_to = "prop"
# MAGIC   )
# MAGIC
# MAGIC # ==============================
# MAGIC # 4) PLOT
# MAGIC # ==============================
# MAGIC
# MAGIC ggplot(spread_long,
# MAGIC        aes(prop, spread_rate, colour = year)) +
# MAGIC   
# MAGIC   geom_point(alpha = 0.25) +
# MAGIC   geom_line() +
# MAGIC   
# MAGIC   facet_grid(Year ~ habitat) +
# MAGIC   
# MAGIC   labs(
# MAGIC     x = "Proportion of habitat type (baseline)",
# MAGIC     y = "Spread rate (m/year)",
# MAGIC     colour = "Climate decade"
# MAGIC   ) +
# MAGIC   
# MAGIC   theme_minimal() +
# MAGIC   theme(
# MAGIC     panel.grid = element_blank(),
# MAGIC     axis.line = element_line(colour = "black"),
# MAGIC     axis.ticks = element_line(colour = "black"),
# MAGIC     axis.text = element_text(colour = "black"),
# MAGIC     axis.title = element_text(colour = "black")
# MAGIC   )
# MAGIC
