
library(terra)
library(dplyr)
library(readr)
library(tidyr)

species_init_map_fast <- function(
    Output_path,
    start_years,
    NumberCells = 1,          # keep as 1 for "one point per rep"
    reps = 1000,              # requested points per year
    buffer_dist = 500,
    seed = 42,
    write_inds = FALSE,       # only write inds later for selected medoids
    reps_policy = c("cap","replace","error"),  # behaviour when reps > available
    place_on_nonhabitat = FALSE                 # TRUE = allow NA habitat for placement
) {
  set.seed(seed)
  reps_policy <- match.arg(reps_policy)
  all_summary <- dplyr::tibble()
  
  for (start_year in start_years) {
    message("=== Year ", start_year, " ===")
    
    quality_file  <- file.path(Output_path, paste0(start_year, "_quality_10yWindFell.tif"))
    habitat_file  <- file.path(Output_path, paste0(start_year, "_habcodes_fell.asc"))
    
    if (!file.exists(quality_file) || !file.exists(habitat_file)) {
      warning("Missing files for year: ", start_year)
      next
    }
    
    Map_rast     <- terra::rast(quality_file)
    habitat_rast <- terra::rast(habitat_file)
    
    # Align geometries if needed
    if (!terra::compareGeom(habitat_rast, Map_rast, stopOnError = FALSE)) {
      habitat_rast <- terra::project(habitat_rast, Map_rast)
    }
    
    # Candidate cells for placing points
    qv <- terra::values(Map_rast, mat = FALSE)
    hv <- terra::values(habitat_rast, mat = FALSE)
    
    if (isTRUE(place_on_nonhabitat)) {
      # Any cell with quality > 0 is allowed (habitat NA considered valid)
      valid <- which(qv > 0 & !is.na(qv))
    } else {
      # Only cells with quality > 0 AND a known habitat code (1..6)
      valid <- which(qv > 0 & !is.na(hv))
    }
    
    n_avail <- length(valid)
    message("Available candidate cells in ", start_year, ": ", n_avail)
    
    if (n_avail == 0) {
      warning("No valid cells for year ", start_year, ". Skipping.")
      next
    }
    
    # Decide how many to take and whether to allow replacement
    if (reps_policy == "error" && n_avail < reps) {
      stop("Not enough valid cells (", n_avail, ") for reps = ", reps, " in year ", start_year)
    }
    
    if (reps_policy == "cap" && n_avail < reps) {
      warning("Year ", start_year, ": requested reps = ", reps, 
              " but only ", n_avail, " available. Capping to ", n_avail, ".")
      reps_eff <- n_avail
      replace_eff <- FALSE
    } else if (reps_policy == "replace" && n_avail < reps) {
      warning("Year ", start_year, ": requested reps = ", reps, 
              " > available ", n_avail, ". Sampling with replacement.")
      reps_eff <- reps
      replace_eff <- TRUE
    } else {
      reps_eff <- reps
      replace_eff <- FALSE
    }
    
    # Sample cells
    chosen <- sample(valid, reps_eff, replace = replace_eff)
    
    # Coordinates for chosen cells
    xy <- terra::xyFromCell(Map_rast, chosen)
    
    # Points
    init_points <- terra::vect(xy, type = "points", crs = terra::crs(Map_rast))
    
    # Buffers once for all points
    buffers <- terra::buffer(init_points, width = buffer_dist)
    
    # ---- Extract habitat within buffers (single call) ----
    ext <- terra::extract(habitat_rast, buffers, cells = TRUE)  # cols: ID, cell, <hab layer>
    hab_col <- setdiff(names(ext), c("ID", "cell"))[1]
    if (is.na(hab_col)) stop("Could not identify habitat column from extract() output.")
    
    ext2 <- ext %>% dplyr::mutate(hab = .data[[hab_col]])
    
    # Totals & NonHabitat
    tot_df <- ext2 %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(
        TotalCells = dplyr::n(),
        NonHabitat = sum(is.na(.data$hab)),
        .groups = "drop"
      )
    
    # Habitat code counts (1..6)
    freq_df <- ext2 %>%
      dplyr::filter(!is.na(.data$hab)) %>%
      dplyr::count(.data$ID, .data$hab, name = "n_cells")
    
    code_counts <- tidyr::pivot_wider(
      freq_df,
      id_cols = ID,
      names_from = hab,
      values_from = n_cells,
      values_fill = 0,
      names_prefix = "Code"
    )
    for (k in 1:6) {
      nm <- paste0("Code", k)
      if (!nm %in% names(code_counts)) code_counts[[nm]] <- 0L
    }
    
    coords_df <- tibble::tibble(
      ID = seq_len(nrow(xy)),
      X  = xy[, 1],
      Y  = xy[, 2]
    )
    
    sum_y <- coords_df %>%
      dplyr::left_join(code_counts, by = "ID") %>%
      dplyr::left_join(tot_df,      by = "ID") %>%
      dplyr::mutate(
        StartYear = start_year,
        Rep       = seq_len(dplyr::n()),
        InitSite  = 1L
      ) %>%
      dplyr::select(
        StartYear, Rep, InitSite, X, Y,
        dplyr::all_of(paste0("Code", 1:6)),
        NonHabitat, TotalCells
      )
    
    all_summary <- dplyr::bind_rows(all_summary, sum_y)
    
    # We keep write_inds = FALSE here; you'll write inds only for selected medoids later
  }
  
  out_csv <- file.path(Output_path, "init_summary_all.csv")
  readr::write_csv(all_summary, out_csv)
  message("Saved: ", out_csv)
  
  invisible(all_summary)
}




library(dplyr)
library(readr)
library(cluster)
library(ggplot2)
library(uwot)
library(tidyr)

# --------- Selection: per-year ASW-optimal PAM medoids ----------


select_representatives_by_year <- function(
    summary_csv,
    out_dir = dirname(summary_csv),
    use_xy   = TRUE,
    use_props = TRUE,
    include_nonhab = TRUE,
    w_codes = rep(1, 6),
    w_nonhab = 1,
    w_xy = c(0.5, 0.5),
    k_min = 2, k_max = 20,
    asw_tol = 0.02,
    make_plots = TRUE,
    seed = 42
) {
  set.seed(seed)
  df <- readr::read_csv(summary_csv, show_col_types = FALSE)
  
  years <- sort(unique(df$StartYear))
  all_reps <- list()
  year_summ <- list()
  all_membership <- list()  # NEW: store cluster assignment for ALL points
  
  for (yr in years) {
    message("---- Year ", yr, " ----")
    df_y <- df %>% dplyr::filter(StartYear == yr)
    n_y <- nrow(df_y)
    if (n_y < 2) {
      reps <- df_y
      reps$ChosenCluster <- 1L
      all_reps[[as.character(yr)]] <- reps
      year_summ[[as.character(yr)]] <- data.frame(StartYear = yr, n = n_y, k_star = 1, asw_max = NA, coverage = NA)
      
      # trivial membership
      all_membership[[as.character(yr)]] <- df_y %>%
        dplyr::transmute(StartYear = yr, Rep, ChosenCluster = 1L)
      next
    }
    
    # Build features (counts OR proportions + optional NonHab + optional XY)
    if (!"NonHabitat" %in% names(df_y)) {
      df_y <- df_y %>%
        dplyr::mutate(NonHabitat = pmax(0, TotalCells - (Code1+Code2+Code3+Code4+Code5+Code6)))
    }
    
    if (isTRUE(use_props)) {
      eps <- 1e-9
      for (k in 1:6) df_y[[paste0("Code", k, "_p")]] <- df_y[[paste0("Code", k)]] / pmax(df_y$TotalCells, eps)
      df_y$NonHabitat_p <- df_y$NonHabitat / pmax(df_y$TotalCells, eps)
      code_cols <- paste0("Code", 1:6, "_p")
      nonhab_col <- "NonHabitat_p"
    } else {
      code_cols <- paste0("Code", 1:6)
      nonhab_col <- "NonHabitat"
    }
    
    feature_cols <- code_cols
    if (isTRUE(include_nonhab)) feature_cols <- c(feature_cols, nonhab_col)
    if (isTRUE(use_xy)) feature_cols <- c(feature_cols, "X","Y")
    
    X <- df_y %>%
      dplyr::select(dplyr::all_of(feature_cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
    
    # Drop zero-variance columns
    zero_var <- sapply(X, function(z) length(unique(z[!is.na(z)])) <= 1)
    X <- X[, !zero_var, drop = FALSE]
    
    # Weights
    w_map <- c(setNames(w_codes, paste0("Code", 1:6)))
    if (isTRUE(include_nonhab)) {
      if (isTRUE(use_props)) {
        w_map["NonHabitat_p"] <- w_nonhab
      } else {
        w_map["NonHabitat"] <- w_nonhab
      }
    }
    if (isTRUE(use_xy)) {
      w_map[c("X", "Y")] <- w_xy
    }
    w_vec <- as.numeric(w_map[colnames(X)])
    w_vec[is.na(w_vec)] <- 1
    W <- diag(w_vec, nrow = length(w_vec))
    
    # Impute NAs (median), weight + scale
    for (j in seq_along(X)) {
      v <- X[[j]]
      if (anyNA(v)) {
        med <- stats::median(v, na.rm = TRUE)
        if (is.na(med)) med <- 0
        v[is.na(v)] <- med
        X[[j]] <- v
      }
    }
    M_w <- as.matrix(X) %*% W
    M_scaled <- scale(M_w)
    stopifnot(!any(is.na(M_scaled)), nrow(M_scaled) >= 2)
    
    n <- nrow(M_scaled)
    k_max_y <- min(k_max, n - 1)
    if (k_max_y < k_min) {
      D <- as.matrix(dist(M_scaled))
      center_idx <- which.min(rowMeans(D))
      reps <- df_y[center_idx, ]
      reps$ChosenCluster <- 1L
      all_reps[[as.character(yr)]] <- reps
      year_summ[[as.character(yr)]] <- data.frame(StartYear = yr, n = n, k_star = 1, asw_max = NA, coverage = NA)
      # membership (all points -> cluster 1)
      all_membership[[as.character(yr)]] <- df_y %>% dplyr::transmute(StartYear = yr, Rep, ChosenCluster = 1L)
      next
    }
    
    # Evaluate k by ASW
    coverage_metric <- function(Ms, idx) {
      D <- as.matrix(dist(Ms))
      mean(apply(D[, idx, drop = FALSE], 1, min))
    }
    
    eval_k <- data.frame(k = integer(), asw = numeric(), coverage = numeric())
    for (k in k_min:k_max_y) {
      fit <- cluster::pam(M_scaled, k = k, metric = "euclidean")
      eval_k <- rbind(eval_k, data.frame(
        k = k,
        asw = fit$silinfo$avg.width,
        coverage = coverage_metric(M_scaled, fit$id.med)
      ))
    }
    asw_max <- max(eval_k$asw, na.rm = TRUE)
    cand_k  <- eval_k$k[eval_k$asw >= (asw_max - asw_tol)]
    k_star  <- min(cand_k)
    
    pam_star <- cluster::pam(M_scaled, k = k_star, metric = "euclidean")
    idx_star <- pam_star$id.med
    clusters <- pam_star$clustering
    cov_star <- coverage_metric(M_scaled, idx_star)
    
    # Representatives (medoids)
    reps <- df_y[idx_star, ] %>% dplyr::mutate(ChosenCluster = seq_len(nrow(.)))
    all_reps[[as.character(yr)]] <- reps
    year_summ[[as.character(yr)]] <- data.frame(StartYear = yr, n = n, k_star = k_star, asw_max = asw_max, coverage = cov_star)
    
    # NEW: Save membership for all points for this StartYear
    memb_y <- df_y %>%
      dplyr::mutate(ChosenCluster = as.integer(clusters)) %>%
      dplyr::select(StartYear, Rep, ChosenCluster)
    all_membership[[as.character(yr)]] <- memb_y
    
    # Plots (unchanged)
    if (isTRUE(make_plots)) {
      p_asw <- ggplot2::ggplot(eval_k, ggplot2::aes(k, asw)) +
        ggplot2::geom_line(color = "#0072B2") + ggplot2::geom_point(color = "#0072B2") +
        ggplot2::geom_vline(xintercept = k_star, linetype = 2, color = "#D55E00") +
        ggplot2::labs(title = paste0("ASW vs k — ", yr),
                      subtitle = paste0("Chosen k* = ", k_star, " (tol=", asw_tol, ")")) +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(file.path(out_dir, paste0("ASW_vs_k_", yr, ".png")), p_asw, width = 7, height = 5, dpi = 150)
      
      p_cov <- ggplot2::ggplot(eval_k, ggplot2::aes(k, coverage)) +
        ggplot2::geom_line(color = "#009E73") + ggplot2::geom_point(color = "#009E73") +
        ggplot2::geom_vline(xintercept = k_star, linetype = 2, color = "#D55E00") +
        ggplot2::labs(title = paste0("Coverage vs k — ", yr), y = "mean min distance") +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(file.path(out_dir, paste0("Coverage_vs_k_", yr, ".png")), p_cov, width = 7, height = 5, dpi = 150)
      
      n_neighbors <- min(max(10, round(sqrt(n) + 10)), max(2, n - 1))
      emb <- uwot::umap(M_scaled, n_neighbors = n_neighbors, min_dist = 0.2,
                        metric = "euclidean", n_components = 2, verbose = FALSE)
      umap_df <- as.data.frame(emb)
      if (ncol(umap_df) == 2) names(umap_df) <- c("UMAP1","UMAP2")
      p_umap <- ggplot2::ggplot(cbind(umap_df, cluster = factor(clusters)),
                                ggplot2::aes(UMAP1, UMAP2, color = cluster)) +
        ggplot2::geom_point(size = 1.6, alpha = 0.9) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(title = paste0("UMAP — ", yr, " (k*=", k_star, ")"))
      ggplot2::ggsave(file.path(out_dir, paste0("UMAP_", yr, ".png")), p_umap, width = 7, height = 5, dpi = 150)
    }
  }
  
  reps_combined   <- dplyr::bind_rows(all_reps)
  summ            <- dplyr::bind_rows(year_summ)
  membership_all  <- dplyr::bind_rows(all_membership)
  
  # Write your existing files
  readr::write_csv(reps_combined, file.path(out_dir, "subset_PAM_medoids_by_year.csv"))
  readr::write_csv(summ,          file.path(out_dir, "summary_by_year.csv"))
  
  # NEW: write membership (all points)
  readr::write_csv(membership_all, file.path(out_dir, "cluster_membership_all_points.csv"))
  
  # NEW: also write the full summary file augmented with ChosenCluster (handy!)
  df_with_clusters <- df %>%
    dplyr::left_join(membership_all, by = c("StartYear","Rep"))
  readr::write_csv(df_with_clusters, file.path(out_dir, "init_summary_all_with_clusters.csv"))
  
  invisible(list(
    reps = reps_combined,
    summary = summ,
    membership = membership_all,
    full_with_clusters = df_with_clusters
  ))
}




out_path <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/Kielder/"


summary_df <- species_init_map_fast(
  Output_path         = out_path,
  start_years         = c(2020, 2030, 2040),
  NumberCells         = 1,
  reps                = 1000,
  buffer_dist         = 500,
  seed                = 42,
  write_inds          = FALSE,
  reps_policy         = "cap",           # or "replace" if you MUST hit 1000
  place_on_nonhabitat = FALSE            # set TRUE to allow NA habitat placement
)


# Step 2: select a per-year shortlist (ASW-optimal k, 1 medoid per cluster)
sel <- select_representatives_by_year(
  summary_csv     = file.path(out_path, "init_summary_all.csv"),
  out_dir         = out_path,
  use_xy          = TRUE,     # include spatial spread
  use_props       = TRUE,     # composition (proportions) recommended
  include_nonhab  = TRUE,     # include NonHabitat if it varies
  w_codes         = c(1, 1, 1, 1, 1, 1), # tweak if some codes should be down/up-weighted
  w_nonhab        = 1,
  w_xy            = c(0.5, 0.5),
  k_min           = 2,
  k_max           = 20,
  asw_tol         = 0.02,
  make_plots      = TRUE
)





# Writes inds_*.txt (and optionally initial_dist_*.tif) for the selected medoids only
# Assumes the sampling/summary used the same raster stack (same CRS/res) and saved X/Y as map coordinates.
write_inds_for_selected <- function(
    selected_csv,               # e.g., file.path(out_path, "subset_PAM_medoids_by_year.csv")
    Output_path,                # folder where year rasters live AND where to write outputs
    quality_suffix = "_quality_10yWindFell.tif",  # to load the raster per year
    pop_size = 2,               # total individuals per selected site (split F/M equally)
    include_cluster_in_name = TRUE, # name inds files with cluster id for clarity
    also_write_rasters = FALSE       # write initial_dist_*.tif for selected only (optional)
) {
  # Dependencies
  requireNamespace("terra", quietly = TRUE)
  requireNamespace("readr", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  reps <- readr::read_csv(selected_csv, show_col_types = FALSE)
  
  # Check fields we rely on
  needed <- c("StartYear", "X", "Y")
  missing <- setdiff(needed, names(reps))
  if (length(missing) > 0) stop("Selected CSV missing columns: ", paste(missing, collapse = ", "))
  
  # Optional columns
  has_rep_col <- "Rep" %in% names(reps)
  has_cluster <- "ChosenCluster" %in% names(reps)
  
  # Split by year
  years <- sort(unique(reps$StartYear))
  
  for (yr in years) {
    message("Writing inds for StartYear = ", yr, " ...")
    
    # Load the year’s base raster to get grid geometry
    quality_file <- file.path(Output_path, paste0(yr, quality_suffix))
    if (!file.exists(quality_file)) {
      warning("Quality raster not found for year ", yr, ": ", quality_file, 
              ". Skipping inds for this year.")
      next
    }
    Map_rast <- terra::rast(quality_file)
    nr <- terra::nrow(Map_rast); nc <- terra::ncol(Map_rast)
    
    # Filter rows for this year
    reps_y <- reps[reps$StartYear == yr, , drop = FALSE]
    
    # Convert map coordinates back to cell numbers
    cells <- terra::cellFromXY(Map_rast, reps_y[, c("X","Y")])
    if (any(is.na(cells))) {
      # Fallback: try snapping to nearest cell by small jitter (rarely needed)
      warning("Some selected points returned NA from cellFromXY; check CRS/coordinates.")
    }
    
    # Compute 0-based grid indices (to match your previous inds format)
    col0 <- (cells - 1) %% nc
    row0 <- (cells - 1) %/% nc
    X_ind <- col0
    Y_ind <- (nr - 1) - row0   # origin at bottom-left, as you used before
    
    # Prepare filenames and write
    for (i in seq_len(nrow(reps_y))) {
      rid <- if (has_rep_col) reps_y$Rep[i] else i
      cl  <- if (has_cluster) reps_y$ChosenCluster[i] else NA_integer_
      
      # Build file name
      if (include_cluster_in_name && has_cluster) {
        inds_name <- sprintf("inds_%d_cluster%d_rep%d.txt", yr, cl, rid)
      } else {
        inds_name <- sprintf("inds_%d_rep%d.txt", yr, rid)
      }
      
      # Female & Male rows (split pop_size evenly)
      inds <- data.frame(
        Year    = c(0, 0),
        Species = c(0, 0),
        X       = c(X_ind[i], X_ind[i]),
        Y       = c(Y_ind[i], Y_ind[i]),
        Ninds   = c(pop_size/2, pop_size/2),
        Sex     = c(1, 0)
      )
      
      inds_file <- file.path(Output_path, inds_name)
      utils::write.table(inds, inds_file, row.names = FALSE, quote = FALSE, sep = "\t")
      
      # Optional: write single-cell raster for the selected site
      if (isTRUE(also_write_rasters)) {
        rst <- Map_rast
        terra::values(rst) <- NA
        terra::values(rst)[cells[i]] <- 1
        names(rst) <- "pop"
        # Include cluster in the name if available
        if (include_cluster_in_name && has_cluster) {
          rst_name <- sprintf("initial_dist_%d_cluster%d_rep%d.tif", yr, cl, rid)
        } else {
          rst_name <- sprintf("initial_dist_%d_rep%d.tif", yr, rid)
        }
        terra::writeRaster(rst, file.path(Output_path, rst_name), overwrite = TRUE, datatype = "INT2S")
      }
    }
  }
  
  message("Done writing inds (and rasters if requested).")
}


selected_csv <- paste(out_path, "subset_PAM_medoids_by_year.csv", sep = "")

# Write inds (and optionally rasters) ONLY for selected medoids
write_inds_for_selected(
  selected_csv        = selected_csv,
  Output_path         = out_path,
  quality_suffix      = "_quality_10yWindFell.tif",
  pop_size            = 2,
  include_cluster_in_name = TRUE,
  also_write_rasters  = FALSE   # set TRUE if you also want the .tif for selected points
)


