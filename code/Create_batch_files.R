library("writexl")

# --- edit this to your folder ---
folder <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/Kielder"

# Regex pattern:
# inds_ + 4 digits + _cluster + number + _rep + number
pattern <- "^inds_[0-9]{4}_cluster[0-9]+_rep[0-9]+" # need to locat the inds files generated in hab model (define init pop)

# List all files in folder (recursive in case of subfolders)
files <- list.files(folder, full.names = TRUE, recursive = TRUE)

# Keep only those matching the pattern
matched <- files[str_detect(basename(files), pattern)]

# Create a dataframe for export
df <- data.frame(file_path = matched)

# Export to Excel
write_xlsx(df, "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/Kielder/inds_file_paths.xlsx")

df # this will be a table with all ind files available and their metadata






# Packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(fs)
library/writexl)

# ---------------------------------------------------------------------------
# build_batch_table()  -- updated per your request
# ---------------------------------------------------------------------------
# Arguments:
# - batches_path: path to your batch template (TSV/CSV)
# - inds_dir:     root folder to scan for inds_* files (recursive)
# - start_year:   minimum year to include (e.g., 2020)
# - site:         optional site filter (substring in path, case-insensitive)
# - dynamic_path: dynamic landscapes path; either literal or template using {year}
# - init_rel_dir: relative directory where init_inds_file should point (default "../data/processed")
# - phenips_template: optional template for phenips_file using {year}
# - year_max:     max year to include (default Inf). Set = start_year to get only that year.
# - extensions:   allowed inds file extensions (default c("txt","txtZ"))
# - recursive:    recurse in inds_dir (default TRUE)
# - out_xlsx/csv: optional output files to write
# - one_row_per_batch: TRUE for batches × inds cartesian; FALSE = use only first batch row as template
# - make_batches_consecutive: TRUE to renumber Batch1..BatchN in the final output
# - consecutive_strategy: "global" (default) or "per_year" (reset numbering each year)
# - consecutive_sort: "none" (keep join order) or "stable" (sort by year/cluster/rep then original batch)
# ---------------------------------------------------------------------------
build_batch_table <- function(
    batches_path,
    inds_dir,
    start_year,
    site,
    dynamic_path,
    init_rel_dir = "../data/processed",
    phenips_template = NULL,
    year_max = Inf,
    extensions = c("txt","txtZ"),
    recursive = TRUE,
    out_xlsx = NULL,
    out_csv  = NULL,
    one_row_per_batch = TRUE,
    make_batches_consecutive = TRUE,
    consecutive_strategy = c("global", "per_year"),
    consecutive_sort = c("none", "stable")
) {
  consecutive_strategy <- match.arg(consecutive_strategy)
  consecutive_sort     <- match.arg(consecutive_sort)
  
  # --- Local helpers
  normalize_names <- function(df) {
    names(df) <- gsub("\\s+", "_", names(df))
    df
  }
  read_batches <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("tsv", "tab", "txt")) {
      readr::read_delim(path, delim = "\t", col_types = readr::cols())
    } else if (ext == "csv") {
      readr::read_csv(path, show_col_types = FALSE)
    } else {
      first <- readLines(path, n = 1, warn = FALSE)
      if (grepl("\t", first)) {
        readr::read_delim(path, delim = "\t", col_types = readr::cols())
      } else {
        readr::read_csv(path, show_col_types = FALSE)
      }
    }
  }
  fill_year <- function(x, yr) {
    if (stringr::str_detect(x, stringr::fixed("{year}"))) {
      stringr::str_replace_all(x, stringr::fixed("{year}"), as.character(yr))
    } else if (stringr::str_detect(x, "\\d{4}")) {
      stringr::str_replace(x, "\\d{4}", as.character(yr))
    } else {
      x
    }
  }
  
  # --- Read template and normalize names
  batches <- read_batches(batches_path)
  batches <- normalize_names(batches)
  
  # Loosen/adjust if your template differs
  required_cols <- c("batch_name", "phenips_file", "init_inds_file",
                     "dynamic_landscapes", "years_for_destruction",
                     "apply_destruction_annually", "keep_input_files",
                     "pop_impacts_quality", "replicates", "years")
  missing_cols <- setdiff(required_cols, names(batches))
  if (length(missing_cols)) {
    stop("Missing required columns in batches file: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # --- Find inds files
  ext_rx <- paste0("\\.(", paste(extensions, collapse = "|"), ")$")
  inds_rx <- paste0("(?i)inds_\\d{4}_cluster\\d+_rep\\d+.*", ext_rx)
  
  inds_files <- fs::dir_ls(
    path = inds_dir,
    recurse = recursive,
    type = "file",
    regexp = inds_rx
  )
  if (!length(inds_files)) stop("No inds_* files found in: ", inds_dir)
  
  # Optional site filter: keep only paths containing 'site'
  inds_files2 <- inds_files
  if (!is.null(site) && nzchar(site)) {
    inds_files2 <- inds_files[stringr::str_detect(tolower(inds_files), tolower(site))]
  }
  if (!length(inds_files2)) {
    stop("After filtering by site='", site, "', no inds files remain.")
  }
  
  # --- Parse parts and filter by year
  inds_tbl <- tibble::tibble(
    full_path = as.character(inds_files2),
    file_name = basename(as.character(inds_files2))
  ) %>%
    mutate(
      year    = as.integer(stringr::str_match(file_name, "^inds_(\\d{4})_")[, 2]),
      cluster = as.integer(stringr::str_match(file_name, "_cluster(\\d+)")[, 2]),
      rep     = as.integer(stringr::str_match(file_name, "_rep(\\d+)")[, 2])
    ) %>%
    filter(!is.na(year), year >= start_year, year <= year_max) %>%
    # Use *_new names to avoid collisions in the join
    mutate(
      init_inds_file_new     = file.path(init_rel_dir, file_name, fsep = "/"),
      dynamic_landscapes_new = fill_year(dynamic_path, year),
      phenips_new            = if (!is.null(phenips_template))
        fill_year(phenips_template, year)
      else NA_character_
    )
  
  if (nrow(inds_tbl) == 0) {
    stop("No inds files left after filtering by year range [",
         start_year, ", ", year_max, "].")
  }
  
  # --- Combine with batches
  base_tbl  <- if (isTRUE(one_row_per_batch)) batches else batches[1, ]
  combined  <- tidyr::crossing(base_tbl, inds_tbl, .name_repair = "unique")
  
  # --- Overwrite template columns with inds-derived values
  expanded <- combined %>%
    mutate(
      init_inds_file     = init_inds_file_new,
      dynamic_landscapes = dynamic_landscapes_new,
      phenips_file       = if (!is.null(phenips_template)) phenips_new else phenips_file,
      site               = site,
      .orig_batch_idx    = match(batch_name, unique(batch_name))
    ) %>%
    select(
      batch_name, phenips_file, init_inds_file, dynamic_landscapes,
      years_for_destruction, apply_destruction_annually, keep_input_files,
      pop_impacts_quality, replicates, years,
      site, year, cluster, rep, .orig_batch_idx
    )
  
  # --- Make batch_name consecutive (optional)
  if (isTRUE(make_batches_consecutive)) {
    
    # Optional deterministic ordering
    if (consecutive_sort == "stable") {
      expanded <- expanded %>%
        arrange(year, cluster, rep, .orig_batch_idx)
    }
    
    if (consecutive_strategy == "global") {
      expanded <- expanded %>%
        mutate(batch_name = paste0("Batch", dplyr::row_number()))
    } else if (consecutive_strategy == "per_year") {
      expanded <- expanded %>%
        group_by(year) %>%
        mutate(batch_name = paste0("Batch", dplyr::row_number())) %>%
        ungroup()
    }
  }
  
  # Drop helper
  expanded <- expanded %>% select(-.orig_batch_idx)
  
  # --- Write outputs (optional)
  if (!is.null(out_xlsx)) writexl::write_xlsx(expanded, out_xlsx)
  if (!is.null(out_csv))  readr::write_csv(expanded, out_csv)
  
  expanded
}


expanded <- build_batch_table(
  batches_path  = "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/Batches_FoD_2020.csv",
  inds_dir      = "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed/FoD",
  start_year    = 2020,
  year_max      = 2020,  # only that year
  site          = "FoD",
  dynamic_path  = "../data/dynamic_landscapes_{year}.csv",
  init_rel_dir  = "../data/processed",
  # phenips_template = "../data/phenipsFoD{year}.csv",   # optional
  out_xlsx      = "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/Batches_FoD_2020_expanded.xlsx",
  out_csv       = "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data/Batches_FoD_2020_expanded.csv",
  one_row_per_batch = TRUE,
  make_batches_consecutive = TRUE,   # renumber Batch1..BatchN
  consecutive_strategy     = "global",
  consecutive_sort         = "stable"
)





# -----------------------------------------------
# Run build_batch_table() for 3 sites × 2 years
# -----------------------------------------------

sites <- c("FoD", "NYM", "Kielder")
years <- c(2020, 2030)

# Base paths
base_data <- "Z:/Common/Fiona Plenderleith/IpsShift2/Ips_RangeShiftR/data"
base_inds <- "Z:/Common/Fiona Plenderleith/IpsShift2/Landscapes/processed"

# Your existing batch template (works for all sites)
batches_template <- file.path(base_data, "Batches_FoD_2020.csv")

for (site in sites) {
  for (yr in years) {
    
    message("\n----------------------------------------")
    message("Running: site = ", site, " | year = ", yr)
    message("----------------------------------------")
    
    # Input folder for inds_* files
    inds_dir <- file.path(base_inds, site)
    
    # Output file names
    out_xlsx <- file.path(base_data, sprintf("Batches_%s_%d_expanded.xlsx", site, yr))
    out_csv  <- file.path(base_data, sprintf("Batches_%s_%d_expanded.csv",   site, yr))
    
    # Phenips file pattern for the site
    phenips_tmpl <- sprintf("../data/phenips%s{year}.csv", site)
    
    # Run your function
    res <- tryCatch(
      build_batch_table(
        batches_path  = batches_template,
        inds_dir      = inds_dir,
        start_year    = yr,
        year_max      = yr,  # only that year
        site          = site,
        dynamic_path  = "../data/dynamic_landscapes_{year}.csv",
        init_rel_dir  = "../data/processed",
        phenips_template = phenips_tmpl,
        out_xlsx      = out_xlsx,
        out_csv       = out_csv,
        one_row_per_batch = TRUE,
        make_batches_consecutive = TRUE,
        consecutive_strategy     = "global",
        consecutive_sort         = "stable"
      ),
      error = function(e) {
        message("❌ Failed for site=", site, ", year=", yr)
        message("   Error: ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(res)) {
      message("✔ SUCCESS: ", nrow(res), " rows written")
      message("   Saved: ", out_csv)
      message("   Saved: ", out_xlsx)
    }
  }
}

message("\nDone — all 6 batch files attempted.")
