#!/usr/bin/env Rscript
## Build the canonical breeding-season detection history used by every model.
## y dimensions: species x grid cell x 5-year period x year within period.

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(readr)
  library(dplyr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_spatial.R"))
ensure_v3_dirs()

log_time("02", "Starting canonical four-dimensional survey-history build")

events_path <- file.path(PROJECT_ROOT, "data", "derived_v3", "combined_events_dedup_v3.rds")
if (!file.exists(events_path)) {
  events_path <- file.path(DIRS$v2_derived, "combined_events_merged_dedup_2000_2025.rds")
}
if (!file.exists(events_path)) stop("Deduplicated events are missing. Run stage 01 first.")

events <- as.data.table(readRDS(events_path))
needed <- c("species", "longitude", "latitude", "year", "month", "source")
missing_cols <- setdiff(needed, names(events))
if (length(missing_cols)) stop("Event data lack: ", paste(missing_cols, collapse = ", "))
events <- events[
  year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI &
    month %in% BREEDING_MONTHS & is.finite(longitude) & is.finite(latitude) &
    !is.na(species) & nzchar(species)
]
if (!nrow(events)) stop("No breeding-season events remain after filtering.")

period_starts <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH)
n_periods <- length(period_starts)
events[, period_index := ((year - ANALYSIS_YR_LO) %/% PERIOD_LENGTH) + 1L]
events[, year_in_period := ((year - ANALYSIS_YR_LO) %% PERIOD_LENGTH) + 1L]

audit_all <- as.data.table(readRDS(events_path))[
  year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI,
  .(n_total = .N), by = .(period_index = ((year - ANALYSIS_YR_LO) %/% PERIOD_LENGTH) + 1L)
]
audit_breeding <- events[, .(n_breeding = .N), by = period_index]
audit <- merge(data.table(period_index = seq_len(n_periods)), audit_all,
               by = "period_index", all.x = TRUE)
audit <- merge(audit, audit_breeding, by = "period_index", all.x = TRUE)
audit[is.na(n_total), n_total := 0L]
audit[is.na(n_breeding), n_breeding := 0L]
audit[, `:=`(
  period = paste0("P", period_index),
  pct_retained = fifelse(n_total > 0, 100 * n_breeding / n_total, NA_real_)
)]
write_csv(as.data.frame(audit), v3_file("results", "table_breeding_filter_audit"))

grid_path <- file.path(PROJECT_ROOT, "data", "derived_v2", "china_grid_100km_v2.rds")
tagged_grid <- v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds")
grid_sf <- NULL
if (GRID_SIZE_KM == 100L && file.exists(grid_path)) grid_sf <- readRDS(grid_path)
if (is.null(grid_sf) && file.exists(tagged_grid)) grid_sf <- readRDS(tagged_grid)
if (is.null(grid_sf)) {
  pts <- st_as_sf(events, coords = c("longitude", "latitude"), crs = 4326,
                  remove = FALSE)
  pts_aea <- project_china_albers(pts)
  grid_sf <- st_sf(geometry = st_make_grid(pts_aea,
    cellsize = GRID_SIZE_KM * 1000, square = TRUE))
  grid_sf$grid_cell <- seq_len(nrow(grid_sf))
  grid_sf <- st_transform(grid_sf, 4326)
  saveRDS(grid_sf, tagged_grid)
}
if (!"grid_cell" %in% names(grid_sf)) grid_sf$grid_cell <- seq_len(nrow(grid_sf))

pts <- st_as_sf(events, coords = c("longitude", "latitude"), crs = 4326,
                remove = FALSE)
sf_use_s2(FALSE)
matched <- st_join(pts, grid_sf[, "grid_cell"], join = st_intersects, left = FALSE)
dt <- as.data.table(st_drop_geometry(matched))
rm(events, pts, matched); gc()

species_counts <- dt[, .(n_records = .N,
  n_grid_years = uniqueN(paste(grid_cell, year, sep = "_"))), by = species]
setorder(species_counts, -n_records, species)
species_counts <- species_counts[n_records >= MIN_VISITS]
candidate_species <- species_counts$species
sites <- sort(unique(dt$grid_cell))
n_sp <- length(candidate_species); n_sites <- length(sites)
if (!n_sp || !n_sites) stop("No eligible species or grid cells after matching.")
write_csv(as.data.frame(species_counts),
          v3_file("results", paste0("table_candidate_species", GRID_TAG, "_v3")))

dt <- dt[species %chin% candidate_species]
dt[, site_index := match(grid_cell, sites)]
species_lookup <- data.table(species = candidate_species,
                             species_index = seq_along(candidate_species))
dt <- species_lookup[dt, on = "species"]

dims_y <- c(n_sp, n_sites, n_periods, PERIOD_LENGTH)
dims_det <- c(n_sites, n_periods, PERIOD_LENGTH)
y <- array(NA_integer_, dim = dims_y,
           dimnames = list(candidate_species, as.character(sites), paste0("P", seq_len(n_periods)),
                           as.character(seq_len(PERIOD_LENGTH))))
visited <- unique(dt[, .(site_index, period_index, year_in_period)])
visited_mask <- array(FALSE, dim = dims_det)
visited_mask[as.matrix(visited)] <- TRUE

# Fill surveyed grid-years with non-detections, then overwrite observed species with 1.
for (k in seq_len(nrow(visited))) {
  y[, visited$site_index[k], visited$period_index[k], visited$year_in_period[k]] <- 0L
}
detections <- unique(dt[, .(species_index, site_index, period_index, year_in_period)])
y[as.matrix(detections)] <- 1L

effort <- dt[, .(
  n_events = .N,
  mean_duration_min = if (any(is.finite(duration_min) & duration_min > 0))
    mean(duration_min[is.finite(duration_min) & duration_min > 0]) else NA_real_,
  source_ebird_prop = mean(tolower(source) == "ebird", na.rm = TRUE),
  n_sources = uniqueN(source)
), by = .(site_index, grid_cell, period_index, year_in_period, year)]

make_det_array <- function(value, default = NA_real_) {
  out <- array(default, dim = dims_det)
  out[as.matrix(effort[, .(site_index, period_index, year_in_period)])] <- value
  out
}
has_duration <- as.numeric(is.finite(effort$mean_duration_min) & effort$mean_duration_min > 0)
det_covs <- list(
  log_events = make_det_array(log1p(effort$n_events)),
  log_duration = make_det_array(ifelse(has_duration == 1,
    log1p(effort$mean_duration_min), 0)),
  has_duration = make_det_array(has_duration),
  source_ebird_prop = make_det_array(effort$source_ebird_prop),
  source_mixed = make_det_array(as.numeric(effort$n_sources > 1))
)

period_labels <- sprintf("%d-%d", period_starts,
  pmin(period_starts + PERIOD_LENGTH - 1L, ANALYSIS_YR_HI))
survey <- list(
  y = y,
  det_covs = det_covs,
  visited_mask = visited_mask,
  sites = sites,
  species = candidate_species,
  periods = period_starts,
  period_labels = period_labels,
  replicate_years = seq_len(PERIOD_LENGTH),
  breeding_months = BREEDING_MONTHS,
  grid_size_km = GRID_SIZE_KM,
  effort_summary = effort,
  input_file = normalizePath(events_path),
  built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
survey_path <- v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds")
saveRDS(survey, survey_path, compress = "xz")
saveRDS(det_covs,
        v3_file("derived", paste0("detection_covariates", GRID_TAG, "_v3"), "rds"),
        compress = "xz")

summary_out <- data.frame(
  n_species = n_sp, n_sites = n_sites, n_periods = n_periods,
  n_replicate_years = PERIOD_LENGTH, n_events_breeding = nrow(dt),
  n_visited_grid_years = nrow(visited),
  observed_fraction = mean(visited_mask),
  source_ebird_fraction = mean(tolower(dt$source) == "ebird", na.rm = TRUE),
  grid_size_km = GRID_SIZE_KM,
  breeding_months = paste(BREEDING_MONTHS, collapse = ",")
)
write_csv(summary_out,
  v3_file("results", paste0("table_02_survey_summary", GRID_TAG, "_v3")))
log_time("02", sprintf("DONE: %d species x %d sites x %d periods x %d years",
                        n_sp, n_sites, n_periods, PERIOD_LENGTH))
