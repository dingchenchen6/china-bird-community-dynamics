#!/usr/bin/env Rscript
## Canonical 200-species spatial multi-species occupancy model.
## Both the base and source-aware detection models consume the same 4D survey.

gc()
suppressPackageStartupMessages({
  library(spOccupancy)
  library(readr)
  library(dplyr)
  library(sf)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
ensure_v3_dirs()

include_source <- Sys.getenv("V3_INCLUDE_SOURCE", "0") == "1"
max_species <- as.integer(Sys.getenv("V3_MAX_SPECIES", "200"))
is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
chain_text <- Sys.getenv("V3_CHAIN_ID", "")
chain_id <- if (nzchar(chain_text)) as.integer(chain_text) else NA_integer_
n_omp <- as.integer(Sys.getenv("V3_N_OMP_THREADS", "1"))
log_time("04", sprintf("Starting %s (source covariate: %s, chain: %s)",
  run_label, include_source, ifelse(is.na(chain_id), "all", chain_id)))

survey_path <- v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds")
grid_env_path <- v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds")
if (!file.exists(survey_path)) stop("Canonical survey history missing: run stage 02.")
if (!file.exists(grid_env_path)) stop("Grid environment missing: run stage 03.")
survey <- readRDS(survey_path)
grid_env <- readRDS(grid_env_path)
if (length(dim(survey$y)) != 4L) {
  stop("Submission models require y[species, site, period, year-within-period].")
}
if (is.null(survey$det_covs) || is.null(survey$visited_mask)) {
  stop("Canonical detection covariates/visited mask are missing; rebuild stage 02.")
}

candidate_species <- head(survey$species, max_species)
sp_idx <- match(candidate_species, survey$species)
y <- survey$y[sp_idx, , , , drop = FALSE]
sites <- survey$sites
n_sp <- dim(y)[1]; n_sites <- dim(y)[2]; n_periods <- dim(y)[3]
if (n_sp < max_species) warning("Only ", n_sp, " eligible species are available.")
if (!any(is.na(y)) || !any(y == 0, na.rm = TRUE) || !any(y == 1, na.rm = TRUE)) {
  stop("Invalid y structure: expected NA unvisited cells plus observed 0/1 values.")
}

grid_env <- grid_env[match(sites, grid_env$grid_cell), , drop = FALSE]
if (nrow(grid_env) != n_sites || any(is.na(grid_env$grid_cell))) {
  stop("Grid environment cannot be aligned one-to-one with survey sites.")
}
standardize <- function(x) {
  x <- as.numeric(x); ok <- is.finite(x)
  if (!any(ok) || stats::sd(x[ok]) == 0) return(rep(0, length(x)))
  x[!ok] <- stats::median(x[ok]); as.numeric(scale(x))
}
occ_vars <- intersect(c("bio4", "bio7", "bio11", "bio13", "elev_mean", "elev_sd",
  "texture_shannon", "habitat_diversity_shannon", "hfi_mean",
  "landcover_built", "landcover_cropland", "centroid_lon", "centroid_lat"),
  names(grid_env))
if (!length(occ_vars)) stop("No occupancy covariates found in grid environment.")
occ_covs <- setNames(lapply(occ_vars, function(v) standardize(grid_env[[v]])), occ_vars)
occ_covs$year_scaled <- matrix(rep(as.numeric(scale(seq_len(n_periods))), each = n_sites),
                               nrow = n_sites, ncol = n_periods)

scale_det <- function(a) {
  ok <- is.finite(a)
  if (!any(ok) || stats::sd(a[ok]) == 0) return(a)
  a[ok] <- as.numeric(scale(a[ok])); a
}
det_names <- c("log_events", "log_duration", "has_duration")
if (include_source) det_names <- c(det_names, "source_ebird_prop")
missing_det <- setdiff(det_names, names(survey$det_covs))
if (length(missing_det)) stop("Missing detection covariates: ", paste(missing_det, collapse = ", "))
det_covs <- survey$det_covs[det_names]
det_covs[c("log_events", "log_duration")] <-
  lapply(det_covs[c("log_events", "log_duration")], scale_det)

coords_sf <- st_as_sf(grid_env, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
coords <- st_coordinates(st_transform(coords_sf,
  "+proj=aea +lat_1=25 +lat_2=47 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"))[, 1:2, drop = FALSE]
if (any(!is.finite(coords))) stop("Projected site coordinates contain non-finite values.")

data_list <- list(y = y, occ.covs = occ_covs, det.covs = det_covs,
                  coords = coords, grid.index = seq_len(n_sites))
occ_formula <- reformulate(c(occ_vars, "year_scaled"))
det_formula <- reformulate(det_names)

n_batch <- if (is_pilot) as.integer(Sys.getenv("V3_PILOT_N_BATCH", "10")) else FULL_N_BATCH
n_burn <- if (is_pilot) as.integer(Sys.getenv("V3_PILOT_N_BURN", "100")) else FULL_N_BURN
n_thin <- if (is_pilot) as.integer(Sys.getenv("V3_PILOT_N_THIN", "1")) else FULL_N_THIN
n_chains <- if (is.na(chain_id)) FULL_N_CHAINS else 1L
if (is_pilot) {
  pilot_n <- min(as.integer(Sys.getenv("V3_PILOT_SPECIES", "10")), n_sp)
  data_list$y <- data_list$y[seq_len(pilot_n), , , , drop = FALSE]
  n_sp <- pilot_n
}
n_factors <- min(as.integer(Sys.getenv("V3_N_FACTORS", "10")), max(1L, n_sp - 1L))

z_init <- apply(data_list$y, c(1, 2, 3), function(v) as.integer(any(v == 1, na.rm = TRUE)))
set.seed(if (is.na(chain_id)) 2024L else 2024L + 1000L * chain_id)
dmat <- as.matrix(dist(coords)); pairwise <- dmat[upper.tri(dmat)]
nearest <- apply(replace(dmat, row(dmat) == col(dmat), Inf), 1, min)
min_range <- max(GRID_SIZE_KM, median(nearest[is.finite(nearest)]))
max_range <- max(pairwise[is.finite(pairwise)])
phi_bounds <- sort(c(3 / max_range, 3 / min_range))
phi_init <- min(max(3 / median(pairwise), phi_bounds[1] * 1.01), phi_bounds[2] * 0.99)
rm(dmat, pairwise, nearest); gc()

inits <- list(alpha.comm = 0, beta.comm = 0, beta = 0, alpha = 0,
  tau.sq.beta = 1, tau.sq.alpha = 1, z = z_init, sigma.sq = 1, phi = phi_init)
priors <- list(
  beta.comm.normal = list(mean = 0, var = 2.72),
  alpha.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
  sigma.sq.ig = list(a = 2, b = 1),
  phi.unif = list(a = phi_bounds[1], b = phi_bounds[2]))
tuning <- list(phi = 0.5, sigma.sq = 1, rho = 0.2)

message(sprintf("[04] %s: %d species x %d sites x %d periods x %d years",
  run_label, dim(data_list$y)[1], n_sites, n_periods, dim(data_list$y)[4]))
message(sprintf("[04] MCMC %d x 25; burn=%d; thin=%d; chains=%d; OMP=%d",
  n_batch, n_burn, n_thin, n_chains, n_omp))
t0 <- Sys.time()
fit <- stMsPGOcc(occ.formula = occ_formula, det.formula = det_formula,
  data = data_list, inits = inits, priors = priors, tuning = tuning,
  cov.model = COV_MODEL, NNGP = TRUE, n.neighbors = N_NEIGHBORS,
  n.factors = n_factors, n.batch = n_batch, batch.length = 25,
  accept.rate = 0.43, n.omp.threads = n_omp, verbose = TRUE, ar1 = TRUE,
  n.report = max(1, floor(n_batch / 10)), n.burn = n_burn,
  n.thin = n_thin, n.chains = n_chains)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

stem <- paste0("stMsPGOcc_fit_", run_label,
               if (!is.na(chain_id)) paste0("_chain", chain_id) else "")
fit_path <- v3_file("derived", stem, "rds")
saveRDS(fit, fit_path, compress = "xz")
if (!file.exists(fit_path) || file.info(fit_path)$size < 1000) stop("Fit checkpoint failed validation.")
write_csv(data.frame(run_label = run_label, chain_id = chain_id,
  include_source = include_source, n_species = dim(data_list$y)[1], n_sites = n_sites,
  n_periods = n_periods, n_replicates = dim(data_list$y)[4], n_batch = n_batch,
  n_burn = n_burn, n_thin = n_thin, elapsed_mins = elapsed,
  fit_bytes = file.info(fit_path)$size),
  v3_file("results", paste0("table_model_summary_", run_label,
    if (!is.na(chain_id)) paste0("_chain", chain_id) else "")))
log_time("04", sprintf("DONE: %s (%.1f min)", basename(fit_path), elapsed))
