#!/usr/bin/env Rscript
## Canonical 500-species temporal sensitivity model (non-spatial, AR1).

suppressPackageStartupMessages({
  library(spOccupancy)
  library(readr)
  library(dplyr)
})
CODE_V3 <- Sys.getenv("V3_CODE_DIR", file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
ensure_v3_dirs()

run_label <- Sys.getenv("V3_RUN_LABEL", "submission_500sp_temporal")
max_species <- as.integer(Sys.getenv("V3_MAX_SPECIES", "500"))
chain_text <- Sys.getenv("V3_CHAIN_ID", "")
chain_id <- if (nzchar(chain_text)) as.integer(chain_text) else NA_integer_
n_omp <- as.integer(Sys.getenv("V3_N_OMP_THREADS", Sys.getenv("V3_N_OMP", "1")))
include_source <- Sys.getenv("V3_INCLUDE_SOURCE", "0") == "1"

survey_path <- v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds")
env_path <- v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds")
if (!file.exists(survey_path) || !file.exists(env_path)) {
  stop("Canonical survey/environment inputs are missing; run stages 02-03.")
}
survey <- readRDS(survey_path); grid_env <- readRDS(env_path)
if (length(dim(survey$y)) != 4L) stop("500-species model requires canonical 4D y.")
sp_idx <- seq_len(min(max_species, length(survey$species)))
y <- survey$y[sp_idx, , , , drop = FALSE]
n_sp <- dim(y)[1]; n_sites <- dim(y)[2]; n_periods <- dim(y)[3]
sites <- survey$sites

grid_env <- grid_env[match(sites, grid_env$grid_cell), , drop = FALSE]
if (nrow(grid_env) != n_sites || any(is.na(grid_env$grid_cell))) stop("Environment/site mismatch.")
standardize <- function(x) {
  x <- as.numeric(x); ok <- is.finite(x)
  if (!any(ok) || stats::sd(x[ok]) == 0) return(rep(0, length(x)))
  x[!ok] <- median(x[ok]); as.numeric(scale(x))
}
occ_vars <- intersect(c("bio4", "bio7", "bio11", "bio13", "elev_mean", "elev_sd",
  "texture_shannon", "habitat_diversity_shannon", "hfi_mean",
  "landcover_built", "landcover_cropland", "centroid_lon", "centroid_lat"), names(grid_env))
occ_covs <- setNames(lapply(occ_vars, function(v) matrix(
  rep(standardize(grid_env[[v]]), n_periods), nrow = n_sites, ncol = n_periods)), occ_vars)
occ_covs$year_scaled <- matrix(rep(as.numeric(scale(seq_len(n_periods))), each = n_sites),
                               nrow = n_sites, ncol = n_periods)

det_names <- c("log_events", "log_duration", "has_duration")
if (include_source) det_names <- c(det_names, "source_ebird_prop")
if (!all(det_names %in% names(survey$det_covs))) stop("Required detection covariates missing.")
det_covs <- survey$det_covs[det_names]
for (nm in intersect(c("log_events", "log_duration"), names(det_covs))) {
  ok <- is.finite(det_covs[[nm]])
  if (sum(ok) > 1 && sd(det_covs[[nm]][ok]) > 0) {
    det_covs[[nm]][ok] <- as.numeric(scale(det_covs[[nm]][ok]))
  }
}

site_has_data <- apply(y, 2, function(a) any(!is.na(a)))
if (!all(site_has_data)) {
  y <- y[, site_has_data, , , drop = FALSE]
  occ_covs <- lapply(occ_covs, function(a) a[site_has_data, , drop = FALSE])
  det_covs <- lapply(det_covs, function(a) a[site_has_data, , , drop = FALSE])
  n_sites <- sum(site_has_data)
}
data_list <- list(y = y, occ.covs = occ_covs, det.covs = det_covs)
occ_formula <- reformulate(c(occ_vars, "year_scaled"))
det_formula <- reformulate(det_names)

n_batch <- FULL_N_BATCH; n_burn <- FULL_N_BURN; n_thin <- FULL_N_THIN
n_chains <- if (is.na(chain_id)) FULL_N_CHAINS else 1L
n_factors <- min(as.integer(Sys.getenv("V3_N_FACTORS", "10")), max(1L, n_sp - 1L))
set.seed(if (is.na(chain_id)) 5024L else 5024L + chain_id * 1000L)
priors <- list(
  beta.comm.normal = list(mean = 0, var = 2.72),
  alpha.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
  sigma.sq.t.ig = list(a = 2, b = 0.5), rho.unif = list(a = -1, b = 1))
inits <- list(beta.comm = 0, alpha.comm = 0, tau.sq.beta = 1,
  tau.sq.alpha = 1, sigma.sq.t = 0.5, rho = 0)

message(sprintf("[04_t500] %s: %d species x %d sites x %d periods x %d years",
  run_label, n_sp, n_sites, n_periods, dim(y)[4]))
t0 <- Sys.time()
fit <- tMsPGOcc(occ.formula = occ_formula, det.formula = det_formula,
  data = data_list, inits = inits, priors = priors, n.factors = n_factors,
  n.batch = n_batch, batch.length = 25, accept.rate = 0.43,
  n.omp.threads = n_omp, verbose = TRUE, ar1 = TRUE,
  n.report = max(1, floor(n_batch / 10)), n.burn = n_burn,
  n.thin = n_thin, n.chains = n_chains)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
stem <- paste0("tMsPGOcc_fit_", run_label,
               if (!is.na(chain_id)) paste0("_chain", chain_id) else "")
fit_path <- v3_file("derived", stem, "rds")
saveRDS(fit, fit_path, compress = "xz")
if (!file.exists(fit_path) || file.info(fit_path)$size < 1000) stop("Fit checkpoint failed validation.")
write_csv(data.frame(run_label = run_label, chain_id = chain_id,
  n_species = n_sp, n_sites = n_sites, n_periods = n_periods,
  n_replicates = dim(y)[4], n_batch = n_batch, n_burn = n_burn,
  n_thin = n_thin, elapsed_mins = elapsed, fit_bytes = file.info(fit_path)$size),
  v3_file("results", paste0("table_model_summary_", run_label,
    if (!is.na(chain_id)) paste0("_chain", chain_id) else "")))
log_time("04_t500", sprintf("DONE: %s (%.1f min)", basename(fit_path), elapsed))
