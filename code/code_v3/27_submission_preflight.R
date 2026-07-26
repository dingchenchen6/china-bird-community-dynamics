#!/usr/bin/env Rscript
## Fail-fast offline preflight for the submission rerun.

CODE_V3 <- Sys.getenv("V3_CODE_DIR", file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
ensure_v3_dirs()

required_packages <- c("data.table", "sf", "terra", "spOccupancy", "readr",
  "dplyr", "tidyr", "abind", "coda", "posterior", "brms", "cmdstanr",
  "ape", "vegan", "betapart", "fundiversity", "Kendall", "ggplot2",
  "patchwork", "scales", "forcats", "stringr", "officer", "rvg")
pkg_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing offline R packages: ", paste(names(pkg_ok)[!pkg_ok], collapse = ", "))

input_paths <- c(
  events = file.path(PROJECT_ROOT, "data", "derived_v3", "combined_events_dedup_v3.rds"),
  grid = file.path(PROJECT_ROOT, "data", "derived_v2", "china_grid_100km_v2.rds"),
  environment = file.path(PROJECT_ROOT, "data", "derived_v2", "grid_environment_dynamic_occupancy.rds"))
if (any(!file.exists(input_paths))) stop("Missing local inputs: ",
  paste(names(input_paths)[!file.exists(input_paths)], collapse = ", "))

survey_path <- v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds")
if (file.exists(survey_path)) {
  survey <- readRDS(survey_path)
  d <- dim(survey$y)
  if (length(d) != 4L || d[3] != 5L || d[4] != PERIOD_LENGTH) stop("Survey dimensions invalid.")
  if (!identical(dim(survey$visited_mask), d[-1])) stop("Visited-mask dimensions invalid.")
  needed_det <- c("log_events", "log_duration", "has_duration", "source_ebird_prop")
  if (!all(needed_det %in% names(survey$det_covs))) stop("Detection covariates invalid.")
  if (any(vapply(survey$det_covs[needed_det], function(a) !identical(dim(a), d[-1]), logical(1)))) {
    stop("Detection-covariate dimensions do not match y.")
  }
  if (stats::sd(survey$det_covs$source_ebird_prop, na.rm = TRUE) == 0) {
    stop("Source detection covariate has no variation.")
  }
  message("[27] Canonical survey valid: ", paste(d, collapse = " x "))
} else {
  message("[27] Canonical survey not built yet; stage 02 is required.")
}

checks <- data.frame(check = c(paste0("package:", names(pkg_ok)), paste0("input:", names(input_paths))),
  passed = c(unname(pkg_ok), file.exists(input_paths)),
  path_or_version = c(vapply(names(pkg_ok), function(p) as.character(packageVersion(p)), ""), input_paths))
write.csv(checks, v3_file("results", "table_submission_preflight"), row.names = FALSE)
message("[27] Offline preflight passed; no network resources are required by the model run.")
