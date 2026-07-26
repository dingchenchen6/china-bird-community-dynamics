#!/usr/bin/env Rscript
## Final machine-readable gate for the server rerun.
suppressPackageStartupMessages(library(readr))
CODE_V3 <- Sys.getenv("V3_CODE_DIR", file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R")); source(file.path(CODE_V3, "utils_paths.R"))
label <- Sys.getenv("V3_RUN_LABEL", "submission_20260722_500sp_temporal")
ext <- paste0(label, "_extended")
required <- c(
  combined_fit = v3_file("derived", paste0("tMsPGOcc_fit_", label), "rds"),
  convergence = v3_file("results", paste0("table_convergence_gate_", label)),
  metrics = v3_file("results", paste0("table_community_metrics_with_cri_", ext)),
  trends = v3_file("results", paste0("table_trend_summary_", ext)),
  pdecline = v3_file("results", paste0("table_functional_trend_pdecline_", label)),
  homogenization = v3_file("results", paste0("table_homogenization_trend_", label)))
present <- file.exists(required)
if (all(present)) {
  metrics <- read_csv(required[["metrics"]], show_col_types = FALSE)
  n_metrics <- length(unique(metrics$metric))
  finite_metrics <- sum(vapply(split(metrics$mean, metrics$metric), function(x) any(is.finite(x)), logical(1)))
} else {
  n_metrics <- finite_metrics <- 0L
}
gate <- data.frame(item = names(required), path = unname(required), present = present)
gate <- rbind(gate, data.frame(item = c("metric_count_at_least_20", "all_metrics_finite"),
  path = c(as.character(n_metrics), as.character(finite_metrics)),
  present = c(n_metrics >= 20L, n_metrics >= 20L && finite_metrics == n_metrics)))
write_csv(gate, v3_file("results", "table_submission_rerun_postflight"))
if (!all(gate$present)) stop("Submission postflight failed: ", paste(gate$item[!gate$present], collapse = ", "))
message("[28] Postflight passed: ", n_metrics, " metrics with finite results.")
