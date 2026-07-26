#!/usr/bin/env Rscript
## Compare base and data-source-aware 200-species occupancy trajectories.
suppressPackageStartupMessages(library(readr))
CODE_V3 <- Sys.getenv("V3_CODE_DIR", file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R")); source(file.path(CODE_V3, "utils_paths.R"))
base_label <- Sys.getenv("V3_BASE_LABEL", "submission_20260722_200sp_spatial_base")
source_label <- Sys.getenv("V3_SOURCE_LABEL", "submission_20260722_200sp_spatial_source")
base <- readRDS(v3_file("derived", paste0("stMsPGOcc_fit_", base_label), "rds"))
src <- readRDS(v3_file("derived", paste0("stMsPGOcc_fit_", source_label), "rds"))
if (!identical(dim(base$psi.samples)[2:4], dim(src$psi.samples)[2:4])) stop("Model dimensions differ.")
trajectory <- function(x) apply(x$psi.samples, c(1, 4), mean, na.rm = TRUE)
b <- trajectory(base); s <- trajectory(src)
n <- min(nrow(b), nrow(s)); b <- b[seq_len(n), , drop = FALSE]; s <- s[seq_len(n), , drop = FALSE]
out <- do.call(rbind, lapply(seq_len(ncol(b)), function(t) {
  delta <- s[, t] - b[, t]
  data.frame(period = paste0("P", t), base_mean = mean(b[, t]), source_mean = mean(s[, t]),
    delta_mean = mean(delta), delta_q025 = quantile(delta, .025), delta_q975 = quantile(delta, .975),
    p_delta_positive = mean(delta > 0))
}))
write_csv(out, v3_file("results", "table_source_detection_trajectory_sensitivity"))

a <- as.matrix(src$alpha.comm.samples)
cn <- colnames(a); idx <- grep("source_ebird_prop", cn, fixed = TRUE)
if (length(idx)) {
  vals <- as.numeric(a[, idx[1]])
  write_csv(data.frame(term = cn[idx[1]], mean = mean(vals), q025 = quantile(vals, .025),
    q975 = quantile(vals, .975), p_positive = mean(vals > 0)),
    v3_file("results", "table_source_detection_coefficient"))
}
message("[29] Source-detection sensitivity comparison written.")
