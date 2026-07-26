#!/usr/bin/env Rscript
## Memory-aware chain validation, convergence diagnostics and posterior merge.

suppressPackageStartupMessages({
  library(abind)
  library(readr)
  library(posterior)
})
CODE_V3 <- Sys.getenv("V3_CODE_DIR", file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
ensure_v3_dirs()

run_label <- RUN_LABEL
model_prefix <- Sys.getenv("V3_MODEL_PREFIX", "stMsPGOcc")
n_expected <- as.integer(Sys.getenv("V3_N_CHAINS", "4"))
psi_max <- as.integer(Sys.getenv("V3_PSI_MERGE_DRAWS", as.character(PSI_MAX_DRAWS)))
paths <- file.path(DIRS$derived,
  paste0(model_prefix, "_fit_", run_label, "_chain", seq_len(n_expected), ".rds"))
missing <- paths[!file.exists(paths)]
if (length(missing)) stop("Missing exact chain files: ", paste(basename(missing), collapse = ", "))

fits <- vector("list", n_expected)
chain_rows <- vector("list", n_expected)
for (ch in seq_len(n_expected)) {
  message("[26] Reading chain ", ch, ": ", basename(paths[ch]))
  fits[[ch]] <- tryCatch(readRDS(paths[ch]), error = function(e) {
    stop("Unreadable chain ", ch, ": ", conditionMessage(e))
  })
  d <- dim(fits[[ch]]$psi.samples)
  if (length(d) != 4L || any(d < 1)) stop("Invalid psi.samples in chain ", ch)
  chain_rows[[ch]] <- data.frame(chain = ch, file = basename(paths[ch]),
    bytes = file.info(paths[ch])$size, n_draws = d[1], n_species = d[2],
    n_sites = d[3], n_periods = d[4], readable = TRUE)
}
shape_key <- vapply(fits, function(f) paste(dim(f$psi.samples)[-1], collapse = "x"), "")
if (length(unique(shape_key)) != 1L) stop("Chain posterior dimensions disagree.")
write_csv(do.call(rbind, chain_rows),
          v3_file("results", paste0("table_chain_integrity_", run_label)))

diag_slots <- intersect(c("beta.comm.samples", "alpha.comm.samples",
  "tau.sq.beta.samples", "tau.sq.alpha.samples", "phi.samples",
  "sigma.sq.samples", "rho.samples"), names(fits[[1]]))
diag_out <- list()
for (slot in diag_slots) {
  mats <- lapply(fits, function(f) as.matrix(f[[slot]]))
  if (any(vapply(mats, is.null, logical(1)))) next
  n_iter <- min(vapply(mats, nrow, integer(1)))
  n_param <- min(vapply(mats, ncol, integer(1)))
  for (j in seq_len(n_param)) {
    x <- sapply(mats, function(m) m[seq_len(n_iter), j])
    diag_out[[length(diag_out) + 1L]] <- data.frame(
      group = sub("\\.samples$", "", slot),
      parameter = colnames(mats[[1]])[j] %||% paste0("V", j),
      rhat = posterior::rhat(x),
      ess_bulk = posterior::ess_bulk(x),
      ess_tail = posterior::ess_tail(x))
  }
}
diag_df <- if (length(diag_out)) do.call(rbind, diag_out) else data.frame()
write_csv(diag_df, v3_file("results", paste0("table_convergence_diagnostics_", run_label)))
if (!nrow(diag_df) || any(!is.finite(diag_df$rhat))) stop("Convergence diagnostics unavailable.")

merge_first_dimension <- function(xs) {
  if (any(vapply(xs, is.null, logical(1)))) return(NULL)
  ds <- lapply(xs, dim)
  if (any(vapply(ds, is.null, logical(1)))) return(unlist(xs, use.names = FALSE))
  trailing <- vapply(ds, function(d) paste(d[-1], collapse = "x"), "")
  if (length(unique(trailing)) != 1L) return(xs[[1]])
  do.call(abind::abind, c(xs, along = 1))
}

combined <- fits[[1]]
# Retain an even, chain-balanced psi subset; omit latent arrays not used downstream.
per_chain <- max(1L, floor(psi_max / n_expected))
psi_parts <- lapply(fits, function(f) {
  idx <- unique(round(seq(1, dim(f$psi.samples)[1], length.out = per_chain)))
  f$psi.samples[idx, , , , drop = FALSE]
})
combined$psi.samples <- do.call(abind::abind, c(psi_parts, along = 1))
combined$z.samples <- NULL
combined$like.samples <- NULL

sample_slots <- setdiff(names(combined)[grepl("\\.samples$", names(combined))],
                        c("psi.samples", "z.samples", "like.samples"))
for (slot in sample_slots) {
  merged <- merge_first_dimension(lapply(fits, `[[`, slot))
  if (!is.null(merged)) combined[[slot]] <- merged
}
combined$n.chains <- n_expected
combined$chain_ids <- seq_len(n_expected)
combined$merge_metadata <- list(run_label = run_label, model_prefix = model_prefix,
  psi_draws = dim(combined$psi.samples)[1], merged_at = Sys.time(),
  chain_files = basename(paths))
out <- v3_file("derived", paste0(model_prefix, "_fit_", run_label), "rds")
saveRDS(combined, out, compress = "gzip")
if (!file.exists(out) || file.info(out)$size < 1000) stop("Combined fit save failed.")

gate <- data.frame(run_label = run_label, model_prefix = model_prefix,
  n_chains = n_expected, max_rhat = max(diag_df$rhat, na.rm = TRUE),
  min_ess_bulk = min(diag_df$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(diag_df$ess_tail, na.rm = TRUE),
  rhat_pass = all(diag_df$rhat <= RHAT_THRESHOLD),
  ess_pass = all(diag_df$ess_bulk >= ESS_THRESHOLD),
  combined_file = basename(out), combined_bytes = file.info(out)$size)
write_csv(gate, v3_file("results", paste0("table_convergence_gate_", run_label)))
message(sprintf("[26] Combined %d chains; max Rhat %.3f, min bulk ESS %.0f",
                n_expected, gate$max_rhat, gate$min_ess_bulk))
if (!gate$rhat_pass || !gate$ess_pass) {
  stop("Convergence gate failed. Do not continue to submission post-processing.")
}
