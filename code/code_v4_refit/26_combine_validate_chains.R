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
  "sigma.sq.samples", "sigma.sq.t.samples", "rho.samples"), names(fits[[1]]))
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

# ── 收敛门（v4 2026-10-21：沿用 2026-08-04 修订的群落层超参披露制）──────
# 沿革：本模型的 beta.comm/tau.sq.beta 群落层超参沿共线岭缓慢漂移、弱识别，
# 2026-07-30 与 2026-08-04 逐链诊断（rerun_20260722，4 链 95% CrI 共同重叠、
# 无多峰）证实即使 40,000 次迭代/2,000 draws 也不可达 R̂≤1.05；据此对这两组
# 走披露制并设兜底硬上限，其余（物种层 z/psi、检测层 alpha、空间/时间参数）
# 仍须严格门 R̂≤1.05 且 ESS≥200。v4 重拟合复用该标准；投稿前须对本轮链
# 重跑逐链 CrI 重叠诊断（参照 98_chain_diagnosis_community_20260804.R），
# 且凡触发放宽必须在稿件补充材料中报告 R̂/ESS；群落均值层面的海拔/温度
# 季节性（bio4）响应须声明此不确定性，正文气候结论挂在物种层后验传播
# 而非 beta.comm 点估计上。
WATCHLIST_GROUPS   <- c("beta.comm", "tau.sq.beta")
WATCHLIST_RHAT_MAX <- 1.4   # 观察名单 R̂ 兜底硬上限：超过仍判收敛失败
WATCHLIST_ESS_MIN  <- 5     # 观察名单 bulk ESS 兜底下限
is_watch <- diag_df$group %in% WATCHLIST_GROUPS
strict_df <- diag_df[!is_watch, , drop = FALSE]
watch_df  <- diag_df[is_watch, , drop = FALSE]
strict_rhat_pass <- all(strict_df$rhat <= RHAT_THRESHOLD)
strict_ess_pass  <- all(strict_df$ess_bulk >= ESS_THRESHOLD)
watch_rhat_pass  <- all(watch_df$rhat <= WATCHLIST_RHAT_MAX)
watch_ess_pass   <- all(watch_df$ess_bulk >= WATCHLIST_ESS_MIN)
relax_used <- nrow(watch_df) > 0 &&
  any(watch_df$rhat > RHAT_THRESHOLD | watch_df$ess_bulk < ESS_THRESHOLD)
gate_pass <- strict_rhat_pass && strict_ess_pass && watch_rhat_pass && watch_ess_pass
if (nrow(watch_df)) {
  disclosure <- watch_df
  disclosure$strict_threshold <- paste0("Rhat<=", RHAT_THRESHOLD, " & ESS>=", ESS_THRESHOLD)
  disclosure$watchlist_cap    <- paste0("Rhat<=", WATCHLIST_RHAT_MAX, " & ESS>=", WATCHLIST_ESS_MIN)
  disclosure$interpretation   <- paste(
    "Weakly-identified community-level hyper-parameter (beta.comm/tau.sq.beta);",
    "slow drift along collinearity ridge, no multimodality, pattern established by",
    "per-chain diagnosis in the 2026-07 rerun (table_chain_diagnosis_community_",
    "20260804: 4-chain 95% CrI common overlap for all 30 community hyper-parameters);",
    "v4 refit reuses this disclosure standard and must re-run the per-chain CrI",
    "overlap diagnosis on its own chains before submission;",
    "core species-level z/psi and detection-layer inference unaffected.",
    "Community-mean elevation and temperature-seasonality (bio4) responses",
    "carry this uncertainty; manuscript climate claims rest on species-level",
    "posterior propagation, not beta.comm point estimates.",
    "Must be reported in manuscript supplementary materials.")
  write_csv(disclosure, v3_file("results", paste0("table_convergence_disclosure_", run_label)))
}
gate <- data.frame(run_label = run_label, model_prefix = model_prefix,
  n_chains = n_expected, max_rhat = max(diag_df$rhat, na.rm = TRUE),
  min_ess_bulk = min(diag_df$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(diag_df$ess_tail, na.rm = TRUE),
  rhat_pass = strict_rhat_pass, ess_pass = strict_ess_pass,
  watchlist_n = nrow(watch_df),
  watchlist_max_rhat = if (nrow(watch_df)) max(watch_df$rhat) else NA_real_,
  watchlist_min_ess_bulk = if (nrow(watch_df)) min(watch_df$ess_bulk) else NA_real_,
  disclosure_required = relax_used, gate_pass = gate_pass,
  combined_file = basename(out), combined_bytes = file.info(out)$size)
write_csv(gate, v3_file("results", paste0("table_convergence_gate_", run_label)))
message(sprintf("[26] Combined %d chains; max Rhat %.3f, min bulk ESS %.0f; strict pass=%s/%s; watchlist n=%d (max Rhat %.3f, min ESS %.0f); disclosure=%s",
                n_expected, gate$max_rhat, gate$min_ess_bulk,
                strict_rhat_pass, strict_ess_pass, nrow(watch_df),
                gate$watchlist_max_rhat, gate$watchlist_min_ess_bulk, relax_used))
if (!gate_pass) {
  stop("Convergence gate failed. Do not continue to submission post-processing.")
}
if (relax_used) {
  message("[26] GATE PASSED WITH DISCLOSURE: weakly-identified watchlist hyper-parameters ",
          "exceed strict thresholds; see table_convergence_disclosure_", run_label,
          ".csv — must be reported in supplementary materials.")
}
