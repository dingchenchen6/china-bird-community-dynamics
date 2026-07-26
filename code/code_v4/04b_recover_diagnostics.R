#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   合并多链 stMsPGOcc 后验，验证模型收敛是否满足 R-hat ≤ 1.05、
#   ESS ≥ 200 的发表门槛；同步把 psi.samples 4D chain 维 flatten
#   进 draw 维，让下游 05 不再因维度错位崩溃（C3 修复）
#
# Objective / 分析目标:
#   - 多链合并：abind 沿 chain 维 → 单 fit 对象
#   - 单链 R-hat：返回 NA + warning（C2 修复）
#   - psi.samples chain flatten：[draw, sp, site, period, chain] → [draw*chain, sp, site, period]
#   - WAIC / 收敛诊断 / 群落 & 物种系数表 / 空间参数（phi, sigma.sq, rho）
#
# Input data / 输入数据:
#   data/derived_v4/stMsPGOcc_fit_<run_label>[_chainK].qs
#   data/derived_v4/psi_samples_thinned_<run_label>[_chainK].qs
#
# Main workflow / 主要流程:
#   1. 扫描链文件，分链/单文件 两种模式
#   2. 合并 fit；flatten psi.samples chain 维
#   3. 收敛诊断（utils_diagnostics::compute_rhat_safe / ESS）
#   4. 写出群落系数 / 物种系数 / 空间参数表
#   5. 输出 mcmc_diagnostics_summary CSV
#
# Key assumptions / 关键假设:
#   - spOccupancy single-chain fit slots: matrix [param, sample] / vector [sample]
#   - psi.samples 4D: [draws, species, sites, periods]，多链时是 5D
#
# Main packages / 主要包:
#   abind, qs, coda, dplyr, tidyr, tibble
#
# Output directory / 输出路径:
#   data/derived_v4/stMsPGOcc_fit_<run_label>_combined.qs
#   data/derived_v4/psi_samples_thinned_<run_label>_combined.qs
#   results_v4/table_beta_community_<run_label>.csv
#   results_v4/table_alpha_community_<run_label>.csv
#   results_v4/table_beta_species_<run_label>.csv
#   results_v4/table_spatial_params_<run_label>.csv
#   results_v4/table_convergence_diagnostics_<run_label>.csv
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(qs); library(abind)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_diversity.R"))   # summarise_post
source(file.path(CODE_V4, "utils_diagnostics.R"))
P <- ensure_v4_dirs()

log_time("04b", "Starting diagnostics recovery (v4)")

# ── 1. 加载 fit ──────────────────────────────────────────────────────
is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

chain_pattern <- paste0("stMsPGOcc_fit_", run_label, "_chain[0-9]+\\.qs$")
chain_files <- list.files(DIRS$derived, pattern = chain_pattern, full.names = TRUE)
if (length(chain_files) == 0) {
  # fallback: rds
  chain_pattern_rds <- paste0("stMsPGOcc_fit_", run_label, "_chain[0-9]+\\.rds$")
  chain_files <- list.files(DIRS$derived, pattern = chain_pattern_rds, full.names = TRUE)
}

single_qs  <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, ".qs"))
single_rds <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, ".rds"))

if (length(chain_files) >= 2) {
  chain_ids <- as.integer(gsub(".*chain([0-9]+)\\.(qs|rds)$", "\\1", basename(chain_files)))
  chain_files <- chain_files[order(chain_ids)]
  message(sprintf("[04b] Multi-chain mode: combining %d files", length(chain_files)))
  chain_fits <- lapply(chain_files, function(f) {
    message(sprintf("[04b]   Loading: %s", basename(f)))
    safe_read(f)
  })
  n_chains_actual <- length(chain_fits)
} else if (file.exists(single_qs)) {
  message(sprintf("[04b] Single fit (qs): %s", single_qs))
  chain_fits <- list(safe_read(single_qs))
  n_chains_actual <- 1L
} else if (file.exists(single_rds)) {
  message(sprintf("[04b] Single fit (rds): %s", single_rds))
  chain_fits <- list(safe_read(single_rds))
  n_chains_actual <- 1L
} else {
  stop("[04b] Model fit not found. Run 04 first.")
}

# ── 2. 合并多链（v4 改进：psi.samples 立刻 chain flatten） ──────────
combine_chains_v4 <- function(chain_fits) {
  if (length(chain_fits) == 0) stop("No chain fits")
  if (length(chain_fits) == 1) {
    fit <- chain_fits[[1]]
    fit$n.chains <- 1L
    return(fit)
  }

  combined <- chain_fits[[1]]
  sample_slots_matrix <- c(
    "beta.comm.samples", "alpha.comm.samples",
    "tau.sq.beta.samples", "tau.sq.alpha.samples",
    "tau.sq.samples", "phi.samples", "sigma.sq.samples", "rho.samples",
    "theta.samples", "lambda.samples"
  )

  for (slot in sample_slots_matrix) {
    if (slot %in% names(combined) && !is.null(combined[[slot]])) {
      combined[[slot]] <- .abind_chain_axis(
        lapply(chain_fits, function(f) f[[slot]]))
    }
  }

  # Species-level slots (list of matrix per species)
  for (slot in c("beta.samples", "alpha.samples")) {
    if (slot %in% names(combined) && is.list(combined[[slot]])) {
      n_sp <- length(combined[[slot]])
      for (i in seq_len(n_sp)) {
        combined[[slot]][[i]] <- .abind_chain_axis(
          lapply(chain_fits, function(f) f[[slot]][[i]]))
      }
    }
  }

  # psi.samples / z.samples：4D 后验
  for (slot in c("psi.samples", "z.samples")) {
    if (slot %in% names(combined) && !is.null(combined[[slot]])) {
      # 5D 合并（沿新维 chain），再立刻 flatten
      arr5d <- abind::abind(
        lapply(chain_fits, function(f) f[[slot]]),
        along = length(dim(combined[[slot]])) + 1L
      )
      combined[[slot]] <- flatten_chain_dim(arr5d)
      message(sprintf("[04b]   %s flattened: %s",
                      slot, paste(dim(combined[[slot]]), collapse = "x")))
    }
  }

  combined$n.chains <- length(chain_fits)
  combined$chain_ids <- seq_len(length(chain_fits))
  combined
}

.abind_chain_axis <- function(chain_samples) {
  chain_samples <- chain_samples[!sapply(chain_samples, is.null)]
  if (length(chain_samples) == 0) return(NULL)
  if (length(chain_samples) == 1) return(chain_samples[[1]])

  first <- chain_samples[[1]]
  if (is.vector(first)) {
    # vectors → cbind into [sample, chain]
    return(do.call(cbind, chain_samples))
  }
  if (is.matrix(first)) {
    return(abind::abind(chain_samples, along = 3))
  }
  if (length(dim(first)) >= 3) {
    return(abind::abind(chain_samples, along = length(dim(first)) + 1L))
  }
  first
}

fit <- combine_chains_v4(chain_fits)
rm(chain_fits); gc()

combined_qs <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, "_combined.qs"))
qs_save_safe(fit, combined_qs)

# ── 3. 同步合并 thinned psi.samples（已经在 04 写过单链版本） ────────
thin_pattern <- paste0("psi_samples_thinned_", run_label, "_chain[0-9]+\\.qs$")
thin_files <- list.files(DIRS$derived, pattern = thin_pattern, full.names = TRUE)
if (length(thin_files) >= 2) {
  thin_ids <- as.integer(gsub(".*chain([0-9]+)\\.qs$", "\\1", basename(thin_files)))
  thin_files <- thin_files[order(thin_ids)]
  thin_list <- lapply(thin_files, safe_read)
  # 检查 species/sites 一致性
  species <- thin_list[[1]]$species
  sites   <- thin_list[[1]]$sites
  # 沿 chain 维 abind，再 flatten
  psi_arrs <- lapply(thin_list, function(x) x$psi_samples_thinned)
  psi_5d <- abind::abind(psi_arrs, along = length(dim(psi_arrs[[1]])) + 1L)
  psi_flat <- flatten_chain_dim(psi_5d)
  qs_save_safe(list(
    psi_samples_thinned = psi_flat,
    species = species, sites = sites,
    psi_dim = dim(psi_flat),
    n_periods = thin_list[[1]]$n_periods,
    run_label = run_label,
    chain_id  = NA_integer_,
    n_chains_combined = length(thin_files)
  ), file.path(DIRS$derived, paste0("psi_samples_thinned_", run_label, "_combined.qs")))
  message(sprintf("[04b] Thinned psi combined: %s",
                  paste(dim(psi_flat), collapse = "x")))
}

# ── 4. 群落 / 物种系数表 ─────────────────────────────────────────────
beta_comm  <- summarise_post(fit$beta.comm.samples)
alpha_comm <- summarise_post(fit$alpha.comm.samples)

write_csv(beta_comm  |> mutate(param = "beta.comm"),
          v4_file("results", paste0("table_beta_community_", run_label)))
write_csv(alpha_comm |> mutate(param = "alpha.comm"),
          v4_file("results", paste0("table_alpha_community_", run_label)))

# 空间参数（如有）
sp_rows <- list()
for (slot in c("phi.samples", "sigma.sq.samples", "rho.samples", "tau.sq.samples")) {
  if (slot %in% names(fit) && !is.null(fit[[slot]])) {
    smry <- summarise_post(fit[[slot]])
    smry$param <- gsub("\\.samples$", "", slot)
    sp_rows[[slot]] <- smry
  }
}
if (length(sp_rows) > 0) {
  write_csv(bind_rows(sp_rows),
            v4_file("results", paste0("table_spatial_params_", run_label)))
}

# 种级系数
n_sp <- nrow(fit$y)
species <- rownames(fit$y) %||% paste0("sp", seq_len(n_sp))
beta_sp_list <- list()
for (i in seq_len(n_sp)) {
  bi <- summarise_post(fit$beta.samples[[i]])
  bi$species <- species[i]
  beta_sp_list[[i]] <- bi
}
write_csv(bind_rows(beta_sp_list),
          v4_file("results", paste0("table_beta_species_", run_label)))

# ── 5. 收敛诊断 ──────────────────────────────────────────────────────
diagnose_param_array_v4 <- function(samples, param_name, n_chains) {
  if (is.null(samples)) return(NULL)
  d <- dim(samples)

  if (is.null(d)) {
    # 纯 vector
    ess <- compute_ess_safe(samples)
    return(tibble(parameter = param_name, group = param_name,
                  rhat = NA_real_, ess = ess))
  }

  if (length(d) == 2) {
    # matrix: [param, sample] (单链) 或 [sample, chain] (vector 拼成)
    # 假设 [param, sample]
    n_params <- d[1]
    out <- vector("list", n_params)
    for (i in seq_len(n_params)) {
      vec <- samples[i, ]
      ess <- compute_ess_safe(vec)
      # 单链：R-hat = NA
      rhat <- if (n_chains > 1) NA_real_ else NA_real_
      lbl <- rownames(samples)[i] %||% paste0(param_name, "[", i, "]")
      out[[i]] <- tibble(parameter = lbl, group = param_name,
                          rhat = rhat, ess = ess)
    }
    return(bind_rows(out))
  }

  if (length(d) == 3) {
    # [param, sample, chain] —— 多链矩阵
    n_params <- d[1]; n_ch <- d[3]
    out <- vector("list", n_params)
    for (i in seq_len(n_params)) {
      chain_list <- lapply(seq_len(n_ch), function(ch) samples[i, , ch])
      rhat <- compute_rhat_safe(chain_list, warn_single = FALSE)
      ess  <- sum(sapply(chain_list, compute_ess_safe), na.rm = TRUE)
      lbl <- dimnames(samples)[[1]][i] %||% paste0(param_name, "[", i, "]")
      out[[i]] <- tibble(parameter = lbl, group = param_name,
                          rhat = rhat, ess = ess)
    }
    return(bind_rows(out))
  }

  # 4D / 5D: 简化只诊断第一维前 20 项
  n_dim1 <- d[1]
  out <- vector("list", min(n_dim1, 20))
  for (i in seq_len(min(n_dim1, 20))) {
    # 假设最后一维是 chain
    chain_list <- lapply(seq_len(d[length(d)]), function(ch) {
      idx <- as.list(rep(TRUE, length(d)))
      idx[[1]] <- i; idx[[length(d)]] <- ch
      as.vector(do.call(`[`, c(list(samples), idx)))
    })
    rhat <- compute_rhat_safe(chain_list, warn_single = FALSE)
    ess  <- sum(sapply(chain_list, compute_ess_safe), na.rm = TRUE)
    out[[i]] <- tibble(parameter = paste0(param_name, "[", i, "]"),
                       group = param_name,
                       rhat = rhat, ess = ess)
  }
  bind_rows(out)
}

diag_list <- list()
for (slot in c("beta.comm.samples", "alpha.comm.samples",
               "phi.samples", "sigma.sq.samples", "rho.samples")) {
  if (slot %in% names(fit)) {
    d <- diagnose_param_array_v4(fit[[slot]], gsub("\\.samples$", "", slot),
                                  fit$n.chains %||% n_chains_actual)
    if (!is.null(d) && nrow(d) > 0) diag_list[[slot]] <- d
  }
}

if (length(diag_list) > 0) {
  diag_df <- bind_rows(diag_list)
  write_csv(diag_df, v4_file("results",
                              paste0("table_convergence_diagnostics_", run_label)))

  # 汇总
  conv_summary <- summarise_convergence(diag_df,
                                         rhat_thresh = RHAT_THRESHOLD,
                                         ess_thresh  = ESS_THRESHOLD)
  write_csv(conv_summary, v4_file("results",
                                    paste0("table_convergence_summary_", run_label)))

  if (n_chains_actual >= 2) {
    max_rhat <- max(diag_df$rhat, na.rm = TRUE)
    min_ess  <- min(diag_df$ess, na.rm = TRUE)
    n_bad_rhat <- sum(diag_df$rhat > RHAT_THRESHOLD, na.rm = TRUE)
    n_bad_ess  <- sum(diag_df$ess < ESS_THRESHOLD, na.rm = TRUE)
    message(sprintf("[04b] Convergence: max R-hat=%.3f (%d > %.2f), min ESS=%g (%d < %d)",
                    max_rhat, n_bad_rhat, RHAT_THRESHOLD, min_ess, n_bad_ess, ESS_THRESHOLD))
    if (max_rhat > 1.10) {
      warning("[04b] Some R-hat > 1.10 — consider longer chains")
    } else if (max_rhat > RHAT_THRESHOLD) {
      warning("[04b] Some R-hat > ", RHAT_THRESHOLD, " — may need more iterations")
    } else {
      message("[04b] All R-hat OK")
    }
  } else {
    message("[04b] Single chain → R-hat = NA (use multi-chain for convergence assessment)")
  }
}

# WAIC
waic_val <- tryCatch(spOccupancy::waicOcc(fit)$waic, error = function(e) NA_real_)
message(sprintf("[04b] WAIC = %.2f", waic_val))

log_time("04b", sprintf("DONE: %s, WAIC=%.2f, n_chains=%d",
                        run_label, waic_val, n_chains_actual))
