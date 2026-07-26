#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   修复 v3 04b_recover_diagnostics.R 中：
#     (a) 单链 R-hat 假报 1.0（C2）
#     (b) psi.samples 4D chain 维未 flatten（C3 在 04b 中已处理；这里
#         提供工具函数）
#     (c) DHARMa gate 对 brms 模型的统一接口
#
# Objective / 分析目标:
#   - compute_rhat_safe(chain_list)         : 单链返回 NA + warning
#   - compute_ess_safe(samples)              : coda::effectiveSize 安全包装
#   - flatten_chain_dim(samples_4d)          : 把 chain 维 fold 进 draw 维
#   - check_dharma_gate(brms_fit, ...)       : 沿用 v3，加强容错
#   - summarise_convergence(diag_df, ...)    : 阈值汇总
#
# Input data / 输入数据:
#   chain_list / array / brms fit
#
# Main workflow / 主要流程:
#   纯函数定义
#
# Key assumptions / 关键假设:
#   coda / DHARMa 可用（按需 require）
#
# Main packages / 主要包:
#   coda, DHARMa
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

# ── compute_rhat_safe（C2 修复：单链返回 NA） ────────────────────────
compute_rhat_safe <- function(chain_list, warn_single = TRUE) {
  chain_list <- chain_list[!sapply(chain_list, is.null)]
  if (length(chain_list) < 2) {
    if (warn_single) {
      warning("[rhat] single chain detected; R-hat undefined → returning NA",
              call. = FALSE)
    }
    return(NA_real_)
  }
  min_len <- min(sapply(chain_list, length))
  chains_matrix <- sapply(chain_list, function(ch) ch[1:min_len])
  n <- min_len; m <- length(chain_list)
  chain_means <- colMeans(chains_matrix)
  chain_vars  <- apply(chains_matrix, 2, var)
  B <- n * var(chain_means)
  W <- mean(chain_vars)
  if (W == 0) return(NA_real_)
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  max(sqrt(var_hat / W), 1)
}

# ── compute_ess_safe ─────────────────────────────────────────────────
compute_ess_safe <- function(samples) {
  if (!requireNamespace("coda", quietly = TRUE)) return(NA_real_)
  tryCatch(as.numeric(coda::effectiveSize(samples)),
           error = function(e) NA_real_)
}

# ── flatten_chain_dim（C3 修复：把 chain 维 fold 进 draw 维） ────────
#' 把 4D / 5D 后验样本（多了一维 chain）合并到 draw 维。
#' 输入：array(dim = c(draws_per_chain, ..., chain))
#' 输出：array(dim = c(draws_per_chain * chain, ...))
#'
#' 例：psi.samples 输入维度 [draws, species, sites, periods, chain]
#'     → 输出 [draws*chain, species, sites, periods]
flatten_chain_dim <- function(arr, chain_axis = NULL) {
  d <- dim(arr)
  if (is.null(d) || length(d) < 2) return(arr)

  # 默认最后一维是 chain
  if (is.null(chain_axis)) chain_axis <- length(d)

  if (chain_axis == length(d) && d[chain_axis] == 1) {
    # 只有 1 chain → 直接 drop chain 维
    return(array(arr, dim = d[-chain_axis]))
  }

  # 移动 chain 维到最前
  perm <- c(chain_axis, setdiff(seq_along(d), chain_axis))
  arr_p <- aperm(arr, perm)
  # 现在第一维 = chain，第二维 = draw_per_chain
  d_p <- dim(arr_p)
  # reshape: flatten 前两维
  out_dim <- c(d_p[1] * d_p[2], d_p[-c(1, 2)])
  array(as.numeric(arr_p), dim = out_dim)
}

# ── summarise_convergence（基于诊断 df 做阈值汇总） ─────────────────
summarise_convergence <- function(diag_df, rhat_thresh = 1.05, ess_thresh = 200) {
  diag_df |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n_params       = dplyr::n(),
      rhat_max       = suppressWarnings(max(rhat, na.rm = TRUE)),
      rhat_gt_thresh = sum(rhat > rhat_thresh, na.rm = TRUE),
      ess_min        = suppressWarnings(min(ess, na.rm = TRUE)),
      ess_lt_thresh  = sum(ess < ess_thresh, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      rhat_ok = rhat_gt_thresh == 0 & !is.na(rhat_max) & is.finite(rhat_max),
      ess_ok  = ess_lt_thresh  == 0 & !is.na(ess_min) & is.finite(ess_min),
      overall = rhat_ok & ess_ok
    )
}

# ── extract_species_diag（沿用 v3 接口） ──────────────────────────────
extract_species_diag <- function(fit, n_species = NULL, species_names = NULL) {
  rhat_out <- fit$rhat
  ess_out  <- fit$ESS
  if (is.null(rhat_out) || is.null(ess_out)) {
    warning("Model fit missing rhat/ESS slots. Skip.")
    return(NULL)
  }
  if (is.null(n_species)) {
    n_species <- if (!is.null(fit$y)) nrow(fit$y) else NROW(rhat_out)
  }
  if (is.null(species_names)) {
    species_names <- if (!is.null(rownames(fit$y))) rownames(fit$y) else paste0("sp", seq_len(n_species))
  }
  param_names <- rownames(rhat_out)
  if (is.null(param_names)) param_names <- paste0("param", seq_len(nrow(rhat_out)))
  diag_list <- list()
  for (p in seq_along(param_names)) {
    diag_list[[p]] <- tibble::tibble(
      species = species_names,
      param   = param_names[p],
      rhat    = as.numeric(rhat_out[p, ]),
      ess     = as.numeric(ess_out[p, ])
    )
  }
  dplyr::bind_rows(diag_list)
}

# ── check_dharma_gate（v3 沿用 + 容错强化） ──────────────────────────
check_dharma_gate <- function(brms_fit, save_dir = NULL, stem = "dharma",
                                n_sim = 250) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    warning("DHARMa package not installed. Skipping gate check.")
    return(NULL)
  }
  sim_res <- tryCatch(
    DHARMa::simulateResiduals(brms_fit, n = n_sim, plot = FALSE),
    error = function(e) {
      message("[dharma] simulateResiduals failed: ", e$message)
      NULL
    }
  )
  if (is.null(sim_res)) return(NULL)

  od_test  <- tryCatch(DHARMa::testDispersion(sim_res, plot = FALSE), error = function(e) NULL)
  unif_test <- tryCatch(DHARMa::testUniformity(sim_res, plot = FALSE), error = function(e) NULL)
  out_test  <- tryCatch(DHARMa::testOutliers(sim_res, plot = FALSE), error = function(e) NULL)

  result <- tibble::tibble(
    test = c("overdispersion", "uniformity", "outliers"),
    statistic = c(od_test$statistic %||% NA_real_,
                   unif_test$statistic %||% NA_real_,
                   out_test$statistic %||% NA_real_),
    p_value = c(od_test$p.value %||% NA_real_,
                 unif_test$p.value %||% NA_real_,
                 out_test$p.value %||% NA_real_)
  ) |>
    dplyr::mutate(passed = !is.na(p_value) & p_value > 0.05)

  # 保存图
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    png_path <- file.path(save_dir, paste0(stem, ".png"))
    tryCatch({
      png(png_path, width = 12, height = 5, units = "in", res = 300)
      par(mfrow = c(1, 2))
      plot(sim_res)
      dev.off()
      message("[dharma] saved → ", png_path)
    }, error = function(e) message("[dharma] plot save failed: ", e$message))
  }

  all_passed <- all(result$passed, na.rm = TRUE)
  if (!all_passed) {
    failed <- result$test[!result$passed]
    warning("[dharma_gate] FAILED: ", paste(failed, collapse = ", "), call. = FALSE)
  } else {
    message("[dharma_gate] PASSED all tests")
  }
  result
}

message("[utils_diagnostics_v4] loaded — single-chain R-hat = NA, chain flatten ready")
