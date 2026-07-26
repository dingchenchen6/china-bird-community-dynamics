## utils_diagnostics.R
## v2 模型诊断工具：MCMC 收敛 + 后验预测检验 + DHARMa 残差。
##
## 设计：
##   - spOccupancy 主模型（tMsPGOcc）：用 spOccupancy::ppcOcc + 自写 R-hat/ESS 提取（基于 coda）。
##   - brms 驱动回归：DHARMa::createDHARMa(posterior_predict, observed) 走标准残差诊断。
##   - 输出统一图（PNG + PDF）+ 表（CSV），所有结果落在 figures_v2/ 与 results_v2/。

suppressPackageStartupMessages({
  library(coda)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
})

if (!exists("v2_script_dir", mode = "function")) {
  source(Sys.getenv("V2_PATHS_R",
    "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2/utils_paths.R"))
}

## --- 1. spOccupancy MCMC 诊断 ----------------------------------------------

## extract_mcmc_diag: 给定 spOccupancy 拟合对象 + 关心的样本名，返回 R-hat / ESS。
## spOccupancy 把 n.chains 条链的样本沿行拼成一个 mcmc 矩阵；
## 这里按 fit$n.chains 切回 mcmc.list 才能算 R-hat。
extract_mcmc_diag <- function(fit, sample_slot) {
  s <- fit[[sample_slot]]
  if (is.null(s)) return(NULL)
  M <- as.matrix(s)
  n_chains <- if (!is.null(fit$n.chains)) fit$n.chains else 1L
  per_chain <- nrow(M) %/% max(1L, n_chains)
  mc <- if (n_chains > 1 && per_chain > 1) {
    coda::mcmc.list(lapply(seq_len(n_chains), function(c) {
      i0 <- (c - 1L) * per_chain + 1L
      i1 <- c * per_chain
      coda::mcmc(M[i0:i1, , drop = FALSE])
    }))
  } else {
    coda::mcmc.list(coda::mcmc(M))
  }
  rhat <- tryCatch(
    coda::gelman.diag(mc, autoburnin = FALSE, multivariate = FALSE)$psrf,
    error = function(e) NULL
  )
  ess <- tryCatch(coda::effectiveSize(mc), error = function(e) NULL)
  tibble(
    sample_slot = sample_slot,
    parameter   = colnames(M),
    rhat        = if (!is.null(rhat)) rhat[, "Point est."] else NA_real_,
    ess         = if (!is.null(ess)) as.numeric(ess) else NA_real_
  )
}

write_mcmc_diag_report <- function(fit, run_label,
                                   slots = c("beta.comm.samples", "alpha.comm.samples",
                                             "tau.sq.beta.samples", "tau.sq.alpha.samples",
                                             "sigma.sq.t.samples", "rho.samples")) {
  rows <- lapply(slots, function(s) extract_mcmc_diag(fit, s))
  diag_tbl <- bind_rows(Filter(Negate(is.null), rows))
  out_csv <- v2_file("results", paste0("mcmc_diagnostics_", run_label), "csv")
  write_csv(diag_tbl, out_csv)
  diag_tbl
}

## plot_rhat_ess: 把诊断表变成两面板的 R-hat / ESS 直方图。
plot_rhat_ess <- function(diag_tbl, run_label) {
  source(file.path(v2_paths()$code_v2, "utils_mapping.R"))
  p_rhat <- ggplot(diag_tbl, aes(rhat)) +
    geom_histogram(bins = 40, fill = "#0E5A78", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = 1.05, linetype = 2, colour = "#8B2E1E") +
    labs(title = "R-hat distribution", subtitle = "Threshold 1.05 (red dashed)",
         x = "R-hat", y = "# parameters") +
    theme_v2_pub(10.5)
  p_ess <- ggplot(diag_tbl, aes(ess)) +
    geom_histogram(bins = 40, fill = "#3C8C5A", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = 400, linetype = 2, colour = "#8B2E1E") +
    labs(title = "Effective sample size", subtitle = "Threshold 400 (red dashed)",
         x = "ESS", y = "# parameters") +
    theme_v2_pub(10.5)
  patchwork::wrap_plots(p_rhat, p_ess, ncol = 2)
}

## --- 2. spOccupancy 后验预测检验 -------------------------------------------

run_ppc_summary <- function(fit, group = 1, type = "freeman-tukey", run_label) {
  ppc <- spOccupancy::ppcOcc(fit, fit.stat = type, group = group)
  bayes_p <- if (is.list(ppc) && !is.null(ppc$fit.y) && !is.null(ppc$fit.y.rep)) {
    mean(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)
  } else NA_real_
  res <- list(ppc = ppc, bayes_p = bayes_p)
  saveRDS(res, v2_file("derived", paste0("ppc_", run_label), "rds"))
  res
}

## --- 3. brms + DHARMa 包装 -------------------------------------------------

dharma_from_brms <- function(brms_fit, n_sim = 1000, integer_response = FALSE,
                             run_label = "model") {
  if (!requireNamespace("DHARMa", quietly = TRUE))
    stop("DHARMa is required for residual diagnostics.")
  pp  <- brms::posterior_predict(brms_fit, ndraws = n_sim)   # [draws x n_obs]
  mf  <- stats::model.frame(brms_fit)
  obs <- mf[[1]]                                              # response vec
  if (length(obs) != ncol(pp))
    stop(sprintf("DHARMa size mismatch: obs=%d, pp ncol=%d",
                 length(obs), ncol(pp)))
  fit_med <- apply(pp, 2, median)
  sim <- DHARMa::createDHARMa(
    simulatedResponse  = t(pp),
    observedResponse   = obs,
    fittedPredictedResponse = fit_med,
    integerResponse    = integer_response
  )
  png_path <- v2_file("figure", paste0("dharma_", run_label), "png")
  grDevices::png(png_path, width = 1700, height = 900, res = 220)
  plot(sim)
  grDevices::dev.off()
  saveRDS(sim, v2_file("derived", paste0("dharma_sim_", run_label), "rds"))
  list(sim = sim, png = png_path)
}

## --- 4. brms 收敛 + LOO ---------------------------------------------------

brms_diag_report <- function(brms_fit, run_label) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms required.")
  sm <- as.data.frame(brms::posterior_summary(brms_fit))
  sm$parameter <- rownames(sm)
  rhat_tbl <- tibble(
    parameter = sm$parameter,
    rhat      = brms::rhat(brms_fit),
    ess_bulk  = brms::neff_ratio(brms_fit) * sum(brms::ndraws(brms_fit))
  )
  loo_obj <- tryCatch(brms::loo(brms_fit, moment_match = FALSE), error = function(e) NULL)
  write_csv(rhat_tbl, v2_file("results", paste0("brms_diag_", run_label), "csv"))
  if (!is.null(loo_obj)) saveRDS(loo_obj, v2_file("derived", paste0("brms_loo_", run_label), "rds"))
  list(rhat = rhat_tbl, loo = loo_obj)
}
