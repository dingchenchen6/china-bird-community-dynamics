#!/usr/bin/env Rscript
## 04_run_tMsPGOcc_main.R
##
## 阶段 4：spOccupancy::tMsPGOcc 多物种动态占域主拟合（合并+去重数据 v2）。
##
## 关键设计：
##   - 4 chains
##   - 默认 PILOT (60 species, 短链) 用于端到端验证；正式 run 用 env var 触发：
##         V2_RUN_MODE=full V2_SPECIES_LIMIT=200 V2_N_BATCH=200 V2_N_BURN=2500 ...
##   - 保存 thinned posterior samples（5x thin 后的 psi.samples 子集 → 给 stage-5 算 CRI）
##   - 输出：MCMC 收敛诊断（R-hat / ESS）、trace / density 图、PPC、coefficient summary
##
## 输入：data/derived_v2/{visit_effort_2000_2024.rds, species_visit_2000_2024.rds,
##                       grid_environment_dynamic_occupancy.rds}
##       results_v2/table_dynamic_occupancy_candidate_species_all.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(spOccupancy)
  library(coda)
  library(ggplot2)
  library(patchwork)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_diagnostics.R"))

P <- ensure_v2_dirs()
set.seed(20260502)

## --- 0. 运行参数（PILOT vs FULL） ------------------------------------------

RUN_MODE <- Sys.getenv("V2_RUN_MODE", "pilot")          # pilot | full
SPECIES_LIMIT <- as.integer(Sys.getenv("V2_SPECIES_LIMIT",
                                        if (RUN_MODE == "pilot") "60" else "200"))
N_BATCH       <- as.integer(Sys.getenv("V2_N_BATCH",
                                        if (RUN_MODE == "pilot") "30" else "200"))
BATCH_LENGTH  <- as.integer(Sys.getenv("V2_BATCH_LENGTH", "25"))
N_BURN        <- as.integer(Sys.getenv("V2_N_BURN",
                                        if (RUN_MODE == "pilot") "300" else "2500"))
N_THIN        <- as.integer(Sys.getenv("V2_N_THIN",
                                        if (RUN_MODE == "pilot") "2" else "5"))
N_CHAINS      <- as.integer(Sys.getenv("V2_N_CHAINS", "4"))
AR1_FLAG      <- identical(Sys.getenv("V2_AR1", "true"), "true")
RUN_LABEL     <- Sys.getenv("V2_RUN_LABEL",
  sprintf("v2_%s_%dsp_%s", RUN_MODE, SPECIES_LIMIT,
          if (AR1_FLAG) "ar1" else "noar1"))

message(sprintf("[stage-4] RUN_MODE=%s | label=%s", RUN_MODE, RUN_LABEL))
message(sprintf("  species=%d | n.batch=%d | batch.length=%d | n.burn=%d | n.thin=%d | n.chains=%d | ar1=%s",
                SPECIES_LIMIT, N_BATCH, BATCH_LENGTH, N_BURN, N_THIN,
                N_CHAINS, AR1_FLAG))

## --- 1. 输入加载 -----------------------------------------------------------

visit_effort   <- readRDS(file.path(P$derived_v2, "visit_effort_2000_2024.rds"))
species_visit  <- readRDS(file.path(P$derived_v2, "species_visit_2000_2024.rds"))
grid_env       <- readRDS(file.path(P$derived_v2, "grid_environment_dynamic_occupancy.rds"))
candidate_all  <- read_csv(v2_file("results", "table_dynamic_occupancy_candidate_species_all"),
                            show_col_types = FALSE)
primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)

candidate_species <- head(candidate_all$species, SPECIES_LIMIT)
sites <- sort(unique(visit_effort$grid_cell))
site_ids <- as.character(sites)
n_species  <- length(candidate_species)
n_sites    <- length(sites)
n_primary  <- nrow(primary_blocks)
n_secondary <- max(visit_effort$year_in_block, na.rm = TRUE)

message(sprintf("  building y array: species=%d × sites=%d × periods=%d × reps=%d",
                n_species, n_sites, n_primary, n_secondary))

## --- 2. 构建 y / det.covs / occ.covs --------------------------------------

y <- array(NA_integer_,
            dim = c(n_species, n_sites, n_primary, n_secondary),
            dimnames = list(candidate_species, site_ids,
                            primary_blocks$block_label,
                            paste0("rep", seq_len(n_secondary))))
visited_mask     <- array(FALSE, dim = c(n_sites, n_primary, n_secondary))
log_events_arr   <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
log_observers_arr <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
duration_arr     <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))

eff <- visit_effort |> mutate(site_index = match(grid_cell, sites))
for (ii in seq_len(nrow(eff))) {
  rr <- eff[ii, ]
  visited_mask[rr$site_index, rr$block_id, rr$year_in_block] <- TRUE
  log_events_arr[rr$site_index, rr$block_id, rr$year_in_block]    <- rr$log_events
  log_observers_arr[rr$site_index, rr$block_id, rr$year_in_block] <- rr$log_observers
  duration_arr[rr$site_index, rr$block_id, rr$year_in_block]      <- rr$mean_duration_min
}

det <- species_visit |>
  filter(species %in% candidate_species) |>
  mutate(species_index = match(species, candidate_species),
         site_index    = match(grid_cell, sites))
for (ii in seq_len(nrow(det))) {
  rr <- det[ii, ]
  y[rr$species_index, rr$site_index, rr$block_id, rr$year_in_block] <- rr$detected
}

# 在所有访问过但未检测到该物种的位置补 0
visited_idx <- which(visited_mask, arr.ind = TRUE)
for (ii in seq_len(nrow(visited_idx))) {
  j <- visited_idx[ii, 1]; t <- visited_idx[ii, 2]; k <- visited_idx[ii, 3]
  miss <- is.na(y[, j, t, k])
  if (any(miss)) y[miss, j, t, k] <- 0L
}

# det.covs 缺失填中位数（仅在 visited 位置）
fill_arr_na <- function(arr, mask) {
  obs <- is.finite(arr)
  fill_val <- if (any(obs)) median(arr[obs], na.rm = TRUE) else 0
  arr[mask & !obs] <- fill_val
  arr
}
log_events_arr     <- fill_arr_na(log_events_arr,     visited_mask)
log_observers_arr  <- fill_arr_na(log_observers_arr,  visited_mask)
duration_arr       <- fill_arr_na(duration_arr,       visited_mask)
duration_arr[!is.finite(duration_arr)] <- 0      # 极少数全 NA 网格

## --- 3. occ.covs 与协变量筛选（vif + cor） --------------------------------

standardize <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(0, length(x)))
  x[!is.finite(x)] <- median(x[is.finite(x)], na.rm = TRUE)
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

grid_env_aligned <- grid_env |>
  filter(grid_cell %in% sites) |>
  arrange(match(grid_cell, sites))

candidate_vars <- c(
  paste0("bio", 1:19),
  "elev_mean", "elev_sd",
  "texture_shannon", "texture_entropy", "texture_contrast",
  "npp_mean",
  "landcover_trees", "landcover_cropland", "landcover_built",
  "landcover_shrubs", "landcover_grassland", "landcover_water",
  "habitat_diversity_shannon",
  "natural_landcover_fraction", "human_modified_fraction",
  "hfi_mean", "hfi_sd",
  "centroid_lat", "centroid_lon"
) |> intersect(names(grid_env_aligned))
mandatory <- intersect(c("elev_mean", "elev_sd", "centroid_lat",
                         "npp_mean", "hfi_mean"),
                       candidate_vars)

screen <- {
  df <- grid_env_aligned |> mutate(across(all_of(candidate_vars), standardize))
  keep <- candidate_vars
  cor_cutoff <- 0.7
  vif_cutoff <- 5
  log_drops <- tibble(step = character(), variable = character(),
                      reason = character())
  repeat {
    cm <- suppressWarnings(cor(df[, keep, drop = FALSE],
                               use = "pairwise.complete.obs",
                               method = "spearman"))
    diag(cm) <- 0
    mx <- suppressWarnings(max(abs(cm), na.rm = TRUE))
    if (!is.finite(mx) || mx <= cor_cutoff) break
    pair <- which(abs(cm) == mx, arr.ind = TRUE)[1, ]
    pair_vars <- rownames(cm)[pair]
    drop_cands <- setdiff(pair_vars, mandatory)
    if (length(drop_cands) == 0) break
    mab <- vapply(drop_cands, function(v)
      mean(abs(cm[v, setdiff(keep, v)]), na.rm = TRUE), numeric(1))
    drop <- drop_cands[which.max(mab)]
    keep <- setdiff(keep, drop)
    log_drops <- bind_rows(log_drops,
      tibble(step = "correlation", variable = drop,
             reason = sprintf("|spearman r|>%.2f", cor_cutoff)))
    if (length(keep) <= 2) break
  }
  repeat {
    if (length(keep) <= 2) break
    vif_vec <- vapply(keep, function(v) {
      others <- setdiff(keep, v)
      fit <- lm(reformulate(others, response = v), data = df)
      r2 <- summary(fit)$r.squared
      if (r2 >= 0.999) Inf else 1 / (1 - r2)
    }, numeric(1))
    cand <- vif_vec[!names(vif_vec) %in% mandatory]
    if (!length(cand) || max(cand) <= vif_cutoff) break
    drop <- names(cand)[which.max(cand)]
    keep <- setdiff(keep, drop)
    log_drops <- bind_rows(log_drops,
      tibble(step = "vif", variable = drop,
             reason = sprintf("VIF=%.2f > %.0f", max(cand), vif_cutoff)))
  }
  list(keep = keep, log = log_drops)
}
occ_keep <- screen$keep
write_csv(
  bind_rows(
    screen$log |> mutate(scope = "occupancy"),
    tibble(step = "retained", variable = occ_keep,
           reason = "kept after cor/vif", scope = "occupancy")
  ),
  v2_file("results", paste0("table_environment_screening_", RUN_LABEL))
)
message("  retained occupancy predictors: ", paste(occ_keep, collapse = ", "))

site_cov <- function(v) standardize(grid_env_aligned[[v]])
occ_covs <- setNames(lapply(occ_keep, site_cov), occ_keep)
occ_covs$year_scaled <- matrix(
  rep(as.numeric(scale(primary_blocks$block_start)), each = n_sites),
  nrow = n_sites, ncol = n_primary
)
occ_formula <- reformulate(c(occ_keep, "year_scaled"))
det_covs <- list(log_events    = log_events_arr,
                 log_observers = log_observers_arr,
                 duration_min  = duration_arr)

## --- 4. tMsPGOcc 拟合 ------------------------------------------------------

z_init <- apply(y, c(1, 2, 3), function(a) as.integer(any(a == 1, na.rm = TRUE)))
inits <- list(alpha.comm = 0, beta.comm = 0, beta = 0, alpha = 0,
              tau.sq.beta = 1, tau.sq.alpha = 1, z = z_init)
priors <- list(beta.comm.normal  = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig    = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig   = list(a = 0.1, b = 0.1))
tuning <- list()
if (AR1_FLAG) {
  inits$sigma.sq.t  <- rep(1, n_species)
  inits$rho         <- rep(0.2, n_species)
  priors$sigma.sq.t.ig <- list(a = 0.1, b = 0.1)
  priors$rho.unif      <- list(a = -1, b = 1)
  tuning$rho           <- 0.2
} else {
  tuning$phi <- 1
}

fit_path <- v2_file("derived", paste0("tMsPGOcc_fit_", RUN_LABEL), "rds")
if (file.exists(fit_path) && identical(Sys.getenv("V2_REUSE_FIT", "false"), "true")) {
  message("  reusing existing fit at ", fit_path)
  fit <- readRDS(fit_path)
} else {
  t0 <- Sys.time()
  fit <- spOccupancy::tMsPGOcc(
    occ.formula  = occ_formula,
    det.formula  = ~ log_events + log_observers + duration_min,
    data         = list(y = y, occ.covs = occ_covs, det.covs = det_covs),
    inits        = inits,
    priors       = priors,
    tuning       = tuning,
    n.batch      = N_BATCH,
    batch.length = BATCH_LENGTH,
    accept.rate  = 0.43,
    n.omp.threads = 1,
    verbose      = TRUE,
    ar1          = AR1_FLAG,
    n.report     = max(1, floor(N_BATCH / 10)),
    n.burn       = N_BURN,
    n.thin       = N_THIN,
    n.chains     = N_CHAINS
  )
  t1 <- Sys.time()
  attr(fit, "v2_elapsed_min") <- as.numeric(difftime(t1, t0, units = "mins"))
  saveRDS(fit, fit_path, compress = FALSE)
  message(sprintf("  fit done in %.2f min, written to %s",
                  attr(fit, "v2_elapsed_min"), fit_path))
}

## --- 5. 收敛诊断（R-hat / ESS） -------------------------------------------

message("[stage-4] running convergence diagnostics")
diag_tbl <- write_mcmc_diag_report(fit, RUN_LABEL,
  slots = c("beta.comm.samples", "alpha.comm.samples",
            "tau.sq.beta.samples", "tau.sq.alpha.samples"))
if (AR1_FLAG) {
  diag_tbl <- bind_rows(
    diag_tbl,
    extract_mcmc_diag(fit, "sigma.sq.t.samples"),
    extract_mcmc_diag(fit, "rho.samples")
  )
  write_csv(diag_tbl,
            v2_file("results", paste0("mcmc_diagnostics_", RUN_LABEL)))
}

p_diag <- plot_rhat_ess(diag_tbl, RUN_LABEL) +
  patchwork::plot_annotation(
    title = sprintf("MCMC convergence diagnostics - %s", RUN_LABEL),
    subtitle = sprintf("4 chains | n.batch=%d × batch.length=%d, burn=%d, thin=%d",
                       N_BATCH, BATCH_LENGTH, N_BURN, N_THIN),
    theme = theme_v2_pub(11)
  )
save_dual(p_diag, paste0("fig_mcmc_diag_", RUN_LABEL),
          width = 11, height = 4.6)

## --- 6. trace 图（采样核心参数子集） --------------------------------------

trace_params <- c("(Intercept)", occ_keep, "year_scaled")
trace_df <- as.matrix(fit$beta.comm.samples) |>
  as.data.frame() |>
  mutate(iter = row_number()) |>
  pivot_longer(-iter, names_to = "parameter", values_to = "value") |>
  filter(parameter %in% trace_params)

trace_plot <- ggplot(trace_df, aes(iter, value)) +
  geom_line(colour = "#0E5A78", linewidth = 0.25, alpha = 0.85) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  labs(title = sprintf("Community-level beta posterior trace - %s", RUN_LABEL),
       x = "Iteration (post-burn, post-thin)", y = "beta value") +
  theme_v2_pub(10)
save_dual(trace_plot, paste0("fig_mcmc_trace_betacomm_", RUN_LABEL),
          width = 11, height = max(5, ceiling(length(trace_params) / 3) * 1.6))

## --- 7. 群落系数 caterpillar plot -----------------------------------------

beta_comm <- as.matrix(fit$beta.comm.samples)
alpha_comm <- as.matrix(fit$alpha.comm.samples)

summarise_post <- function(M, kind) {
  tibble(
    parameter = colnames(M),
    mean      = apply(M, 2, mean),
    median    = apply(M, 2, median),
    l95       = apply(M, 2, quantile, 0.025),
    u95       = apply(M, 2, quantile, 0.975),
    p_pos     = apply(M, 2, function(x) mean(x > 0)),
    kind      = kind
  )
}

post_summary <- bind_rows(
  summarise_post(beta_comm, "occupancy"),
  summarise_post(alpha_comm, "detection")
) |>
  mutate(sig95 = (l95 > 0) | (u95 < 0))
write_csv(post_summary,
          v2_file("results", paste0("table_community_coefficients_", RUN_LABEL)))

caterpillar_plot <- post_summary |>
  mutate(parameter = forcats::fct_reorder(parameter, mean)) |>
  ggplot(aes(mean, parameter, xmin = l95, xmax = u95, colour = sig95)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#8B2E1E"),
                       labels = c("CI incl. 0", "CI excl. 0"),
                       name = NULL) +
  facet_wrap(~ kind, scales = "free", ncol = 2) +
  labs(title = sprintf("Community-level coefficients - %s", RUN_LABEL),
       subtitle = "Posterior mean +/- 95% credible interval (4 chains pooled)",
       x = "Coefficient (logit scale)", y = NULL) +
  theme_v2_pub(11)
save_dual(caterpillar_plot,
          paste0("fig_community_caterpillar_", RUN_LABEL),
          width = 11, height = 5.6)

## --- 8. 物种 × 协变量系数 heatmap （hierarchical shrinkage 视图） ---------

beta_sp <- as.matrix(fit$beta.samples)
sp_long <- as.data.frame(beta_sp) |>
  mutate(iter = row_number()) |>
  pivot_longer(-iter, names_to = "term", values_to = "value") |>
  separate(term, into = c("covariate", "species"), sep = "-",
           extra = "merge", remove = FALSE) |>
  group_by(covariate, species) |>
  summarise(mean = mean(value), l95 = quantile(value, 0.025),
            u95 = quantile(value, 0.975),
            .groups = "drop") |>
  mutate(sig95 = (l95 > 0) | (u95 < 0))
write_csv(sp_long,
          v2_file("results", paste0("table_species_coefficients_", RUN_LABEL)))

heat_plot <- sp_long |>
  ggplot(aes(covariate, species, fill = mean)) +
  geom_tile(colour = "white", linewidth = 0.05) +
  geom_text(aes(label = ifelse(sig95, "*", "")),
            size = 2.6, colour = "grey15") +
  scale_fill_v2_diverging(name = "Posterior mean beta",
                          limits = c(-2, 2)) +
  labs(title = sprintf("Species × covariate occupancy effects - %s", RUN_LABEL),
       subtitle = "Color = posterior mean beta; * = 95% CRI excludes 0",
       x = NULL, y = NULL) +
  theme_v2_pub(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = max(4, 9 - n_species / 30)))
save_dual(heat_plot,
          paste0("fig_species_coef_heatmap_", RUN_LABEL),
          width = 11, height = max(6, n_species * 0.10))

## --- 9. 后验预测检验 ------------------------------------------------------

message("[stage-4] running posterior predictive check (Freeman-Tukey)")
# 内存防护：先用 group=2（聚合到主期）+ 子采样限制 draws；失败再 fallback。
do_ppc <- function(grp) {
  spOccupancy::ppcOcc(fit, fit.stat = "freeman-tukey", group = grp)
}
ppc_res <- tryCatch(do_ppc(2),
  error = function(e) { message("  ppc group=2 failed: ", conditionMessage(e))
                         tryCatch(do_ppc(1),
                                   error = function(e2) {
                                     message("  ppc group=1 also failed: ",
                                             conditionMessage(e2)); NULL })
  })
if (!is.null(ppc_res)) {
  saveRDS(ppc_res, v2_file("derived", paste0("ppc_", RUN_LABEL), "rds"))
  bp <- mean(rowMeans(ppc_res$fit.y.rep, na.rm = TRUE) >
             rowMeans(ppc_res$fit.y, na.rm = TRUE), na.rm = TRUE)
  message(sprintf("  Bayesian p-value (Freeman-Tukey) = %.3f", bp))
  ppc_df <- tibble(
    fit_y = rowMeans(ppc_res$fit.y, na.rm = TRUE),
    fit_y_rep = rowMeans(ppc_res$fit.y.rep, na.rm = TRUE)
  )
  ppc_plot <- ggplot(ppc_df, aes(fit_y, fit_y_rep)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(alpha = 0.3, size = 1.2, colour = "#0E5A78") +
    labs(title = sprintf("Posterior predictive check - %s", RUN_LABEL),
         subtitle = sprintf("Freeman-Tukey | Bayesian p-value = %.3f", bp),
         x = "fit.y (observed)", y = "fit.y.rep (replicated)") +
    theme_v2_pub(11)
  save_dual(ppc_plot, paste0("fig_ppc_", RUN_LABEL),
            width = 6.6, height = 5.6)
}

## --- 10. 落盘 thinned posterior（给 stage-5 算 CRI） ----------------------

# psi.samples 维度 [chain*draws_per_chain, species, sites, periods]
# 强制把 thinning 上限到 ≤ 200 draws（足以算 95% CRI），控制 RDS 体积。
PSI_MAX_DRAWS <- as.integer(Sys.getenv("V2_PSI_MAX_DRAWS", "200"))
psi_samps <- fit$psi.samples
n_draws <- dim(psi_samps)[1]
keep_draws <- if (n_draws > PSI_MAX_DRAWS)
  round(seq(1, n_draws, length.out = PSI_MAX_DRAWS)) else
  seq_len(n_draws)
psi_thin <- psi_samps[keep_draws, , , , drop = FALSE]
saveRDS(list(psi_samples_thinned = psi_thin,
              keep_draws = keep_draws,
              dims = dim(psi_samps)),
        v2_file("derived", paste0("psi_samples_thinned_", RUN_LABEL), "rds"))
message(sprintf("  psi.samples thinned: %d draws kept of %d (cap=%d)",
                length(keep_draws), n_draws, PSI_MAX_DRAWS))

## --- 11. 写运行说明 -------------------------------------------------------

writeLines(c(
  sprintf("# v2 阶段 4 主拟合说明（%s）", RUN_LABEL),
  "",
  sprintf("- run_mode: %s", RUN_MODE),
  sprintf("- 物种数: %d", n_species),
  sprintf("- 站点数: %d", n_sites),
  sprintf("- 主期 / 重复: %d / %d", n_primary, n_secondary),
  sprintf("- MCMC: n.batch=%d, batch.length=%d, n.burn=%d, n.thin=%d, n.chains=%d",
          N_BATCH, BATCH_LENGTH, N_BURN, N_THIN, N_CHAINS),
  sprintf("- AR1: %s", AR1_FLAG),
  sprintf("- 协变量: %s + year_scaled", paste(occ_keep, collapse = ", ")),
  sprintf("- 拟合用时: %.2f min", attr(fit, "v2_elapsed_min") %||% NA_real_),
  if (exists("ppc_res") && !is.null(ppc_res))
    sprintf("- Bayesian p-value (Freeman-Tukey): %.3f", bp) else
    "- PPC: 未成功（见日志）"
), v2_file("results", paste0("multispecies_dynamic_main_summary_", RUN_LABEL), "md"))

message("[stage-4] Done.")
