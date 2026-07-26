#!/usr/bin/env Rscript
## 05_postprocess_diversity.R
##
## 阶段 5：把 stage-4 主拟合的 psi 后验（thinned）传播到群落多样性、
## temporal beta、driver 回归（brms+cmdstanr+4chain+DHARMa）。
##
## 关键改进相对 v1:
##   1) 真正"逐 MCMC draw"算 richness/Shannon/CWM/CWSD/Rao's Q/PD/MPD，
##      落盘 mean + 95% CRI，而不是只在后验均值上算一次。
##   2) Faith's PD 用 probability-weighted（Reese-Tucker），不二值化 psi。
##   3) Temporal beta 同时输出 Bray、Sørensen、Baselga turnover/nestedness。
##   4) Driver 回归用 brms + cmdstanr + 4 chains，DHARMa 残差诊断。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(sf)
  library(ape)
  library(phangorn)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_diversity.R"))
source(file.path(CODE_V2, "utils_diagnostics.R"))

P <- ensure_v2_dirs()
set.seed(20260502)

RUN_LABEL  <- Sys.getenv("V2_RUN_LABEL", "v2_pilot_60sp_ar1")
N_DRAWS_USE <- as.integer(Sys.getenv("V2_POST_DRAWS", "150"))   # 上限节流

message(sprintf("[stage-5] postprocess for %s | draws cap = %d",
                RUN_LABEL, N_DRAWS_USE))

## --- 1. 加载 stage-4 输出 + 元数据 -----------------------------------------

psi_thin_obj <- readRDS(v2_file("derived",
                                 paste0("psi_samples_thinned_", RUN_LABEL),
                                 "rds"))
psi_arr <- psi_thin_obj$psi_samples_thinned   # [draws, species, sites, periods]
n_draws_avail <- dim(psi_arr)[1]
sel_draws <- if (n_draws_avail > N_DRAWS_USE)
  seq(1, n_draws_avail, length.out = N_DRAWS_USE) |> round() else
  seq_len(n_draws_avail)
psi_arr <- psi_arr[sel_draws, , , , drop = FALSE]
n_draws  <- dim(psi_arr)[1]

# 物种 / 站点 / 主期标签
candidate_all <- read_csv(v2_file("results",
                            "table_dynamic_occupancy_candidate_species_all"),
                            show_col_types = FALSE)
primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
visit_effort <- readRDS(file.path(P$derived_v2, "visit_effort_2000_2024.rds"))
sites <- sort(unique(visit_effort$grid_cell))
n_species <- dim(psi_arr)[2]
candidate_species <- head(candidate_all$species, n_species)
n_sites <- dim(psi_arr)[3]
n_periods <- dim(psi_arr)[4]

dimnames(psi_arr) <- list(
  draw   = paste0("d", seq_len(n_draws)),
  species = candidate_species,
  grid    = as.character(sites),
  period  = primary_blocks$block_label
)

message(sprintf("  psi tensor: %d draws × %d sp × %d sites × %d periods",
                n_draws, n_species, n_sites, n_periods))

## --- 2. 性状 / 系统发育 ---------------------------------------------------

trait_path <- file.path(P$project_root, "bird_grid_community_analysis",
                        "results", "table_species_traits_imputed.csv")
trait_tbl <- read_csv(trait_path, show_col_types = FALSE)
trait_vars <- c("body_mass_g", "clutch_size", "longevity_y",
                 "maturity_y", "avonet_hwi", "avonet_range_size")
trait_vars <- intersect(trait_vars, names(trait_tbl))

# 严格按 candidate_species 顺序对齐，缺失物种补 NA 行（后面用列中位数 fill）。
trait_aligned <- tibble(species = candidate_species) |>
  left_join(
    trait_tbl |>
      select(species, all_of(trait_vars)) |>
      mutate(across(where(is.numeric),
                    ~ if_else(.x > 0, log10(.x), .x))),
    by = "species"
  )
trait_mat <- trait_aligned |> select(-species) |> as.matrix()
rownames(trait_mat) <- trait_aligned$species

# 缺失 trait 物种用列中位数填补，避免 NA 传染
for (j in seq_len(ncol(trait_mat))) {
  miss <- is.na(trait_mat[, j])
  if (any(miss)) trait_mat[miss, j] <- median(trait_mat[, j], na.rm = TRUE)
}
# 标准化
trait_mat <- scale(trait_mat)
attr(trait_mat, "scaled:center") <- NULL
attr(trait_mat, "scaled:scale") <- NULL

# Gower-like Euclidean distance on standardized traits
dist_mat <- as.matrix(stats::dist(trait_mat))

# 系统发育：clootl 需要先显式加载内部数据到全局 env
suppressMessages(library(clootl))
data("clootl_data", package = "clootl", envir = globalenv())
phy <- NULL; phy_dist <- NULL; species_match <- NULL
phy_obj <- tryCatch(
  clootl::extractTree(species = candidate_species,
                       label_type = "scientific",
                       taxonomy_year = 2025,
                       version = "1.6",
                       force = TRUE),
  error = function(e) { message("  clootl::extractTree failed: ",
                                conditionMessage(e)); NULL }
)
if (!is.null(phy_obj)) {
  phy <- phy_obj
  phy_dist <- ape::cophenetic.phylo(phy)
  species_match <- match(candidate_species, phy$tip.label)
  message(sprintf("  phylogeny: %d / %d species matched",
                  sum(!is.na(species_match)), length(candidate_species)))
} else {
  message("  no phylogeny - PD/MPD will be NA")
}

## --- 3. 逐 draw 计算群落指标 ----------------------------------------------

message("[stage-5] computing per-draw community metrics")
# 输出维度：[draws, sites, periods, metric]
metric_names_taxon <- c("corrected_richness", "shannon",
                         "inverse_simpson", "effective_species")
cwm_names  <- paste0("cwm_",  colnames(trait_mat))
cwsd_names <- paste0("cwsd_", colnames(trait_mat))
metric_names_func <- c(cwm_names, cwsd_names,
                        "trait_volume", "trait_dispersion", "rao_q")
metric_names_phy  <- c("pd_prob", "mpd_prob")
metric_names <- c(metric_names_taxon, metric_names_func, metric_names_phy)
n_metrics <- length(metric_names)

metric_arr <- array(NA_real_,
                     dim = c(n_draws, n_sites, n_periods, n_metrics),
                     dimnames = list(draw = NULL,
                                     grid = as.character(sites),
                                     period = primary_blocks$block_label,
                                     metric = metric_names))

t0 <- Sys.time()
for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    for (t in seq_len(n_periods)) {
      psi_v <- psi_arr[d, , s, t]
      tax <- div_taxonomic(psi_v)
      fn  <- div_functional(psi_v, trait_mat)
      raq <- div_rao_q(psi_v, dist_mat)
      pd  <- div_pd_prob(psi_v, phy, species_match)
      mpv <- div_mpd_prob(psi_v, phy_dist, species_match)
      metric_arr[d, s, t, "corrected_richness"] <- tax["corrected_richness"]
      metric_arr[d, s, t, "shannon"]            <- tax["shannon"]
      metric_arr[d, s, t, "inverse_simpson"]    <- tax["inverse_simpson"]
      metric_arr[d, s, t, "effective_species"]  <- tax["effective_species"]
      metric_arr[d, s, t, cwm_names]  <- fn$cwm
      metric_arr[d, s, t, cwsd_names] <- fn$cwsd
      metric_arr[d, s, t, "trait_volume"]     <- fn$trait_volume
      metric_arr[d, s, t, "trait_dispersion"] <- fn$trait_dispersion
      metric_arr[d, s, t, "rao_q"]            <- raq
      metric_arr[d, s, t, "pd_prob"]          <- pd
      metric_arr[d, s, t, "mpd_prob"]         <- mpv
    }
  }
  if (d %% max(1, floor(n_draws / 10)) == 0)
    message(sprintf("  draw %d / %d (%.1f min)", d, n_draws,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## --- 4. 后验 mean + 95% CRI 汇总 ------------------------------------------

agg_post <- function(arr, FUN, ...) apply(arr, c(2, 3, 4), FUN, ...)
metric_mean   <- agg_post(metric_arr, mean, na.rm = TRUE)
metric_l95    <- agg_post(metric_arr, quantile, probs = 0.025, na.rm = TRUE)
metric_u95    <- agg_post(metric_arr, quantile, probs = 0.975, na.rm = TRUE)
metric_sd     <- agg_post(metric_arr, sd, na.rm = TRUE)

# 整成 long table
to_long <- function(arr, suffix) {
  ar <- as.data.frame.table(arr, responseName = paste0("value", suffix),
                             stringsAsFactors = FALSE)
  names(ar)[1:3] <- c("grid_cell", "block_label", "metric")
  ar
}
metrics_long <- list(
  to_long(metric_mean, "_mean"),
  to_long(metric_l95,  "_l95"),
  to_long(metric_u95,  "_u95"),
  to_long(metric_sd,   "_sd")
) |>
  reduce(\(a, b) full_join(a, b, by = c("grid_cell", "block_label", "metric"))) |>
  mutate(grid_cell = as.integer(as.character(grid_cell)))

write_csv(metrics_long,
          v2_file("results",
                  paste0("table_community_metrics_with_cri_", RUN_LABEL)))

## --- 5. Temporal beta（per draw → mean + CRI）-----------------------------

message("[stage-5] computing temporal beta with Baselga decomposition")
beta_metric_names <- c("bray", "soerensen", "turnover", "nestedness",
                        "delta_richness")
n_pairs <- n_periods - 1
beta_arr <- array(NA_real_,
                   dim = c(n_draws, n_sites, n_pairs, length(beta_metric_names)),
                   dimnames = list(draw = NULL,
                                   grid = as.character(sites),
                                   pair = paste0(primary_blocks$block_label[-1],
                                                 " vs ",
                                                 primary_blocks$block_label[-n_periods]),
                                   metric = beta_metric_names))
for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    for (k in seq_len(n_pairs)) {
      v <- div_temporal_beta_pair(psi_arr[d, , s, k], psi_arr[d, , s, k + 1])
      beta_arr[d, s, k, ] <- v
    }
  }
}
beta_mean <- apply(beta_arr, c(2, 3, 4), mean, na.rm = TRUE)
beta_l95  <- apply(beta_arr, c(2, 3, 4), quantile, 0.025, na.rm = TRUE)
beta_u95  <- apply(beta_arr, c(2, 3, 4), quantile, 0.975, na.rm = TRUE)

beta_long <- list(
  to_long(beta_mean, "_mean"),
  to_long(beta_l95,  "_l95"),
  to_long(beta_u95,  "_u95")
) |>
  reduce(\(a, b) full_join(a, b, by = c("grid_cell", "block_label", "metric"))) |>
  rename(period_pair = block_label) |>
  mutate(grid_cell = as.integer(as.character(grid_cell)))
write_csv(beta_long,
          v2_file("results", paste0("table_temporal_beta_with_cri_", RUN_LABEL)))

## --- 6. Per-grid linear trend across 5 periods (per draw) -----------------

message("[stage-5] computing per-grid linear trends with CRI")
period_idx <- seq_len(n_periods) - mean(seq_len(n_periods))
denom <- sum(period_idx^2)

# slope = sum_t (t - tbar) * y[t] / denom
trend_per_draw <- function(y_mat) {
  # y_mat: [n_periods, K]
  apply(y_mat, 2, function(y) sum(period_idx * y, na.rm = TRUE) / denom)
}

# Output: [n_draws, n_sites, n_metrics] of slope
trend_arr <- array(NA_real_,
                    dim = c(n_draws, n_sites, n_metrics),
                    dimnames = list(draw = NULL,
                                    grid = as.character(sites),
                                    metric = metric_names))
for (d in seq_len(n_draws)) {
  for (m in seq_len(n_metrics)) {
    yt <- t(metric_arr[d, , , m])     # [periods, sites]
    trend_arr[d, , m] <- trend_per_draw(yt)
  }
}
trend_mean <- apply(trend_arr, c(2, 3), mean, na.rm = TRUE)
trend_l95  <- apply(trend_arr, c(2, 3), quantile, 0.025, na.rm = TRUE)
trend_u95  <- apply(trend_arr, c(2, 3), quantile, 0.975, na.rm = TRUE)
trend_sig  <- (apply(trend_arr, c(2, 3), function(x) mean(x > 0, na.rm = TRUE)))

to_long2 <- function(mat, val_name) {
  d <- as.data.frame.table(mat, responseName = val_name, stringsAsFactors = FALSE)
  names(d)[1:2] <- c("grid_cell", "metric")
  d$grid_cell <- as.integer(as.character(d$grid_cell))
  as_tibble(d)
}
trend_long <- to_long2(trend_mean, "trend_mean") |>
  left_join(to_long2(trend_l95, "trend_l95"),
            by = c("grid_cell", "metric")) |>
  left_join(to_long2(trend_u95, "trend_u95"),
            by = c("grid_cell", "metric")) |>
  left_join(to_long2(trend_sig, "p_pos"),
            by = c("grid_cell", "metric")) |>
  mutate(sig95 = (trend_l95 > 0) | (trend_u95 < 0))
write_csv(trend_long,
          v2_file("results", paste0("table_grid_trends_with_cri_", RUN_LABEL)))

## --- 7. Driver 回归：brms + cmdstanr + DHARMa -----------------------------

message("[stage-5] driver regressions (brms + cmdstanr + 4 chains)")
grid_env <- readRDS(file.path(P$derived_v2, "grid_environment_dynamic_occupancy.rds"))

# 把 trend 表 pivot 到 wide：每个 metric 一列趋势值
trend_wide_mean <- trend_long |>
  select(grid_cell, metric, trend_mean) |>
  pivot_wider(names_from = metric, values_from = trend_mean,
              names_prefix = "trend_")

driver_input <- trend_wide_mean |>
  inner_join(grid_env, by = "grid_cell")

# 标准化驱动变量
drv_vars <- c("bio4", "elev_sd", "texture_shannon", "npp_mean",
              "landcover_built", "habitat_diversity_shannon", "hfi_mean",
              "centroid_lat") |>
  intersect(names(driver_input))
for (v in drv_vars)
  driver_input[[paste0("z_", v)]] <- as.numeric(scale(driver_input[[v]]))

resp_vars <- c("trend_corrected_richness", "trend_shannon",
               "trend_trait_volume", "trend_pd_prob")
resp_vars <- intersect(resp_vars, names(driver_input))

driver_results <- list()
for (resp in resp_vars) {
  message("  fitting brms model for ", resp)
  use <- driver_input |>
    select(all_of(c(resp, paste0("z_", drv_vars),
                     "centroid_lon", "centroid_lat"))) |>
    drop_na()
  if (nrow(use) < 100) {
    message("    skip (n<100)"); next
  }
  rhs <- paste0("z_", drv_vars, collapse = " + ")
  form <- as.formula(sprintf("%s ~ %s", resp, rhs))
  fit_b <- tryCatch(
    brms::brm(form, data = use, family = gaussian(),
              chains = 4, iter = 2000, warmup = 1000,
              cores = 4, backend = "cmdstanr",
              refresh = 0,
              save_pars = brms::save_pars(all = TRUE)),
    error = function(e) {
      message("    brms failed: ", conditionMessage(e)); NULL })
  if (is.null(fit_b)) next

  # diagnostics
  diag_b <- brms_diag_report(fit_b, paste0("driver_", resp, "_", RUN_LABEL))
  dh <- tryCatch(dharma_from_brms(fit_b, n_sim = 500,
                                   run_label = paste0("driver_", resp, "_", RUN_LABEL)),
                 error = function(e) {
                   message("    DHARMa failed: ", conditionMessage(e)); NULL })

  # forest plot of standardized coefficients
  post <- as.data.frame(brms::fixef(fit_b, probs = c(0.025, 0.975)))
  post$term <- rownames(post)
  fp <- post |>
    filter(term != "Intercept") |>
    mutate(term = forcats::fct_reorder(term, Estimate))
  fp_plot <- ggplot(fp, aes(Estimate, term, xmin = `Q2.5`, xmax = `Q97.5`)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbar(orientation = "y", width = 0, linewidth = 0.5,
                  colour = "#0E5A78") +
    geom_point(size = 2.4, colour = "#0E5A78") +
    labs(title = sprintf("Drivers of %s", resp),
         subtitle = sprintf("brms + cmdstanr | n=%d grids", nrow(use)),
         x = "Standardised coefficient (95% CRI)", y = NULL) +
    theme_v2_pub(11)
  save_dual(fp_plot,
            paste0("fig_driver_forest_", resp, "_", RUN_LABEL),
            width = 8.4, height = 5.2)
  driver_results[[resp]] <- list(fit = fit_b, fixef = post)
  saveRDS(fit_b,
          v2_file("derived",
                  paste0("brms_driver_fit_", resp, "_", RUN_LABEL),
                  "rds"))
}

## --- 8. 输出说明 -----------------------------------------------------------

writeLines(c(
  sprintf("# v2 阶段 5 后处理说明（%s）", RUN_LABEL),
  "",
  sprintf("- thinned posterior draws used: %d", n_draws),
  sprintf("- 物种 / 站点 / 主期 / 配对: %d / %d / %d / %d",
          n_species, n_sites, n_periods, n_pairs),
  sprintf("- trait variables used: %s", paste(colnames(trait_mat), collapse = ", ")),
  if (!is.null(phy)) sprintf("- phylogeny matched species: %d / %d",
                              sum(!is.na(species_match)), n_species)
                     else "- phylogeny: not available",
  sprintf("- driver responses fitted (brms+cmdstanr): %s",
          paste(names(driver_results), collapse = ", "))
), v2_file("results",
           paste0("postprocess_summary_", RUN_LABEL), "md"))

message("[stage-5] Done.")
