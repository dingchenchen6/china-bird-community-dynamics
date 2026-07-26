#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把 stMsPGOcc 后验 psi 一致地传播到群落（分类、系统发育、功能、
#   时间动态、Baselga 周转）多样性，并对 5 主期趋势做稳健估计与
#   驱动因子贝叶斯回归（Q1/Q2/Q3/Q5 核心）
#
# Objective / 分析目标:
#   v4 关键改进相对 v3：
#     (C1) library(loo)；loo::loo(...) 替代 lo(...)
#     (C4) 用 div_phylogenetic_prob（真正概率加权 PD）
#     (C6) brms 驱动回归用 horseshoe 稀疏先验 + GP k 三档敏感性
#     (C7) 只读 psi_samples_thinned_*_combined.qs，不读完整 fit
#     新增：CWM 空间表（feed 06 的 F 段）
#     新增：每个驱动回归同时报告 marginal_R² vs conditional_R²
#
# Input data / 输入数据:
#   data/derived_v4/psi_samples_thinned_<run_label>_combined.qs
#   data/derived_v4/grid_environment{GRID_TAG}_v4.rds  (or v3 fallback)
#   data/derived_v4/trait_extended_v4.rds              (or v3 fallback)
#   data/derived_v4/phylogeny_matched_v4.rds           (or v3 fallback)
#
# Main workflow / 主要流程:
#   1. 加载 thinned psi（[draws, sp, sites, periods]）
#   2. 逐 draw × period × site 计算 10+ 个多样性指标
#   3. 计算时间动态：synchrony / variance ratio / turnover
#   4. 概率加权 Baselga 分解（相邻 period 对）
#   5. 趋势：OLS + Theil-Sen + Mann-Kendall
#   6. Naive vs corrected richness 对比
#   7. brms 驱动回归（horseshoe + GP k 三档）
#   8. 物种层趋势 + 性状关联 + 环境关联 + 热点
#   9. CWM 空间表（每网格 × period × trait）
#
# Key assumptions / 关键假设:
#   psi.samples 已经 chain-flattened（由 04b 输出 _combined.qs）
#
# Main packages / 主要包:
#   brms, cmdstanr, loo, qs, dplyr, tidyr, abind
#
# Output directory / 输出路径:
#   results_v4/table_diversity_summary_<run_label>.csv
#   results_v4/table_temporal_dynamics_summary_<run_label>.csv
#   results_v4/table_baselga_summary_<run_label>.csv
#   results_v4/table_baselga_global_<run_label>.csv
#   results_v4/table_trend_summary_<run_label>.csv
#   results_v4/table_mann_kendall_<run_label>.csv
#   results_v4/table_naive_vs_corrected_<run_label>.csv
#   results_v4/table_species_trend_<run_label>.csv
#   results_v4/table_species_trend_classify_<run_label>.csv
#   results_v4/table_species_hotspot_<run_label>.csv
#   results_v4/table_species_env_trend_<run_label>.csv
#   results_v4/table_trait_trend_correlation_<run_label>.csv
#   results_v4/table_env_trend_correlation_<run_label>.csv
#   results_v4/table_cwm_spatial_pattern_<run_label>.csv
#   results_v4/table_brms_driver_coefs_<metric>_<gp_k>_<run_label>.csv
#   results_v4/table_brms_driver_R2_<run_label>.csv
#   results_v4/table_loo_driver_<run_label>.csv
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(brms); library(cmdstanr); library(loo); library(abind); library(qs)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_seeds.R"))
source(file.path(CODE_V4, "utils_diversity.R"))
source(file.path(CODE_V4, "utils_diagnostics.R"))
P <- ensure_v4_dirs()

set_seeds("05_postprocess")
log_time("05", "Starting postprocess diversity (v4)")

# ── 1. 加载 thinned psi + 环境 + 性状 + 系统发育 ─────────────────────
is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

thin_paths <- c(
  file.path(DIRS$derived, paste0("psi_samples_thinned_", run_label, "_combined.qs")),
  file.path(DIRS$derived, paste0("psi_samples_thinned_", run_label, ".qs"))
)
thin <- NULL
for (p in thin_paths) if (is.null(thin) && file.exists(p)) thin <- safe_read(p)
if (is.null(thin)) stop("[05] Thinned psi not found. Run 04 (and 04b for multi-chain) first.")

psi_samples <- thin$psi_samples_thinned
species <- thin$species
sites   <- thin$sites
psi_dim <- dim(psi_samples)
n_draws  <- psi_dim[1]
n_sp     <- psi_dim[2]
n_sites  <- psi_dim[3]
n_periods <- if (length(psi_dim) >= 4) psi_dim[4] else 1L
message(sprintf("[05] psi: %d draws × %d sp × %d sites × %d periods",
                n_draws, n_sp, n_sites, n_periods))

grid_env <- safe_read(file.path(DIRS$derived, paste0("grid_environment", GRID_TAG, "_v4.rds")), quiet = TRUE) %||%
            safe_read(file.path(DIRS$v3_derived, "grid_environment_v3.rds"), quiet = TRUE) %||%
            safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"), quiet = TRUE)
if (is.null(grid_env)) warning("[05] grid_environment not found; driver regression will be skipped")

trait_ext <- safe_read(file.path(DIRS$derived, "trait_extended_v4.rds"), quiet = TRUE) %||%
             safe_read(file.path(DIRS$v3_derived, "trait_extended_v3.rds"), quiet = TRUE)
if (is.null(trait_ext)) {
  warning("[05] trait_extended not found, using basic trait fallback")
  trait_ext <- read_csv_safe(TRAIT_IMPUTED_PATH)
  if (!is.null(trait_ext)) {
    trait_ext$diet_specialization <- NA_real_
    trait_ext$habitat_breadth    <- NA_real_
  }
}

phylo <- safe_read(file.path(DIRS$derived, "phylogeny_matched_v4.rds"), quiet = TRUE) %||%
         safe_read(file.path(DIRS$v3_derived, "phylogeny_matched_v3.rds"), quiet = TRUE) %||%
         safe_read(file.path(DIRS$v2_derived, "phylogeny_matched.rds"), quiet = TRUE)

# ── 2. 性状矩阵 ──────────────────────────────────────────────────────
trait_aligned <- tibble(species = species) |> left_join(trait_ext, by = "species")
all_trait_vars <- intersect(TRAIT_VARS_ALL, names(trait_aligned))
trait_mat <- prepare_trait_matrix(trait_aligned, all_trait_vars)

# ── 3. 多样性指标（10 类） ──────────────────────────────────────────
metrics <- c("corrected_richness", "shannon", "inv_simpson",
             "trait_volume", "trait_dispersion", "rao_q",
             "feve", "fdiv",
             "pd_prob", "mpd_prob")
metric_arr <- array(NA_real_,
                    dim = c(n_draws, n_sites, n_periods, length(metrics)),
                    dimnames = list(NULL, sites, NULL, metrics))

message(sprintf("[05] Computing diversity: %d × %d × %d cells",
                n_draws, n_sites, n_periods))

for (d in seq_len(n_draws)) {
  for (t in seq_len(n_periods)) {
    psi_d <- if (length(psi_dim) >= 4) psi_samples[d, , , t] else psi_samples[d, , ]
    for (s in seq_len(n_sites)) {
      psi_v <- psi_d[, s]
      tax <- div_taxonomic(psi_v)
      metric_arr[d, s, t, "corrected_richness"] <- tax$richness
      metric_arr[d, s, t, "shannon"]            <- tax$shannon
      metric_arr[d, s, t, "inv_simpson"]        <- tax$inv_simpson

      fn <- div_functional(psi_v, trait_mat)
      metric_arr[d, s, t, "trait_volume"]     <- fn$trait_volume
      metric_arr[d, s, t, "trait_dispersion"] <- fn$trait_dispersion
      metric_arr[d, s, t, "rao_q"]            <- fn$rao_q

      metric_arr[d, s, t, "feve"] <- div_feve(psi_v, trait_mat)
      metric_arr[d, s, t, "fdiv"] <- div_fdiv(psi_v, trait_mat)

      if (!is.null(phylo)) {
        pd <- div_phylogenetic_prob(psi_v, phylo, tip_order = species)
        metric_arr[d, s, t, "pd_prob"]  <- pd$pd_prob
        metric_arr[d, s, t, "mpd_prob"] <- pd$mpd_prob
      }
    }
  }
  if (d %% 50 == 0) message(sprintf("  [diversity] draw %d/%d", d, n_draws))
}

# ── 4. 多样性后验汇总 ────────────────────────────────────────────────
summary_long <- list()
for (m in metrics) {
  for (t in seq_len(n_periods)) {
    smry <- summarise_post(metric_arr[, , t, m])
    smry$metric    <- m
    smry$period    <- paste0("P", t)
    smry$grid_cell <- sites
    summary_long[[length(summary_long) + 1]] <- smry
  }
}
metric_summary_df <- bind_rows(summary_long)
write_csv(metric_summary_df,
          v4_file("results", paste0("table_diversity_summary_", run_label)))

# 06 兼容长格式
metric_cri_long <- metric_summary_df |>
  select(grid_cell, period, metric, value_mean = mean,
         value_l95 = q025, value_u95 = q975) |>
  mutate(block_label = period)
write_csv(metric_cri_long,
          v4_file("results", paste0("table_community_metrics_with_cri_", run_label)))

# ── 5. 时间动态（synchrony / variance ratio / turnover） ────────────
temporal_metrics <- c("synchrony", "variance_ratio",
                      "turnover_total", "turnover_gain", "turnover_loss")
temporal_arr <- array(NA_real_,
                      dim = c(n_draws, n_sites, length(temporal_metrics)),
                      dimnames = list(NULL, sites, temporal_metrics))

if (n_periods >= 2 && length(psi_dim) >= 4) {
  for (d in seq_len(n_draws)) {
    for (s in seq_len(n_sites)) {
      comm_mat <- psi_samples[d, , s, ]
      if (!is.matrix(comm_mat)) comm_mat <- matrix(comm_mat, nrow = n_sp, ncol = n_periods)
      if (ncol(comm_mat) < 2) next
      temporal_arr[d, s, "synchrony"]       <- div_synchrony(comm_mat)
      temporal_arr[d, s, "variance_ratio"]  <- div_variance_ratio(comm_mat)
      to <- div_temporal_turnover(comm_mat, threshold = 0.1)
      temporal_arr[d, s, "turnover_total"]  <- to$turnover_total
      temporal_arr[d, s, "turnover_gain"]   <- to$turnover_gain
      temporal_arr[d, s, "turnover_loss"]   <- to$turnover_loss
    }
    if (d %% 50 == 0) message(sprintf("  [temporal] draw %d/%d", d, n_draws))
  }
  temporal_summary <- list()
  for (m in temporal_metrics) {
    smry <- summarise_post(temporal_arr[, , m])
    smry$metric    <- m
    smry$grid_cell <- sites
    temporal_summary[[length(temporal_summary) + 1]] <- smry
  }
  write_csv(bind_rows(temporal_summary),
            v4_file("results", paste0("table_temporal_dynamics_summary_", run_label)))
}

# ── 6. Baselga 分解 ─────────────────────────────────────────────────
baselga_metrics <- c("beta_sor", "beta_sim", "beta_sne", "prop_turnover")
n_pp <- n_periods - 1
if (n_pp >= 1 && length(psi_dim) >= 4) {
  baselga_arr <- array(NA_real_,
                       dim = c(n_draws, n_sites, n_pp, length(baselga_metrics)),
                       dimnames = list(NULL, sites, NULL, baselga_metrics))
  for (d in seq_len(n_draws)) {
    for (s in seq_len(n_sites)) {
      comm_mat <- psi_samples[d, , s, ]
      if (!is.matrix(comm_mat)) comm_mat <- matrix(comm_mat, nrow = n_sp, ncol = n_periods)
      if (ncol(comm_mat) < 2) next
      for (pp in seq_len(n_pp)) {
        bg <- div_baselga(comm_mat[, pp], comm_mat[, pp + 1])
        baselga_arr[d, s, pp, "beta_sor"]      <- bg$beta_sor
        baselga_arr[d, s, pp, "beta_sim"]      <- bg$beta_sim
        baselga_arr[d, s, pp, "beta_sne"]      <- bg$beta_sne
        baselga_arr[d, s, pp, "prop_turnover"] <- bg$prop_turnover
      }
    }
    if (d %% 50 == 0) message(sprintf("  [baselga] draw %d/%d", d, n_draws))
  }
  baselga_rows <- list()
  for (pp in seq_len(n_pp)) {
    for (m in baselga_metrics) {
      smry <- summarise_post(baselga_arr[, , pp, m])
      smry$metric      <- m
      smry$period_pair <- paste0("P", pp, "_P", pp + 1)
      smry$grid_cell   <- sites
      baselga_rows[[length(baselga_rows) + 1]] <- smry
    }
  }
  baselga_summary_df <- bind_rows(baselga_rows)
  write_csv(baselga_summary_df,
            v4_file("results", paste0("table_baselga_summary_", run_label)))

  global_baselga <- baselga_summary_df |>
    filter(metric == "prop_turnover") |>
    group_by(period_pair) |>
    summarise(
      prop_turnover_mean = mean(mean, na.rm = TRUE),
      prop_turnover_q025 = mean(q025, na.rm = TRUE),
      prop_turnover_q975 = mean(q975, na.rm = TRUE),
      n_grids = sum(!is.na(mean)),
      .groups = "drop"
    )
  write_csv(global_baselga,
            v4_file("results", paste0("table_baselga_global_", run_label)))
}

# ── 7. 趋势（OLS + Theil-Sen） ──────────────────────────────────────
trend_methods <- c("ols", "theil_sen")
trend_arr <- array(NA_real_,
                   dim = c(n_draws, n_sites, length(metrics), length(trend_methods)),
                   dimnames = list(NULL, sites, metrics, trend_methods))

for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    for (m in metrics) {
      vals <- metric_arr[d, s, , m]
      if (sum(!is.na(vals)) >= 3) {
        trend_arr[d, s, m, "ols"] <- coef(lm(vals ~ seq_len(n_periods)))[2]
        trend_arr[d, s, m, "theil_sen"] <- theil_sen_slope(vals)
      }
    }
  }
}

trend_rows <- list()
for (method in trend_methods) {
  for (m in metrics) {
    smry <- summarise_post(trend_arr[, , m, method])
    smry$metric    <- m
    smry$method    <- method
    smry$grid_cell <- sites
    trend_rows[[length(trend_rows) + 1]] <- smry
  }
}
trend_summary_df <- bind_rows(trend_rows)
write_csv(trend_summary_df,
          v4_file("results", paste0("table_trend_summary_", run_label)))

# Mann-Kendall
mk_rows <- list()
for (s in seq_len(n_sites)) {
  for (m in metrics) {
    vals_mean <- apply(metric_arr[, s, , m], 2, mean, na.rm = TRUE)
    if (sum(!is.na(vals_mean)) >= 4) {
      mk <- mann_kendall_test(vals_mean)
      mk_rows[[length(mk_rows) + 1]] <- tibble(
        grid_cell = sites[s], metric = m,
        mk_S = mk$S, mk_tau = mk$tau, mk_p = mk$p_value
      )
    }
  }
}
write_csv(bind_rows(mk_rows),
          v4_file("results", paste0("table_mann_kendall_", run_label)))

# ── 8. CWM 空间表（v4 新增，供 06 F 段） ────────────────────────────
message("[05] Computing CWM spatial pattern")
cwm_rows <- list()
trait_focus <- intersect(c("diet_specialization", "habitat_breadth"),
                          colnames(trait_mat))
if (length(trait_focus) > 0) {
  psi_mean_arr <- if (length(psi_dim) >= 4) {
    apply(psi_samples, c(2, 3, 4), mean, na.rm = TRUE)
  } else {
    array(apply(psi_samples, c(2, 3), mean, na.rm = TRUE),
          dim = c(n_sp, n_sites, 1))
  }

  for (t in seq_len(n_periods)) {
    for (s in seq_len(n_sites)) {
      psi_v <- psi_mean_arr[, s, t]
      fn <- div_functional(psi_v, trait_mat)
      row <- list(grid_cell = sites[s], period = paste0("P", t))
      for (tr in trait_focus) {
        row[[paste0("cwm_", tr)]] <- fn$cwm[[tr]] %||% NA_real_
      }
      cwm_rows[[length(cwm_rows) + 1]] <- as_tibble(row)
    }
  }
  cwm_df <- bind_rows(cwm_rows)
  write_csv(cwm_df,
            v4_file("results", paste0("table_cwm_spatial_pattern_", run_label)))
}

# ── 9. Naive vs corrected richness ──────────────────────────────────
# 从 v3 survey_history 重建朴素丰富度（基于原始 y 检测频次）
survey <- safe_read(file.path(DIRS$derived, "survey_history_v4.rds"), quiet = TRUE) %||%
          safe_read(file.path(DIRS$v3_derived, "survey_history_v3.rds"), quiet = TRUE)
if (!is.null(survey)) {
  y_mat <- survey$y
  sp_ok <- intersect(species, rownames(y_mat))
  if (length(sp_ok) > 0) {
    y_use <- y_mat[sp_ok, , drop = FALSE]
    naive_rich <- matrix(NA_integer_, nrow = n_sites, ncol = n_periods)
    for (t in seq_len(n_periods)) {
      cs <- (t - 1) * n_sites + 1
      ce <- t * n_sites
      if (ce <= ncol(y_use)) {
        naive_rich[, t] <- colSums(y_use[, cs:ce, drop = FALSE] > 0, na.rm = TRUE)
      }
    }

    corrected_mean <- apply(metric_arr[, , , "corrected_richness"], c(2, 3), mean, na.rm = TRUE)

    naive_trends <- apply(naive_rich, 1, function(x)
      if (sum(!is.na(x)) >= 3) coef(lm(x ~ seq_len(n_periods)))[2] else NA_real_)
    corrected_trends <- apply(corrected_mean, 1, function(x)
      if (sum(!is.na(x)) >= 3) coef(lm(x ~ seq_len(n_periods)))[2] else NA_real_)

    direction_flipped <- sign(naive_trends) != sign(corrected_trends) &
                          !is.na(naive_trends) & !is.na(corrected_trends) &
                          naive_trends != 0 & corrected_trends != 0
    flip_rate <- mean(direction_flipped, na.rm = TRUE)
    abs_diff_med <- median(abs(corrected_trends - naive_trends), na.rm = TRUE)
    message(sprintf("[05] Naive vs Corrected: flip_rate=%.1f%%, median |diff|=%.3f",
                    flip_rate * 100, abs_diff_med))

    write_csv(tibble(
      grid_cell        = sites,
      naive_trend      = naive_trends,
      corrected_trend  = corrected_trends,
      trend_diff       = corrected_trends - naive_trends,
      direction_flipped = direction_flipped
    ), v4_file("results", paste0("table_naive_vs_corrected_", run_label)))
  }
}

# ── 10. brms 驱动回归（horseshoe + GP k 三档） ──────────────────────
if (!is.null(grid_env)) {
  driver_metrics <- c("corrected_richness", "shannon", "trait_volume", "pd_prob")
  driver_vars <- c("bio4", "bio7", "bio11", "bio13",
                   "elev_mean", "elev_sd", "texture_shannon",
                   "habitat_diversity_shannon",
                   "hfi_mean", "landcover_built", "landcover_cropland",
                   "centroid_lon", "centroid_lat")

  brms_results <- list()
  loo_results  <- list()
  r2_results   <- list()

  for (m in driver_metrics) {
    trend_m <- trend_summary_df |>
      filter(metric == m, method == "theil_sen") |>
      inner_join(grid_env, by = "grid_cell") |>
      drop_na()

    if (nrow(trend_m) < 50) {
      message(sprintf("[05] Skip brms %s: only %d rows", m, nrow(trend_m)))
      next
    }

    for (v in driver_vars) {
      if (v %in% names(trend_m)) trend_m[[paste0("z_", v)]] <- as.numeric(scale(trend_m[[v]]))
    }
    z_vars <- intersect(paste0("z_", driver_vars), names(trend_m))
    z_vars_gp <- setdiff(z_vars, c("z_centroid_lon", "z_centroid_lat"))

    # 三档 GP k
    for (k in BRMS_GP_K_VEC) {
      stem <- sprintf("brms_driver_%s_k%d_%s", m, k, run_label)
      message(sprintf("[05] brms: metric=%s, GP k=%d, n=%d", m, k, nrow(trend_m)))

      f <- as.formula(sprintf(
        "mean ~ %s + gp(centroid_lon, centroid_lat, k = %d)",
        paste(z_vars_gp, collapse = " + "), k
      ))

      # v4：horseshoe 稀疏先验（C6 修复）
      prior_v <- c(
        prior(horseshoe(df = 1, par_ratio = 0.1), class = "b"),
        prior(exponential(1), class = "sd"),
        prior(normal(0, 1), class = "Intercept")
      )

      fit <- tryCatch(
        brm(formula = f, data = trend_m, family = gaussian(),
            prior = prior_v,
            iter = BRMS_ITER, warmup = BRMS_WARMUP,
            chains = BRMS_CHAINS, cores = max(1, parallel::detectCores() - 1),
            control = list(adapt_delta = BRMS_ADAPT_DELTA,
                            max_treedepth = BRMS_MAX_TREED),
            seed = BRMS_SEED + k,
            backend = "cmdstanr",
            silent = 2, refresh = 0),
        error = function(e) {
          message(sprintf("[05] brms %s k=%d failed: %s", m, k, e$message))
          NULL
        }
      )
      if (is.null(fit)) next

      checkpoint_save(fit, stem, subdir = "derived")

      # 系数
      coef_df <- as_tibble(fixef(fit), rownames = "term") |>
        rename(estimate = Estimate, q025 = `Q2.5`, q975 = `Q97.5`)
      coef_df$metric <- m; coef_df$gp_k <- k
      brms_results[[length(brms_results) + 1]] <- coef_df

      # R² marginal vs conditional（C6 修复要求报告）
      r2 <- tryCatch(bayes_R2(fit), error = function(e) NULL)
      if (!is.null(r2)) {
        # bayes_R2 默认 conditional；marginal 需手动算
        r2_marg <- tryCatch(bayes_R2(fit, re_formula = NA), error = function(e) NULL)
        r2_results[[length(r2_results) + 1]] <- tibble(
          metric = m, gp_k = k,
          R2_conditional = r2[1, "Estimate"] %||% NA_real_,
          R2_conditional_q025 = r2[1, "Q2.5"] %||% NA_real_,
          R2_conditional_q975 = r2[1, "Q97.5"] %||% NA_real_,
          R2_marginal    = r2_marg[1, "Estimate"] %||% NA_real_,
          R2_marginal_q025 = r2_marg[1, "Q2.5"] %||% NA_real_,
          R2_marginal_q975 = r2_marg[1, "Q97.5"] %||% NA_real_
        )
      }

      # LOO（C1 修复：loo::loo 替代 lo）
      loo_obj <- tryCatch(loo::loo(fit), error = function(e) {
        message(sprintf("[05] LOO %s k=%d failed: %s", m, k, e$message))
        NULL
      })
      if (!is.null(loo_obj)) {
        loo_results[[length(loo_results) + 1]] <- tibble(
          metric = m, gp_k = k,
          elpd_loo  = loo_obj$estimates["elpd_loo", "Estimate"],
          elpd_se   = loo_obj$estimates["elpd_loo", "SE"],
          p_loo     = loo_obj$estimates["p_loo", "Estimate"],
          looic     = loo_obj$estimates["looic", "Estimate"]
        )
        checkpoint_save(loo_obj, paste0("loo_", stem), subdir = "derived")
      }

      # DHARMa
      check_dharma_gate(fit, save_dir = DIRS$figures,
                         stem = sprintf("dharma_driver_%s_k%d_%s", m, k, run_label))
    }
  }

  if (length(brms_results) > 0) {
    write_csv(bind_rows(brms_results),
              v4_file("results", paste0("table_brms_driver_coefs_", run_label)))
  }
  if (length(r2_results) > 0) {
    write_csv(bind_rows(r2_results),
              v4_file("results", paste0("table_brms_driver_R2_", run_label)))
  }
  if (length(loo_results) > 0) {
    write_csv(bind_rows(loo_results),
              v4_file("results", paste0("table_loo_driver_", run_label)))
  }
}

# ── 11. 物种层趋势 + 性状/环境关联 + 热点 ───────────────────────────
message("[05] Species-level trend analysis")
sp_trend_arr <- array(NA_real_,
                       dim = c(n_draws, n_sp, 2),
                       dimnames = list(NULL, species, c("theil_sen", "ols")))
sp_psi_mean <- array(NA_real_,
                      dim = c(n_draws, n_sp, n_periods),
                      dimnames = list(NULL, species, NULL))

if (length(psi_dim) >= 4) {
  for (d in seq_len(n_draws)) {
    for (sp in seq_len(n_sp)) {
      psi_sp <- psi_samples[d, sp, , ]
      psi_mean_t <- colMeans(psi_sp, na.rm = TRUE)
      sp_psi_mean[d, sp, ] <- psi_mean_t
      if (sum(!is.na(psi_mean_t)) >= 3) {
        sp_trend_arr[d, sp, "theil_sen"] <- theil_sen_slope(psi_mean_t)
        sp_trend_arr[d, sp, "ols"]       <- coef(lm(psi_mean_t ~ seq_len(n_periods)))[2]
      }
    }
    if (d %% 50 == 0) message(sprintf("  [sp trend] draw %d/%d", d, n_draws))
  }
}

sp_trend_rows <- list()
for (method in c("theil_sen", "ols")) {
  for (sp in seq_len(n_sp)) {
    vals <- sp_trend_arr[, sp, method]; vals <- vals[!is.na(vals)]
    if (length(vals) < 10) next
    q <- quantile(vals, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    sp_trend_rows[[length(sp_trend_rows) + 1]] <- tibble(
      species = species[sp], method = method,
      mean = mean(vals), sd = sd(vals),
      q025 = q[1], q25 = q[2], median = q[3], q75 = q[4], q975 = q[5],
      p_positive = mean(vals > 0), p_negative = mean(vals < 0)
    )
  }
}
sp_trend_df <- bind_rows(sp_trend_rows)
write_csv(sp_trend_df, v4_file("results", paste0("table_species_trend_", run_label)))

sp_classify <- sp_trend_df |>
  filter(method == "theil_sen") |>
  mutate(
    trend_class = case_when(
      p_positive > 0.95 ~ "expanding",
      p_negative > 0.95 ~ "contracting",
      TRUE              ~ "stable"
    )
  )
write_csv(sp_classify, v4_file("results", paste0("table_species_trend_classify_", run_label)))
message(sprintf("[05] Species: %d expanding, %d stable, %d contracting",
                sum(sp_classify$trend_class == "expanding"),
                sum(sp_classify$trend_class == "stable"),
                sum(sp_classify$trend_class == "contracting")))

# 热点地图（每网格的扩张/收缩物种 psi 之和）
expanding_sp  <- sp_classify$species[sp_classify$trend_class == "expanding"]
contracting_sp <- sp_classify$species[sp_classify$trend_class == "contracting"]
if (length(psi_dim) >= 4 && (length(expanding_sp) + length(contracting_sp)) > 0) {
  psi_mean_arr <- apply(psi_samples, c(2, 3, 4), mean, na.rm = TRUE)
  hotspot_rows <- list()
  for (t in seq_len(n_periods)) {
    exp_rich <- if (length(expanding_sp) > 0) {
      colSums(psi_mean_arr[species %in% expanding_sp, , t, drop = FALSE], na.rm = TRUE)
    } else rep(0, n_sites)
    con_rich <- if (length(contracting_sp) > 0) {
      colSums(psi_mean_arr[species %in% contracting_sp, , t, drop = FALSE], na.rm = TRUE)
    } else rep(0, n_sites)
    hotspot_rows[[length(hotspot_rows) + 1]] <- tibble(
      grid_cell = sites, period = paste0("P", t),
      expanding_richness = exp_rich,
      contracting_richness = con_rich,
      net_trend_index = exp_rich - con_rich
    )
  }
  write_csv(bind_rows(hotspot_rows),
            v4_file("results", paste0("table_species_hotspot_", run_label)))
}

log_time("05", "DONE")
