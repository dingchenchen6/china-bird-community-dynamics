#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   v3 的多样性函数中 eps=1e-12 让"几乎所有物种"被计入加权，
#   PSI_EPS_DEFAULT 默认值改 0.05 后多样性指标会改变多少？
#   （C12 修复对应的敏感性测试）
#
# Objective / 分析目标:
#   扫描 eps ∈ {1e-12, 1e-6, 0.01, 0.05, 0.10}：
#     - 重新计算每网格每 period 的 trait_volume / rao_q / feve / fdiv
#     - 计算各 eps 之间的 Spearman 相关
#     - 各 eps 下 5-period 趋势方向一致性
#
# Input data / 输入数据:
#   data/derived_v4/psi_samples_thinned_<run_label>_combined.qs
#   data/derived_v4/trait_extended_v4.rds
#
# Main workflow / 主要流程:
#   1. 加载 thinned psi（只用 50 draws 节省时间）
#   2. 5 档 eps，每档对每网格每 period 算 trait_volume / rao_q
#   3. 5 档之间的相关 + 趋势方向一致性
#
# Key assumptions / 关键假设:
#   utils_diversity_v4 已加载（div_functional 支持 eps 参数）
#
# Main packages / 主要包:
#   dplyr, tidyr, qs
#
# Output directory / 输出路径:
#   results_v4/table_eps_sensitivity_correlation.csv
#   results_v4/table_eps_sensitivity_trend_direction.csv
#   figures_v4/fig_v4_eps_sensitivity.{png,pdf}
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(qs)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_diversity.R"))
source(file.path(CODE_V4, "utils_mapping.R"))
P <- ensure_v4_dirs()

log_time("17", "Eps threshold sensitivity")

is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

# ── 加载 thinned psi ────────────────────────────────────────────────
thin <- safe_read(file.path(DIRS$derived,
  paste0("psi_samples_thinned_", run_label, "_combined.qs"))) %||%
        safe_read(file.path(DIRS$derived,
  paste0("psi_samples_thinned_", run_label, ".qs")))
if (is.null(thin)) stop("[17] thinned psi not found")

psi <- thin$psi_samples_thinned
species <- thin$species
sites   <- thin$sites
psi_dim <- dim(psi)
n_draws <- min(50, psi_dim[1])   # 17 只用 50 draws
draw_idx <- round(seq(1, psi_dim[1], length.out = n_draws))
n_sp <- psi_dim[2]
n_sites <- psi_dim[3]
n_periods <- if (length(psi_dim) >= 4) psi_dim[4] else 1L

# 性状矩阵
trait_ext <- safe_read(file.path(DIRS$derived, "trait_extended_v4.rds"), quiet = TRUE) %||%
             safe_read(file.path(DIRS$v3_derived, "trait_extended_v3.rds"), quiet = TRUE)
if (is.null(trait_ext)) stop("[17] trait_extended not found")

trait_aligned <- tibble(species = species) |> left_join(trait_ext, by = "species")
all_trait_vars <- intersect(TRAIT_VARS_ALL, names(trait_aligned))
trait_mat <- prepare_trait_matrix(trait_aligned, all_trait_vars)

# ── 扫描 eps ────────────────────────────────────────────────────────
eps_vec <- c(1e-12, 1e-6, 0.01, 0.05, 0.10)
metrics <- c("trait_volume", "rao_q", "feve", "fdiv")

result_arr <- array(NA_real_,
                     dim = c(length(eps_vec), n_sites, n_periods, length(metrics)),
                     dimnames = list(as.character(eps_vec), sites, NULL, metrics))

# 用 draws 的后验均值 psi 即可（不必每个 draw 都算）
psi_mean <- if (length(psi_dim) >= 4) {
  apply(psi[draw_idx, , , , drop = FALSE], c(2, 3, 4), mean)
} else {
  apply(psi[draw_idx, , , drop = FALSE], c(2, 3), mean)
}

for (ie in seq_along(eps_vec)) {
  eps <- eps_vec[ie]
  message(sprintf("[17] eps = %g", eps))
  for (t in seq_len(n_periods)) {
    for (s in seq_len(n_sites)) {
      psi_v <- if (length(psi_dim) >= 4) psi_mean[, s, t] else psi_mean[, s]
      fn <- div_functional(psi_v, trait_mat, eps = eps)
      result_arr[ie, s, t, "trait_volume"] <- fn$trait_volume
      result_arr[ie, s, t, "rao_q"]        <- fn$rao_q
      result_arr[ie, s, t, "feve"] <- div_feve(psi_v, trait_mat, eps = eps)
      result_arr[ie, s, t, "fdiv"] <- div_fdiv(psi_v, trait_mat, eps = eps)
    }
  }
}

# ── Spearman 相关：每两档 eps 之间 ──────────────────────────────────
cor_rows <- list()
for (m in metrics) {
  for (i in seq_len(length(eps_vec) - 1)) {
    for (j in (i + 1):length(eps_vec)) {
      a <- as.vector(result_arr[i, , , m])
      b <- as.vector(result_arr[j, , , m])
      ok <- !is.na(a) & !is.na(b)
      if (sum(ok) > 50) {
        rho <- cor(a[ok], b[ok], method = "spearman")
        cor_rows[[length(cor_rows) + 1]] <- tibble(
          metric = m,
          eps_i = eps_vec[i], eps_j = eps_vec[j],
          spearman_rho = rho, n = sum(ok)
        )
      }
    }
  }
}
cor_df <- bind_rows(cor_rows)
write_csv(cor_df, v4_file("results", "table_eps_sensitivity_correlation"))

# ── 趋势方向一致性 ──────────────────────────────────────────────────
trd_rows <- list()
for (m in metrics) {
  for (ie in seq_along(eps_vec)) {
    trd <- apply(result_arr[ie, , , m], 1, function(x) {
      if (sum(!is.na(x)) >= 3) coef(lm(x ~ seq_len(n_periods)))[2] else NA_real_
    })
    trd_rows[[length(trd_rows) + 1]] <- tibble(
      metric = m, eps = eps_vec[ie],
      grid_cell = sites, trend = trd, sign = sign(trd)
    )
  }
}
trd_df <- bind_rows(trd_rows)

# 各 metric 在 eps_default(0.05) 与其他 eps 的方向一致性
default_eps <- 0.05
direction_rows <- list()
for (m in metrics) {
  ref <- trd_df |> filter(metric == m, eps == default_eps) |>
    select(grid_cell, ref_sign = sign)
  for (ie in eps_vec) {
    if (ie == default_eps) next
    cmp <- trd_df |> filter(metric == m, eps == ie) |>
      select(grid_cell, this_sign = sign) |>
      inner_join(ref, by = "grid_cell") |>
      filter(!is.na(this_sign), !is.na(ref_sign), this_sign != 0, ref_sign != 0)
    direction_rows[[length(direction_rows) + 1]] <- tibble(
      metric = m, eps = ie, default_eps = default_eps,
      consistency = mean(cmp$this_sign == cmp$ref_sign, na.rm = TRUE),
      n_compared = nrow(cmp)
    )
  }
}
write_csv(bind_rows(direction_rows),
          v4_file("results", "table_eps_sensitivity_trend_direction"))

# ── 简单可视化 ──────────────────────────────────────────────────────
p_cor <- ggplot(cor_df, aes(factor(eps_i), factor(eps_j), fill = spearman_rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(spearman_rho, 2)), size = 2.5) +
  scale_fill_gradient(low = "white", high = NATURE_ACCENT,
                      name = "Spearman ρ", limits = c(0.5, 1)) +
  facet_wrap(~ metric) +
  labs(x = "eps_i", y = "eps_j",
       title = "Eps sensitivity: pairwise Spearman correlation") +
  theme_nature_pub()
save_nature(p_cor, "fig_v4_eps_sensitivity",
            width_mm = NATURE_WIDTH_L, height_mm = 100)

log_time("17", "DONE")
