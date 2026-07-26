#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   stMsPGOcc 模型对真实观测数据的拟合优度如何？
#   后验预测检验（PPC）是审稿人对贝叶斯占域模型的标准要求
#
# Objective / 分析目标:
#   - spOccupancy::ppcOcc(fit, fit.stat="freeman-tukey", group=1)
#     按 site 分组计算 Bayesian p-value，理想 0.1 < p < 0.9
#   - 同样用 chi-square 统计量做一遍 robustness
#   - fit.y vs fit.y.rep 散点图，按 primary period 分面
#
# Input data / 输入数据:
#   data/derived_v4/stMsPGOcc_fit_<run_label>_combined.qs
#
# Main workflow / 主要流程:
#   1. 加载 fit
#   2. ppcOcc 两种统计量
#   3. 输出 Bayesian p-value 表 + 散点图
#
# Key assumptions / 关键假设:
#   - spOccupancy >= 0.7（含 ppcOcc）
#   - 在 4 chain 合并后 fit.y.rep 内存足够
#
# Main packages / 主要包:
#   spOccupancy, ggplot2, dplyr, tidyr, qs
#
# Output directory / 输出路径:
#   results_v4/table_ppc_bayesian_pvalue_<run_label>.csv
#   figures_v4/fig_v4_ppc_fit_y_vs_rep_<run_label>.{png,pdf}
# ============================================================

suppressPackageStartupMessages({
  library(spOccupancy); library(readr); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(qs)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
P <- ensure_v4_dirs()

log_time("05c", "PPC Bayesian p-value")

is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

fit_path <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, "_combined.qs"))
if (!file.exists(fit_path)) {
  fit_path <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, ".qs"))
}
fit <- safe_read(fit_path)
if (is.null(fit)) stop("[05c] fit not found")

# ── PPC: Freeman-Tukey 与 Chi-square ────────────────────────────────
run_ppc <- function(fit, stat, group) {
  tryCatch({
    spOccupancy::ppcOcc(fit, fit.stat = stat, group = group)
  }, error = function(e) {
    message(sprintf("[05c] ppcOcc %s/group=%d failed: %s",
                    stat, group, e$message))
    NULL
  })
}

results <- list()
plot_data <- list()

for (stat in c("freeman-tukey", "chi-squared")) {
  for (grp in c(1)) {   # group=1 按 site，group=2 按 replicate（spOccupancy）
    ppc <- run_ppc(fit, stat, grp)
    if (is.null(ppc)) next

    # Bayesian p-value：按 species 求 P(fit.y.rep >= fit.y)
    fy <- ppc$fit.y      # draws × species × sites (3D)
    fy_rep <- ppc$fit.y.rep
    if (length(dim(fy)) >= 2) {
      # by species
      n_sp <- dim(fy)[2]
      sp_names <- rownames(fit$y) %||% paste0("sp", seq_len(n_sp))
      bp <- numeric(n_sp)
      for (i in seq_len(n_sp)) {
        if (length(dim(fy)) == 3) {
          bp[i] <- mean(fy_rep[, i, ] > fy[, i, ], na.rm = TRUE)
        } else {
          bp[i] <- mean(fy_rep[, i] > fy[, i], na.rm = TRUE)
        }
      }
      results[[length(results) + 1]] <- tibble(
        species = sp_names,
        statistic = stat,
        group = grp,
        bayesian_p = bp,
        flag = case_when(
          bp < 0.05 | bp > 0.95 ~ "EXTREME",
          bp < 0.1  | bp > 0.9  ~ "marginal",
          TRUE                   ~ "ok"
        )
      )

      # 散点图：fit.y vs fit.y.rep 的后验均值（每物种）
      if (length(dim(fy)) == 3) {
        fy_mean    <- apply(fy,    c(2, 3), mean, na.rm = TRUE)
        fyrep_mean <- apply(fy_rep, c(2, 3), mean, na.rm = TRUE)
        # flatten by species
        for (i in seq_len(n_sp)) {
          plot_data[[length(plot_data) + 1]] <- tibble(
            species   = sp_names[i],
            statistic = stat,
            fit_y     = fy_mean[i, ],
            fit_y_rep = fyrep_mean[i, ]
          )
        }
      }
    }
  }
}

if (length(results) > 0) {
  res_df <- bind_rows(results)
  write_csv(res_df, v4_file("results", paste0("table_ppc_bayesian_pvalue_", run_label)))

  # 全局摘要
  global <- res_df |>
    group_by(statistic) |>
    summarise(
      median_bp = median(bayesian_p, na.rm = TRUE),
      pct_extreme = mean(flag == "EXTREME", na.rm = TRUE) * 100,
      pct_marginal = mean(flag == "marginal", na.rm = TRUE) * 100,
      .groups = "drop"
    )
  write_csv(global, v4_file("results", paste0("table_ppc_global_summary_", run_label)))
  message(sprintf("[05c] PPC global: median bp = %s",
                  paste(global$statistic, round(global$median_bp, 3), collapse = "; ")))
}

# ── 散点图：fit.y vs fit.y.rep（前 12 物种） ─────────────────────────
if (length(plot_data) > 0) {
  pd_df <- bind_rows(plot_data)
  top_sp <- pd_df |>
    group_by(species) |>
    summarise(s = sum(fit_y, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(s)) |> slice_head(n = 12) |> pull(species)
  pd_show <- pd_df |> filter(species %in% top_sp)

  p <- ggplot(pd_show, aes(fit_y, fit_y_rep, colour = statistic)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2,
                colour = "grey50", linewidth = 0.2) +
    geom_point(alpha = 0.6, size = 0.6) +
    facet_wrap(~ species, scales = "free", ncol = 4) +
    scale_colour_brewer(palette = "Dark2", name = "Statistic") +
    labs(x = expression(T(y) ~ "(observed)"),
         y = expression(T(y[rep]) ~ "(replicated)"),
         title = "Posterior predictive check: fit.y vs fit.y.rep (top-12 species by detection)") +
    theme_minimal(base_size = NATURE_PT, base_family = NATURE_FONT) +
    theme(legend.position = "top",
          strip.text = element_text(size = 6))

  ggsave(file.path(DIRS$figures,
                    paste0("fig_v4_ppc_fit_y_vs_rep_", run_label, ".png")),
         plot = p, width = NATURE_WIDTH_L, height = 130, units = "mm",
         dpi = NATURE_DPI)
  ggsave(file.path(DIRS$figures,
                    paste0("fig_v4_ppc_fit_y_vs_rep_", run_label, ".pdf")),
         plot = p, width = NATURE_WIDTH_L, height = 130, units = "mm",
         device = cairo_pdf)
}

log_time("05c", "DONE")
