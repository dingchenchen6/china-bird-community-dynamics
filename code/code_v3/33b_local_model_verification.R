#!/usr/bin/env Rscript
# ============================================================
# 33b_local_model_verification.R
#
# Scientific question / 科学问题:
#   首轮"保护年限剂量-反应显著为正"与"period×保护交互极显著"是真实效应，
#   还是随机截距模型造成的伪重复假显著？本脚本在本地用结果表直接复核。
#   Are the first-round protection effects real, or artefacts of a
#   random-intercept-only specification (pseudoreplication)?
#
# Objective / 分析目标:
#   把【旧设定：随机截距】与【新设定：随机斜率】并排拟合，并与三个
#   不依赖模型的量对照：逐格斜率与保护量的相关、原始分组均值、分层轨迹。
#   四者一致才可下结论。
#
# Input / 输入（均为已产出的结果表，不需要 psi 后验数组）:
#   results/results_v3/table_community_metrics_with_cri_<run>_extended.csv
#   results/results_v3/table_pa_grid_coverage_<run>.csv
# Output / 输出:
#   results/results_v3/table_model_specification_check_<run>.csv
#
# Key assumption / 关键假设:
#   保护覆盖率与保护年限都是网格层常量，因此"随时间的斜率如何随保护变化"
#   属于网格间比较，有效样本量是网格数而非网格×时期数。
# Packages: lme4, dplyr, readr
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(lme4)
})

REPO <- Sys.getenv("REPO_DIR",
  file.path("~", "Documents", "New project", "china-bird-community-dynamics"))
RES  <- file.path(path.expand(REPO), "results", "results_v3")
RUN  <- "v3_full_200sp_ar1_spatial"
PA_LOW <- 0.05; PA_HIGH <- 0.30

cm_p  <- file.path(RES, sprintf("table_community_metrics_with_cri_%s_extended.csv", RUN))
cov_p <- file.path(RES, sprintf("table_pa_grid_coverage_%s.csv", RUN))
stopifnot(file.exists(cm_p), file.exists(cov_p))

cm  <- read_csv(cm_p,  show_col_types = FALSE)
cov <- read_csv(cov_p, show_col_types = FALSE) |>
  mutate(pa_frac = if ("pa_frac" %in% names(cur_data_all())) pa_frac else pa_frac_all,
         protect_age = ifelse(is.na(pa_year_first_major), NA_real_,
                              2024 - pa_year_first_major),
         pa_stratum = cut(pa_frac, c(-Inf, PA_LOW, PA_HIGH, Inf),
                          labels = c("unprotected", "partial", "protected")))

FOCAL <- c("corrected_richness", "trait_volume", "rao_q", "shannon")
cat(sprintf("\n网格数 %d；高保护 %d，部分 %d，未保护 %d\n", nrow(cov),
            sum(cov$pa_stratum == "protected",  na.rm = TRUE),
            sum(cov$pa_stratum == "partial",    na.rm = TRUE),
            sum(cov$pa_stratum == "unprotected", na.rm = TRUE)))

# ── 工具：拟合两种设定并抽取交互项 ───────────────────────────────────
fit_pair <- function(d, pred) {
  f_int <- as.formula(sprintf("z ~ period_num * %s + (1 | grid_cell)", pred))
  f_slp <- as.formula(sprintf("z ~ period_num * %s + (period_num | grid_cell)", pred))
  grab <- function(fit, tag) {
    if (inherits(fit, "try-error")) return(NULL)
    cf <- summary(fit)$coefficients
    r  <- cf[paste0("period_num:", pred), ]
    data.frame(spec = tag, estimate = r[1], se = r[2], t = r[3],
               p_approx = 2 * pnorm(-abs(r[3])), stringsAsFactors = FALSE)
  }
  rbind(grab(try(lmer(f_int, data = d), silent = TRUE), "旧: 随机截距"),
        grab(try(lmer(f_slp, data = d), silent = TRUE), "新: 随机斜率"))
}

# 逐格线性斜率（不依赖任何混合模型）
per_grid_slope <- function(d) {
  d |> group_by(grid_cell) |>
    summarise(sl = if (sum(!is.na(value_mean)) >= 3)
                     unname(coef(lm(value_mean ~ period_num))[2]) else NA_real_,
              .groups = "drop") |>
    filter(!is.na(sl))
}

out <- list()

# ── A. period × 保护覆盖率 交互 ──────────────────────────────────────
cat("\n══════ A. period × 保护覆盖率（全部网格）══════\n")
for (m in FOCAL) {
  d <- cm |> filter(metric == m) |>
    inner_join(cov |> select(grid_cell, pa_frac), by = "grid_cell") |>
    mutate(period_num = as.integer(sub("P", "", period)),
           z = as.numeric(scale(value_mean)))
  if (!nrow(d)) next
  res <- fit_pair(d, "pa_frac")
  if (is.null(res)) next
  sl <- per_grid_slope(d) |> inner_join(cov |> select(grid_cell, pa_frac), by = "grid_cell")
  r  <- suppressWarnings(cor(sl$pa_frac, sl$sl, use = "complete.obs"))
  res$analysis <- "period_x_pa_frac"; res$metric <- m
  res$effect_size_r <- r; res$n_grid <- nrow(sl); res$n_obs <- nrow(d)
  out[[length(out) + 1]] <- res
  cat(sprintf("\n%-20s  逐格斜率 vs 保护覆盖率 r = %+.3f  (n_grid=%d, n_obs=%d)\n",
              m, r, nrow(sl), nrow(d)))
  print(res[, c("spec", "estimate", "se", "t", "p_approx")], row.names = FALSE, digits = 3)
}

# ── B. period × 保护年限 剂量-反应（仅已保护网格）────────────────────
cat("\n══════ B. period × 保护年限（仅有设立年份的网格）══════\n")
for (m in FOCAL) {
  d <- cm |> filter(metric == m) |>
    inner_join(cov |> select(grid_cell, protect_age), by = "grid_cell") |>
    filter(!is.na(protect_age), protect_age > 0) |>
    mutate(period_num = as.integer(sub("P", "", period)),
           z = as.numeric(scale(value_mean)),
           age_s = as.numeric(scale(protect_age)))
  if (nrow(d) < 50) next
  res <- fit_pair(d, "age_s")
  if (is.null(res)) next
  sl <- per_grid_slope(d) |> inner_join(cov |> select(grid_cell, protect_age), by = "grid_cell")
  r  <- suppressWarnings(cor(sl$protect_age, sl$sl, use = "complete.obs"))
  res$analysis <- "period_x_protect_age"; res$metric <- m
  res$effect_size_r <- r; res$n_grid <- nrow(sl); res$n_obs <- nrow(d)
  out[[length(out) + 1]] <- res
  cat(sprintf("\n%-20s  逐格斜率 vs 保护年限 r = %+.3f  (n_grid=%d, n_obs=%d)\n",
              m, r, nrow(sl), nrow(d)))
  print(res[, c("spec", "estimate", "se", "t", "p_approx")], row.names = FALSE, digits = 3)
}

# ── C. 原始分层轨迹（完全不依赖模型的对照）──────────────────────────
cat("\n══════ C. 原始分层：P1→P5 变化（不依赖模型）══════\n")
raw <- cm |> filter(metric %in% FOCAL, period %in% c("P1", "P5")) |>
  inner_join(cov |> select(grid_cell, pa_stratum), by = "grid_cell") |>
  filter(!is.na(pa_stratum)) |>
  group_by(metric, pa_stratum, period) |>
  summarise(mean = mean(value_mean, na.rm = TRUE), n = n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = period, values_from = c(mean, n)) |>
  mutate(change = mean_P5 - mean_P1, pct = 100 * change / mean_P1)
print(as.data.frame(raw |> select(metric, pa_stratum, n = n_P1, mean_P1, mean_P5, pct)),
      row.names = FALSE, digits = 4)

# ── 汇总落盘 ─────────────────────────────────────────────────────────
if (length(out)) {
  final <- bind_rows(out) |>
    select(analysis, metric, spec, estimate, se, t, p_approx, effect_size_r, n_grid, n_obs)
  write_csv(final, file.path(RES, sprintf("table_model_specification_check_%s.csv", RUN)))
  cat("\n已写出: table_model_specification_check_", RUN, ".csv\n", sep = "")

  cat("\n══════ 判读 ══════\n")
  chg <- final |> filter(analysis == "period_x_protect_age") |>
    group_by(metric) |>
    summarise(p_old = p_approx[spec == "旧: 随机截距"][1],
              p_new = p_approx[spec == "新: 随机斜率"][1],
              r = effect_size_r[1], .groups = "drop")
  for (i in seq_len(nrow(chg))) {
    w <- chg[i, ]
    verdict <- if (!is.na(w$r) && abs(w$r) < 0.1) "效应量近零 → 无实质剂量-反应"
               else if (!is.na(w$p_new) && w$p_new < 0.05) "随机斜率下仍显著 → 结论稳健"
               else "随机斜率下不再显著 → 原结论为伪重复假象"
    cat(sprintf("  %-20s p旧=%.4g  p新=%.4g  r=%+.3f  → %s\n",
                w$metric, w$p_old, w$p_new, w$r, verdict))
  }
}
