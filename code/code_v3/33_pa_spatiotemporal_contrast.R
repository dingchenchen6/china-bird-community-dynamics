#!/usr/bin/env Rscript
# ============================================================
# 33_pa_spatiotemporal_contrast.R
#
# Scientific question / 科学问题:
#   保护地内外，鸟类群落 25 年的时空演变轨迹是否分化？具体地：
#   (1) 物种占域：保护区是否维持了特化种、抑制了泛化种的相对优势？
#   (2) 群落组成：保护区内的时间 beta 与空间 beta 演变是否不同？
#   (3) 功能同质化：保护区是否延缓了功能空间的收缩与群落趋同？
#   Do community trajectories diverge inside vs outside protected areas?
#
#   与脚本 31 的分工：31 比较"端点趋势"的处理效应（ATT/DiD）；
#   本脚本比较"逐期轨迹"，并把对比推进到物种占域与空间同质化速率层面。
#
# Objective / 分析目标:
#   1) 保护梯度分层的多维多样性逐期轨迹（含后验可信区间）
#   2) period × PA 交互的正式统计检验（轨迹是否真的分化）
#   3) 物种层面：特化种 vs 泛化种在保护地内外的占域轨迹分化
#   4) 群落组成：时间 beta（turnover/nestedness）与空间 beta（同质化速率）内外对比
#   5) 事件研究：以保护区设立年份为事件时间的动态双重差分
#
# Input / 输入:
#   derived: psi_samples_thinned_<run_label>.rds   (4D: draws × species × sites × periods)
#   results: table_pa_grid_coverage_<run_label>.csv (脚本 31 产出)
#   results: table_community_metrics_with_cri_<run_label>_extended.csv
#   results: table_species_trend_traits_<run_label>_extended.csv (性状)
#   derived: grid_environment<GRID_TAG>_v3.rds
# Output / 输出:
#   results: table_pa_trajectory_<run_label>.csv          逐期 × 保护层 的多维轨迹
#   results: table_pa_interaction_tests_<run_label>.csv   period × PA 交互检验
#   results: table_pa_species_occupancy_<run_label>.csv   物种占域内外对比
#   results: table_pa_beta_contrast_<run_label>.csv       时间/空间 beta 内外对比
#   results: table_pa_event_study_<run_label>.csv         事件研究系数
#
# Key assumptions / 关键假设:
#   - 保护地"处理"以网格保护区面积占比定义，非二元；同时给出二元分层结果
#   - 空间 beta 的内外对比在匹配子样本上做，以缓解选址偏倚
#   - 事件研究要求保护区设立年份可得（字段 始建时间/years）
#   - 所有结论表述为"与保护相关的差异"，非纯因果
# Packages: sf, dplyr, readr, tidyr, lme4/lmerTest（可选 brms）
# ============================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(tidyr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

ensure_v3_dirs()
is_pilot  <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
run_ext   <- paste0(run_label, "_extended")
PA_HIGH   <- 0.30   # 高保护阈值（覆盖率≥30%，对齐昆明-蒙特利尔框架目标）
PA_LOW    <- 0.05   # 低保护阈值
log_time("33", "PA inside-outside spatiotemporal contrast")

res_path <- function(stem) paste0(v3_file("results", paste0(stem, "_", run_label)), ".csv")

# ── 0. 载入保护区覆盖率（脚本 31 产出）───────────────────────────────
cov_path <- res_path("table_pa_grid_coverage")
if (!file.exists(cov_path))
  stop("[33] 需要先运行 31_protected_area_effectiveness.R 生成 ", basename(cov_path))
cov <- read_csv(cov_path, show_col_types = FALSE)
frac_col <- intersect(c("pa_fraction", "pa_frac", "pa_cover"), names(cov))[1]
if (is.na(frac_col)) stop("[33] 覆盖率表中未找到保护区面积占比列")
cov <- cov |> rename(pa_frac = all_of(frac_col))
cov$pa_stratum <- cut(cov$pa_frac, breaks = c(-Inf, PA_LOW, PA_HIGH, Inf),
                      labels = c("unprotected", "partial", "protected"))
message(sprintf("[33] 保护分层: %s",
        paste(names(table(cov$pa_stratum)), table(cov$pa_stratum), sep = "=", collapse = ", ")))

# ── 1. 多维多样性的逐期轨迹 × 保护分层 ───────────────────────────────
cm_path <- paste0(v3_file("results", paste0("table_community_metrics_with_cri_", run_ext)), ".csv")
if (!file.exists(cm_path)) cm_path <- res_path("table_community_metrics_with_cri")
cm <- read_csv(cm_path, show_col_types = FALSE)

traj <- cm |>
  inner_join(cov |> select(grid_cell, pa_frac, pa_stratum), by = "grid_cell") |>
  group_by(metric, period, pa_stratum) |>
  summarise(n_grid = n(),
            mean = mean(value_mean, na.rm = TRUE),
            se   = sd(value_mean, na.rm = TRUE) / sqrt(n()),
            lwr  = mean(value_l95, na.rm = TRUE),
            upr  = mean(value_u95, na.rm = TRUE), .groups = "drop")
write_csv(traj, res_path("table_pa_trajectory"))
message("[33] 轨迹表已输出：", basename(res_path("table_pa_trajectory")))

# 轨迹分化的直观量：末期与首期的组间差值变化
div <- traj |>
  filter(period %in% c("P1", "P5")) |>
  select(metric, period, pa_stratum, mean) |>
  pivot_wider(names_from = period, values_from = mean) |>
  mutate(change = P5 - P1) |>
  select(metric, pa_stratum, change) |>
  pivot_wider(names_from = pa_stratum, values_from = change) |>
  mutate(protected_minus_unprotected = protected - unprotected)
print(as.data.frame(div), digits = 4)
message("[33] >>> protected_minus_unprotected > 0 的功能指标 = 保护区延缓了该维度的衰退")

# ── 2. period × PA 交互检验（轨迹是否统计上分化）─────────────────────
have_lmer <- requireNamespace("lme4", quietly = TRUE)
if (have_lmer) {
  suppressPackageStartupMessages(library(lme4))
  focal <- intersect(unique(cm$metric),
    c("corrected_richness", "shannon", "trait_volume", "rao_q", "fric_prob", "fdis_prob"))
  int_res <- do.call(rbind, lapply(focal, function(m) {
    d <- cm |> filter(metric == m) |>
      inner_join(cov |> select(grid_cell, pa_frac), by = "grid_cell") |>
      mutate(period_num = as.integer(sub("P", "", period)),
             z = as.numeric(scale(value_mean)))
    fit <- try(lmer(z ~ period_num * pa_frac + (1 | grid_cell), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    cf <- summary(fit)$coefficients
    r <- cf["period_num:pa_frac", ]
    data.frame(metric = m, interaction_estimate = r[1], se = r[2], t = r[3],
               n_obs = nrow(d), stringsAsFactors = FALSE)
  }))
  int_res$p_approx <- 2 * pnorm(-abs(int_res$t))
  print(int_res, digits = 3)
  write_csv(int_res, res_path("table_pa_interaction_tests"))
  message("[33] >>> 交互项显著为正 = 保护区内该指标随时间的表现优于区外。")
} else message("[33] lme4 未安装，跳过交互检验（install.packages('lme4')）")

# ── 3. 物种占域：特化种 vs 泛化种在保护地内外的分化 ───────────────────
psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))
if (is.null(psi_obj)) stop("[33] psi_samples_thinned 未找到，请先跑 05。")
psi <- psi_obj$psi_samples_thinned
if (length(dim(psi)) < 4)
  stop(sprintf("[33] psi 为 %dD，需 4D（draws×species×sites×periods）。", length(dim(psi))))
species <- psi_obj$species; sites <- psi_obj$sites
psi_mean <- apply(psi, c(2, 3, 4), mean)        # species × sites × periods
np <- dim(psi_mean)[3]

# 网格是否属于高保护层
idx_cov <- match(sites, cov$grid_cell)
in_pa  <- !is.na(idx_cov) & cov$pa_frac[idx_cov] >= PA_HIGH
out_pa <- !is.na(idx_cov) & cov$pa_frac[idx_cov] <  PA_LOW
message(sprintf("[33] 高保护格点 %d 个；低/无保护格点 %d 个", sum(in_pa), sum(out_pa)))

sp_occ <- do.call(rbind, lapply(seq_len(np), function(t) {
  data.frame(period = paste0("P", t), species = species,
             mean_psi_inside  = rowMeans(psi_mean[, in_pa,  t, drop = FALSE]),
             mean_psi_outside = rowMeans(psi_mean[, out_pa, t, drop = FALSE]),
             stringsAsFactors = FALSE)
}))
sp_occ$psi_diff <- sp_occ$mean_psi_inside - sp_occ$mean_psi_outside

# 关联性状：栖息地宽度（低=特化种）
tr_path <- paste0(v3_file("results", paste0("table_species_trend_traits_", run_ext)), ".csv")
if (!file.exists(tr_path)) tr_path <- res_path("table_species_trend_traits")
if (file.exists(tr_path)) {
  tr <- read_csv(tr_path, show_col_types = FALSE)
  hb <- intersect(c("habitat_breadth", "Habitat.Breadth"), names(tr))[1]
  if (!is.na(hb)) {
    tr2 <- tr |> select(species, habitat_breadth = all_of(hb)) |>
      mutate(sp_group = ifelse(habitat_breadth <= median(habitat_breadth, na.rm = TRUE),
                               "specialist", "generalist"))
    sp_occ <- sp_occ |> left_join(tr2, by = "species")
    grp <- sp_occ |> filter(!is.na(sp_group)) |>
      group_by(sp_group, period) |>
      summarise(mean_inside = mean(mean_psi_inside, na.rm = TRUE),
                mean_outside = mean(mean_psi_outside, na.rm = TRUE),
                mean_diff = mean(psi_diff, na.rm = TRUE), .groups = "drop")
    print(as.data.frame(grp), digits = 3)
    message("[33] >>> 若特化种的 mean_diff 随时间扩大而泛化种不变，")
    message("[33]     即为保护区选择性维持特化种的直接证据。")
  }
}
write_csv(sp_occ, res_path("table_pa_species_occupancy"))

# ── 4. 群落组成：时间 beta 与空间 beta 的内外对比 ─────────────────────
# 空间 beta = 同一期内格点两两 Sørensen 相异度的均值；其随时间下降即同质化。
spatial_beta <- function(mat) {                    # mat: species × sites（占域概率）
  if (ncol(mat) < 2) return(NA_real_)
  a <- crossprod(mat)                              # 期望共享占域（sites × sites）
  b <- colSums(mat)
  den <- outer(b, b, "+")
  d <- 1 - 2 * a / den
  d[den == 0] <- NA_real_
  mean(d[upper.tri(d)], na.rm = TRUE)
}
# 为控制格点数差异带来的偏倚，对两组做等样本量重抽样
n_sub <- min(sum(in_pa), sum(out_pa), 120)
set.seed(if (exists("SEED")) SEED else 20260728L)
beta_tbl <- do.call(rbind, lapply(seq_len(np), function(t) {
  reps <- t(sapply(1:30, function(r) {
    i_in  <- sample(which(in_pa),  n_sub)
    i_out <- sample(which(out_pa), n_sub)
    c(spatial_beta(psi_mean[, i_in,  t]), spatial_beta(psi_mean[, i_out, t]))
  }))
  data.frame(period = paste0("P", t), n_grid_each = n_sub,
             spatial_beta_inside  = mean(reps[, 1], na.rm = TRUE),
             spatial_beta_inside_sd = sd(reps[, 1], na.rm = TRUE),
             spatial_beta_outside = mean(reps[, 2], na.rm = TRUE),
             spatial_beta_outside_sd = sd(reps[, 2], na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
beta_tbl$beta_gap <- beta_tbl$spatial_beta_inside - beta_tbl$spatial_beta_outside
print(beta_tbl, digits = 4)
write_csv(beta_tbl, res_path("table_pa_beta_contrast"))
message("[33] >>> 若区外空间 beta 下降更快（同质化更快），保护区即起到了抑制趋同的作用。")

# ── 5. 事件研究：以保护区设立年份为事件时间的动态 DiD ─────────────────
yr_col <- intersect(c("pa_year_min", "establish_year", "years"), names(cov))[1]
if (!is.na(yr_col) && requireNamespace("lme4", quietly = TRUE)) {
  ev <- cm |>
    filter(metric %in% c("corrected_richness", "trait_volume", "rao_q")) |>
    inner_join(cov |> select(grid_cell, pa_frac, est_year = all_of(yr_col)), by = "grid_cell") |>
    mutate(period_num = as.integer(sub("P", "", period)),
           period_mid = 2002 + (period_num - 1) * 5,
           treated = pa_frac >= PA_HIGH & !is.na(est_year),
           event_time = ifelse(treated, floor((period_mid - est_year) / 5), NA_integer_)) |>
    filter(is.na(event_time) | abs(event_time) <= 3)
  es <- do.call(rbind, lapply(unique(ev$metric), function(m) {
    d <- ev |> filter(metric == m) |> mutate(z = as.numeric(scale(value_mean)))
    d$et <- factor(ifelse(is.na(d$event_time), "control", paste0("t", d$event_time)),
                   levels = c("control", paste0("t", sort(unique(na.omit(d$event_time))))))
    fit <- try(lme4::lmer(z ~ et + (1 | grid_cell), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    cf <- summary(fit)$coefficients
    data.frame(metric = m, term = rownames(cf), estimate = cf[, 1],
               se = cf[, 2], t = cf[, 3], stringsAsFactors = FALSE)
  }))
  if (!is.null(es)) {
    write_csv(es, res_path("table_pa_event_study"))
    message("[33] 事件研究系数已输出（事件前系数应接近 0 = 平行趋势假设成立）")
  }
} else message("[33] 无保护区设立年份列或缺 lme4，跳过事件研究。")

log_time("33", "DONE")
