#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题: 观鸟努力增长最快的格点，是否正是校正物种数增加最多的格点？
#   Do the cells with the fastest effort growth gain the most corrected species?
# Objective / 分析目标: 为讨论中"努力增长与增益在空间上重合"一句提供量化依据或将其删除。
#   逐格计算努力增长斜率（log 记录数）与校正物种数趋势，求 Spearman 相关，并用空间分块 bootstrap 估计不确定性。
# Input / 输入: survey_history_v3.rds（det_covs$log_events）；table_community_metrics_with_cri_<RUN_LABEL>_extended.csv
# Output / 输出: results_v3/table_effort_gain_test_<RUN_LABEL>.csv
# Key assumptions / 关键假设: 空间分块约 500 km（4.5°）；用分块 bootstrap 给出相关系数的 95% 区间。
# ============================================================
CODE <- Sys.getenv("V3_CODE_DIR", "code_v4_refit")
source(file.path(CODE, "00_config.R")); source(file.path(CODE, "utils_paths.R")); set.seed(20261020)
s <- readRDS(v3_file("derived", "survey_history_v3", "rds"))
eff <- apply(s$det_covs$log_events, c(1, 2), mean, na.rm = TRUE)               # 格点 × 时期
slope <- function(y) { x <- seq_along(y); k <- is.finite(y); if (sum(k) < 3) NA else coef(lm(y[k] ~ x[k]))[2] }
eff_sl <- apply(eff, 1, slope)
cm <- read.csv(paste0(v3_file("results", paste0("table_community_metrics_with_cri_", RUN_LABEL, "_extended")), ""))
if (!nrow(cm)) stop("[S2_04] 找不到群落指标表")
cm <- cm[cm$metric == "corrected_richness", ]
cm$p <- as.integer(sub("P", "", cm$period))
rich_sl <- sapply(split(cm, cm$grid_cell), function(d) slope(d$value_mean[order(d$p)]))
env <- readRDS(v3_file("derived", "grid_environment_v3", "rds"))
d <- data.frame(grid_cell = s$sites, effort_slope = eff_sl)
d$richness_slope <- rich_sl[as.character(d$grid_cell)]
d <- merge(d, env[, c("grid_cell", "centroid_lon", "centroid_lat")], by = "grid_cell")
d <- d[is.finite(d$effort_slope) & is.finite(d$richness_slope), ]
rho <- cor(d$effort_slope, d$richness_slope, method = "spearman")
blk <- interaction(floor(d$centroid_lon / 4.5), floor(d$centroid_lat / 4.5), drop = TRUE)   # 约 500 km 分块
# 空间分块 bootstrap：按块有放回重抽，保留块内空间自相关 / block bootstrap keeps within-block autocorrelation
ids <- split(seq_len(nrow(d)), blk)
boot <- replicate(999, { i <- unlist(ids[sample(length(ids), replace = TRUE)])
  cor(d$effort_slope[i], d$richness_slope[i], method = "spearman") })
ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
out <- data.frame(n_cells = nrow(d), n_blocks = length(ids), spearman_rho = round(rho, 3),
                  block_boot_lwr = round(ci[1], 3), block_boot_upr = round(ci[2], 3),
                  note = "区间全部 > 0 才可保留讨论中的空间重合表述；否则删除该句")
write.csv(out, v3_file("results", paste0("table_effort_gain_test_", RUN_LABEL)), row.names = FALSE); print(out)
