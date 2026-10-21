#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题: 占域模型能否复现观测到的检测数据？拟合不足是否集中在某些时期或努力水平？
#   Does the model reproduce the observed detections, and does misfit concentrate by period or effort?
# Objective / 分析目标: 用 spOccupancy::ppcOcc（Freeman–Tukey）计算后验预测检验的 Bayesian p 值，
#   按物种 × 时期汇总，并按调查努力五分位汇总格点层面的拟合偏差。（BLOCKERS Gate 2；SI S4）
# Input / 输入: <OUT>/data/derived_v3/stMsPGOcc_fit_<RUN_LABEL>.rds；survey_history_v3.rds
# Output / 输出: results_v3/table_ppc_by_species_period_<RUN_LABEL>.csv、table_ppc_by_effort_<RUN_LABEL>.csv
# Key assumptions / 关键假设: Bayesian p 在 0.1–0.9 视为无明显拟合不足；先确认所装 spOccupancy 版本
#   的 ppcOcc 支持 stMsPGOcc（≥ 0.7）。
# Packages / 主要包: spOccupancy
# ============================================================
suppressPackageStartupMessages(library(spOccupancy))
source(file.path(Sys.getenv("V3_CODE_DIR", "code_v4_refit"), "00_config.R"))
source(file.path(Sys.getenv("V3_CODE_DIR", "code_v4_refit"), "utils_paths.R"))
fit_path <- v3_file("derived", paste0("stMsPGOcc_fit_", RUN_LABEL), "rds")
if (!file.exists(fit_path)) stop("[05c] 找不到拟合对象：", fit_path)
fit <- readRDS(fit_path)
set.seed(20261020)

# group = 1：按格点汇总，得到每个后验抽样 × 物种 × 格点 × 时期的偏差 / site-level discrepancy
ppc <- tryCatch(ppcOcc(fit, fit.stat = "freeman-tukey", group = 1),
                error = function(e) stop("[05c] ppcOcc 不支持该模型或版本：", conditionMessage(e)))
fy <- ppc$fit.y; fr <- ppc$fit.y.rep            # samples × species × periods（多季节多物种）
bp <- apply(fr > fy, c(2, 3), mean)
sp <- if (!is.null(fit$sp.names)) fit$sp.names else seq_len(dim(bp)[1])
tab <- data.frame(species = rep(sp, times = ncol(bp)), period = rep(seq_len(ncol(bp)), each = nrow(bp)),
                  bayes_p = as.vector(bp))
tab$flag <- tab$bayes_p < 0.1 | tab$bayes_p > 0.9
write.csv(tab, v3_file("results", paste0("table_ppc_by_species_period_", RUN_LABEL)), row.names = FALSE)
message(sprintf("[05c] 物种 × 时期 Bayesian p 落在 0.1–0.9 之外的比例：%.1f%%", 100 * mean(tab$flag)))

# 按努力五分位汇总格点偏差（若 ppcOcc 返回了格点层分位数）/ by effort quintile, if site-level output exists
gq <- ppc$fit.y.group.quants; rq <- ppc$fit.y.rep.group.quants
if (!is.null(gq) && !is.null(rq)) {
  s <- readRDS(v3_file("derived", "survey_history_v3", "rds"))
  eff <- apply(s$det_covs$log_events, c(1, 2), mean, na.rm = TRUE)            # 格点 × 时期 / site × period
  med_o <- gq[3, , , , drop = TRUE]; med_r <- rq[3, , , , drop = TRUE]          # 中位数：物种 × 格点 × 时期
  excess <- apply(med_o - med_r, c(2, 3), mean, na.rm = TRUE)                   # 格点 × 时期的平均超额偏差
  q <- cut(eff, quantile(eff, 0:5 / 5, na.rm = TRUE), include.lowest = TRUE, labels = paste0("Q", 1:5))
  out <- aggregate(as.vector(excess), list(effort_quintile = as.vector(q), period = rep(seq_len(ncol(eff)), each = nrow(eff))),
                   function(x) mean(x, na.rm = TRUE))
  names(out)[3] <- "mean_obs_minus_rep_discrepancy"
  write.csv(out, v3_file("results", paste0("table_ppc_by_effort_", RUN_LABEL)), row.names = FALSE)
} else warning("[05c] ppcOcc 未返回格点层分位数，未生成按努力汇总的表")
