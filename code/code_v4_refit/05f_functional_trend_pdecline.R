#!/usr/bin/env Rscript
# ============================================================
# 05f_functional_trend_pdecline.R
#
# Scientific question / 科学问题:
#   校正后功能多样性(功能体积 trait_volume、Rao's Q)是否在全国尺度上
#   随时间显著下降? —— 为 GCB 稿件 "functional homogenization" 主张提供
#   斜率后验概率 P(slope<0),闭合当前 [PENDING: P(decline)] 占位。
#   Do functional trait volume and Rao's Q decline at the national scale?
#
# Objective / 分析目标:
#   从逐 draw 的多样性立方(metric_arrays_checkpoint)出发,对每个后验 draw
#   计算"全国均值时间序列"的 Theil-Sen 斜率,得到斜率的后验分布,
#   报告 posterior mean / median / 95% CrI / P(decline) / P(increase)。
#
# Input / 输入:
#   data/derived_v3/metric_arrays_checkpoint_<run_label>_extended.rds
#     (ck$metric_arr: [n_draws, n_sites, n_periods, n_metrics])
# Output / 输出:
#   results_v3/table_functional_trend_pdecline_<run_label>.csv + 控制台打印
#
# Key assumption / 关键假设:
#   "全国功能同质化"以全国 grid 均值的时间斜率为准(与稿件全国均值口径一致);
#   Theil-Sen 对 5 期短序列稳健。
# Packages: parallel, readr; utils_diversity::theil_sen_slope
# Reuse: 完全复用 05e 的 checkpoint 与工具函数,不需重跑模型。
# ============================================================

suppressPackageStartupMessages({
  library(parallel)
  library(readr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_diversity.R"))  # theil_sen_slope

ensure_v3_dirs()

is_pilot      <- Sys.getenv("V3_PILOT", "0") == "1"
run_label     <- if (is_pilot) PILOT_LABEL else RUN_LABEL
run_label_ext <- paste0(run_label, "_extended")
log_time("05f", "Functional trend P(decline) from metric checkpoint")

# ── 1. 加载 metric checkpoint（与 05e 同源）───────────────────────────
ck_path <- v3_file("derived", paste0("metric_arrays_checkpoint_", run_label_ext), "rds")
if (!file.exists(ck_path)) {
  alt <- v3_file("derived", paste0("metric_arrays_", run_label_ext), "rds")
  if (file.exists(alt)) ck_path <- alt else
    stop(sprintf("[05f] metric checkpoint not found: %s", ck_path))
}
message(sprintf("[05f] loading: %s", ck_path))
ck <- readRDS(ck_path)

metric_arr <- ck$metric_arr            # [draws, sites, periods, metrics]
metrics    <- ck$metrics
n_draws    <- dim(metric_arr)[1]
n_sites    <- dim(metric_arr)[2]
n_periods  <- dim(metric_arr)[3]
message(sprintf("[05f] metric_arr: %d draws x %d sites x %d periods x %d metrics",
                n_draws, n_sites, n_periods, length(metrics)))

# ── 2. 目标指标：功能同质化核心 + 对照（若存在则纳入）─────────────────
target <- c("trait_volume", "rao_q", "feve_fund", "fdiv_fund",
            "fric_prob", "fdis_prob", "corrected_richness", "shannon")
target <- intersect(target, metrics)
if (!length(target)) stop("[05f] none of the target metrics present in checkpoint")
message(sprintf("[05f] metrics analysed: %s", paste(target, collapse = ", ")))

mc_cores <- min(if (exists("MC_CORES")) MC_CORES else parallel::detectCores() - 2L, 16L)
mc_cores <- max(1L, mc_cores)

# ── 3. 逐指标：每个 draw 的全国均值时间序列 → Theil-Sen 斜率后验 ───────
summ <- lapply(target, function(m) {
  arr_m <- metric_arr[, , , m]         # [draws, sites, periods]
  slopes <- unlist(mclapply(seq_len(n_draws), function(d) {
    # 该 draw 下各期的全国 grid 均值（对 sites 求均值）
    mat  <- matrix(arr_m[d, , ], nrow = n_sites, ncol = n_periods)
    natl <- colMeans(mat, na.rm = TRUE)
    if (sum(is.finite(natl)) >= 3) theil_sen_slope(natl) else NA_real_
  }, mc.cores = mc_cores))

  data.frame(
    metric       = m,
    slope_mean   = mean(slopes, na.rm = TRUE),
    slope_median = median(slopes, na.rm = TRUE),
    slope_l95    = unname(quantile(slopes, 0.025, na.rm = TRUE)),
    slope_u95    = unname(quantile(slopes, 0.975, na.rm = TRUE)),
    P_decline    = mean(slopes < 0, na.rm = TRUE),
    P_increase   = mean(slopes > 0, na.rm = TRUE),
    n_draws_used = sum(is.finite(slopes)),
    stringsAsFactors = FALSE)
})
res <- do.call(rbind, summ)

message("\n[05f] ===== Functional / diversity trend posterior (national mean, per period) =====")
print(res, row.names = FALSE, digits = 4)

# ── 4. 落盘（供回填稿件 [PENDING] 处）────────────────────────────────
out_path <- v3_file("results", paste0("table_functional_trend_pdecline_", run_label))
write_csv(res, if (grepl("\\.csv$", out_path)) out_path else paste0(out_path, ".csv"))
message(sprintf("[05f] saved: %s", out_path))

# ── 5. 直接给出可粘贴进稿件的措辞（trait_volume / rao_q）─────────────
for (m in intersect(c("trait_volume", "rao_q"), res$metric)) {
  r <- res[res$metric == m, ]
  message(sprintf(
    "[05f] %s: slope = %.4f/period (95%% CrI %.4f to %.4f), P(decline) = %.3f",
    m, r$slope_mean, r$slope_l95, r$slope_u95, r$P_decline))
}
log_time("05f", "DONE")
