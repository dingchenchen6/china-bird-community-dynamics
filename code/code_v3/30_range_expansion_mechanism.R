#!/usr/bin/env Rscript
# ============================================================
# 30_range_expansion_mechanism.R
#
# Scientific question / 科学问题:
#   校正后每格物种数为何上升？其数学等价物是"平均物种分布面积扩张"，
#   但扩张可能来自两种截然不同的机制：
#     (A) 边缘扩张 range-edge expansion —— 占据以前没有的新区域（气候驱动位移的特征）
#     (B) 范围内填充 range infilling   —— 原分布区内占域概率上升（种群恢复的特征）
#   本脚本把 ΔAOO 分解为这两个分量，并进一步检验气候位移与类群归因。
#   Why is per-grid richness rising: range-edge expansion or in-range infilling?
#
# Objective / 分析目标:
#   1) 范围扩张 vs 填充分解（全国与逐物种）
#   2) 定殖-灭绝分解（期望定殖数 / 灭绝数）
#   3) 分布重心位移：纬度与海拔（气候驱动的方向性检验）
#   4) 群落温度指数 CTI 与热适应性变化（Devictor 等的气候债框架）
#   5) 生境类群归因：谁贡献了每格 +21.6 种的增量
#
# Input / 输入:
#   derived: psi_samples_thinned_<run_label>.rds  (4D: draws × species × sites × periods)
#   derived: grid_environment<GRID_TAG>_v3.rds    (lat/lon/elevation/bio1 等)
#   results: table_traits_extended_v3.csv 或 trait 表（habitat / trophic 类群）
# Output / 输出:
#   results: table_range_decomposition_<run_label>.csv      逐物种扩张/填充分量
#   results: table_range_decomposition_summary_<run_label>.csv  全国汇总
#   results: table_colonization_extinction_<run_label>.csv  期间定殖/灭绝
#   results: table_centroid_shift_<run_label>.csv           重心纬度/海拔位移
#   results: table_cti_<run_label>.csv                      格点×期 CTI
#   results: table_guild_contribution_<run_label>.csv       类群对增量的贡献
#
# Key assumptions / 关键假设:
#   - "核心分布区"以 P1 的 psi > CORE_THRESH 定义（默认 0.5），敏感性见 CORE_THRESH_ALT
#   - 定殖/灭绝按格点独立近似由相邻期 psi 计算期望值
#   - STI 用 P1 占域加权的年均温，避免用全期数据造成循环论证
# Packages: matrixStats, readr, dplyr
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(matrixStats)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

ensure_v3_dirs()
is_pilot   <- Sys.getenv("V3_PILOT", "0") == "1"
run_label  <- if (is_pilot) PILOT_LABEL else RUN_LABEL
CORE_THRESH     <- 0.5    # 核心分布区阈值 / core-range threshold
CORE_THRESH_ALT <- c(0.3, 0.7)   # 敏感性 / sensitivity
log_time("30", "Range-expansion mechanism decomposition")

# ── 1. 载入后验占域数组 ──────────────────────────────────────────────
psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))
if (is.null(psi_obj)) stop("[30] psi_samples_thinned not found; run 05 first.")
psi <- psi_obj$psi_samples_thinned
if (length(dim(psi)) < 4)
  stop(sprintf("[30] psi is %dD; need 4D (draws x species x sites x periods). Re-thin via 05.",
               length(dim(psi))))
species <- psi_obj$species; sites <- psi_obj$sites
nd <- dim(psi)[1]; ns <- dim(psi)[2]; ng <- dim(psi)[3]; np <- dim(psi)[4]
message(sprintf("[30] psi: %d draws x %d species x %d grids x %d periods", nd, ns, ng, np))

# 后验均值占域面（多数分解在均值上做，CI 由 draw 循环给出）
psi_mean <- apply(psi, c(2, 3, 4), mean)   # species × sites × periods

# ── 2. 恒等式核对：每格物种数均值 == (S/N) × 平均 AOO ────────────────
Rbar <- sapply(seq_len(np), function(t) mean(colSums(psi_mean[, , t])))
AOO  <- sapply(seq_len(np), function(t) rowSums(psi_mean[, , t]))   # species × periods
message("[30] 每格平均物种数: ", paste(round(Rbar, 2), collapse = " -> "))
message("[30] 平均物种 AOO(格): ", paste(round(colMeans(AOO), 1), collapse = " -> "))
message(sprintf("[30] 恒等式校验 |Rbar - (S/N)*meanAOO| max = %.6f",
                max(abs(Rbar - (ns / ng) * colMeans(AOO)))))

# ── 3. 范围扩张 vs 范围内填充分解（P1 -> P5）────────────────────────
decomp_one <- function(thresh) {
  core <- psi_mean[, , 1] > thresh                    # species × sites 逻辑核心区
  dpsi <- psi_mean[, , np] - psi_mean[, , 1]
  infill <- rowSums(dpsi * core)                      # 核心区内的变化
  expand <- rowSums(dpsi * !core)                     # 核心区外的变化
  data.frame(species = species, threshold = thresh,
             d_AOO = rowSums(dpsi), infilling = infill, expansion = expand,
             stringsAsFactors = FALSE)
}
dec <- do.call(rbind, lapply(c(CORE_THRESH, CORE_THRESH_ALT), decomp_one))
write_csv(dec, paste0(v3_file("results", paste0("table_range_decomposition_", run_label)), ".csv"))

summ <- dec |>
  group_by(threshold) |>
  summarise(total_dAOO = sum(d_AOO),
            total_infilling = sum(infilling),
            total_expansion = sum(expansion),
            pct_infilling = 100 * sum(infilling) / sum(d_AOO),
            pct_expansion = 100 * sum(expansion) / sum(d_AOO),
            n_species_net_gain = sum(d_AOO > 0), .groups = "drop")
print(as.data.frame(summ))
write_csv(summ, paste0(v3_file("results", paste0("table_range_decomposition_summary_", run_label)), ".csv"))
message("[30] >>> 若 pct_expansion 占优 = 边缘扩张主导（气候位移特征）；")
message("[30] >>> 若 pct_infilling 占优 = 范围内填充主导（种群恢复特征）。")

# ── 4. 定殖-灭绝分解（相邻期，格点独立近似）─────────────────────────
ce <- do.call(rbind, lapply(seq_len(np - 1), function(t) {
  a <- psi_mean[, , t]; b <- psi_mean[, , t + 1]
  data.frame(period_pair = sprintf("P%d_P%d", t, t + 1),
             expected_colonizations = sum((1 - a) * b),
             expected_extinctions   = sum(a * (1 - b)),
             net = sum(b) - sum(a), stringsAsFactors = FALSE)
}))
ce$colonization_extinction_ratio <- ce$expected_colonizations / ce$expected_extinctions
print(ce); write_csv(ce, paste0(v3_file("results", paste0("table_colonization_extinction_", run_label)), ".csv"))

# ── 5. 分布重心位移（纬度 / 海拔）——气候驱动的方向性检验 ────────────
genv <- safe_read(v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds"))
if (is.null(genv)) genv <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
if (is.null(genv)) stop("[30] grid_environment not found.")
idx <- match(sites, genv$grid_cell)
lat <- genv$lat[idx]; lon <- genv$lon[idx]
elev <- if ("elevation_mean" %in% names(genv)) genv$elevation_mean[idx] else genv$elevation[idx]
bio1 <- if ("bio1" %in% names(genv)) genv$bio1[idx] else genv$BIO1[idx]

wmean <- function(w, x) sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
cent <- do.call(rbind, lapply(seq_len(ns), function(s) {
  data.frame(species = species[s],
             lat_P1  = wmean(psi_mean[s, , 1], lat),
             lat_P5  = wmean(psi_mean[s, , np], lat),
             elev_P1 = wmean(psi_mean[s, , 1], elev),
             elev_P5 = wmean(psi_mean[s, , np], elev),
             stringsAsFactors = FALSE)
}))
cent$d_lat_deg <- cent$lat_P5 - cent$lat_P1
cent$d_elev_m  <- cent$elev_P5 - cent$elev_P1
write_csv(cent, paste0(v3_file("results", paste0("table_centroid_shift_", run_label)), ".csv"))
message(sprintf("[30] 重心纬度位移: 中位 %+.3f deg（北移为正），上移物种比例 %.1f%%",
                median(cent$d_lat_deg, na.rm = TRUE), 100 * mean(cent$d_lat_deg > 0, na.rm = TRUE)))
message(sprintf("[30] 重心海拔位移: 中位 %+.1f m，上移物种比例 %.1f%%",
                median(cent$d_elev_m, na.rm = TRUE), 100 * mean(cent$d_elev_m > 0, na.rm = TRUE)))
print(t.test(cent$d_lat_deg)); print(t.test(cent$d_elev_m))

# ── 6. 群落温度指数 CTI（气候指纹）──────────────────────────────────
# STI 用 P1 占域加权年均温定义，避免用全期数据造成循环论证
STI <- sapply(seq_len(ns), function(s) wmean(psi_mean[s, , 1], bio1))
cti <- do.call(rbind, lapply(seq_len(np), function(t) {
  w <- psi_mean[, , t]                       # species × sites
  data.frame(period = paste0("P", t), grid_cell = sites,
             CTI = colSums(w * STI, na.rm = TRUE) / colSums(w, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
write_csv(cti, paste0(v3_file("results", paste0("table_cti_", run_label)), ".csv"))
cti_nat <- tapply(cti$CTI, cti$period, mean, na.rm = TRUE)
message("[30] 全国 CTI 轨迹: ", paste(sprintf("%.3f", cti_nat), collapse = " -> "))
message("[30] >>> CTI 持续上升 = 群落热适应种占比上升 = 气候驱动的经典指纹。")

# ── 7. 生境类群归因：谁贡献了每格物种数的增量 ───────────────────────
# 修复：原写法把 run_label 传成了 v3_file 的扩展名参数，拼出的路径不存在，
# 只能靠下一行回退。Fixed: run_label was passed as v3_file's extension argument.
tr_path <- paste0(v3_file("results", paste0("table_species_trend_traits_", run_label, "_extended")), ".csv")
if (!file.exists(tr_path)) tr_path <- paste0(v3_file("results", paste0("table_species_trend_traits_", run_label)), ".csv")
if (file.exists(tr_path)) {
  tr <- read_csv(tr_path, show_col_types = FALSE)
  gcol <- intersect(c("Habitat", "habitat", "Primary.Lifestyle", "Trophic.Level"), names(tr))[1]
  if (!is.na(gcol)) {
    contrib <- data.frame(species = species, d_AOO = AOO[, np] - AOO[, 1], stringsAsFactors = FALSE) |>
      left_join(tr |> select(species, guild = all_of(gcol)), by = "species") |>
      mutate(guild = ifelse(is.na(guild), "unknown", guild)) |>
      group_by(guild) |>
      summarise(n_species = n(),
                d_richness_per_grid = sum(d_AOO) / ng,     # 该类群对每格增量的贡献
                .groups = "drop") |>
      mutate(pct_of_total = 100 * d_richness_per_grid / sum(d_richness_per_grid)) |>
      arrange(desc(d_richness_per_grid))
    print(as.data.frame(contrib))
    write_csv(contrib, paste0(v3_file("results", paste0("table_guild_contribution_", run_label)), ".csv"))
    message("[30] >>> 若增量高度集中于少数类群（如水鸟），即可解释功能多样性为何未同步上升。")
  } else message("[30] 未找到类群列，跳过类群归因。")
} else message("[30] 未找到性状表，跳过类群归因。")

log_time("30", "DONE")
