#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   100 km 网格是否过粗？50 km 网格的占域趋势方向与 100 km 是否一致？
#   （D3 / 审稿人必问的尺度敏感性）
#
# Objective / 分析目标:
#   - 在 50 km 网格上跑 top-60 sp 的 stMsPGOcc pilot
#   - 把 100 km 结果按空间聚合到 50 km 子网格再比较
#   - 输出 table_grid_size_sensitivity_50km_vs_100km.csv：
#       grid_cell_100km, trend_100km, trend_50km_median, direction_consistent
#
# Input data / 输入数据:
#   results_v4/table_trend_summary_<run_label>.csv（100 km）
#   data/derived_v4/stMsPGOcc_fit_<v4_pilot_60sp_ar1_spatial_50km>.qs（50 km）
#
# Main workflow / 主要流程:
#   1. 用 V4_GRID_SIZE_KM=50 + V4_PILOT=1 跑一次 02/03/04（外部执行）
#   2. 本脚本读两套 trend，按 50km 网格质心落到 100km 网格做聚合
#   3. 计算方向一致性 + 相关性
#
# Key assumptions / 关键假设:
#   50 km 网格 pilot 已跑完
#   100 km full / pilot 也已跑完
#
# Main packages / 主要包:
#   sf, dplyr, tidyr, readr
#
# Output directory / 输出路径:
#   results_v4/table_grid_size_sensitivity_50km_vs_100km.csv
#   figures_v4/fig_v4_grid_size_sensitivity.{png,pdf}
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(sf); library(ggplot2)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_mapping.R"))
P <- ensure_v4_dirs()

log_time("16", "Grid size sensitivity 50km vs 100km")

# ── 加载两套 trend（用 OLS richness trend）──────────────────────────
trd_100 <- read_csv_safe(v4_file("results",
                                   paste0("table_trend_summary_", RUN_LABEL))) %||%
            read_csv_safe(v4_file("results",
                                   paste0("table_trend_summary_", PILOT_LABEL)))

trd_50_label <- gsub("_100km$|$", "_50km", RUN_LABEL)
trd_50_label <- gsub("(_v4_full_200sp_ar1_spatial)(_50km)?", "\\1_50km", trd_50_label)
trd_50 <- read_csv_safe(v4_file("results", paste0("table_trend_summary_", trd_50_label)))
if (is.null(trd_50)) {
  trd_50 <- read_csv_safe(v4_file("results",
                                    "table_trend_summary_v4_pilot_60sp_ar1_spatial_50km"))
}

if (is.null(trd_100) || is.null(trd_50)) {
  stop("[16] Need both 100km and 50km trend tables. Run 04+05 with V4_GRID_SIZE_KM=50 first.")
}

# ── 准备：100 km / 50 km 网格 sf ────────────────────────────────────
grid_100 <- safe_read(v4_file("derived", "china_grid_100km_v4", "rds"), quiet = TRUE) %||%
             safe_read(v3_file("derived", "china_grid_100km_v3", "rds"), quiet = TRUE) %||%
             safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"), quiet = TRUE)
grid_50 <- safe_read(v4_file("derived", "china_grid_50km_v4", "rds"), quiet = TRUE) %||%
            safe_read(v3_file("derived", "china_grid_50km_v3", "rds"), quiet = TRUE)
if (is.null(grid_100) || is.null(grid_50)) {
  stop("[16] Missing grid sf for 100km or 50km. Run 02 with both grid sizes first.")
}

# ── 把 50 km 网格质心落到 100 km 网格 ───────────────────────────────
sf_use_s2(FALSE)
ctr_50 <- suppressWarnings(st_centroid(grid_50))
joined <- st_join(ctr_50, grid_100 |> select(grid_cell_100 = grid_cell),
                  join = st_intersects)
mapping <- tibble(
  grid_cell = grid_50$grid_cell,
  grid_cell_100 = joined$grid_cell_100
) |> drop_na()

# ── 聚合 50 km trend 到 100 km ──────────────────────────────────────
metric_focus <- "corrected_richness"
trd_50_rich <- trd_50 |>
  filter(metric == metric_focus, method == "theil_sen") |>
  select(grid_cell, trend_50 = mean) |>
  left_join(mapping, by = "grid_cell") |>
  drop_na(grid_cell_100) |>
  group_by(grid_cell_100) |>
  summarise(
    trend_50_median = median(trend_50, na.rm = TRUE),
    trend_50_n = n(),
    .groups = "drop"
  )

trd_100_rich <- trd_100 |>
  filter(metric == metric_focus, method == "theil_sen") |>
  select(grid_cell_100 = grid_cell, trend_100 = mean)

comp <- trd_100_rich |>
  inner_join(trd_50_rich, by = "grid_cell_100") |>
  drop_na()

comp <- comp |>
  mutate(
    sign_100 = sign(trend_100),
    sign_50  = sign(trend_50_median),
    direction_consistent = sign_100 == sign_50 & sign_100 != 0
  )

n_total <- nrow(comp)
n_consistent <- sum(comp$direction_consistent, na.rm = TRUE)
rho <- cor(comp$trend_100, comp$trend_50_median,
            use = "complete.obs", method = "spearman")
message(sprintf("[16] Direction consistency: %d/%d (%.1f%%); Spearman rho = %.3f",
                n_consistent, n_total, 100 * n_consistent / n_total, rho))

write_csv(comp, v4_file("results", "table_grid_size_sensitivity_50km_vs_100km"))

# ── 配对图 ──────────────────────────────────────────────────────────
p <- ggplot(comp, aes(trend_100, trend_50_median, colour = direction_consistent)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = "grey60", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.15) +
  geom_point(alpha = 0.55, size = 0.8) +
  scale_colour_manual(values = c("FALSE" = "#B2182B", "TRUE" = "grey40"),
                      labels = c("FALSE" = "direction flipped",
                                 "TRUE"  = "consistent"),
                      name = NULL) +
  labs(x = "Richness trend (100 km)",
       y = "Richness trend (50 km, median over child grids)",
       title = sprintf("Grid-size sensitivity: %d/%d consistent, ρ = %.3f",
                       n_consistent, n_total, rho)) +
  theme_nature_pub()
save_nature(p, "fig_v4_grid_size_sensitivity",
            width_mm = NATURE_WIDTH_M, height_mm = 90)

log_time("16", "DONE")
