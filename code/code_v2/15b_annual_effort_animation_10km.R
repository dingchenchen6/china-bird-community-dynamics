#!/usr/bin/env Rscript
## 15b_annual_effort_animation_10km.R
##
## 独立脚本：仅处理 10km 分辨率年度动态地图。
## 假设 100km 部分已由 15_annual_effort_animation.R 完成。
##
## 输入：data/derived_v2/combined_events_merged_dedup_2000_2025.rds
##       data/derived_v2/china_grid_10km_v2.rds（若无则现场构建）
## 输出：results_v2/table_effort_by_grid_year_source_10km.csv
##       figures_v2/anim_effort_events_10km_combined.gif
##       figures_v2/anim_effort_events_10km_by_source.gif

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(sf)
  library(ggplot2)
  library(gifski)
  library(scales)
  library(scico)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_spatial.R"))
source(file.path(CODE_V2, "utils_mapping.R"))

P <- ensure_v2_dirs()

ANALYSIS_YR_LO <- 2000L
ANALYSIS_YR_HI <- 2024L

## --- 辅助函数 ---------------------------------------------------------------

save_gif_from_plots <- function(plot_list, out_path, width = 1200, height = 800,
                                 res = 150, delay = 0.5) {
  tmpdir <- tempfile(pattern = "gif_frames_")
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  frame_files <- file.path(tmpdir, sprintf("frame_%03d.png", seq_along(plot_list)))
  for (i in seq_along(plot_list)) {
    png(frame_files[i], width = width, height = height, res = res)
    print(plot_list[[i]])
    dev.off()
  }
  gifski::gifski(frame_files, gif_file = out_path,
                 width = width, height = height, delay = delay)
  message("  GIF saved: ", out_path)
}

## --- 1. 读取事件数据 --------------------------------------------------------

message("[stage-15b] Loading dedup'd events...")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                                "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]
message(sprintf("  events %d-%d: %s",
                ANALYSIS_YR_LO, ANALYSIS_YR_HI,
                format(nrow(events), big.mark = ",")))

events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]
events[, day_key := sprintf("%04d-%02d-%02d", year, month, day)]
events[, observer_key := tolower(coalesce(username, "anon"))]
events[, visit_key := paste(observer_key, day_key, sep = "|")]
events[, has_dur := !is.na(duration_min)]

## --- 2. 10km 网格 -----------------------------------------------------------

message("[stage-15b] Loading/building 10km grid...")
grid_10_cache <- file.path(P$derived_v2, "china_grid_10km_v2.rds")
if (file.exists(grid_10_cache)) {
  grid_10 <- readRDS(grid_10_cache)
  message(sprintf("  10km grid loaded from cache: %d cells", nrow(grid_10)))
} else {
  grid_10 <- load_or_build_v2_grid(P$china_boundary, resolution_m = 1e4,
                                   cache_path = grid_10_cache)
  message(sprintf("  10km grid built: %d cells", nrow(grid_10)))
}

## --- 3. 使用 terra 快速分配事件到 10km 网格 ----------------------------------

message("[stage-15b] Fast spatial assignment with terra...")

## 将 10km 网格转为 terra::vect，并创建覆盖全中国的 raster template
## 策略：利用 grid_10 的 bbox + 分辨率，直接用 cellFromXY 计算单元格ID

## 获取 grid_10 的 3857 投影 bbox 和分辨率（grid_10 当前是 4326）
grid_10_proj <- sf::st_transform(grid_10, 3857)
bb_3857 <- sf::st_bbox(grid_10_proj)
res_m <- 1e4  # 10km

## 将事件点转为 3857 坐标
pts_3857 <- sf::st_as_sf(events[, .(longitude, latitude)],
                         coords = c("longitude", "latitude"),
                         crs = 4326, remove = FALSE) |>
  sf::st_transform(3857)
coords_3857 <- sf::st_coordinates(pts_3857)

## 计算行列号（基于 grid_10 的实际 bbox）
col_idx <- floor((coords_3857[, 1] - bb_3857$xmin) / res_m) + 1
row_idx <- floor((coords_3857[, 2] - bb_3857$ymin) / res_m) + 1

## 构建一个从 (col, row) -> grid_cell 的查找表
## 先计算 grid_10_proj 中每个 cell 的 (col, row)
centroids_3857 <- sf::st_coordinates(sf::st_centroid(grid_10_proj))
cols <- floor((centroids_3857[, 1] - bb_3857$xmin) / res_m) + 1
rows <- floor((centroids_3857[, 2] - bb_3857$ymin) / res_m) + 1

## 使用 match + 唯一字符串键，避免 data.table cartesian join
lookup <- unique(data.table(grid_cell = grid_10$grid_cell, col = cols, row = rows))
lookup[, cr_key := paste(col, row, sep = "_")]
event_cr_key <- paste(col_idx, row_idx, sep = "_")
match_idx <- match(event_cr_key, lookup$cr_key)
events[, grid_cell_10km := lookup$grid_cell[match_idx]]

events_10 <- events[!is.na(grid_cell_10km)]
message(sprintf("  events matched to 10km grid: %s (%0.1f%%)",
                format(nrow(events_10), big.mark = ","),
                100 * nrow(events_10) / nrow(events)))

rm(pts_3857, coords_3857, col_idx, row_idx, lookup, matched)
gc()

## --- 4. 按 (grid × year × source) 汇总 --------------------------------------

message("[stage-15b] Aggregating effort by grid-year-source (10km)...")

eff_10_source <- events_10[, .(
  n_events    = uniqueN(visit_key),
  n_observers = uniqueN(observer_key),
  n_dates     = uniqueN(day_key),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE)
), by = .(grid_cell_10km, year, source_short)]
setnames(eff_10_source, "grid_cell_10km", "grid_cell")

eff_10_combined <- events_10[, .(
  n_events    = uniqueN(visit_key),
  n_observers = uniqueN(observer_key),
  n_dates     = uniqueN(day_key),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE)
), by = .(grid_cell_10km, year)][, source_short := "Combined"]
setnames(eff_10_combined, "grid_cell_10km", "grid_cell")

eff_10_all <- rbindlist(list(eff_10_source, eff_10_combined), use.names = TRUE)
eff_10_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
setorder(eff_10_all, grid_cell, year, source_short)

write_csv(as_tibble(eff_10_all),
          v2_file("results", "table_effort_by_grid_year_source_10km"))
message("  Saved: results_v2/table_effort_by_grid_year_source_10km.csv")

## 清理
gc()

## --- 5. 10km 动态地图（GIF）--------------------------------------------------

message("[stage-15b] Rendering 10km annual animated maps...")

china_layers <- load_china_layers(P$china_boundary, P$province_line)
grid_10_lean <- grid_10 |> select(grid_cell)
years <- sort(unique(eff_10_all$year))

## 5a. Combined
message("  -> Combined n_events animation (10km)...")
global_max_10_combined <- log10(max(eff_10_all[source_short == "Combined"]$n_events, na.rm = TRUE) + 1)

plots_10_combined <- list()
for (i in seq_along(years)) {
  yr <- years[i]
  yr_dt <- eff_10_all[year == yr & source_short == "Combined", .(grid_cell, n_events)]
  map_dt <- grid_10_lean |>
    left_join(yr_dt, by = "grid_cell") |>
    mutate(fill_val = log10(n_events + 1))

  plots_10_combined[[i]] <- ggplot(map_dt) +
    geom_sf(aes(fill = fill_val),
            colour = alpha("white", 0.03), linewidth = 0.01) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_sequential(
      name = "log10(events+1)", palette = "lajolla", direction = 1,
      limits = c(0, global_max_10_combined), oob = squish, na.value = "grey92"
    ) +
    v2_china_coord() +
    theme_v2_map(11) +
    labs(
      title = "Annual survey effort — Combined (merged + deduplicated)",
      subtitle = sprintf("Year: %d  |  10 km grid  |  Events per grid cell", yr),
      caption = sprintf("Max log10(events+1) across all years = %.2f", global_max_10_combined)
    )
}

save_gif_from_plots(
  plots_10_combined,
  file.path(P$figures_v2, "anim_effort_events_10km_combined.gif"),
  width = 1200, height = 900, delay = 0.6
)

## 5b. By-source
message("  -> By-source n_events animation (10km)...")
global_max_10_src <- log10(max(eff_10_all$n_events, na.rm = TRUE) + 1)

plots_10_bysource <- list()
idx <- 1
for (yr in years) {
  yr_dt <- eff_10_all[year == yr, .(grid_cell, source_short, n_events)]
  map_dt <- grid_10_lean |>
    left_join(yr_dt, by = "grid_cell") |>
    mutate(fill_val = log10(n_events + 1))

  plots_10_bysource[[idx]] <- ggplot(map_dt) +
    geom_sf(aes(fill = fill_val),
            colour = alpha("white", 0.02), linewidth = 0.005) +
    china_map_layers(china_layers, country_lwd = 0.20, province_lwd = 0.10) +
    scale_fill_v2_sequential(
      name = "log10(events+1)", palette = "lajolla", direction = 1,
      limits = c(0, global_max_10_src), oob = squish, na.value = "grey92"
    ) +
    facet_wrap(~ source_short, ncol = 3) +
    v2_china_coord() +
    theme_v2_map(10) +
    labs(
      title = sprintf("Annual survey effort by data source — Year %d", yr),
      subtitle = "10 km grid  |  log10(events + 1)",
      caption = "Grey = no data."
    )
  idx <- idx + 1
}

save_gif_from_plots(
  plots_10_bysource,
  file.path(P$figures_v2, "anim_effort_events_10km_by_source.gif"),
  width = 1800, height = 700, delay = 0.6
)

message("[stage-15b] 10km animated maps completed.")
