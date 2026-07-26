#!/usr/bin/env Rscript
## 02e_effort_annual_tables_animaps.R
##
## 不同数据源、不同调查努力指标的年度统计表 + 年度动态地图（GIF）。
## 产出：
##   results_v2/table_effort_annual_by_source.csv      — 年度 × 数据源 × 指标宽表
##   results_v2/table_effort_annual_by_source_long.csv  — 长表（适合可视化）
##   figures_v2/fig_effort_annual_100km_*.gif           — 100km 网格年度动态地图
##   figures_v2/fig_effort_annual_100km_*_by_source.gif — 100km 分源版
##   figures_v2/fig_effort_annual_10km_*.gif            — 10km 网格年度动态地图
##   figures_v2/fig_effort_annual_10km_*_by_source.gif  — 10km 分源版

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(sf)
  library(terra)
  library(ggplot2)
  library(gifski)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_spatial.R"))
source(file.path(CODE_V2, "utils_mapping.R"))

P <- ensure_v2_dirs()

ANALYSIS_YR_LO <- 2000L
ANALYSIS_YR_HI <- 2024L
years <- ANALYSIS_YR_LO:ANALYSIS_YR_HI

## =============================================================================
## 0. 加载去重事件 + 预处理
## =============================================================================

message("[02e] Loading dedup'd events")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                            "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]
message(sprintf("  events 2000-2024: %s", format(nrow(events), big.mark = ",")))

events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]
events[, day_key := sprintf("%04d-%02d-%02d", year, month, day)]
events[, observer_key := tolower(coalesce(username, "anon"))]
events[, visit_key := paste(observer_key, day_key, sep = "|")]
events[, has_dur := !is.na(duration_min)]

china_layers <- load_china_layers(P$china_boundary, P$province_line)

## =============================================================================
## 1. 100km 网格：加载 + 事件分配 + 年度统计 + 动画
## =============================================================================

message("[02e] === 100 km resolution ===")

grid_cache_100 <- file.path(P$derived_v2, "china_grid_100km_v2.rds")
grid_100 <- if (file.exists(grid_cache_100)) readRDS(grid_cache_100) else
  load_or_build_v2_grid(P$china_boundary, 1e5, grid_cache_100)
message(sprintf("  100km grid: %d cells", nrow(grid_100)))

# 如果事件已有 grid_cell（来自 stage-2 缓存），跳过空间 join
if (!"grid_cell" %in% names(events) || any(is.na(events$grid_cell))) {
  message("  assigning events to 100km grid (sf planar)")
  old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
  pts <- sf::st_as_sf(events[, .(longitude, latitude)],
                      coords = c("longitude", "latitude"),
                      crs = 4326, remove = FALSE)
  joined <- sf::st_join(pts, grid_100[, "grid_cell"],
                        join = sf::st_intersects, left = TRUE)
  events[, grid_cell := joined$grid_cell]
  sf::sf_use_s2(old_s2)
}
events_100 <- events[!is.na(grid_cell)]

## --- 1a. 年度统计表 ---------------------------------------------------------

message("[02e] Building annual effort tables by source (100km)")

annual_by_source <- events_100[, .(
  n_records         = .N,
  n_visits          = uniqueN(visit_key),
  n_observers       = uniqueN(observer_key),
  n_birding_days    = uniqueN(day_key),
  n_grids_100km     = uniqueN(grid_cell),
  n_species         = uniqueN(species),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE),
  total_duration_min = sum(duration_min[has_dur], na.rm = TRUE)
), by = .(year, source_short)]

annual_combined <- events_100[, .(
  n_records         = .N,
  n_visits          = uniqueN(visit_key),
  n_observers       = uniqueN(observer_key),
  n_birding_days    = uniqueN(day_key),
  n_grids_100km     = uniqueN(grid_cell),
  n_species         = uniqueN(species),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE),
  total_duration_min = sum(duration_min[has_dur], na.rm = TRUE)
), by = .(year)][, source_short := "Combined"]

annual_all <- rbindlist(list(annual_by_source, annual_combined), use.names = TRUE)
annual_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
setorder(annual_all, source_short, year)

write_csv(as_tibble(annual_all),
          file.path(P$results_v2, "table_effort_annual_by_source.csv"))

annual_long <- as_tibble(annual_all) |>
  pivot_longer(cols = c(n_records, n_visits, n_observers, n_birding_days,
                        n_grids_100km, n_species, mean_duration_min, total_duration_min),
               names_to = "metric", values_to = "value")
write_csv(annual_long,
          file.path(P$results_v2, "table_effort_annual_by_source_long.csv"))

message(sprintf("  Annual table: %d rows, %d sources", nrow(annual_all), uniqueN(annual_all$source_short)))

## --- 1b. 100km 网格 × 年 × 数据源（空间统计） ------------------------------

message("[02e] Computing 100km grid x year effort")

grid_year_src <- events_100[, .(
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_species      = uniqueN(species),
  n_birding_days = uniqueN(day_key)
), by = .(grid_cell, year, source_short)]

grid_year_comb <- events_100[, .(
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_species      = uniqueN(species),
  n_birding_days = uniqueN(day_key)
), by = .(grid_cell, year)][, source_short := "Combined"]

grid_year_all <- rbindlist(list(grid_year_src, grid_year_comb), use.names = TRUE)
grid_year_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]

## --- 1c. 100km 动画帧 -------------------------------------------------------

message("[02e] Rendering 100km annual animation frames")

anim_dir_100 <- file.path(P$figures_v2, "anim_100km")
dir.create(anim_dir_100, recursive = TRUE, showWarnings = FALSE)

effort_100_comb <- grid_year_all[source_short == "Combined"]
visits_max_100   <- quantile(effort_100_comb$n_visits, 0.99, na.rm = TRUE)
obs_max_100      <- quantile(effort_100_comb$n_observers, 0.99, na.rm = TRUE)
species_max_100  <- quantile(effort_100_comb$n_species, 0.99, na.rm = TRUE)
days_max_100     <- quantile(effort_100_comb$n_birding_days, 0.99, na.rm = TRUE)

render_100km_frame <- function(yr, metric_col, fill_label, fill_max, stem) {
  sub <- grid_year_all[year == yr & source_short == "Combined"]
  if (nrow(sub) == 0) return(NULL)

  sf_yr <- grid_100 |> left_join(as_tibble(sub), by = "grid_cell")
  fv <- rlang::sym(metric_col)

  p <- ggplot(sf_yr) +
    geom_sf(aes(fill = !!fv),
            colour = scales::alpha("white", 0.08), linewidth = 0.03) +
    china_map_layers(china_layers) +
    scale_fill_v2_sequential(name = fill_label, palette = "lajolla", direction = 1,
                             limits = c(0, fill_max), oob = scales::squish,
                             na.value = "grey95") +
    v2_china_coord() + theme_v2_map(11) +
    labs(title = sprintf("Survey effort — %d", yr),
         subtitle = sprintf("%s per 100 km grid (Combined)", fill_label))

  out <- file.path(anim_dir_100, sprintf("frame_%s_%04d.png", stem, yr))
  ggsave(out, p, width = 9.4, height = 6.6, dpi = 200, bg = "white")
  out
}

render_100km_bysrc <- function(yr, metric_col, fill_label, fill_max, stem) {
  sub <- grid_year_all[year == yr]
  if (nrow(sub) == 0) return(NULL)

  sf_yr <- grid_100 |> left_join(as_tibble(sub), by = "grid_cell")
  fv <- rlang::sym(metric_col)

  p <- ggplot(sf_yr) +
    geom_sf(aes(fill = !!fv),
            colour = scales::alpha("white", 0.06), linewidth = 0.03) +
    china_map_layers(china_layers) +
    scale_fill_v2_sequential(name = fill_label, palette = "lajolla", direction = 1,
                             limits = c(0, fill_max), oob = scales::squish,
                             na.value = "grey95") +
    facet_wrap(~ source_short, ncol = 3) +
    v2_china_coord() + theme_v2_map(9) +
    labs(title = sprintf("Survey effort by source — %d", yr),
         subtitle = sprintf("%s per 100 km grid", fill_label))

  out <- file.path(anim_dir_100, sprintf("frame_bysrc_%s_%04d.png", stem, yr))
  ggsave(out, p, width = 16, height = 6, dpi = 200, bg = "white")
  out
}

metrics_100 <- list(
  list(col = "n_visits",       label = "Visit events",     max = visits_max_100,  stem = "visits"),
  list(col = "n_observers",    label = "Unique observers",  max = obs_max_100,     stem = "observers"),
  list(col = "n_species",      label = "Species detected",  max = species_max_100, stem = "species"),
  list(col = "n_birding_days", label = "Birding days",     max = days_max_100,    stem = "days")
)

for (m in metrics_100) {
  message(sprintf("  100km Combined: %s", m$stem))
  frames <- unlist(lapply(years, function(yr)
    render_100km_frame(yr, m$col, m$label, m$max, m$stem)))
  gifski::gifski(frames,
                 file.path(P$figures_v2, sprintf("fig_effort_annual_100km_%s.gif", m$stem)),
                 width = 1880, height = 1320, delay = 0.8, loop = TRUE)

  message(sprintf("  100km by-source: %s", m$stem))
  frames_src <- unlist(lapply(years, function(yr)
    render_100km_bysrc(yr, m$col, m$label, m$max, m$stem)))
  gifski::gifski(frames_src,
                 file.path(P$figures_v2, sprintf("fig_effort_annual_100km_%s_by_source.gif", m$stem)),
                 width = 3200, height = 1200, delay = 0.8, loop = TRUE)
}
message("[02e] 100km GIFs done.")

## =============================================================================
## 2. 10km 网格：构建 + 快速分配 + 年度统计 + 动画
## =============================================================================

message("[02e] === 10 km resolution ===")

grid_cache_10 <- file.path(P$derived_v2, "china_grid_10km_v2.rds")
grid_10 <- if (file.exists(grid_cache_10)) {
  readRDS(grid_cache_10)
} else {
  load_or_build_v2_grid(P$china_boundary, 1e4, grid_cache_10)
}
message(sprintf("  10km grid: %d cells", nrow(grid_10)))

## --- 2a. 用 terra::extract 快速分配点到 10km 网格 -------------------------
## 构建与 create_china_grid 一致的 raster 模板，用 terra::extract 向量化分配。

message("[02e] Assigning events to 10km grid via terra::extract")

# 重构原始 raster 模板参数（与 create_china_grid 一致）
china_proj <- sf::st_transform(
  suppressWarnings(sf::st_read(P$china_boundary, quiet = TRUE)) |> sf::st_make_valid(),
  3857
)
bb <- sf::st_bbox(china_proj)
aligned_ext <- align_extent(bb$xmin, bb$xmax, bb$ymin, bb$ymax, 1e4)
raster_10k <- terra::rast(ext = aligned_ext, resolution = 1e4,
                           crs = sf::st_crs(china_proj)$wkt)
terra::values(raster_10k) <- seq_len(terra::ncell(raster_10k))

# 将事件坐标投影到 3857，再用 terra::extract 分配 cell
pts_3857 <- sf::st_transform(
  sf::st_as_sf(events_100[, .(longitude, latitude)],
               coords = c("longitude", "latitude"), crs = 4326),
  3857
)
pts_vect <- terra::vect(pts_3857)

# 分批 extract 防止内存问题
batch <- 2e6
n_ev <- nrow(events_100)
grid_cell_10km <- integer(n_ev)

for (i in seq(1, n_ev, by = batch)) {
  end_i <- min(i + batch - 1, n_ev)
  idx <- i:end_i
  ext_batch <- terra::extract(raster_10k, pts_vect[idx], fun = NULL, na.rm = FALSE)
  grid_cell_10km[idx] <- ext_batch[[2]]  # 第 2 列是 lyr.1 = grid_cell
  if (i + batch > n_ev) message(sprintf("  10km extract: %d / %d", end_i, n_ev))
}
events_100[, grid_cell_10km := grid_cell_10km]

# 只保留落在 grid_10 范围内的事件
valid_cells_10 <- grid_10$grid_cell
events_10 <- events_100[grid_cell_10km %in% valid_cells_10]
message(sprintf("  events with valid 10km grid_cell: %s",
                format(nrow(events_10), big.mark = ",")))

rm(pts_3857, pts_vect, raster_10k, china_proj)
gc()

## --- 2b. 10km 网格 × 年 × 数据源 ------------------------------------------

message("[02e] Computing 10km grid x year effort")

grid10_year_src <- events_10[, .(
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_species      = uniqueN(species),
  n_birding_days = uniqueN(day_key)
), by = .(grid_cell = grid_cell_10km, year, source_short)]

grid10_year_comb <- events_10[, .(
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_species      = uniqueN(species),
  n_birding_days = uniqueN(day_key)
), by = .(grid_cell = grid_cell_10km, year)][, source_short := "Combined"]

grid10_year_all <- rbindlist(list(grid10_year_src, grid10_year_comb), use.names = TRUE)
grid10_year_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]

# 保存 10km 网格年度努力（方便后续使用）
write_csv(as_tibble(grid10_year_all),
          file.path(P$results_v2, "table_effort_10km_grid_year_source.csv"))

## --- 2c. 10km 动画帧 -------------------------------------------------------

message("[02e] Rendering 10km annual animation frames")

anim_dir_10 <- file.path(P$figures_v2, "anim_10km")
dir.create(anim_dir_10, recursive = TRUE, showWarnings = FALSE)

effort_10_comb <- grid10_year_all[source_short == "Combined"]
visits_max_10   <- quantile(effort_10_comb$n_visits, 0.99, na.rm = TRUE)
obs_max_10      <- quantile(effort_10_comb$n_observers, 0.99, na.rm = TRUE)
species_max_10  <- quantile(effort_10_comb$n_species, 0.99, na.rm = TRUE)
days_max_10     <- quantile(effort_10_comb$n_birding_days, 0.99, na.rm = TRUE)

render_10km_frame <- function(yr, metric_col, fill_label, fill_max, stem) {
  sub <- grid10_year_all[year == yr & source_short == "Combined"]
  if (nrow(sub) == 0) return(NULL)

  sf_yr <- grid_10 |> left_join(as_tibble(sub), by = "grid_cell")
  fv <- rlang::sym(metric_col)

  p <- ggplot(sf_yr) +
    geom_sf(aes(fill = !!fv), colour = NA, linewidth = 0) +
    china_map_layers(china_layers) +
    scale_fill_v2_sequential(name = fill_label, palette = "lajolla", direction = 1,
                             limits = c(0, fill_max), oob = scales::squish,
                             na.value = "grey95") +
    v2_china_coord() + theme_v2_map(11) +
    labs(title = sprintf("Survey effort — %d (10 km)", yr),
         subtitle = sprintf("%s per 10 km grid (Combined)", fill_label))

  out <- file.path(anim_dir_10, sprintf("frame_%s_%04d.png", stem, yr))
  ggsave(out, p, width = 9.4, height = 6.6, dpi = 200, bg = "white")
  out
}

render_10km_bysrc <- function(yr, metric_col, fill_label, fill_max, stem) {
  sub <- grid10_year_all[year == yr]
  if (nrow(sub) == 0) return(NULL)

  sf_yr <- grid_10 |> left_join(as_tibble(sub), by = "grid_cell")
  fv <- rlang::sym(metric_col)

  p <- ggplot(sf_yr) +
    geom_sf(aes(fill = !!fv), colour = NA, linewidth = 0) +
    china_map_layers(china_layers) +
    scale_fill_v2_sequential(name = fill_label, palette = "lajolla", direction = 1,
                             limits = c(0, fill_max), oob = scales::squish,
                             na.value = "grey95") +
    facet_wrap(~ source_short, ncol = 3) +
    v2_china_coord() + theme_v2_map(9) +
    labs(title = sprintf("Survey effort by source — %d (10 km)", yr),
         subtitle = sprintf("%s per 10 km grid", fill_label))

  out <- file.path(anim_dir_10, sprintf("frame_bysrc_%s_%04d.png", stem, yr))
  ggsave(out, p, width = 16, height = 6, dpi = 200, bg = "white")
  out
}

metrics_10 <- list(
  list(col = "n_visits",       label = "Visit events",     max = visits_max_10,  stem = "visits"),
  list(col = "n_observers",    label = "Unique observers",  max = obs_max_10,     stem = "observers"),
  list(col = "n_species",      label = "Species detected",  max = species_max_10, stem = "species"),
  list(col = "n_birding_days", label = "Birding days",     max = days_max_10,    stem = "days")
)

for (m in metrics_10) {
  message(sprintf("  10km Combined: %s", m$stem))
  frames <- unlist(lapply(years, function(yr)
    render_10km_frame(yr, m$col, m$label, m$max, m$stem)))
  gifski::gifski(frames,
                 file.path(P$figures_v2, sprintf("fig_effort_annual_10km_%s.gif", m$stem)),
                 width = 1880, height = 1320, delay = 0.8, loop = TRUE)

  message(sprintf("  10km by-source: %s", m$stem))
  frames_src <- unlist(lapply(years, function(yr)
    render_10km_bysrc(yr, m$col, m$label, m$max, m$stem)))
  gifski::gifski(frames_src,
                 file.path(P$figures_v2, sprintf("fig_effort_annual_10km_%s_by_source.gif", m$stem)),
                 width = 3200, height = 1200, delay = 0.8, loop = TRUE)
}
message("[02e] 10km GIFs done.")

## =============================================================================
## 3. 汇总
## =============================================================================

message("[02e] === All done ===")
message("  Tables:")
message("    table_effort_annual_by_source.csv")
message("    table_effort_annual_by_source_long.csv")
message("    table_effort_10km_grid_year_source.csv")
message("  100km GIFs:")
for (m in metrics_100)
  message(sprintf("    fig_effort_annual_100km_%s.gif + _by_source", m$stem))
message("  10km GIFs:")
for (m in metrics_10)
  message(sprintf("    fig_effort_annual_10km_%s.gif + _by_source", m$stem))
