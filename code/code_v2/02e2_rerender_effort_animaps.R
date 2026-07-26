#!/usr/bin/env Rscript
## 02e2_rerender_effort_animaps.R
##
## 重新渲染年度调查努力动态地图：修复底图问题。
## 修改：NA 网格不填灰色（白色）、省份无填充、确保用 data/中国shp/ 底图。
##
## 产出：figures_v2/fig_effort_annual_100km_*.gif  (覆盖旧版)
##       figures_v2/fig_effort_annual_10km_*.gif   (覆盖旧版)

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(sf)
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
## 0. 加载底图（确保来自 data/中国shp/）
## =============================================================================

message("[rerender] Loading China base layers from data/中国shp/")

# 显式路径，避免回退到其它目录
shp_dir <- file.path(P$task_root, "data", "中国shp")
china_province_shp  <- file.path(shp_dir, "省.shp")
china_provline_shp  <- file.path(shp_dir, "省_境界线.shp")
china_ten_dash_shp  <- file.path(shp_dir, "十段线.shp")

stopifnot("省.shp not found" = file.exists(china_province_shp))

china_layers <- load_china_layers(
  china_boundary_path   = china_province_shp,
  province_line_path    = china_provline_shp,
  ten_dash_path         = china_ten_dash_shp,
  strip_inset           = TRUE,
  province_lines_from_polygons = FALSE
)
message(sprintf("  provinces loaded: %d features", nrow(china_layers$china)))
message(sprintf("  province lines:  %s",
                if (!is.null(china_layers$province)) nrow(china_layers$province) else "NULL"))
message(sprintf("  ten-dash line:   %s",
                if (!is.null(china_layers$ten_dash)) "present" else "NULL"))

## =============================================================================
## 1. 加载去重事件，计算 grid × year × source 聚合（100km + 10km）
## =============================================================================

message("[rerender] Loading events for aggregation")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                            "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]

events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]
events[, day_key := sprintf("%04d-%02d-%02d", year, month, day)]
events[, observer_key := tolower(coalesce(username, "anon"))]
events[, visit_key := paste(observer_key, day_key, sep = "|")]

## --- 100km 聚合 ---------------------------------------------------------------

grid_cache_100 <- file.path(P$derived_v2, "china_grid_100km_v2.rds")
grid_100 <- readRDS(grid_cache_100)

if (!"grid_cell" %in% names(events) || any(is.na(events$grid_cell))) {
  message("  assigning events to 100km grid")
  old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
  pts <- sf::st_as_sf(events[, .(longitude, latitude)],
                      coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  joined <- sf::st_join(pts, grid_100[, "grid_cell"],
                        join = sf::st_intersects, left = TRUE)
  events[, grid_cell := joined$grid_cell]
  sf::sf_use_s2(old_s2)
}
events_100 <- events[!is.na(grid_cell)]

message("[rerender] Aggregating 100km grid x year x source")
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

## --- 10km 聚合 ----------------------------------------------------------------

grid_cache_10 <- file.path(P$derived_v2, "china_grid_10km_v2.rds")
grid_10 <- readRDS(grid_cache_10)

message("[rerender] Loading pre-computed 10km grid x year effort from CSV")
grid10_year_all <- as.data.table(read_csv(
  file.path(P$results_v2, "table_effort_10km_grid_year_source.csv"),
  show_col_types = FALSE))
grid10_year_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]

## =============================================================================
## 2. 统一色标 + 渲染函数（修复底图 + NA=白色）
## =============================================================================

# 100km 色标
eff100 <- grid_year_all[source_short == "Combined"]
vmax_100 <- list(
  n_visits       = quantile(eff100$n_visits, 0.99, na.rm = TRUE),
  n_observers    = quantile(eff100$n_observers, 0.99, na.rm = TRUE),
  n_species      = quantile(eff100$n_species, 0.99, na.rm = TRUE),
  n_birding_days = quantile(eff100$n_birding_days, 0.99, na.rm = TRUE)
)

# 10km 色标
eff10 <- grid10_year_all[source_short == "Combined"]
vmax_10 <- list(
  n_visits       = quantile(eff10$n_visits, 0.99, na.rm = TRUE),
  n_observers    = quantile(eff10$n_observers, 0.99, na.rm = TRUE),
  n_species      = quantile(eff10$n_species, 0.99, na.rm = TRUE),
  n_birding_days = quantile(eff10$n_birding_days, 0.99, na.rm = TRUE)
)

# 通用渲染函数：Combined 版
render_combined <- function(yr, grid_sf, effort_dt, metric_col, fill_label,
                            fill_max, out_dir, stem, resolution_label) {
  sub <- effort_dt[year == yr & source_short == "Combined"]
  if (nrow(sub) == 0) return(NULL)

  sf_yr <- grid_sf |> left_join(as_tibble(sub), by = "grid_cell")
  fv <- rlang::sym(metric_col)

  # 100km 用淡白描边，10km 不描边
  grid_colour <- if (resolution_label == "100 km")
    scales::alpha("white", 0.08) else NA
  grid_lwd <- if (resolution_label == "100 km") 0.03 else 0

  p <- ggplot(sf_yr) +
    geom_sf(aes(fill = !!fv),
            colour = grid_colour, linewidth = grid_lwd) +
    china_map_layers(china_layers,
                     country_colour = "#1F1F1F", country_lwd = 0.4,
                     province_colour = "#888888", province_lwd = 0.15,
                     province_alpha = 0.5,
                     with_ten_dash = TRUE, ten_dash_colour = "#1F1F1F",
                     ten_dash_lwd = 0.45) +
    scale_fill_v2_sequential(name = fill_label, palette = "lajolla", direction = 1,
                             limits = c(0, fill_max), oob = scales::squish,
                             na.value = "white") +
    v2_china_coord() + theme_v2_map(11) +
    labs(title = sprintf("Survey effort — %d (%s)", yr, resolution_label),
         subtitle = sprintf("%s per %s grid (Combined)", fill_label, resolution_label))

  out <- file.path(out_dir, sprintf("frame_%s_%04d.png", stem, yr))
  ggsave(out, p, width = 9.4, height = 6.6, dpi = 200, bg = "white")
  out
}

# 通用渲染函数：分源版
render_bysrc <- function(yr, grid_sf, effort_dt, metric_col, fill_label,
                         fill_max, out_dir, stem, resolution_label) {
  sub <- effort_dt[year == yr]
  if (nrow(sub) == 0) return(NULL)

  sf_yr <- grid_sf |> left_join(as_tibble(sub), by = "grid_cell")
  fv <- rlang::sym(metric_col)

  grid_colour <- if (resolution_label == "100 km")
    scales::alpha("white", 0.06) else NA
  grid_lwd <- if (resolution_label == "100 km") 0.03 else 0

  p <- ggplot(sf_yr) +
    geom_sf(aes(fill = !!fv),
            colour = grid_colour, linewidth = grid_lwd) +
    china_map_layers(china_layers,
                     country_colour = "#1F1F1F", country_lwd = 0.35,
                     province_colour = "#888888", province_lwd = 0.12,
                     province_alpha = 0.45,
                     with_ten_dash = TRUE, ten_dash_colour = "#1F1F1F",
                     ten_dash_lwd = 0.4) +
    scale_fill_v2_sequential(name = fill_label, palette = "lajolla", direction = 1,
                             limits = c(0, fill_max), oob = scales::squish,
                             na.value = "white") +
    facet_wrap(~ source_short, ncol = 3) +
    v2_china_coord() + theme_v2_map(9) +
    labs(title = sprintf("Survey effort by source — %d (%s)", yr, resolution_label),
         subtitle = sprintf("%s per %s grid", fill_label, resolution_label))

  out <- file.path(out_dir, sprintf("frame_bysrc_%s_%04d.png", stem, yr))
  ggsave(out, p, width = 16, height = 6, dpi = 200, bg = "white")
  out
}

## =============================================================================
## 3. 100km GIF
## =============================================================================

message("[rerender] === 100 km GIFs ===")

anim_dir_100 <- file.path(P$figures_v2, "anim_100km")
dir.create(anim_dir_100, recursive = TRUE, showWarnings = FALSE)

metrics <- list(
  list(col = "n_visits",       label = "Visit events",    stem = "visits"),
  list(col = "n_observers",    label = "Unique observers", stem = "observers"),
  list(col = "n_species",      label = "Species detected", stem = "species"),
  list(col = "n_birding_days", label = "Birding days",    stem = "days")
)

for (m in metrics) {
  fill_max <- vmax_100[[m$col]]

  message(sprintf("  100km Combined: %s", m$stem))
  frames <- unlist(lapply(years, function(yr)
    render_combined(yr, grid_100, grid_year_all, m$col, m$label,
                    fill_max, anim_dir_100, m$stem, "100 km")))
  gif_path <- file.path(P$figures_v2,
                        sprintf("fig_effort_annual_100km_%s.gif", m$stem))
  gifski::gifski(frames, gif_path, width = 1880, height = 1320,
                 delay = 0.8, loop = TRUE)

  message(sprintf("  100km by-source: %s", m$stem))
  frames_src <- unlist(lapply(years, function(yr)
    render_bysrc(yr, grid_100, grid_year_all, m$col, m$label,
                 fill_max, anim_dir_100, m$stem, "100 km")))
  gif_src_path <- file.path(P$figures_v2,
                            sprintf("fig_effort_annual_100km_%s_by_source.gif", m$stem))
  gifski::gifski(frames_src, gif_src_path, width = 3200, height = 1200,
                 delay = 0.8, loop = TRUE)
}

## =============================================================================
## 4. 10km GIF
## =============================================================================

message("[rerender] === 10 km GIFs ===")

anim_dir_10 <- file.path(P$figures_v2, "anim_10km")
dir.create(anim_dir_10, recursive = TRUE, showWarnings = FALSE)

for (m in metrics) {
  fill_max <- vmax_10[[m$col]]

  message(sprintf("  10km Combined: %s", m$stem))
  frames <- unlist(lapply(years, function(yr)
    render_combined(yr, grid_10, grid10_year_all, m$col, m$label,
                    fill_max, anim_dir_10, m$stem, "10 km")))
  gif_path <- file.path(P$figures_v2,
                        sprintf("fig_effort_annual_10km_%s.gif", m$stem))
  gifski::gifski(frames, gif_path, width = 1880, height = 1320,
                 delay = 0.8, loop = TRUE)

  message(sprintf("  10km by-source: %s", m$stem))
  frames_src <- unlist(lapply(years, function(yr)
    render_bysrc(yr, grid_10, grid10_year_all, m$col, m$label,
                 fill_max, anim_dir_10, m$stem, "10 km")))
  gif_src_path <- file.path(P$figures_v2,
                            sprintf("fig_effort_annual_10km_%s_by_source.gif", m$stem))
  gifski::gifski(frames_src, gif_src_path, width = 3200, height = 1200,
                 delay = 0.8, loop = TRUE)
}

message("[rerender] === All done ===")
