#!/usr/bin/env Rscript
## 15_annual_effort_animation.R
##
## 目标：
##   1. 输出不同数据源 × 不同调查努力指标的年度统计表（2000-2024 每年）
##   2. 生成 100km 分辨率年度动态地图（GIF）
##   3. 生成 10km 分辨率年度动态地图（GIF）
##
## 输入：data/derived_v2/combined_events_merged_dedup_2000_2025.rds
##       data/derived_v2/china_grid_100km_v2.rds
## 输出：results_v2/table_effort_metrics_year_source_detailed.csv
##       results_v2/table_effort_by_grid_year_source_100km.csv
##       results_v2/table_effort_by_grid_year_source_10km.csv
##       figures_v2/anim_effort_*_{100km,10km}_*.gif

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(sf)
  library(ggplot2)
  library(gganimate)
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

#' 安全保存GIF：先写出PNG帧序列，再调用gifski合成
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

#' 生成统一色标的网格动态地图（返回 ggplot 列表）
make_annual_grid_maps <- function(eff_dt, grid_sf, metric_col = "n_events",
                                  metric_label = "log10(events+1)",
                                  by_source = FALSE,
                                  global_max = NULL,
                                  mainland_only = FALSE) {
  ## 与网格合并
  eff_dt <- as.data.table(eff_dt)
  setnames(eff_dt, metric_col, "metric_val", skip_absent = TRUE)

  ## 计算全局最大值（统一色标）
  if (is.null(global_max)) {
    global_max <- log10(max(eff_dt$metric_val, na.rm = TRUE) + 1)
  }

  years <- sort(unique(eff_dt$year))
  sources <- if (by_source) {
    c("China_Birdwatch", "eBird_GBIF", "Combined")
  } else {
    "Combined"
  }

  ## 准备空网格（所有年份/源都显示完整网格底图）
  grid_all <- grid_sf |>
    select(grid_cell) |>
    mutate(fill_val = NA_real_)

  coord_args <- if (mainland_only) {
    list(xlim = V2_MAINLAND_BBOX$xlim, ylim = V2_MAINLAND_BBOX$ylim)
  } else {
    list(xlim = V2_CHINA_BBOX$xlim, ylim = V2_CHINA_BBOX$ylim)
  }

  plot_list <- list()
  idx <- 1

  for (yr in years) {
    for (src in sources) {
      yr_dt <- eff_dt[year == yr & source_short == src, .(grid_cell, metric_val)]
      map_dt <- grid_all |>
        left_join(yr_dt, by = "grid_cell") |>
        mutate(fill_val = log10(metric_val + 1))

      src_label <- switch(src,
        "China_Birdwatch" = "China Birdwatch",
        "eBird_GBIF"      = "eBird/GBIF",
        "Combined"        = "Combined")

      p <- ggplot(map_dt) +
        geom_sf(aes(fill = fill_val),
                colour = alpha("white", 0.08), linewidth = 0.04) +
        china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
        scale_fill_v2_sequential(
          name = metric_label,
          palette = "lajolla", direction = 1,
          limits = c(0, global_max),
          oob = squish,
          na.value = "grey92"
        ) +
        do.call(v2_china_coord, coord_args) +
        theme_v2_map(10.5) +
        labs(
          title = sprintf("Survey effort — %s", src_label),
          subtitle = sprintf("Year: %d  |  %s", yr, metric_label),
          caption = sprintf("Cells with no data shown in grey. Max %s = %.2f",
                            metric_label, global_max)
        )

      if (by_source) {
        p <- p + facet_wrap(~ "", ncol = 1)  # 占位，后续按源分面
      }

      plot_list[[idx]] <- p
      names(plot_list)[idx] <- sprintf("%d_%s", yr, src)
      idx <- idx + 1
    }
  }
  plot_list
}

## --- 1. 读取事件数据 --------------------------------------------------------

message("[stage-15] Loading dedup'd events (this may take a minute)...")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                                "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]
message(sprintf("  events %d-%d after filter: %s",
                ANALYSIS_YR_LO, ANALYSIS_YR_HI,
                format(nrow(events), big.mark = ",")))

## 数据源分类（与 02b 保持一致）
events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]
events[, day_key := sprintf("%04d-%02d-%02d", year, month, day)]
events[, observer_key := tolower(coalesce(username, "anon"))]
events[, visit_key := paste(observer_key, day_key, sep = "|")]
events[, has_dur := !is.na(duration_min)]

## --- 2. 年度统计表（增强版）--------------------------------------------------

message("[stage-15] Building annual effort metrics table by source...")

annual_per_source <- events[, .(
  n_events       = .N,
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_birding_days = uniqueN(day_key),
  n_grids_round1 = uniqueN(paste(round(longitude, 1), round(latitude, 1))),
  n_species      = uniqueN(species),
  total_duration = sum(duration_min[has_dur], na.rm = TRUE),
  mean_duration  = mean(duration_min[has_dur], na.rm = TRUE),
  median_duration = median(duration_min[has_dur], na.rm = TRUE)
), by = .(year, source_short)]

annual_combined <- events[, .(
  n_events       = .N,
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_birding_days = uniqueN(day_key),
  n_grids_round1 = uniqueN(paste(round(longitude, 1), round(latitude, 1))),
  n_species      = uniqueN(species),
  total_duration = sum(duration_min[has_dur], na.rm = TRUE),
  mean_duration  = mean(duration_min[has_dur], na.rm = TRUE),
  median_duration = median(duration_min[has_dur], na.rm = TRUE)
), by = year][, source_short := "Combined"]

annual_all <- rbindlist(list(annual_per_source, annual_combined), use.names = TRUE)
annual_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
setorder(annual_all, year, source_short)

write_csv(as_tibble(annual_all),
          v2_file("results", "table_effort_metrics_year_source_detailed"))
message("  Saved: results_v2/table_effort_metrics_year_source_detailed.csv")

## --- 3. 100km 网格处理 -------------------------------------------------------

message("[stage-15] Processing 100km grid...")

grid_100 <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
message(sprintf("  100km grid cells: %d", nrow(grid_100)))

## 若事件已有 grid_cell 且非全NA，复用；否则重新分配
if (!"grid_cell" %in% names(events) || all(is.na(events$grid_cell))) {
  message("  assigning events to 100km grid...")
  old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  pts <- st_as_sf(events[, .(longitude, latitude)],
                  coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  joined <- st_join(pts, grid_100[, "grid_cell"],
                    join = st_intersects, left = TRUE)
  events[, grid_cell := joined$grid_cell]
}
events_100 <- events[!is.na(grid_cell)]
message(sprintf("  events with 100km grid_cell: %s",
                format(nrow(events_100), big.mark = ",")))

## 按 (grid × year × source) 汇总
eff_100_source <- events_100[, .(
  n_events    = uniqueN(visit_key),
  n_observers = uniqueN(observer_key),
  n_dates     = uniqueN(day_key),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE)
), by = .(grid_cell, year, source_short)]

eff_100_combined <- events_100[, .(
  n_events    = uniqueN(visit_key),
  n_observers = uniqueN(observer_key),
  n_dates     = uniqueN(day_key),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE)
), by = .(grid_cell, year)][, source_short := "Combined"]

eff_100_all <- rbindlist(list(eff_100_source, eff_100_combined), use.names = TRUE)
eff_100_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
setorder(eff_100_all, grid_cell, year, source_short)

write_csv(as_tibble(eff_100_all),
          v2_file("results", "table_effort_by_grid_year_source_100km"))
message("  Saved: results_v2/table_effort_by_grid_year_source_100km.csv")

## --- 4. 100km 动态地图（GIF）-------------------------------------------------

message("[stage-15] Rendering 100km annual animated maps...")

china_layers <- load_china_layers(P$china_boundary, P$province_line)

## 4a. Combined 单一源动态地图
message("  -> Combined n_events animation (100km)...")
global_max_100_combined <- log10(max(eff_100_all[source_short == "Combined"]$n_events, na.rm = TRUE) + 1)

years <- sort(unique(eff_100_all$year))
plots_100_combined <- list()
for (i in seq_along(years)) {
  yr <- years[i]
  yr_dt <- eff_100_all[year == yr & source_short == "Combined", .(grid_cell, n_events)]
  map_dt <- grid_100 |>
    select(grid_cell) |>
    left_join(yr_dt, by = "grid_cell") |>
    mutate(fill_val = log10(n_events + 1))

  plots_100_combined[[i]] <- ggplot(map_dt) +
    geom_sf(aes(fill = fill_val),
            colour = alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_sequential(
      name = "log10(events+1)", palette = "lajolla", direction = 1,
      limits = c(0, global_max_100_combined), oob = squish, na.value = "grey92"
    ) +
    v2_china_coord() +
    theme_v2_map(11) +
    labs(
      title = "Annual survey effort — Combined (merged + deduplicated)",
      subtitle = sprintf("Year: %d  |  100 km grid  |  Events per grid cell", yr),
      caption = sprintf("Max log10(events+1) across all years = %.2f", global_max_100_combined)
    )
}

save_gif_from_plots(
  plots_100_combined,
  file.path(P$figures_v2, "anim_effort_events_100km_combined.gif"),
  width = 1200, height = 900, delay = 0.6
)

## 4b. 分源动态地图（三源并列）
message("  -> By-source n_events animation (100km)...")
global_max_100_src <- log10(max(eff_100_all$n_events, na.rm = TRUE) + 1)

plots_100_bysource <- list()
idx <- 1
for (yr in years) {
  yr_dt <- eff_100_all[year == yr, .(grid_cell, source_short, n_events)]
  map_dt <- grid_100 |>
    select(grid_cell) |>
    left_join(yr_dt, by = "grid_cell") |>
    mutate(fill_val = log10(n_events + 1))

  plots_100_bysource[[idx]] <- ggplot(map_dt) +
    geom_sf(aes(fill = fill_val),
            colour = alpha("white", 0.06), linewidth = 0.03) +
    china_map_layers(china_layers, country_lwd = 0.20, province_lwd = 0.10) +
    scale_fill_v2_sequential(
      name = "log10(events+1)", palette = "lajolla", direction = 1,
      limits = c(0, global_max_100_src), oob = squish, na.value = "grey92"
    ) +
    facet_wrap(~ source_short, ncol = 3) +
    v2_china_coord() +
    theme_v2_map(10) +
    labs(
      title = sprintf("Annual survey effort by data source — Year %d", yr),
      subtitle = "100 km grid  |  log10(events + 1)",
      caption = "Grey = no data."
    )
  idx <- idx + 1
}

save_gif_from_plots(
  plots_100_bysource,
  file.path(P$figures_v2, "anim_effort_events_100km_by_source.gif"),
  width = 1800, height = 700, delay = 0.6
)

## 清理内存
gc()

## --- 5. 10km 网格处理 --------------------------------------------------------

message("[stage-15] Building 10km grid (this may take a few minutes)...")

grid_10_cache <- file.path(P$derived_v2, "china_grid_10km_v2.rds")
if (file.exists(grid_10_cache)) {
  grid_10 <- readRDS(grid_10_cache)
  message(sprintf("  10km grid loaded from cache: %d cells", nrow(grid_10)))
} else {
  grid_10 <- load_or_build_v2_grid(P$china_boundary, resolution_m = 1e4,
                                   cache_path = grid_10_cache)
  message(sprintf("  10km grid built: %d cells", nrow(grid_10)))
}

## 分配事件到 10km 网格（数据量大，使用 sf planar）
message("[stage-15] Assigning events to 10km grid (may take several minutes)...")
old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
on.exit(sf::sf_use_s2(old_s2), add = TRUE)

## 分块处理以避免内存爆炸：每 50 万行一批
batch_size <- 500000L
events_10_list <- list()
n_events <- nrow(events)
n_batches <- ceiling(n_events / batch_size)

for (b in seq_len(n_batches)) {
  start_idx <- (b - 1) * batch_size + 1
  end_idx   <- min(b * batch_size, n_events)
  batch <- events[start_idx:end_idx]

  pts <- st_as_sf(batch[, .(longitude, latitude)],
                  coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  joined <- st_join(pts, grid_10[, "grid_cell"],
                    join = st_intersects, left = TRUE)
  batch[, grid_cell_10km := joined$grid_cell]
  events_10_list[[b]] <- batch[!is.na(grid_cell_10km)]

  message(sprintf("  batch %d/%d: rows %s-%s -> assigned %s",
                  b, n_batches, format(start_idx, big.mark = ","),
                  format(end_idx, big.mark = ","),
                  format(nrow(events_10_list[[b]]), big.mark = ",")))
  rm(batch, pts, joined)
  gc()
}

events_10 <- rbindlist(events_10_list, use.names = TRUE)
message(sprintf("  total events with 10km grid_cell: %s",
                format(nrow(events_10), big.mark = ",")))
rm(events_10_list)
gc()

## 按 (grid × year × source) 汇总
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

## 清理事件数据释放内存
rm(events, events_10)
gc()

## --- 6. 10km 动态地图（GIF）--------------------------------------------------

message("[stage-15] Rendering 10km annual animated maps...")

## 为了渲染效率，仅渲染 Combined 源（10km分源数据量巨大，渲染很慢）
## 若用户需要10km分源，可后续单独运行

message("  -> Combined n_events animation (10km)...")
global_max_10_combined <- log10(max(eff_10_all[source_short == "Combined"]$n_events, na.rm = TRUE) + 1)

## 10km 网格也仅需必要列
grid_10_lean <- grid_10 |> select(grid_cell)

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

## 10km 分源（可选，数据量大，可能很慢）—— 默认启用，若内存/时间不足可注释掉
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

message("[stage-15] All annual effort tables and animated maps completed.")
