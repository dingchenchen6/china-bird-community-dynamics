#!/usr/bin/env Rscript
## 00_remake_maps_from_v1_results.R
##
## 阶段 0：用 v1 已有的 result CSV，按 v2 出图工具重画所有地图。
## 解决两件事：
##   1) v1 地图"南北截断 + 鹰眼图错觉"——根因是没显式 coord_sf + 画了十段线把 bbox 拉太大。
##   2) 统一字体 / 调色 / 主题，做出"专业美观"基线。
##
## 输入：results/ 下 v1 表（top200sp_ar1）
## 输出：figures_v2/*_v2.{png,pdf}
##
## 注意：本脚本不重跑模型、不动 v1 文件。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(patchwork)
  library(forcats)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_spatial.R"))
source(file.path(CODE_V2, "utils_mapping.R"))

P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "top200sp_ar1")

message(sprintf("[stage-0] Remaking maps from v1 run=%s", RUN_LABEL))

## --- 1. 载入空间 base 层 + v1 网格 -----------------------------------------

china_layers <- load_china_layers(P$china_boundary, P$province_line)
grid_cache_path <- file.path(P$derived_v2, sprintf("china_grid_100km_v2.rds"))
grid_sf <- load_or_build_v2_grid(P$china_boundary, resolution_m = 1e5,
                                 cache_path = grid_cache_path)

## --- 2. 载入 v1 result 表 --------------------------------------------------

read_v1 <- function(stem) {
  path <- file.path(P$results_v1, sprintf("%s_%s.csv", stem, RUN_LABEL))
  if (!file.exists(path)) {
    message("  [skip] missing v1 file: ", basename(path))
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}

trends_tbl     <- read_v1("table_multispecies_dynamic_grid_metric_trends")
community_tbl  <- read_v1("table_multispecies_dynamic_community_metrics_by_grid_block")
beta_tbl       <- read_v1("table_multispecies_dynamic_temporal_beta_by_grid")
occ_grid_tbl   <- read_v1("table_multispecies_dynamic_community_occupancy_by_grid_block")
det_grid_tbl   <- read_v1("table_multispecies_dynamic_community_detection_by_grid_block")
blocks_tbl     <- read_csv(file.path(P$results_v1, "table_primary_5year_blocks.csv"),
                           show_col_types = FALSE)

## --- 3. 网格趋势地图（4 panel：richness / phylo MPD / trait volume / temporal beta） --

if (!is.null(trends_tbl) && !is.null(beta_tbl)) {
  beta_grid_summary <- beta_tbl |>
    group_by(grid_cell) |>
    summarise(mean_temporal_beta = mean(temporal_beta_bray, na.rm = TRUE),
              .groups = "drop")

  metric_levels <- c("Richness trend", "Phylogenetic MPD trend",
                     "Trait-volume trend", "Temporal beta (mean)")

  ## --- 4a. 趋势：逐指标 z-score + 同一 diverging 色盘 -----------------------
  trend_long <- trends_tbl |>
    select(grid_cell, richness_trend, phylo_mpd_trend, trait_volume_trend) |>
    left_join(beta_grid_summary, by = "grid_cell") |>
    pivot_longer(-grid_cell, names_to = "metric_raw", values_to = "value") |>
    mutate(metric = recode(metric_raw,
      richness_trend     = "Richness trend",
      phylo_mpd_trend    = "Phylogenetic MPD trend",
      trait_volume_trend = "Trait-volume trend",
      mean_temporal_beta = "Temporal beta (mean)"
    )) |>
    filter(!is.na(metric)) |>
    group_by(metric) |>
    mutate(value_z = as.numeric(scale(value)),
           value_z = pmin(pmax(value_z, -2.5), 2.5)) |>
    ungroup() |>
    mutate(metric = factor(metric, levels = metric_levels))

  trend_sf <- grid_sf |> inner_join(trend_long, by = "grid_cell") |>
    filter(!is.na(metric))

  trend_map <- ggplot(trend_sf) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.10), linewidth = 0.05) +
    china_map_layers(china_layers) +
    scale_fill_v2_diverging(
      name = "Within-metric z-score (clipped at +/-2.5)",
      limits = c(-2.5, 2.5)
    ) +
    facet_wrap(~ metric, ncol = 2, drop = TRUE) +
    v2_china_coord() +
    theme_v2_map(11) +
    labs(
      title = "Community dynamic trends across China (v2 cleaned)",
      subtitle = sprintf("Top-200-species multi-species dynamic occupancy run | label=%s", RUN_LABEL),
      caption = "Each panel z-scaled within metric so cross-panel intensity is comparable; central white = near-zero trend."
    )
  save_dual(trend_map, sprintf("fig_community_trends_4panel_%s_v2", RUN_LABEL),
            width = 11, height = 9.2)
  message("  [ok] community 4-panel trend map written (z-scored, NA dropped)")

  ## --- 4b. Temporal beta 单图（非负量，sequential 色盘） -------------------
  beta_sf <- grid_sf |> inner_join(beta_grid_summary, by = "grid_cell")
  beta_lim <- quantile(beta_sf$mean_temporal_beta, c(0.02, 0.98), na.rm = TRUE)
  beta_map <- ggplot(beta_sf) +
    geom_sf(aes(fill = mean_temporal_beta),
            colour = scales::alpha("white", 0.10), linewidth = 0.05) +
    china_map_layers(china_layers) +
    scale_fill_v2_sequential(
      name = "Mean Bray temporal beta",
      palette = "lajolla", direction = 1,
      limits = unname(beta_lim)
    ) +
    v2_china_coord() +
    theme_v2_map(11) +
    labs(
      title = "Mean temporal turnover (Bray) across primary periods",
      subtitle = sprintf("Top-200-species multi-species dynamic occupancy | label=%s", RUN_LABEL),
      caption = "Color clipped to robust 2-98% range; higher = stronger period-to-period community change."
    )
  save_dual(beta_map, sprintf("fig_temporal_beta_solo_%s_v2", RUN_LABEL),
            width = 7.4, height = 5.6)
  message("  [ok] temporal-beta solo map written")
}

## --- 4. 五时间切片多样性面板 ----------------------------------------------

if (!is.null(community_tbl)) {
  # 包含 v1 实际有的 6 个多样性维度：分类 / 系统发育 / 功能
  metric_map <- c(
    corrected_richness        = "Taxonomic richness",
    taxon_shannon             = "Taxonomic Shannon",
    phylo_pd                  = "Faith's PD",
    phylo_mpd_weighted        = "Phylogenetic MPD (weighted)",
    life_history_trait_volume = "Functional trait volume",
    trait_rao_q               = "Functional Rao's Q"
  )
  metric_levels <- unname(metric_map)
  show_cols <- intersect(names(metric_map), names(community_tbl))

  community_long <- community_tbl |>
    select(grid_cell, block_label, all_of(show_cols)) |>
    pivot_longer(-c(grid_cell, block_label), names_to = "metric_raw",
                 values_to = "value") |>
    mutate(
      metric = factor(metric_map[metric_raw], levels = metric_levels),
      block_label = factor(block_label, levels = blocks_tbl$block_label)
    ) |>
    filter(!is.na(block_label), !is.na(metric)) |>
    group_by(metric) |>
    mutate(value_z = as.numeric(scale(value)),
           value_z = pmin(pmax(value_z, -2.5), 2.5)) |>
    ungroup()
  community_sf <- grid_sf |>
    inner_join(community_long, by = "grid_cell")

  ts_map <- ggplot(community_sf) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_diverging(name = "Within-metric z-score (clipped at +/-2.5)") +
    facet_grid(metric ~ block_label) +
    v2_china_coord() +
    theme_v2_map(9.4) +
    labs(
      title = "Multidiversity time-slice maps (taxonomic + phylogenetic + functional)",
      subtitle = sprintf("Six diversity dimensions across five 5-year periods | run=%s",
                         RUN_LABEL),
      caption = paste0("Each row z-scaled within its metric so cross-panel intensity is comparable. ",
                       "Note: low-variance metrics (e.g. Rao's Q) get amplified — small absolute changes look strong.")
    )
  save_dual(ts_map, sprintf("fig_multidiversity_timeslices_%s_v2", RUN_LABEL),
            width = 14, height = 12)
  message(sprintf("  [ok] multidiversity time-slice map written (%d metrics)",
                  length(show_cols)))
}

## --- 5. 占域 / 探测概率时空图（如可用） ------------------------------------

if (!is.null(occ_grid_tbl)) {
  occ_long <- occ_grid_tbl |>
    transmute(grid_cell, block_label,
              metric = "Community occupancy",
              value  = mean_site_species_occupancy)
  if (!is.null(det_grid_tbl)) {
    det_long <- det_grid_tbl |>
      transmute(grid_cell, block_label,
                metric = "Community detection",
                value  = detection_prob_comm_mean)
    occ_long <- bind_rows(occ_long, det_long)
  }
  occ_long <- occ_long |>
    mutate(block_label = factor(block_label, levels = blocks_tbl$block_label)) |>
    filter(!is.na(block_label))
  occ_sf <- grid_sf |> inner_join(occ_long, by = "grid_cell")

  occ_map <- ggplot(occ_sf) +
    geom_sf(aes(fill = value),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_sequential(name = "Probability", limits = c(0, 1),
                             palette = "lajolla", direction = -1) +
    facet_grid(metric ~ block_label) +
    v2_china_coord() +
    theme_v2_map(9.4) +
    labs(
      title = "Community occupancy and detection across primary periods (v2 cleaned)",
      subtitle = sprintf("Posterior mean per 100 km grid | run=%s", RUN_LABEL)
    )
  save_dual(occ_map, sprintf("fig_occ_det_timeslices_%s_v2", RUN_LABEL),
            width = 14, height = 6)
  message("  [ok] occupancy/detection time-slice map written")
}

message("[stage-0] All v2 maps written to figures_v2/.")
