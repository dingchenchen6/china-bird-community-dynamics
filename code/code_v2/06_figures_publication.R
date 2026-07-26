#!/usr/bin/env Rscript
## 06_figures_publication.R
##
## 阶段 6：把 stage-4/5 输出的所有结果整合为出版级图集。
## 全部地图遵守 v2 地图规则：coord_sf 强制 bbox、不画十段线、无 NA 子图、
## 跨指标 z-score 对齐、专业字体与调色。
##
## 输出（figures_v2/）：
##   A. 群落系数 caterpillar（已在 stage-4 出，这里复制为标准命名）
##   B. 群落多样性时间切片（CRI 带宽副图）
##   C. 群落动态趋势 4-panel + CRI 不确定度地图
##   D. Bray + Sørensen + Baselga turnover/nestedness 双面板地图
##   E. naive vs corrected richness 对比
##   F. 物种 occupancy / 探测概率轨迹图
##   G. 多样性时序轨迹（带 95% CRI 阴影）

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(sf)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_spatial.R"))
source(file.path(CODE_V2, "utils_mapping.R"))

P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_pilot_60sp_ar1")

message(sprintf("[stage-6] assembling publication figures for %s", RUN_LABEL))

china_layers <- load_china_layers(P$china_boundary, P$province_line)
grid_sf <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)

## --- B. 时间切片多样性地图（多指标 × 5 主期） ----------------------------

metrics_path <- v2_file("results",
                         paste0("table_community_metrics_with_cri_", RUN_LABEL))
if (file.exists(metrics_path)) {
  m_long <- read_csv(metrics_path, show_col_types = FALSE)
  show_metrics <- c("corrected_richness", "shannon", "pd_prob", "mpd_prob",
                    "trait_volume", "rao_q")
  m_show <- m_long |>
    filter(metric %in% show_metrics) |>
    mutate(
      metric_label = recode(metric,
        corrected_richness = "Taxonomic richness",
        shannon            = "Taxonomic Shannon",
        pd_prob            = "Faith's PD (prob-weighted)",
        mpd_prob           = "Phylogenetic MPD",
        trait_volume       = "Functional trait volume",
        rao_q              = "Functional Rao's Q"),
      metric_label = factor(metric_label, levels = c(
        "Taxonomic richness", "Taxonomic Shannon",
        "Faith's PD (prob-weighted)", "Phylogenetic MPD",
        "Functional trait volume", "Functional Rao's Q")),
      block_label = factor(block_label, levels = primary_blocks$block_label)
    ) |>
    filter(!is.na(metric_label), !is.na(block_label)) |>
    group_by(metric_label) |>
    mutate(value_z = as.numeric(scale(value_mean)),
           value_z = pmin(pmax(value_z, -2.5), 2.5)) |>
    ungroup()
  m_sf <- grid_sf |> inner_join(m_show, by = "grid_cell")

  ts_plot <- ggplot(m_sf) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_diverging(name = "Within-metric z-score (clipped at +/-2.5)") +
    facet_grid(metric_label ~ block_label) +
    v2_china_coord() +
    theme_v2_map(9.4) +
    labs(
      title = "Occupancy-corrected multidiversity time-slice maps",
      subtitle = sprintf("Posterior mean per 100 km grid | run=%s", RUN_LABEL),
      caption = "Each row is z-scaled within metric; central-white = near-mean grid value."
    )
  save_dual(ts_plot,
            paste0("fig_v2_multidiversity_timeslices_", RUN_LABEL),
            width = 14, height = 12)

  ## --- B'. 不确定度（CRI 宽度）地图 ---------------------------------------
  m_sd <- m_show |>
    mutate(cri_width = value_u95 - value_l95) |>
    group_by(metric_label) |>
    mutate(width_z = as.numeric(scale(cri_width)),
           width_z = pmin(pmax(width_z, -2.5), 2.5)) |>
    ungroup()
  m_sd_sf <- grid_sf |> inner_join(m_sd, by = "grid_cell")
  uncertainty_plot <- ggplot(m_sd_sf) +
    geom_sf(aes(fill = cri_width),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_sequential(name = "95% CRI width",
                             palette = "lajolla", direction = 1) +
    facet_grid(metric_label ~ block_label, scales = "fixed") +
    v2_china_coord() +
    theme_v2_map(9.4) +
    labs(title = "Posterior uncertainty (95% CRI width) of community metrics",
         subtitle = sprintf("run=%s | wider = more uncertain", RUN_LABEL))
  save_dual(uncertainty_plot,
            paste0("fig_v2_multidiversity_uncertainty_", RUN_LABEL),
            width = 14, height = 12)
}

## --- C. 群落动态趋势 4-panel + 对应不确定度 -------------------------------

trend_path <- v2_file("results",
                       paste0("table_grid_trends_with_cri_", RUN_LABEL))
if (file.exists(trend_path)) {
  trd <- read_csv(trend_path, show_col_types = FALSE)
  panel_metrics <- c("corrected_richness", "shannon", "trait_volume", "pd_prob")
  trd_show <- trd |>
    filter(metric %in% panel_metrics) |>
    mutate(metric_label = recode(metric,
      corrected_richness = "Richness trend",
      shannon            = "Shannon trend",
      trait_volume       = "Trait-volume trend",
      pd_prob            = "Faith's PD trend")) |>
    group_by(metric_label) |>
    mutate(z = as.numeric(scale(trend_mean)),
           z = pmin(pmax(z, -2.5), 2.5)) |>
    ungroup()
  trd_sf <- grid_sf |> inner_join(trd_show, by = "grid_cell")

  trd_plot <- ggplot(trd_sf) +
    geom_sf(aes(fill = z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = trd_sf |> filter(sig95),
            aes(geometry = geometry),
            inherit.aes = FALSE, fill = NA,
            colour = "#222222", linewidth = 0.10) +
    china_map_layers(china_layers) +
    scale_fill_v2_diverging(name = "Within-metric z-score (95% CRI excluding 0 outlined)") +
    facet_wrap(~ metric_label, ncol = 2) +
    v2_china_coord() +
    theme_v2_map(11) +
    labs(title = "Community dynamic trends with posterior uncertainty",
         subtitle = sprintf("run=%s | per-grid linear slope across 5 periods, posterior mean (z-scaled)", RUN_LABEL))
  save_dual(trd_plot,
            paste0("fig_v2_community_trends_with_cri_", RUN_LABEL),
            width = 11, height = 9.2)
}

## --- D. Temporal beta：Baselga 双面板地图 ---------------------------------

beta_path <- v2_file("results",
                      paste0("table_temporal_beta_with_cri_", RUN_LABEL))
if (file.exists(beta_path)) {
  bt <- read_csv(beta_path, show_col_types = FALSE)
  bt_grid <- bt |>
    filter(metric %in% c("turnover", "nestedness", "bray")) |>
    group_by(grid_cell, metric) |>
    summarise(value = mean(value_mean, na.rm = TRUE), .groups = "drop") |>
    mutate(metric_label = recode(metric,
      turnover   = "Baselga turnover",
      nestedness = "Baselga nestedness",
      bray       = "Bray-Curtis"))
  bt_sf <- grid_sf |> inner_join(bt_grid, by = "grid_cell")
  bt_plot <- ggplot(bt_sf) +
    geom_sf(aes(fill = value),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers) +
    scale_fill_v2_sequential(name = "Mean dissimilarity",
                             palette = "lajolla", direction = 1) +
    facet_wrap(~ metric_label, ncol = 3) +
    v2_china_coord() +
    theme_v2_map(11) +
    labs(title = "Temporal beta diversity decomposition",
         subtitle = sprintf("Baselga partition + Bray, mean across period pairs | run=%s",
                            RUN_LABEL),
         caption = "Bray = total community change; turnover = species replacement; nestedness = richness difference component.")
  save_dual(bt_plot,
            paste0("fig_v2_temporal_beta_baselga_", RUN_LABEL),
            width = 14, height = 5.6)
}

## --- E. naive vs occupancy-corrected richness ----------------------------

naive_grid <- read_csv(v2_file("results", "table_species_detection_coverage"),
                        show_col_types = FALSE)
# 用 visit_effort + species_visit 直接算 naive richness per grid×block
sv <- readRDS(file.path(P$derived_v2, "species_visit_2000_2024.rds"))
naive_block <- as_tibble(sv) |>
  group_by(grid_cell, block_id) |>
  summarise(naive_richness = n_distinct(species), .groups = "drop") |>
  left_join(primary_blocks, by = "block_id")

if (file.exists(metrics_path)) {
  corr_block <- m_long |>
    filter(metric == "corrected_richness") |>
    select(grid_cell, block_label, corrected_richness = value_mean)
  cmp <- naive_block |>
    inner_join(corr_block, by = c("grid_cell", "block_label"))
  cmp_plot <- ggplot(cmp, aes(naive_richness, corrected_richness,
                                 colour = block_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(alpha = 0.55, size = 1.4) +
    scale_colour_manual(values = V2_PALETTES$qualitative[1:5],
                         name = "Period") +
    facet_wrap(~ block_label, ncol = 5) +
    labs(title = "Naive vs occupancy-corrected richness across primary periods",
         subtitle = sprintf("Each point = one 100 km grid | run=%s", RUN_LABEL),
         x = "Naive species richness (observed)",
         y = "Occupancy-corrected richness (sum of psi)") +
    theme_v2_pub(11) +
    theme(legend.position = "none")
  save_dual(cmp_plot,
            paste0("fig_v2_naive_vs_corrected_richness_", RUN_LABEL),
            width = 13, height = 4.2)
}

## --- G. 多样性时序轨迹（CRI 阴影） ---------------------------------------

if (file.exists(metrics_path)) {
  traj <- m_long |>
    filter(metric %in% c("corrected_richness", "shannon",
                          "pd_prob", "trait_volume", "rao_q")) |>
    group_by(metric, block_label) |>
    summarise(
      mean_med = median(value_mean, na.rm = TRUE),
      l95 = quantile(value_mean, 0.025, na.rm = TRUE),
      u95 = quantile(value_mean, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(block_label = factor(block_label, levels = primary_blocks$block_label),
           metric_label = recode(metric,
             corrected_richness = "Taxonomic richness",
             shannon            = "Shannon",
             pd_prob            = "Faith's PD",
             trait_volume       = "Trait volume",
             rao_q              = "Rao's Q"))
  traj_plot <- ggplot(traj, aes(block_label, mean_med, group = metric_label)) +
    geom_ribbon(aes(ymin = l95, ymax = u95, fill = metric_label),
                alpha = 0.20, colour = NA) +
    geom_line(aes(colour = metric_label), linewidth = 0.7) +
    geom_point(aes(colour = metric_label), size = 2.2) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = V2_PALETTES$qualitative[1:5], guide = "none") +
    scale_colour_manual(values = V2_PALETTES$qualitative[1:5], guide = "none") +
    labs(title = "National multidiversity trajectories across 5-year periods",
         subtitle = sprintf("Median across grids; band = 2.5%%-97.5%% across-grid quantiles | run=%s",
                            RUN_LABEL),
         x = NULL, y = NULL) +
    theme_v2_pub(11) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))
  save_dual(traj_plot,
            paste0("fig_v2_multidiversity_trajectories_", RUN_LABEL),
            width = 12, height = 6.4)
}

message("[stage-6] Publication figures assembled.")
