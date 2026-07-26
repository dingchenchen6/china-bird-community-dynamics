#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把 stage-4/5/05b/05c 输出整合为 Nature/Science 级图集（覆盖
#   standard_occupancy_figures.md 的 A–J 10 组）
#
# Objective / 分析目标:
#   v4 比 v3 新增 / 修复：
#     - H. 物种 × 协变量 系数 heatmap（× 标 95% CRI 不跨 0）
#     - I. AR1 ρ 后验密度（由 05b 出，但这里补一张 forest）
#     - F+. CWM diet_specialization / habitat_breadth 空间图
#         （依赖 05 产出 table_cwm_spatial_pattern_*）
#     - K. 朴素 vs 校正趋势配对图（颜色标方向翻转）
#     - 地图全部用 utils_mapping_v4 的 v4_china_coord（固定 bbox）
#     - 缺值灰色保留；within-metric z-score
#
# Input data / 输入数据:
#   results_v4/table_community_metrics_with_cri_<run_label>.csv
#   results_v4/table_trend_summary_<run_label>.csv
#   results_v4/table_baselga_summary_<run_label>.csv
#   results_v4/table_baselga_global_<run_label>.csv
#   results_v4/table_naive_vs_corrected_<run_label>.csv
#   results_v4/table_species_trend_<run_label>.csv
#   results_v4/table_species_trend_classify_<run_label>.csv
#   results_v4/table_species_hotspot_<run_label>.csv
#   results_v4/table_cwm_spatial_pattern_<run_label>.csv
#   results_v4/table_beta_species_<run_label>.csv
#   results_v4/table_brms_driver_coefs_<run_label>.csv
#   data/derived_v4/china_grid_<size>km_v4.rds (or v3 fallback)
#
# Main workflow / 主要流程:
#   按 A–K 段落依次产出
#
# Key assumptions / 关键假设:
#   - 04→05→05b→05c 已完成
#   - grid_sf 与 grid_cell 编码一致
#
# Main packages / 主要包:
#   ggplot2, patchwork, sf, dplyr, tidyr, scales, scico（可选）
#
# Output directory / 输出路径:
#   figures_v4/fig_v4_*_<run_label>.{png,pdf}
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(forcats); library(ggplot2)
  library(patchwork); library(sf)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_mapping.R"))
source(file.path(CODE_V4, "utils_spatial.R"))
P <- ensure_v4_dirs()

is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
RUN_LABEL <- if (is_pilot) PILOT_LABEL else RUN_LABEL
log_time("06", sprintf("Assembling v4 publication figures: %s", RUN_LABEL))

# ── 共享 ─────────────────────────────────────────────────────────────
grid_sf <- safe_read(v4_file("derived",
                              paste0("china_grid_", GRID_SIZE_KM, "km_v4"), "rds"),
                      quiet = TRUE) %||%
            safe_read(v3_file("derived",
                              paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds"),
                      quiet = TRUE) %||%
            safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"),
                      quiet = TRUE)
if (is.null(grid_sf)) stop("[06] grid_sf not found")

china_layers <- china_province_basemap()

METRIC_LABELS <- c(
  corrected_richness = "Taxonomic richness",
  shannon            = "Shannon diversity",
  pd_prob            = "Faith's PD (prob-weighted)",
  mpd_prob           = "MPD (prob-weighted)",
  trait_volume       = "Functional trait volume",
  rao_q              = "Rao's Q",
  feve               = "Functional evenness",
  fdiv               = "Functional divergence",
  cwm_diet_specialization = "CWM diet specialization",
  cwm_habitat_breadth     = "CWM habitat breadth"
)

PANEL_METRICS <- c("corrected_richness", "shannon", "trait_volume", "pd_prob")

metrics_long <- read_csv_safe(v4_file("results",
                                       paste0("table_community_metrics_with_cri_", RUN_LABEL)))

# ── A. 多样性时间切片地图 ───────────────────────────────────────────
if (!is.null(metrics_long)) {
  show_metrics <- intersect(c("corrected_richness", "shannon",
                               "pd_prob", "mpd_prob",
                               "trait_volume", "rao_q", "feve", "fdiv"),
                             unique(metrics_long$metric))
  m_show <- metrics_long |>
    filter(metric %in% show_metrics) |>
    mutate(
      metric_label = recode(metric, !!!METRIC_LABELS),
      metric_label = factor(metric_label,
                            levels = METRIC_LABELS[show_metrics]),
      block_label  = factor(block_label,
                            levels = sort(unique(block_label)))
    ) |>
    filter(!is.na(metric_label), !is.na(block_label)) |>
    group_by(metric_label) |>
    mutate(value_z = pmin(pmax(as.numeric(scale(value_mean)), -2.5), 2.5)) |>
    ungroup()

  # 重要：left_join 保留所有网格，缺值 NA 显示为灰色（map_quality_rules.md #4）
  m_sf <- grid_sf |>
    left_join(m_show, by = "grid_cell") |>
    filter(!is.na(metric_label))

  ts_plot <- ggplot(m_sf) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "z-score", palette = "diverging",
                        na.value = "grey92") +
    facet_grid(metric_label ~ block_label) +
    v4_china_coord() +
    theme_nature_map() +
    theme(strip.text = element_text(size = 6))
  save_nature(ts_plot,
              paste0("fig_v4_multidiversity_timeslices_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 160)

  # A'. 不确定度（CRI 宽度）
  m_sd <- m_show |>
    mutate(cri_width = value_u95 - value_l95) |>
    group_by(metric_label) |>
    mutate(width_z = pmin(pmax(as.numeric(scale(cri_width)), -2.5), 2.5)) |>
    ungroup()
  m_sd_sf <- grid_sf |> left_join(m_sd, by = "grid_cell") |>
    filter(!is.na(metric_label))
  unc_plot <- ggplot(m_sd_sf) +
    geom_sf(aes(fill = cri_width),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "95% CRI width", palette = "accent",
                        na.value = "grey92") +
    facet_grid(metric_label ~ block_label) +
    v4_china_coord() + theme_nature_map() +
    theme(strip.text = element_text(size = 6))
  save_nature(unc_plot,
              paste0("fig_v4_multidiversity_uncertainty_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 160)
}

# ── B. 群落趋势 4-panel ─────────────────────────────────────────────
trd <- read_csv_safe(v4_file("results", paste0("table_trend_summary_", RUN_LABEL)))
if (!is.null(trd)) {
  trd <- trd |> filter(method == "theil_sen")
  if (!"sig95" %in% names(trd) && all(c("q025", "q975") %in% names(trd))) {
    trd$sig95 <- (trd$q025 > 0) | (trd$q975 < 0)
  }
  if (!"trend_mean" %in% names(trd) && "mean" %in% names(trd)) {
    trd$trend_mean <- trd$mean
  }
  trd_show <- trd |>
    filter(metric %in% PANEL_METRICS) |>
    mutate(metric_label = recode(metric,
                                  corrected_richness = "Richness trend",
                                  shannon            = "Shannon trend",
                                  trait_volume       = "Trait-volume trend",
                                  pd_prob            = "PD trend")) |>
    group_by(metric_label) |>
    mutate(z = pmin(pmax(as.numeric(scale(trend_mean)), -2.5), 2.5)) |>
    ungroup()
  trd_sf <- grid_sf |> left_join(trd_show, by = "grid_cell") |>
    filter(!is.na(metric_label))

  p_trd <- ggplot(trd_sf) +
    geom_sf(aes(fill = z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = trd_sf |> filter(sig95),
            aes(geometry = geometry),
            inherit.aes = FALSE, fill = NA,
            colour = "#222222", linewidth = 0.10) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "Trend z-score", palette = "diverging",
                        na.value = "grey92") +
    facet_wrap(~ metric_label, ncol = 2) +
    v4_china_coord() + theme_nature_map()
  save_nature(p_trd, paste0("fig_v4_community_trends_with_cri_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 140)
}

# ── C. Baselga 双面板 ───────────────────────────────────────────────
bt <- read_csv_safe(v4_file("results", paste0("table_baselga_summary_", RUN_LABEL)))
if (!is.null(bt)) {
  bt_g <- bt |>
    filter(metric %in% c("beta_sim", "beta_sne", "beta_sor")) |>
    mutate(metric = recode(metric,
                           beta_sim = "turnover",
                           beta_sne = "nestedness",
                           beta_sor = "soerensen")) |>
    group_by(grid_cell, metric) |>
    summarise(value = mean(mean, na.rm = TRUE), .groups = "drop") |>
    mutate(metric_label = recode(metric,
                                  turnover   = "Baselga turnover",
                                  nestedness = "Baselga nestedness",
                                  soerensen  = "Sørensen total"))
  bt_sf <- grid_sf |> left_join(bt_g, by = "grid_cell") |>
    filter(!is.na(metric_label))
  p_bt <- ggplot(bt_sf) +
    geom_sf(aes(fill = value),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "Mean dissimilarity",
                        palette = "sequential",
                        na.value = "grey92") +
    facet_wrap(~ metric_label, ncol = 3) +
    v4_china_coord() + theme_nature_map()
  save_nature(p_bt, paste0("fig_v4_temporal_beta_baselga_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 60)
}

# ── D. Naive vs Corrected richness trend 配对图 ─────────────────────
nvc <- read_csv_safe(v4_file("results", paste0("table_naive_vs_corrected_", RUN_LABEL)))
if (!is.null(nvc)) {
  p_nvc <- ggplot(nvc, aes(naive_trend, corrected_trend,
                            colour = direction_flipped)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2,
                colour = "grey60", linewidth = 0.2) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.15) +
    geom_vline(xintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.15) +
    geom_point(alpha = 0.6, size = 0.7) +
    scale_colour_manual(values = c("FALSE" = "grey40",
                                    "TRUE"  = "#B2182B"),
                        labels = c("FALSE" = "same direction",
                                   "TRUE"  = "direction flipped"),
                        name = NULL) +
    labs(x = "Naive richness trend (raw)",
         y = "Occupancy-corrected richness trend") +
    theme_nature_pub()
  save_nature(p_nvc, paste0("fig_v4_naive_vs_corrected_trend_", RUN_LABEL),
              width_mm = NATURE_WIDTH_M, height_mm = 80)
}

# ── E. 多样性时序轨迹（CRI 阴影） ────────────────────────────────────
if (!is.null(metrics_long)) {
  traj_metrics <- c("corrected_richness", "shannon", "pd_prob",
                    "trait_volume", "rao_q", "feve", "fdiv",
                    "cwm_diet_specialization", "cwm_habitat_breadth")
  traj <- metrics_long |>
    filter(metric %in% traj_metrics) |>
    group_by(metric, block_label) |>
    summarise(
      mean_med = median(value_mean, na.rm = TRUE),
      l95 = quantile(value_mean, 0.025, na.rm = TRUE),
      u95 = quantile(value_mean, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(metric_label = recode(metric, !!!METRIC_LABELS)) |>
    filter(!is.na(metric_label))

  p_traj <- ggplot(traj, aes(block_label, mean_med, group = metric_label)) +
    geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.15, fill = NATURE_ACCENT) +
    geom_line(colour = NATURE_ACCENT, linewidth = 0.4) +
    geom_point(colour = NATURE_ACCENT, size = 1.2) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = NULL) +
    theme_nature_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_nature(p_traj, paste0("fig_v4_multidiversity_trajectories_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 100)
}

# ── F. CWM 空间（含 diet_specialization + habitat_breadth） ─────────
cwm <- read_csv_safe(v4_file("results", paste0("table_cwm_spatial_pattern_", RUN_LABEL)))
if (!is.null(cwm)) {
  cwm_long <- cwm |>
    pivot_longer(matches("^cwm_"), names_to = "trait", values_to = "value") |>
    filter(trait %in% c("cwm_diet_specialization", "cwm_habitat_breadth")) |>
    mutate(trait_label = recode(trait,
                                 cwm_diet_specialization = "CWM diet specialization",
                                 cwm_habitat_breadth     = "CWM habitat breadth"))
  # 取最后 period
  if ("period" %in% names(cwm_long)) {
    last_p <- max(as.numeric(gsub("P", "", cwm_long$period)), na.rm = TRUE)
    cwm_long <- cwm_long |> filter(period == paste0("P", last_p))
  }
  cwm_sf <- grid_sf |> left_join(cwm_long, by = "grid_cell") |>
    filter(!is.na(trait_label))
  if (nrow(cwm_sf) > 0) {
    p_cwm <- ggplot(cwm_sf) +
      geom_sf(aes(fill = value),
              colour = scales::alpha("white", 0.08), linewidth = 0.04) +
      geom_sf(data = china_layers, fill = NA, colour = "grey50",
              linewidth = NATURE_LINE_MAP) +
      scale_fill_nature_c(name = "CWM (last period)",
                          palette = "accent",
                          na.value = "grey92") +
      facet_wrap(~ trait_label, ncol = 2) +
      v4_china_coord() + theme_nature_map()
    save_nature(p_cwm, paste0("fig_v4_cwm_spatial_traits_", RUN_LABEL),
                width_mm = NATURE_WIDTH_M, height_mm = 70)
  }
}

# ── G. 物种趋势分类柱状 + 密度 + 热点地图 ───────────────────────────
sp_cls <- read_csv_safe(v4_file("results", paste0("table_species_trend_classify_", RUN_LABEL)))
sp_trd <- read_csv_safe(v4_file("results", paste0("table_species_trend_", RUN_LABEL)))
hotspot <- read_csv_safe(v4_file("results", paste0("table_species_hotspot_", RUN_LABEL)))

if (!is.null(sp_cls)) {
  cls_counts <- sp_cls |>
    count(trend_class, name = "n") |>
    mutate(trend_class = factor(trend_class,
                                 levels = c("contracting", "stable", "expanding"),
                                 labels = c("Contracting", "Stable", "Expanding")))
  p_cls <- ggplot(cls_counts, aes(trend_class, n, fill = trend_class)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n), vjust = -0.5, size = 7.5 / 2.835) +
    scale_fill_manual(values = c("Contracting" = "#B2182B",
                                  "Stable"     = "#969696",
                                  "Expanding"  = "#2166AC"),
                      name = NULL) +
    labs(x = NULL, y = "Number of species") +
    theme_nature_pub() + theme(legend.position = "none")
  save_nature(p_cls, paste0("fig_v4_species_trend_classification_", RUN_LABEL),
              width_mm = NATURE_WIDTH_S, height_mm = 70)
}

if (!is.null(sp_trd)) {
  sp_t <- sp_trd |> filter(method == "theil_sen")
  p_d <- ggplot(sp_t, aes(mean)) +
    geom_density(fill = NATURE_ACCENT, alpha = 0.3, colour = NATURE_ACCENT) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.2) +
    labs(x = "Species occupancy trend (Theil-Sen slope)", y = "Density") +
    theme_nature_pub()
  save_nature(p_d, paste0("fig_v4_species_trend_density_", RUN_LABEL),
              width_mm = NATURE_WIDTH_S, height_mm = 60)
}

if (!is.null(hotspot)) {
  hs_period_ok <- hotspot |> filter(!is.na(period))
  if (nrow(hs_period_ok) > 0) {
    max_p <- max(as.numeric(gsub("P", "", hs_period_ok$period)), na.rm = TRUE)
    hs_final <- hs_period_ok |> filter(period == paste0("P", max_p))
    hs_sf <- grid_sf |> left_join(hs_final, by = "grid_cell") |>
      filter(!is.na(net_trend_index))
    if (nrow(hs_sf) > 0) {
      p_hs <- ggplot(hs_sf) +
        geom_sf(aes(fill = net_trend_index),
                colour = scales::alpha("white", 0.08), linewidth = 0.04) +
        geom_sf(data = china_layers, fill = NA, colour = "grey50",
                linewidth = NATURE_LINE_MAP) +
        scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                             midpoint = 0, name = "Net trend index",
                             na.value = "grey92") +
        v4_china_coord() + theme_nature_map()
      save_nature(p_hs, paste0("fig_v4_species_hotspot_map_", RUN_LABEL),
                  width_mm = NATURE_WIDTH_M, height_mm = 80)
    }
  }
}

# ── H. 物种 × 协变量系数 heatmap（v4 新增） ──────────────────────────
beta_sp <- read_csv_safe(v4_file("results", paste0("table_beta_species_", RUN_LABEL)))
if (!is.null(beta_sp) && nrow(beta_sp) > 0) {
  # 把 param（如 "beta.1", "beta.2"...）映射到协变量名
  occ_terms <- c("(Intercept)", "bio4", "bio7", "bio11", "bio13",
                 "elev_mean", "elev_sd", "texture_shannon",
                 "habitat_diversity_shannon",
                 "hfi_mean", "landcover_built", "landcover_cropland",
                 "centroid_lon", "centroid_lat", "year_scaled")
  if (length(unique(beta_sp$param)) <= length(occ_terms)) {
    param_map <- setNames(head(occ_terms, length(unique(beta_sp$param))),
                          sort(unique(beta_sp$param)))
    beta_sp$covariate <- recode(beta_sp$param, !!!param_map)
  } else {
    beta_sp$covariate <- beta_sp$param
  }

  beta_sp <- beta_sp |>
    mutate(sig95 = (q025 > 0) | (q975 < 0))

  # 只保留前 60 物种（避免图过密）
  top_sp <- beta_sp |>
    group_by(species) |>
    summarise(s = sum(abs(mean), na.rm = TRUE), .groups = "drop") |>
    arrange(desc(s)) |> slice_head(n = 60) |> pull(species)
  beta_show <- beta_sp |> filter(species %in% top_sp,
                                  covariate != "(Intercept)")

  p_hm <- ggplot(beta_show, aes(covariate,
                                 fct_reorder(species, mean, .fun = sum),
                                 fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.05) +
    geom_text(data = beta_show |> filter(sig95),
              aes(label = "*"), size = 2.5, colour = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, name = "Posterior mean") +
    labs(x = NULL, y = NULL,
         title = "Species × covariate occupancy coefficients (* = 95% CRI excludes 0)") +
    theme_nature_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5))
  save_nature(p_hm, paste0("fig_v4_species_covariate_heatmap_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 160)
}

# ── I. brms 驱动回归系数 forest plot ─────────────────────────────────
brms_coefs <- read_csv_safe(v4_file("results",
                                     paste0("table_brms_driver_coefs_", RUN_LABEL)))
if (!is.null(brms_coefs)) {
  fp <- brms_coefs |>
    filter(grepl("^z_", term) | term %in% c("Intercept"))  |>
    mutate(sig95 = (q025 > 0) | (q975 < 0),
           term = gsub("^z_", "", term))
  p_fp <- ggplot(fp, aes(estimate, fct_reorder(term, estimate),
                          colour = factor(gp_k))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.2) +
    geom_errorbarh(aes(xmin = q025, xmax = q975),
                   height = 0.2, position = position_dodge(0.5),
                   linewidth = 0.3) +
    geom_point(position = position_dodge(0.5), size = 1.5) +
    facet_wrap(~ metric, scales = "free_x", ncol = 2) +
    scale_colour_brewer(palette = "Set2", name = "GP k") +
    labs(x = "Standardized coefficient (95% CRI)", y = NULL,
         title = "brms driver regression coefficients (horseshoe + GP)") +
    theme_nature_pub()
  save_nature(p_fp, paste0("fig_v4_brms_driver_forest_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 130)
}

# ── J. Baselga 全局比例柱状（turnover vs nestedness） ──────────────
bg_g <- read_csv_safe(v4_file("results", paste0("table_baselga_global_", RUN_LABEL)))
if (!is.null(bg_g)) {
  bg_long <- bg_g |>
    select(period_pair, prop_turnover_mean) |>
    mutate(prop_nestedness = 1 - prop_turnover_mean) |>
    pivot_longer(cols = c(prop_turnover_mean, prop_nestedness),
                  names_to = "component", values_to = "proportion") |>
    mutate(component = recode(component,
                              prop_turnover_mean = "Turnover",
                              prop_nestedness = "Nestedness"))
  p_bg <- ggplot(bg_long, aes(period_pair, proportion, fill = component)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("Turnover" = NATURE_ACCENT,
                                  "Nestedness" = "#D4DEE4"),
                      name = NULL) +
    labs(x = "Period pair", y = "Proportion of β diversity") +
    theme_nature_pub() + theme(legend.position = "top")
  save_nature(p_bg, paste0("fig_v4_baselga_proportion_", RUN_LABEL),
              width_mm = NATURE_WIDTH_M, height_mm = 70)
}

log_time("06", "v4 publication figures assembled")
