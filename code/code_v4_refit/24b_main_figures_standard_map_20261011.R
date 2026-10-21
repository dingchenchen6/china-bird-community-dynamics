#!/usr/bin/env Rscript
# ============================================================
# 24b_main_figures_standard_map_20261011.R
#
# Scientific question / 科学问题:
#   主文 Figure 1–6 展示中国鸟类群落 2000–2024 年的多维多样性格局、趋势与观测偏差。
#   本脚本只改呈现方式：地图改用中国标准地图，修复 Figure 1a 文字溢出，统一导出可编辑 PPTX。
#   Main Figures 1-6 of the community-dynamics manuscript; presentation-only revision.
#
# Objective / 分析目标:
#   以 24_make_nature_submission_figures.R 的同一数据处理代码重绘 Figure 1–6：
#   Fig 2/3 用标准地图（省面.shp 拆分、十段线仅在南海插图中、固定 bbox、缺数据格显灰并加图例键）；
#   Fig 1a 按实测宽度断行；Fig 6a 注释加白底，避免被 1:1 虚线穿过。数值、模型、筛选一律不变。
#   Redraw Figs 1-6 with the identical data code; only maps, text fit and export change.
#
# Input data / 输入数据:
#   results_v3/table_*.csv（与 24 号脚本相同的 14 张表）; data/derived_v2/china_grid_100km_v2.rds;
#   data/china_boundary/{省面,十段线}.shp
# Main workflow / 主要流程:
#   读表 → 逐图构建（与 24 号脚本相同的数据代码）→ 标准地图图层 → 导出 PDF/PNG/TIFF/PPTX + 源数据
#   read tables -> build figures (same data code) -> standard map -> export
# Expected output / 预期输出:
#   main/Figure_1..6_*.{pdf,png,tiff,pptx}; source_data/SourceData_Figure_1..6.csv;
#   qa/Figure_1a_text_fit.csv; qa/manifest_main.csv
# Key assumptions / 关键假设:
#   - 源数据 CSV 必须与 Shared/Figures/source_data 逐字节一致（QA 脚本核验）。
#   - 地图 PPTX 中格子填色层为嵌入栅格（Fig 2 有 30 个面板 × 约 1,700 个格子），
#     其余（省界、国界、十段线、插图框、文字、图例）均为可编辑矢量。
#     In map PPTX files the cell-fill layer is an embedded raster; everything else is vector.
# Main packages / 主要包: sf, dplyr, tidyr, readr, ggplot2, patchwork, scales, officer, rvg, terra
# Output directory / 输出路径: figures_community_submission_20261011/{main,source_data,qa}
# ============================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(stringr); library(forcats); library(grid)
})
source(file.path(getwd(), "code_v3", "utils_figs_export_20261011.R"))

RESULTS_DIR <- file.path(PROJECT_ROOT, "results_v3")
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUT_MAIN <- file.path(FIG_OUT_ROOT, "main")
SRC_DIR <- file.path(FIG_OUT_ROOT, "source_data")
QA_DIR <- file.path(FIG_OUT_ROOT, "qa")
for (d in c(OUT_MAIN, SRC_DIR, QA_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── 以下辅助函数与样式逐字沿用 24 号脚本 / helpers and style copied verbatim from script 24 ──
read_csv_req <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}
safe_make_valid <- function(x) {
  out <- try(sf::st_make_valid(x), silent = TRUE)
  if (inherits(out, "try-error")) suppressWarnings(sf::st_buffer(x, 0)) else out
}
safe_rescale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) return(rep(NA_real_, length(x)))
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  scales::rescale(x, to = c(0, 1), from = rng)
}
COL <- c(ink = "#202A33", neutral_dark = "#4D4D4D", neutral = "#8C96A3", neutral_light = "#D9DEE3",
         panel = "#F3F5F7", blue = "#0072B2", sky = "#56B4E9", orange = "#D55E00", gold = "#E69F00",
         purple = "#CC79A7", teal = "#009E73")
pal_dataset <- c("200-species spatial" = COL[["blue"]], "500-species temporal" = COL[["purple"]])
pal_class <- c("Increasing" = COL[["blue"]], "Stable" = COL[["neutral"]], "Declining" = COL[["orange"]])
pal_driver <- c("Spatial baseline" = COL[["neutral_dark"]], "Climate change" = COL[["blue"]],
                "Land-use change" = COL[["gold"]], "Human pressure change" = COL[["purple"]])
theme_nature <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(text = element_text(colour = COL[["ink"]]),
          axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 0.4, colour = COL[["ink"]]),
          legend.title = element_text(size = base_size - 0.2, face = "bold"),
          legend.text = element_text(size = base_size - 0.5),
          legend.key.height = unit(2.5, "mm"),
          strip.background = element_rect(fill = "#EEF1F4", colour = NA),
          strip.text = element_text(size = base_size - 0.1, face = "bold"),
          plot.title = element_text(size = base_size + 0.4, face = "bold", margin = margin(b = 2)),
          plot.subtitle = element_text(size = base_size - 0.2, colour = "#52616B", margin = margin(b = 3)),
          plot.caption = element_text(size = base_size - 1, colour = "#606D76", hjust = 0),
          plot.tag = element_text(size = 8, face = "bold", colour = "black"),
          panel.grid = element_blank(), plot.margin = margin(3, 3, 3, 3))
}
theme_map_nature <- function(base_size = 6.2) {
  theme_void(base_size = base_size, base_family = "Helvetica") +
    theme(strip.background = element_rect(fill = "#EEF1F4", colour = NA),
          strip.text = element_text(size = base_size, face = "bold", colour = COL[["ink"]]),
          strip.placement = "outside", legend.position = "bottom",
          legend.title = element_text(size = base_size - 0.1, face = "bold"),
          legend.text = element_text(size = base_size - 0.4),
          legend.key.width = unit(15, "mm"), panel.spacing = unit(1.0, "mm"),
          plot.tag = element_text(size = 8, face = "bold"), plot.margin = margin(2, 2, 2, 2))
}
period_label <- c(P1 = "2000-2004", P2 = "2005-2009", P3 = "2010-2014", P4 = "2015-2019", P5 = "2020-2024")
metric_label <- c(corrected_richness = "Richness", shannon = "Shannon diversity", inv_simpson = "Inverse Simpson",
                  pd_prob_mctavish = "Phylogenetic diversity", mpd_prob_mctavish = "Mean phylogenetic distance",
                  trait_volume = "Functional volume", rao_q = "Rao's Q", feve_fund = "Functional evenness")
metric_group <- c(corrected_richness = "Taxonomic", shannon = "Taxonomic", inv_simpson = "Taxonomic",
                  pd_prob_mctavish = "Phylogenetic", mpd_prob_mctavish = "Phylogenetic",
                  trait_volume = "Functional", rao_q = "Functional", feve_fund = "Functional")

MANIFEST <- list()
save_main <- function(plot, stem, width_mm, height_mm, ppt_plot = plot, mode = "native_vector") {
  # Nature 尺寸约束沿用 24 号脚本 / Nature size contract as in script 24
  if (width_mm > 183 || height_mm > 170) stop("Figure exceeds Nature size contract: ", stem)
  MANIFEST[[length(MANIFEST) + 1L]] <<- save_figure_set(plot, OUT_MAIN, stem, width_mm, height_mm,
                                                        ppt_plot = ppt_plot, pptx_mode = mode)
}

# ── 输入 / Inputs (identical to script 24) ───────────────────────────
grid_path <- c(file.path(DATA_DIR, "derived_v3", "china_grid_100km_v3.rds"),
               file.path(DATA_DIR, "derived_v2", "china_grid_100km_v2.rds"))
grid_path <- grid_path[file.exists(grid_path)][1]
grid_sf <- readRDS(grid_path) |> safe_make_valid() |> st_transform(4326) |> select(grid_cell, geometry)
metrics_200 <- read_csv_req(file.path(RESULTS_DIR, "table_community_metrics_with_cri_v3_full_200sp_ar1_spatial_extended.csv"))
metrics_500 <- read_csv_req(file.path(RESULTS_DIR, "table_community_metrics_with_cri_v3_full_500sp_ar1_temporal_extended.csv"))
trend_200 <- read_csv_req(file.path(RESULTS_DIR, "table_trend_summary_v3_full_200sp_ar1_spatial_extended.csv"))
beta_200 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_200sp_ar1_spatial_extended.csv"))
beta_500 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_500sp_ar1_temporal_extended.csv"))
class_500 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_500sp_ar1_temporal_extended.csv"))
trait_corr <- read_csv_req(file.path(RESULTS_DIR, "table_trait_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"))
naive_500 <- read_csv_req(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_500sp_ar1_temporal_extended.csv"))
rf_group <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_group_summary_trend.csv"))
varpart_500 <- read_csv_req(file.path(RESULTS_DIR, "table_varpart_corrected_richness_v3_full_500sp_ar1_temporal_extended.csv"))
conv200 <- read_csv_req(file.path(RESULTS_DIR, "table_convergence_diagnostics_v3_full_200sp_ar1_spatial_4chain.csv"))
conv500 <- read_csv_req(file.path(RESULTS_DIR, "table_convergence_diagnostics_v3_full_500sp_ar1_temporal_4chain.csv"))

# 标准地图图层与绘图网格 / standard-map layers and plotting grid
CHN <- china_layers_pub()
grid_plot <- grid_for_maps(grid_sf, model_cells = unique(metrics_200$grid_cell))
SEAM_LW <- 0.08   # 格子同色描边宽度（mm），消除抗锯齿白缝 / same-colour stroke closing seams
INSET_FACET <- scs_inset_pub(CHN, line_col = COL[["neutral"]], line_lw = 0.1, dash_col = "grey25",
                             dash_lw = 0.25, frame_col = "grey45", frame_lw = 0.2)

source(file.path(PROJECT_ROOT, "code_v3", "24c_figure1_layout_20261011.R"))   # Figure 1 (layout fix)

# ── Figure 2：五个时段的空间多维多样性 / spatial multidiversity across five periods ──
map_metrics <- c("corrected_richness", "shannon", "trait_volume", "pd_prob_mctavish", "mpd_prob_mctavish", "rao_q")
fig2_data <- expand_grid(grid_cell = grid_sf$grid_cell, metric = map_metrics, period = names(period_label)) |>
  left_join(metrics_200 |> select(grid_cell, metric, period, value_mean, value_l95, value_u95),
            by = c("grid_cell", "metric", "period")) |>
  group_by(metric) |>
  mutate(relative_value = safe_rescale(value_mean)) |>
  ungroup() |>
  mutate(metric_plot = factor(unname(metric_label[metric]), levels = unname(metric_label[map_metrics])),
         period_plot = factor(unname(period_label[period]), levels = unname(period_label)))
write_csv(fig2_data |> select(grid_cell, period, metric, value_mean, value_l95, value_u95, relative_value),
          file.path(SRC_DIR, "SourceData_Figure_2.csv"))

map_frame <- function(p, base_size) {
  # 标准地图公共图层：省界、国界、南海插图（含十段线）、无数据图例键、固定 bbox
  # Shared standard-map layers: provinces, outline, SCS inset with dash line, no-data key, fixed bbox
  p + geom_sf(data = CHN$prov, fill = NA, colour = "#D8DEE3", linewidth = 0.055) +
    geom_sf(data = CHN$outline, fill = NA, colour = COL[["ink"]], linewidth = 0.24) +
    scs_inset_layer_pub(INSET_FACET) +
    nodata_key() +
    coord_sf(xlim = MAP_BB[1:2], ylim = MAP_BB[3:4], expand = FALSE, datum = NA) +
    theme_map_nature(base_size)
}
fig2_sf <- grid_plot |> inner_join(fig2_data, by = "grid_cell")
fig2_px <- cell_pixels(grid_plot, res_deg = 0.125) |>   # 1/8 度：二进制精确，像元不错位 / exact binary fraction
  inner_join(fig2_data |> select(grid_cell, metric_plot, period_plot, relative_value),
             by = "grid_cell", relationship = "many-to-many")
fig2_build <- function(fill_layer) {
  map_frame(ggplot() + fill_layer, 6.0) +
    facet_grid(rows = vars(metric_plot), cols = vars(period_plot), switch = "y") +
    scale_fill_gradientn(colours = c("#F7FBFF", "#D8EAF4", "#8FC3DD", "#3B8BC2", "#075A8C"),
                         na.value = NA_GREY, breaks = c(0, 0.5, 1), labels = c("Low", "Mid", "High"),
                         name = "Relative value within metric") +
    guides(fill = guide_colourbar(order = 1, title.position = "top", barwidth = unit(42, "mm"),
                                  barheight = unit(2.3, "mm"))) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1), legend.position = "bottom")
}
fig2 <- fig2_build(geom_sf(data = fig2_sf, aes(fill = relative_value, colour = after_scale(fill)), linewidth = SEAM_LW))
fig2_ppt <- fig2_build(geom_raster(data = fig2_px, aes(x, y, fill = relative_value)))
FIG2_H <- as.numeric(Sys.getenv("FIG2_H", "150"))
save_main(fig2, "Figure_2_spatial_multidiversity", 183, FIG2_H, ppt_plot = fig2_ppt,
          mode = "vector_with_raster_cell_fill")

# ── Figure 3：方向性变化的空间指纹 / spatial fingerprints of directional change ──
fig3_data <- trend_200 |>
  filter(metric %in% map_metrics, method %in% c("ols", "theil_sen", "theilsen") | is.na(method)) |>
  group_by(grid_cell, metric) |>
  slice(1) |>
  ungroup() |>
  mutate(credible = q025 > 0 | q975 < 0) |>
  group_by(metric) |>
  mutate(standardized_trend = mean / max(abs(mean), na.rm = TRUE)) |>
  ungroup() |>
  mutate(metric_plot = factor(unname(metric_label[metric]), levels = unname(metric_label[map_metrics])))
fig3_full <- expand_grid(grid_cell = grid_sf$grid_cell, metric = map_metrics) |>
  left_join(fig3_data |> select(grid_cell, metric, mean, q025, q975, credible, standardized_trend),
            by = c("grid_cell", "metric")) |>
  mutate(metric_plot = factor(unname(metric_label[metric]), levels = unname(metric_label[map_metrics])))
write_csv(fig3_full |> select(grid_cell, metric, mean, q025, q975, credible, standardized_trend),
          file.path(SRC_DIR, "SourceData_Figure_3.csv"))

fig3_sf <- grid_plot |> inner_join(fig3_full, by = "grid_cell")
fig3_px <- cell_pixels(grid_plot, res_deg = 0.0625) |>  # 1/16 度 / 1/16 degree
  inner_join(fig3_full |> select(grid_cell, metric_plot, standardized_trend, credible),
             by = "grid_cell", relationship = "many-to-many")
# 透明度（区间跨 0 的格子变浅）换算为白底等效色，保证同色描边不产生缝或深色边
# Alpha (paler cells = interval overlaps zero) converted to the equivalent opaque colour.
fig3_aes <- aes(fill = stage(standardized_trend, after_scale = blend_white(fill, alpha)),
                colour = after_scale(blend_white(fill, alpha)),
                alpha = stage(credible, after_scale = rep(1, length(alpha))))
fig3_build <- function(fill_layer) {
  map_frame(ggplot() + fill_layer, 6.2) +
    facet_wrap(~metric_plot, ncol = 3) +
    scale_fill_gradientn(colours = c(COL[["blue"]], "#B8D8EB", "#F7F7F7", "#F4C27A", COL[["orange"]]),
                         limits = c(-1, 1), oob = squish, na.value = NA_GREY,
                         breaks = c(-1, 0, 1), labels = c("Decrease", "0", "Increase"), name = "Standardized trend") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.38), guide = "none", na.value = 1) +
    guides(fill = guide_colourbar(order = 1, title.position = "top", barwidth = unit(42, "mm"),
                                  barheight = unit(2.3, "mm")))
}
fig3 <- fig3_build(geom_sf(data = fig3_sf, fig3_aes, linewidth = SEAM_LW))
fig3_ppt <- fig3_build(geom_raster(data = fig3_px, aes(x, y, fill = stage(standardized_trend,
  after_scale = blend_white(fill, alpha)), alpha = stage(credible, after_scale = rep(1, length(alpha))))))
FIG3_H <- as.numeric(Sys.getenv("FIG3_H", "104"))
save_main(fig3, "Figure_3_spatial_trend_fingerprints", 183, FIG3_H, ppt_plot = fig3_ppt,
          mode = "vector_with_raster_cell_fill")

source(file.path(PROJECT_ROOT, "code_v3", "24d_figures4to6_20261011.R"))   # Figures 4-6 (unchanged code)

manifest <- bind_rows(MANIFEST)
write_csv(manifest, file.path(QA_DIR, "manifest_main.csv"))
message("[24b] main figures written: ", nrow(manifest), " -> ", OUT_MAIN)
message("[24b] dropped artefact grid cells (baked-in inset zone): ",
        length(attr(grid_plot, "dropped_artefact_cells")))
