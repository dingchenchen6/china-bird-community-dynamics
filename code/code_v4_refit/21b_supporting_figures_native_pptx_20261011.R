#!/usr/bin/env Rscript
# ============================================================
# 21b_supporting_figures_native_pptx_20261011.R
#
# Scientific question / 科学问题:
#   支撑材料图 fig05–fig09：物种趋势类别、朴素与校正趋势对比、全国 Baselga 分解、
#   随机森林驱动因子重要性、性状与环境关联。
#   Supporting figures fig05-fig09 (species trend classes, naive vs corrected trends,
#   national Baselga partition, random-forest importance, trait/environment associations).
#
# Objective / 分析目标:
#   2026-07-14 的 PPTX（figures_top_journal_20260714_v2/）以 12.75 × 6.75 英寸重排后放在
#   10 × 7.5 英寸的 4:3 幻灯片上，约 3 英寸画面落在幻灯片外（fig07 的 P4–P5 区间整段不可见），
#   并多出文件名标题；因此不复用，改为用 21 号脚本的同一作图代码、按原始物理尺寸重新导出
#   原生可编辑 PPTX（幻灯片尺寸 = 图件尺寸），并同时重出 PDF/PNG/TIFF 以核对与投稿版逐像素一致。
#   Re-export the five figures from the identical plot code at native size (slide = figure).
#
# Input data / 输入数据: results_v3/ 下 21 号脚本使用的 10 张表（物种分类、物种趋势、朴素对比、
#   Baselga 全国、随机森林、性状与环境相关）
# Main workflow / 主要流程: 读表 → 逐字沿用 21 号脚本的作图函数 → 导出四种格式
#   read tables -> plot functions copied verbatim from script 21 -> export four formats
# Expected output / 预期输出: supporting/fig05..fig09_*_top_journal_editable.{pdf,png,tiff,pptx};
#   qa/manifest_supporting.csv
# Key assumptions / 关键假设: 作图代码、尺寸（英寸）与 21 号脚本完全相同；PNG 应与 07-14 版逐像素一致。
# Main packages / 主要包: dplyr, tidyr, readr, ggplot2, patchwork, scales, forcats, stringr, officer, rvg
# Output directory / 输出路径: figures_community_submission_20261011/supporting
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
  library(scales); library(forcats); library(stringr); library(grid)
})
source(file.path(getwd(), "code_v3", "utils_figs_export_20261011.R"))
RESULTS_DIR <- file.path(PROJECT_ROOT, "results_v3")
OUT_SUP <- file.path(FIG_OUT_ROOT, "supporting")
QA_DIR <- file.path(FIG_OUT_ROOT, "qa")
for (d in c(OUT_SUP, QA_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── 以下逐字沿用 21_make_top_journal_figures.R / verbatim from script 21 ─────────────
read_csv_req <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}
class_label <- c(expanding = "Increasing", stable = "Stable", contracting = "Declining",
                 increase = "Increasing", no_change = "Stable", decrease = "Declining")
group_label <- c(baseline = "Spatial baseline", climate_change = "Climate change",
                 landuse_change = "Land-use change", human_change = "Human pressure change")
term_label <- c(delta_t_mean = "Mean temperature change", delta_t_extreme = "Thermal-extreme change",
                delta_precip = "Precipitation change", delta_forest = "Forest-cover change",
                delta_cropland = "Cropland change", delta_crop = "Cropland change", delta_grass = "Grassland change",
                delta_urban = "Urban-cover change", delta_hfi = "Human footprint change",
                delta_hfp = "Human pressure change", centroid_lon = "Longitude", centroid_lat = "Latitude",
                elev_mean = "Elevation", bio4 = "Temperature seasonality", bio7 = "Temperature annual range",
                bio11 = "Mean temperature of coldest quarter", bio13 = "Precipitation of wettest month")
trait_label <- c(habitat_breadth = "Habitat breadth", diet_specialization = "Diet specialization",
                 migration_score = "Migratory tendency")
pal_classes <- c(Increasing = "#0072B2", Stable = "#8C96A3", Declining = "#D55E00")
pal_500 <- c("200 spatial species" = "#0072B2", "500 temporal species" = "#CC79A7")
pal_driver <- c("Spatial baseline" = "#4D4D4D", "Climate change" = "#0072B2", "Land-use change" = "#E69F00",
                "Human pressure change" = "#CC79A7", Other = "#7F8790")
theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(text = element_text(colour = "#1E2A32"),
          plot.title = element_text(face = "bold", size = base_size + 2.8, margin = margin(b = 4)),
          plot.subtitle = element_text(size = base_size, colour = "#52616B", margin = margin(b = 8)),
          plot.caption = element_text(size = base_size - 1.5, colour = "#6B7780", hjust = 0),
          axis.title = element_text(face = "bold", size = base_size),
          axis.text = element_text(size = base_size - 0.8, colour = "#35424A"),
          axis.line = element_line(colour = "#26343D", linewidth = 0.35),
          axis.ticks = element_line(colour = "#26343D", linewidth = 0.35),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = base_size, colour = "#23313A"),
          strip.background = element_rect(fill = "#EEF3F7", colour = NA),
          legend.title = element_text(face = "bold", size = base_size - 0.4),
          legend.text = element_text(size = base_size - 0.8), legend.position = "bottom",
          plot.margin = margin(8, 10, 8, 10))
}
wrap_label <- function(x, width = 24) stringr::str_wrap(x, width = width)
nice_term <- function(x) {
  out <- unname(term_label[x]); out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " ")); out
}
nice_group <- function(x) {
  out <- unname(group_label[x]); out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " ")); out
}
nice_trait <- function(x) {
  out <- unname(trait_label[x]); out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " ")); out
}

class_200 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_200sp_ar1_spatial_extended.csv"))
class_500 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_500sp_ar1_temporal_extended.csv"))
species_500 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_v3_full_500sp_ar1_temporal_extended.csv"))
naive_200 <- read_csv_req(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_200sp_ar1_spatial_extended.csv"))
naive_500 <- read_csv_req(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_500sp_ar1_temporal_extended.csv"))
beta_global_200 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_200sp_ar1_spatial_extended.csv"))
beta_global_500 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_500sp_ar1_temporal_extended.csv"))
rf_var <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_variable_summary_trend.csv"))
rf_group <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_group_summary_trend.csv"))
trait_500 <- read_csv_req(file.path(RESULTS_DIR, "table_trait_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"))
env_500 <- read_csv_req(file.path(RESULTS_DIR, "table_env_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"))

source(file.path(PROJECT_ROOT, "code_v3", "21c_supporting_plot_functions_20261011.R"))   # verbatim plot functions

# ── 导出：尺寸（英寸）与 21 号脚本相同 / export at the script-21 sizes (inches) ──────
SUP <- list(
  list(p = make_species_class_plot(class_200, class_500, species_500), stem = "fig05_species_trend_classes_500sp_top_journal_editable", w = 12.2, h = 6.6),
  list(p = make_naive_corrected_plot(naive_200, naive_500), stem = "fig06_naive_vs_corrected_trends_top_journal_editable", w = 10.6, h = 5.7),
  list(p = make_beta_global_plot(beta_global_200, beta_global_500), stem = "fig07_global_baselga_turnover_nestedness_top_journal_editable", w = 9.8, h = 6.6),
  list(p = make_rf_plot(rf_var, rf_group), stem = "fig08_random_forest_driver_importance_top_journal_editable", w = 12.1, h = 6.8),
  list(p = make_trait_env_plot(trait_500, env_500), stem = "fig09_trait_environment_mechanisms_500sp_top_journal_editable", w = 11.2, h = 5.8)
)
manifest <- bind_rows(lapply(SUP, function(f) save_figure_set(f$p, OUT_SUP, f$stem, f$w * 25.4, f$h * 25.4)))
write_csv(manifest, file.path(QA_DIR, "manifest_supporting.csv"))
message("[21b] supporting figures written: ", nrow(manifest), " -> ", OUT_SUP)
