#!/usr/bin/env Rscript
# ============================================================
# 34b_pa_figures_english_standard_map_20261011.R
#
# Scientific question / 科学问题:
#   保护地内外的鸟类群落轨迹、空间同质化速率、特化/泛化种占域、事件研究、保护年限剂量-反应、
#   匹配平衡与保护规划方案（脚本 31–33 的结果）。
#   Protected-area (PA) contrasts, matching balance and conservation prioritization (scripts 31-33).
#
# Objective / 分析目标:
#   按 34_pa_publication_figures.R（仓库 china-bird-community-dynamics 提交 3f52951 版本，即生成
#   现有 SI 图的版本）的同一数据处理重绘 PA1–PA3、PA5、PA5b、PA6–PA9，只改呈现：
#   (1) 所有中文标题、坐标轴、图例、注释改为简洁英文；(2) PA3 只画真实的匹配前/后 SMD
#   （原脚本 matches("smd|...|after") 把逻辑列 balanced_after 当作 SMD，在 |SMD| = 1.00 处画出伪点）；
#   (3) PA1、PA8 改为中国标准地图（真实网格多边形，十段线仅在南海插图，缺数据格显灰）；
#   (4) PA5b 图例如实标注（原图把保护地内/外两条斜率误标为"区间含 0"）；
#   (5) PA7 面板顺序与 PA2/PA9 一致。数值、筛选、统计量一律不变。
#   Same data handling as the script-34 version that produced the SI figures; presentation only.
#
# Input data / 输入数据:
#   submission_packages_20261010/Shared/Supplementary/S12_protected_areas/table_pa_*_v3_full_200sp_ar1_spatial.csv,
#   S13_prioritization/table_priority_solutions_v3_full_200sp_ar1_spatial.csv（与 china-bird-community-dynamics/
#   results/results_v3 逐字节一致，已核验）; data/derived_v2/china_grid_100km_v2.rds; data/china_boundary/*.shp
# Main workflow / 主要流程: 读表 → 英文标签 → 作图（同一数据代码）→ 标准地图 → 导出四种格式
#   read tables -> English labels -> plots (same data code) -> standard map -> export
# Expected output / 预期输出: protected_area/fig_pa{1,2,3,5,5b,6,7,8,9}_*_v3_full_200sp_ar1_spatial.{pdf,png,tiff,pptx};
#   qa/fig_pa3_plotted_smd.csv; qa/manifest_protected_area.csv
# Key assumptions / 关键假设:
#   - 保护分层沿用脚本 33：outside ≤ 5%，partial 5–30%，inside > 30% 网格面积为保护区（图例沿用原图写法）。
#   - PA1/PA8 中模型之外的网格（无 PA 结果）显示为灰色"No data"，不再被省略。
#   - PA8 的 PPTX 中格子填色层为嵌入栅格（6 个面板 × 约 1,700 格），其余为可编辑矢量。
# Main packages / 主要包: ggplot2, sf, dplyr, readr, tidyr, officer, rvg
# Output directory / 输出路径: figures_community_submission_20261011/protected_area
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(dplyr); library(readr); library(tidyr)
})
source(file.path(getwd(), "code_v3", "utils_figs_export_20261011.R"))

RUN_LABEL <- "v3_full_200sp_ar1_spatial"
SI_DIR <- file.path(PROJECT_ROOT, "submission_packages_20261010", "Shared", "Supplementary")
OUT_PA <- file.path(FIG_OUT_ROOT, "protected_area")
QA_DIR <- file.path(FIG_OUT_ROOT, "qa")
for (d in c(OUT_PA, QA_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

read_si <- function(sub, stem) {
  # 只读 SI 中已核验的结果表 / read the verified SI copies of the result tables
  p <- file.path(SI_DIR, sub, paste0(stem, "_", RUN_LABEL, ".csv"))
  if (!file.exists(p)) stop("Missing input: ", p, call. = FALSE)
  read_csv(p, show_col_types = FALSE)
}

# ── 主题：逐字复制自 utils_mapping.R（不 source 该文件，以免触发配置与建目录副作用）──
# Themes copied verbatim from utils_mapping.R (not sourced, to avoid config side effects).
theme_nature_pub <- function(base_size = 7.5, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.border = element_rect(colour = "black", linewidth = 0.25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black", linewidth = 0.25),
          axis.ticks = element_line(colour = "black", linewidth = 0.15),
          axis.text = element_text(colour = "black", size = base_size),
          axis.title = element_text(colour = "black", size = base_size, face = "plain"),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.key.size = unit(3, "pt"), legend.text = element_text(size = base_size - 1),
          legend.title = element_text(size = base_size, face = "plain"), legend.position = "right",
          strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.15),
          strip.text = element_text(size = base_size, face = "plain", margin = margin(2, 2, 2, 2)),
          plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = base_size, hjust = 0, colour = "grey30"),
          plot.margin = margin(4, 4, 4, 4, unit = "pt"))
}
theme_nature_map <- function(base_size = 7.5, base_family = "Arial") {
  theme_nature_pub(base_size, base_family) %+replace%
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_text(size = base_size - 1.5, colour = "grey40"), axis.title = element_blank(),
          panel.border = element_rect(colour = "grey60", linewidth = 0.15),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm"))
}
# 副标题从整幅图左缘起排，并按实测宽度断行，避免原图的右侧截断
# Subtitles start at the plot edge and are wrapped on measured widths (the originals were cut off).
theme_pa_extra <- theme(plot.title.position = "plot")
wrap_sub <- function(text, fig_width_mm, size_pt = 7.5) {
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 8, height = 6)
  on.exit(grDevices::dev.off())
  paste(wrap_to_width(text, fig_width_mm - 6, size_pt, family = "Arial"), collapse = "\n")
}

# ── 配色与英文标签 / palette (unchanged) and English labels ─────────────
COL_IN <- "#2E75B6"; COL_OUT <- "#D97A34"; COL_MID <- "#9AA7B4"
PA_COLS <- c(protected = COL_IN, partial = COL_MID, unprotected = COL_OUT)
PA_LABS <- c(protected = "Inside PAs (≥30% cover)", partial = "Partially protected",
             unprotected = "Outside PAs (<5% cover)")
PERIOD_LABS <- c("2000-04", "2005-09", "2010-14", "2015-19", "2020-24")
metric_label <- c(corrected_richness = "Corrected richness", shannon = "Shannon diversity",
                  inv_simpson = "Inverse Simpson", trait_volume = "Functional volume", rao_q = "Rao's Q",
                  fric_prob = "Functional richness", fdis_prob = "Functional dispersion",
                  feve_fund = "Functional evenness", pd_prob = "Phylogenetic diversity")
lab_of <- function(x) unname(ifelse(x %in% names(metric_label), metric_label[x], x))
Y_TRAJ <- "Posterior mean (shading, ±1 s.e.)"

MANIFEST <- list()
save_pa <- function(plot, stem, w, h, ppt_plot = plot, mode = "native_vector") {
  MANIFEST[[length(MANIFEST) + 1L]] <<- save_figure_set(plot, OUT_PA, paste0(stem, "_", RUN_LABEL), w, h,
                                                        ppt_plot = ppt_plot, pptx_mode = mode)
}

# ── PA2 多维多样性逐期轨迹：保护地内 vs 外 / trajectories inside vs outside ──
traj <- read_si("S12_protected_areas", "table_pa_trajectory")
focal <- intersect(unique(traj$metric), c("corrected_richness", "shannon", "trait_volume", "rao_q"))
d2 <- traj |> filter(metric %in% focal, !is.na(pa_stratum)) |>
  mutate(metric_lab = factor(lab_of(metric), levels = lab_of(focal)),
         period_num = as.integer(sub("P", "", period)),
         pa_stratum = factor(pa_stratum, levels = names(PA_COLS)))
p2 <- ggplot(d2, aes(period_num, mean, colour = pa_stratum, fill = pa_stratum)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
  facet_wrap(~ metric_lab, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = PA_COLS, labels = PA_LABS, name = NULL) +
  scale_fill_manual(values = PA_COLS, labels = PA_LABS, name = NULL) +
  scale_x_continuous(breaks = 1:5, labels = PERIOD_LABS) +
  labs(x = NULL, y = Y_TRAJ) +
  theme_nature_pub() + theme_pa_extra +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
save_pa(p2, "fig_pa2_trajectories", 183, 68)

# ── PA5 空间 beta 内外对比 / spatial beta inside vs outside ───────────────
bet <- read_si("S12_protected_areas", "table_pa_beta_contrast")
has_ci <- all(c("spatial_beta_inside_lwr", "spatial_beta_outside_lwr") %in% names(bet))
d5 <- if (has_ci) {
  bet |> select(period, inside = spatial_beta_inside, outside = spatial_beta_outside,
                in_l = spatial_beta_inside_lwr, in_u = spatial_beta_inside_upr,
                out_l = spatial_beta_outside_lwr, out_u = spatial_beta_outside_upr) |>
    pivot_longer(c(inside, outside), names_to = "grp", values_to = "beta") |>
    mutate(lwr = ifelse(grp == "inside", in_l, out_l), upr = ifelse(grp == "inside", in_u, out_u),
           period_num = as.integer(sub("P", "", period)))
} else {
  bet |> select(period, inside = spatial_beta_inside, outside = spatial_beta_outside,
                sd_in = spatial_beta_inside_sd, sd_out = spatial_beta_outside_sd) |>
    pivot_longer(c(inside, outside), names_to = "grp", values_to = "beta") |>
    mutate(sd = ifelse(grp == "inside", sd_in, sd_out), lwr = beta - sd, upr = beta + sd,
           period_num = as.integer(sub("P", "", period)))
}
sub5 <- if (has_ci) "Shading: bootstrap 95% intervals; a decline indicates community convergence (homogenization)" else
  "A decline indicates community convergence (homogenization); the steeper group homogenizes faster"
p5 <- ggplot(d5, aes(period_num, beta, colour = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(inside = COL_IN, outside = COL_OUT), labels = c("Inside PAs", "Outside PAs"), name = NULL) +
  scale_fill_manual(values = c(inside = COL_IN, outside = COL_OUT), labels = c("Inside PAs", "Outside PAs"), name = NULL) +
  scale_x_continuous(breaks = 1:5, labels = PERIOD_LABS) +
  labs(x = NULL, y = "Spatial beta diversity among cells", subtitle = wrap_sub(sub5, 89)) +
  theme_nature_pub() + theme_pa_extra +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
save_pa(p5, "fig_pa5_spatial_beta", 89, 72)

# ── PA5b 同质化速率的内外差异 / homogenization-rate contrast ───────────────
# 图例改为如实标注：只有"差值"一行检验区间是否含 0；内/外两条斜率是组内估计，不是"区间含 0"。
# Legend now states what the colour encodes; inside/outside rows were mislabelled "interval includes 0".
rate <- read_si("S12_protected_areas", "table_pa_beta_rate")
SIG_COLS <- c("Difference: 95% interval excludes 0" = COL_IN, "Difference: 95% interval includes 0" = "grey60",
              "Inside or outside slope" = "grey60")
d5b <- rate |>
  mutate(grp = factor(group, levels = c("inside", "outside", "difference"),
                      labels = c("Inside PAs", "Outside PAs", "Difference (inside − outside)")),
         sig = case_when(group == "difference" & (lwr > 0 | upr < 0) ~ names(SIG_COLS)[1],
                         group == "difference" ~ names(SIG_COLS)[2], TRUE ~ names(SIG_COLS)[3]))
p5b <- ggplot(d5b, aes(slope, grp, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lwr, xmax = upr), orientation = "y", width = 0.16, linewidth = 0.5) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = SIG_COLS, name = NULL, guide = guide_legend(nrow = 2)) +
  # 窄图（89 mm）中 x 轴标题分两行，避免右端被截断 / two-line title so it is not clipped in an 89-mm figure
  labs(x = "Slope of spatial beta across periods\n(negative = convergence)", y = NULL,
       subtitle = wrap_sub("Error bars: bootstrap 95% intervals; inside and outside rates differ only if the difference interval excludes 0", 89)) +
  theme_nature_pub() + theme_pa_extra + theme(legend.position = "top")
save_pa(p5b, "fig_pa5b_homogenization_rate", 89, 64)

# ── PA6 特化种 vs 泛化种的占域差 / specialists vs generalists ─────────────
sp <- read_si("S12_protected_areas", "table_pa_species_occupancy")
d6 <- sp |> filter(!is.na(sp_group)) |> group_by(sp_group, period) |>
  summarise(diff = mean(psi_diff, na.rm = TRUE), se = sd(psi_diff, na.rm = TRUE) / sqrt(n()), .groups = "drop") |>
  mutate(period_num = as.integer(sub("P", "", period)), sp_group = factor(sp_group, levels = c("specialist", "generalist")))
SP_LABS <- c("Specialists (narrow habitat breadth)", "Generalists (broad habitat breadth)")
p6 <- ggplot(d6, aes(period_num, diff, colour = sp_group, fill = sp_group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = diff - se, ymax = diff + se), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(specialist = "#7B3294", generalist = COL_OUT), labels = SP_LABS, name = NULL,
                      guide = guide_legend(nrow = 2)) +
  scale_fill_manual(values = c(specialist = "#7B3294", generalist = COL_OUT), labels = SP_LABS, name = NULL,
                    guide = guide_legend(nrow = 2)) +
  scale_x_continuous(breaks = 1:5, labels = PERIOD_LABS) +
  labs(x = NULL, y = "Occupancy difference (inside − outside PAs)",
       subtitle = wrap_sub("Positive values: higher occupancy inside PAs; a widening specialist gap would indicate selective retention of specialists", 89)) +
  theme_nature_pub() + theme_pa_extra +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
save_pa(p6, "fig_pa6_specialist_generalist", 89, 76)

# ── PA7 事件研究 / event study ─────────────────────────────────────────
es <- read_si("S12_protected_areas", "table_pa_event_study")
PA7_ORDER <- c("corrected_richness", "rao_q", "trait_volume")       # 与 PA2/PA9 一致 / as in PA2 and PA9
d7 <- es |> filter(grepl("^ett-?[0-9]|^et", term), term != "(Intercept)") |>
  mutate(event_time = suppressWarnings(as.integer(gsub("[^0-9-]", "", term)))) |>
  filter(!is.na(event_time)) |>
  mutate(lwr = estimate - 1.96 * se, upr = estimate + 1.96 * se,
         metric_lab = factor(lab_of(metric), levels = lab_of(intersect(PA7_ORDER, unique(metric)))))
p7 <- ggplot(d7, aes(event_time, estimate)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = -0.5, linetype = 3, colour = COL_OUT, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.12, linewidth = 0.4, colour = COL_IN) +
  geom_point(size = 1.8, colour = COL_IN) +
  facet_wrap(~ metric_lab, scales = "free_y") +
  labs(x = "Period relative to reserve establishment (0 = establishment period)", y = "Standardized effect",
       subtitle = wrap_sub("Error bars: 95% CI (±1.96 s.e.). Pre-establishment coefficients should be near 0 (parallel trends); post-establishment departures from 0 indicate a protection effect", 183)) +
  theme_nature_pub() + theme_pa_extra
save_pa(p7, "fig_pa7_event_study", 183, 70)

# ── PA9 保护年限剂量-反应 / protection-age dose-response ──────────────────
age <- read_si("S12_protected_areas", "table_pa_age_trajectory")
focal9 <- intersect(unique(age$metric), c("corrected_richness", "trait_volume", "rao_q"))
AGE_COLS <- c("unprotected" = COL_OUT, "young(<15a)" = "#9AC3E0", "mid(15-30a)" = COL_IN, "old(>30a)" = "#12395E")
AGE_LABS <- c("unprotected" = "Unprotected (<5% cover)", "young(<15a)" = "<15 years",
              "mid(15-30a)" = "15–30 years", "old(>30a)" = ">30 years")
d9 <- age |> filter(metric %in% focal9, age_class %in% names(AGE_COLS)) |>
  mutate(metric_lab = factor(lab_of(metric), levels = lab_of(focal9)),
         period_num = as.integer(sub("P", "", period)), age_class = factor(age_class, levels = names(AGE_COLS)))
p9 <- ggplot(d9, aes(period_num, mean, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
  facet_wrap(~ metric_lab, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = AGE_COLS, labels = function(b) unname(AGE_LABS[b]), name = "Protection age") +
  scale_fill_manual(values = AGE_COLS, labels = function(b) unname(AGE_LABS[b]), name = "Protection age") +
  scale_x_continuous(breaks = 1:5, labels = PERIOD_LABS) +
  labs(x = NULL, y = Y_TRAJ,
       subtitle = wrap_sub("Dose–response: if protection is effective, cells protected for longer should show better functional trajectories", 183)) +
  theme_nature_pub() + theme_pa_extra +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
save_pa(p9, "fig_pa9_protection_age", 183, 74)

source(file.path(PROJECT_ROOT, "code_v3", "34c_pa_balance_and_maps_20261011.R"))   # PA3, PA1, PA8

manifest <- bind_rows(MANIFEST)
write_csv(manifest, file.path(QA_DIR, "manifest_protected_area.csv"))
message("[34b] protected-area figures written: ", nrow(manifest), " -> ", OUT_PA)
