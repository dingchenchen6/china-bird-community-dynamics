#!/usr/bin/env Rscript
# ============================================================
# 34_pa_publication_figures.R
#
# Scientific question / 科学问题:
#   把保护地内外的时空对比、成效评估与保护规划结果，做成可直接投稿的图件。
#   Publication-quality figures for the protected-area contrast, effectiveness
#   and prioritization analyses (scripts 31-33).
#
# Objective / 分析目标:
#   PA1 保护梯度地图（保护区覆盖率 + 保护区边界）
#   PA2 多维多样性逐期轨迹：保护地内 vs 外（核心图）
#   PA3 匹配协变量平衡图（love plot）—— 证明准实验设计有效，顶刊必备
#   PA4 处理效应森林图（各维度 ATT）
#   PA5 空间 beta（同质化速率）内外对比
#   PA6 特化种 vs 泛化种的占域分化
#   PA7 事件研究图（动态 DiD，检验平行趋势）
#   PA8 保护规划三方案对比地图 + 方案重叠度
#
# Input / 输入: 脚本 31/32/33 产出的结果表 + grid_environment + 保护区边界
# Output / 输出: figures_v3/fig_pa_*.{pdf,png}（Nature 尺寸，300+ dpi）
#
# Key conventions / 出图硬规则（沿用项目既定标准）:
#   - 统一 bbox：coord_sf(xlim=c(73,135), ylim=c(18,54), expand=FALSE)
#   - 不画十段线、不加鹰眼图/inset
#   - facet 中不得出现 NA 子图；缺值格点显式显示为灰色，不伪装成"无网格"
#   - 跨指标比较一律 within-metric z-score
#   - 色盲友好配色；PDF（矢量投稿）+ PNG（预览）双出
# Packages: ggplot2, sf, dplyr, readr, tidyr, patchwork, scico
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(dplyr); library(readr)
  library(tidyr); library(patchwork)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_mapping.R"))        # theme_nature_map / save_nature / 调色
source(file.path(CODE_V3, "utils_plots_advanced.R")) # forest_plot_nature 等

ensure_v3_dirs()
is_pilot  <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
log_time("34", "Publication figures for PA analyses")

res_path <- function(stem) paste0(v3_file("results", paste0(stem, "_", run_label)), ".csv")
read_if <- function(stem) { p <- res_path(stem); if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL }

# 配色：保护地内=蓝，外=橙（与全篇一致，色盲友好）
COL_IN  <- "#2E75B6"; COL_OUT <- "#D97A34"; COL_MID <- "#9AA7B4"
PA_COLS <- c(protected = COL_IN, partial = COL_MID, unprotected = COL_OUT)
PA_LABS <- c(protected = "保护地内 (≥30%)", partial = "部分保护", unprotected = "保护地外 (<5%)")

# 统一的指标显示名（避免图上出现变量名）
metric_label <- c(corrected_richness = "校正后物种数", shannon = "Shannon 多样性",
                  inv_simpson = "inverse Simpson", trait_volume = "功能体积",
                  rao_q = "Rao's Q", fric_prob = "功能丰富度", fdis_prob = "功能离散度",
                  feve_fund = "功能均匀度", pd_prob = "系统发育多样性")
lab_of <- function(x) ifelse(x %in% names(metric_label), metric_label[x], x)

# ── PA2. 多维多样性逐期轨迹：保护地内 vs 外（核心图）──────────────────
traj <- read_if("table_pa_trajectory")
if (!is.null(traj)) {
  focal <- intersect(unique(traj$metric),
    c("corrected_richness", "shannon", "trait_volume", "rao_q"))
  d <- traj |> filter(metric %in% focal, !is.na(pa_stratum)) |>
    mutate(metric_lab = factor(lab_of(metric), levels = lab_of(focal)),
           period_num = as.integer(sub("P", "", period)),
           pa_stratum = factor(pa_stratum, levels = names(PA_COLS)))

  p2 <- ggplot(d, aes(period_num, mean, colour = pa_stratum, fill = pa_stratum)) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    facet_wrap(~ metric_lab, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = PA_COLS, labels = PA_LABS, name = NULL) +
    scale_fill_manual(values = PA_COLS, labels = PA_LABS, name = NULL) +
    scale_x_continuous(breaks = 1:5,
      labels = c("2000-04", "2005-09", "2010-14", "2015-19", "2020-24")) +
    labs(x = NULL, y = "后验均值（阴影为标准误）") +
    theme_nature_pub() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
  save_nature(p2, paste0("fig_pa2_trajectories_", run_label), width_mm = 183, height_mm = 68)
  message("[34] PA2 轨迹图已输出")
}

# ── PA5. 空间 beta（同质化速率）内外对比 ─────────────────────────────
bet <- read_if("table_pa_beta_contrast")
if (!is.null(bet)) {
  d <- bet |>
    select(period, inside = spatial_beta_inside, outside = spatial_beta_outside,
           sd_in = spatial_beta_inside_sd, sd_out = spatial_beta_outside_sd) |>
    pivot_longer(c(inside, outside), names_to = "grp", values_to = "beta") |>
    mutate(sd = ifelse(grp == "inside", sd_in, sd_out),
           period_num = as.integer(sub("P", "", period)))

  p5 <- ggplot(d, aes(period_num, beta, colour = grp, fill = grp)) +
    geom_ribbon(aes(ymin = beta - sd, ymax = beta + sd), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    scale_colour_manual(values = c(inside = COL_IN, outside = COL_OUT),
                        labels = c("保护地内", "保护地外"), name = NULL) +
    scale_fill_manual(values = c(inside = COL_IN, outside = COL_OUT),
                      labels = c("保护地内", "保护地外"), name = NULL) +
    scale_x_continuous(breaks = 1:5,
      labels = c("2000-04", "2005-09", "2010-14", "2015-19", "2020-24")) +
    labs(x = NULL, y = "格点间空间 beta 多样性",
         subtitle = "下降即群落趋同（同质化）；斜率更陡的一组同质化更快") +
    theme_nature_pub() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
  save_nature(p5, paste0("fig_pa5_spatial_beta_", run_label), width_mm = 89, height_mm = 72)
  message("[34] PA5 空间 beta 对比图已输出")
}

# ── PA6. 特化种 vs 泛化种的占域分化 ──────────────────────────────────
sp <- read_if("table_pa_species_occupancy")
if (!is.null(sp) && "sp_group" %in% names(sp)) {
  d <- sp |> filter(!is.na(sp_group)) |>
    group_by(sp_group, period) |>
    summarise(diff = mean(psi_diff, na.rm = TRUE),
              se = sd(psi_diff, na.rm = TRUE) / sqrt(n()), .groups = "drop") |>
    mutate(period_num = as.integer(sub("P", "", period)),
           sp_group = factor(sp_group, levels = c("specialist", "generalist")))

  p6 <- ggplot(d, aes(period_num, diff, colour = sp_group, fill = sp_group)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = diff - se, ymax = diff + se), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    scale_colour_manual(values = c(specialist = "#7B3294", generalist = COL_OUT),
                        labels = c("特化种（窄栖息地）", "泛化种（广栖息地）"), name = NULL) +
    scale_fill_manual(values = c(specialist = "#7B3294", generalist = COL_OUT),
                      labels = c("特化种（窄栖息地）", "泛化种（广栖息地）"), name = NULL) +
    scale_x_continuous(breaks = 1:5,
      labels = c("2000-04", "2005-09", "2010-14", "2015-19", "2020-24")) +
    labs(x = NULL, y = "占域概率差（保护地内 − 外）",
         subtitle = "正值=保护地内占域更高；特化种差值扩大即保护区选择性维持特化种") +
    theme_nature_pub() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
  save_nature(p6, paste0("fig_pa6_specialist_generalist_", run_label), width_mm = 89, height_mm = 72)
  message("[34] PA6 特化/泛化种分化图已输出")
}

# ── PA4. 处理效应森林图（各维度 ATT）────────────────────────────────
att <- read_if("table_pa_matching_att")
if (!is.null(att)) {
  est <- intersect(c("att", "estimate", "ATT"), names(att))[1]
  se  <- intersect(c("se", "std_error", "SE"), names(att))[1]
  out <- intersect(c("outcome", "metric"), names(att))[1]
  if (!is.na(est) && !is.na(se) && !is.na(out)) {
    d <- att |> rename(estimate = all_of(est), se = all_of(se), outcome = all_of(out)) |>
      mutate(lwr = estimate - 1.96 * se, upr = estimate + 1.96 * se,
             sig = ifelse(lwr > 0 | upr < 0, "显著", "不显著"),
             lab = lab_of(sub("^trend_", "", outcome)))
    p4 <- ggplot(d, aes(estimate, reorder(lab, estimate), colour = sig)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.18, linewidth = 0.45) +
      geom_point(size = 2) +
      scale_colour_manual(values = c("显著" = COL_IN, "不显著" = "grey65"), name = NULL) +
      labs(x = "处理效应 ATT（保护地内 − 匹配对照）", y = NULL,
           subtitle = "误差棒为 95% 置信区间；正值=保护地内该维度表现更好") +
      theme_nature_pub() + theme(legend.position = "top")
    save_nature(p4, paste0("fig_pa4_att_forest_", run_label), width_mm = 89, height_mm = 70)
    message("[34] PA4 ATT 森林图已输出")
  }
}

# ── PA7. 事件研究图（动态 DiD，检验平行趋势）─────────────────────────
es <- read_if("table_pa_event_study")
if (!is.null(es)) {
  d <- es |> filter(grepl("^ett-?[0-9]|^et", term), term != "(Intercept)") |>
    mutate(event_time = suppressWarnings(as.integer(gsub("[^0-9-]", "", term)))) |>
    filter(!is.na(event_time)) |>
    mutate(lwr = estimate - 1.96 * se, upr = estimate + 1.96 * se,
           metric_lab = lab_of(metric))
  if (nrow(d) > 0) {
    p7 <- ggplot(d, aes(event_time, estimate)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
      geom_vline(xintercept = -0.5, linetype = 3, colour = COL_OUT, linewidth = 0.4) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.12, linewidth = 0.4, colour = COL_IN) +
      geom_point(size = 1.8, colour = COL_IN) +
      facet_wrap(~ metric_lab, scales = "free_y") +
      labs(x = "相对保护区设立的时期（0 = 设立当期）", y = "标准化效应",
           subtitle = "设立前系数应接近 0（平行趋势）；设立后偏离 0 即为保护效应") +
      theme_nature_pub()
    save_nature(p7, paste0("fig_pa7_event_study_", run_label), width_mm = 183, height_mm = 68)
    message("[34] PA7 事件研究图已输出")
  }
}

# ── PA1 & PA8. 地图类：保护梯度 与 保护规划方案 ──────────────────────
genv <- safe_read(v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds"))
if (is.null(genv)) genv <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
cov <- read_if("table_pa_grid_coverage")

if (!is.null(genv) && !is.null(cov)) {
  frac_col <- intersect(c("pa_fraction", "pa_frac", "pa_cover"), names(cov))[1]
  gmap <- genv |>
    select(grid_cell, lon = centroid_lon, lat = centroid_lat) |>
    left_join(cov |> select(grid_cell, pa_frac = all_of(frac_col)), by = "grid_cell")

  # PA1：保护区覆盖率（缺值显式灰色，不隐藏）
  p1 <- ggplot(gmap, aes(lon, lat, fill = pa_frac)) +
    geom_tile(width = 1.0, height = 1.0) +
    scale_fill_viridis_c(option = "mako", direction = -1, na.value = "grey92",
                         name = "保护区\n面积占比", labels = scales::percent) +
    coord_sf(xlim = c(73, 135), ylim = c(18, 54), expand = FALSE) +
    labs(x = NULL, y = NULL, subtitle = "灰色=无数据网格（显式标出，非缺省隐藏）") +
    theme_nature_map()
  save_nature(p1, paste0("fig_pa1_coverage_map_", run_label), width_mm = 120, height_mm = 96)
  message("[34] PA1 保护梯度地图已输出")

  # PA8：保护规划三方案对比（脚本 32 产出）
  sol <- read_if("table_priority_solutions")
  if (!is.null(sol)) {
    sc_cols <- setdiff(names(sol), "grid_cell")
    d <- sol |> pivot_longer(all_of(sc_cols), names_to = "scenario", values_to = "selected") |>
      filter(!is.na(selected)) |>
      left_join(gmap |> select(grid_cell, lon, lat), by = "grid_cell") |>
      mutate(selected = factor(ifelse(selected > 0, "选中", "未选中"),
                               levels = c("选中", "未选中")))
    p8 <- ggplot(d, aes(lon, lat, fill = selected)) +
      geom_tile(width = 1.0, height = 1.0) +
      facet_wrap(~ scenario) +
      scale_fill_manual(values = c("选中" = COL_IN, "未选中" = "grey90"),
                        na.value = "grey92", name = NULL) +
      coord_sf(xlim = c(73, 135), ylim = c(18, 54), expand = FALSE) +
      labs(x = NULL, y = NULL,
           subtitle = "三套优先区方案；若空间重叠低，说明按物种数规划会错过阻止同质化的关键区") +
      theme_nature_map() + theme(legend.position = "top")
    save_nature(p8, paste0("fig_pa8_prioritization_", run_label), width_mm = 183, height_mm = 84)
    message("[34] PA8 保护规划方案图已输出")
  }
}

# ── PA9. 保护年限梯度轨迹（剂量-反应）────────────────────────────────
age <- read_if("table_pa_age_trajectory")
if (!is.null(age)) {
  focal <- intersect(unique(age$metric), c("corrected_richness", "trait_volume", "rao_q"))
  AGE_COLS <- c("unprotected" = COL_OUT, "young(<15a)" = "#9AC3E0",
                "mid(15-30a)" = COL_IN, "old(>30a)" = "#12395E")
  d <- age |> filter(metric %in% focal, age_class %in% names(AGE_COLS)) |>
    mutate(metric_lab = factor(lab_of(metric), levels = lab_of(focal)),
           period_num = as.integer(sub("P", "", period)),
           age_class = factor(age_class, levels = names(AGE_COLS)))
  p9 <- ggplot(d, aes(period_num, mean, colour = age_class, fill = age_class)) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
    facet_wrap(~ metric_lab, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = AGE_COLS, name = "保护年限") +
    scale_fill_manual(values = AGE_COLS, name = "保护年限") +
    scale_x_continuous(breaks = 1:5,
      labels = c("2000-04", "2005-09", "2010-14", "2015-19", "2020-24")) +
    labs(x = NULL, y = "后验均值（阴影为标准误）",
         subtitle = "剂量-反应：若保护有效，保护年限越长的网格功能维度轨迹应越好") +
    theme_nature_pub() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
  save_nature(p9, paste0("fig_pa9_protection_age_", run_label), width_mm = 183, height_mm = 70)
  message("[34] PA9 保护年限梯度图已输出")
}

# ── PA3. 匹配协变量平衡图（love plot）──────────────────────────────
bal <- read_if("table_pa_matching_balance")
if (!is.null(bal)) {
  vcol <- intersect(c("variable", "covariate", "term"), names(bal))[1]
  d <- bal |> rename(variable = all_of(vcol)) |>
    pivot_longer(matches("smd|std_diff|before|after"),
                 names_to = "stage", values_to = "smd") |>
    filter(!is.na(smd)) |>
    mutate(stage = ifelse(grepl("before|un|raw", stage, ignore.case = TRUE), "匹配前", "匹配后"))
  p3 <- ggplot(d, aes(abs(smd), reorder(variable, abs(smd)), colour = stage)) +
    geom_vline(xintercept = 0.1, linetype = 2, colour = COL_OUT, linewidth = 0.35) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("匹配前" = "grey60", "匹配后" = COL_IN), name = NULL) +
    labs(x = "标准化均值差绝对值 |SMD|", y = NULL,
         subtitle = "虚线=0.1 平衡阈值；匹配后应全部落在阈值左侧") +
    theme_nature_pub() + theme(legend.position = "top")
  save_nature(p3, paste0("fig_pa3_balance_love_", run_label), width_mm = 89, height_mm = 80)
  message("[34] PA3 匹配平衡图已输出")
} else message("[34] 未找到匹配平衡表，跳过 PA3（需脚本 31 输出 table_pa_matching_balance）")

log_time("34", "DONE")
