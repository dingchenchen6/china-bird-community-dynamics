#!/usr/bin/env Rscript
# ============================================================
# 34c_pa_balance_and_maps_20261011.R  (sourced by 34b_pa_figures_english_standard_map_20261011.R)
#
# Scientific question / 科学问题:
#   PA3：倾向值匹配是否使保护区内外网格在基线协变量上可比（love plot）；
#   PA1：各网格的保护区覆盖率；PA8：三种规划目标在两种面积预算下选中的优先网格。
#   PA3 covariate balance; PA1 protected-area cover; PA8 priority networks under three objectives.
#
# Objective / 分析目标:
#   PA3：只画真实的匹配前/后 |SMD|（修复原脚本把逻辑列 balanced_after 画成 |SMD| = 1.00 的 bug），
#        协变量名改为英文全称，按匹配前 |SMD| 排序；
#   PA1、PA8：改用真实 100 km 网格多边形（原图以经纬度质心画 1°×1° 方块，留下条状空隙），
#        叠加省界、国界与南海诸岛插图（十段线），模型之外的网格显示为灰色 No data。
#   PA3: plot only the real SMD columns; PA1/PA8: standard map with true grid polygons.
#
# Input data / 输入数据: table_pa_matching_balance, table_pa_grid_coverage（S12），
#   table_priority_solutions（S13）; data/derived_v2/china_grid_100km_v2.rds
# Main workflow / 主要流程: 选列 → 作图 → 导出；地图：网格 + 标准地图图层 → 导出
# Expected output / 预期输出: protected_area/fig_pa{3,1,8}_*; qa/fig_pa3_plotted_smd.csv
# Key assumptions / 关键假设: "distance" 为 MatchIt 的倾向值（distance = "glm"，见脚本 31）。
# Main packages / 主要包: ggplot2, sf, dplyr, tidyr
# Output directory / 输出路径: figures_community_submission_20261011/{protected_area,qa}
# ============================================================

# ── PA3 匹配协变量平衡 / covariate balance (love plot) ───────────────────
COV_LABS <- c(distance = "Propensity score", elev_mean = "Elevation (mean)", elev_sd = "Elevation (s.d.)",
              bio1 = "Annual mean temperature (BIO1)", bio12 = "Annual precipitation (BIO12)",
              hfi_2000 = "Human footprint, 2000", landcover_trees = "Tree cover",
              landcover_cropland = "Cropland cover", landcover_built = "Built-up cover",
              npp_mean = "Net primary productivity", centroid_lon = "Longitude", centroid_lat = "Latitude")
cov_lab <- function(v) unname(ifelse(v %in% names(COV_LABS), COV_LABS[v], v))
bal <- read_si("S12_protected_areas", "table_pa_matching_balance")
vcol <- intersect(c("variable", "covariate", "term"), names(bal))[1]
bal <- bal |> rename(variable = all_of(vcol))
# 修复：只取 smd_before / smd_after 两列；原 matches("smd|std_diff|before|after") 还匹配到 balanced_after
# Fix: only the two real SMD columns; the old regex also matched the logical column balanced_after.
d3 <- bal |> select(variable, smd_before, smd_after) |>
  pivot_longer(c(smd_before, smd_after), names_to = "stage", values_to = "smd") |>
  filter(!is.na(smd)) |>
  mutate(stage = factor(ifelse(stage == "smd_before", "Before matching", "After matching"),
                        levels = c("Before matching", "After matching")),
         var_lab = factor(cov_lab(variable), levels = cov_lab(bal$variable[order(abs(bal$smd_before))])))
write_csv(d3 |> transmute(variable, covariate = as.character(var_lab), stage, smd, abs_smd = abs(smd)),
          file.path(QA_DIR, "fig_pa3_plotted_smd.csv"))
p3 <- ggplot(d3, aes(abs(smd), var_lab, colour = stage)) +
  geom_vline(xintercept = 0.1, linetype = 2, colour = COL_OUT, linewidth = 0.35) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Before matching" = "grey60", "After matching" = COL_IN), name = NULL) +
  labs(x = "Absolute standardized mean\ndifference, |SMD|", y = NULL,   # 两行，避免截断 / two lines, no clipping
       subtitle = "Dashed line: |SMD| = 0.1 balance threshold") +
  theme_nature_pub() + theme_pa_extra + theme(legend.position = "top")
save_pa(p3, "fig_pa3_balance_love", 89, 80)

# ── 地图公共部分 / shared map layers ─────────────────────────────────
safe_make_valid <- function(x) {
  out <- try(sf::st_make_valid(x), silent = TRUE)
  if (inherits(out, "try-error")) suppressWarnings(sf::st_buffer(x, 0)) else out
}
cov <- read_si("S12_protected_areas", "table_pa_grid_coverage")
frac_col <- intersect(c("pa_fraction", "pa_frac", "pa_cover"), names(cov))[1]   # 同脚本 34 / as in script 34
grid_sf <- readRDS(file.path(PROJECT_ROOT, "data", "derived_v2", "china_grid_100km_v2.rds")) |>
  safe_make_valid() |> st_transform(4326) |> select(grid_cell, geometry)
CHN <- china_layers_pub()
grid_plot <- grid_for_maps(grid_sf, model_cells = cov$grid_cell)
if (!all(cov$grid_cell %in% grid_plot$grid_cell)) stop("[34c] coverage cells missing from the grid")
SEAM_LW <- 0.08
pa_map_frame <- function(p, inset) {
  # 省界、国界、南海插图、固定 bbox；坐标轴标注沿用原 PA 地图风格
  # Provinces, outline, SCS inset, fixed bbox; axis labels as in the original PA maps.
  p + geom_sf(data = CHN$prov, fill = NA, colour = V5$MAP_LINE, linewidth = 0.1) +
    geom_sf(data = CHN$outline, fill = NA, colour = "grey20", linewidth = 0.25) +
    scs_inset_layer_pub(inset) +
    coord_sf(xlim = MAP_BB[1:2], ylim = MAP_BB[3:4], expand = FALSE)
}

# ── PA1 保护区覆盖率 / protected-area cover ─────────────────────────────
gmap <- grid_plot |> left_join(cov |> select(grid_cell, pa_frac = all_of(frac_col)), by = "grid_cell")
pa1_fill <- geom_sf(data = gmap, aes(fill = pa_frac, colour = after_scale(fill)), linewidth = SEAM_LW)
pa1_scale <- scale_fill_viridis_c(option = "mako", direction = -1, na.value = NA_GREY,
                                  name = "Protected-area\ncover", labels = scales::percent)
inset1 <- scs_inset_pub(CHN, line_col = "grey45", line_lw = 0.1, fill_layers = list(pa1_fill, pa1_scale))
p1 <- pa_map_frame(ggplot() + pa1_fill, inset1) + pa1_scale +
  labs(x = NULL, y = NULL, subtitle = "Grey cells: no data (shown explicitly, not hidden)") +
  theme_nature_map() + theme_pa_extra
save_pa(p1, "fig_pa1_coverage_map", 120, 100)

# ── PA8 三种规划目标 × 两种预算 / priority networks ──────────────────────
sol <- read_si("S13_prioritization", "table_priority_solutions")
PLAN_LABS <- c(A_richness_led = "Plan A: richness-led", B_functional_distinctiveness = "Plan B: functional distinctiveness",
               C_homogenization_risk = "Plan C: homogenization risk")
plan_lab <- function(p) unname(ifelse(p %in% names(PLAN_LABS), PLAN_LABS[p], p))
budget_lab <- function(b) paste0("Budget: ", round(100 * b), "% of cells")
plans <- unique(sol$plan); budgets <- sort(unique(sol$budget))
STATUS_COLS <- c("Selected" = COL_IN, "Not selected" = "grey90", "No data" = NA_GREY)
d8 <- tidyr::expand_grid(plan = plans, budget = budgets, grid_cell = grid_plot$grid_cell) |>
  left_join(sol |> filter(!is.na(selected)), by = c("plan", "budget", "grid_cell")) |>
  mutate(status = factor(case_when(is.na(selected) ~ "No data", selected > 0 ~ "Selected", TRUE ~ "Not selected"),
                         levels = names(STATUS_COLS)),
         plan_lab = factor(plan_lab(plan), levels = plan_lab(intersect(names(PLAN_LABS), plans))),
         budget_lab = factor(budget_lab(budget), levels = budget_lab(budgets)))
# 核对：每个方案×预算的选中格数与输入表一致 / check selected counts against the input table
chk <- d8 |> count(plan, budget, status) |> filter(status == "Selected") |>
  full_join(sol |> group_by(plan, budget) |> summarise(n_in = sum(selected > 0, na.rm = TRUE), .groups = "drop"),
            by = c("plan", "budget"))
if (any(chk$n != chk$n_in)) stop("[34c] PA8 selected-cell counts differ from the input table")
d8_sf <- grid_plot |> inner_join(d8, by = "grid_cell")
d8_px <- cell_pixels(grid_plot, res_deg = 0.0625) |>   # 1/16 度，二进制精确 / exact binary fraction
  inner_join(d8 |> select(grid_cell, plan_lab, budget_lab, status), by = "grid_cell", relationship = "many-to-many")
inset8 <- scs_inset_pub(CHN, line_col = "grey45", line_lw = 0.1)
p8_build <- function(fill_layer) {
  pa_map_frame(ggplot() + fill_layer, inset8) +
    facet_grid(budget_lab ~ plan_lab) +
    scale_fill_manual(values = STATUS_COLS, name = NULL, drop = FALSE) +
    labs(x = NULL, y = NULL,
         subtitle = wrap_sub("Three priority-area plans; low spatial overlap would mean that richness-led planning misses areas key to preventing homogenization", 183)) +
    theme_nature_map() + theme_pa_extra + theme(legend.position = "top")
}
p8 <- p8_build(geom_sf(data = d8_sf, aes(fill = status, colour = after_scale(fill)), linewidth = SEAM_LW))
p8_ppt <- p8_build(geom_raster(data = d8_px, aes(x, y, fill = status)))
save_pa(p8, "fig_pa8_prioritization", 183, 101, ppt_plot = p8_ppt, mode = "vector_with_raster_cell_fill")
message("[34c] PA8 selected-cell counts verified against the input table: ", nrow(chk), " plan x budget combinations")
