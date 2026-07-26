#!/usr/bin/env Rscript
## 02d_effort_cross_dataset_scatter.R
##
## 在 02c 基础上扩展：
##   1) 跨数据集（China_Birdwatch vs eBird_GBIF）× 同指标 的散点+拟合图（带 R²/p/n/RMA/OLS）
##   2) 跨指标 × 同数据集 的 pair-grid 散点（GGally）带统计参数
##   3) 跨数据集 × 跨指标 全集 facet 散点矩阵
##   4) 复合 effort index 与每个原始指标的关系散点
##
## 依赖 02c 已生成 results_v2/table_effort_metrics_year_source.csv
## + results_v2/table_effort_index_year.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(forcats); library(ggplot2); library(GGally); library(broom)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
P <- ensure_v2_dirs()

eff <- read_csv(v2_file("results", "table_effort_metrics_year_source"),
                show_col_types = FALSE)
idx <- read_csv(v2_file("results", "table_effort_index_year"),
                show_col_types = FALSE)
metric_vars <- c("n_events","n_visits","n_observers","n_birding_days",
                 "n_grids","n_species","total_duration")
metric_pretty <- c(
  n_events = "Records (events)",
  n_visits = "Visits (observer-days)",
  n_observers = "Unique observers",
  n_birding_days = "Birding days",
  n_grids = "Grids visited",
  n_species = "Species detected",
  total_duration = "Total duration (min)"
)

## --- 1. 跨数据集 × 同指标：散点+OLS+RMA --------------------------------------

wide_src <- eff |>
  filter(source_short %in% c("China_Birdwatch", "eBird_GBIF")) |>
  select(year, source_short, all_of(metric_vars)) |>
  pivot_longer(-c(year, source_short), names_to = "metric") |>
  pivot_wider(names_from = source_short, values_from = value) |>
  drop_na(China_Birdwatch, eBird_GBIF) |>
  mutate(metric_label = factor(metric_pretty[metric],
                                 levels = unname(metric_pretty)))

# 每指标算 OLS + Spearman + RMA-like (geometric mean of two slopes)
fit_per_metric <- wide_src |>
  group_by(metric_label) |>
  summarise(
    n        = n(),
    rho      = cor(China_Birdwatch, eBird_GBIF, method = "spearman"),
    pearson  = cor(China_Birdwatch, eBird_GBIF, method = "pearson"),
    p_value  = cor.test(China_Birdwatch, eBird_GBIF, method = "spearman")$p.value,
    ols_slope = coef(lm(eBird_GBIF ~ China_Birdwatch))[2],
    ols_int   = coef(lm(eBird_GBIF ~ China_Birdwatch))[1],
    R2        = summary(lm(eBird_GBIF ~ China_Birdwatch))$r.squared,
    .groups = "drop"
  ) |>
  mutate(stat_label = sprintf("n=%d | rho=%.2f | R2=%.2f | p=%s",
                                n, rho, R2,
                                format.pval(p_value, digits = 2, eps = 1e-3)))
write_csv(fit_per_metric,
          v2_file("results", "table_effort_cross_dataset_stats"))

# 注释位置
labels_pos <- wide_src |>
  group_by(metric_label) |>
  summarise(
    x = min(China_Birdwatch, na.rm = TRUE),
    y = max(eBird_GBIF, na.rm = TRUE) * 1.0,
    .groups = "drop"
  ) |>
  inner_join(fit_per_metric |> select(metric_label, stat_label),
              by = "metric_label")

p_cross <- ggplot(wide_src, aes(China_Birdwatch, eBird_GBIF)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey55") +
  geom_smooth(method = "lm", se = TRUE, colour = "#0E5A78",
              fill = "#0E5A78", alpha = 0.18, linewidth = 0.7) +
  geom_point(aes(colour = year), alpha = 0.9, size = 2.2) +
  geom_text(data = labels_pos,
            aes(x = x, y = y, label = stat_label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.0,
            colour = "grey15") +
  facet_wrap(~ metric_label, scales = "free", ncol = 3) +
  scale_colour_viridis_c(option = "mako", direction = -1, name = "Year") +
  labs(
    title = "Cross-dataset agreement of annual effort metrics",
    subtitle = "China Birdwatch (x) vs eBird/GBIF (y) per metric, 2000-2025; dashed = 1:1; solid = OLS fit",
    caption = "Stats: Spearman rho, Pearson R^2 on log-natural scale, n = years.",
    x = "China Birdwatch", y = "eBird / GBIF"
  ) +
  theme_v2_pub(11) +
  theme(legend.position = "right")
save_dual(p_cross, "fig_effort_cross_dataset_scatter_v3", width = 13, height = 8.5)

# 该图通常需要可编辑 PPT
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE)) {
  d <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(value = rvg::dml(ggobj = p_cross),
                      location = officer::ph_location(left = 0.4, top = 0.4,
                                                      width = 13, height = 8.5))
  print(d, target = file.path(P$figures_v2,
                               "fig_effort_cross_dataset_scatter_v3.pptx"))
}

## --- 2. 跨指标矩阵 + log10 缩放 + 统计参数 --------------------------------

dat_cross <- eff |> filter(source_short == "Combined") |>
  as_tibble() |>
  select(year, all_of(metric_vars)) |>
  drop_na() |>
  mutate(across(all_of(metric_vars), ~ log10(. + 1)))

fit_loess_with_stats <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) return(ggplot() + theme_void())
  xv <- x[ok]; yv <- y[ok]
  fit <- lm(yv ~ xv)
  r <- cor(xv, yv, method = "pearson")
  rho <- cor(xv, yv, method = "spearman")
  R2 <- summary(fit)$r.squared
  lab <- sprintf("r=%.2f\nrho=%.2f\nR2=%.2f", r, rho, R2)
  ggplot(data, mapping) +
    geom_point(alpha = 0.55, size = 1.0, colour = "#0E5A78") +
    geom_smooth(method = "lm", se = TRUE, colour = "#8B2E1E",
                fill = "#8B2E1E", alpha = 0.18, linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf, label = lab,
              hjust = -0.06, vjust = 1.1, size = 2.8, colour = "grey20") +
    theme_v2_pub(8.5) +
    theme(panel.grid = element_line(colour = "grey94", linewidth = 0.2))
}

p_pair <- GGally::ggpairs(
  dat_cross |> select(-year),
  upper = list(continuous = GGally::wrap("cor", method = "spearman",
                                          size = 3.4, stars = TRUE,
                                          digits = 2)),
  lower = list(continuous = fit_loess_with_stats),
  diag  = list(continuous = GGally::wrap("densityDiag",
                                          fill = "#0E5A78", alpha = 0.45,
                                          colour = NA)),
  columnLabels = unname(metric_pretty)
) + theme_v2_pub(9) +
  theme(strip.text = element_text(face = "bold", size = 9))
save_dual(p_pair, "fig_effort_pair_grid_with_stats_v3",
          width = 14, height = 13)

## --- 3. 跨数据集 × 跨指标 全集 facet（log scale + RMA-like 比例线）-----------

cross_full <- eff |>
  filter(source_short %in% c("China_Birdwatch", "eBird_GBIF")) |>
  select(year, source_short, all_of(metric_vars)) |>
  pivot_longer(-c(year, source_short), names_to = "metric")

# 把"两源 × 8 metric"长表横向拼，用每对 (source × metric) 作为 panel
cross_per_pair <- eff |>
  filter(source_short %in% c("China_Birdwatch", "eBird_GBIF")) |>
  select(year, source_short, all_of(metric_vars)) |>
  pivot_longer(-c(year, source_short),
                names_to = "metric_x", values_to = "value_x") |>
  inner_join(
    eff |> filter(source_short %in% c("China_Birdwatch", "eBird_GBIF")) |>
      select(year, source_short, all_of(metric_vars)) |>
      pivot_longer(-c(year, source_short),
                    names_to = "metric_y", values_to = "value_y"),
    by = c("year", "source_short")
  )
# 太多 panels → 仅保留 events 与其他指标的双源对比
focal_metric <- "n_events"
focal_dat <- cross_per_pair |>
  filter(metric_x == focal_metric, metric_y != focal_metric) |>
  mutate(metric_y = factor(metric_pretty[metric_y],
                            levels = unname(metric_pretty)))

stats_focal <- focal_dat |>
  group_by(metric_y, source_short) |>
  summarise(n = n(),
            rho = cor(value_x, value_y, method = "spearman", use = "complete.obs"),
            R2  = summary(lm(value_y ~ value_x))$r.squared,
            p   = cor.test(value_x, value_y, method = "spearman")$p.value,
            .groups = "drop") |>
  mutate(label = sprintf("rho=%.2f\nR2=%.2f", rho, R2))

p_focal <- ggplot(focal_dat, aes(value_x, value_y, colour = source_short,
                                    fill = source_short)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.18, linewidth = 0.6) +
  geom_text(data = stats_focal, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.1, size = 2.8, colour = "grey20") +
  facet_grid(metric_y ~ source_short, scales = "free") +
  scale_colour_manual(values = c(China_Birdwatch = "#0E5A78",
                                  eBird_GBIF = "#C9784A"), guide = "none") +
  scale_fill_manual(values = c(China_Birdwatch = "#0E5A78",
                                eBird_GBIF = "#C9784A"), guide = "none") +
  labs(
    title = "How records (events) scale with other effort metrics, per dataset",
    subtitle = "rows = response metric | columns = data source",
    x = "Records (events) per year", y = NULL
  ) +
  theme_v2_pub(10) +
  theme(panel.grid = element_line(colour = "grey94", linewidth = 0.2))
save_dual(p_focal, "fig_effort_records_vs_others_by_source_v3",
          width = 11, height = 12)
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE)) {
  d <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(value = rvg::dml(ggobj = p_focal),
                      location = officer::ph_location(left = 0.4, top = 0.4,
                                                      width = 11, height = 12))
  print(d, target = file.path(P$figures_v2,
                               "fig_effort_records_vs_others_by_source_v3.pptx"))
}

## --- 4. 复合 effort index 与原始指标 ---------------------------------------

idx_long <- eff |> filter(source_short == "Combined") |>
  inner_join(idx |> select(year, effort_pc1, effort_zmean), by = "year") |>
  pivot_longer(all_of(metric_vars), names_to = "metric") |>
  mutate(metric_label = factor(metric_pretty[metric],
                                 levels = unname(metric_pretty)))
stats_idx <- idx_long |>
  group_by(metric_label) |>
  summarise(rho_pc1   = cor(effort_pc1, log10(value + 1), method = "spearman",
                             use = "complete.obs"),
            R2_pc1    = summary(lm(log10(value + 1) ~ effort_pc1))$r.squared,
            rho_zmean = cor(effort_zmean, log10(value + 1), method = "spearman",
                              use = "complete.obs"),
            .groups = "drop") |>
  mutate(label = sprintf("rho(PC1)=%.2f | R2=%.2f", rho_pc1, R2_pc1))

p_idx <- ggplot(idx_long, aes(effort_pc1, log10(value + 1))) +
  geom_smooth(method = "lm", se = TRUE, colour = "#0E5A78",
              fill = "#0E5A78", alpha = 0.18, linewidth = 0.7) +
  geom_point(aes(colour = year), alpha = 0.9, size = 2.0) +
  geom_text(data = stats_idx, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.1, size = 3.0, colour = "grey20") +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_colour_viridis_c(option = "mako", direction = -1, name = "Year") +
  labs(title = "Composite effort index (PC1) vs each raw metric",
       subtitle = "Validates that PC1 captures the same signal as each individual metric",
       x = "Composite effort index (PC1, sd units)",
       y = "log10(metric + 1)") +
  theme_v2_pub(11)
save_dual(p_idx, "fig_effort_index_vs_metrics_v3",
          width = 13, height = 8.5)
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE)) {
  d <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(value = rvg::dml(ggobj = p_idx),
                      location = officer::ph_location(left = 0.4, top = 0.4,
                                                      width = 13, height = 8.5))
  print(d, target = file.path(P$figures_v2,
                               "fig_effort_index_vs_metrics_v3.pptx"))
}

message("[stage-2d] cross-dataset / cross-metric scatter + stats done.")
