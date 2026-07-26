#!/usr/bin/env Rscript
## 02c_effort_index_correlation.R
##
## 调查努力指标相关性分析 + 复合 effort index 构建。
## 输入：data/derived_v2/combined_events_merged_dedup_2000_2025.rds
## 输出：results_v2/table_effort_metrics_year_source.csv（年度多指标）
##         results_v2/table_effort_correlation_matrix.csv（各 metric × source 的 Spearman/Pearson）
##         results_v2/table_effort_index_year.csv（PCA + z-mean 复合指数）
##         figures_v2/fig_effort_corr_heatmap_v2.{png,pdf}
##         figures_v2/fig_effort_pairwise_scatter_v2.{png,pdf}
##         figures_v2/fig_effort_pca_biplot_v2.{png,pdf}
##         figures_v2/fig_effort_composite_index_timeseries_v2.{png,pdf}

suppressPackageStartupMessages({
  library(data.table); library(readr); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(GGally); library(forcats)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))

P <- ensure_v2_dirs()

events <- as.data.table(readRDS(file.path(P$derived_v2,
                          "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= 2000 & year <= 2025]
events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]
events[, day_key := sprintf("%04d-%02d-%02d", year, month, day)]
events[, observer_key := tolower(coalesce(username, "anon"))]
events[, visit_key := paste(observer_key, day_key, sep = "|")]
events[, has_dur := !is.na(duration_min)]

## --- 1. 年度 × 源 多指标 ---------------------------------------------------

annual_per_source <- events[, .(
  n_events       = .N,                         # 原始记录数
  n_visits       = uniqueN(visit_key),         # 观测次数（unique observer-day）
  n_observers    = uniqueN(observer_key),      # 观测人数
  n_birding_days = uniqueN(day_key),           # 观鸟天数（独立日历日）
  n_grids        = uniqueN(paste(round(longitude,1), round(latitude,1))),
  n_species      = uniqueN(species),
  total_duration = sum(duration_min[has_dur], na.rm = TRUE)
), by = .(year, source_short)]

annual_combined <- events[, .(
  n_events       = .N,
  n_visits       = uniqueN(visit_key),
  n_observers    = uniqueN(observer_key),
  n_birding_days = uniqueN(day_key),
  n_grids        = uniqueN(paste(round(longitude,1), round(latitude,1))),
  n_species      = uniqueN(species),
  total_duration = sum(duration_min[has_dur], na.rm = TRUE)
), by = year][, source_short := "Combined"]

annual_all <- rbindlist(list(annual_per_source, annual_combined),
                         use.names = TRUE)
annual_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
write_csv(as_tibble(annual_all),
          v2_file("results", "table_effort_metrics_year_source"))

## --- 2. 相关性矩阵 ---------------------------------------------------------

metric_vars <- c("n_events", "n_visits", "n_observers", "n_birding_days",
                  "n_grids", "n_species", "total_duration")

corr_long <- annual_all |>
  as_tibble() |>
  group_by(source_short) |>
  group_modify(~ {
    M <- as.matrix(.x[, metric_vars])
    sp <- cor(M, use = "pairwise.complete.obs", method = "spearman")
    pe <- cor(M, use = "pairwise.complete.obs", method = "pearson")
    tibble(
      metric_x = rep(metric_vars, each = length(metric_vars)),
      metric_y = rep(metric_vars, times = length(metric_vars)),
      spearman = as.vector(sp),
      pearson  = as.vector(pe)
    )
  }) |>
  ungroup()
write_csv(corr_long, v2_file("results", "table_effort_correlation_matrix"))

heat_dat <- corr_long |>
  filter(source_short == "Combined") |>
  mutate(metric_x = factor(metric_x, levels = metric_vars),
         metric_y = factor(metric_y, levels = rev(metric_vars)))

heat_plot <- ggplot(heat_dat, aes(metric_x, metric_y, fill = spearman)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman)),
            size = 3.4, colour = "grey15") +
  scale_fill_v2_diverging(name = "Spearman rho",
                          limits = c(-1, 1), midpoint = 0) +
  facet_wrap(~ source_short) +
  labs(title = "Effort metric correlations (Spearman)",
       subtitle = "Annual time series 2000-2025 in the merged + deduplicated dataset",
       x = NULL, y = NULL,
       caption = "Cells show pairwise correlation across years.") +
  theme_v2_pub(11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid = element_blank())
save_dual(heat_plot, "fig_effort_corr_heatmap_v2",
          width = 8, height = 7)

# 三源拼一张
heat_all <- corr_long |>
  mutate(metric_x = factor(metric_x, levels = metric_vars),
         metric_y = factor(metric_y, levels = rev(metric_vars)))
heat3 <- ggplot(heat_all, aes(metric_x, metric_y, fill = spearman)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", spearman)),
            size = 2.8, colour = "grey15") +
  scale_fill_v2_diverging(name = "Spearman rho",
                          limits = c(-1, 1), midpoint = 0) +
  facet_wrap(~ source_short, ncol = 3) +
  labs(title = "Effort metric correlations across sources",
       subtitle = "Spearman correlation among annual metrics, per data source",
       x = NULL, y = NULL) +
  theme_v2_pub(10.4) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid = element_blank())
save_dual(heat3, "fig_effort_corr_heatmap_3sources_v2",
          width = 16, height = 6)

## --- 3. 配对散点 + 拟合 ---------------------------------------------------

pair_dat <- annual_all |>
  filter(source_short == "Combined") |>
  as_tibble() |>
  select(year, all_of(metric_vars))

pair_plot <- GGally::ggpairs(
  pair_dat |> select(-year),
  upper = list(continuous = GGally::wrap("cor", method = "spearman",
                                          size = 3.6, stars = FALSE)),
  lower = list(continuous = GGally::wrap("smooth", method = "lm",
                                          se = TRUE,
                                          colour = "#0E5A78",
                                          alpha = 0.6, size = 0.8)),
  diag  = list(continuous = GGally::wrap("densityDiag",
                                          fill = "#0E5A78", alpha = 0.45,
                                          colour = NA))
) + theme_v2_pub(9.5) +
  theme(panel.grid = element_line(colour = "grey94", linewidth = 0.2),
        strip.text = element_text(face = "bold"))
save_dual(pair_plot, "fig_effort_pairwise_scatter_v2",
          width = 12, height = 11)

## --- 4. PCA 复合 effort index --------------------------------------------

# 在 Combined 行上 fit PCA（用 log1p 缩平偏态）
mat_combined <- annual_all |>
  filter(source_short == "Combined") |>
  as_tibble() |>
  arrange(year)
X <- mat_combined |> select(all_of(metric_vars)) |>
  mutate(across(everything(), ~ log1p(as.numeric(.x))))
X_clean <- X[apply(X, 1, function(r) all(is.finite(r))), ]
year_clean <- mat_combined$year[apply(X, 1, function(r) all(is.finite(r)))]

pca <- prcomp(X_clean, center = TRUE, scale. = TRUE)
expl <- summary(pca)$importance[2, ]
loadings <- as_tibble(pca$rotation, rownames = "metric")
write_csv(loadings, v2_file("results", "table_effort_pca_loadings"))

# PC1 作为复合 effort index（保证主要量随效率单调上升）
pc1 <- pca$x[, 1]
if (cor(pc1, X_clean$n_events, use = "complete.obs") < 0) pc1 <- -pc1
zmean <- rowMeans(scale(X_clean))      # 备用：标准化均值法

idx_tbl <- tibble(year = year_clean,
                  effort_pc1 = as.numeric(pc1),
                  effort_zmean = as.numeric(zmean),
                  pc1_var_explained = expl["PC1"],
                  pc2_var_explained = expl["PC2"]) |>
  mutate(across(c(effort_pc1, effort_zmean),
                ~ as.numeric(scale(.x))))
write_csv(idx_tbl, v2_file("results", "table_effort_index_year"))

# Biplot
bp_dat <- as_tibble(pca$x[, 1:2]) |> mutate(year = year_clean)
load_dat <- loadings |> select(metric, PC1, PC2) |>
  mutate(PC1 = PC1 * max(abs(bp_dat$PC1)),
         PC2 = PC2 * max(abs(bp_dat$PC2)))
biplot <- ggplot(bp_dat, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey75") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey75") +
  geom_path(aes(group = 1), colour = "grey55", linewidth = 0.4,
             arrow = grid::arrow(length = grid::unit(0.18, "cm"),
                                  type = "closed")) +
  geom_point(aes(colour = year), size = 3.0) +
  geom_text(aes(label = year), nudge_y = 0.18, size = 3.0, colour = "grey25") +
  geom_segment(data = load_dat,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               colour = "#8B2E1E", linewidth = 0.55,
               arrow = grid::arrow(length = grid::unit(0.18, "cm"))) +
  geom_text(data = load_dat,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = metric),
            inherit.aes = FALSE,
            colour = "#8B2E1E", size = 3.4, fontface = "italic") +
  scale_colour_viridis_c(option = "mako", direction = -1, name = "Year") +
  labs(title = "PCA of annual effort metrics (merged + deduplicated)",
       subtitle = sprintf("PC1 explains %.1f%%, PC2 %.1f%% of variance",
                          100 * expl["PC1"], 100 * expl["PC2"]),
       caption = "Each point = one year, trajectory shows how effort changed over time. Red arrows = metric loadings.",
       x = sprintf("PC1 (%.1f%%)", 100 * expl["PC1"]),
       y = sprintf("PC2 (%.1f%%)", 100 * expl["PC2"])) +
  theme_v2_pub(11)
save_dual(biplot, "fig_effort_pca_biplot_v2",
          width = 9, height = 7)

# 复合 index 时序
idx_long <- idx_tbl |>
  pivot_longer(c(effort_pc1, effort_zmean),
                names_to = "method", values_to = "value") |>
  mutate(method = recode(method,
    effort_pc1   = "PC1 (PCA on log1p)",
    effort_zmean = "z-mean (mean of standardised metrics)"))
idx_plot <- ggplot(idx_long, aes(year, value, colour = method)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey75") +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("PC1 (PCA on log1p)" = "#0E5A78",
                                  "z-mean (mean of standardised metrics)" = "#C9784A"),
                       name = NULL) +
  scale_x_continuous(breaks = seq(2000, 2024, 4)) +
  labs(title = "Composite annual birding-effort index for the merged dataset",
       subtitle = "PCA-based PC1 vs z-score mean — both standardised; 0 = mean across years",
       x = NULL, y = "Standardised effort (sd units, year mean = 0)",
       caption = "Higher = more cumulative birding effort that year.") +
  theme_v2_pub(11)
save_dual(idx_plot, "fig_effort_composite_index_timeseries_v2",
          width = 10, height = 5.4)

message("[stage-2c] Effort correlation + composite index done.")
