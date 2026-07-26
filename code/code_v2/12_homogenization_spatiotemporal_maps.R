#!/usr/bin/env Rscript
## 12_homogenization_spatiotemporal_maps.R
##
## 把生物同质化的空间-时间动态以地图展示：
##   - 每格 × 每期 的"局部组成异质性"= 该格 psi 向量与最近 K 个邻居 psi 向量
##     的平均 probability-Sørensen 距离。值低 → 与邻居相似（已同质化）；
##     值高 → 与邻居差异大（仍异质）。
##   - 每格的 5 期线性趋势：负 = 越来越像邻居（同质化加深）；正 = 异质化。
##
## 默认 V2_RUN_LABEL=v2_full_200sp_ar1，K_NEIGHBORS=20。
## 输出：
##   results_v2/table_homogenization_grid_period_<LABEL>.csv
##   figures_v2/fig_homogenization_timeslices_<LABEL>.{png,pdf}
##   figures_v2/fig_homogenization_trend_<LABEL>.{png,pdf}
##   figures_v2/fig_homogenization_combined_<LABEL>.{png,pdf}

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(forcats); library(ggplot2); library(patchwork); library(sf)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")
K_NEIGHBORS <- as.integer(Sys.getenv("V2_K_NEIGHBORS", "20"))

message(sprintf("[stage-12] spatiotemporal homogenization maps | %s | K=%d",
                RUN_LABEL, K_NEIGHBORS))

## --- 1. 加载 + 用 posterior mean psi（per-draw 太重） ---------------------

psi_obj <- readRDS(v2_file("derived",
                paste0("psi_samples_thinned_", RUN_LABEL), "rds"))
psi_arr <- psi_obj$psi_samples_thinned
n_sp    <- dim(psi_arr)[2]
n_site  <- dim(psi_arr)[3]
n_per   <- dim(psi_arr)[4]
psi_mean <- apply(psi_arr, c(2, 3, 4), mean)   # [sp × site × period]

primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
visit_effort <- readRDS(file.path(P$derived_v2, "visit_effort_2000_2024.rds"))
sites <- sort(unique(visit_effort$grid_cell))
grid_sf <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
grid_centroids <- grid_sf |>
  filter(grid_cell %in% sites) |>
  arrange(match(grid_cell, sites)) |>
  st_drop_geometry() |>
  select(grid_cell, centroid_lon, centroid_lat)
stopifnot(nrow(grid_centroids) == n_site)

## --- 2. 邻居索引：基于经纬度欧氏距离 -------------------------------------

coords <- as.matrix(grid_centroids[, c("centroid_lon", "centroid_lat")])
dist_xy <- as.matrix(dist(coords))
diag(dist_xy) <- Inf
neighbor_idx <- t(apply(dist_xy, 1, function(d) order(d)[seq_len(K_NEIGHBORS)]))

## --- 3. probability-Sørensen 距离（向量化 helper） -----------------------

sor_prob <- function(p1, p2) {
  ok <- is.finite(p1) & is.finite(p2)
  p1 <- p1[ok]; p2 <- p2[ok]
  A <- sum(pmin(p1, p2))
  B <- sum(pmax(p1 - p2, 0))
  C <- sum(pmax(p2 - p1, 0))
  if (2 * A + B + C <= 1e-12) NA_real_ else (B + C) / (2 * A + B + C)
}

## --- 4. 每格每期：本格 vs K 邻居 平均 Sørensen ---------------------------

res_arr <- matrix(NA_real_, n_site, n_per,
                   dimnames = list(as.character(sites),
                                   primary_blocks$block_label))
for (t in seq_len(n_per)) {
  P_t <- psi_mean[, , t]            # [sp × site]
  for (i in seq_len(n_site)) {
    nb <- neighbor_idx[i, ]
    p_i <- P_t[, i]
    vals <- vapply(nb, function(j) sor_prob(p_i, P_t[, j]), numeric(1))
    res_arr[i, t] <- mean(vals, na.rm = TRUE)
  }
  message(sprintf("  period %d/%d done", t, n_per))
}

## --- 5. per-grid 线性趋势 -------------------------------------------------

period_idx <- seq_len(n_per) - mean(seq_len(n_per))
denom <- sum(period_idx^2)
slope_vec <- apply(res_arr, 1, function(y) {
  ok <- is.finite(y); if (sum(ok) < 3) return(NA_real_)
  sum(period_idx[ok] * (y[ok] - mean(y[ok]))) / sum(period_idx[ok]^2)
})

mean_local <- rowMeans(res_arr, na.rm = TRUE)

## --- 6. 落表 -------------------------------------------------------------

long_tbl <- as.data.frame.table(res_arr, responseName = "mean_local_sorensen",
                                  stringsAsFactors = FALSE)
names(long_tbl)[1:2] <- c("grid_cell", "block_label")
long_tbl$grid_cell <- as.integer(as.character(long_tbl$grid_cell))
trend_tbl <- tibble(grid_cell = sites,
                     slope_per_period = slope_vec,
                     mean_local_sorensen = mean_local)
write_csv(long_tbl,
          v2_file("results",
                  paste0("table_homogenization_grid_period_", RUN_LABEL)))
write_csv(trend_tbl,
          v2_file("results",
                  paste0("table_homogenization_grid_trend_", RUN_LABEL)))

## --- 7. 出图 ------------------------------------------------------------

china_layers <- load_china_layers(P$china_boundary, P$province_line)

# 7a. time-slice
sf_long <- grid_sf |>
  inner_join(long_tbl |>
              mutate(block_label = factor(block_label,
                                            levels = primary_blocks$block_label)),
              by = "grid_cell")
p_ts <- ggplot(sf_long) +
  geom_sf(aes(fill = mean_local_sorensen),
          colour = scales::alpha("white", 0.08), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_v2_sequential(name = "Mean Sorensen to K-nearest neighbours",
                           palette = "lajolla", direction = 1,
                           limits = quantile(sf_long$mean_local_sorensen,
                                              c(0.02, 0.98), na.rm = TRUE),
                           oob = scales::squish) +
  facet_wrap(~ block_label, ncol = 5) +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = "Local biotic heterogeneity through time",
       subtitle = sprintf("Each grid: mean probability-Sorensen to its %d nearest neighbours per 5-year period",
                          K_NEIGHBORS),
       caption = "Lower (lighter) = more similar to neighbours (homogenised).  Higher (darker) = locally distinct.")
save_dual(p_ts, paste0("fig_homogenization_timeslices_", RUN_LABEL),
          width = 16, height = 5)

# 7b. trend
sf_trend <- grid_sf |> inner_join(trend_tbl, by = "grid_cell")
slope_lim <- quantile(abs(sf_trend$slope_per_period), 0.97, na.rm = TRUE)
p_trend <- ggplot(sf_trend) +
  geom_sf(aes(fill = slope_per_period),
          colour = scales::alpha("white", 0.10), linewidth = 0.04) +
  china_map_layers(china_layers) +
  scale_fill_v2_diverging(name = "Slope of local Sorensen vs period",
                          limits = c(-slope_lim, slope_lim),
                          oob = scales::squish) +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = "Where biotic homogenization is strongest",
       subtitle = "Per-grid linear trend of local Sorensen across the five 5-year periods",
       caption = "Blue = becoming more similar to neighbours (homogenization).  Red = becoming more distinct.")
save_dual(p_trend, paste0("fig_homogenization_trend_", RUN_LABEL),
          width = 9, height = 7)

# 7c. 组合图：左趋势 + 右 5 期 facet
combined <- p_trend / p_ts +
  patchwork::plot_layout(heights = c(1.05, 0.85)) +
  patchwork::plot_annotation(
    title = "Spatiotemporal dynamics of biotic homogenization across China",
    subtitle = sprintf("Top: per-grid 5-period trend slope.  Bottom: time-slice maps of local heterogeneity (K=%d nearest neighbours).",
                       K_NEIGHBORS),
    theme = theme_v2_pub(12)
  )
save_dual(combined, paste0("fig_homogenization_combined_", RUN_LABEL),
          width = 16, height = 12.5)

# 7d. PPT vector for the trend map (single map = small enough for dml)
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE)) {
  d <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(value = rvg::dml(ggobj = p_trend),
                      location = officer::ph_location(left = 0.4, top = 0.4,
                                                      width = 9, height = 7))
  print(d, target = file.path(P$figures_v2,
                               paste0("fig_homogenization_trend_",
                                       RUN_LABEL, ".pptx")))
}

message("[stage-12] Done.")
