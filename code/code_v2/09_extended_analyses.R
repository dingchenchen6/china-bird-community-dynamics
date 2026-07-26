#!/usr/bin/env Rscript
## 09_extended_analyses.R
##
## 基于 v2 stage-4/5 输出做的 6 类扩展分析（参考 He 2025 GCB / Liang 2024 ZR /
## Sun 2022 SciAdv / Lu 2026 NatCities）。
##
## 输出（results_v2/ + figures_v2/）：
##   A. Spatial pairwise beta + Baselga 分解（每期）+ biotic homogenization 时序
##   B. 物种 gain/loss 分解（per grid × period pair）→ 地图
##   C. Variance partitioning：占域趋势 ~ 气候 / 地形+栖息地 / HFI / 空间（vegan::varpart）
##   D. Urbanization stratification：HFI 四分位上的趋势对比
##   E. 群落类型聚类（k-means on multivariate biodiversity profile）+ 地图
##   F. 纬度梯度强度逐期 quantification（gradient steepness over time）

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(forcats); library(ggplot2); library(sf)
  library(vegan); library(broom); library(patchwork)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_diversity.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")
set.seed(20260503)

message(sprintf("[stage-9] extended analyses for %s", RUN_LABEL))

## --- 0. 加载 -----------------------------------------------------------------

psi_obj  <- readRDS(v2_file("derived",
                paste0("psi_samples_thinned_", RUN_LABEL), "rds"))
psi_arr  <- psi_obj$psi_samples_thinned     # [draws × sp × site × period]
n_draws_avail <- dim(psi_arr)[1]
N_DRAWS_USE <- min(40L, n_draws_avail)       # 仅 40 draws 给 spatial pairwise（够算 CRI）
sel <- if (n_draws_avail > N_DRAWS_USE)
  round(seq(1, n_draws_avail, length.out = N_DRAWS_USE)) else seq_len(n_draws_avail)
psi_arr <- psi_arr[sel, , , , drop = FALSE]
n_draws <- dim(psi_arr)[1]
n_sp    <- dim(psi_arr)[2]
n_site  <- dim(psi_arr)[3]
n_per   <- dim(psi_arr)[4]

primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
visit_effort <- readRDS(file.path(P$derived_v2, "visit_effort_2000_2024.rds"))
sites <- sort(unique(visit_effort$grid_cell))
grid_env <- readRDS(file.path(P$derived_v2,
                    "grid_environment_dynamic_occupancy.rds")) |>
  filter(grid_cell %in% sites) |>
  arrange(match(grid_cell, sites))
grid_sf  <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
metrics  <- read_csv(v2_file("results",
                    paste0("table_community_metrics_with_cri_", RUN_LABEL)),
                    show_col_types = FALSE)
trends   <- read_csv(v2_file("results",
                    paste0("table_grid_trends_with_cri_", RUN_LABEL)),
                    show_col_types = FALSE)

china_layers <- load_china_layers(P$china_boundary, P$province_line)

## --- A. Spatial pairwise beta + Baselga + biotic homogenization -------------
##
##  对每期：把 psi 矩阵当作"行=grid, 列=species, 值=detection prob"用，
##  构造 Sørensen-style probability beta；Baselga 拆 turnover/nestedness。
##  再统计每期内 grid 间相似度的 median + IQR（升高 → 同质化）。

message("[A] spatial pairwise beta per period (posterior subsample)")

beta_pairwise_one_draw <- function(P_mat) {
  # P_mat: site × species; rows are sites, columns species
  # Sørensen-prob: 2A / (2A + B + C); A=sum(min), B=sum(p1-p2)+, C=sum(p2-p1)+
  n_s <- nrow(P_mat)
  out <- matrix(NA_real_, n_s, n_s)
  rs  <- rowSums(P_mat)
  for (i in seq_len(n_s - 1)) {
    pi <- P_mat[i, ]
    for (j in (i + 1):n_s) {
      pj <- P_mat[j, ]
      A <- sum(pmin(pi, pj))
      B <- sum(pmax(pi - pj, 0))
      C <- sum(pmax(pj - pi, 0))
      sor <- (B + C) / (2 * A + B + C + 1e-12)
      out[i, j] <- sor; out[j, i] <- sor
    }
  }
  out
}

# 对 grid 数量太大 (1308) 时全 pairwise 太重；按 stratified sample（200 grids）算
n_sample <- min(180L, n_site)
samp_idx <- sort(sample(seq_len(n_site), n_sample))

agg_per_period <- map_dfr(seq_len(n_per), function(t) {
  # 在 N_DRAWS_USE 个 draws 上分别算，再聚合
  vals <- numeric(0)
  iqrs <- numeric(0)
  for (d in seq_len(n_draws)) {
    P_mat <- t(psi_arr[d, , samp_idx, t])     # site × species
    bm <- beta_pairwise_one_draw(P_mat)
    vals <- c(vals, median(bm[upper.tri(bm)], na.rm = TRUE))
    iqrs <- c(iqrs, IQR(bm[upper.tri(bm)], na.rm = TRUE))
  }
  tibble(block_id = t,
         block_label = primary_blocks$block_label[t],
         median_pairwise_sorensen = mean(vals),
         median_pairwise_sorensen_l95 = quantile(vals, 0.025),
         median_pairwise_sorensen_u95 = quantile(vals, 0.975),
         iqr_pairwise_sorensen = mean(iqrs))
})
write_csv(agg_per_period,
          v2_file("results", paste0("table_spatial_homogenization_", RUN_LABEL)))

homog_plot <- ggplot(agg_per_period, aes(x = block_id)) +
  geom_ribbon(aes(ymin = median_pairwise_sorensen_l95,
                   ymax = median_pairwise_sorensen_u95),
              fill = "#0E5A78", alpha = 0.25) +
  geom_line(aes(y = median_pairwise_sorensen), colour = "#0E5A78",
            linewidth = 0.8) +
  geom_point(aes(y = median_pairwise_sorensen), colour = "#0E5A78", size = 2.6) +
  scale_x_continuous(breaks = seq_len(n_per),
                      labels = primary_blocks$block_label) +
  labs(
    title = "Biotic homogenization across China (occupancy-corrected)",
    subtitle = "Median pairwise Sørensen dissimilarity among 100 km grids per 5-year period",
    caption = "Decline = homogenization (grids becoming more similar). 95% CRI shaded.",
    x = NULL, y = "Median pairwise Sørensen (probability-based)"
  ) +
  theme_v2_pub(11)
save_dual(homog_plot,
          paste0("fig_v2_biotic_homogenization_trend_", RUN_LABEL),
          width = 9, height = 5.4)

## --- B. Species gain/loss decomposition --------------------------------------

message("[B] per-grid species gain/loss across period pairs")

gainloss_one <- function(p1, p2, eps = 1e-6) {
  c(gain = sum(pmax(p2 - p1, 0)), loss = sum(pmax(p1 - p2, 0)),
    net  = sum(p2 - p1))
}

# 用 psi posterior mean（per draw 太重）
psi_mean <- apply(psi_arr, c(2, 3, 4), mean)   # sp × site × period
gl_long <- map_dfr(seq_len(n_per - 1), function(k) {
  res <- map_dfr(seq_len(n_site), function(s) {
    p1 <- psi_mean[, s, k]; p2 <- psi_mean[, s, k + 1]
    gl <- gainloss_one(p1, p2)
    tibble(grid_cell = sites[s], pair_id = k,
           pair_label = sprintf("%s -> %s",
                                  primary_blocks$block_label[k],
                                  primary_blocks$block_label[k + 1]),
           gain = unname(gl["gain"]), loss = unname(gl["loss"]),
           net  = unname(gl["net"]))
  })
})
write_csv(gl_long,
          v2_file("results", paste0("table_species_gain_loss_", RUN_LABEL)))

# 平均 gain/loss/net per grid（跨所有相邻期）
gl_avg <- gl_long |>
  group_by(grid_cell) |>
  summarise(mean_gain = mean(gain),
            mean_loss = mean(loss),
            mean_net  = mean(net),
            mean_total_turnover = mean(gain + loss),
            .groups = "drop")
gl_sf <- grid_sf |> inner_join(gl_avg, by = "grid_cell")

# 4-panel：gain / loss / net / total turnover
gl_long_for_map <- gl_sf |>
  pivot_longer(c(mean_gain, mean_loss, mean_net, mean_total_turnover),
               names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
    mean_gain = "Mean species gain (per period)",
    mean_loss = "Mean species loss (per period)",
    mean_net  = "Mean net change",
    mean_total_turnover = "Mean total turnover"
  ))
gl_map <- ggplot(gl_long_for_map) +
  geom_sf(aes(fill = value), colour = scales::alpha("white", 0.10),
          linewidth = 0.04) +
  china_map_layers(china_layers) +
  scale_fill_v2_diverging(name = "Sum of psi change") +
  facet_wrap(~ metric, ncol = 2) +
  v2_china_coord() +
  theme_v2_map(10.5) +
  labs(
    title = "Per-grid species gain / loss / net change averaged over period pairs",
    subtitle = sprintf("Probability-weighted (sum of psi differences) | run=%s", RUN_LABEL),
    caption = "Gain = sum of psi increases between consecutive periods; loss = sum of psi decreases."
  )
save_dual(gl_map, paste0("fig_v2_species_gain_loss_", RUN_LABEL),
          width = 11, height = 9)

## --- C. Variance partitioning：trend ~ climate + topo + human + space --------

message("[C] variance partitioning of richness trend across driver groups")

vp_dat <- trends |>
  filter(metric == "corrected_richness") |>
  inner_join(grid_env, by = "grid_cell") |>
  drop_na()

# 4 个驱动组
groups <- list(
  climate = c("bio4", "bio7", "bio11", "bio13"),
  topo_habitat = c("elev_mean", "elev_sd", "texture_shannon",
                    "habitat_diversity_shannon"),
  human = c("hfi_mean", "landcover_built", "landcover_cropland"),
  space = c("centroid_lon", "centroid_lat")
)
groups <- lapply(groups, intersect, names(vp_dat))
make_X <- function(vars) {
  out <- vp_dat |> select(all_of(vars)) |>
    mutate(across(everything(), as.numeric)) |>
    mutate(across(everything(), ~ as.numeric(scale(.x))))
  out
}

vp <- vegan::varpart(vp_dat$trend_mean,
                      make_X(groups$climate),
                      make_X(groups$topo_habitat),
                      make_X(groups$human),
                      make_X(groups$space))
saveRDS(vp, v2_file("derived",
                     paste0("varpart_richness_trend_", RUN_LABEL), "rds"))

vp_part <- as_tibble(vp$part$indfract, rownames = "fraction")
adjcol <- grep("^[Aa]dj", names(vp_part), value = TRUE)[1]
vp_part <- vp_part |>
  rename(adj_r2 = !!adjcol) |>
  mutate(fraction = factor(fraction, levels = unique(fraction)))
write_csv(vp_part,
          v2_file("results",
                  paste0("table_varpart_richness_trend_", RUN_LABEL)))

vp_pure <- vp_part |>
  filter(grepl("^\\[\\w+\\]$", fraction)) |>
  mutate(component = recode(as.character(fraction),
    "[a]" = "Climate (pure)", "[b]" = "Topo+Habitat (pure)",
    "[c]" = "Human (pure)", "[d]" = "Space (pure)")) |>
  mutate(adj_R2 = pmax(adj_r2, 0))
vp_plot <- ggplot(vp_pure, aes(reorder(component, adj_R2), adj_R2)) +
  geom_col(fill = "#0E5A78", alpha = 0.85, width = 0.55) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * adj_R2)),
            hjust = -0.15, size = 3.6) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                      labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Variance partitioning of grid-level richness trend",
       subtitle = sprintf("Independent (pure) fractions of adj-R2 | run=%s", RUN_LABEL),
       x = NULL, y = "Pure fraction of variance explained") +
  theme_v2_pub(11)
save_dual(vp_plot, paste0("fig_v2_varpart_richness_trend_", RUN_LABEL),
          width = 8, height = 5)

## --- D. Urbanization stratification (HFI quartiles) -------------------------

message("[D] HFI-quartile stratified trend comparison")

q_breaks <- quantile(grid_env$hfi_mean, c(0, .25, .5, .75, 1), na.rm = TRUE)
hfi_class <- tibble(grid_cell = grid_env$grid_cell,
                     hfi_q = cut(grid_env$hfi_mean, breaks = unique(q_breaks),
                                 include.lowest = TRUE,
                                 labels = c("Q1 (low HFI)", "Q2", "Q3",
                                             "Q4 (high HFI)")))

hfi_trends <- trends |>
  filter(metric %in% c("corrected_richness", "shannon",
                        "trait_volume", "pd_prob")) |>
  inner_join(hfi_class, by = "grid_cell") |>
  mutate(metric_label = recode(metric,
    corrected_richness = "Richness trend",
    shannon = "Shannon trend",
    trait_volume = "Trait-volume trend",
    pd_prob = "Faith's PD trend"))

hfi_summary <- hfi_trends |>
  group_by(metric_label, hfi_q) |>
  summarise(median = median(trend_mean, na.rm = TRUE),
            l95 = quantile(trend_mean, 0.025, na.rm = TRUE),
            u95 = quantile(trend_mean, 0.975, na.rm = TRUE),
            n = n(), .groups = "drop")
write_csv(hfi_summary,
          v2_file("results", paste0("table_hfi_stratified_trends_", RUN_LABEL)))

hfi_plot <- ggplot(hfi_trends, aes(hfi_q, trend_mean, fill = hfi_q,
                                    colour = hfi_q)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_violin(alpha = 0.25, colour = NA, width = 0.85, trim = FALSE) +
  ggbeeswarm::geom_quasirandom(alpha = 0.18, size = 0.4, width = 0.20,
                                show.legend = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85,
               colour = "#3A3A3A", linewidth = 0.2) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_fill_brewer(palette = "YlOrRd", guide = "none") +
  scale_colour_brewer(palette = "YlOrRd", guide = "none") +
  labs(title = "Community dynamic trends across human-footprint quartiles",
       subtitle = sprintf("Q1 = least urbanised; Q4 = most urbanised | run=%s",
                          RUN_LABEL),
       x = NULL, y = "Trend (per period)") +
  theme_v2_pub(11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
save_dual(hfi_plot, paste0("fig_v2_hfi_stratified_trends_", RUN_LABEL),
          width = 11, height = 8)

## --- E. Community typology via k-means clustering ---------------------------

message("[E] community typology (k-means on multivariate biodiversity profile)")

prof_wide <- metrics |>
  filter(metric %in% c("corrected_richness", "shannon",
                        "trait_volume", "rao_q", "pd_prob", "mpd_prob")) |>
  group_by(grid_cell, metric) |>
  summarise(mean_value = mean(value_mean, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = metric, values_from = mean_value) |>
  drop_na()
mat <- prof_wide |> select(-grid_cell) |> as.matrix()
mat_z <- scale(mat)

ks <- 2:8
wss <- sapply(ks, function(k)
  kmeans(mat_z, centers = k, nstart = 25, iter.max = 100)$tot.withinss)
elbow <- ks[which.min(diff(wss) / wss[-length(wss)])] + 1
elbow <- min(max(elbow, 3), 6)        # 限制在 3-6
km <- kmeans(mat_z, centers = elbow, nstart = 50, iter.max = 200)
prof_wide$cluster <- factor(km$cluster,
                              levels = seq_len(elbow),
                              labels = paste0("Type-", LETTERS[seq_len(elbow)]))

cluster_centers <- as_tibble(km$centers, rownames = "cluster_id") |>
  mutate(cluster_id = factor(cluster_id,
                              labels = paste0("Type-", LETTERS[seq_len(elbow)])))
write_csv(cluster_centers,
          v2_file("results", paste0("table_community_types_centers_", RUN_LABEL)))
write_csv(prof_wide |> select(grid_cell, cluster),
          v2_file("results", paste0("table_community_types_grids_", RUN_LABEL)))

cluster_sf <- grid_sf |> inner_join(prof_wide |> select(grid_cell, cluster),
                                      by = "grid_cell")
type_map <- ggplot(cluster_sf) +
  geom_sf(aes(fill = cluster), colour = scales::alpha("white", 0.10),
          linewidth = 0.04) +
  china_map_layers(china_layers) +
  scale_fill_manual(values = V2_PALETTES$qualitative[seq_len(elbow)],
                     name = "Community type") +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = sprintf("China bird community typology — %d-type k-means",
                        elbow),
       subtitle = sprintf("Clustering on z-scaled occupancy-corrected diversity profile | run=%s",
                          RUN_LABEL),
       caption = "Each grid assigned to nearest centroid in standardised diversity space.")

# 类型 z-profile 雷达式 bar
profile_long <- as_tibble(km$centers, rownames = "type") |>
  mutate(type = paste0("Type-", LETTERS[as.integer(type)])) |>
  pivot_longer(-type, names_to = "metric", values_to = "z") |>
  mutate(metric = recode(metric,
    corrected_richness = "Richness",
    shannon = "Shannon",
    trait_volume = "Trait vol",
    rao_q = "Rao Q",
    pd_prob = "PD",
    mpd_prob = "MPD"))
profile_plot <- ggplot(profile_long,
                        aes(metric, z, fill = z > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey55") +
  facet_wrap(~ type, ncol = elbow) +
  scale_fill_manual(values = c(`TRUE` = "#8B2E1E", `FALSE` = "#0E5A78"),
                     guide = "none") +
  labs(title = "Z-profile of each community type",
       subtitle = "Positive (red) = higher than national mean; Negative (blue) = lower",
       x = NULL, y = "z (sd units)") +
  theme_v2_pub(10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
combined <- type_map / profile_plot + patchwork::plot_layout(heights = c(2.6, 1))
save_dual(combined, paste0("fig_v2_community_types_", RUN_LABEL),
          width = 12, height = 11)

## --- F. Latitudinal gradient strength per period ----------------------------

message("[F] latitudinal gradient steepness per period")

lat_dat <- metrics |>
  filter(metric %in% c("corrected_richness", "shannon", "trait_volume",
                        "pd_prob", "rao_q")) |>
  inner_join(grid_env |> select(grid_cell, centroid_lat), by = "grid_cell")

lat_fit <- lat_dat |>
  group_by(metric, block_label) |>
  summarise(
    fit = list(lm(value_mean ~ centroid_lat, data = pick(everything()))),
    .groups = "drop") |>
  mutate(slope = map_dbl(fit, ~ coef(.x)[["centroid_lat"]]),
         se    = map_dbl(fit, ~ summary(.x)$coefficients["centroid_lat", "Std. Error"]),
         r2    = map_dbl(fit, ~ summary(.x)$r.squared)) |>
  select(-fit) |>
  mutate(block_label = factor(block_label, levels = primary_blocks$block_label),
         metric = factor(metric))
write_csv(lat_fit,
          v2_file("results", paste0("table_latitudinal_gradient_strength_", RUN_LABEL)))

lat_plot <- ggplot(lat_fit, aes(block_label, slope, group = metric,
                                  colour = metric)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se),
                width = 0.15, linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = V2_PALETTES$qualitative[1:5], guide = "none") +
  labs(title = "Latitudinal gradient steepness across primary periods",
       subtitle = sprintf("Slope of metric ~ centroid_lat per 5-year period (95%% CI) | run=%s",
                          RUN_LABEL),
       caption = "Negative = metric declines with latitude (typical N->S pattern).",
       x = NULL, y = "Slope (per degree latitude)") +
  theme_v2_pub(11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
save_dual(lat_plot, paste0("fig_v2_latitudinal_gradient_", RUN_LABEL),
          width = 12, height = 6)

message("[stage-9] Extended analyses done.")
