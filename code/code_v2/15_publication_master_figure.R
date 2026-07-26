#!/usr/bin/env Rscript
## 15_publication_master_figure.R
##
## 顶刊级综合时空动态主图：
##  - 上：每指标 5 期时序曲线 + 95% CRI 带（多 panel）
##  - 中：每指标 5 期空间切片地图（z-score, 跨指标 facet_grid）
##  - 下：年度变化率 bar + 95% CRI（小型 summary）
## 全部用 v2 工具栈：scico + theme_v2_pub + theme_v2_map + 十段线 + 无鹰眼。

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(forcats); library(scales); library(sf)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")

primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
m_long <- read_csv(v2_file("results",
                  paste0("table_community_metrics_with_cri_", RUN_LABEL)),
                  show_col_types = FALSE)
grid_sf <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
china_layers <- load_china_layers(P$china_boundary, P$province_line)

metric_set <- tibble::tribble(
  ~metric,               ~label,                          ~kind,
  "corrected_richness",  "Taxonomic richness",            "Taxonomic",
  "shannon",             "Shannon diversity",             "Taxonomic",
  "pd_prob",             "Faith's PD (prob-weighted)",    "Phylogenetic",
  "mpd_prob",            "Phylogenetic MPD",              "Phylogenetic",
  "trait_volume",        "Functional trait volume",       "Functional",
  "rao_q",               "Functional Rao's Q",            "Functional"
)
metric_set <- metric_set |>
  mutate(label = factor(label, levels = label),
         kind  = factor(kind, levels = c("Taxonomic","Phylogenetic","Functional")))

m_show <- m_long |>
  inner_join(metric_set, by = "metric") |>
  mutate(block_label = factor(block_label, levels = primary_blocks$block_label)) |>
  filter(!is.na(block_label))

kind_pal <- c(Taxonomic = "#0E5A78", Phylogenetic = "#6A4C93",
              Functional = "#C9784A")

## ----- Panel A: trajectory curves with 95% CRI bands ----------------------

traj <- m_show |>
  group_by(label, kind, block_label) |>
  summarise(
    med = median(value_mean, na.rm = TRUE),
    l   = quantile(value_mean, 0.025, na.rm = TRUE),
    u   = quantile(value_mean, 0.975, na.rm = TRUE),
    .groups = "drop")

pA <- ggplot(traj, aes(block_label, med, group = label, colour = kind, fill = kind)) +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = kind_pal, name = NULL) +
  scale_fill_manual(values = kind_pal, guide = "none") +
  labs(
    title    = "a  Detection-corrected diversity trajectories across China (2000-2024)",
    subtitle = "Median across 1,308 grids; shaded band = 2.5-97.5% across-grid quantiles",
    x = NULL, y = NULL
  ) +
  theme_v2_pub(11) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text = element_text(face = "bold"))

## ----- Panel B: spatial time slices (z-scored per metric) -----------------

m_z <- m_show |>
  group_by(label) |>
  mutate(value_z = pmin(pmax(as.numeric(scale(value_mean)), -2.5), 2.5)) |>
  ungroup()
m_sf <- grid_sf |> inner_join(m_z, by = "grid_cell")

pB <- ggplot(m_sf) +
  geom_sf(aes(fill = value_z),
          colour = scales::alpha("white", 0.08), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.10) +
  scale_fill_v2_diverging(
    name = "Within-metric z-score (+/-2.5 clipped)") +
  facet_grid(label ~ block_label) +
  v2_china_coord() +
  theme_v2_map(8.8) +
  labs(title = "b  Spatial time-slice maps: detection-corrected diversity across 5 primary periods",
       subtitle = "Six diversity dimensions across taxonomic, phylogenetic and functional layers; each row z-scaled to expose patterns regardless of native scale")

## ----- Panel C: change-rate bar (2020-24 vs 2000-04) ----------------------

change <- traj |>
  group_by(label, kind) |>
  summarise(start_med = med[block_label == "2000-2004"],
            end_med   = med[block_label == "2020-2024"],
            .groups = "drop") |>
  mutate(pct_change = 100 * (end_med - start_med) / start_med)

# CRI 来自跨网格 quantile 在 start/end 的差（简化）
change_ci <- m_show |>
  filter(block_label %in% c("2000-2004", "2020-2024")) |>
  group_by(label, kind, grid_cell) |>
  summarise(start = value_mean[block_label == "2000-2004"],
            end   = value_mean[block_label == "2020-2024"],
            .groups = "drop") |>
  mutate(delta_pct = 100 * (end - start) / (start + 1e-6)) |>
  group_by(label, kind) |>
  summarise(med = median(delta_pct, na.rm = TRUE),
            l   = quantile(delta_pct, 0.05, na.rm = TRUE),
            u   = quantile(delta_pct, 0.95, na.rm = TRUE),
            .groups = "drop")

pC <- ggplot(change_ci,
             aes(label, med, colour = kind, fill = kind)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.18, linewidth = 0.5,
                colour = "grey25") +
  geom_text(aes(label = sprintf("%+.1f%%", med),
                y = ifelse(med >= 0, u + 1, l - 1.5)),
            colour = "grey15", size = 3.2, fontface = "bold") +
  scale_fill_manual(values = kind_pal, guide = "none") +
  scale_colour_manual(values = kind_pal, guide = "none") +
  labs(title = "c  Per-grid percent change 2020-2024 vs 2000-2004",
       subtitle = "Bar = grid-wise median; whisker = 5-95% across grids",
       x = NULL, y = "% change") +
  theme_v2_pub(11) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

## ----- Assemble + save ----------------------------------------------------

master <- (pA / pB / pC) +
  patchwork::plot_layout(heights = c(1.0, 1.8, 0.85)) +
  patchwork::plot_annotation(
    title = "Detection-corrected community dynamics of Chinese birds (2000-2024)",
    subtitle = "Bayesian multi-species dynamic occupancy on merged & deduplicated citizen-science data | run = v2_full_200sp_ar1",
    caption  = "Source: China Bird Watching Records Platform + eBird/GBIF China (deduplicated). Model: spOccupancy::tMsPGOcc, 4 chains, R-hat <=1.09. Map base: data/中国shp/ with ten-dash line.",
    theme = theme_v2_pub(12)
  )
save_dual(master, paste0("fig_master_diversity_dynamics_", RUN_LABEL),
          width = 16, height = 18)

## 单独保存子 panel 方便 PPT/Word 插入
save_dual(pA, paste0("fig_master_panel_A_trajectories_", RUN_LABEL),
          width = 13, height = 6.5)
save_dual(pC, paste0("fig_master_panel_C_pctchange_", RUN_LABEL),
          width = 11, height = 5.4)

# 编辑级 PPTX（曲线 + bar，地图保存为 PNG 嵌入版）
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg",     quietly = TRUE)) {
  doc <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(value = rvg::dml(ggobj = pA),
                      location = officer::ph_location(left = 0.3, top = 0.3,
                                                      width = 13, height = 6.5)) |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(value = rvg::dml(ggobj = pC),
                      location = officer::ph_location(left = 0.3, top = 0.3,
                                                      width = 11, height = 5.4))
  print(doc, target = file.path(P$figures_v2,
                paste0("fig_master_panels_AC_", RUN_LABEL, ".pptx")))
}

message(sprintf("[stage-15] master publication figure: figures_v2/fig_master_diversity_dynamics_%s.{png,pdf}",
                RUN_LABEL))
