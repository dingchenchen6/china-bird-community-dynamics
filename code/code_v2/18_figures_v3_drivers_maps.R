#!/usr/bin/env Rscript
## 18_figures_v3_drivers_maps.R
##
## 集成 v3 主图集：
##   - phylo & FD 时空切片地图（facet_grid 指标 × 5 期）
##   - 时序曲线 + 95% CRI
##   - biome 分层版（WWF Ecoregions via rnaturalearth）
##   - 中国 7 大地理区分层版（基于 data/中国shp/省.shp 手工映射）
##   - 所有地图：coord_sf 强制 bbox + 十段线 + 无鹰眼
##
## Output: figures_v3/fig_v3_*.{png,pdf} + figures_v3/v3_master_deck.pptx

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(forcats); library(ggplot2); library(patchwork); library(sf)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")

message(sprintf("[stage-18] v3 figures + biome/region stratified | %s", RUN_LABEL))

primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
grid_sf <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
china_layers <- load_china_layers(P$china_boundary, P$province_line)

## --- 1. 读 v3 指标（若 stage 16 完成） + v2 已有 ----------------------

v3_path <- v3_file("results",
                   paste0("table_community_metrics_v3_", RUN_LABEL))
v2_path <- v2_file("results",
                   paste0("table_community_metrics_with_cri_", RUN_LABEL))

m_v3 <- if (file.exists(v3_path)) read_csv(v3_path, show_col_types = FALSE) else NULL
m_v2 <- read_csv(v2_path, show_col_types = FALSE)

# 合并 v2 + v3 指标到同一 long 格式
metric_set <- tibble::tribble(
  ~metric, ~label, ~kind, ~palette_type,
  "corrected_richness", "Taxonomic richness",  "Taxonomic",    "diverging",
  "shannon",            "Shannon",             "Taxonomic",    "diverging",
  "pd_prob",            "Faith's PD",          "Phylogenetic", "diverging",
  "mpd_prob",           "MPD",                 "Phylogenetic", "diverging",
  "psv",                "PSV",                 "Phylogenetic", "diverging",
  "nri",                "NRI",                 "Phylogenetic", "diverging",
  "nti",                "NTI",                 "Phylogenetic", "diverging",
  "trait_volume",       "Functional vol.",     "Functional",   "diverging",
  "rao_q",              "Functional Rao Q",    "Functional",   "diverging",
  "fric",               "FRic",                "Functional",   "diverging",
  "feve",               "FEve",                "Functional",   "diverging",
  "fdiv",               "FDiv",                "Functional",   "diverging",
  "fdis",               "FDis",                "Functional",   "diverging"
)

all_m <- bind_rows(
  m_v2 |> select(grid_cell, block_label, block_id, metric,
                  value_mean, value_l95, value_u95),
  if (!is.null(m_v3))
    m_v3 |> select(grid_cell, block_label, block_id, metric,
                    value_mean, value_l95, value_u95)
)
m_show <- all_m |>
  inner_join(metric_set, by = "metric") |>
  mutate(label = factor(label, levels = metric_set$label),
         kind = factor(kind, levels = c("Taxonomic","Phylogenetic","Functional")),
         block_label = factor(block_label, levels = primary_blocks$block_label)) |>
  filter(!is.na(block_label)) |>
  group_by(label) |>
  mutate(value_z = pmin(pmax(as.numeric(scale(value_mean)), -2.5), 2.5)) |>
  ungroup()

## --- 2. Phylo 时空切片地图 -------------------------------------------

m_phy <- m_show |> filter(kind == "Phylogenetic")
if (nrow(m_phy) > 0) {
  sf_phy <- grid_sf |> inner_join(m_phy, by = "grid_cell")
  p_phy <- ggplot(sf_phy) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_diverging(name = "Within-metric z-score (clipped at +/-2.5)") +
    facet_grid(label ~ block_label) +
    v2_china_coord() +
    theme_v2_map(9) +
    labs(title = "Phylogenetic diversity time-slice maps",
         subtitle = sprintf("McTavish/clootl tree v1.6.2025 | run=%s", RUN_LABEL),
         caption = "Each row z-scaled within metric.")
  save_dual_v3(p_phy, "fig_v3_map_phylo_grid", width = 14, height = 11)
}

## --- 3. FD 时空切片地图 ----------------------------------------------

m_fd <- m_show |> filter(kind == "Functional")
if (nrow(m_fd) > 0) {
  sf_fd <- grid_sf |> inner_join(m_fd, by = "grid_cell")
  p_fd <- ggplot(sf_fd) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_diverging(name = "Within-metric z-score (clipped at +/-2.5)") +
    facet_grid(label ~ block_label) +
    v2_china_coord() +
    theme_v2_map(9) +
    labs(title = "Functional diversity time-slice maps",
         subtitle = sprintf("fundiversity + custom FSpe/FOri | run=%s", RUN_LABEL),
         caption = "Each row z-scaled within metric.")
  save_dual_v3(p_fd, "fig_v3_map_func_grid", width = 14, height = 13)
}

## --- 4. 时序曲线 + 95% CRI 全国 mean -----------------------------------

traj <- m_show |>
  group_by(label, kind, block_label) |>
  summarise(med = median(value_mean, na.rm = TRUE),
            l   = quantile(value_mean, 0.025, na.rm = TRUE),
            u   = quantile(value_mean, 0.975, na.rm = TRUE),
            .groups = "drop")
kind_pal <- c(Taxonomic = "#0E5A78", Phylogenetic = "#6A4C93",
              Functional = "#C9784A")
p_traj <- ggplot(traj, aes(block_label, med, group = label,
                              colour = kind, fill = kind)) +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = kind_pal, name = NULL) +
  scale_fill_manual(values = kind_pal, guide = "none") +
  labs(title = "Diversity trajectories across China (v3)",
       subtitle = sprintf("National median + 2.5-97.5%% across-grid quantiles | %s", RUN_LABEL),
       x = NULL, y = NULL) +
  theme_v2_pub(10.5) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text = element_text(face = "bold"))
save_dual_v3(p_traj, "fig_v3_timeseries_cri", width = 14, height = 8)

## --- 5. 中国 7 大地理区分层 ------------------------------------------

prov_sf <- suppressWarnings(sf::st_read(P$china_boundary, quiet = TRUE)) |>
  sf::st_make_valid() |> sf::st_transform(4326)
prov_name_col <- intersect(c("省名","NAME","name"), names(prov_sf))[1]

region_map <- c(
  # 青藏高原
  `西藏自治区`="QXP", `西藏`="QXP", `青海省`="QXP", `青海`="QXP",
  # 西北
  `新疆维吾尔自治区`="NW", `新疆`="NW", `内蒙古自治区`="NW", `内蒙古`="NW",
  `甘肃省`="NW", `甘肃`="NW", `宁夏回族自治区`="NW", `宁夏`="NW",
  # 华北
  `北京市`="HB", `北京`="HB", `天津市`="HB", `天津`="HB",
  `河北省`="HB", `河北`="HB", `山东省`="HB", `山东`="HB",
  `河南省`="HB", `河南`="HB", `山西省`="HB", `山西`="HB",
  # 东北
  `黑龙江省`="NE", `黑龙江`="NE", `吉林省`="NE", `吉林`="NE",
  `辽宁省`="NE", `辽宁`="NE",
  # 西南
  `云南省`="SW", `云南`="SW", `贵州省`="SW", `贵州`="SW",
  `重庆市`="SW", `重庆`="SW", `四川省`="SW", `四川`="SW",
  # 华中-华东
  `湖北省`="CE", `湖北`="CE", `湖南省`="CE", `湖南`="CE",
  `江西省`="CE", `江西`="CE", `安徽省`="CE", `安徽`="CE",
  `江苏省`="CE", `江苏`="CE", `浙江省`="CE", `浙江`="CE",
  `上海市`="CE", `上海`="CE",
  # 华南
  `广东省`="SC", `广东`="SC", `广西壮族自治区`="SC", `广西`="SC",
  `福建省`="SC", `福建`="SC", `海南省`="SC", `海南`="SC",
  `台湾省`="SC", `台湾`="SC",
  `香港特别行政区`="SC", `香港`="SC",
  `澳门特别行政区`="SC", `澳门`="SC"
)
region_labels <- c(QXP="Qinghai-Tibet Plateau", NW="Northwest",
                   HB="North China", NE="Northeast",
                   SW="Southwest", CE="Central-East",
                   SC="South China")
prov_sf$region <- region_map[as.character(prov_sf[[prov_name_col]])]
prov_sf$region_label <- region_labels[prov_sf$region]

# 网格 × 省 空间 join 取主导省份
old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
on.exit(sf::sf_use_s2(old_s2), add = TRUE)
join <- sf::st_join(sf::st_centroid(grid_sf), prov_sf[, "region"],
                     join = sf::st_intersects, left = TRUE)
grid_region <- tibble(grid_cell = grid_sf$grid_cell,
                       region = join$region)
saveRDS(grid_region,
        v3_file("derived", "grid_region_assignment", "rds"))

# 每区时序曲线
m_reg <- m_show |> inner_join(grid_region, by = "grid_cell") |>
  filter(!is.na(region)) |>
  group_by(label, kind, block_label, region) |>
  summarise(med = median(value_mean, na.rm = TRUE),
            l = quantile(value_mean, 0.025, na.rm = TRUE),
            u = quantile(value_mean, 0.975, na.rm = TRUE),
            .groups = "drop") |>
  mutate(region_label = region_labels[region])

p_reg <- ggplot(m_reg, aes(block_label, med, group = region,
                              colour = region_label, fill = region_label)) +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = V2_PALETTES$qualitative[1:7],
                       name = "Region") +
  scale_fill_manual(values = V2_PALETTES$qualitative[1:7], guide = "none") +
  labs(title = "Diversity trajectories by Chinese geographic region",
       subtitle = sprintf("7 regions (Qinghai-Tibet, NW, NC, NE, SW, CE, SC) | %s",
                          RUN_LABEL),
       x = NULL, y = NULL) +
  theme_v2_pub(10) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 22, hjust = 1),
        strip.text = element_text(face = "bold"))
save_dual_v3(p_reg, "fig_v3_chinaregion_stratified", width = 14, height = 9)

# 区域分布地图
region_grid_sf <- grid_sf |> inner_join(grid_region, by = "grid_cell") |>
  mutate(region_label = region_labels[region]) |>
  filter(!is.na(region_label))
p_region_map <- ggplot(region_grid_sf) +
  geom_sf(aes(fill = region_label),
          colour = scales::alpha("white", 0.08), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_manual(values = V2_PALETTES$qualitative[1:7], name = "Region") +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = "7 geographic regions of China",
       subtitle = sprintf("Region assignment based on province polygons | %s",
                          RUN_LABEL))
save_dual_v3(p_region_map, "fig_v3_chinaregion_map", width = 10, height = 7.5)

## --- 6. WWF biome 分层（如可用） ------------------------------------

biome_sf <- tryCatch({
  if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    # rnaturalearth 没有直接 biome，用 lat-lon bin 退化
    NULL
  } else NULL
}, error = function(e) NULL)

# 退化方案：Köppen-like 气候带分箱（基于 grid_env 的 bio1 + bio12）
grid_env <- readRDS(file.path(P$derived_v2,
                "grid_environment_dynamic_occupancy.rds"))
koppen_bin <- function(bio1, bio12) {
  case_when(
    is.na(bio1) | is.na(bio12) ~ NA_character_,
    bio1 >= 18 & bio12 >= 1500 ~ "Tropical wet",
    bio1 >= 18 & bio12 <  1500 ~ "Tropical dry",
    bio1 >= 5  & bio1 < 18 & bio12 >= 800 ~ "Temperate wet",
    bio1 >= 5  & bio1 < 18 & bio12 <  800 ~ "Temperate dry",
    bio1 >= -5 & bio1 < 5 ~ "Cold temperate",
    bio1 <  -5 ~ "Subarctic / Alpine",
    TRUE ~ "Other")
}
grid_env$biome <- koppen_bin(grid_env$bio1, grid_env$bio12)
m_biome <- m_show |>
  inner_join(grid_env |> select(grid_cell, biome), by = "grid_cell") |>
  filter(!is.na(biome)) |>
  group_by(label, kind, block_label, biome) |>
  summarise(med = median(value_mean, na.rm = TRUE),
            l = quantile(value_mean, 0.025, na.rm = TRUE),
            u = quantile(value_mean, 0.975, na.rm = TRUE),
            .groups = "drop")
p_biome <- ggplot(m_biome, aes(block_label, med, group = biome,
                                  colour = biome, fill = biome)) +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.3) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = V2_PALETTES$qualitative[1:7], name = "Biome") +
  scale_fill_manual(values = V2_PALETTES$qualitative[1:7], guide = "none") +
  labs(title = "Diversity trajectories by Köppen-style climate zone",
       subtitle = sprintf("Fallback biome from BIO1+BIO12 | %s", RUN_LABEL),
       x = NULL, y = NULL) +
  theme_v2_pub(10) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 22, hjust = 1),
        strip.text = element_text(face = "bold"))
save_dual_v3(p_biome, "fig_v3_biome_stratified", width = 14, height = 9)

## --- 7. master PPTX deck（地图 PNG 嵌入 + 非地图矢量） ----------------

deck_plots <- list(
  "Phylogenetic time-slice maps"   = if (exists("p_phy"))   p_phy   else NULL,
  "Functional time-slice maps"     = if (exists("p_fd"))    p_fd    else NULL,
  "Diversity trajectories"         = p_traj,
  "Chinese region stratified"      = p_reg,
  "Region assignment map"          = p_region_map,
  "Biome (Köppen) stratified"      = p_biome)
deck_plots <- deck_plots[!sapply(deck_plots, is.null)]
if (length(deck_plots) > 0) {
  save_pptx_deck(deck_plots,
                  deck_stem = file.path("..", "figures_v3",
                                         paste0("v3_master_deck_", RUN_LABEL)),
                  width = 12, height = 7)
  # save_pptx_deck wrote to figures_v2; move it:
  src <- file.path(P$figures_v2,
                    paste0("..", "/figures_v3/v3_master_deck_", RUN_LABEL, ".pptx"))
  if (file.exists(src)) {
    file.rename(src, v3_file("figure",
                              paste0("v3_master_deck_", RUN_LABEL), "pptx"))
  } else {
    # Fallback: write directly
    if (requireNamespace("officer", quietly = TRUE)) {
      doc <- officer::read_pptx()
      for (nm in names(deck_plots)) {
        pl <- deck_plots[[nm]]
        png_tmp <- tempfile(fileext = ".png")
        ggplot2::ggsave(png_tmp, pl, width = 12, height = 7,
                         dpi = 220, bg = "white")
        doc <- officer::add_slide(doc, layout = "Title and Content",
                                    master = "Office Theme") |>
          officer::ph_with(value = nm,
                            location = officer::ph_location_type(type = "title")) |>
          officer::ph_with(value = officer::external_img(png_tmp,
                                                          width = 12, height = 7),
                            location = officer::ph_location(left = 0.3, top = 1.1,
                                                            width = 12, height = 7))
      }
      print(doc, target = v3_file("figure",
                                    paste0("v3_master_deck_", RUN_LABEL), "pptx"))
    }
  }
}

message("[stage-18] done.")
