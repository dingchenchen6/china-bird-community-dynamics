#!/usr/bin/env Rscript
## 03_prepare_environment.R
##
## 阶段 3：把 v1 已抽取好的网格 × 环境驱动表对齐到 v2 网格。
## v2 与 v1 grid_cell 编号 100% 一致（同一边界、同一分辨率、同一算法），
## 因此本阶段不重新跑栅格抽取，仅做：
##   1) 对齐 + 校验（grid_cell × area_km2 × centroid 对得上）
##   2) 落盘到 data/derived_v2/grid_environment_dynamic_occupancy.rds
##   3) 输出层覆盖审计 manifest

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))

P <- ensure_v2_dirs()

message("[stage-3] Aligning v1 environmental layer table to v2 grid")

v1_env_path <- file.path(P$derived_v1, "grid_environment_dynamic_occupancy.rds")
v2_grid     <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))

stopifnot(file.exists(v1_env_path))
v1_env <- readRDS(v1_env_path)

# 校验：grid_cell 集合必须一致
miss_v1 <- setdiff(v2_grid$grid_cell, v1_env$grid_cell)
miss_v2 <- setdiff(v1_env$grid_cell, v2_grid$grid_cell)
stopifnot(length(miss_v1) == 0L, length(miss_v2) == 0L)

# 对齐顺序、补 area_km2 / centroid（如不一致以 v2 网格为准）
v2_meta <- sf::st_drop_geometry(v2_grid) |>
  as_tibble() |>
  select(grid_cell, grid_id, centroid_lon, centroid_lat, area_km2)

env_v2 <- v2_meta |>
  left_join(v1_env |> select(-any_of(c("grid_id", "centroid_lon",
                                        "centroid_lat", "area_km2"))),
            by = "grid_cell")

saveRDS(env_v2, file.path(P$derived_v2, "grid_environment_dynamic_occupancy.rds"))
write_csv(env_v2, v2_file("results", "table_grid_environment_dynamic_occupancy"))

# layer manifest
layer_groups <- list(
  climate_bioclim   = paste0("bio", 1:19),
  topography        = c("elev_mean", "elev_sd"),
  texture           = c("texture_shannon", "texture_entropy", "texture_contrast"),
  productivity      = c("npp_mean", "ndvi_mean"),
  landcover         = c("landcover_trees", "landcover_cropland", "landcover_built",
                         "landcover_shrubs", "landcover_grassland", "landcover_water"),
  habitat_diversity = c("habitat_diversity_shannon", "habitat_diversity_richness",
                         "habitat_diversity_evenness", "habitat_dominance",
                         "habitat_diversity_simpson",
                         "natural_landcover_fraction", "human_modified_fraction"),
  human_footprint   = c("hfi_2000","hfi_2005","hfi_2010","hfi_2015","hfi_2020",
                         "hfi_2024","hfi_mean","hfi_sd"),
  geography         = c("centroid_lon", "centroid_lat", "area_km2")
)
manifest <- bind_rows(lapply(names(layer_groups), function(g) {
  vars <- layer_groups[[g]]
  tibble(
    layer_group = g,
    variable    = vars,
    available   = vars %in% names(env_v2),
    n_finite    = vapply(vars, function(v)
      if (v %in% names(env_v2)) sum(is.finite(env_v2[[v]])) else 0L, integer(1)),
    coverage_pct = vapply(vars, function(v)
      if (v %in% names(env_v2))
        100 * mean(is.finite(env_v2[[v]])) else 0, numeric(1))
  )
}))
write_csv(manifest, v2_file("results", "table_dynamic_environment_layer_manifest"))

message(sprintf("  v2 grid env: %d cells × %d variables",
                nrow(env_v2), ncol(env_v2)))
message(sprintf("  available variables: %d / %d",
                sum(manifest$available), nrow(manifest)))
print(manifest |> filter(!available))
message("[stage-3] Environment table aligned.")
