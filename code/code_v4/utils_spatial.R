#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把 v3 的 extract_env_from_raster 从 03_prepare_environment.R 抽离，
#   消除 "调用先于定义" 的潜在 source-order bug（C5 修复）
#
# Objective / 分析目标:
#   提供 project_china_albers / extract_env_from_raster / centroids_lonlat
#   等空间工具，所有 03 / 03c-e 共用
#
# Input data / 输入数据:
#   00_config.R + utils_paths.R 已 source
#
# Main workflow / 主要流程:
#   1. project_china_albers：sf 对象转 China Albers 等积投影
#   2. extract_env_from_raster：从外部栅格按网格质心提取
#   3. nearest_neighbor_impute：缺值用空间最近邻填补
#
# Key assumptions / 关键假设:
#   sf / terra 可用；外部栅格按 data/external/ 目录组织
#
# Main packages / 主要包:
#   sf, terra
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

suppressPackageStartupMessages({
  library(sf)
})

# ── project_china_albers ─────────────────────────────────────────────
#' 投影到中国 Albers 等积坐标系
#' 中央经线 105°E，标准纬线 25°N/47°N，原点 (0, 0)
project_china_albers <- function(sf_obj) {
  albers_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  st_transform(sf_obj, albers_crs)
}

# ── centroids_lonlat ─────────────────────────────────────────────────
#' 计算 sf 多边形质心并返回 lon/lat 表
centroids_lonlat <- function(grid_sf) {
  ctr <- suppressWarnings(st_centroid(grid_sf))
  coords <- st_coordinates(st_transform(ctr, 4326))
  data.frame(grid_cell = grid_sf$grid_cell,
             centroid_lon = coords[, 1],
             centroid_lat = coords[, 2],
             stringsAsFactors = FALSE)
}

# ── extract_env_from_raster（v4：抽离自 v3 03 脚本） ────────────────
#' 从外部栅格按网格质心提取环境协变量
#' @param grid_sf sf 多边形，含 grid_cell 列
#' @param sites integer：要保留的 grid_cell 子集
#' @param vars character：环境变量名（bio4/7/11/13、elev_mean 等）
#' @param raster_dirs named list：每类栅格的目录候选（按优先级）
#' @return tibble: grid_cell + 各环境变量列
extract_env_from_raster <- function(grid_sf, sites,
                                     vars = c("bio4", "bio7", "bio11", "bio13",
                                              "elev_mean", "elev_sd",
                                              "hfi_mean",
                                              "landcover_built", "landcover_cropland",
                                              "texture_shannon",
                                              "habitat_diversity_shannon"),
                                     raster_dirs = list(
                                       worldclim = c(WORLDCLIM_DIR,
                                                      file.path(DIRS$external, "worldclim")),
                                       elevation = c(file.path(DIRS$external, "elevation"),
                                                     file.path(DIRS$data_raw, "elevation")),
                                       hfi       = c(HFI_DIR),
                                       landcover = c(CLCD_DIR),
                                       texture   = c(file.path(DIRS$external, "earthenv_texture")),
                                       habitat   = c(file.path(DIRS$external, "habitat_diversity"))
                                     ),
                                     v3_fallback = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    warning("terra package required for raster extraction")
    return(NULL)
  }

  if (is.null(grid_sf) || nrow(grid_sf) == 0) {
    warning("[extract_env] grid_sf empty")
    return(NULL)
  }

  grid_subset <- grid_sf[grid_sf$grid_cell %in% sites, ]
  if (nrow(grid_subset) == 0) {
    warning("[extract_env] no grids match sites")
    return(NULL)
  }

  cent <- centroids_lonlat(grid_subset)
  pts <- st_as_sf(cent, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
  result <- tibble::tibble(grid_cell    = cent$grid_cell,
                            centroid_lon = cent$centroid_lon,
                            centroid_lat = cent$centroid_lat)

  # 占位：所有 vars 默认 NA，后续逐组填
  for (v in vars) result[[v]] <- NA_real_

  # ── WorldClim BIO ──────────────────────────────────────────────────
  bio_vars <- intersect(vars, c("bio4", "bio7", "bio11", "bio13",
                                "bio1", "bio12"))
  if (length(bio_vars) > 0) {
    bio_stack <- .find_and_stack(raster_dirs$worldclim,
                                  pattern = "^bio[0-9]+.*\\.tif$",
                                  min_layers = 4)
    if (!is.null(bio_stack)) {
      bio_vals <- terra::extract(bio_stack, pts, ID = FALSE)
      for (bv in bio_vars) {
        if (bv %in% names(bio_vals)) result[[bv]] <- bio_vals[[bv]]
      }
      message(sprintf("[extract_env] WorldClim BIO extracted: %d vars",
                      sum(bio_vars %in% names(bio_vals))))
    }
  }

  # ── Elevation（mean + sd） ─────────────────────────────────────────
  if (any(c("elev_mean", "elev_sd") %in% vars)) {
    elev_rast <- .find_first_raster(raster_dirs$elevation,
                                     candidate_names = c("dem_china.tif",
                                                          "srtm_china.tif",
                                                          "elev_china.tif"))
    if (!is.null(elev_rast)) {
      ev <- terra::extract(elev_rast, pts, ID = FALSE)
      if ("elev_mean" %in% vars) result$elev_mean <- ev[, 1]
      message("[extract_env] elevation extracted")
    }
  }

  # ── HFI ───────────────────────────────────────────────────────────
  if ("hfi_mean" %in% vars) {
    hfi_rast <- .find_first_raster(raster_dirs$hfi,
                                    pattern = "^hfp.*\\.tif$")
    if (!is.null(hfi_rast)) {
      hv <- terra::extract(hfi_rast, pts, ID = FALSE)
      result$hfi_mean <- hv[, 1]
      message("[extract_env] HFI extracted")
    }
  }

  # ── v3 fallback：缺失列用 v3 grid_env 最近邻填补 ───────────────────
  if (v3_fallback) {
    v3_env <- safe_read(file.path(DIRS$v3_derived, "grid_environment_v3.rds"),
                         quiet = TRUE)
    if (is.null(v3_env)) {
      v3_env <- safe_read(file.path(DIRS$v2_derived,
                                       "grid_environment_dynamic_occupancy.rds"),
                           quiet = TRUE)
    }
    if (!is.null(v3_env)) {
      result <- .nn_impute_from_grid_env(result, v3_env, vars)
    }
  }

  tibble::as_tibble(result)
}

# ── 内部辅助：找到第一个匹配的栅格 ────────────────────────────────
.find_first_raster <- function(dirs, candidate_names = NULL, pattern = NULL) {
  for (d in dirs) {
    if (!dir.exists(d)) next
    if (!is.null(candidate_names)) {
      for (nm in candidate_names) {
        p <- file.path(d, nm)
        if (file.exists(p)) return(tryCatch(terra::rast(p), error = function(e) NULL))
      }
    }
    if (!is.null(pattern)) {
      files <- list.files(d, pattern = pattern, full.names = TRUE)
      if (length(files) > 0) {
        return(tryCatch(terra::rast(files[1]), error = function(e) NULL))
      }
    }
  }
  NULL
}

.find_and_stack <- function(dirs, pattern, min_layers = 1) {
  for (d in dirs) {
    if (!dir.exists(d)) next
    files <- list.files(d, pattern = pattern, full.names = TRUE)
    if (length(files) >= min_layers) {
      return(tryCatch(terra::rast(files), error = function(e) NULL))
    }
  }
  NULL
}

# ── 空间最近邻填补 ────────────────────────────────────────────────
.nn_impute_from_grid_env <- function(result, ref_env, vars) {
  if (!all(c("centroid_lon", "centroid_lat") %in% names(ref_env))) return(result)

  na_cols <- vars[sapply(vars, function(v) {
    v %in% names(result) && sum(is.na(result[[v]])) > nrow(result) * 0.5
  })]
  if (length(na_cols) == 0) return(result)

  message(sprintf("[extract_env] NN-imputing %d high-NA columns from v3/v2 grid_env",
                  length(na_cols)))

  for (nc in na_cols) {
    if (!nc %in% names(ref_env)) next
    for (i in seq_len(nrow(result))) {
      if (!is.na(result[[nc]][i])) next
      dlat <- abs(ref_env$centroid_lat - result$centroid_lat[i])
      dlon <- abs(ref_env$centroid_lon - result$centroid_lon[i])
      nearest <- which.min(dlat^2 + dlon^2)
      if (length(nearest) > 0 && !is.na(ref_env[[nc]][nearest])) {
        result[[nc]][i] <- ref_env[[nc]][nearest]
      }
    }
  }
  result
}

message("[utils_spatial_v4] loaded")
