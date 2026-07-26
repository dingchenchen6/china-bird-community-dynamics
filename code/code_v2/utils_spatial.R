## utils_spatial.R
## v2 共享空间工具：网格构建、栅格抽取。
## 抽自旧脚本 prepare_dynamic_occupancy_environment.R 与 run_bird_dynamic_occupancy_analysis.R
## 中重复出现的 helper，改用 sf 0.9+ / terra 1.7+ 安全调用。

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tibble)
  library(purrr)
})

align_extent <- function(xmin, xmax, ymin, ymax, resolution) {
  terra::ext(
    floor(xmin / resolution) * resolution,
    ceiling(xmax / resolution) * resolution,
    floor(ymin / resolution) * resolution,
    ceiling(ymax / resolution) * resolution
  )
}

## create_china_grid: 在等积投影下生成 100km 网格，clip 到中国边界。
## 返回 sf（WGS84），含字段 grid_cell, grid_id, area_km2, centroid_lon, centroid_lat。
create_china_grid <- function(china_boundary, resolution_m = 1e5) {
  if (sf::st_crs(china_boundary)$epsg %in% c(4326, NA)) {
    china_proj <- sf::st_transform(china_boundary, 3857)
  } else {
    china_proj <- china_boundary
  }
  bb <- sf::st_bbox(china_proj)
  template <- terra::rast(
    ext = align_extent(bb$xmin, bb$xmax, bb$ymin, bb$ymax, resolution_m),
    resolution = resolution_m,
    crs = sf::st_crs(china_proj)$wkt
  )
  terra::values(template) <- seq_len(terra::ncell(template))
  poly <- terra::as.polygons(template, values = TRUE, na.rm = FALSE) |>
    sf::st_as_sf() |>
    dplyr::rename(grid_cell = lyr.1)
  grid_sf <- suppressWarnings(sf::st_intersection(poly, china_proj)) |>
    sf::st_make_valid()
  if (any(sf::st_geometry_type(grid_sf) == "GEOMETRYCOLLECTION")) {
    grid_sf <- sf::st_collection_extract(grid_sf, "POLYGON", warn = FALSE)
  }
  grid_sf <- grid_sf |>
    dplyr::mutate(grid_id = sprintf("G%05d", grid_cell)) |>
    dplyr::group_by(grid_cell, grid_id) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid() |>
    sf::st_cast("MULTIPOLYGON", warn = FALSE) |>
    dplyr::mutate(area_km2 = as.numeric(sf::st_area(geometry) / 1e6)) |>
    dplyr::filter(area_km2 > 1)
  centroids <- sf::st_centroid(grid_sf) |>
    sf::st_transform(4326) |>
    sf::st_coordinates() |>
    tibble::as_tibble() |>
    dplyr::rename(centroid_lon = X, centroid_lat = Y) |>
    dplyr::mutate(grid_cell = grid_sf$grid_cell)
  grid_sf |>
    dplyr::left_join(centroids, by = "grid_cell") |>
    sf::st_transform(4326)
}

extract_mean_tbl <- function(r, grid_vect, prefix = NULL) {
  if (is.null(r) || terra::nlyr(r) == 0) return(tibble::tibble(grid_cell = integer()))
  out <- terra::extract(r, grid_vect, fun = mean, na.rm = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(grid_cell = ID)
  if (!is.null(prefix)) {
    value_cols <- setdiff(names(out), "grid_cell")
    names(out)[match(value_cols, names(out))] <- paste0(prefix, value_cols)
  }
  out
}

first_existing_path <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  normalizePath(hit[1], mustWork = FALSE)
}

## load_or_build_v2_grid: 优先从 derived_v2 缓存读取；否则现场构建并保存。
load_or_build_v2_grid <- function(china_boundary_path, resolution_m = 1e5,
                                  cache_path = NULL, force = FALSE) {
  if (!is.null(cache_path) && file.exists(cache_path) && !force) {
    return(readRDS(cache_path))
  }
  china <- suppressWarnings(sf::st_read(china_boundary_path, quiet = TRUE)) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  grid_sf <- create_china_grid(china, resolution_m = resolution_m)
  if (!is.null(cache_path)) saveRDS(grid_sf, cache_path)
  grid_sf
}
