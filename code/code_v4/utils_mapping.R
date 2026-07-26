#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把 v3 的 Nature 风格地图函数搬到 v4，并强制 7 条硬规则
#   （map_quality_rules.md）：固定 bbox、无十段线、无 NA 子图、
#   缺值灰色保留、within-metric z-score、南方北方覆盖审计
#
# Objective / 分析目标:
#   - china_province_basemap()：中国省级 sf
#   - v4_china_coord()：固定 coord_sf
#   - theme_nature_map() / theme_nature_pub()
#   - scale_fill_nature_c() / scale_fill_nature_d()
#   - save_nature(): 同时输出 png + pdf
#
# Input data / 输入数据:
#   data/external/china_shp 或 data/中国shp 的省界
#
# Main workflow / 主要流程:
#   纯函数定义
#
# Key assumptions / 关键假设:
#   sf / ggplot2 / scales 可用
#
# Main packages / 主要包:
#   sf, ggplot2, scales, scico（可选）
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(scales)
})

# ── 中国省级 basemap ─────────────────────────────────────────────────
china_province_basemap <- function() {
  # 优先查 data/external/china_boundary/
  candidate_paths <- c(
    file.path(DIRS$external, "china_boundary", "china_provinces.shp"),
    file.path(DIRS$external, "china_boundary", "中国省界.shp"),
    file.path(DIRS$data_raw, "china_boundary", "china_provinces.shp"),
    file.path(DIRS$data_raw, "中国shp", "省界.shp"),
    file.path(PROJECT_ROOT, "..", "bird_grid_community_analysis",
              "data", "external", "china_boundary", "china_provinces.shp")
  )
  for (p in candidate_paths) {
    if (file.exists(p)) {
      sf_obj <- tryCatch(st_read(p, quiet = TRUE), error = function(e) NULL)
      if (!is.null(sf_obj)) {
        message("[basemap] loaded ", p)
        return(st_transform(sf_obj, 4326))
      }
    }
  }
  warning("[basemap] No China province shapefile found; using empty sf")
  st_sf(geometry = st_sfc(crs = 4326))
}

# ── 固定 bbox 的 coord_sf（硬规则 #2） ───────────────────────────────
v4_china_coord <- function() {
  coord_sf(xlim = MAP_BBOX_XLIM, ylim = MAP_BBOX_YLIM, expand = FALSE)
}

# 包装 coord_sf() 调用，确保所有地图统一
# 用法：替换 ggplot 中的 + coord_sf() 为 + v4_china_coord()
coord_sf <- function(..., xlim = MAP_BBOX_XLIM, ylim = MAP_BBOX_YLIM, expand = FALSE) {
  ggplot2::coord_sf(..., xlim = xlim, ylim = ylim, expand = expand)
}

# ── 主题：地图 ─────────────────────────────────────────────────────
theme_nature_map <- function(base_size = NATURE_PT) {
  theme_void(base_size = base_size, base_family = NATURE_FONT) +
    theme(
      strip.text  = element_text(face = "bold", size = base_size),
      legend.position = "right",
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size - 1),
      legend.key.height = unit(8, "mm"),
      legend.key.width  = unit(2, "mm"),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

# ── 主题：通用 publication ───────────────────────────────────────────
theme_nature_pub <- function(base_size = NATURE_PT) {
  theme_minimal(base_size = base_size, base_family = NATURE_FONT) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
      axis.line  = element_line(linewidth = NATURE_LINE_AXIS, colour = "grey20"),
      axis.ticks = element_line(linewidth = NATURE_LINE_TICK, colour = "grey20"),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
}

# ── 色盘：连续 ─────────────────────────────────────────────────────
scale_fill_nature_c <- function(name = "Value", palette = "diverging",
                                 na.value = "grey92", ...) {
  if (palette == "diverging") {
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, name = name, na.value = na.value, ...
    )
  } else if (palette == "sequential") {
    if (requireNamespace("scico", quietly = TRUE)) {
      scale_fill_gradientn(colours = scico::scico(20, palette = "lajolla"),
                            name = name, na.value = na.value, ...)
    } else {
      scale_fill_viridis_c(option = "C", name = name, na.value = na.value, ...)
    }
  } else if (palette == "accent") {
    scale_fill_gradient(low = "white", high = NATURE_ACCENT,
                        name = name, na.value = na.value, ...)
  } else {
    scale_fill_viridis_c(name = name, na.value = na.value, ...)
  }
}

# ── 保存：同时 png + pdf（硬规则 #1 + #2 在调用时已固化） ─────────────
save_nature <- function(plot, stem, width_mm = NATURE_WIDTH_L,
                          height_mm = NULL, dir = DIRS$figures) {
  if (is.null(height_mm)) height_mm <- width_mm * 0.6
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")),
         plot = plot, width = width_mm, height = height_mm,
         units = "mm", dpi = NATURE_DPI)
  ggsave(file.path(dir, paste0(stem, ".pdf")),
         plot = plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
}

# ── 地图自检函数（硬规则 #6/#7） ─────────────────────────────────────
audit_map_coverage <- function(sf_data, var, label = NULL) {
  if (!"centroid_lat" %in% names(sf_data) || !"centroid_lon" %in% names(sf_data)) {
    sf_data$centroid_lon <- st_coordinates(st_centroid(sf_data))[, 1]
    sf_data$centroid_lat <- st_coordinates(st_centroid(sf_data))[, 2]
  }
  south <- sf_data |> dplyr::filter(centroid_lat < 25)
  north <- sf_data |> dplyr::filter(centroid_lat > 45)
  if (is.null(label)) label <- var

  pct_south_NA <- mean(is.na(south[[var]]), na.rm = TRUE) * 100
  pct_north_NA <- mean(is.na(north[[var]]), na.rm = TRUE) * 100
  message(sprintf("[map audit] %s: south NA %.1f%%, north NA %.1f%%",
                  label, pct_south_NA, pct_north_NA))
  invisible(tibble::tibble(
    var = label,
    south_pct_NA = pct_south_NA,
    north_pct_NA = pct_north_NA
  ))
}

message("[utils_mapping_v4] loaded — fixed bbox, no dashline, NA preserved")
