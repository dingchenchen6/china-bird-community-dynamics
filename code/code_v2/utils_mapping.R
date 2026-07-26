## utils_mapping.R
## v2 统一出图工具：地图 + 主题 + 调色板。
## 设计目标：
##   1. 修复 v1 地图"南北截断 + 鹰眼图错觉"的根因——显式 coord_sf 限定 bbox。
##   2. 永远不画十段线 / 鹰眼图（用户要求）。
##   3. 所有 v2 图统一字体、统一调色板、统一边距。

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(scales)
  library(scico)
})

## --- 1. 全局视觉常量 -------------------------------------------------------

V2_FONT_FAMILY <- "Helvetica"

# 中国全境 bbox：含南海诸岛 + 十段线（标准制图惯例）。鹰眼图禁用。
V2_CHINA_BBOX <- list(
  xlim = c(73, 135),
  ylim = c(3,  54)
)
# 仅大陆主体的 bbox（如确需聚焦时手动传入）
V2_MAINLAND_BBOX <- list(
  xlim = c(73, 135),
  ylim = c(18, 54)
)

V2_PALETTES <- list(
  diverging   = scico::scico(11, palette = "vik"),
  sequential  = scico::scico(11, palette = "lajolla", direction = 1),
  rich        = scico::scico(11, palette = "batlow"),
  qualitative = c("#0E5A78", "#C9784A", "#3C8C5A", "#8B2E1E", "#5A4E96",
                  "#D8A23A", "#3A6F8F", "#A85B8B", "#4F6B3F", "#7C7C7C")
)

## --- 2. 主题 ---------------------------------------------------------------

theme_v2_pub <- function(base_size = 11, base_family = V2_FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title          = element_text(face = "bold", size = base_size + 2,
                                         hjust = 0, margin = margin(b = 4)),
      plot.subtitle       = element_text(size = base_size, colour = "grey25",
                                         hjust = 0, margin = margin(b = 8)),
      plot.title.position = "plot",
      plot.caption        = element_text(size = base_size - 2, colour = "grey45",
                                         hjust = 1, margin = margin(t = 6)),
      panel.grid.minor    = element_blank(),
      panel.grid.major    = element_line(colour = "grey92", linewidth = 0.25),
      strip.text          = element_text(face = "bold", size = base_size - 0.5,
                                         margin = margin(t = 4, b = 4)),
      strip.background    = element_blank(),
      legend.position     = "top",
      legend.title        = element_text(size = base_size - 1, face = "plain"),
      legend.text         = element_text(size = base_size - 1.5),
      legend.key.height   = grid::unit(0.32, "cm"),
      legend.key.width    = grid::unit(1.0,  "cm"),
      axis.title          = element_text(size = base_size - 0.5),
      axis.text           = element_text(size = base_size - 1.5, colour = "grey25"),
      plot.margin         = margin(8, 10, 8, 10)
    )
}

theme_v2_map <- function(base_size = 11, base_family = V2_FONT_FAMILY) {
  theme_void(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title          = element_text(face = "bold", size = base_size + 2,
                                         hjust = 0, margin = margin(b = 4)),
      plot.subtitle       = element_text(size = base_size, colour = "grey25",
                                         hjust = 0, margin = margin(b = 8)),
      plot.title.position = "plot",
      plot.caption        = element_text(size = base_size - 2, colour = "grey45",
                                         hjust = 1, margin = margin(t = 6)),
      strip.text          = element_text(face = "bold", size = base_size - 0.5,
                                         margin = margin(t = 4, b = 4)),
      strip.background    = element_blank(),
      legend.position     = "top",
      legend.title        = element_text(size = base_size - 1, face = "plain"),
      legend.text         = element_text(size = base_size - 1.5),
      legend.key.height   = grid::unit(0.30, "cm"),
      legend.key.width    = grid::unit(1.4,  "cm"),
      plot.margin         = margin(6, 10, 6, 10)
    )
}

## --- 3. 调色 / scale 包装 --------------------------------------------------

scale_fill_v2_diverging <- function(name = NULL, limits = NULL, midpoint = 0,
                                    na.value = "grey92", oob = scales::squish, ...) {
  ggplot2::scale_fill_gradientn(
    name     = name,
    colours  = V2_PALETTES$diverging,
    values   = scales::rescale(c(-1, -0.4, 0, 0.4, 1)),
    rescaler = function(x, to = c(0, 1), from = NULL) {
      lim <- if (is.null(limits)) range(x, na.rm = TRUE) else limits
      hw  <- max(abs(lim - midpoint))
      scales::rescale(x, to = to, from = c(midpoint - hw, midpoint + hw))
    },
    limits   = limits,
    na.value = na.value,
    oob      = oob,
    ...
  )
}

scale_fill_v2_sequential <- function(name = NULL, limits = NULL,
                                     na.value = "grey92", oob = scales::squish,
                                     direction = 1, palette = "lajolla", ...) {
  ggplot2::scale_fill_gradientn(
    name     = name,
    colours  = scico::scico(11, palette = palette, direction = direction),
    limits   = limits,
    na.value = na.value,
    oob      = oob,
    ...
  )
}

## --- 4. 中国底图图层（共享，统一所有地图） ---------------------------------

## strip_nansha_inset: 中国官方 boundary shapefile 在福建/广东/海南/台湾等
## 南方省份的 MULTIPOLYGON 内通常内嵌一个"南海诸岛缩略图（鹰眼图）"
## 子多边形，位于经度 ~124-137°E、纬度 ~18-30°N 的东海空白区。
## 这个函数把所有 centroid 落在那个矩形里的子多边形剔除，保留真实国土。
strip_nansha_inset <- function(sf_obj,
                                excl_xmin = 124, excl_xmax = 137,
                                excl_ymin = 18,  excl_ymax = 30) {
  if (is.null(sf_obj) || nrow(sf_obj) == 0) return(sf_obj)
  geom_types <- as.character(sf::st_geometry_type(sf_obj))
  if (!any(grepl("POLYGON", geom_types))) return(sf_obj)
  parts <- suppressWarnings(sf::st_cast(sf_obj, "POLYGON", warn = FALSE))
  ctr <- suppressWarnings(sf::st_centroid(sf::st_geometry(parts)))
  cc  <- sf::st_coordinates(ctr)
  in_inset <- cc[, 1] >= excl_xmin & cc[, 1] <= excl_xmax &
              cc[, 2] >= excl_ymin & cc[, 2] <= excl_ymax
  parts_keep <- parts[!in_inset, , drop = FALSE]
  # 重组回原 feature（按属性 group_by 后 union）
  if (nrow(parts_keep) == 0) return(sf_obj)
  attr_cols <- setdiff(names(parts_keep), attr(parts_keep, "sf_column"))
  if (length(attr_cols) > 0) {
    out <- parts_keep |>
      dplyr::group_by(dplyr::across(dplyr::all_of(attr_cols))) |>
      dplyr::summarise(.groups = "drop")
  } else {
    out <- parts_keep
  }
  sf::st_make_valid(out)
}

load_china_layers <- function(china_boundary_path, province_line_path = NULL,
                              ten_dash_path = NULL,
                              strip_inset = TRUE,
                              province_lines_from_polygons = FALSE) {
  out <- list()
  out$china <- suppressWarnings(sf::st_read(china_boundary_path, quiet = TRUE)) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  if (isTRUE(strip_inset)) out$china <- strip_nansha_inset(out$china)
  # 主路径：直接从 省.shp 多边形派生 province boundaries（保证全部省-省界齐全）
  if (isTRUE(province_lines_from_polygons)) {
    out$province <- suppressWarnings(sf::st_boundary(out$china))
  } else if (!is.null(province_line_path) && file.exists(province_line_path)) {
    out$province <- suppressWarnings(sf::st_read(province_line_path, quiet = TRUE)) |>
      sf::st_make_valid() |>
      sf::st_transform(4326)
    if (isTRUE(strip_inset)) {
      bb <- sf::st_as_sfc(sf::st_bbox(c(xmin = 124, xmax = 137,
                                          ymin = 18,  ymax = 30),
                                       crs = sf::st_crs(out$province)))
      out$province <- out$province[lengths(sf::st_intersects(out$province, bb)) == 0, ]
    }
  } else {
    out$province <- NULL
  }
  # 十段线（南海主权标识；标准中国地图必加；与"鹰眼图"是两回事）
  if (is.null(ten_dash_path)) {
    P <- v2_paths()
    ten_dash_path <- P$ten_dash_line
  }
  if (!is.null(ten_dash_path) && file.exists(ten_dash_path)) {
    out$ten_dash <- suppressWarnings(sf::st_read(ten_dash_path, quiet = TRUE)) |>
      sf::st_make_valid() |>
      sf::st_transform(4326)
  } else {
    out$ten_dash <- NULL
  }
  out
}

## china_map_layers: 一组叠加在数据 geom_sf 之上的省界 + 国界图层。
## NOTE: 永远不画十段线 / inset。需要在论文里声明南海诸岛归属时用 caption / 文字注。
china_map_layers <- function(china_layers, country_colour = "#1F1F1F",
                             province_colour = "#666666",
                             country_lwd = 0.32, province_lwd = 0.18,
                             province_alpha = 0.55,
                             with_ten_dash = TRUE,
                             ten_dash_colour = "#1F1F1F",
                             ten_dash_lwd = 0.45) {
  layers <- list()
  if (!is.null(china_layers$province)) {
    layers <- c(layers, list(
      ggplot2::geom_sf(
        data = china_layers$province, inherit.aes = FALSE,
        colour = scales::alpha(province_colour, province_alpha),
        fill = NA, linewidth = province_lwd
      )
    ))
  }
  layers <- c(layers, list(
    ggplot2::geom_sf(
      data = china_layers$china, inherit.aes = FALSE,
      colour = country_colour, fill = NA, linewidth = country_lwd
    )
  ))
  if (isTRUE(with_ten_dash) && !is.null(china_layers$ten_dash)) {
    layers <- c(layers, list(
      ggplot2::geom_sf(
        data = china_layers$ten_dash, inherit.aes = FALSE,
        colour = scales::alpha(ten_dash_colour, 0.85),
        fill = NA, linewidth = ten_dash_lwd, linetype = "solid"
      )
    ))
  }
  layers
}

## v2_china_coord: 强制 coord_sf 到统一 bbox，禁用 expand。
##   默认 ylim=(3,54) 含十段线整段；如需仅大陆，传 ylim=V2_MAINLAND_BBOX$ylim。
v2_china_coord <- function(xlim = V2_CHINA_BBOX$xlim, ylim = V2_CHINA_BBOX$ylim,
                           crs = 4326) {
  ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = crs)
}

## --- 5. 标准网格地图骨架 ---------------------------------------------------

## build_grid_map: 给定 grid_sf（sf 数据，含字段 fill_var）+ 中国底图层，
## 生成一张干净的中国底图 + 网格 fill 图。
build_grid_map <- function(grid_sf, fill_var, china_layers,
                           palette_type = c("diverging", "sequential"),
                           palette = NULL,
                           legend_title = NULL,
                           limits = NULL,
                           midpoint = 0,
                           na_value = "grey92",
                           grid_outline_alpha = 0.10,
                           base_size = 11) {
  palette_type <- match.arg(palette_type)
  p <- ggplot2::ggplot(grid_sf) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = .data[[fill_var]]),
      colour = scales::alpha("white", grid_outline_alpha),
      linewidth = 0.05
    )
  p <- p + china_map_layers(china_layers)
  if (palette_type == "diverging") {
    p <- p + scale_fill_v2_diverging(
      name = legend_title, limits = limits, midpoint = midpoint, na.value = na_value
    )
  } else {
    p <- p + scale_fill_v2_sequential(
      name = legend_title, limits = limits, na.value = na_value,
      palette = palette %||% "lajolla"
    )
  }
  p + v2_china_coord() + theme_v2_map(base_size = base_size)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## --- 6. 统一保存（PNG + PDF） ---------------------------------------------

save_dual <- function(plot, stem, width = 9, height = 6, dpi = 360) {
  P <- v2_paths()
  ggplot2::ggsave(file.path(P$figures_v2, paste0(stem, ".png")),
                  plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(file.path(P$figures_v2, paste0(stem, ".pdf")),
                  plot, width = width, height = height, bg = "white")
  invisible(NULL)
}

## save_dual_v3: 同 save_dual 但写到 figures_v3/
save_dual_v3 <- function(plot, stem, width = 9, height = 6, dpi = 360) {
  P <- v2_paths()
  ggplot2::ggsave(file.path(P$figures_v3, paste0(stem, ".png")),
                  plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(file.path(P$figures_v3, paste0(stem, ".pdf")),
                  plot, width = width, height = height, bg = "white")
  invisible(NULL)
}

## save_triple: PNG + PDF + 单页可编辑 PPTX（rvg/officer，矢量层在 PPT 中可拖动）
save_triple <- function(plot, stem, width = 9, height = 6, dpi = 360) {
  save_dual(plot, stem, width = width, height = height, dpi = dpi)
  P <- v2_paths()
  if (requireNamespace("officer", quietly = TRUE) &&
      requireNamespace("rvg",     quietly = TRUE)) {
    pptx_path <- file.path(P$figures_v2, paste0(stem, ".pptx"))
    doc <- officer::read_pptx()
    doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
    doc <- officer::ph_with(
      doc, value = rvg::dml(ggobj = plot),
      location = officer::ph_location(left = 0.4, top = 0.4,
                                       width = width, height = height)
    )
    print(doc, target = pptx_path)
  }
  invisible(NULL)
}

## save_pptx_deck: 把多张图合并到同一个 PPTX deck。
##  embed = "png": 嵌入 PNG（快、稳定，不可在 PPT 内编辑）；
##  embed = "vector": 用 rvg::dml 嵌入矢量（可在 PPT 内编辑，复杂 facet 易崩）；
##  embed = "auto": 先尝试 vector，崩则 fallback PNG。
save_pptx_deck <- function(plots_named, deck_stem, width = 12, height = 7,
                            embed = c("auto", "png", "vector"), dpi = 220) {
  if (!requireNamespace("officer", quietly = TRUE)) return(invisible(NULL))
  embed <- match.arg(embed)
  P <- v2_paths()
  doc <- officer::read_pptx()
  for (nm in names(plots_named)) {
    pl <- plots_named[[nm]]
    if (is.null(pl)) next
    doc <- officer::add_slide(doc, layout = "Title and Content",
                               master = "Office Theme")
    doc <- officer::ph_with(doc, value = nm,
                             location = officer::ph_location_type(type = "title"))
    inserted <- FALSE
    if (embed %in% c("auto", "vector") &&
        requireNamespace("rvg", quietly = TRUE)) {
      ok <- tryCatch({
        doc <- officer::ph_with(doc, value = rvg::dml(ggobj = pl),
                                 location = officer::ph_location(left = 0.3, top = 1.1,
                                                                 width = width,
                                                                 height = height))
        TRUE
      }, error = function(e) { message("  vector embed fail (", nm, "): ",
                                       conditionMessage(e)); FALSE })
      if (ok) { inserted <- TRUE }
    }
    if (!inserted) {
      png_tmp <- tempfile(fileext = ".png")
      ggplot2::ggsave(png_tmp, pl, width = width, height = height,
                       dpi = dpi, bg = "white")
      doc <- officer::ph_with(doc, value = officer::external_img(png_tmp,
                                                                  width = width,
                                                                  height = height),
                               location = officer::ph_location(left = 0.3, top = 1.1,
                                                               width = width,
                                                               height = height))
    }
  }
  pptx_path <- file.path(P$figures_v2, paste0(deck_stem, ".pptx"))
  print(doc, target = pptx_path)
  invisible(pptx_path)
}
