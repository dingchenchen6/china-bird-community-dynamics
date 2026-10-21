#!/usr/bin/env Rscript
# ============================================================
# utils_figs_export_20261011.R
#
# Scientific question / 科学问题:
#   群落动态稿件的全部图件须在不改变任何作图数值的前提下，统一满足中国标准地图、
#   英文标注与原生可编辑 PowerPoint 的交付要求。
#   The community-dynamics figures must meet the China standard-map rules and be
#   delivered as natively editable PowerPoint, without changing any plotted value.
#
# Objective / 分析目标:
#   (1) 中国标准地图图层：省面.shp 按位置拆分为主图层与南海插图层（直接复用
#       code_v5_urban/utils_figs_pub.R 的 load_china_layers() 与 BAKED_INSET），
#       完整十段线只画在右下角南海诸岛插图中；
#   (2) 统一导出：PDF（矢量）+ 450 dpi PNG + 450 dpi LZW TIFF + 原生可编辑 PPTX
#       （rvg DrawingML，幻灯片尺寸等于图件尺寸）。
#   (1) China standard-map layers; (2) one export routine for PDF, PNG, TIFF, PPTX.
#
# Input data / 输入数据:
#   data/china_boundary/省面.shp, 十段线.shp; data/derived_v2/china_grid_100km_v2.rds
# Main workflow / 主要流程:
#   读取并拆分边界 → 主图裁切到固定 bbox → 显示用简化 → 国界外轮廓（去除省际缝隙）
#   → 南海插图 → 统一导出
#   load and split boundaries -> crop -> display simplification -> outline -> inset -> export
# Expected output / 预期输出: 函数与常量，不写文件 / functions only, no files written
# Key assumptions / 关键假设:
#   - 显示用简化容差 0.005°（约 0.5 km），在 30–110 mm 宽的地图上远小于 0.01 mm，肉眼不可见；
#     十段线不做简化。Display simplification of 0.005 deg is far below print resolution;
#     the ten-dash line is never simplified.
#   - 100 km 网格中落在烘焙插图副本上的格子（128.4°E 以东）是伪格子，不属于任何模型，
#     绘图时剔除（先断言其中没有模型格子）。
#     Grid cells on the baked-in inset copy (east of 128.4 E) are artefacts and are not drawn.
#   - 插图范围与位置使用协调人 2026-10-11 更正后的 utils_figs_pub.R：MAP_BB_SCS = c(105.5, 124, 2.5, 25.6)
#     （含台湾全岛与完整第十段线）、INSET_BOX = c(125.4, 135.8, 17.7, 28.56)（不遮挡 125.4°E 以西）。
#     Inset extent and box are the corrected values in utils_figs_pub.R (2026-10-11).
#   - v5 的拆分（BAKED_INSET 自 124.5°E 起）会把真实位置上的赤尾屿（台湾省，124.55°E, 25.87°N）
#     一并删除；本文件在拆分后把这一个真实岛屿加回主图层，并保留其所在网格。
#     The v5 split also deletes Chiwei Yu (a true-position island at 124.55 E); it is added back.
# Main packages / 主要包: sf, ggplot2, officer, rvg, ragg, xml2, terra
# Output directory / 输出路径: set by the calling scripts (figures_community_submission_20261011/)
# ============================================================

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(grid)
})
sf::sf_use_s2(FALSE)

# ── 项目根目录与输出根目录 / project and output roots ────────────────
PROJECT_ROOT <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(PROJECT_ROOT, "code_v3")))
  stop("Run from the project root (bird_dynamic_occupancy_analysis).", call. = FALSE)
FIG_OUT_ROOT <- file.path(PROJECT_ROOT, "figures_community_submission_20261011")

# ── 复用 v5 标准地图实现（独立环境，避免与各脚本的 theme 名冲突）──────
# Reuse the v5 standard-map implementation in its own environment.
V5 <- new.env()
sys.source(file.path(PROJECT_ROOT, "code_v5_urban", "utils_figs_pub.R"), envir = V5)
MAP_BB      <- V5$MAP_BB        # c(73, 136, 17.5, 54)，主图固定范围 / fixed main extent
MAP_BB_SCS  <- V5$MAP_BB_SCS    # c(105.5, 124, 2.5, 25.6)，南海插图范围（已更正）/ inset extent (corrected)
BAKED_INSET <- V5$BAKED_INSET   # 烘焙插图区 / baked-in inset zone
# 2026-10-13：共享工具已把起点改为 125.0°E，赤尾屿留在 prov 中；下面的"补回"步骤因此不再选中任何部分，不会重复绘制。
# 2026-10-13: the shared utility now starts the zone at 125.0 E, so Chiwei Yu stays in prov and the add-back below selects nothing.
INSET_BOX_PUB <- V5$INSET_BOX   # c(125.4, 135.8, 17.7, 28.56)，主图中的插图位置（已更正）/ corrected box
# 确认使用的是更正后的插图参数（台湾全岛与第十段线在插图内；125.4°E 以西不被遮挡）
# Guard: the corrected inset parameters must be in place.
if (MAP_BB_SCS[4] < 24.7 || MAP_BB_SCS[2] < 122.9 || INSET_BOX_PUB[1] < 125.4)
  stop("utils_figs_pub.R still has the old inset extent/box; expected the 2026-10-11 correction.")
TRUE_PART_MAX_LON <- 125.0      # 烘焙副本始于 128.4°E；此经度以西的部分是真实岛屿 / west of this = true position
NA_GREY     <- "grey80"         # 缺数据格统一灰色 / one grey for cells without data
SIMPLIFY_TOL <- 0.005           # 显示用简化容差（度）/ display simplification (degrees)

bbox_sfc <- function(b) {
  # 由 c(xmin, xmax, ymin, ymax) 构造 bbox 多边形 / bbox polygon from a numeric extent
  sf::st_as_sfc(sf::st_bbox(c(xmin = b[1], xmax = b[2], ymin = b[3], ymax = b[4]), crs = 4326))
}

simplify_display <- function(x, tol = SIMPLIFY_TOL) {
  # 仅为显示简化几何（保持拓扑，不丢小岛）/ display-only simplification, topology kept
  suppressWarnings(sf::st_simplify(x, preserveTopology = TRUE, dTolerance = tol))
}

drop_holes <- function(geom) {
  # 去除并集多边形内部的环：省面各省边界不完全重合，并集会留下上千个省际缝隙环。
  # Drop interior rings: the province union leaves >1,500 sliver rings along borders.
  polys <- suppressWarnings(sf::st_cast(sf::st_sfc(geom, crs = 4326), "POLYGON"))
  sf::st_sfc(lapply(polys, function(p) sf::st_polygon(list(p[[1]]))), crs = 4326)
}

restored_true_parts <- function(root = PROJECT_ROOT) {
  # 读原始省面，取落在 BAKED_INSET 内、但位于烘焙副本以西（< 125°E）的真实岛屿（即赤尾屿）。
  # True-position parts inside BAKED_INSET but west of the baked copy (i.e. Chiwei Yu).
  raw <- sf::st_transform(sf::st_make_valid(sf::st_read(
    file.path(root, "data", "china_boundary", "省面.shp"), quiet = TRUE)), 4326)
  parts <- suppressWarnings(sf::st_cast(sf::st_cast(raw, "MULTIPOLYGON"), "POLYGON"))
  bb <- t(vapply(sf::st_geometry(parts), function(g) as.numeric(sf::st_bbox(g)), numeric(4)))
  keep <- bb[, 1] >= BAKED_INSET[1] & bb[, 3] < TRUE_PART_MAX_LON &
    bb[, 2] >= BAKED_INSET[3] & bb[, 4] <= BAKED_INSET[4]
  parts[keep, ]
}

china_layers_pub <- function(root = PROJECT_ROOT, tol = SIMPLIFY_TOL) {
  # 1) 与 v5 完全相同的拆分 / identical split to v5
  L <- suppressMessages(suppressWarnings(V5$load_china_layers(root)))
  # 1b) 加回被拆分误删的真实岛屿（赤尾屿）/ add back the true-position island removed by the split
  extra <- restored_true_parts(root)
  if (nrow(extra)) L$prov <- rbind(L$prov, extra[, names(L$prov)])
  # 2) 主图层裁切到固定范围：省面.shp 中 17.5°N 以南的南海岛礁只在插图中出现
  #    Crop the main layer to the fixed extent; SCS islands appear in the inset only.
  main <- suppressWarnings(sf::st_crop(L$prov, sf::st_bbox(bbox_sfc(MAP_BB))))
  outline <- drop_holes(suppressMessages(sf::st_union(main)))
  # 3) 插图层：十段线裁切到插图范围（等同于绘图时的面板裁切），不做简化
  #    Inset: crop the dash line to the inset extent (same as panel clipping); no simplification.
  dash <- suppressWarnings(sf::st_crop(L$dash, sf::st_bbox(bbox_sfc(MAP_BB_SCS))))
  list(prov = simplify_display(main, tol), outline = simplify_display(outline, tol),
       prov_scs = simplify_display(L$prov_scs, tol), dash = dash, restored = extra)
}

grid_for_maps <- function(grid_sf, model_cells, tol = SIMPLIFY_TOL) {
  # 剔除烘焙副本上的伪格子（BAKED_INSET 内且位于 125°E 以东）；赤尾屿所在真实格子保留。
  # 若其中有任何模型格子则报错（不允许静默丢数据）。
  # Drop artefact cells on the baked copy; keep the true Chiwei Yu cell; never drop model cells.
  zone <- BAKED_INSET; zone[1] <- TRUE_PART_MAX_LON
  hit <- lengths(sf::st_intersects(grid_sf, bbox_sfc(zone))) > 0
  bad <- intersect(grid_sf$grid_cell[hit], model_cells)
  if (length(bad)) stop("Model grid cells inside the baked-in inset zone: ", paste(bad, collapse = ", "))
  out <- grid_sf[!hit, ]
  sf::st_geometry(out) <- simplify_display(sf::st_geometry(out), tol)
  attr(out, "dropped_artefact_cells") <- grid_sf$grid_cell[hit]
  out
}

# ── 南海诸岛插图 / South China Sea inset ─────────────────────────────
scs_inset_pub <- function(L, line_col = V5$MAP_LINE, line_lw = 0.08, dash_col = V5$MAP_DASH,
                          dash_lw = 0.3, frame_col = "grey45", frame_lw = 0.25, fill_layers = NULL) {
  # 与 v5 scs_inset() 同范围、同图层顺序；另加白底，保证插图下方不透出任何要素。
  # Same extent and layer order as v5 scs_inset(), plus a white background.
  p <- ggplot()
  for (l in fill_layers) p <- p + l
  p + geom_sf(data = L$prov_scs, fill = NA, colour = line_col, linewidth = line_lw) +
    geom_sf(data = L$dash, colour = dash_col, linewidth = dash_lw) +
    coord_sf(xlim = MAP_BB_SCS[1:2], ylim = MAP_BB_SCS[3:4], expand = FALSE, crs = sf::st_crs(4326), datum = NA) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = frame_col, linewidth = frame_lw))
}

scs_inset_layer_pub <- function(inset, box = INSET_BOX_PUB) {
  # annotation_custom 嵌入主图坐标系；不用 patchwork::inset_element（会破坏面板编号）。
  # annotation_custom keeps panel tagging intact (inset_element would not).
  annotation_custom(ggplotGrob(inset), xmin = box[1], xmax = box[2], ymin = box[3], ymax = box[4])
}

nodata_key <- function(label = "No data", fill = NA_GREY, x = 100, y = 35, order = 2) {
  # 用一个不可见的点层生成"无数据"图例键（连续色标无法显示 na.value）。
  # An invisible point layer that only produces a 'No data' legend key.
  list(
    geom_point(data = data.frame(x = x, y = y, key = label), aes(x = x, y = y, shape = key),
               size = 0.01, alpha = 0, inherit.aes = FALSE),
    scale_shape_manual(values = stats::setNames(22, label), name = NULL,
                       guide = guide_legend(order = order, override.aes = list(
                         size = 2.6, alpha = 1, fill = fill, colour = NA, stroke = 0),
                         theme = theme(legend.key.width = unit(3, "mm"), legend.key.height = unit(3, "mm"))))
  )
}

blend_white <- function(col, a) {
  # 把透明度换算成白底上的等效不透明色：相邻格子可描同色边，消除抗锯齿白缝。
  # Convert alpha to the equivalent opaque colour on white, so cells can be stroked
  # in their own colour (removes anti-aliasing seams) without double-alpha borders.
  a[is.na(a)] <- 1
  m <- grDevices::col2rgb(col) / 255
  out <- m * rep(a, each = 3) + (1 - rep(a, each = 3))
  grDevices::rgb(out[1, ], out[2, ], out[3, ])
}

cell_pixels <- function(grid_sf, res_deg, bb = MAP_BB) {
  # PPTX 用：把格子多边形栅格化为规则经纬度像元（像元中心落在哪个格子即取该格子编号）。
  # For PPTX only: rasterise cell polygons to a regular lon/lat grid of cell ids.
  # 不带 CRS 运算（只需格子编号；本机 terra 找不到 proj.db），坐标仍为经纬度。
  # CRS-free on purpose (only cell ids are needed; terra's PROJ database is unavailable here).
  g <- sf::st_set_crs(grid_sf[, "grid_cell"], NA)
  r <- terra::rast(xmin = bb[1], xmax = bb[2], ymin = bb[3], ymax = bb[4],
                   resolution = res_deg, crs = "")
  rr <- terra::rasterize(terra::vect(g), r, field = "grid_cell")
  df <- terra::as.data.frame(rr, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "grid_cell")
  # 像元中心坐标取整到分辨率网格，避免浮点误差触发 geom_raster 的"间隔不均"位移
  # Snap pixel centres to the resolution grid (avoids geom_raster uneven-interval shifts).
  df$x <- round(df$x, 8); df$y <- round(df$y, 8)
  df
}

# ── 导出 / Export ────────────────────────────────────────────────────
pptx_native <- function(plot, path, width_mm, height_mm) {
  # rvg DrawingML：每条线、每个点、每段文字都是 PowerPoint 原生对象；
  # 幻灯片尺寸设为图件尺寸，图件以原始物理尺寸放在 (0, 0)，字号与投稿版一致。
  # Native DrawingML; slide size equals the figure size, figure placed at true size.
  w <- round(width_mm / 25.4, 9); h <- round(height_mm / 25.4, 9)
  doc <- officer::read_pptx()
  doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
  doc <- officer::ph_with(doc, rvg::dml(ggobj = plot, bg = "white"),
                          location = officer::ph_location(left = 0, top = 0, width = w, height = h))
  sz <- xml2::xml_find_first(doc$presentation$get(), "//p:sldSz")
  xml2::xml_set_attr(sz, "cx", as.character(round(w * 914400)))
  xml2::xml_set_attr(sz, "cy", as.character(round(h * 914400)))
  xml2::xml_set_attr(sz, "type", NULL)
  print(doc, target = path)
  invisible(path)
}

save_figure_set <- function(plot, dir, stem, width_mm, height_mm, ppt_plot = plot, dpi = 450,
                            pptx_mode = "native_vector") {
  # 四种格式一次写出 / write PDF, PNG, TIFF and PPTX in one call
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # 英寸取 9 位小数，消除 mm↔in 往返的浮点误差（否则 12.2 in 会变成 5489 而非 5490 像素）
  # Round inches to 9 decimals so mm->in round trips do not lose a pixel column.
  base <- file.path(dir, stem); w <- round(width_mm / 25.4, 9); h <- round(height_mm / 25.4, 9)
  ggsave(paste0(base, ".pdf"), plot, width = w, height = h, units = "in", device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".png"), plot, width = w, height = h, units = "in", dpi = dpi, bg = "white",
         device = ragg::agg_png)
  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = dpi,
                 background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()
  pptx_native(ppt_plot, paste0(base, ".pptx"), width_mm, height_mm)
  message(sprintf("[export] %s (%.0f x %.0f mm; PPTX %s)", stem, width_mm, height_mm, pptx_mode))
  data.frame(stem = stem, folder = basename(dir), width_mm = width_mm, height_mm = height_mm,
             pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"), tiff = paste0(stem, ".tiff"),
             pptx = paste0(stem, ".pptx"), pptx_mode = pptx_mode, stringsAsFactors = FALSE)
}

measure_text_mm <- function(x, size_pt, family = "Helvetica", face = 1) {
  # 在当前 Cairo 设备上量取文字宽度（mm），与 cairo_pdf 输出一致。
  # Text width in mm on the current Cairo device (matches the cairo_pdf output).
  vapply(x, function(s) grid::convertWidth(grid::grobWidth(grid::textGrob(
    s, gp = grid::gpar(fontfamily = family, fontsize = size_pt, fontface = face))), "mm", valueOnly = TRUE),
    numeric(1))
}

wrap_to_width <- function(text, width_mm, size_pt, family = "Helvetica", face = 1) {
  # 按实测宽度贪心断行（只在空格处断行，不改任何字词）。
  # Greedy word wrap on measured widths; breaks only at spaces, words unchanged.
  words <- strsplit(text, " ", fixed = TRUE)[[1]]
  lines <- character(0); cur <- ""
  for (wd in words) {
    cand <- if (nzchar(cur)) paste(cur, wd) else wd
    if (!nzchar(cur) || measure_text_mm(cand, size_pt, family, face) <= width_mm) cur <- cand
    else { lines <- c(lines, cur); cur <- wd }
  }
  lines <- c(lines, cur)
  if (length(lines) == 2L && length(words) > 2L) {
    # 两行时取最均衡的断点，避免孤词 / balanced two-line break (no orphan word)
    splits <- lapply(seq_len(length(words) - 1L), function(k)
      c(paste(words[seq_len(k)], collapse = " "), paste(words[-seq_len(k)], collapse = " ")))
    widest <- vapply(splits, function(s) max(measure_text_mm(s, size_pt, family, face)), numeric(1))
    ok <- widest <= width_mm
    if (any(ok)) lines <- splits[[which(ok)[which.min(widest[ok])]]]
  }
  lines
}

message("[utils_figs_export_20261011] loaded; inset box: ",
        paste(round(INSET_BOX_PUB, 2), collapse = ", "))
