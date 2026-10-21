#!/usr/bin/env Rscript
# ============================================================
# 24c_figure1_layout_20261011.R  (sourced by 24b_main_figures_standard_map_20261011.R)
#
# Scientific question / 科学问题:
#   Figure 1 说明证据链：观测数据 → 调查设计 → 主模型/广度扩展 → 后验生态证据，
#   以及两套模型各自支持哪些结论（推断边界）。
#   Figure 1: the evidence architecture and the inference boundary of the two models.
#
# Objective / 分析目标:
#   修复 Figure 1a 方框文字溢出：(1) 用 patchwork::free() 使 a 面板不再与 b 面板的
#   y 轴标签对齐（原先因此被压窄约 27 mm）；(2) 按实测文字宽度在空格处断行，
#   不改任何字词与数字；(3) 自动核查每一行都落在方框内（结果写入 qa/Figure_1a_text_fit.csv）。
#   Fix the text overflow in panel a: free the panel from axis alignment, wrap on measured
#   widths (words unchanged) and verify every line fits its box.
#
# Input data / 输入数据: 24b 中已读入的收敛诊断表 conv200 / conv500（仅用于源数据）
# Main workflow / 主要流程: 校准面板物理尺寸 → 实测断行 → 排版 → 核查 → 导出
#   calibrate panel size -> measured wrapping -> layout -> fit check -> export
# Expected output / 预期输出: main/Figure_1_inference_architecture.*; source_data/SourceData_Figure_1.csv;
#   qa/Figure_1a_text_fit.csv
# Key assumptions / 关键假设: 方框内文字与 24 号脚本完全相同（硬编码数值未改）；
#   b 面板代码逐字沿用。Box texts and panel b are identical to script 24.
# Main packages / 主要包: ggplot2, patchwork, ragg, png, grid
# Output directory / 输出路径: figures_community_submission_20261011/{main,source_data,qa}
# ============================================================

FIG1_W <- 183; FIG1_H <- 118
BODY_SIZE <- 2.05; TITLE_SIZE <- 2.45                 # ggplot 文字尺寸（mm），同 24 号脚本
BODY_PT <- BODY_SIZE * .pt; TITLE_PT <- TITLE_SIZE * .pt
PAD_MM <- 1.6; LINE_H <- 1.05; ITEM_GAP_MM <- 0.8      # 内边距、行距、条目间距 / padding, leading, gap
X_RANGE <- c(0, 100); Y_RANGE <- c(0, 38)

# 方框内容与 24 号脚本逐字一致；位置按 free() 后的整宽面板重新排布
# Box texts verbatim from script 24; positions re-spaced for the full-width panel.
fig1_boxes <- list(
  list(id = "Occurrence evidence", xmin = 0.6, xmax = 18.6, ymin = 9, ymax = 30, fill = "#E7F1F7",
       items = c("7.49 million records", "1,426 recorded species", "3.15 million breeding-season records")),
  list(id = "Survey design", xmin = 22.2, xmax = 40.2, ymin = 9, ymax = 30, fill = "#F0ECF7",
       items = c("1,247 100-km grid cells", "Five periods, 2000-2024", "Effort and duration covariates")),
  list(id = "Primary spatial model", xmin = 44.4, xmax = 66.4, ymin = 21, ymax = 36, fill = "#F6E8DD",
       items = c("200-species stMsPGOcc", "NNGP + temporal AR(1)", "4 chains; maximum R-hat 1.039")),
  list(id = "Breadth extension", xmin = 44.4, xmax = 66.4, ymin = 2, ymax = 17, fill = "#FBEEE6",
       items = c("500-species tMsPGOcc", "Temporal AR(1), non-spatial", "R-hat unavailable")),
  list(id = "Posterior ecological evidence", xmin = 70.6, xmax = 99.4, ymin = 9, ymax = 30, fill = "#E7F2EA",
       items = c("Taxonomic, functional and phylogenetic diversity",
                 "Temporal beta decomposition and species trends",
                 "Detection bias, traits and driver associations"))
)
fig1_arrows <- data.frame(x = c(19.1, 40.7, 40.7, 66.9, 66.9), xend = c(21.7, 43.9, 43.9, 70.1, 70.1),
                          y = c(19.5, 23, 16, 28, 9.5), yend = c(19.5, 28, 9.5, 23, 16))

# b 面板逐字沿用 24 号脚本 / panel b verbatim from script 24
scope_data <- tribble(
  ~model, ~target, ~status,
  "200-species spatial", "Spatial patterns", "Primary",
  "200-species spatial", "National trajectories", "Primary",
  "200-species spatial", "Species trends", "Primary",
  "200-species spatial", "Trait associations", "Complete-case",
  "500-species temporal", "Spatial patterns", "Descriptive only",
  "500-species temporal", "National trajectories", "Breadth support",
  "500-species temporal", "Species trends", "Breadth support",
  "500-species temporal", "Trait associations", "n = 200"
) |>
  mutate(
    target = factor(target, levels = c("Spatial patterns", "National trajectories", "Species trends", "Trait associations")),
    model = factor(model, levels = c("500-species temporal", "200-species spatial")),
    status_group = case_when(status == "Primary" ~ "Primary", str_detect(status, "Breadth") ~ "Breadth",
                             str_detect(status, "Descriptive") ~ "Descriptive", TRUE ~ "Limited")
  )
p1b <- ggplot(scope_data, aes(target, model, fill = status_group)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = status), size = 2.05, colour = "black", lineheight = 0.95) +
  scale_fill_manual(values = c(Primary = "#B8D8EB", Breadth = "#E7D5E2", Descriptive = "#E5E7E9", Limited = "#F6DFC2"), guide = "none") +
  labs(title = "Inference boundary", x = NULL, y = NULL) +
  theme_nature(6.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), plot.margin = margin(1, 2, 2, 2))

p1a_frame <- function(layers, calibrate = FALSE) {
  # a 面板骨架；calibrate = TRUE 时画一个品红边框用于量取面板物理尺寸
  # Panel-a skeleton; a magenta frame is drawn only for calibration.
  p <- ggplot() + coord_cartesian(xlim = X_RANGE, ylim = Y_RANGE, expand = FALSE, clip = "off")
  for (l in layers) p <- p + l
  if (calibrate) p <- p + annotate("rect", xmin = X_RANGE[1], xmax = X_RANGE[2], ymin = Y_RANGE[1],
                                   ymax = Y_RANGE[2], fill = NA, colour = "#FF00FF", linewidth = 0.3)
  p + labs(title = "Evidence architecture") +
    theme_void(base_family = "Helvetica", base_size = 6.5) +
    theme(plot.title = element_text(size = 7, face = "bold"), plot.tag = element_text(size = 8, face = "bold"),
          plot.margin = margin(2, 2, 1, 2))
}
fig1_compose <- function(p1a) {
  # free()：a 面板不再与 b 面板的坐标轴对齐，获得整幅宽度 / panel a no longer aligned to b's axis
  free(p1a) / p1b + plot_layout(heights = c(2.25, 1)) + plot_annotation(tag_levels = "a")
}

# 1) 校准：渲染一次，量取品红框像素范围 → 每个数据单位对应的毫米数
#    Calibrate: render once and convert the magenta frame to mm per data unit.
cal_png <- tempfile(fileext = ".png")
ragg::agg_png(cal_png, width = FIG1_W, height = FIG1_H, units = "mm", res = 300)
print(fig1_compose(p1a_frame(list(), calibrate = TRUE)))
invisible(grDevices::dev.off())
cal <- png::readPNG(cal_png)
mag <- cal[, , 1] > 0.9 & cal[, , 2] < 0.2 & cal[, , 3] > 0.9
MM_X <- diff(range(which(colSums(mag) > 0))) / 300 * 25.4 / diff(X_RANGE)
MM_Y <- diff(range(which(rowSums(mag) > 0))) / 300 * 25.4 / diff(Y_RANGE)
message(sprintf("[24c] panel a: %.2f mm per x unit, %.2f mm per y unit", MM_X, MM_Y))

# 2) 实测断行与排版（在 Cairo 设备上量宽）/ measured wrapping and layout on a Cairo device
grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 7, height = 5)
line_mm <- BODY_PT * LINE_H * 25.4 / 72
title_line_mm <- TITLE_PT * LINE_H * 25.4 / 72
box_layers <- list(); fit_rows <- list()
for (b in fig1_boxes) {
  usable <- (b$xmax - b$xmin) * MM_X - 2 * PAD_MM
  mid <- (b$xmin + b$xmax) / 2
  t_lines <- wrap_to_width(b$id, usable, TITLE_PT, face = 2)
  i_lines <- lapply(b$items, wrap_to_width, width_mm = usable, size_pt = BODY_PT)
  title_y <- b$ymax - 2.1                                     # 同 24 号脚本 / as in script 24
  title_bottom <- title_y * MM_Y - length(t_lines) * title_line_mm / 2
  n_l <- lengths(i_lines)
  block <- sum(n_l) * line_mm + (length(n_l) - 1) * ITEM_GAP_MM
  centre <- ((b$ymin + b$ymax) / 2 - 0.8) * MM_Y              # 原正文中心 / original body centre
  top <- min(centre + block / 2, title_bottom - 1.2)          # 与标题保持间距 / keep clear of title
  y_items <- numeric(length(n_l)); cur <- top
  for (i in seq_along(n_l)) { y_items[i] <- cur - n_l[i] * line_mm / 2; cur <- cur - n_l[i] * line_mm - ITEM_GAP_MM }
  bottom <- cur + ITEM_GAP_MM
  box_layers <- c(box_layers, list(
    annotate("rect", xmin = b$xmin, xmax = b$xmax, ymin = b$ymin, ymax = b$ymax,
             fill = b$fill, colour = "#6F7B84", linewidth = 0.35),
    annotate("text", x = mid, y = title_y, label = paste(t_lines, collapse = "\n"), family = "Helvetica",
             fontface = "bold", size = TITLE_SIZE, lineheight = LINE_H, colour = COL[["ink"]])),
    lapply(seq_along(i_lines), function(i) annotate(
      "text", x = mid, y = y_items[i] / MM_Y, label = paste(i_lines[[i]], collapse = "\n"),
      family = "Helvetica", size = BODY_SIZE, lineheight = LINE_H, colour = COL[["ink"]])))
  all_lines <- c(t_lines, unlist(i_lines))
  widths <- c(measure_text_mm(t_lines, TITLE_PT, face = 2), measure_text_mm(unlist(i_lines), BODY_PT))
  fit_rows[[b$id]] <- data.frame(
    box = b$id, line = all_lines, width_mm = round(widths, 2), usable_width_mm = round(usable, 2),
    fits_width = widths <= usable + 1e-6, body_top_mm = round(top, 2), body_bottom_mm = round(bottom, 2),
    box_inner_bottom_mm = round(b$ymin * MM_Y + PAD_MM, 2), title_bottom_mm = round(title_bottom, 2),
    fits_height = bottom >= b$ymin * MM_Y + PAD_MM & top <= title_bottom - 1.0, row.names = NULL)
}
invisible(grDevices::dev.off())
fit <- bind_rows(fit_rows)
write_csv(fit, file.path(QA_DIR, "Figure_1a_text_fit.csv"))
if (!all(fit$fits_width & fit$fits_height)) stop("[24c] Figure 1a text does not fit; see qa/Figure_1a_text_fit.csv")

arrow_layer <- annotate("segment", x = fig1_arrows$x, xend = fig1_arrows$xend, y = fig1_arrows$y,
                        yend = fig1_arrows$yend, colour = "#66737C", linewidth = 0.55,
                        arrow = arrow(length = unit(1.6, "mm")))
fig1 <- fig1_compose(p1a_frame(c(box_layers, list(arrow_layer))))

# 源数据逐字沿用 24 号脚本 / source data identical to script 24
fig1_source <- tibble(
  item = c("occurrence_records", "recorded_species", "breeding_season_records", "grid_cells", "periods",
           "maximum_rhat_200sp", "minimum_ess_200sp", "available_rhat_500sp", "minimum_ess_500sp"),
  value = c(7489201, 1426, 3149635, 1247, 5, max(conv200$rhat, na.rm = TRUE), min(conv200$ess, na.rm = TRUE),
            sum(!is.na(conv500$rhat)), min(conv500$ess, na.rm = TRUE)),
  note = c("all events", "integrated database", "breeding-season filter", "100-km grid", "2000-2024",
           "4-chain diagnostic table", "4-chain diagnostic table", "R-hat unavailable", "4-chain metadata")
)
write_csv(fig1_source, file.path(SRC_DIR, "SourceData_Figure_1.csv"))
save_main(fig1, "Figure_1_inference_architecture", FIG1_W, FIG1_H)
message("[24c] Figure 1a: ", nrow(fit), " text lines checked, all inside their boxes")
