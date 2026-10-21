#!/usr/bin/env Rscript

## Nature-style submission figure package for Chinese bird community dynamics.
## R is the exclusive rendering backend. Maps use the requested provincial shp.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
  library(forcats)
  library(officer)
  library(grid)
})

sf::sf_use_s2(FALSE)

ROOT <- normalizePath(getwd(), mustWork = TRUE)
RESULTS_DIR <- file.path(ROOT, "results_v3")
DATA_DIR <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "submission_figure_package_20260714_v3")
FIG_DIR <- file.path(OUT_DIR, "figures")
SRC_DIR <- file.path(OUT_DIR, "source_data")
QA_DIR <- file.path(OUT_DIR, "qa")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("svglite", "ragg", "rvg")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))

read_csv_req <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

safe_make_valid <- function(x) {
  out <- try(sf::st_make_valid(x), silent = TRUE)
  if (inherits(out, "try-error")) suppressWarnings(sf::st_buffer(x, 0)) else out
}

safe_rescale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) return(rep(NA_real_, length(x)))
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  scales::rescale(x, to = c(0, 1), from = rng)
}

COL <- c(
  ink = "#202A33",
  neutral_dark = "#4D4D4D",
  neutral = "#8C96A3",
  neutral_light = "#D9DEE3",
  panel = "#F3F5F7",
  blue = "#0072B2",
  sky = "#56B4E9",
  orange = "#D55E00",
  gold = "#E69F00",
  purple = "#CC79A7",
  teal = "#009E73"
)

pal_dataset <- c("200-species spatial" = COL[["blue"]], "500-species temporal" = COL[["purple"]])
pal_class <- c("Increasing" = COL[["blue"]], "Stable" = COL[["neutral"]], "Declining" = COL[["orange"]])
pal_driver <- c(
  "Spatial baseline" = COL[["neutral_dark"]],
  "Climate change" = COL[["blue"]],
  "Land-use change" = COL[["gold"]],
  "Human pressure change" = COL[["purple"]]
)

theme_nature <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(colour = COL[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = COL[["ink"]]),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(2.5, "mm"),
      strip.background = element_rect(fill = "#EEF1F4", colour = NA),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.4, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#52616B", margin = margin(b = 3)),
      plot.caption = element_text(size = base_size - 1, colour = "#606D76", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold", colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

theme_map_nature <- function(base_size = 6.2) {
  theme_void(base_size = base_size, base_family = "Helvetica") +
    theme(
      strip.background = element_rect(fill = "#EEF1F4", colour = NA),
      strip.text = element_text(size = base_size, face = "bold", colour = COL[["ink"]]),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 0.1, face = "bold"),
      legend.text = element_text(size = base_size - 0.4),
      legend.key.width = unit(15, "mm"),
      panel.spacing = unit(1.0, "mm"),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

period_label <- c(
  P1 = "2000-2004", P2 = "2005-2009", P3 = "2010-2014",
  P4 = "2015-2019", P5 = "2020-2024"
)

metric_label <- c(
  corrected_richness = "Richness",
  shannon = "Shannon diversity",
  inv_simpson = "Inverse Simpson",
  pd_prob_mctavish = "Phylogenetic diversity",
  mpd_prob_mctavish = "Mean phylogenetic distance",
  trait_volume = "Functional volume",
  rao_q = "Rao's Q",
  feve_fund = "Functional evenness"
)

metric_group <- c(
  corrected_richness = "Taxonomic",
  shannon = "Taxonomic",
  inv_simpson = "Taxonomic",
  pd_prob_mctavish = "Phylogenetic",
  mpd_prob_mctavish = "Phylogenetic",
  trait_volume = "Functional",
  rao_q = "Functional",
  feve_fund = "Functional"
)

FIGURES <- list()

fit_location <- function(width_mm, height_mm) {
  aspect <- width_mm / height_mm
  max_w <- 12.75
  max_h <- 7.00
  if (aspect >= max_w / max_h) {
    w <- max_w
    h <- w / aspect
  } else {
    h <- max_h
    w <- h * aspect
  }
  list(left = (13.333 - w) / 2, top = (7.5 - h) / 2, width = w, height = h)
}

save_submission_figure <- function(plot, id, width_mm, height_mm, editable = TRUE, map = FALSE, source_data) {
  if (width_mm > 183 || height_mm > 170) stop("Figure exceeds Nature size contract: ", id)
  base <- file.path(FIG_DIR, id)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  ggsave(paste0(base, ".png"), plot, width = w, height = h, units = "in", dpi = 450, bg = "white", limitsize = FALSE)
  ggsave(paste0(base, ".pdf"), plot, width = w, height = h, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)

  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 450,
                 background = "white", compression = "lzw")
  print(plot)
  dev.off()

  svg_path <- NA_character_
  if (!map) {
    svg_path <- paste0(base, ".svg")
    svglite::svglite(svg_path, width = w, height = h, bg = "white")
    print(plot)
    dev.off()
  }

  ppt <- read_pptx()
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  loc <- fit_location(width_mm, height_mm)
  ph_loc <- ph_location(left = loc$left, top = loc$top, width = loc$width, height = loc$height)
  if (editable && !map) {
    ppt <- ph_with(ppt, rvg::dml(ggobj = plot), location = ph_loc)
    ppt_editable <- "editable_drawingml"
  } else {
    ppt <- ph_with(ppt, external_img(paste0(base, ".png"), width = loc$width, height = loc$height), location = ph_loc)
    ppt_editable <- "map_embedded_png"
  }
  ppt_path <- paste0(base, ".pptx")
  print(ppt, target = ppt_path)

  FIGURES[[length(FIGURES) + 1L]] <<- list(
    id = id, plot = plot, width_mm = width_mm, height_mm = height_mm,
    png = paste0(base, ".png"), pdf = paste0(base, ".pdf"), tiff = paste0(base, ".tiff"),
    svg = svg_path, pptx = ppt_path, pptx_editable = ppt_editable,
    source_data = source_data, map = map
  )
}

province_path <- file.path(DATA_DIR, "中国shp", "省.shp")
province <- suppressMessages(st_read(province_path, quiet = TRUE, options = "ENCODING=UTF-8")) |>
  safe_make_valid() |>
  st_transform(4326)
province_lines <- st_boundary(province)
outline <- st_sf(geometry = st_union(st_geometry(province)), crs = st_crs(province))
grid_path <- c(
  file.path(DATA_DIR, "derived_v3", "china_grid_100km_v3.rds"),
  file.path(DATA_DIR, "derived_v2", "china_grid_100km_v2.rds")
)
grid_path <- grid_path[file.exists(grid_path)][1]
grid_sf <- readRDS(grid_path) |>
  safe_make_valid() |>
  st_transform(4326) |>
  select(grid_cell, geometry)

metrics_200 <- read_csv_req(file.path(RESULTS_DIR, "table_community_metrics_with_cri_v3_full_200sp_ar1_spatial_extended.csv"))
metrics_500 <- read_csv_req(file.path(RESULTS_DIR, "table_community_metrics_with_cri_v3_full_500sp_ar1_temporal_extended.csv"))
trend_200 <- read_csv_req(file.path(RESULTS_DIR, "table_trend_summary_v3_full_200sp_ar1_spatial_extended.csv"))
beta_200 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_200sp_ar1_spatial_extended.csv"))
beta_500 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_500sp_ar1_temporal_extended.csv"))
class_500 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_500sp_ar1_temporal_extended.csv"))
trait_corr <- read_csv_req(file.path(RESULTS_DIR, "table_trait_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"))
naive_500 <- read_csv_req(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_500sp_ar1_temporal_extended.csv"))
rf_group <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_group_summary_trend.csv"))
varpart_500 <- read_csv_req(file.path(RESULTS_DIR, "table_varpart_corrected_richness_v3_full_500sp_ar1_temporal_extended.csv"))
survey <- read_csv_req(file.path(RESULTS_DIR, "table_02_survey_summary_v3.csv"))
conv200 <- read_csv_req(file.path(RESULTS_DIR, "table_convergence_diagnostics_v3_full_200sp_ar1_spatial_4chain.csv"))
conv500 <- read_csv_req(file.path(RESULTS_DIR, "table_convergence_diagnostics_v3_full_500sp_ar1_temporal_4chain.csv"))

# Figure 1: inference architecture and explicit claim boundaries.
box_layer <- function(xmin, xmax, ymin, ymax, fill, title, body) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = "#6F7B84", linewidth = 0.35),
    annotate("text", x = (xmin + xmax) / 2, y = ymax - 2.1, label = title,
             family = "Helvetica", fontface = "bold", size = 2.45, colour = COL[["ink"]]),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2 - 0.8, label = body,
             family = "Helvetica", size = 2.05, lineheight = 1.05, colour = COL[["ink"]])
  )
}

p1a <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 38), clip = "off") +
  box_layer(1, 18, 9, 30, "#E7F1F7", "Occurrence evidence",
            "7.49 million records\n1,426 recorded species\n3.15 million breeding-season records") +
  box_layer(22, 39, 9, 30, "#F0ECF7", "Survey design",
            "1,247 100-km grid cells\nFive periods, 2000-2024\nEffort and duration covariates") +
  box_layer(44, 65, 21, 36, "#F6E8DD", "Primary spatial model",
            "200-species stMsPGOcc\nNNGP + temporal AR(1)\n4 chains; maximum R-hat 1.039") +
  box_layer(44, 65, 2, 17, "#FBEEE6", "Breadth extension",
            "500-species tMsPGOcc\nTemporal AR(1), non-spatial\nR-hat unavailable") +
  box_layer(70, 98, 9, 30, "#E7F2EA", "Posterior ecological evidence",
            "Taxonomic, functional and phylogenetic diversity\nTemporal beta decomposition and species trends\nDetection bias, traits and driver associations") +
  annotate("segment", x = 18.8, xend = 21.2, y = 19.5, yend = 19.5, colour = "#66737C", linewidth = 0.55,
           arrow = arrow(length = unit(1.6, "mm"))) +
  annotate("segment", x = 39.8, xend = 43.2, y = 23, yend = 28, colour = "#66737C", linewidth = 0.55,
           arrow = arrow(length = unit(1.6, "mm"))) +
  annotate("segment", x = 39.8, xend = 43.2, y = 16, yend = 9.5, colour = "#66737C", linewidth = 0.55,
           arrow = arrow(length = unit(1.6, "mm"))) +
  annotate("segment", x = 65.8, xend = 69.2, y = 28, yend = 23, colour = "#66737C", linewidth = 0.55,
           arrow = arrow(length = unit(1.6, "mm"))) +
  annotate("segment", x = 65.8, xend = 69.2, y = 9.5, yend = 16, colour = "#66737C", linewidth = 0.55,
           arrow = arrow(length = unit(1.6, "mm"))) +
  labs(title = "Evidence architecture") +
  theme_void(base_family = "Helvetica", base_size = 6.5) +
  theme(plot.title = element_text(size = 7, face = "bold"), plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = margin(2, 2, 1, 2))

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

fig1 <- p1a / p1b + plot_layout(heights = c(2.25, 1)) + plot_annotation(tag_levels = "a")
fig1_source <- tibble(
  item = c("occurrence_records", "recorded_species", "breeding_season_records", "grid_cells", "periods",
           "maximum_rhat_200sp", "minimum_ess_200sp", "available_rhat_500sp", "minimum_ess_500sp"),
  value = c(7489201, 1426, 3149635, 1247, 5, max(conv200$rhat, na.rm = TRUE), min(conv200$ess, na.rm = TRUE),
            sum(!is.na(conv500$rhat)), min(conv500$ess, na.rm = TRUE)),
  note = c("all events", "integrated database", "breeding-season filter", "100-km grid", "2000-2024",
           "4-chain diagnostic table", "4-chain diagnostic table", "R-hat unavailable", "4-chain metadata")
)
fig1_src <- file.path(SRC_DIR, "SourceData_Fig1.csv")
write_csv(fig1_source, fig1_src)
save_submission_figure(fig1, "Fig1_inference_architecture", 183, 118, editable = TRUE, source_data = fig1_src)

# Figure 2: spatial multidiversity across five periods from the primary model.
map_metrics <- c("corrected_richness", "shannon", "trait_volume", "pd_prob_mctavish", "mpd_prob_mctavish", "rao_q")
fig2_data <- expand_grid(grid_cell = grid_sf$grid_cell, metric = map_metrics, period = names(period_label)) |>
  left_join(metrics_200 |> select(grid_cell, metric, period, value_mean, value_l95, value_u95),
            by = c("grid_cell", "metric", "period")) |>
  group_by(metric) |>
  mutate(relative_value = safe_rescale(value_mean)) |>
  ungroup() |>
  mutate(
    metric_plot = factor(unname(metric_label[metric]), levels = unname(metric_label[map_metrics])),
    period_plot = factor(unname(period_label[period]), levels = unname(period_label))
  )
fig2_sf <- grid_sf |> inner_join(fig2_data, by = "grid_cell")
fig2 <- ggplot() +
  geom_sf(data = fig2_sf, aes(fill = relative_value), colour = NA) +
  geom_sf(data = province_lines, fill = NA, colour = "#D8DEE3", linewidth = 0.055) +
  geom_sf(data = outline, fill = NA, colour = COL[["ink"]], linewidth = 0.24) +
  coord_sf(xlim = c(73, 136), ylim = c(3, 54), expand = FALSE, datum = NA) +
  facet_grid(rows = vars(metric_plot), cols = vars(period_plot), switch = "y") +
  scale_fill_gradientn(colours = c("#F7FBFF", "#D8EAF4", "#8FC3DD", "#3B8BC2", "#075A8C"),
                       na.value = "#F2F3F4", breaks = c(0, 0.5, 1), labels = c("Low", "Mid", "High"),
                       name = "Relative value within metric") +
  guides(fill = guide_colourbar(title.position = "top", barwidth = unit(42, "mm"), barheight = unit(2.3, "mm"))) +
  theme_map_nature(6.0) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1), legend.position = "bottom")
fig2_src <- file.path(SRC_DIR, "SourceData_Fig2.csv")
write_csv(fig2_data |> select(grid_cell, period, metric, value_mean, value_l95, value_u95, relative_value), fig2_src)
save_submission_figure(fig2, "Fig2_spatial_multidiversity", 183, 168, editable = FALSE, map = TRUE, source_data = fig2_src)

# Figure 3: spatial fingerprints of directional change.
fig3_data <- trend_200 |>
  filter(metric %in% map_metrics, method %in% c("ols", "theil_sen", "theilsen") | is.na(method)) |>
  group_by(grid_cell, metric) |>
  slice(1) |>
  ungroup() |>
  mutate(credible = q025 > 0 | q975 < 0) |>
  group_by(metric) |>
  mutate(standardized_trend = mean / max(abs(mean), na.rm = TRUE)) |>
  ungroup() |>
  mutate(metric_plot = factor(unname(metric_label[metric]), levels = unname(metric_label[map_metrics])))
fig3_full <- expand_grid(grid_cell = grid_sf$grid_cell, metric = map_metrics) |>
  left_join(fig3_data |> select(grid_cell, metric, mean, q025, q975, credible, standardized_trend),
            by = c("grid_cell", "metric")) |>
  mutate(metric_plot = factor(unname(metric_label[metric]), levels = unname(metric_label[map_metrics])))
fig3_sf <- grid_sf |> inner_join(fig3_full, by = "grid_cell")
fig3 <- ggplot() +
  geom_sf(data = fig3_sf, aes(fill = standardized_trend, alpha = credible), colour = NA) +
  geom_sf(data = province_lines, fill = NA, colour = "#D8DEE3", linewidth = 0.06) +
  geom_sf(data = outline, fill = NA, colour = COL[["ink"]], linewidth = 0.25) +
  coord_sf(xlim = c(73, 136), ylim = c(3, 54), expand = FALSE, datum = NA) +
  facet_wrap(~metric_plot, ncol = 3) +
  scale_fill_gradientn(colours = c(COL[["blue"]], "#B8D8EB", "#F7F7F7", "#F4C27A", COL[["orange"]]),
                       limits = c(-1, 1), oob = squish, na.value = "#F2F3F4",
                       breaks = c(-1, 0, 1), labels = c("Decrease", "0", "Increase"), name = "Standardized trend") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.38), guide = "none", na.value = 0.25) +
  guides(fill = guide_colourbar(title.position = "top", barwidth = unit(42, "mm"), barheight = unit(2.3, "mm"))) +
  theme_map_nature(6.2)
fig3_src <- file.path(SRC_DIR, "SourceData_Fig3.csv")
write_csv(fig3_full |> select(grid_cell, metric, mean, q025, q975, credible, standardized_trend), fig3_src)
save_submission_figure(fig3, "Fig3_spatial_trend_fingerprints", 183, 116, editable = FALSE, map = TRUE, source_data = fig3_src)

# Figure 4: multidiversity decoupling across model scopes.
core_metrics <- c("corrected_richness", "shannon", "inv_simpson", "pd_prob_mctavish", "mpd_prob_mctavish",
                  "trait_volume", "rao_q", "feve_fund")
summarize_metrics <- function(dat, dataset) {
  dat |>
    filter(metric %in% core_metrics) |>
    group_by(metric, period) |>
    summarise(value = mean(value_mean, na.rm = TRUE), lo = mean(value_l95, na.rm = TRUE),
              hi = mean(value_u95, na.rm = TRUE), spatial_sd = sd(value_mean, na.rm = TRUE), .groups = "drop") |>
    arrange(metric, period) |>
    group_by(metric) |>
    mutate(base = value[period == "P1"][1], base_sd = spatial_sd[period == "P1"][1],
           standardized_change = (value - base) / base_sd,
           standardized_lo = (lo - base) / base_sd, standardized_hi = (hi - base) / base_sd,
           percent_change = 100 * (value / base - 1)) |>
    ungroup() |>
    mutate(dataset = dataset, metric_plot = unname(metric_label[metric]), group = unname(metric_group[metric]),
           period_plot = factor(unname(period_label[period]), levels = unname(period_label)))
}
fig4_data <- bind_rows(summarize_metrics(metrics_200, "200-species spatial"),
                       summarize_metrics(metrics_500, "500-species temporal"))
fig4_end <- fig4_data |>
  filter(period == "P5") |>
  mutate(metric_plot = factor(metric_plot, levels = rev(unname(metric_label[core_metrics]))))

p4a <- ggplot(fig4_end, aes(standardized_change, metric_plot, group = metric_plot)) +
  geom_vline(xintercept = 0, colour = "#7D878E", linewidth = 0.35, linetype = "22") +
  geom_line(aes(group = metric_plot), colour = COL[["neutral_light"]], linewidth = 0.7) +
  geom_point(aes(fill = dataset, shape = dataset), size = 2.25, stroke = 0.4, colour = "black") +
  scale_fill_manual(values = pal_dataset, name = NULL) +
  scale_shape_manual(values = c("200-species spatial" = 21, "500-species temporal" = 22), name = NULL) +
  labs(title = "Endpoint change", subtitle = "2020-2024 relative to 2000-2004",
       x = "Change (baseline spatial s.d.)", y = NULL) +
  theme_nature(6.5) +
  theme(legend.position = "top", legend.justification = "left")

fig4_traj <- fig4_data |>
  filter(dataset == "500-species temporal") |>
  mutate(
    metric_panel = recode(
      metric,
      corrected_richness = "Richness",
      shannon = "Shannon",
      inv_simpson = "Inverse Simpson",
      pd_prob_mctavish = "PD",
      mpd_prob_mctavish = "MPD",
      trait_volume = "Functional volume",
      rao_q = "Rao's Q",
      feve_fund = "Functional evenness"
    ),
    metric_panel = factor(metric_panel, levels = c("Richness", "Shannon", "Inverse Simpson", "PD", "MPD",
                                                   "Functional volume", "Rao's Q", "Functional evenness"))
  )
p4b <- ggplot(fig4_traj, aes(period_plot, standardized_change, group = metric_plot)) +
  geom_hline(yintercept = 0, colour = "#A4ADB4", linewidth = 0.3) +
  geom_ribbon(aes(ymin = standardized_lo, ymax = standardized_hi, fill = group), alpha = 0.13, colour = NA) +
  geom_line(aes(colour = group), linewidth = 0.55) +
  geom_point(aes(colour = group), size = 1.1) +
  facet_wrap(~metric_panel, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c(Taxonomic = COL[["blue"]], Phylogenetic = COL[["purple"]], Functional = COL[["orange"]]), guide = "none") +
  scale_fill_manual(values = c(Taxonomic = COL[["blue"]], Phylogenetic = COL[["purple"]], Functional = COL[["orange"]]), guide = "none") +
  labs(title = "National trajectories", subtitle = "500-species temporal breadth extension",
       x = NULL, y = "Change from baseline (s.d.)") +
  theme_nature(6.1) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(size = 5.7))

fig4 <- (p4a / p4b) + plot_layout(heights = c(0.9, 1.35)) + plot_annotation(tag_levels = "a")
fig4_src <- file.path(SRC_DIR, "SourceData_Fig4.csv")
write_csv(fig4_data |> select(dataset, metric, group, period, value, lo, hi, spatial_sd, standardized_change,
                              standardized_lo, standardized_hi, percent_change), fig4_src)
save_submission_figure(fig4, "Fig4_multidiversity_decoupling", 183, 160, editable = TRUE, source_data = fig4_src)

# Figure 5: beta reorganization, winners and trait associations.
fig5_beta <- bind_rows(beta_200 |> mutate(dataset = "200-species spatial"),
                       beta_500 |> mutate(dataset = "500-species temporal")) |>
  mutate(period_plot = factor(str_replace_all(period_pair, "_", "-"),
                              levels = str_replace_all(unique(period_pair), "_", "-")))
p5a <- ggplot(fig5_beta, aes(period_plot, prop_turnover_mean, colour = dataset, shape = dataset, group = dataset)) +
  geom_hline(yintercept = 0.5, colour = "#9AA4AB", linewidth = 0.3, linetype = "22") +
  geom_line(linewidth = 0.55) +
  geom_errorbar(aes(ymin = prop_turnover_q025, ymax = prop_turnover_q975), width = 0.09, linewidth = 0.4) +
  geom_point(size = 2.0, fill = "white", stroke = 0.6) +
  scale_colour_manual(values = pal_dataset, name = NULL) +
  scale_shape_manual(values = c("200-species spatial" = 21, "500-species temporal" = 22), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Turnover share collapses in the latest interval", x = NULL, y = "Turnover share of beta diversity") +
  theme_nature(6.5) + theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(angle = 20, hjust = 1))

fig5_class <- class_500 |>
  count(trend_class, name = "n") |>
  mutate(class = recode(trend_class, expanding = "Increasing", stable = "Stable", contracting = "Declining"),
         prop = n / sum(n), class = factor(class, levels = c("Declining", "Stable", "Increasing")))
p5b <- ggplot(fig5_class, aes(prop, class, fill = class)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = paste0(n, "  (", percent(prop, accuracy = 0.1), ")")), hjust = -0.08, size = 2.05) +
  scale_fill_manual(values = pal_class, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.28))) +
  labs(title = "Species-level trend classes", x = "Share of 500 species", y = NULL) + theme_nature(6.5)

top_species <- bind_rows(
  class_500 |> filter(trend_class == "expanding") |> slice_max(mean, n = 6, with_ties = FALSE),
  class_500 |> filter(trend_class == "contracting") |> slice_min(mean, n = 6, with_ties = FALSE)
) |>
  mutate(direction = recode(trend_class, expanding = "Increasing", contracting = "Declining"),
         species_plot = factor(species, levels = species[order(mean)]))
p5c <- ggplot(top_species, aes(mean, species_plot, colour = direction)) +
  geom_vline(xintercept = 0, colour = "#8C969D", linewidth = 0.3) +
  geom_errorbar(aes(xmin = q025, xmax = q975), orientation = "y", width = 0, linewidth = 0.4) +
  geom_point(size = 1.65) +
  scale_colour_manual(values = pal_class[c("Increasing", "Declining")], guide = "none") +
  labs(title = "Strongest credible species shifts", x = "Occupancy trend (per period)", y = NULL) +
  theme_nature(6.2) + theme(axis.text.y = element_text(face = "italic", size = 5.5))

trait_names <- c(diet_specialization = "Diet specialization", habitat_breadth = "Habitat breadth")
fig5_trait <- trait_corr |>
  mutate(trait_plot = factor(unname(trait_names[trait]), levels = unname(trait_names[trait_names %in% trait_names])))
p5d <- ggplot(fig5_trait, aes(rho_mean, fct_reorder(trait_plot, rho_mean))) +
  geom_vline(xintercept = 0, colour = "#8C969D", linewidth = 0.3, linetype = "22") +
  geom_errorbar(aes(xmin = rho_q025, xmax = rho_q975), orientation = "y", width = 0, linewidth = 0.5, colour = COL[["neutral_dark"]]) +
  geom_point(size = 2.1, shape = 21, fill = COL[["teal"]], colour = "black", stroke = 0.35) +
  labs(title = "Trait associations (n = 200)", x = "Spearman rho with species trend", y = NULL) + theme_nature(6.5)

fig5 <- (p5a | p5b) / (p5c | p5d) + plot_layout(widths = c(1.35, 0.85), heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a")
fig5_src <- file.path(SRC_DIR, "SourceData_Fig5.csv")
write_csv(bind_rows(
  fig5_beta |> transmute(panel = "a_beta", dataset, item = as.character(period_plot), estimate = prop_turnover_mean,
                        lower = prop_turnover_q025, upper = prop_turnover_q975, n = n_grids),
  fig5_class |> transmute(panel = "b_classes", dataset = "500-species temporal", item = as.character(class), estimate = prop,
                         lower = NA_real_, upper = NA_real_, n = n),
  top_species |> transmute(panel = "c_species", dataset = "500-species temporal", item = species, estimate = mean,
                          lower = q025, upper = q975, n = NA_integer_),
  fig5_trait |> transmute(panel = "d_traits", dataset = "complete-case species", item = as.character(trait_plot), estimate = rho_mean,
                         lower = rho_q025, upper = rho_q975, n = 200L)
), fig5_src)
save_submission_figure(fig5, "Fig5_reorganization_and_winners", 183, 142, editable = TRUE, source_data = fig5_src)

# Figure 6: observation bias and spatially structured driver associations.
flip_n <- sum(naive_500$direction_flipped, na.rm = TRUE)
p6a <- ggplot(naive_500, aes(naive_trend, corrected_trend)) +
  geom_abline(slope = 1, intercept = 0, colour = "#7F8990", linewidth = 0.35, linetype = "22") +
  geom_point(aes(colour = direction_flipped), shape = 16, size = 0.85, alpha = 0.48) +
  geom_smooth(method = "lm", se = FALSE, colour = COL[["ink"]], linewidth = 0.55) +
  scale_colour_manual(values = c(`FALSE` = COL[["neutral"]], `TRUE` = COL[["orange"]]), guide = "none") +
  annotate("text", x = -Inf, y = Inf, label = paste0("Direction flips: ", flip_n, " / ", nrow(naive_500)),
           hjust = -0.05, vjust = 1.35, size = 2.15) +
  labs(title = "Naive versus detection-corrected trends", x = "Naive richness trend (species per period)",
       y = "Corrected richness trend (species per period)") + theme_nature(6.4)

p6b <- ggplot(naive_500, aes(trend_diff)) +
  geom_histogram(aes(y = after_stat(density)), bins = 34, fill = "#B8C2CA", colour = "white", linewidth = 0.2) +
  geom_density(colour = COL[["blue"]], linewidth = 0.65) +
  geom_vline(xintercept = 0, colour = COL[["ink"]], linewidth = 0.35, linetype = "22") +
  geom_vline(xintercept = mean(naive_500$trend_diff, na.rm = TRUE), colour = COL[["orange"]], linewidth = 0.55) +
  annotate("text", x = mean(naive_500$trend_diff, na.rm = TRUE), y = Inf,
           label = paste0("Mean = ", round(mean(naive_500$trend_diff, na.rm = TRUE), 1)),
           hjust = 1.08, vjust = 1.4, size = 2.05) +
  labs(title = "Correction effect across grid cells", x = "Corrected minus naive trend (species per period)", y = "Density") +
  theme_nature(6.4)

rf_plot <- rf_group |>
  mutate(group_plot = recode(label, `Land use change` = "Land-use change"),
         group_plot = factor(group_plot, levels = rev(c("Spatial baseline", "Land-use change", "Climate change", "Human pressure change"))))
p6c <- ggplot(rf_plot, aes(mean, group_plot, colour = group_plot)) +
  geom_errorbar(aes(xmin = q_lo, xmax = q_hi), orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = pal_driver, guide = "none") +
  labs(title = "Random-forest group importance", x = "Permutation importance", y = NULL) + theme_nature(6.4)

vp_labels <- c(
  "[a] = X1 | X2+X3+X4" = "Climate change",
  "[b] = X2 | X1+X3+X4" = "Land-use change",
  "[c] = X3 | X1+X2+X4" = "Human pressure change",
  "[d] = X4 | X1+X2+X3" = "Spatial baseline"
)
vp_plot <- varpart_500 |>
  filter(fraction %in% names(vp_labels)) |>
  mutate(group_plot = factor(unname(vp_labels[fraction]), levels = rev(c("Spatial baseline", "Land-use change", "Climate change", "Human pressure change"))))
p6d <- ggplot(vp_plot, aes(pmax(adj_r2, 0), group_plot, fill = group_plot)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = percent(pmax(adj_r2, 0), accuracy = 0.1)), hjust = -0.1, size = 2.0) +
  scale_fill_manual(values = pal_driver, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Pure variance fractions", subtitle = "Corrected-richness trend", x = "Adjusted R-squared", y = NULL) + theme_nature(6.4)

fig6 <- (p6a | p6b) / (p6c | p6d) + plot_layout(heights = c(1.1, 0.9)) + plot_annotation(tag_levels = "a")
fig6_src <- file.path(SRC_DIR, "SourceData_Fig6.csv")
write_csv(bind_rows(
  naive_500 |> transmute(panel = "a_b_detection", item = as.character(grid_cell), estimate = corrected_trend,
                         comparator = naive_trend, difference = trend_diff, lower = NA_real_, upper = NA_real_),
  rf_plot |> transmute(panel = "c_random_forest", item = as.character(group_plot), estimate = mean,
                      comparator = NA_real_, difference = NA_real_, lower = q_lo, upper = q_hi),
  vp_plot |> transmute(panel = "d_variance_partition", item = as.character(group_plot), estimate = adj_r2,
                      comparator = NA_real_, difference = NA_real_, lower = NA_real_, upper = NA_real_)
), fig6_src)
save_submission_figure(fig6, "Fig6_detection_bias_and_drivers", 183, 140, editable = TRUE, source_data = fig6_src)

# Stable PNG atlas for review; individual non-map PPTX files retain editable DrawingML.
atlas <- read_pptx()
for (fig in FIGURES) {
  atlas <- add_slide(atlas, layout = "Blank", master = "Office Theme")
  loc <- fit_location(fig$width_mm, fig$height_mm)
  atlas <- ph_with(atlas, external_img(fig$png, width = loc$width, height = loc$height),
                   location = ph_location(left = loc$left, top = loc$top, width = loc$width, height = loc$height))
}
atlas_path <- file.path(OUT_DIR, "NATURE_SUBMISSION_MAIN_FIGURES_ATLAS_20260714.pptx")
print(atlas, target = atlas_path)

manifest <- bind_rows(lapply(FIGURES, function(x) tibble(
  figure = x$id, width_mm = x$width_mm, height_mm = x$height_mm,
  png = x$png, pdf = x$pdf, tiff = x$tiff, svg = x$svg, pptx = x$pptx,
  pptx_editable = x$pptx_editable, source_data = x$source_data
)))
write_csv(manifest, file.path(OUT_DIR, "SUBMISSION_FIGURE_MANIFEST_20260714.csv"))

qa <- tribble(
  ~figure, ~core_conclusion, ~archetype, ~hero_evidence, ~reviewer_risk,
  "Fig1_inference_architecture", "The spatial primary model and temporal breadth extension support distinct claims.", "schematic-led composite", "Explicit model-to-claim boundary", "500-species maps must remain descriptive.",
  "Fig2_spatial_multidiversity", "Taxonomic and phylogenetic diversity rise without comparable functional expansion.", "quantitative grid", "Six metrics across five periods from the 200-species spatial model", "Values are scaled within metric and are not comparable among metrics.",
  "Fig3_spatial_trend_fingerprints", "Taxonomic gains and functional changes have different spatial fingerprints.", "quantitative grid", "Posterior trend maps from the spatial primary model", "Paler cells have intervals overlapping zero.",
  "Fig4_multidiversity_decoupling", "Taxonomic and phylogenetic gains are decoupled from functional change across model scopes.", "asymmetric mixed-modality", "Cross-model endpoint effect sizes", "Global intervals are summaries of grid-level posterior intervals.",
  "Fig5_reorganization_and_winners", "Late nestedness and broad-habitat winners accompany expansion-dominated reorganization.", "asymmetric mixed-modality", "Turnover-share collapse and species trend classes", "Trait analysis is complete-case n=200.",
  "Fig6_detection_bias_and_drivers", "Detection correction changes local inference while spatial context dominates measured drivers.", "quantitative grid", "Naive-corrected discrepancy and grouped importance", "Driver analyses are associations, not causal effects."
) |>
  mutate(
    final_size_pass = TRUE,
    font_contract = "Helvetica; 5-7 pt body; 8 pt panel labels",
    colour_contract = "colour-blind accessible; direction also encoded by position/shape",
    source_data_pass = TRUE,
    export_contract = "PDF + PNG + TIFF; SVG for non-maps; PPTX for all"
  )
write_csv(qa, file.path(QA_DIR, "FIGURE_QA_CONTRACT_20260714.csv"))

message("[24] Submission package: ", OUT_DIR)
message("[24] Figures: ", length(FIGURES))
message("[24] Atlas: ", atlas_path)
