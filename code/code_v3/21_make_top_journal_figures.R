#!/usr/bin/env Rscript

## Top-journal figure rebuild for China bird community dynamics.
## Outputs: high-resolution PNG, vector PDF, and a PPTX atlas with editable
## DrawingML ggplots when rvg is available. Complex sf maps are embedded as
## high-resolution images in PPTX to keep the deck usable.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forcats)
  library(stringr)
  library(grid)
  library(officer)
})

sf::sf_use_s2(FALSE)

ROOT <- normalizePath(getwd(), mustWork = TRUE)
RESULTS_DIR <- file.path(ROOT, "results_v3")
DATA_DIR <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "figures_top_journal_20260714_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("[21] Output directory: ", OUT_DIR)

has_rvg <- requireNamespace("rvg", quietly = TRUE)
message("[21] Editable PPTX backend rvg available: ", has_rvg)

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  path
}

read_csv_req <- function(path) {
  readr::read_csv(stop_if_missing(path), show_col_types = FALSE)
}

read_csv_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

safe_make_valid <- function(x) {
  out <- try(sf::st_make_valid(x), silent = TRUE)
  if (inherits(out, "try-error")) {
    suppressWarnings(sf::st_buffer(x, 0))
  } else {
    out
  }
}

label_metric <- c(
  corrected_richness = "Occupancy-corrected richness",
  shannon = "Shannon diversity",
  inv_simpson = "Inverse Simpson diversity",
  trait_volume = "Functional volume",
  rao_q = "Rao's Q",
  fdis_prob = "Functional dispersion",
  fric_prob = "Functional richness",
  feve_fund = "Functional evenness",
  pd_prob = "Faith's PD (legacy)",
  mpd_prob = "MPD (legacy)",
  pd_prob_mctavish = "Phylogenetic diversity",
  mpd_prob_mctavish = "Mean phylogenetic distance"
)

label_metric_short <- c(
  corrected_richness = "Richness",
  shannon = "Shannon",
  inv_simpson = "Inv. Simpson",
  trait_volume = "Functional volume",
  rao_q = "Rao's Q",
  fdis_prob = "Functional dispersion",
  fric_prob = "Functional richness",
  feve_fund = "Functional evenness",
  pd_prob = "Faith's PD",
  mpd_prob = "MPD",
  pd_prob_mctavish = "Phylogenetic diversity",
  mpd_prob_mctavish = "MPD"
)

period_label <- c(
  P1 = "2000-2004",
  P2 = "2005-2009",
  P3 = "2010-2014",
  P4 = "2015-2019",
  P5 = "2020-2024"
)

class_label <- c(
  expanding = "Increasing",
  stable = "Stable",
  contracting = "Declining",
  increase = "Increasing",
  no_change = "Stable",
  decrease = "Declining"
)

group_label <- c(
  baseline = "Spatial baseline",
  climate_change = "Climate change",
  landuse_change = "Land-use change",
  human_change = "Human pressure change"
)

term_label <- c(
  delta_t_mean = "Mean temperature change",
  delta_t_extreme = "Thermal-extreme change",
  delta_precip = "Precipitation change",
  delta_forest = "Forest-cover change",
  delta_cropland = "Cropland change",
  delta_crop = "Cropland change",
  delta_grass = "Grassland change",
  delta_urban = "Urban-cover change",
  delta_hfi = "Human footprint change",
  delta_hfp = "Human pressure change",
  centroid_lon = "Longitude",
  centroid_lat = "Latitude",
  elev_mean = "Elevation",
  bio4 = "Temperature seasonality",
  bio7 = "Temperature annual range",
  bio11 = "Mean temperature of coldest quarter",
  bio13 = "Precipitation of wettest month"
)

trait_label <- c(
  habitat_breadth = "Habitat breadth",
  diet_specialization = "Diet specialization",
  migration_score = "Migratory tendency"
)

response_label <- c(
  corrected_richness = "Richness trend",
  shannon = "Shannon trend",
  inv_simpson = "Inverse Simpson trend",
  trait_volume = "Functional volume trend",
  rao_q = "Rao's Q trend",
  fdis_prob = "Functional dispersion trend",
  fric_prob = "Functional richness trend",
  feve_fund = "Functional evenness trend",
  pd_prob_mctavish = "Phylogenetic diversity trend",
  mpd_prob_mctavish = "Mean phylogenetic distance trend"
)

pal_seq <- c("#F7FCF0", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B")
pal_seq_blue <- c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C")
pal_div <- c("#0072B2", "#B8D8EB", "#F7F7F7", "#F4C27A", "#D55E00")
pal_classes <- c(Increasing = "#0072B2", Stable = "#8C96A3", Declining = "#D55E00")
pal_500 <- c("200 spatial species" = "#0072B2", "500 temporal species" = "#CC79A7")
pal_driver <- c(
  "Spatial baseline" = "#4D4D4D",
  "Climate change" = "#0072B2",
  "Land-use change" = "#E69F00",
  "Human pressure change" = "#CC79A7",
  Other = "#7F8790"
)

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(colour = "#1E2A32"),
      plot.title = element_text(face = "bold", size = base_size + 2.8, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size, colour = "#52616B", margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 1.5, colour = "#6B7780", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 0.8, colour = "#35424A"),
      axis.line = element_line(colour = "#26343D", linewidth = 0.35),
      axis.ticks = element_line(colour = "#26343D", linewidth = 0.35),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = base_size, colour = "#23313A"),
      strip.background = element_rect(fill = "#EEF3F7", colour = NA),
      legend.title = element_text(face = "bold", size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.8),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 10)
    )
}

theme_map <- function(base_size = 7.5) {
  theme_void(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2.2, colour = "#1E2A32"),
      plot.subtitle = element_text(size = base_size, colour = "#52616B", margin = margin(b = 7)),
      strip.text = element_text(face = "bold", size = base_size, colour = "#23313A"),
      strip.background = element_rect(fill = "#EEF3F7", colour = NA),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.8),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

wrap_label <- function(x, width = 24) stringr::str_wrap(x, width = width)

nice_metric <- function(x) {
  out <- unname(label_metric[x])
  out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " "))
  out
}

nice_metric_short <- function(x) {
  out <- unname(label_metric_short[x])
  out[is.na(out)] <- nice_metric(x[is.na(out)])
  out
}

nice_term <- function(x) {
  out <- unname(term_label[x])
  out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " "))
  out
}

nice_group <- function(x) {
  out <- unname(group_label[x])
  out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " "))
  out
}

nice_trait <- function(x) {
  out <- unname(trait_label[x])
  out[is.na(out)] <- str_to_sentence(str_replace_all(x[is.na(out)], "_", " "))
  out
}

save_fig <- function(plot, name, width, height, editable = TRUE, map = FALSE) {
  png_path <- file.path(OUT_DIR, paste0(name, ".png"))
  pdf_path <- file.path(OUT_DIR, paste0(name, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 450, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  FIGURES[[length(FIGURES) + 1L]] <<- list(
    name = name,
    title = str_replace_all(name, "_", " "),
    png = png_path,
    pdf = pdf_path,
    plot = plot,
    width = width,
    height = height,
    editable = editable && !map && has_rvg,
    map = map
  )
  invisible(list(png = png_path, pdf = pdf_path))
}

add_slide_with_figure <- function(ppt, fig) {
  ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- officer::ph_with(
    ppt,
    value = str_to_sentence(fig$title),
    location = officer::ph_location(left = 0.35, top = 0.14, width = 12.6, height = 0.34)
  )
  loc <- officer::ph_location(left = 0.3, top = 0.55, width = 12.75, height = 6.75)
  if (isTRUE(fig$editable)) {
    ppt <- officer::ph_with(ppt, value = rvg::dml(ggobj = fig$plot), location = loc)
  } else {
    ppt <- officer::ph_with(
      ppt,
      value = officer::external_img(fig$png, width = 12.75, height = 6.75),
      location = loc
    )
  }
  ppt
}

write_pptx_atlas <- function(figures) {
  ppt_path <- file.path(OUT_DIR, "bird_community_top_journal_figure_atlas_20260714_v2.pptx")
  ppt <- officer::read_pptx()
  ppt <- officer::layout_summary(ppt) |> invisible()
  ppt <- officer::read_pptx()
  for (fig in figures) {
    ppt <- add_slide_with_figure(ppt, fig)
  }
  print(ppt, target = ppt_path)
  ppt_path
}

write_individual_pptx <- function(figures) {
  paths <- character(length(figures))
  for (i in seq_along(figures)) {
    fig <- figures[[i]]
    ppt <- officer::read_pptx()
    ppt <- add_slide_with_figure(ppt, fig)
    ppt_path <- file.path(OUT_DIR, paste0(fig$name, ".pptx"))
    print(ppt, target = ppt_path)
    paths[i] <- ppt_path
  }
  paths
}

load_china_layers <- function() {
  poly_path <- file.path(DATA_DIR, "中国shp", "省.shp")

  province <- suppressMessages(sf::st_read(stop_if_missing(poly_path), quiet = TRUE, options = "ENCODING=UTF-8")) |>
    safe_make_valid() |>
    sf::st_transform(4326)
  province_lines <- sf::st_boundary(province)

  outline <- sf::st_sf(geometry = sf::st_union(sf::st_geometry(province)), crs = sf::st_crs(province))

  list(province = province, province_lines = province_lines, outline = outline)
}

load_grid <- function() {
  candidates <- c(
    file.path(DATA_DIR, "derived_v3", "china_grid_100km_v3.rds"),
    file.path(DATA_DIR, "derived_v2", "china_grid_100km_v2.rds")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) stop("No 100-km grid RDS found.", call. = FALSE)
  readRDS(path) |>
    safe_make_valid() |>
    sf::st_transform(4326) |>
    dplyr::select(grid_cell, geometry)
}

china_layers <- load_china_layers()
grid_sf <- load_grid()
FIGURES <- list()

metrics_200 <- read_csv_req(file.path(RESULTS_DIR, "table_community_metrics_with_cri_v3_full_200sp_ar1_spatial_extended.csv"))
metrics_500 <- read_csv_req(file.path(RESULTS_DIR, "table_community_metrics_with_cri_v3_full_500sp_ar1_temporal_extended.csv"))
trend_200 <- read_csv_req(file.path(RESULTS_DIR, "table_trend_summary_v3_full_200sp_ar1_spatial_extended.csv"))
trend_500 <- read_csv_req(file.path(RESULTS_DIR, "table_trend_summary_v3_full_500sp_ar1_temporal_extended.csv"))
beta_200 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_summary_v3_full_200sp_ar1_spatial_extended.csv"))
beta_global_200 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_200sp_ar1_spatial_extended.csv"))
beta_global_500 <- read_csv_req(file.path(RESULTS_DIR, "table_baselga_global_v3_full_500sp_ar1_temporal_extended.csv"))
class_200 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_200sp_ar1_spatial_extended.csv"))
class_500 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_500sp_ar1_temporal_extended.csv"))
species_500 <- read_csv_req(file.path(RESULTS_DIR, "table_species_trend_v3_full_500sp_ar1_temporal_extended.csv"))
naive_200 <- read_csv_req(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_200sp_ar1_spatial_extended.csv"))
naive_500 <- read_csv_req(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_500sp_ar1_temporal_extended.csv"))
rf_var <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_variable_summary_trend.csv"))
rf_group <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_group_summary_trend.csv"))
trait_500 <- read_csv_req(file.path(RESULTS_DIR, "table_trait_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"))
env_500 <- read_csv_req(file.path(RESULTS_DIR, "table_env_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"))

available_main_metrics <- intersect(
  c("corrected_richness", "shannon", "trait_volume", "pd_prob_mctavish", "mpd_prob_mctavish", "rao_q"),
  unique(metrics_200$metric)
)
if (length(available_main_metrics) < 4) {
  available_main_metrics <- intersect(
    c("corrected_richness", "shannon", "trait_volume", "pd_prob", "mpd_prob", "rao_q"),
    unique(metrics_200$metric)
  )
}
periods_main <- intersect(names(period_label), unique(metrics_200$period))

make_timeslice_map <- function(metric_df, dataset_label, metrics_use, periods_use) {
  d <- metric_df |>
    filter(metric %in% metrics_use, period %in% periods_use) |>
    select(grid_cell, period, metric, value_mean)

  combos <- tidyr::expand_grid(
    grid_cell = grid_sf$grid_cell,
    metric = metrics_use,
    period = periods_use
  )

  d_full <- combos |>
    left_join(d, by = c("grid_cell", "metric", "period")) |>
    group_by(metric) |>
    mutate(value_scaled = scales::rescale(value_mean, to = c(0, 1), from = range(value_mean, na.rm = TRUE))) |>
    ungroup() |>
    mutate(
      metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(metrics_use)),
      period_label = factor(unname(period_label[period]), levels = unname(period_label[periods_use]))
    )

  d_sf <- grid_sf |>
    inner_join(d_full, by = "grid_cell")

  ggplot() +
    geom_sf(data = d_sf, aes(fill = value_scaled), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#D7DEE5", linewidth = 0.08) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    coord_sf(xlim = c(73, 136), ylim = c(3, 54), expand = FALSE, datum = NA) +
    facet_grid(metric_label ~ period_label) +
    scale_fill_gradientn(
      colours = pal_seq,
      na.value = "#F2F4F7",
      name = "Relative value\nwithin metric",
      breaks = c(0, 0.5, 1),
      labels = c("Low", "Mid", "High")
    ) +
    labs(
      title = paste0("Multidiversity geography across five survey periods: ", dataset_label),
      subtitle = "Occupancy-corrected diversity across five survey periods",
      caption = "Values are scaled within each diversity dimension to support cross-period spatial comparison."
    ) +
    theme_map(7.2)
}

p1 <- make_timeslice_map(
  metrics_200,
  "200-species spatial model",
  available_main_metrics,
  periods_main
)
save_fig(p1, "fig01_multidiversity_timeslices_200sp_top_journal", 13.2, 9.2, editable = FALSE, map = TRUE)

make_trend_map <- function(trend_df, dataset_label, metrics_use) {
  d <- trend_df |>
    filter(metric %in% metrics_use, method %in% c("ols", "theil_sen", "theilsen") | is.na(method)) |>
    group_by(grid_cell, metric) |>
    slice(1) |>
    ungroup() |>
    mutate(sig = q025 > 0 | q975 < 0) |>
    group_by(metric) |>
    mutate(trend_scaled = mean / max(abs(mean), na.rm = TRUE)) |>
    ungroup() |>
    mutate(metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(metrics_use)))

  combos <- tidyr::expand_grid(grid_cell = grid_sf$grid_cell, metric = metrics_use) |>
    mutate(metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(metrics_use)))

  d_full <- combos |>
    left_join(select(d, grid_cell, metric, mean, q025, q975, sig, trend_scaled), by = c("grid_cell", "metric"))

  d_sf <- grid_sf |>
    inner_join(d_full, by = "grid_cell")

  ggplot() +
    geom_sf(data = d_sf, aes(fill = trend_scaled, alpha = sig), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#D7DEE5", linewidth = 0.08) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    coord_sf(xlim = c(73, 136), ylim = c(3, 54), expand = FALSE, datum = NA) +
    facet_wrap(~ metric_label, ncol = 3) +
    scale_fill_gradientn(
      colours = pal_div,
      limits = c(-1, 1),
      oob = scales::squish,
      na.value = "#F2F4F7",
      name = "Standardized\ntrend"
    ) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.52), guide = "none", na.value = 0.35) +
    labs(
      title = paste0("Spatial fingerprints of community change: ", dataset_label),
      subtitle = "Blue indicates relative decreases; red indicates relative increases; paler cells have 95% CrI overlapping zero"
    ) +
    theme_map(7.6)
}

p2 <- make_trend_map(
  trend_200,
  "200-species spatial model",
  intersect(available_main_metrics, unique(trend_200$metric))
)
save_fig(p2, "fig02_community_trends_200sp_top_journal", 12.5, 7.9, editable = FALSE, map = TRUE)

trend_metrics_500 <- intersect(
  c("corrected_richness", "shannon", "trait_volume", "pd_prob_mctavish", "mpd_prob_mctavish", "rao_q"),
  unique(trend_500$metric)
)
if (length(trend_metrics_500) < 4) {
  trend_metrics_500 <- intersect(
    c("corrected_richness", "shannon", "trait_volume", "pd_prob", "mpd_prob", "rao_q"),
    unique(trend_500$metric)
  )
}

p2b <- make_trend_map(
  trend_500,
  "500-species temporal model",
  trend_metrics_500
)
save_fig(p2b, "fig02b_community_trends_500sp_temporal_top_journal", 12.5, 7.9, editable = FALSE, map = TRUE)

make_beta_map <- function(beta_df) {
  metrics_use <- c("beta_sim", "beta_sne", "beta_sor")
  metric_names <- c(beta_sim = "Turnover", beta_sne = "Nestedness", beta_sor = "Total beta diversity")
  period_pairs <- unique(beta_df$period_pair)

  d <- beta_df |>
    filter(metric %in% metrics_use) |>
    mutate(
      metric_label = factor(unname(metric_names[metric]), levels = unname(metric_names[metrics_use])),
      period_label = factor(str_replace_all(period_pair, "_", " to "), levels = str_replace_all(period_pairs, "_", " to "))
    )

  combos <- tidyr::expand_grid(
    grid_cell = grid_sf$grid_cell,
    metric = metrics_use,
    period_pair = period_pairs
  ) |>
    mutate(
      metric_label = factor(unname(metric_names[metric]), levels = unname(metric_names[metrics_use])),
      period_label = factor(str_replace_all(period_pair, "_", " to "), levels = str_replace_all(period_pairs, "_", " to "))
    )

  d_full <- combos |>
    left_join(select(d, grid_cell, metric, period_pair, mean), by = c("grid_cell", "metric", "period_pair"))

  d_sf <- grid_sf |>
    inner_join(d_full, by = "grid_cell")

  ggplot() +
    geom_sf(data = d_sf, aes(fill = mean), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#D7DEE5", linewidth = 0.08) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    coord_sf(xlim = c(73, 136), ylim = c(3, 54), expand = FALSE, datum = NA) +
    facet_grid(metric_label ~ period_label) +
    scale_fill_gradientn(colours = pal_seq_blue, na.value = "#F2F4F7", name = "Beta diversity") +
    labs(
      title = "Temporal beta diversity decomposition across China",
      subtitle = "Baselga turnover and nestedness components from occupancy-corrected communities"
    ) +
    theme_map(7.1)
}

p3 <- make_beta_map(beta_200)
save_fig(p3, "fig03_temporal_beta_baselga_200sp_top_journal", 13.2, 8.5, editable = FALSE, map = TRUE)

make_global_trajectory <- function(m200, m500) {
  summarize_global <- function(x, label) {
    x |>
      filter(metric %in% c("corrected_richness", "shannon", "trait_volume", "rao_q", "pd_prob_mctavish", "mpd_prob_mctavish")) |>
      group_by(period, metric) |>
      summarise(
        value = mean(value_mean, na.rm = TRUE),
        lo = mean(value_l95, na.rm = TRUE),
        hi = mean(value_u95, na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(metric, period) |>
      group_by(metric) |>
      mutate(
        base = value[match("P1", period)],
        index = 100 * value / base,
        index_lo = 100 * lo / base,
        index_hi = 100 * hi / base
      ) |>
      ungroup() |>
      mutate(
        dataset = label,
        metric_label = factor(nice_metric(metric), levels = nice_metric(unique(metric))),
        period_label = factor(unname(period_label[period]), levels = unname(period_label[names(period_label) %in% unique(period)]))
      )
  }

  d <- bind_rows(
    summarize_global(m200, "200 spatial species"),
    summarize_global(m500, "500 temporal species")
  )

  ggplot(d, aes(period_label, index, colour = dataset, group = dataset, fill = dataset)) +
    geom_hline(yintercept = 100, linewidth = 0.35, colour = "#9EA8B1") +
    geom_ribbon(aes(ymin = index_lo, ymax = index_hi), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8, stroke = 0.25) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = pal_500, name = NULL) +
    scale_fill_manual(values = pal_500, name = NULL) +
    labs(
      title = "National multidiversity trajectories",
      subtitle = "Each dimension is indexed to 100 at 2000-2004 to compare rates of change across metrics",
      x = NULL,
      y = "Index (2000-2004 = 100)"
    ) +
    theme_pub(8.5) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

p4 <- make_global_trajectory(metrics_200, metrics_500)
save_fig(p4, "fig04_global_multidiversity_trajectories_top_journal_editable", 11.8, 7.2, editable = TRUE)

make_species_class_plot <- function(c200, c500, sp500) {
  class_counts <- bind_rows(
    c200 |> mutate(dataset = "200 spatial species"),
    c500 |> mutate(dataset = "500 temporal species")
  ) |>
    mutate(class = unname(class_label[trend_class])) |>
    mutate(class = if_else(is.na(class), str_to_sentence(trend_class), class)) |>
    count(dataset, class, name = "n") |>
    group_by(dataset) |>
    mutate(prop = n / sum(n), label = paste0(n, "\n", percent(prop, accuracy = 0.1))) |>
    ungroup() |>
    mutate(class = factor(class, levels = c("Declining", "Stable", "Increasing")))

  top_sp <- sp500 |>
    filter(method %in% c("theil_sen", "theilsen") | is.na(method)) |>
    mutate(rank_abs = rank(-abs(mean), ties.method = "first")) |>
    arrange(desc(mean)) |>
    slice_head(n = 10) |>
    bind_rows(sp500 |> arrange(mean) |> slice_head(n = 10)) |>
    distinct(species, .keep_all = TRUE) |>
    mutate(
      direction = if_else(mean >= 0, "Increasing", "Declining"),
      species = factor(species, levels = species[order(mean)])
    )

  p_class <- ggplot(class_counts, aes(dataset, prop, fill = class)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.5, lineheight = 0.92, colour = "white") +
    scale_fill_manual(values = pal_classes, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Species share", title = "Trend classes") +
    theme_pub(8.2) +
    theme(axis.text.x = element_text(face = "bold"), legend.position = "bottom")

  p_species <- ggplot(top_sp, aes(mean, species, colour = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#AAB4BC") +
    geom_errorbar(aes(xmin = q025, xmax = q975), orientation = "y", width = 0, linewidth = 0.45, alpha = 0.85) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = pal_classes[c("Increasing", "Declining")], name = NULL) +
    labs(
      x = "Species occupancy trend",
      y = NULL,
      title = "Most shifted species in the 500-species extension"
    ) +
    theme_pub(8.2) +
    theme(axis.text.y = element_text(face = "italic", size = 7.3))

  p_class + p_species + plot_layout(widths = c(0.9, 1.45)) +
    plot_annotation(
      title = "Species-level winners and losers",
      subtitle = "Trend classification and the strongest 500-species temporal shifts",
      theme = theme_pub(9)
    )
}

p5 <- make_species_class_plot(class_200, class_500, species_500)
save_fig(p5, "fig05_species_trend_classes_500sp_top_journal_editable", 12.2, 6.6, editable = TRUE)

make_naive_corrected_plot <- function(n200, n500) {
  d <- bind_rows(
    n200 |> mutate(dataset = "200 spatial species"),
    n500 |> mutate(dataset = "500 temporal species")
  ) |>
    mutate(direction_flipped = if_else(is.na(direction_flipped), FALSE, direction_flipped))

  flip_lab <- d |>
    group_by(dataset) |>
    summarise(
      x = quantile(naive_trend, 0.03, na.rm = TRUE),
      y = quantile(corrected_trend, 0.94, na.rm = TRUE),
      lab = paste0("Direction flips: ", sum(direction_flipped, na.rm = TRUE), " / ", n()),
      .groups = "drop"
    )

  ggplot(d, aes(naive_trend, corrected_trend)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.45, colour = "#8F9BA4", linetype = "22") +
    geom_point(aes(colour = direction_flipped), size = 1.15, alpha = 0.58) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, colour = "#23313A") +
    geom_text(data = flip_lab, aes(x = x, y = y, label = lab), inherit.aes = FALSE, hjust = 0, size = 2.65, colour = "#23313A") +
    facet_wrap(~ dataset, scales = "free", ncol = 2) +
    scale_colour_manual(values = c(`FALSE` = "#61707A", `TRUE` = "#D55E00"), labels = c("No flip", "Direction flip"), name = NULL) +
    labs(
      title = "Observation bias changes ecological conclusions",
      subtitle = "Comparison of naive detection trends and occupancy-corrected richness trends across 100-km grid cells",
      x = "Naive observed richness trend",
      y = "Occupancy-corrected richness trend"
    ) +
    theme_pub(8.6)
}

p6 <- make_naive_corrected_plot(naive_200, naive_500)
save_fig(p6, "fig06_naive_vs_corrected_trends_top_journal_editable", 10.6, 5.7, editable = TRUE)

make_beta_global_plot <- function(g200, g500) {
  d <- bind_rows(
    g200 |> mutate(dataset = "200 spatial species"),
    g500 |> mutate(dataset = "500 temporal species")
  ) |>
    mutate(
      period_pair = factor(str_replace_all(period_pair, "_", " to "), levels = unique(str_replace_all(period_pair, "_", " to "))),
      nested_mean = 1 - prop_turnover_mean,
      nested_q025 = 1 - prop_turnover_q975,
      nested_q975 = 1 - prop_turnover_q025
    ) |>
    select(dataset, period_pair, prop_turnover_mean, prop_turnover_q025, prop_turnover_q975, nested_mean, nested_q025, nested_q975) |>
    pivot_longer(
      cols = c(prop_turnover_mean, nested_mean),
      names_to = "component",
      values_to = "mean"
    ) |>
    mutate(
      q025 = if_else(component == "prop_turnover_mean", prop_turnover_q025, nested_q025),
      q975 = if_else(component == "prop_turnover_mean", prop_turnover_q975, nested_q975),
      component = recode(component, prop_turnover_mean = "Turnover", nested_mean = "Nestedness")
    )

  ggplot(d, aes(period_pair, mean, fill = component)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = "white", linewidth = 0.22) +
    geom_errorbar(aes(ymin = q025, ymax = q975), position = position_dodge(width = 0.7), width = 0.17, linewidth = 0.35) +
    facet_wrap(~ dataset, ncol = 1) +
    scale_fill_manual(values = c(Turnover = "#0072B2", Nestedness = "#E69F00"), name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "Homogenization is governed by turnover collapse",
      subtitle = "Global Baselga decomposition across consecutive five-year periods",
      x = NULL,
      y = "Component share of total beta diversity"
    ) +
    theme_pub(8.6) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

p7 <- make_beta_global_plot(beta_global_200, beta_global_500)
save_fig(p7, "fig07_global_baselga_turnover_nestedness_top_journal_editable", 9.8, 6.6, editable = TRUE)

make_rf_plot <- function(var_tbl, group_tbl) {
  var_d <- var_tbl |>
    mutate(
      group_label = nice_group(group),
      label_plot = if_else(!is.na(label), label, nice_term(variable)),
      label_plot = wrap_label(label_plot, 28)
    ) |>
    arrange(desc(mean)) |>
    slice_head(n = 18) |>
    mutate(label_plot = factor(label_plot, levels = rev(label_plot)))

  group_d <- group_tbl
  group_d$group_label <- if ("label" %in% names(group_d)) group_d$label else nice_group(group_d$group)
  group_d <- group_d |>
    mutate(group_label = wrap_label(group_label, 24)) |>
    arrange(mean) |>
    mutate(group_label = factor(group_label, levels = group_label))

  p_var <- ggplot(var_d, aes(mean, label_plot, colour = group_label)) +
    geom_errorbar(aes(xmin = q_lo, xmax = q_hi), orientation = "y", width = 0, linewidth = 0.45, alpha = 0.88) +
    geom_point(size = 2.1) +
    scale_colour_manual(values = pal_driver, name = NULL, na.value = "#7F8790") +
    labs(
      title = "Top predictors",
      x = "Permutation importance",
      y = NULL
    ) +
    theme_pub(8.1) +
    theme(axis.text.y = element_text(size = 7.1))

  p_group <- ggplot(group_d, aes(mean, group_label, fill = group_label)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
    geom_errorbar(aes(xmin = q_lo, xmax = q_hi), orientation = "y", width = 0.18, linewidth = 0.38, colour = "#23313A") +
    scale_fill_manual(values = pal_driver, guide = "none", na.value = "#7F8790") +
    labs(
      title = "Predictor groups",
      x = "Grouped importance",
      y = NULL
    ) +
    theme_pub(8.1)

  p_var + p_group + plot_layout(widths = c(1.45, 0.85)) +
    plot_annotation(
      title = "Environmental drivers of community change",
      subtitle = "Random-forest importance from posterior draws of diversity trends",
      theme = theme_pub(9)
    )
}

p8 <- make_rf_plot(rf_var, rf_group)
save_fig(p8, "fig08_random_forest_driver_importance_top_journal_editable", 12.1, 6.8, editable = TRUE)

make_trait_env_plot <- function(trait_tbl, env_tbl) {
  trait_d <- trait_tbl |>
    mutate(
      label = factor(nice_trait(trait), levels = rev(nice_trait(trait))),
      significant = rho_q025 > 0 | rho_q975 < 0
    )

  env_d <- env_tbl |>
    mutate(
      label = factor(nice_term(env_var), levels = nice_term(env_var)[order(spearman_rho)]),
      n_lab = paste0("n = ", n_species)
    )

  p_trait <- ggplot(trait_d, aes(rho_mean, label, colour = significant)) +
    geom_vline(xintercept = 0, colour = "#9EA8B1", linewidth = 0.35, linetype = "22") +
    geom_errorbar(aes(xmin = rho_q025, xmax = rho_q975), orientation = "y", width = 0, linewidth = 0.55) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#61707A"), guide = "none") +
    labs(title = "Trait associations", x = "Spearman rho with species trend", y = NULL) +
    theme_pub(8.4)

  p_env <- ggplot(env_d, aes(spearman_rho, label)) +
    geom_vline(xintercept = 0, colour = "#9EA8B1", linewidth = 0.35, linetype = "22") +
    geom_col(width = 0.58, fill = "#0072B2", alpha = 0.88) +
    geom_text(aes(label = n_lab), hjust = -0.1, size = 2.5, colour = "#52616B") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.2))) +
    labs(title = "Climate-gradient associations", x = "Spearman rho with species trend", y = NULL) +
    theme_pub(8.4)

  p_trait + p_env + plot_layout(widths = c(0.95, 1.2)) +
    plot_annotation(
      title = "Ecological mechanisms in the 500-species extension",
      subtitle = "Complete-case trait and environmental screens; interpret as association, not causal effect",
      theme = theme_pub(9)
    )
}

p9 <- make_trait_env_plot(trait_500, env_500)
save_fig(p9, "fig09_trait_environment_mechanisms_500sp_top_journal_editable", 11.2, 5.8, editable = TRUE)

make_coefficients_plot <- function() {
  files <- list.files(
    RESULTS_DIR,
    pattern = "^table_brms_driver_trend_change_coefficients_.*_v3_full_500sp_ar1_temporal_extended\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) return(NULL)

  d <- lapply(files, function(f) {
    response <- basename(f) |>
      str_remove("^table_brms_driver_trend_change_coefficients_") |>
      str_remove("_v3_full_500sp_ar1_temporal_extended\\.csv$")
    read_csv(f, show_col_types = FALSE) |>
      mutate(response = response)
  }) |>
    bind_rows() |>
    filter(term != "Intercept") |>
    mutate(
      term_label = factor(nice_term(term), levels = rev(unique(nice_term(term)))),
      response_label = factor(
        unname(response_label[response]),
        levels = unname(response_label[unique(response)])
      ),
      credible = Q2.5 > 0 | Q97.5 < 0
    )

  ggplot(d, aes(Estimate, term_label, colour = credible)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9EA8B1", linetype = "22") +
    geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), orientation = "y", width = 0, linewidth = 0.43) +
    geom_point(size = 1.7) +
    facet_wrap(~ response_label, scales = "free_x", ncol = 3) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#61707A"), labels = c("CrI overlaps 0", "CrI excludes 0"), name = NULL) +
    labs(
      title = "Posterior driver coefficients across 500-species diversity responses",
      subtitle = "Temporal-change covariates fitted to occupancy-corrected trend estimates",
      x = "Posterior coefficient",
      y = NULL
    ) +
    theme_pub(8.1) +
    theme(axis.text.y = element_text(size = 7.2))
}

p10 <- make_coefficients_plot()
if (!is.null(p10)) {
  save_fig(p10, "fig10_driver_coefficients_500sp_top_journal_editable", 12.2, 7.4, editable = TRUE)
}

individual_pptx_paths <- write_individual_pptx(FIGURES)
pptx_path <- write_pptx_atlas(FIGURES)

manifest <- tibble(
  figure = vapply(FIGURES, `[[`, character(1), "name"),
  png = vapply(FIGURES, `[[`, character(1), "png"),
  pdf = vapply(FIGURES, `[[`, character(1), "pdf"),
  pptx = individual_pptx_paths,
  pptx_editable = vapply(FIGURES, function(x) if (isTRUE(x$editable)) "editable_drawingml" else if (isTRUE(x$map)) "map_embedded_png" else "embedded_png", character(1))
)
readr::write_csv(manifest, file.path(OUT_DIR, "figure_manifest_20260714_v2.csv"))

message("[21] Wrote ", length(FIGURES), " figures.")
message("[21] PPTX atlas: ", pptx_path)
message("[21] Manifest: ", file.path(OUT_DIR, "figure_manifest_20260714_v2.csv"))
