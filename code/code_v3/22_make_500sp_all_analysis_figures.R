#!/usr/bin/env Rscript

## 500-species all-analysis figure suite.
## Basemap is forced to data/中国shp/省.shp, as requested.
## Outputs PNG, PDF, individual PPTX, and one PPTX atlas.

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

INPUT_ROOT <- normalizePath(Sys.getenv("BIRD_PROJECT_ROOT", getwd()), mustWork = TRUE)
OUTPUT_ROOT <- normalizePath(Sys.getenv("BIRD_OUTPUT_ROOT", INPUT_ROOT), mustWork = FALSE)
RUN_500 <- Sys.getenv("V3_RUN_LABEL", "v3_full_500sp_ar1_temporal")
RUN_500_EXT <- paste0(RUN_500, "_extended")
RESULTS_DIR <- file.path(OUTPUT_ROOT, "results_v3")
DATA_DIR <- file.path(INPUT_ROOT, "data")
OUT_DIR <- file.path(OUTPUT_ROOT, paste0("figures_500sp_all_analysis_", RUN_500))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("[22] Output directory: ", OUT_DIR)
has_rvg <- requireNamespace("rvg", quietly = TRUE)
message("[22] Editable PPTX backend rvg available: ", has_rvg)

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  path
}

read_csv_req <- function(path) readr::read_csv(stop_if_missing(path), show_col_types = FALSE)
read_csv_optional <- function(path) if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else NULL

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

safe_standardize <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  denom <- max(abs(x), na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) return(rep(NA_real_, length(x)))
  x / denom
}

metric_label <- c(
  corrected_richness = "Occupancy-corrected richness",
  shannon = "Shannon diversity",
  inv_simpson = "Inverse Simpson diversity",
  pd_prob_mctavish = "Phylogenetic diversity",
  mpd_prob_mctavish = "Mean phylogenetic distance",
  pd_prob = "Legacy Faith's PD",
  mpd_prob = "Legacy MPD",
  fmpd_prob = "Functional MPD",
  trait_volume = "Functional volume",
  trait_dispersion = "Trait dispersion",
  fdis_prob = "Functional dispersion",
  fric_prob = "Functional richness",
  fdiv = "Functional divergence",
  fdiv_fund = "Functional divergence (fund.)",
  feve = "Functional evenness",
  feve_fund = "Functional evenness (fund.)",
  rao_q = "Rao's Q",
  raoq_fund = "Rao's Q (fund.)",
  cwm_pc1 = "CWM PC1",
  cwm_pc2 = "CWM PC2"
)

metric_short <- c(
  corrected_richness = "Richness",
  shannon = "Shannon",
  inv_simpson = "Inv. Simpson",
  pd_prob_mctavish = "PD",
  mpd_prob_mctavish = "MPD",
  pd_prob = "Legacy PD",
  mpd_prob = "Legacy MPD",
  fmpd_prob = "Functional MPD",
  trait_volume = "Trait volume",
  trait_dispersion = "Trait dispersion",
  fdis_prob = "FDis",
  fric_prob = "FRic",
  fdiv = "FDiv",
  fdiv_fund = "FDiv fund.",
  feve = "FEve",
  feve_fund = "FEve fund.",
  rao_q = "Rao's Q",
  raoq_fund = "RaoQ fund.",
  cwm_pc1 = "CWM PC1",
  cwm_pc2 = "CWM PC2"
)

metric_order <- c(
  "corrected_richness", "shannon", "inv_simpson", "pd_prob_mctavish", "mpd_prob_mctavish",
  "pd_prob", "mpd_prob", "fmpd_prob", "trait_volume", "trait_dispersion",
  "fdis_prob", "fric_prob", "fdiv", "fdiv_fund", "feve",
  "feve_fund", "rao_q", "raoq_fund", "cwm_pc1", "cwm_pc2"
)

period_label <- c(
  P1 = "2000-2004",
  P2 = "2005-2009",
  P3 = "2010-2014",
  P4 = "2015-2019",
  P5 = "2020-2024"
)

class_label <- c(expanding = "Increasing", stable = "Stable", contracting = "Declining")

term_label <- c(
  delta_t_mean = "Mean temperature change",
  delta_t_extreme = "Thermal-extreme change",
  delta_t_std = "Standardized warming",
  delta_forest = "Forest-cover change",
  delta_impervious = "Impervious-surface change",
  delta_natural = "Natural land-cover change",
  delta_hfi = "Human footprint change",
  bio11 = "Mean temperature (baseline)",
  elev_mean = "Elevation",
  centroid_lon = "Longitude",
  centroid_lat = "Latitude",
  bio4 = "Temperature seasonality",
  bio7 = "Temperature annual range",
  bio13 = "Precipitation of wettest month",
  landcover_built = "Built-up land",
  landcover_cropland = "Cropland",
  hfi_mean = "Human footprint"
)

group_label <- c(
  climate_change = "Climate change",
  landuse_change = "Land-use change",
  human_change = "Human pressure change",
  baseline = "Spatial baseline"
)

varpart_label <- c(
  "[a] = X1 | X2+X3+X4" = "Climate change",
  "[b] = X2 | X1+X3+X4" = "Land-use change",
  "[c] = X3 | X1+X2+X4" = "Human pressure change",
  "[d] = X4 | X1+X2+X3" = "Spatial baseline"
)

pal_seq <- c("#F7FCF0", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B")
pal_blue <- c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C")
pal_div <- c("#0072B2", "#B8D8EB", "#F7F7F7", "#F4C27A", "#D55E00")
pal_class <- c(Increasing = "#0072B2", Stable = "#8C96A3", Declining = "#D55E00")
pal_driver <- c(
  "Climate change" = "#0072B2",
  "Land-use change" = "#E69F00",
  "Human pressure change" = "#CC79A7",
  "Spatial baseline" = "#4D4D4D"
)

nice_metric <- function(x) {
  out <- unname(metric_label[x])
  out[is.na(out)] <- stringr::str_to_sentence(stringr::str_replace_all(x[is.na(out)], "_", " "))
  out
}

nice_metric_short <- function(x) {
  out <- unname(metric_short[x])
  out[is.na(out)] <- nice_metric(x[is.na(out)])
  out
}

nice_term <- function(x) {
  out <- unname(term_label[x])
  out[is.na(out)] <- stringr::str_to_sentence(stringr::str_replace_all(x[is.na(out)], "_", " "))
  out
}

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(colour = "#1E2A32"),
      plot.title = element_text(face = "bold", size = base_size + 3, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size, colour = "#52616B", margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 1.4, colour = "#6B7780", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 0.7, colour = "#35424A"),
      axis.line = element_line(colour = "#26343D", linewidth = 0.35),
      axis.ticks = element_line(colour = "#26343D", linewidth = 0.35),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = base_size, colour = "#23313A"),
      strip.background = element_rect(fill = "#EEF3F7", colour = NA),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      plot.margin = margin(8, 10, 8, 10)
    )
}

theme_map <- function(base_size = 7.3) {
  theme_void(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2.5, colour = "#1E2A32"),
      plot.subtitle = element_text(size = base_size, colour = "#52616B", margin = margin(b = 7)),
      plot.caption = element_text(size = base_size - 1.2, colour = "#6B7780", hjust = 0),
      strip.text = element_text(face = "bold", size = base_size, colour = "#23313A"),
      strip.background = element_rect(fill = "#EEF3F7", colour = NA),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.8),
      panel.spacing = unit(0.55, "lines"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

FIGURES <- list()

save_fig <- function(plot, name, width, height, editable = TRUE, map = FALSE) {
  png_path <- file.path(OUT_DIR, paste0(name, ".png"))
  pdf_path <- file.path(OUT_DIR, paste0(name, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 450, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  FIGURES[[length(FIGURES) + 1L]] <<- list(
    name = name,
    title = stringr::str_replace_all(name, "_", " "),
    png = png_path,
    pdf = pdf_path,
    plot = plot,
    width = width,
    height = height,
    editable = editable && !map && has_rvg,
    map = map
  )
}

add_slide_with_figure <- function(ppt, fig) {
  ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- officer::ph_with(
    ppt,
    value = stringr::str_to_sentence(fig$title),
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

write_pptx_outputs <- function(figures) {
  individual <- character(length(figures))
  for (i in seq_along(figures)) {
    ppt <- officer::read_pptx()
    ppt <- add_slide_with_figure(ppt, figures[[i]])
    target <- file.path(OUT_DIR, paste0(figures[[i]]$name, ".pptx"))
    print(ppt, target = target)
    individual[i] <- target
  }

  atlas <- officer::read_pptx()
  for (fig in figures) atlas <- add_slide_with_figure(atlas, fig)
  atlas_path <- file.path(OUT_DIR, "bird_community_500sp_all_analysis_atlas_20260714_v2.pptx")
  print(atlas, target = atlas_path)
  list(individual = individual, atlas = atlas_path)
}

load_china_layers <- function() {
  province_path <- file.path(DATA_DIR, "中国shp", "省.shp")
  province <- suppressMessages(sf::st_read(stop_if_missing(province_path), quiet = TRUE, options = "ENCODING=UTF-8")) |>
    safe_make_valid() |>
    sf::st_transform(4326)
  province_lines <- sf::st_boundary(province)
  outline <- sf::st_sf(geometry = sf::st_union(sf::st_geometry(province)), crs = sf::st_crs(province))
  list(province = province, province_lines = province_lines, outline = outline, source = province_path)
}

load_grid <- function() {
  path <- c(
    file.path(DATA_DIR, "derived_v3", "china_grid_100km_v3.rds"),
    file.path(DATA_DIR, "derived_v2", "china_grid_100km_v2.rds")
  )
  path <- path[file.exists(path)][1]
  if (is.na(path)) stop("No 100-km grid found.", call. = FALSE)
  readRDS(path) |>
    safe_make_valid() |>
    sf::st_transform(4326) |>
    select(grid_cell, geometry)
}

china_layers <- load_china_layers()
grid_sf <- load_grid()
map_xlim <- c(73, 136)
map_ylim <- c(3, 54)

message("[22] Basemap source: ", china_layers$source)

metrics_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_community_metrics_with_cri_", RUN_500_EXT, ".csv")))
trend_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_trend_summary_", RUN_500_EXT, ".csv")))
beta_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_baselga_summary_", RUN_500_EXT, ".csv")))
beta_global_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_baselga_global_", RUN_500_EXT, ".csv")))
dyn_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_temporal_dynamics_summary_", RUN_500_EXT, ".csv")))
hotspot_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_species_hotspot_", RUN_500_EXT, ".csv")))
naive_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_naive_vs_corrected_", RUN_500_EXT, ".csv")))
mk_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_mann_kendall_", RUN_500_EXT, ".csv")))
class_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_species_trend_classify_", RUN_500_EXT, ".csv")))
species_trend_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_species_trend_", RUN_500_EXT, ".csv")))
species_traits_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_species_trend_traits_", RUN_500_EXT, ".csv")))
species_env_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_species_env_trend_", RUN_500_EXT, ".csv")))
trait_corr_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_trait_trend_correlation_", RUN_500_EXT, ".csv")))
env_corr_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_env_trend_correlation_", RUN_500_EXT, ".csv")))
trait_reg_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_trait_regression_coefs_", RUN_500, ".csv")))
hfi_q_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_hfi_quartile_trends_", RUN_500_EXT, ".csv")))
hfi_change_q_500 <- read_csv_req(file.path(RESULTS_DIR, paste0("table_hfi_change_quartile_trends_", RUN_500_EXT, ".csv")))
rf_var <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_variable_summary_trend.csv"))
rf_group <- read_csv_req(file.path(RESULTS_DIR, "table_rf_importance_group_summary_trend.csv"))

available_metrics <- intersect(metric_order, unique(metrics_500$metric))
if (length(available_metrics) != 20) {
  warning("Expected 20 metrics, found ", length(available_metrics), ": ", paste(available_metrics, collapse = ", "))
}
metric_groups <- split(available_metrics, ceiling(seq_along(available_metrics) / 5))
names(metric_groups) <- c("A_alpha_phylogeny_core", "B_phylo_functional_probability", "C_functional_geometry", "D_functional_evenness_cwm")[seq_along(metric_groups)]
group_display <- c(
  A_alpha_phylogeny_core = "taxonomic and phylogenetic diversity",
  B_phylo_functional_probability = "phylogenetic and functional probability metrics",
  C_functional_geometry = "functional richness, dispersion and divergence",
  D_functional_evenness_cwm = "functional evenness, Rao diversity and trait composition"
)

map_base <- list(
  geom_sf(data = china_layers$province, fill = "#F8FAFC", colour = NA),
  geom_sf(data = china_layers$province_lines, fill = NA, colour = "#D7DEE5", linewidth = 0.08),
  geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28),
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE, datum = NA)
)

make_timeslice_panel <- function(metrics_use, group_name) {
  periods <- intersect(names(period_label), unique(metrics_500$period))
  d <- metrics_500 |>
    filter(metric %in% metrics_use, period %in% periods) |>
    select(grid_cell, metric, period, value_mean)

  d_full <- tidyr::expand_grid(grid_cell = grid_sf$grid_cell, metric = metrics_use, period = periods) |>
    left_join(d, by = c("grid_cell", "metric", "period")) |>
    group_by(metric) |>
    mutate(value_scaled = safe_rescale(value_mean)) |>
    ungroup() |>
    mutate(
      metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(metrics_use)),
      period_label = factor(unname(period_label[period]), levels = unname(period_label[periods]))
    )

  d_sf <- grid_sf |> inner_join(d_full, by = "grid_cell")

  ggplot() +
    map_base +
    geom_sf(data = d_sf, aes(fill = value_scaled), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    facet_grid(metric_label ~ period_label) +
    scale_fill_gradientn(
      colours = pal_seq,
      na.value = "#F2F4F7",
      name = "Relative\nvalue",
      breaks = c(0, 0.5, 1),
      labels = c("Low", "Mid", "High")
    ) +
    labs(
      title = paste0("Spatiotemporal patterns of ", unname(group_display[group_name])),
      subtitle = "Occupancy-corrected 500-species communities across five survey periods",
      caption = "Values are scaled within each metric. South China Sea features are shown in the main map without an inset."
    ) +
    theme_map(7.0)
}

for (i in seq_along(metric_groups)) {
  p <- make_timeslice_panel(metric_groups[[i]], names(metric_groups)[i])
  save_fig(p, sprintf("fig01%s_500sp_timeslice_maps_%s", LETTERS[i], names(metric_groups)[i]), 13.6, 9.4, editable = FALSE, map = TRUE)
}

make_trend_panel <- function(metrics_use, group_name) {
  d <- trend_500 |>
    filter(metric %in% metrics_use) |>
    filter(method == "ols" | is.na(method)) |>
    group_by(grid_cell, metric) |>
    slice(1) |>
    ungroup() |>
    mutate(credible = q025 > 0 | q975 < 0) |>
    group_by(metric) |>
    mutate(trend_scaled = safe_standardize(mean)) |>
    ungroup()

  d_full <- tidyr::expand_grid(grid_cell = grid_sf$grid_cell, metric = metrics_use) |>
    left_join(select(d, grid_cell, metric, mean, q025, q975, credible, trend_scaled), by = c("grid_cell", "metric")) |>
    mutate(metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(metrics_use)))

  d_sf <- grid_sf |> inner_join(d_full, by = "grid_cell")

  ggplot() +
    map_base +
    geom_sf(data = d_sf, aes(fill = trend_scaled, alpha = credible), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    facet_wrap(~ metric_label, ncol = 5) +
    scale_fill_gradientn(colours = pal_div, limits = c(-1, 1), oob = squish, na.value = "#F2F4F7", name = "Standardized\ntrend") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.50), guide = "none", na.value = 0.35) +
    labs(
      title = paste0("Spatial fingerprints of change in ", unname(group_display[group_name])),
      subtitle = "Red indicates relative increase; blue indicates relative decrease; pale cells have 95% CrI overlapping zero"
    ) +
    theme_map(7.4)
}

for (i in seq_along(metric_groups)) {
  p <- make_trend_panel(metric_groups[[i]], names(metric_groups)[i])
  save_fig(p, sprintf("fig02%s_500sp_trend_maps_%s", LETTERS[i], names(metric_groups)[i]), 13.6, 4.9, editable = FALSE, map = TRUE)
}

make_global_trajectory <- function() {
  d <- metrics_500 |>
    filter(metric %in% available_metrics) |>
    group_by(metric, period) |>
    summarise(
      value = mean(value_mean, na.rm = TRUE),
      lo = mean(value_l95, na.rm = TRUE),
      hi = mean(value_u95, na.rm = TRUE),
      spatial_sd = sd(value_mean, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(metric, period) |>
    group_by(metric) |>
    mutate(
      base = value[match("P1", period)],
      base_sd = spatial_sd[match("P1", period)],
      index = (value - base) / base_sd,
      index_lo = (lo - base) / base_sd,
      index_hi = (hi - base) / base_sd
    ) |>
    ungroup() |>
    mutate(
      metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(available_metrics)),
      period_label = factor(unname(period_label[period]), levels = unname(period_label[names(period_label) %in% unique(period)]))
    )

  ggplot(d, aes(period_label, index, group = metric_label)) +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "#9EA8B1") +
    geom_ribbon(aes(ymin = index_lo, ymax = index_hi), fill = "#0072B2", alpha = 0.12, colour = NA) +
    geom_line(colour = "#0072B2", linewidth = 0.58) +
    geom_point(colour = "#0072B2", size = 1.35) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 5) +
    labs(
      title = "Standardized national trajectories across 20 diversity indicators",
      subtitle = "Change relative to 2000-2004, expressed in baseline spatial standard-deviation units",
      x = NULL,
      y = "Standardized change from 2000-2004"
    ) +
    theme_pub(8.0) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

save_fig(make_global_trajectory(), "fig03_500sp_all20_global_trajectories_editable", 13.2, 9.0, editable = TRUE)

make_endpoint_change <- function() {
  endpoint <- metrics_500 |>
    filter(metric %in% available_metrics, period %in% c("P1", "P5")) |>
    group_by(metric, period) |>
    summarise(value = mean(value_mean, na.rm = TRUE), spatial_sd = sd(value_mean, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = period, values_from = c(value, spatial_sd)) |>
    mutate(
      standardized_change = (value_P5 - value_P1) / spatial_sd_P1,
      metric_label = factor(nice_metric(metric), levels = nice_metric(metric[order(standardized_change)])),
      direction = if_else(standardized_change >= 0, "Increase", "Decrease")
    )

  ggplot(endpoint, aes(standardized_change, metric_label, colour = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9EA8B1", linetype = "22") +
    geom_segment(aes(x = 0, xend = standardized_change, yend = metric_label), linewidth = 0.65) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(Increase = "#0072B2", Decrease = "#D55E00"), name = NULL) +
    labs(
      title = "Standardized net change across 20 diversity indicators",
      subtitle = "National change from 2000-2004 to 2020-2024, scaled by baseline spatial variation",
      x = "Standardized change (P5 - P1) / baseline SD",
      y = NULL
    ) +
    theme_pub(8.4)
}

save_fig(make_endpoint_change(), "fig04_500sp_all20_endpoint_change_editable", 10.8, 7.6, editable = TRUE)

make_mann_kendall <- function() {
  d <- mk_500 |>
    filter(metric %in% available_metrics) |>
    mutate(sig = mk_p < 0.05, direction = case_when(mk_tau > 0 ~ "Positive", mk_tau < 0 ~ "Negative", TRUE ~ "Flat")) |>
    group_by(metric) |>
    summarise(
      mean_tau = mean(mk_tau, na.rm = TRUE),
      prop_sig_positive = mean(sig & mk_tau > 0, na.rm = TRUE),
      prop_sig_negative = mean(sig & mk_tau < 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(metric_label = factor(nice_metric(metric), levels = nice_metric(metric[order(mean_tau)])))

  p_tau <- ggplot(d, aes(mean_tau, metric_label, colour = mean_tau)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9EA8B1", linetype = "22") +
    geom_point(size = 2.1) +
    scale_colour_gradient2(low = "#0072B2", mid = "#7F8790", high = "#D55E00", midpoint = 0, guide = "none") +
    labs(title = "Mean Mann-Kendall tau", x = "Mean tau across grid cells", y = NULL) +
    theme_pub(8.2)

  sig_d <- d |>
    select(metric_label, prop_sig_positive, prop_sig_negative) |>
    pivot_longer(starts_with("prop"), names_to = "component", values_to = "prop") |>
    mutate(component = recode(component, prop_sig_positive = "Significant positive", prop_sig_negative = "Significant negative"))

  p_sig <- ggplot(sig_d, aes(prop, metric_label, fill = component)) +
    geom_col(position = "stack", width = 0.62, colour = "white", linewidth = 0.2) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("Significant positive" = "#D55E00", "Significant negative" = "#0072B2"), name = NULL) +
    labs(title = "Grid-cell significance share", x = "Share of grid cells", y = NULL) +
    theme_pub(8.2)

  p_tau + p_sig + plot_layout(widths = c(1, 1.1)) +
    plot_annotation(
      title = "Mann-Kendall trend diagnostics for all 20 indicators",
      subtitle = "Non-parametric monotonic trend check across five periods",
      theme = theme_pub(9)
    )
}

save_fig(make_mann_kendall(), "fig05_500sp_all20_mann_kendall_diagnostics_editable", 12.4, 8.2, editable = TRUE)

make_beta_maps <- function() {
  metrics_use <- c("beta_sor", "beta_sim", "beta_sne")
  metric_names <- c(beta_sor = "Total beta diversity", beta_sim = "Turnover", beta_sne = "Nestedness")
  pairs <- unique(beta_500$period_pair)

  d_full <- tidyr::expand_grid(grid_cell = grid_sf$grid_cell, metric = metrics_use, period_pair = pairs) |>
    left_join(beta_500 |> filter(metric %in% metrics_use) |> select(grid_cell, metric, period_pair, mean),
              by = c("grid_cell", "metric", "period_pair")) |>
    mutate(
      metric_label = factor(unname(metric_names[metric]), levels = unname(metric_names[metrics_use])),
      period_label = factor(stringr::str_replace_all(period_pair, "_", " to "), levels = stringr::str_replace_all(pairs, "_", " to "))
    )

  d_sf <- grid_sf |> inner_join(d_full, by = "grid_cell")

  ggplot() +
    map_base +
    geom_sf(data = d_sf, aes(fill = mean), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    facet_grid(metric_label ~ period_label) +
    scale_fill_gradientn(colours = pal_blue, na.value = "#F2F4F7", name = "Beta\ndiversity") +
    labs(
      title = "500-species temporal beta diversity decomposition",
      subtitle = "Baselga total beta, turnover, and nestedness across consecutive periods"
    ) +
    theme_map(7.0)
}

save_fig(make_beta_maps(), "fig06_500sp_temporal_beta_maps", 13.2, 7.8, editable = FALSE, map = TRUE)

make_beta_global <- function() {
  d <- beta_global_500 |>
    mutate(
      period_pair = factor(stringr::str_replace_all(period_pair, "_", " to "), levels = stringr::str_replace_all(period_pair, "_", " to ")),
      nested_mean = 1 - prop_turnover_mean,
      nested_q025 = 1 - prop_turnover_q975,
      nested_q975 = 1 - prop_turnover_q025
    ) |>
    select(period_pair, prop_turnover_mean, prop_turnover_q025, prop_turnover_q975, nested_mean, nested_q025, nested_q975) |>
    pivot_longer(c(prop_turnover_mean, nested_mean), names_to = "component", values_to = "mean") |>
    mutate(
      q025 = if_else(component == "prop_turnover_mean", prop_turnover_q025, nested_q025),
      q975 = if_else(component == "prop_turnover_mean", prop_turnover_q975, nested_q975),
      component = recode(component, prop_turnover_mean = "Turnover", nested_mean = "Nestedness")
    )

  ggplot(d, aes(period_pair, mean, fill = component)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = "white", linewidth = 0.25) +
    geom_errorbar(aes(ymin = q025, ymax = q975), position = position_dodge(width = 0.7), width = 0.17, linewidth = 0.38) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c(Turnover = "#0072B2", Nestedness = "#E69F00"), name = NULL) +
    labs(
      title = "500-species global beta-diversity composition",
      subtitle = "Turnover share versus nestedness share across consecutive periods",
      x = NULL,
      y = "Share of total beta diversity"
    ) +
    theme_pub(8.8)
}

save_fig(make_beta_global(), "fig07_500sp_global_baselga_turnover_nestedness_editable", 9.2, 5.6, editable = TRUE)

make_temporal_dynamics_maps <- function() {
  metrics_use <- sort(unique(dyn_500$metric))
  d <- dyn_500 |>
    filter(metric %in% metrics_use) |>
    group_by(metric) |>
    mutate(value_scaled = safe_rescale(mean)) |>
    ungroup() |>
    mutate(metric_label = factor(nice_term(metric), levels = nice_term(metrics_use)))

  d_sf <- grid_sf |> inner_join(d, by = "grid_cell")

  ggplot() +
    map_base +
    geom_sf(data = d_sf, aes(fill = value_scaled), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    facet_wrap(~ metric_label, ncol = 5) +
    scale_fill_gradientn(colours = pal_seq, na.value = "#F2F4F7", name = "Relative\nvalue") +
    labs(
      title = "500-species temporal-dynamics diagnostics",
      subtitle = "Synchrony, turnover gain/loss/total, and variance ratio"
    ) +
    theme_map(7.4)
}

save_fig(make_temporal_dynamics_maps(), "fig08_500sp_temporal_dynamics_maps", 13.0, 4.7, editable = FALSE, map = TRUE)

make_hotspot_maps <- function() {
  hs_metrics <- c("expanding_richness", "contracting_richness", "net_trend_index", "turnover_balance")
  hs_names <- c(
    expanding_richness = "Increasing-species richness",
    contracting_richness = "Declining-species richness",
    net_trend_index = "Net trend index",
    turnover_balance = "Turnover balance"
  )
  d <- hotspot_500 |>
    filter(period %in% names(period_label)) |>
    pivot_longer(all_of(hs_metrics), names_to = "metric", values_to = "value") |>
    group_by(metric) |>
    mutate(value_scaled = safe_rescale(value)) |>
    ungroup() |>
    mutate(
      metric_label = factor(unname(hs_names[metric]), levels = unname(hs_names[hs_metrics])),
      period_label = factor(unname(period_label[period]), levels = unname(period_label[names(period_label) %in% unique(period)]))
    )

  d_sf <- grid_sf |> inner_join(d, by = "grid_cell")

  ggplot() +
    map_base +
    geom_sf(data = d_sf, aes(fill = value_scaled), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    facet_grid(metric_label ~ period_label) +
    scale_fill_gradientn(colours = pal_seq, na.value = "#F2F4F7", name = "Relative\nvalue") +
    labs(
      title = "500-species species-turnover hotspots",
      subtitle = "Spatial balance of increasing and declining species across five periods"
    ) +
    theme_map(7.0)
}

save_fig(make_hotspot_maps(), "fig09_500sp_species_turnover_hotspot_maps", 13.4, 7.6, editable = FALSE, map = TRUE)

make_naive_corrected <- function() {
  map_d <- grid_sf |>
    left_join(naive_500, by = "grid_cell")

  p_scatter <- ggplot(naive_500, aes(naive_trend, corrected_trend)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.42, colour = "#8F9BA4", linetype = "22") +
    geom_point(aes(colour = direction_flipped), size = 1.2, alpha = 0.58) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, colour = "#23313A") +
    scale_colour_manual(values = c(`FALSE` = "#61707A", `TRUE` = "#D55E00"), labels = c("No flip", "Direction flip"), name = NULL) +
    labs(title = "Naive versus corrected trends", x = "Naive observed richness trend", y = "Occupancy-corrected trend") +
    theme_pub(8.3)

  p_diff <- ggplot() +
    map_base +
    geom_sf(data = map_d, aes(fill = trend_diff), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    scale_fill_gradient2(low = "#0072B2", mid = "#F7F7F7", high = "#D55E00", midpoint = 0, na.value = "#F2F4F7", name = "Corrected - naive") +
    labs(title = "Bias-correction difference") +
    theme_map(7.4)

  p_flip <- ggplot() +
    map_base +
    geom_sf(data = map_d, aes(fill = direction_flipped), colour = NA) +
    geom_sf(data = china_layers$province_lines, fill = NA, colour = "#DFE5EA", linewidth = 0.06) +
    geom_sf(data = china_layers$outline, fill = NA, colour = "#23313A", linewidth = 0.28) +
    scale_fill_manual(values = c(`FALSE` = "#D9E2EA", `TRUE` = "#D55E00"), na.value = "#F2F4F7", name = "Direction flip") +
    labs(title = "Directional reversals") +
    theme_map(7.4)

  p_scatter | (p_diff / p_flip) +
    plot_layout(widths = c(0.95, 1.25)) +
    plot_annotation(
      title = "Detection correction reshapes 500-species richness trends",
      subtitle = paste0("Direction flips: ", sum(naive_500$direction_flipped, na.rm = TRUE), " / ", nrow(naive_500), " grid cells"),
      theme = theme_pub(9)
    )
}

save_fig(make_naive_corrected(), "fig10_500sp_naive_vs_corrected_bias_maps", 12.6, 7.0, editable = FALSE, map = TRUE)

make_species_classes <- function() {
  class_counts <- class_500 |>
    mutate(class = unname(class_label[trend_class]), class = if_else(is.na(class), stringr::str_to_sentence(trend_class), class)) |>
    count(class, name = "n") |>
    mutate(prop = n / sum(n), class = factor(class, levels = c("Declining", "Stable", "Increasing")))

  top_sp <- species_trend_500 |>
    arrange(desc(mean)) |>
    slice_head(n = 12) |>
    bind_rows(species_trend_500 |> arrange(mean) |> slice_head(n = 12)) |>
    distinct(species, .keep_all = TRUE) |>
    mutate(direction = if_else(mean >= 0, "Increasing", "Declining"), species = factor(species, levels = species[order(mean)]))

  p_class <- ggplot(class_counts, aes(class, prop, fill = class)) +
    geom_col(width = 0.64, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = paste0(n, "\n", percent(prop, accuracy = 0.1))), vjust = 0.5, size = 3.0, colour = "white") +
    scale_fill_manual(values = pal_class, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
    labs(title = "Trend classes", x = NULL, y = "Species share") +
    theme_pub(8.6)

  p_species <- ggplot(top_sp, aes(mean, species, colour = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.36, colour = "#9EA8B1") +
    geom_errorbar(aes(xmin = q025, xmax = q975), orientation = "y", width = 0, linewidth = 0.43) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = pal_class[c("Increasing", "Declining")], name = NULL) +
    labs(title = "Strongest species shifts", x = "Species occupancy trend", y = NULL) +
    theme_pub(8.4) +
    theme(axis.text.y = element_text(face = "italic", size = 7.2))

  p_class + p_species + plot_layout(widths = c(0.72, 1.45)) +
    plot_annotation(
      title = "500-species winners, stable taxa, and losers",
      subtitle = "Species-level temporal trend classification",
      theme = theme_pub(9)
    )
}

save_fig(make_species_classes(), "fig11_500sp_species_trend_classes_editable", 12.2, 6.8, editable = TRUE)

make_rf_plot <- function() {
  var_d <- rf_var |>
    mutate(
      group_label = unname(group_label[group]),
      group_label = if_else(is.na(group_label), stringr::str_to_sentence(group), group_label),
      label_plot = if_else(!is.na(label), label, nice_term(variable)),
      label_plot = stringr::str_wrap(label_plot, 28)
    ) |>
    arrange(desc(mean)) |>
    slice_head(n = 20) |>
    mutate(label_plot = factor(label_plot, levels = rev(label_plot)))

  group_d <- rf_group |>
    mutate(group_label = if_else(!is.na(label), label, unname(group_label[group]))) |>
    mutate(group_label = stringr::str_wrap(group_label, 24)) |>
    arrange(mean) |>
    mutate(group_label = factor(group_label, levels = group_label))

  p_var <- ggplot(var_d, aes(mean, label_plot, colour = group_label)) +
    geom_errorbar(aes(xmin = q_lo, xmax = q_hi), orientation = "y", width = 0, linewidth = 0.42, alpha = 0.88) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = pal_driver, name = NULL, na.value = "#7F8790") +
    labs(title = "Top predictors", x = "Permutation importance", y = NULL) +
    theme_pub(8.1) +
    theme(axis.text.y = element_text(size = 7.0))

  p_group <- ggplot(group_d, aes(mean, group_label, fill = group_label)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
    geom_errorbar(aes(xmin = q_lo, xmax = q_hi), orientation = "y", width = 0.16, linewidth = 0.36, colour = "#23313A") +
    scale_fill_manual(values = pal_driver, guide = "none", na.value = "#7F8790") +
    labs(title = "Predictor groups", x = "Grouped importance", y = NULL) +
    theme_pub(8.1)

  p_var + p_group + plot_layout(widths = c(1.45, 0.85)) +
    plot_annotation(
      title = "500-species random-forest driver importance",
      subtitle = "Non-linear predictor ranking across posterior draws",
      theme = theme_pub(9)
    )
}

save_fig(make_rf_plot(), "fig12_500sp_rf_driver_importance_editable", 12.1, 7.0, editable = TRUE)

make_driver_coefficients <- function() {
  files <- list.files(RESULTS_DIR, pattern = paste0("^table_brms_driver_trend_change_coefficients_.*_", RUN_500_EXT, "\\.csv$"), full.names = TRUE)
  d <- lapply(files, function(f) {
    response <- basename(f) |>
      str_remove("^table_brms_driver_trend_change_coefficients_") |>
      str_remove(paste0("_", RUN_500_EXT, "\\.csv$"))
    read_csv(f, show_col_types = FALSE) |> mutate(response = response)
  }) |>
    bind_rows() |>
    filter(term != "Intercept") |>
    mutate(
      term_label = factor(nice_term(term), levels = rev(unique(nice_term(term)))),
      response_label = factor(nice_metric_short(response), levels = nice_metric_short(unique(response))),
      credible = Q2.5 > 0 | Q97.5 < 0
    )

  ggplot(d, aes(Estimate, term_label, colour = credible)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9EA8B1", linetype = "22") +
    geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), orientation = "y", width = 0, linewidth = 0.40) +
    geom_point(size = 1.55) +
    facet_wrap(~ response_label, scales = "free_x", ncol = 4) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#61707A"), labels = c("CrI overlaps 0", "CrI excludes 0"), name = NULL) +
    labs(
      title = "500-species posterior driver coefficients",
      subtitle = "Temporal-change covariates fitted to occupancy-corrected trend estimates",
      x = "Posterior coefficient",
      y = NULL
    ) +
    theme_pub(8.0) +
    theme(axis.text.y = element_text(size = 7.0))
}

save_fig(make_driver_coefficients(), "fig13_500sp_driver_coefficients_editable", 12.8, 7.4, editable = TRUE)

make_varpart <- function() {
  files <- list.files(RESULTS_DIR, pattern = paste0("^table_varpart_.*_", RUN_500_EXT, "\\.csv$"), full.names = TRUE)
  d <- lapply(files, function(f) {
    response <- basename(f) |>
      str_remove("^table_varpart_") |>
      str_remove(paste0("_", RUN_500_EXT, "\\.csv$"))
    read_csv(f, show_col_types = FALSE) |> mutate(response = response)
  }) |>
    bind_rows() |>
    filter(fraction %in% names(varpart_label)) |>
    mutate(
      component = factor(unname(varpart_label[fraction]), levels = c("Climate change", "Land-use change", "Human pressure change", "Spatial baseline")),
      response_label = factor(nice_metric_short(response), levels = nice_metric_short(unique(response))),
      adj_r2_pos = pmax(adj_r2, 0)
    )

  ggplot(d, aes(adj_r2_pos, response_label, fill = component)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.22) +
    scale_fill_manual(values = pal_driver, name = NULL) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "500-species linear variance partitioning",
      subtitle = "Pure adjusted R2 fractions; negative adjusted fractions truncated at zero for display",
      x = "Pure adjusted R2",
      y = NULL
    ) +
    theme_pub(8.4)
}

save_fig(make_varpart(), "fig14_500sp_varpart_driver_groups_editable", 10.8, 6.4, editable = TRUE)

make_trait_environment <- function() {
  trait_d <- trait_corr_500 |>
    mutate(
      label = factor(nice_term(trait), levels = rev(nice_term(trait))),
      credible = rho_q025 > 0 | rho_q975 < 0
    )

  env_d <- env_corr_500 |>
    mutate(label = factor(nice_term(env_var), levels = nice_term(env_var)[order(spearman_rho)]))

  reg_d <- trait_reg_500 |>
    filter(term != "Intercept") |>
    mutate(
      label = factor(nice_term(stringr::str_remove(term, "^z_")), levels = rev(nice_term(stringr::str_remove(term, "^z_")))),
      credible = q025 > 0 | q975 < 0
    )

  p_trait <- ggplot(trait_d, aes(rho_mean, label, colour = credible)) +
    geom_vline(xintercept = 0, colour = "#9EA8B1", linewidth = 0.35, linetype = "22") +
    geom_errorbar(aes(xmin = rho_q025, xmax = rho_q975), orientation = "y", width = 0, linewidth = 0.48) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#61707A"), guide = "none") +
    labs(title = "Trait correlations", x = "Spearman rho", y = NULL) +
    theme_pub(8.3)

  p_env <- ggplot(env_d, aes(spearman_rho, label)) +
    geom_vline(xintercept = 0, colour = "#9EA8B1", linewidth = 0.35, linetype = "22") +
    geom_col(width = 0.58, fill = "#0072B2", alpha = 0.88) +
    geom_text(aes(label = paste0("n=", n_species)), hjust = -0.12, size = 2.35, colour = "#52616B") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(title = "Environmental correlations", x = "Spearman rho", y = NULL) +
    theme_pub(8.3)

  p_reg <- ggplot(reg_d, aes(estimate, label, colour = credible)) +
    geom_vline(xintercept = 0, colour = "#9EA8B1", linewidth = 0.35, linetype = "22") +
    geom_errorbar(aes(xmin = q025, xmax = q975), orientation = "y", width = 0, linewidth = 0.48) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#61707A"), guide = "none") +
    labs(title = "Trait regression", x = "Coefficient", y = NULL) +
    theme_pub(8.3)

  p_trait + p_env + p_reg + plot_layout(widths = c(0.85, 1.15, 0.85)) +
    plot_annotation(
      title = "500-species trait and environmental mechanisms",
      subtitle = "Complete-case association screens for species-level trend variation",
      theme = theme_pub(9)
    )
}

save_fig(make_trait_environment(), "fig15_500sp_trait_environment_mechanisms_editable", 13.0, 5.6, editable = TRUE)

make_hfi_quartile <- function(dat, quartile_col, title, subtitle) {
  metrics_use <- intersect(c("corrected_richness", "shannon", "pd_prob", "trait_volume"), unique(dat$metric))
  metrics_use <- metrics_use[vapply(metrics_use, function(m) any(is.finite(dat$mean[dat$metric == m])), logical(1))]
  d <- dat |>
    filter(metric %in% metrics_use, method == "ols" | is.na(method), !is.na(.data[[quartile_col]])) |>
    group_by(metric, .data[[quartile_col]]) |>
    summarise(
      mean = mean(mean, na.rm = TRUE),
      lo = mean(q025, na.rm = TRUE),
      hi = mean(q975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(quartile = !!quartile_col) |>
    mutate(
      quartile = factor(quartile, levels = unique(quartile[order(quartile)])),
      metric_label = factor(nice_metric_short(metric), levels = nice_metric_short(metrics_use))
    )

  ggplot(d, aes(quartile, mean, colour = metric_label, group = metric_label)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "#9EA8B1", linetype = "22") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.13, linewidth = 0.40) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.8) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = c("#0072B2", "#009E73", "#E69F00", "#CC79A7"), guide = "none") +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean trend") +
    theme_pub(8.2) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

save_fig(
  make_hfi_quartile(hfi_q_500, "hfi_q", "500-species trends across human-footprint quartiles", "Four headline indicators summarized by baseline HFI quartile"),
  "fig16_500sp_hfi_quartile_trends_editable",
  12.0, 5.8, editable = TRUE
)

save_fig(
  make_hfi_quartile(hfi_change_q_500, "hfi_change_q", "500-species trends across human-footprint change quartiles", "Four headline indicators summarized by HFI-change quartile"),
  "fig17_500sp_hfi_change_quartile_trends_editable",
  12.0, 5.8, editable = TRUE
)

make_species_env_profiles <- function() {
  env_vars <- intersect(c("bio4", "bio7", "bio11", "bio13", "elev_mean", "hfi_mean", "landcover_built", "landcover_cropland"), names(species_env_500))
  d <- species_env_500 |>
    filter(trend_class %in% c("expanding", "stable", "contracting")) |>
    mutate(class = factor(unname(class_label[trend_class]), levels = c("Declining", "Stable", "Increasing"))) |>
    pivot_longer(all_of(env_vars), names_to = "env_var", values_to = "value") |>
    filter(is.finite(value)) |>
    mutate(env_label = factor(nice_term(env_var), levels = nice_term(env_vars)))

  ggplot(d, aes(class, value, fill = class)) +
    geom_boxplot(width = 0.62, outlier.size = 0.35, linewidth = 0.28) +
    facet_wrap(~ env_label, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = pal_class, guide = "none") +
    labs(
      title = "Environmental profiles of increasing, stable, and declining species",
      subtitle = "500-species complete-case environmental summaries by trend class",
      x = NULL,
      y = "Species-level mean environment"
    ) +
    theme_pub(8.0) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

save_fig(make_species_env_profiles(), "fig18_500sp_species_environment_profiles_editable", 12.4, 7.4, editable = TRUE)

pptx <- write_pptx_outputs(FIGURES)

manifest <- tibble(
  figure = vapply(FIGURES, `[[`, character(1), "name"),
  png = vapply(FIGURES, `[[`, character(1), "png"),
  pdf = vapply(FIGURES, `[[`, character(1), "pdf"),
  pptx = pptx$individual,
  pptx_editable = vapply(
    FIGURES,
    function(x) if (isTRUE(x$editable)) "editable_drawingml" else if (isTRUE(x$map)) "map_embedded_png" else "embedded_png",
    character(1)
  )
)
readr::write_csv(manifest, file.path(OUT_DIR, "figure_manifest_500sp_all_analysis_20260714_v2.csv"))

message("[22] Wrote ", length(FIGURES), " figures.")
message("[22] Atlas: ", pptx$atlas)
message("[22] Manifest: ", file.path(OUT_DIR, "figure_manifest_500sp_all_analysis_20260714_v2.csv"))
