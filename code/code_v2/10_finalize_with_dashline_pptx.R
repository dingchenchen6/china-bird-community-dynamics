#!/usr/bin/env Rscript
## 10_finalize_with_dashline_pptx.R
##
## 1) 重画 v2 全部核心图，地图统一加十段线（不加鹰眼图）；
## 2) 森林图：补一份 ridgeline 山脊图版本；
## 3) 集成可编辑 PPTX deck（每图一页 + 总目录）。
##
## 默认 V2_RUN_LABEL=v2_full_200sp_ar1。
## 输出：figures_v2/*_v3.{png,pdf,pptx} + figures_v2/v2_master_deck_<LABEL>.pptx

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(forcats); library(ggplot2); library(patchwork); library(sf)
  library(ggridges); library(officer); library(rvg)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_spatial.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_plots_advanced.R"))

P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")

message(sprintf("[stage-10] regenerating maps with 十段线 + PPTX export | %s",
                RUN_LABEL))

china_layers <- load_china_layers(P$china_boundary, P$province_line)
stopifnot(!is.null(china_layers$ten_dash))      # must have dashline
grid_sf <- readRDS(file.path(P$derived_v2, "china_grid_100km_v2.rds"))
primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
# 分两个 deck：地图用 PNG 嵌入（稳定）、非地图用 vector dml（可在 PPT 内编辑）
deck_maps  <- list()
deck_stats <- list()

add_plot <- function(name, plot, w = 12, h = 7.2, is_map = FALSE) {
  save_dual(plot, paste0(name, "_v3"), width = w, height = h)
  if (is_map) {
    deck_maps[[name]] <<- plot           # PNG 嵌入到 master deck
  } else {
    # 单独可编辑 PPTX
    P <- v2_paths()
    if (requireNamespace("officer", quietly = TRUE) &&
        requireNamespace("rvg",     quietly = TRUE)) {
      ok <- tryCatch({
        d <- officer::read_pptx() |>
          officer::add_slide(layout = "Blank", master = "Office Theme") |>
          officer::ph_with(value = rvg::dml(ggobj = plot),
                            location = officer::ph_location(left = 0.4, top = 0.4,
                                                            width = w, height = h))
        print(d, target = file.path(P$figures_v2,
                                     paste0(name, "_v3.pptx")))
        TRUE
      }, error = function(e) { message("  pptx fail (", name, "): ",
                                        conditionMessage(e)); FALSE })
    }
    deck_stats[[name]] <<- plot          # vector deck
  }
}

## --- A. multidiversity time-slice (6 metric × 5 period) -------------------

m_long <- read_csv(v2_file("results",
            paste0("table_community_metrics_with_cri_", RUN_LABEL)),
            show_col_types = FALSE)
show_metrics <- c("corrected_richness","shannon","pd_prob","mpd_prob",
                  "trait_volume","rao_q")
m_show <- m_long |>
  filter(metric %in% show_metrics) |>
  mutate(
    metric_label = factor(recode(metric,
      corrected_richness = "Taxonomic richness",
      shannon            = "Taxonomic Shannon",
      pd_prob            = "Faith's PD",
      mpd_prob           = "Phylogenetic MPD",
      trait_volume       = "Functional trait volume",
      rao_q              = "Functional Rao's Q"),
      levels = c("Taxonomic richness","Taxonomic Shannon",
                  "Faith's PD","Phylogenetic MPD",
                  "Functional trait volume","Functional Rao's Q")),
    block_label = factor(block_label, levels = primary_blocks$block_label)) |>
  filter(!is.na(metric_label), !is.na(block_label)) |>
  group_by(metric_label) |>
  mutate(value_z = pmin(pmax(as.numeric(scale(value_mean)), -2.5), 2.5)) |>
  ungroup()
m_sf <- grid_sf |> inner_join(m_show, by = "grid_cell")

p_ts <- ggplot(m_sf) +
  geom_sf(aes(fill = value_z),
          colour = scales::alpha("white", 0.08), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_v2_diverging(name = "Within-metric z-score (clipped at +/-2.5)") +
  facet_grid(metric_label ~ block_label) +
  v2_china_coord() +
  theme_v2_map(9.4) +
  labs(title = "Occupancy-corrected multidiversity time-slice maps",
       subtitle = sprintf("Posterior mean | run=%s | with 十段线", RUN_LABEL),
       caption = "Each row z-scaled within metric.")
add_plot("multidiversity_timeslices", p_ts, w = 14, h = 12, is_map = TRUE)

## --- B. community trends + sig95 outline ---------------------------------

trd <- read_csv(v2_file("results",
                paste0("table_grid_trends_with_cri_", RUN_LABEL)),
                show_col_types = FALSE)
panel_metrics <- c("corrected_richness","shannon","trait_volume","pd_prob")
trd_show <- trd |>
  filter(metric %in% panel_metrics) |>
  mutate(metric_label = recode(metric,
    corrected_richness = "Richness trend",
    shannon = "Shannon trend",
    trait_volume = "Trait-volume trend",
    pd_prob = "Faith's PD trend")) |>
  group_by(metric_label) |>
  mutate(z = pmin(pmax(as.numeric(scale(trend_mean)), -2.5), 2.5)) |>
  ungroup()
trd_sf <- grid_sf |> inner_join(trd_show, by = "grid_cell")
p_tr <- ggplot(trd_sf) +
  geom_sf(aes(fill = z), colour = scales::alpha("white", 0.08),
          linewidth = 0.04) +
  geom_sf(data = trd_sf |> filter(sig95),
          aes(geometry = geometry), inherit.aes = FALSE, fill = NA,
          colour = "#222222", linewidth = 0.10) +
  china_map_layers(china_layers) +
  scale_fill_v2_diverging(name = "Within-metric z (95% CRI excl. 0 outlined)") +
  facet_wrap(~ metric_label, ncol = 2) +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = "Community dynamic trends with posterior uncertainty",
       subtitle = sprintf("Per-grid linear slope across 5 periods | run=%s",
                          RUN_LABEL))
add_plot("community_trends_with_cri", p_tr, w = 11, h = 9.4, is_map = TRUE)

## --- C. Temporal beta Baselga ---------------------------------------------

bt <- read_csv(v2_file("results",
                paste0("table_temporal_beta_with_cri_", RUN_LABEL)),
                show_col_types = FALSE)
bt_grid <- bt |>
  filter(metric %in% c("turnover","nestedness","bray")) |>
  group_by(grid_cell, metric) |>
  summarise(value = mean(value_mean, na.rm = TRUE), .groups = "drop") |>
  mutate(metric_label = recode(metric,
    turnover = "Baselga turnover",
    nestedness = "Baselga nestedness",
    bray = "Bray-Curtis"))
bt_sf <- grid_sf |> inner_join(bt_grid, by = "grid_cell")
p_bt <- ggplot(bt_sf) +
  geom_sf(aes(fill = value), colour = scales::alpha("white", 0.08),
          linewidth = 0.04) +
  china_map_layers(china_layers) +
  scale_fill_v2_sequential(name = "Mean dissimilarity",
                           palette = "lajolla", direction = 1) +
  facet_wrap(~ metric_label, ncol = 3) +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = "Temporal beta diversity decomposition",
       subtitle = sprintf("Baselga + Bray | run=%s", RUN_LABEL))
add_plot("temporal_beta_baselga", p_bt, w = 14, h = 6, is_map = TRUE)

## --- D. species gain/loss --------------------------------------------------

gl <- read_csv(v2_file("results",
                paste0("table_species_gain_loss_", RUN_LABEL)),
                show_col_types = FALSE) |>
  group_by(grid_cell) |>
  summarise(mean_gain = mean(gain), mean_loss = mean(loss),
            mean_net  = mean(net),
            mean_total = mean(gain + loss), .groups = "drop")
gl_sf <- grid_sf |> inner_join(gl, by = "grid_cell") |>
  pivot_longer(c(mean_gain, mean_loss, mean_net, mean_total),
                names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
    mean_gain = "Mean species gain",
    mean_loss = "Mean species loss",
    mean_net = "Mean net change",
    mean_total = "Mean total turnover"))
p_gl <- ggplot(gl_sf) +
  geom_sf(aes(fill = value), colour = scales::alpha("white", 0.10),
          linewidth = 0.04) +
  china_map_layers(china_layers) +
  scale_fill_v2_diverging(name = "Sum of psi change") +
  facet_wrap(~ metric, ncol = 2) +
  v2_china_coord() +
  theme_v2_map(10.5) +
  labs(title = "Per-grid species gain / loss / net / total turnover",
       subtitle = sprintf("Probability-weighted | run=%s", RUN_LABEL))
add_plot("species_gain_loss", p_gl, w = 11, h = 9, is_map = TRUE)

## --- E. community typology -----------------------------------------------

types_grids <- read_csv(v2_file("results",
                paste0("table_community_types_grids_", RUN_LABEL)),
                show_col_types = FALSE)
types_centers <- read_csv(v2_file("results",
                paste0("table_community_types_centers_", RUN_LABEL)),
                show_col_types = FALSE)
types_sf <- grid_sf |> inner_join(types_grids, by = "grid_cell")
n_types <- length(unique(types_grids$cluster))
type_colours <- V2_PALETTES$qualitative[seq_len(n_types)]

p_types_map <- ggplot(types_sf) +
  geom_sf(aes(fill = cluster), colour = scales::alpha("white", 0.10),
          linewidth = 0.04) +
  china_map_layers(china_layers) +
  scale_fill_manual(values = type_colours, name = "Community type") +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(title = "China bird community typology",
       subtitle = sprintf("k-means on z-scaled diversity profile | run=%s",
                          RUN_LABEL))

profile_long <- types_centers |>
  pivot_longer(-cluster_id, names_to = "metric", values_to = "z") |>
  mutate(metric = recode(metric,
    corrected_richness = "Richness", shannon = "Shannon",
    trait_volume = "Trait vol", rao_q = "Rao Q",
    pd_prob = "PD", mpd_prob = "MPD"))
p_profiles <- ggplot(profile_long, aes(metric, z, fill = z > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey55") +
  facet_wrap(~ cluster_id, ncol = n_types) +
  scale_fill_manual(values = c(`TRUE` = "#8B2E1E", `FALSE` = "#0E5A78"),
                     guide = "none") +
  labs(title = "Z-profile of each community type",
       x = NULL, y = "z (sd units)") +
  theme_v2_pub(10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p_types <- p_types_map / p_profiles + patchwork::plot_layout(heights = c(2.6, 1))
add_plot("community_types", p_types, w = 12, h = 11, is_map = TRUE)

## --- F. driver coefficient: 雨林 + 山脊（升级 forest plot） ---------------

fit_files <- list.files(P$derived_v2,
  pattern = sprintf("^brms_driver_fit_.*_%s\\.rds$", RUN_LABEL),
  full.names = TRUE)
if (length(fit_files) > 0) {
  resp_pretty <- c(
    trend_corrected_richness = "Richness trend",
    trend_shannon            = "Shannon trend",
    trend_pd_prob            = "Faith's PD trend",
    trend_trait_volume       = "Trait-volume trend"
  )
  draws_all <- map_dfr(fit_files, function(p) {
    nm <- sub(sprintf("^brms_driver_fit_(.*)_%s\\.rds$", RUN_LABEL),
              "\\1", basename(p))
    fit <- readRDS(p)
    d <- brms_fixef_draws(fit)
    d$term <- gsub("^z_", "", d$term)
    d$response <- resp_pretty[[nm]] %||% nm
    d
  })

  med_within <- draws_all |>
    group_by(response, term) |>
    summarise(med = median(value), .groups = "drop")
  draws_all <- draws_all |>
    left_join(med_within, by = c("response","term")) |>
    mutate(term = forcats::fct_reorder(term, med))

  # 雨林 multipanel
  p_rain <- ggplot(draws_all, aes(value, term)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55",
               linewidth = 0.4) +
    ggdist::stat_halfeye(.width = c(0.5, 0.95), thickness = 0.55,
                         slab_alpha = 0.55, slab_colour = NA,
                         interval_colour = "grey25", point_size = 1.5,
                         fill = "#0E5A78") +
    facet_wrap(~ response, scales = "free", ncol = 2) +
    labs(title = "Drivers of community dynamic trends — raincloud",
         subtitle = sprintf("brms + cmdstanr | 4 chains | %s", RUN_LABEL),
         x = "Standardised coefficient (95% CRI)", y = NULL) +
    theme_v2_pub(11)
  add_plot("driver_raincloud_multipanel", p_rain, w = 12, h = 8)

  # 山脊 ridgeline multipanel
  p_ridge <- ggplot(draws_all, aes(value, term, fill = after_stat(x))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55",
               linewidth = 0.4) +
    ggridges::geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.005,
                                            colour = "white", linewidth = 0.2,
                                            alpha = 0.95) +
    facet_wrap(~ response, scales = "free", ncol = 2) +
    scale_fill_gradient2(low = "#0E5A78", mid = "#F2E8D8",
                          high = "#8B2E1E", midpoint = 0, guide = "none") +
    labs(title = "Drivers of community dynamic trends — ridgeline",
         subtitle = sprintf("brms + cmdstanr | 4 chains | %s", RUN_LABEL),
         x = "Standardised coefficient", y = NULL) +
    theme_v2_pub(11)
  add_plot("driver_ridgeline_multipanel", p_ridge, w = 12, h = 8)

  # 蜂群 overlay
  p_bee <- ggplot(draws_all, aes(value, term, colour = response)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggbeeswarm::geom_quasirandom(alpha = 0.18, size = 0.4, width = 0.32,
                                  groupOnX = FALSE, show.legend = FALSE) +
    geom_point(data = med_within |>
                 mutate(term = factor(term, levels = levels(draws_all$term))) |>
                 left_join(draws_all |> distinct(response, term),
                            by = c("response","term")),
                aes(x = med, y = term, colour = response),
                inherit.aes = FALSE, size = 2.4) +
    scale_colour_manual(values = V2_PALETTES$qualitative[1:length(unique(draws_all$response))],
                          name = NULL) +
    labs(title = "Drivers across responses — beeswarm overlay",
         subtitle = sprintf("Each color = one response | %s", RUN_LABEL),
         x = "Standardised coefficient", y = NULL) +
    theme_v2_pub(11) +
    theme(legend.position = "top")
  add_plot("driver_beeswarm_overlay", p_bee, w = 11, h = 7)
}

## --- G. effort spatial composite (双源) ----------------------------------

eff_path <- v2_file("results", "table_effort_by_grid_period_source")
if (file.exists(eff_path)) {
  eff <- read_csv(eff_path, show_col_types = FALSE) |>
    mutate(source_short = factor(source_short,
      levels = c("China_Birdwatch", "eBird_GBIF", "Combined")))
  log_max <- log10(quantile(eff$n_events, 0.99, na.rm = TRUE) + 1)
  eff_sf <- grid_sf |> inner_join(eff, by = "grid_cell")
  p_eff <- ggplot(eff_sf) +
    geom_sf(aes(fill = log10(n_events + 1)),
            colour = scales::alpha("white", 0.06), linewidth = 0.04) +
    china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
    scale_fill_v2_sequential(name = "log10(events+1)", palette = "lajolla",
                             direction = 1, limits = c(0, log_max),
                             oob = scales::squish) +
    facet_grid(source_short ~ block_label) +
    v2_china_coord() +
    theme_v2_map(9.4) +
    labs(title = "Survey effort by source × 5-year period",
         subtitle = "Rows: data source; cols: primary period")
  add_plot("effort_events_by_source_period", p_eff, w = 14, h = 7.6, is_map = TRUE)
}

## --- H. naive vs corrected richness --------------------------------------

p_naive_path <- v2_file("results",
                paste0("table_community_metrics_with_cri_", RUN_LABEL))
sv <- readRDS(file.path(P$derived_v2, "species_visit_2000_2024.rds"))
naive_block <- sv |>
  group_by(grid_cell, block_id) |>
  summarise(naive_richness = n_distinct(species), .groups = "drop") |>
  inner_join(primary_blocks, by = "block_id")
corr_block <- m_long |>
  filter(metric == "corrected_richness") |>
  select(grid_cell, block_label, corrected_richness = value_mean)
cmp <- naive_block |> inner_join(corr_block, by = c("grid_cell","block_label"))
p_cmp <- ggplot(cmp, aes(naive_richness, corrected_richness, colour = block_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.55, size = 1.2) +
  scale_colour_manual(values = V2_PALETTES$qualitative[1:5], name = "Period") +
  facet_wrap(~ block_label, ncol = 5) +
  labs(title = "Naive vs occupancy-corrected richness",
       x = "Naive species richness (observed)",
       y = "Occupancy-corrected richness (sum of psi)") +
  theme_v2_pub(11) +
  theme(legend.position = "none")
add_plot("naive_vs_corrected_richness", p_cmp, w = 14, h = 4.4)

## --- I. biotic homogenization trend --------------------------------------

homog_path <- v2_file("results",
                paste0("table_spatial_homogenization_", RUN_LABEL))
if (file.exists(homog_path)) {
  homog <- read_csv(homog_path, show_col_types = FALSE)
  p_homog <- ggplot(homog, aes(x = block_id)) +
    geom_ribbon(aes(ymin = median_pairwise_sorensen_l95,
                     ymax = median_pairwise_sorensen_u95),
                fill = "#0E5A78", alpha = 0.25) +
    geom_line(aes(y = median_pairwise_sorensen),
              colour = "#0E5A78", linewidth = 0.8) +
    geom_point(aes(y = median_pairwise_sorensen),
               colour = "#0E5A78", size = 2.6) +
    scale_x_continuous(breaks = seq_len(nrow(homog)),
                        labels = homog$block_label) +
    labs(title = "Biotic homogenization across China",
         subtitle = "Median pairwise Sorensen | occupancy-corrected",
         x = NULL, y = "Median pairwise Sorensen") +
    theme_v2_pub(11)
  add_plot("biotic_homogenization", p_homog, w = 9, h = 5.4)
}

## --- J. latitudinal gradient steepness over time -------------------------

lat_path <- v2_file("results",
                paste0("table_latitudinal_gradient_strength_", RUN_LABEL))
if (file.exists(lat_path)) {
  lat <- read_csv(lat_path, show_col_types = FALSE) |>
    mutate(block_label = factor(block_label,
                                  levels = primary_blocks$block_label))
  p_lat <- ggplot(lat, aes(block_label, slope, group = metric, colour = metric)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbar(aes(ymin = slope - 1.96*se, ymax = slope + 1.96*se),
                  width = 0.15, linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.3) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = V2_PALETTES$qualitative[1:5], guide = "none") +
    labs(title = "Latitudinal gradient steepness across primary periods",
         x = NULL, y = "Slope (per degree latitude)") +
    theme_v2_pub(11) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))
  add_plot("latitudinal_gradient", p_lat, w = 12, h = 6)
}

## --- 输出 master PPTX deck -----------------------------------------------

maps_deck <- save_pptx_deck(deck_maps,
                             paste0("v2_maps_deck_", RUN_LABEL),
                             width = 12, height = 7,
                             embed = "png", dpi = 240)
stats_deck <- save_pptx_deck(deck_stats,
                              paste0("v2_stats_deck_editable_", RUN_LABEL),
                              width = 12, height = 7,
                              embed = "vector")
message("[stage-10] maps deck (PNG embed): ", maps_deck)
message("[stage-10] stats deck (editable vector): ", stats_deck)
message(sprintf("[stage-10] Done. %d maps + %d non-map figures.",
                length(deck_maps), length(deck_stats)))
