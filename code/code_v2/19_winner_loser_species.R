#!/usr/bin/env Rscript
## 19_winner_loser_species.R
##
## 基于 stage-4 物种 β 后验 + 200 thinned psi draws，识别"赢家/输家"物种。
## 输出：3 张组合图（森林 + 雨林 + 山脊）+ 1 张性状-趋势散点。
##
## Output:
##   results_v3/table_species_winners_losers_v3.csv
##   figures_v3/fig_v3_winner_loser_forest_v3.{png,pdf}
##   figures_v3/fig_v3_winner_loser_raincloud_by_trait_v3.{png,pdf}
##   figures_v3/fig_v3_winner_loser_ridge_by_family_v3.{png,pdf}
##   figures_v3/fig_v3_trait_vs_trend_scatter_v3.{png,pdf}
##   figures_v3/v3_winner_loser_deck.pptx

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(forcats); library(ggplot2); library(patchwork)
  library(ggridges); library(ggdist); library(ggbeeswarm)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_plots_advanced.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")

message(sprintf("[stage-19] winner/loser species analysis for %s", RUN_LABEL))

## --- 1. 读 species coefficient 表 + trait ext ---------------------------

coef_path <- v2_file("results",
                paste0("table_species_coefficients_", RUN_LABEL))
species_coef <- read_csv(coef_path, show_col_types = FALSE)
trait_ext <- readRDS(file.path(P$derived_v2, "trait_extended.rds"))

year_trend <- species_coef |>
  filter(covariate == "year_scaled") |>
  rename(year_slope = mean, year_l95 = l95, year_u95 = u95, year_sig = sig95) |>
  select(species, year_slope, year_l95, year_u95, year_sig)

## --- 2. 合并性状 -------------------------------------------------------

ws <- year_trend |>
  left_join(trait_ext, by = "species") |>
  mutate(
    status = case_when(
      year_sig & year_slope > 0 ~ "Winner",
      year_sig & year_slope < 0 ~ "Loser",
      TRUE                       ~ "Stable"
    ),
    status = factor(status, levels = c("Winner", "Stable", "Loser"))
  )
n_win <- sum(ws$status == "Winner")
n_los <- sum(ws$status == "Loser")
n_stab <- sum(ws$status == "Stable")
message(sprintf("  Winners=%d  Stable=%d  Losers=%d", n_win, n_stab, n_los))
write_csv(ws, v3_file("results", "table_species_winners_losers_v3"))

## --- 3. 森林图：top 30 winners + top 30 losers -------------------------

top_w <- ws |> filter(status == "Winner") |> arrange(desc(year_slope)) |> head(30)
top_l <- ws |> filter(status == "Loser")  |> arrange(year_slope) |> head(30)
top_show <- bind_rows(top_w, top_l) |>
  mutate(species = forcats::fct_reorder(species, year_slope))

p_forest <- ggplot(top_show, aes(year_slope, species,
                                   colour = status, fill = status)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = year_l95, xmax = year_u95),
                  height = 0, linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(Winner = "#8B2E1E", Loser = "#0E5A78",
                                    Stable = "grey55"), guide = "none") +
  scale_fill_manual(values = c(Winner = "#8B2E1E", Loser = "#0E5A78",
                                  Stable = "grey55"), guide = "none") +
  labs(title = "Top 30 winners and top 30 losers",
       subtitle = sprintf("Per-species occupancy time slope (year_scaled) | run=%s",
                          RUN_LABEL),
       caption  = "Posterior mean +/- 95% CRI. Winners (red) = significant increase; Losers (blue) = significant decrease.",
       x = "Year-trend coefficient (logit-scale)",
       y = NULL) +
  theme_v2_pub(10.5) +
  theme(axis.text.y = element_text(size = 7.5))
save_dual_v3(p_forest, "fig_v3_winner_loser_forest_v3", width = 9.4, height = 11)

## --- 4. 雨林图：按性状 group ------------------------------------------

group_cols <- c(Trophic.Niche      = "Trophic niche",
                Habitat            = "AVONET Habitat",
                Primary.Lifestyle  = "Primary lifestyle")
group_cols <- group_cols[names(group_cols) %in% names(ws)]

make_group_raincloud <- function(g_col, g_label) {
  d <- ws |> filter(!is.na(.data[[g_col]])) |>
    mutate(grp = forcats::fct_reorder(.data[[g_col]], year_slope, .fun = median))
  ggplot(d, aes(year_slope, grp, fill = status, colour = status)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggdist::stat_halfeye(thickness = 0.55, .width = c(0.5, 0.95),
                          slab_alpha = 0.45, slab_colour = NA,
                          interval_colour = "grey30",
                          side = "right", justification = -0.15,
                          show.legend = FALSE) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE, alpha = 0.45, size = 0.7,
                                  width = 0.18, show.legend = FALSE) +
    scale_colour_manual(values = c(Winner = "#8B2E1E", Loser = "#0E5A78",
                                      Stable = "grey55"), name = NULL) +
    scale_fill_manual(values = c(Winner = "#8B2E1E", Loser = "#0E5A78",
                                    Stable = "grey55"), name = NULL) +
    labs(title = g_label, x = "Year-trend coefficient", y = NULL) +
    theme_v2_pub(10) +
    theme(legend.position = "top")
}
if (length(group_cols) > 0) {
  plots_g <- lapply(names(group_cols),
                     function(g) make_group_raincloud(g, group_cols[[g]]))
  p_rain <- patchwork::wrap_plots(plots_g, ncol = 1) +
    patchwork::plot_annotation(
      title = "Winner/loser distributions by trait group",
      subtitle = sprintf("AVONET trait categories | run=%s", RUN_LABEL),
      theme = theme_v2_pub(11))
  save_dual_v3(p_rain, "fig_v3_winner_loser_raincloud_by_trait_v3",
            width = 11, height = 11)
} else p_rain <- NULL

## --- 5. 山脊图：按 family（前 15） ------------------------------------

p_ridge <- NULL
if ("family_lat" %in% names(ws)) {
  fam_top <- ws |> filter(!is.na(family_lat)) |>
    count(family_lat, sort = TRUE) |>
    filter(n >= 4) |> head(15) |> pull(family_lat)
  d_fam <- ws |> filter(family_lat %in% fam_top) |>
    mutate(family_lat = forcats::fct_reorder(family_lat, year_slope, .fun = median))
  p_ridge <- ggplot(d_fam, aes(year_slope, family_lat, fill = after_stat(x))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggridges::geom_density_ridges_gradient(scale = 2.2, rel_min_height = 0.008,
                                            colour = "white", linewidth = 0.2,
                                            alpha = 0.95) +
    scale_fill_gradient2(low = "#0E5A78", mid = "#F2E8D8", high = "#8B2E1E",
                          midpoint = 0, guide = "none") +
    labs(title = "Year-trend distributions by family",
         subtitle = sprintf("Top 15 families (>=4 species each) | run=%s", RUN_LABEL),
         x = "Year-trend coefficient", y = NULL) +
    theme_v2_pub(10.5)
  save_dual_v3(p_ridge, "fig_v3_winner_loser_ridge_by_family_v3",
            width = 10, height = 8.4)
}

## --- 6. 性状-趋势散点 -----------------------------------------------

scatter_targets <- intersect(
  c("body_mass_g", "clutch_size", "longevity_y", "avonet_hwi",
    "avonet_range_size", "Habitat.Density", "habitat_breadth_data",
    "diet_specialization", "habitat_openness_avonet"),
  names(ws))

if (length(scatter_targets) > 0) {
  scatter_long <- ws |>
    select(species, year_slope, year_sig, status, all_of(scatter_targets)) |>
    pivot_longer(all_of(scatter_targets), names_to = "trait",
                  values_to = "trait_value") |>
    filter(is.finite(trait_value)) |>
    mutate(trait_label = recode(trait,
      body_mass_g = "log10 body mass (g)",
      clutch_size = "log10 clutch size",
      longevity_y = "log10 longevity (y)",
      avonet_hwi  = "AVONET HWI",
      avonet_range_size = "log10 range size",
      Habitat.Density = "AVONET habitat density",
      habitat_breadth_data = "Habitat breadth (Shannon)",
      diet_specialization = "Diet specialisation",
      habitat_openness_avonet = "Habitat openness"
    ))
  p_sc <- ggplot(scatter_long, aes(trait_value, year_slope)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    geom_point(aes(colour = status), alpha = 0.55, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, colour = "#2A2A2A",
                fill = "#0E5A78", alpha = 0.12, linewidth = 0.6) +
    facet_wrap(~ trait_label, scales = "free_x", ncol = 3) +
    scale_colour_manual(values = c(Winner = "#8B2E1E", Loser = "#0E5A78",
                                      Stable = "grey55"), name = NULL) +
    labs(title = "Species traits vs occupancy time trend",
         subtitle = sprintf("Each point = one species | run=%s", RUN_LABEL),
         x = NULL, y = "Year-trend coefficient") +
    theme_v2_pub(10.5) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = 9))
  save_dual_v3(p_sc, "fig_v3_trait_vs_trend_scatter_v3",
            width = 12, height = 9)
} else p_sc <- NULL

## --- 7. PPTX 可编辑（非地图图） --------------------------------------

if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE)) {
  deck_plots <- list(
    "Top winners & losers (forest plot)" = p_forest)
  if (!is.null(p_rain))  deck_plots[["Raincloud by trait group"]] <- p_rain
  if (!is.null(p_ridge)) deck_plots[["Ridgeline by family"]]       <- p_ridge
  if (!is.null(p_sc))    deck_plots[["Trait vs year-trend"]]       <- p_sc
  doc <- officer::read_pptx()
  for (nm in names(deck_plots)) {
    doc <- officer::add_slide(doc, layout = "Title and Content",
                               master = "Office Theme") |>
      officer::ph_with(value = nm,
                        location = officer::ph_location_type(type = "title")) |>
      officer::ph_with(value = rvg::dml(ggobj = deck_plots[[nm]]),
                        location = officer::ph_location(left = 0.3, top = 1.1,
                                                        width = 12, height = 7))
  }
  print(doc, target = v3_file("figure", "v3_winner_loser_deck", "pptx"))
}

message("[stage-19] done.")
