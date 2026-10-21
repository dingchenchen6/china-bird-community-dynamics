#!/usr/bin/env Rscript

## Assemble the final research package from the reviewed 200sp and 500sp suites.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(officer)
  library(grid)
})

ROOT <- normalizePath(getwd(), mustWork = TRUE)
TOP_DIR <- file.path(ROOT, "figures_top_journal_20260602")
SP500_DIR <- file.path(ROOT, "figures_500sp_all_analysis_20260602")
RESULTS_DIR <- file.path(ROOT, "results_v3")
OUT_DIR <- file.path(ROOT, "final_research_package_20260714")
MAIN_DIR <- file.path(OUT_DIR, "main_figures")
SUPP_DIR <- file.path(OUT_DIR, "supplementary_atlases")
dir.create(MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

has_rvg <- requireNamespace("rvg", quietly = TRUE)

theme_workflow <- theme_void(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 18, colour = "#1E2A32", hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = "#52616B", hjust = 0, margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, colour = "#6B7780", hjust = 0),
    plot.margin = margin(18, 18, 14, 18)
  )

box <- function(xmin, xmax, ymin, ymax, fill, title, body) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = "#52616B", linewidth = 0.45),
    annotate("text", x = (xmin + xmax) / 2, y = ymax - 3.2, label = title,
             fontface = "bold", size = 4.2, colour = "#1E2A32"),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2 - 1.5,
             label = body, size = 3.25, lineheight = 1.08, colour = "#35424A")
  )
}

p_workflow <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 62), clip = "off") +
  box(2, 23, 26, 53, "#E8F2F8", "Data evidence",
      "7.49 million occurrence records\n1,426 recorded species\n3.15 million breeding-season records") +
  box(27, 47, 26, 53, "#EFEAF7", "Study design",
      "1,247 100-km grid cells\nFive periods, 2000-2024\nEffort and duration covariates") +
  box(51, 73, 39.5, 56, "#F6E8DD", "Primary spatial model",
      "200-species stMsPGOcc\nSpatial NNGP + temporal AR(1)\n4 chains; maximum R-hat 1.039") +
  box(51, 73, 20.5, 37, "#FBEFE7", "Breadth extension",
      "500-species tMsPGOcc\nTemporal AR(1), non-spatial\n4-chain metadata; R-hat unavailable") +
  box(77, 98, 26, 53, "#E7F2EA", "Ecological inference",
      "20 diversity indicators\nBaselga beta decomposition\nSpecies trends, traits and drivers\nNaive-versus-corrected comparison") +
  annotate("segment", x = 23.6, xend = 26.4, y = 39.5, yend = 39.5,
           linewidth = 0.75, colour = "#68757E", arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 47.6, xend = 50.4, y = 43.5, yend = 47.5,
           linewidth = 0.75, colour = "#68757E", arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 47.6, xend = 50.4, y = 35.5, yend = 28.5,
           linewidth = 0.75, colour = "#68757E", arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 73.6, xend = 76.4, y = 47.5, yend = 43.5,
           linewidth = 0.75, colour = "#68757E", arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 73.6, xend = 76.4, y = 28.5, yend = 35.5,
           linewidth = 0.75, colour = "#68757E", arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("rect", xmin = 5, xmax = 95, ymin = 4, ymax = 15.5,
           fill = "#F4F6F8", colour = "#AAB4BC", linewidth = 0.45) +
  annotate("text", x = 50, y = 12.2, label = "Core synthesis", fontface = "bold",
           size = 4.3, colour = "#1E2A32") +
  annotate("text", x = 50, y = 8.0,
           label = "Rising taxonomic and phylogenetic diversity | mild functional contraction | late nestedness | broad-habitat winners",
           size = 3.45, colour = "#35424A") +
  labs(
    title = "Inference architecture for Chinese bird community change",
    subtitle = "A spatially explicit primary analysis and a broader temporal extension support complementary claims",
    caption = "Maps from the 200-species model support spatial inference; 500-species maps are descriptive gridded posterior summaries."
  ) +
  theme_workflow

workflow_base <- file.path(MAIN_DIR, "Main_Fig01_research_design_and_inference")
ggsave(paste0(workflow_base, ".png"), p_workflow, width = 13.2, height = 7.2, dpi = 450, bg = "white")
ggsave(paste0(workflow_base, ".pdf"), p_workflow, width = 13.2, height = 7.2, device = cairo_pdf, bg = "white")

ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
loc <- ph_location(left = 0.25, top = 0.2, width = 12.8, height = 7.05)
if (has_rvg) {
  ppt <- ph_with(ppt, rvg::dml(ggobj = p_workflow), location = loc)
} else {
  ppt <- ph_with(ppt, external_img(paste0(workflow_base, ".png"), width = 12.8, height = 7.05), location = loc)
}
print(ppt, target = paste0(workflow_base, ".pptx"))

main_items <- tibble::tribble(
  ~main_order, ~short_name, ~title, ~source_base,
  1, "research_design_and_inference", "Research design and inference hierarchy", workflow_base,
  2, "spatial_multidiversity_200sp", "Spatial multidiversity across five periods (200sp primary model)", file.path(TOP_DIR, "fig01_multidiversity_timeslices_200sp_top_journal"),
  3, "spatial_trends_200sp", "Spatial fingerprints of change (200sp primary model)", file.path(TOP_DIR, "fig02_community_trends_200sp_top_journal"),
  4, "standardized_trajectories_500sp", "Standardized trajectories across 20 indicators (500sp extension)", file.path(SP500_DIR, "fig03_500sp_all20_global_trajectories_editable"),
  5, "beta_turnover_nestedness_500sp", "Turnover-to-nestedness transition", file.path(SP500_DIR, "fig07_500sp_global_baselga_turnover_nestedness_editable"),
  6, "species_winners_losers_500sp", "Species winners, stable taxa and losers", file.path(SP500_DIR, "fig11_500sp_species_trend_classes_editable"),
  7, "naive_vs_corrected", "Observation bias and occupancy correction", file.path(TOP_DIR, "fig06_naive_vs_corrected_trends_top_journal_editable"),
  8, "driver_importance", "Environmental driver importance", file.path(SP500_DIR, "fig12_500sp_rf_driver_importance_editable")
)

for (i in seq_len(nrow(main_items))) {
  if (i == 1) next
  for (ext in c("png", "pdf", "pptx")) {
    src <- paste0(main_items$source_base[i], ".", ext)
    dst <- file.path(MAIN_DIR, sprintf("Main_Fig%02d_%s.%s", main_items$main_order[i], main_items$short_name[i], ext))
    if (!file.exists(src)) stop("Missing main figure source: ", src)
    file.copy(src, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
}

atlas <- read_pptx()
for (i in seq_len(nrow(main_items))) {
  png <- file.path(MAIN_DIR, sprintf("Main_Fig%02d_%s.png", main_items$main_order[i], main_items$short_name[i]))
  atlas <- add_slide(atlas, layout = "Blank", master = "Office Theme")
  atlas <- ph_with(atlas, main_items$title[i], ph_location(left = 0.35, top = 0.12, width = 12.4, height = 0.38))
  atlas <- ph_with(atlas, external_img(png, width = 12.65, height = 6.65),
                   ph_location(left = 0.34, top = 0.55, width = 12.65, height = 6.65))
}
print(atlas, target = file.path(OUT_DIR, "MAIN_FIGURE_ATLAS_20260714.pptx"))

file.copy(file.path(SP500_DIR, "bird_community_500sp_all_analysis_atlas_20260602.pptx"),
          file.path(SUPP_DIR, "SUPPLEMENTARY_500sp_ALL_ANALYSES.pptx"), overwrite = TRUE)
file.copy(file.path(TOP_DIR, "bird_community_top_journal_figure_atlas_20260602.pptx"),
          file.path(SUPP_DIR, "SUPPLEMENTARY_200sp_500sp_SYNTHESIS.pptx"), overwrite = TRUE)

m500 <- read_csv(file.path(SP500_DIR, "figure_manifest_500sp_all_analysis_20260602.csv"), show_col_types = FALSE) |>
  mutate(suite = "500sp temporal extension")
mtop <- read_csv(file.path(TOP_DIR, "figure_manifest_20260602.csv"), show_col_types = FALSE) |>
  mutate(suite = "200sp spatial primary + cross-model synthesis")

main_lookup <- c(
  fig01_multidiversity_timeslices_200sp_top_journal = 2,
  fig02_community_trends_200sp_top_journal = 3,
  fig03_500sp_all20_global_trajectories_editable = 4,
  fig07_500sp_global_baselga_turnover_nestedness_editable = 5,
  fig11_500sp_species_trend_classes_editable = 6,
  fig06_naive_vs_corrected_trends_top_journal_editable = 7,
  fig12_500sp_rf_driver_importance_editable = 8
)
optional_main <- c(
  "fig04_500sp_all20_endpoint_change_editable",
  "fig15_500sp_trait_environment_mechanisms_editable"
)

catalog <- bind_rows(mtop, m500) |>
  mutate(
    main_order = unname(main_lookup[figure]),
    recommended_role = case_when(
      !is.na(main_order) ~ "Main figure",
      figure %in% optional_main ~ "Optional main / high-priority supplement",
      TRUE ~ "Supplementary figure"
    ),
    inference_scope = case_when(
      grepl("200sp", figure) & !grepl("500sp", figure) ~ "200sp spatial primary inference",
      grepl("500sp", figure) & grepl("timeslice|trend_maps|community_trends|beta_maps|dynamics_maps|hotspot_maps|bias_maps", figure, ignore.case = TRUE) ~ "500sp temporal extension; descriptive gridded summaries",
      grepl("500sp", figure) ~ "500sp temporal breadth extension",
      TRUE ~ "Cross-model synthesis"
    ),
    caution = case_when(
      grepl("500sp", figure) & grepl("timeslice|trend_maps|community_trends|beta_maps|dynamics_maps|hotspot_maps|bias_maps", figure, ignore.case = TRUE) ~
        "Do not describe as a 500sp spatial NNGP model.",
      grepl("trait|environment", figure, ignore.case = TRUE) ~
        "Mechanism analysis is complete-case (n=200 species).",
      TRUE ~ ""
    )
  ) |>
  select(main_order, recommended_role, suite, figure, inference_scope, pptx_editable, caution, png, pdf, pptx) |>
  arrange(is.na(main_order), main_order, suite, figure)

catalog <- bind_rows(
  tibble(
    main_order = 1,
    recommended_role = "Main figure",
    suite = "Final synthesis",
    figure = "Main_Fig01_research_design_and_inference",
    inference_scope = "Study design and inference hierarchy",
    pptx_editable = ifelse(has_rvg, "editable_drawingml", "embedded_png"),
    caution = "",
    png = paste0(workflow_base, ".png"),
    pdf = paste0(workflow_base, ".pdf"),
    pptx = paste0(workflow_base, ".pptx")
  ),
  catalog
)
write_csv(catalog, file.path(OUT_DIR, "FINAL_FIGURE_CATALOG_20260714.csv"))

read_metrics <- function(label) {
  read_csv(file.path(RESULTS_DIR, paste0("table_community_metrics_with_cri_v3_full_", label, ".csv")), show_col_types = FALSE)
}

core_metrics <- c("corrected_richness", "shannon", "inv_simpson", "trait_volume", "rao_q", "feve_fund", "pd_prob_mctavish", "mpd_prob_mctavish")
metric_results <- bind_rows(
  read_metrics("200sp_ar1_spatial_extended") |> mutate(dataset = "200sp spatial primary"),
  read_metrics("500sp_ar1_temporal_extended") |> mutate(dataset = "500sp temporal extension")
) |>
  filter(metric %in% core_metrics, period %in% c("P1", "P5")) |>
  group_by(dataset, metric, period) |>
  summarise(value = mean(value_mean, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = value) |>
  mutate(absolute_change = P5 - P1, percent_change = 100 * (P5 / P1 - 1))

classes <- read_csv(file.path(RESULTS_DIR, "table_species_trend_classify_v3_full_500sp_ar1_temporal_extended.csv"), show_col_types = FALSE) |>
  count(trend_class, name = "estimate") |>
  transmute(section = "species trends", dataset = "500sp temporal extension", metric = trend_class,
            estimate = as.numeric(estimate), unit = "species", note = "")

beta <- read_csv(file.path(RESULTS_DIR, "table_baselga_global_v3_full_500sp_ar1_temporal_extended.csv"), show_col_types = FALSE) |>
  transmute(section = "beta diversity", dataset = "500sp temporal extension", metric = paste0(period_pair, " turnover share"),
            estimate = prop_turnover_mean, unit = "proportion", note = paste0("95% CrI ", round(prop_turnover_q025, 3), "-", round(prop_turnover_q975, 3)))

naive <- read_csv(file.path(RESULTS_DIR, "table_naive_vs_corrected_v3_full_500sp_ar1_temporal_extended.csv"), show_col_types = FALSE)
naive_rows <- tibble(
  section = "detection correction",
  dataset = "500sp temporal extension",
  metric = c("direction flips", "median absolute trend difference", "mean corrected-minus-naive trend difference"),
  estimate = c(sum(naive$direction_flipped, na.rm = TRUE), median(abs(naive$trend_diff), na.rm = TRUE), mean(naive$trend_diff, na.rm = TRUE)),
  unit = c("grid cells", "species per period", "species per period"),
  note = c(paste0("out of ", nrow(naive)), "", "")
)

survey <- read_csv(file.path(RESULTS_DIR, "table_02_survey_summary_v3.csv"), show_col_types = FALSE)
survey_rows <- tibble(
  section = "study scope",
  dataset = "integrated occurrence database",
  metric = c("recorded species", "analysed grid cells", "primary periods", "total events", "breeding-season events"),
  estimate = c(survey$n_species[1], survey$n_sites[1], survey$n_periods[1], survey$n_events_total[1], survey$n_events_breeding[1]),
  unit = c("species", "100-km grid cells", "periods", "events", "events"),
  note = c("", "", "2000-2024", "", "")
)

conv200 <- read_csv(file.path(RESULTS_DIR, "table_convergence_diagnostics_v3_full_200sp_ar1_spatial_4chain.csv"), show_col_types = FALSE)
conv500 <- read_csv(file.path(RESULTS_DIR, "table_convergence_diagnostics_v3_full_500sp_ar1_temporal_4chain.csv"), show_col_types = FALSE)
diagnostic_rows <- tibble(
  section = "model diagnostics",
  dataset = c("200sp spatial primary", "200sp spatial primary", "500sp temporal extension", "500sp temporal extension"),
  metric = c("maximum R-hat", "minimum ESS", "available R-hat values", "minimum ESS"),
  estimate = c(max(conv200$rhat, na.rm = TRUE), min(conv200$ess, na.rm = TRUE), sum(!is.na(conv500$rhat)), min(conv500$ess, na.rm = TRUE)),
  unit = c("R-hat", "ESS", "parameters", "ESS"),
  note = c("4-chain diagnostic table", "4-chain diagnostic table", "R-hat is unavailable for all 500sp diagnostic rows", "4-chain metadata")
)

rf_rows <- read_csv(file.path(RESULTS_DIR, "table_rf_importance_group_summary_trend.csv"), show_col_types = FALSE) |>
  transmute(section = "driver analysis", dataset = "posterior random forest", metric = label,
            estimate = mean, unit = "group importance", note = paste0("95% interval ", round(q_lo, 4), "-", round(q_hi, 4)))

trait_rows <- read_csv(file.path(RESULTS_DIR, "table_trait_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"), show_col_types = FALSE) |>
  transmute(section = "trait mechanism", dataset = "500sp-labelled complete-case analysis", metric = trait,
            estimate = rho_mean, unit = "Spearman rho", note = paste0("95% interval ", round(rho_q025, 4), "-", round(rho_q975, 4), "; complete-case n=200"))

env_rows <- read_csv(file.path(RESULTS_DIR, "table_env_trend_correlation_v3_full_500sp_ar1_temporal_extended.csv"), show_col_types = FALSE) |>
  transmute(section = "environment association", dataset = "500sp-labelled complete-case analysis", metric = env_var,
            estimate = spearman_rho, unit = "Spearman rho", note = paste0("complete-case n=", n_species))

vp_names <- c(
  "[a] = X1 | X2+X3+X4" = "Climate change pure fraction",
  "[b] = X2 | X1+X3+X4" = "Land-use change pure fraction",
  "[c] = X3 | X1+X2+X4" = "Human pressure pure fraction",
  "[d] = X4 | X1+X2+X3" = "Spatial baseline pure fraction",
  "[p] = Residuals" = "Residual fraction"
)
varpart_rows <- read_csv(file.path(RESULTS_DIR, "table_varpart_corrected_richness_v3_full_500sp_ar1_temporal_extended.csv"), show_col_types = FALSE) |>
  filter(fraction %in% names(vp_names)) |>
  transmute(section = "variance partitioning", dataset = "500sp corrected-richness trend", metric = unname(vp_names[fraction]),
            estimate = adj_r2, unit = "adjusted R2", note = "")

numeric_summary <- bind_rows(
  survey_rows,
  diagnostic_rows,
  metric_results |>
    transmute(section = "multidiversity", dataset, metric,
              estimate = absolute_change, unit = "P5 minus P1", note = paste0("P1=", signif(P1, 5), "; P5=", signif(P5, 5), "; change=", round(percent_change, 2), "%")),
  classes,
  beta,
  naive_rows,
  rf_rows,
  trait_rows,
  env_rows,
  varpart_rows
)
write_csv(numeric_summary, file.path(OUT_DIR, "FINAL_RESULTS_NUMERIC_SUMMARY_20260714.csv"))

message("[23] Final package: ", OUT_DIR)
message("[23] Main figures: ", nrow(main_items))
message("[23] Catalog rows: ", nrow(catalog))
