#!/usr/bin/env Rscript
# ============================================================
# 24d_figures4to6_20261011.R  (sourced by 24b_main_figures_standard_map_20261011.R)
#
# Scientific question / 科学问题:
#   Fig 4：多维多样性是否解耦；Fig 5：β 重组、赢家物种与性状关联；
#   Fig 6：探测校正如何改变推断，以及驱动因子的关联强度。
#   Fig 4 multidiversity decoupling; Fig 5 reorganization and winners; Fig 6 detection bias and drivers.
#
# Objective / 分析目标:
#   数据处理与作图代码逐字沿用 24_make_nature_submission_figures.R（第 373–568 行），
#   唯一改动：Fig 6a 的 "Direction flips" 注释加白底，避免被 1:1 虚线穿过（纯呈现）。
#   Code copied verbatim from script 24; the only change is a white backing for the Fig 6a note.
#
# Input data / 输入数据: 24b 中读入的 metrics_200/500, beta_200/500, class_500, trait_corr,
#   naive_500, rf_group, varpart_500
# Main workflow / 主要流程: 汇总 → 作图 → 写源数据 → 导出 / summarise -> plot -> source data -> export
# Expected output / 预期输出: main/Figure_4..6_*.{pdf,png,tiff,pptx}; source_data/SourceData_Figure_4..6.csv
# Key assumptions / 关键假设: 源数据须与 Shared/Figures/source_data 逐字节一致。
# Main packages / 主要包: dplyr, ggplot2, patchwork, scales
# Output directory / 输出路径: figures_community_submission_20261011/{main,source_data}
# ============================================================

# ── Figure 4 (verbatim) ─────────────────────────────────────────────
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
    metric_panel = recode(metric, corrected_richness = "Richness", shannon = "Shannon",
                          inv_simpson = "Inverse Simpson", pd_prob_mctavish = "PD", mpd_prob_mctavish = "MPD",
                          trait_volume = "Functional volume", rao_q = "Rao's Q", feve_fund = "Functional evenness"),
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
write_csv(fig4_data |> select(dataset, metric, group, period, value, lo, hi, spatial_sd, standardized_change,
                              standardized_lo, standardized_hi, percent_change),
          file.path(SRC_DIR, "SourceData_Figure_4.csv"))
save_main(fig4, "Figure_4_multidiversity_decoupling", 183, 160)

# ── Figure 5 (verbatim) ─────────────────────────────────────────────
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
write_csv(bind_rows(
  fig5_beta |> transmute(panel = "a_beta", dataset, item = as.character(period_plot), estimate = prop_turnover_mean,
                        lower = prop_turnover_q025, upper = prop_turnover_q975, n = n_grids),
  fig5_class |> transmute(panel = "b_classes", dataset = "500-species temporal", item = as.character(class), estimate = prop,
                         lower = NA_real_, upper = NA_real_, n = n),
  top_species |> transmute(panel = "c_species", dataset = "500-species temporal", item = species, estimate = mean,
                          lower = q025, upper = q975, n = NA_integer_),
  fig5_trait |> transmute(panel = "d_traits", dataset = "complete-case species", item = as.character(trait_plot), estimate = rho_mean,
                         lower = rho_q025, upper = rho_q975, n = 200L)
), file.path(SRC_DIR, "SourceData_Figure_5.csv"))
save_main(fig5, "Figure_5_reorganization_and_winners", 183, 142)

# ── Figure 6 (verbatim except the white-backed note in panel a) ───────
flip_n <- sum(naive_500$direction_flipped, na.rm = TRUE)
p6a <- ggplot(naive_500, aes(naive_trend, corrected_trend)) +
  geom_abline(slope = 1, intercept = 0, colour = "#7F8990", linewidth = 0.35, linetype = "22") +
  geom_point(aes(colour = direction_flipped), shape = 16, size = 0.85, alpha = 0.48) +
  geom_smooth(method = "lm", se = FALSE, colour = COL[["ink"]], linewidth = 0.55) +
  scale_colour_manual(values = c(`FALSE` = COL[["neutral"]], `TRUE` = COL[["orange"]]), guide = "none") +
  # 改动：annotate("text") → 白底 label，1:1 虚线不再穿过文字 / white backing so the 1:1 line no longer crosses the text
  annotate("label", x = -Inf, y = Inf, label = paste0("Direction flips: ", flip_n, " / ", nrow(naive_500)),
           hjust = -0.05, vjust = 1.35, size = 2.15, fill = "white", border.colour = NA,
           label.padding = unit(0.4, "mm"), family = "Helvetica") +
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
vp_labels <- c("[a] = X1 | X2+X3+X4" = "Climate change", "[b] = X2 | X1+X3+X4" = "Land-use change",
               "[c] = X3 | X1+X2+X4" = "Human pressure change", "[d] = X4 | X1+X2+X3" = "Spatial baseline")
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
write_csv(bind_rows(
  naive_500 |> transmute(panel = "a_b_detection", item = as.character(grid_cell), estimate = corrected_trend,
                         comparator = naive_trend, difference = trend_diff, lower = NA_real_, upper = NA_real_),
  rf_plot |> transmute(panel = "c_random_forest", item = as.character(group_plot), estimate = mean,
                      comparator = NA_real_, difference = NA_real_, lower = q_lo, upper = q_hi),
  vp_plot |> transmute(panel = "d_variance_partition", item = as.character(group_plot), estimate = adj_r2,
                      comparator = NA_real_, difference = NA_real_, lower = NA_real_, upper = NA_real_)
), file.path(SRC_DIR, "SourceData_Figure_6.csv"))
save_main(fig6, "Figure_6_detection_bias_and_drivers", 183, 140)
