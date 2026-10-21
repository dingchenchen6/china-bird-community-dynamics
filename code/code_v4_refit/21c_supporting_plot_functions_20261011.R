#!/usr/bin/env Rscript
# ============================================================
# 21c_supporting_plot_functions_20261011.R  (sourced by 21b_supporting_figures_native_pptx_20261011.R)
#
# Scientific question / 科学问题:
#   支撑材料图 fig05–fig09 的作图函数（物种趋势类别、朴素与校正趋势、Baselga 分解、
#   随机森林重要性、性状与环境关联）。
#   Plot functions for supporting figures fig05-fig09.
#
# Objective / 分析目标:
#   从 21_make_top_journal_figures.R 第 573–798 行复制；只有两处标注 [2026-10-11 改动] 的呈现修正：
#   fig05 小段标签移到柱顶（原白字不可读），fig08 "Land use change" 颜色映射（原为缺省灰）。
#   fig06、fig07、fig09 的 PDF/PNG 与 2026-07-14 投稿版逐像素一致。
#   Copied from script 21 (lines 573-798) with two marked presentation fixes (fig05 labels, fig08 colour).
#
# Input data / 输入数据: 由 21b 读入的表 / tables read by 21b
# Main workflow / 主要流程: 定义 make_species_class_plot, make_naive_corrected_plot,
#   make_beta_global_plot, make_rf_plot, make_trait_env_plot / define the five plot functions
# Expected output / 预期输出: 函数，不写文件 / functions only
# Key assumptions / 关键假设: 依赖 21b 中定义的调色板、主题与标签函数 / uses palettes and helpers from 21b
# Main packages / 主要包: dplyr, tidyr, ggplot2, patchwork, scales, forcats, stringr
# Output directory / 输出路径: 无 / none
# ============================================================

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
    # [2026-10-11 改动 / change] 占比 < 6% 的段太窄：白字溢出到相邻段与柱顶白底上（"2"、"15" 不可见）。
    # 这些段的标签改为柱顶上方的深色单行文字，并给 y 轴顶部留 6% 空间；其余标签不变。
    # Labels of segments < 6% (white text was unreadable) move above the bar in dark text.
    geom_text(aes(label = ifelse(prop >= 0.06, label, "")), position = position_stack(vjust = 0.5), size = 2.5, lineheight = 0.92, colour = "white") +
    geom_text(data = dplyr::filter(class_counts, prop < 0.06),
              aes(x = dataset, y = 1.01, label = paste0(n, " (", percent(prop, accuracy = 0.1), ")")),
              inherit.aes = FALSE, vjust = 0, size = 2.5, colour = "#23313A") +
    scale_fill_manual(values = pal_classes, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.06))) +
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
  # [2026-10-11 改动 / change] 表中标签 "Land use change" 与调色板键 "Land-use change" 不一致，
  # 该柱原被画成缺省灰色；统一写法后与左图同为橙色（同主文 Fig 6c 的处理）。
  # Harmonize the label so the bar takes the Land-use colour (as in main Fig 6c).
  group_d$group_label <- dplyr::recode(group_d$group_label, `Land use change` = "Land-use change")
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

