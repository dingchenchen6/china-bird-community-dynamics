## utils_plots_advanced.R
## 高级出图工具：替换简陋的 forest plot 为雨林图 / 蜂群图。
## 适用于：
##   - brms 驱动模型的系数 posterior（每条 cmdstanr 抽出来的 draws）
##   - tMsPGOcc beta.comm / alpha.comm 多链 posterior
##   - 任何 (term, draw, value) long table

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(ggdist)
  library(ggbeeswarm)
})

## ---- 1. 通用：从 brms fit 抽 fixed-effect posterior draws -----------------

brms_fixef_draws <- function(brms_fit, drop_intercept = TRUE) {
  # 兼容多版本 brms
  d <- tryCatch(
    posterior::as_draws_df(brms_fit) |> as.data.frame(),
    error = function(e) {
      m <- as.matrix(brms_fit)
      as.data.frame(m)
    }
  )
  cols <- grep("^b_", names(d), value = TRUE)
  if (drop_intercept) cols <- setdiff(cols, "b_Intercept")
  long <- tidyr::pivot_longer(d[, cols, drop = FALSE],
                               cols = everything(),
                               names_to = "term",
                               values_to = "value")
  long$term <- sub("^b_", "", long$term)
  long
}

## ---- 2. 雨林图（halfeye + jitter + interval） -----------------------------
##  推荐：稳健概览每个系数的整个后验分布

raincloud_posterior <- function(draws_long,
                                 title = NULL, subtitle = NULL,
                                 caption = NULL,
                                 x_lab = "Standardised coefficient",
                                 col_neg = "#0E5A78",
                                 col_pos = "#8B2E1E",
                                 col_zero = "#7C7C7C",
                                 ref_line = 0,
                                 max_points = 600) {
  med_tbl <- draws_long |>
    group_by(term) |>
    summarise(median = median(value),
              l95 = quantile(value, 0.025),
              u95 = quantile(value, 0.975),
              .groups = "drop") |>
    mutate(sign = case_when(l95 > 0 ~ "+",
                             u95 < 0 ~ "-",
                             TRUE ~ "0")) |>
    arrange(median) |>
    mutate(term = factor(term, levels = term))
  draws_long <- draws_long |>
    mutate(term = factor(term, levels = levels(med_tbl$term))) |>
    left_join(med_tbl |> select(term, sign), by = "term")
  # 抽样限点防过度拥挤
  if (nrow(draws_long) > max_points * length(levels(med_tbl$term))) {
    draws_long <- draws_long |>
      group_by(term) |>
      slice_sample(n = max_points) |>
      ungroup()
  }
  ggplot(draws_long, aes(value, term, colour = sign, fill = sign)) +
    geom_vline(xintercept = ref_line, linetype = 2, colour = "grey55",
               linewidth = 0.4) +
    ggdist::stat_halfeye(
      thickness = 0.55, .width = c(0.5, 0.95),
      slab_colour = NA, slab_size = 0,
      slab_alpha = 0.42,
      side = "right", justification = -0.18,
      adjust = 1.1, density = "bounded",
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      groupOnX = FALSE, alpha = 0.18, size = 0.45,
      width = 0.18, show.legend = FALSE
    ) +
    geom_point(
      data = med_tbl, aes(x = median, y = term, colour = sign),
      inherit.aes = FALSE, size = 2.6
    ) +
    geom_errorbar(
      data = med_tbl, aes(xmin = l95, xmax = u95, y = term, colour = sign),
      inherit.aes = FALSE, width = 0, linewidth = 0.55, alpha = 0.9
    ) +
    scale_colour_manual(values = c(`-` = col_neg, `+` = col_pos,
                                    `0` = col_zero), guide = "none") +
    scale_fill_manual(values = c(`-` = col_neg, `+` = col_pos,
                                  `0` = col_zero), guide = "none") +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = x_lab, y = NULL)
}

## ---- 3. 蜂群图（紧凑版，对系数密集场景） ---------------------------------

beeswarm_posterior <- function(draws_long,
                                title = NULL, subtitle = NULL,
                                caption = NULL,
                                x_lab = "Standardised coefficient",
                                ref_line = 0,
                                max_points = 400) {
  med_tbl <- draws_long |>
    group_by(term) |>
    summarise(median = median(value),
              l95 = quantile(value, 0.025),
              u95 = quantile(value, 0.975),
              .groups = "drop") |>
    arrange(median) |>
    mutate(term = factor(term, levels = term),
           sign = case_when(l95 > 0 ~ "+", u95 < 0 ~ "-", TRUE ~ "0"))
  draws_long <- draws_long |>
    mutate(term = factor(term, levels = levels(med_tbl$term))) |>
    left_join(med_tbl |> select(term, sign), by = "term")
  if (nrow(draws_long) > max_points * length(levels(med_tbl$term))) {
    draws_long <- draws_long |>
      group_by(term) |>
      slice_sample(n = max_points) |>
      ungroup()
  }
  ggplot(draws_long, aes(value, term, colour = sign)) +
    geom_vline(xintercept = ref_line, linetype = 2, colour = "grey55",
               linewidth = 0.4) +
    ggbeeswarm::geom_quasirandom(
      groupOnX = FALSE, alpha = 0.32, size = 0.5, width = 0.32,
      show.legend = FALSE
    ) +
    geom_errorbarh(
      data = med_tbl, aes(xmin = l95, xmax = u95, y = term, colour = sign),
      inherit.aes = FALSE, height = 0, linewidth = 0.6, alpha = 0.95
    ) +
    geom_point(
      data = med_tbl, aes(x = median, y = term, colour = sign),
      inherit.aes = FALSE, size = 2.6, fill = "white", shape = 21,
      stroke = 0.7
    ) +
    scale_colour_manual(values = c(`-` = "#0E5A78", `+` = "#8B2E1E",
                                    `0` = "#7C7C7C"), guide = "none") +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = x_lab, y = NULL)
}

## ---- 4. 从一组 brms 模型一次出"多响应面板雨林图" --------------------------

multi_response_raincloud <- function(brms_list, response_labels = NULL,
                                      ncol = 2,
                                      title = "Drivers across responses",
                                      subtitle = "Posterior densities; medians + 95% CRI",
                                      caption = "Coefficients sorted within each response by posterior median.") {
  if (is.null(response_labels)) response_labels <- names(brms_list)
  draws_all <- purrr::imap_dfr(brms_list, function(fit, nm) {
    brms_fixef_draws(fit) |>
      mutate(response = response_labels[[nm]] %||% nm)
  })
  # 每个响应内部按 median 排序
  draws_all <- draws_all |>
    group_by(response, term) |>
    mutate(med = median(value)) |>
    ungroup() |>
    mutate(term = forcats::fct_reorder(term, med))
  ggplot(draws_all, aes(value, term)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggdist::stat_halfeye(.width = c(0.5, 0.95), thickness = 0.55,
                         slab_colour = NA, slab_alpha = 0.55,
                         interval_colour = "grey25", point_size = 1.6,
                         fill = "#0E5A78") +
    facet_wrap(~ response, scales = "free", ncol = ncol) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = "Standardised coefficient (95% CRI)", y = NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- 5. 山脊图（ggridges）：每条系数一条 ridge density --------------------

ridgeline_posterior <- function(draws_long,
                                 title = NULL, subtitle = NULL,
                                 caption = NULL,
                                 x_lab = "Standardised coefficient",
                                 ref_line = 0,
                                 scale_value = 2.6,
                                 fill_low = "#0E5A78",
                                 fill_high = "#8B2E1E") {
  if (!requireNamespace("ggridges", quietly = TRUE))
    stop("ggridges required for ridgeline_posterior()")
  med_tbl <- draws_long |>
    group_by(term) |>
    summarise(median = median(value),
              l95 = quantile(value, 0.025),
              u95 = quantile(value, 0.975),
              .groups = "drop") |>
    arrange(median) |>
    mutate(term = factor(term, levels = term))
  draws_long <- draws_long |>
    mutate(term = factor(term, levels = levels(med_tbl$term)))
  ggplot(draws_long, aes(value, term, fill = after_stat(x))) +
    geom_vline(xintercept = ref_line, linetype = 2, colour = "grey55",
               linewidth = 0.4) +
    ggridges::geom_density_ridges_gradient(
      scale = scale_value, rel_min_height = 0.005,
      colour = "white", linewidth = 0.25, alpha = 0.95
    ) +
    geom_point(data = med_tbl, aes(x = median, y = term),
               inherit.aes = FALSE, size = 1.6, colour = "grey15") +
    geom_errorbarh(data = med_tbl,
                    aes(xmin = l95, xmax = u95, y = term),
                    inherit.aes = FALSE, height = 0,
                    linewidth = 0.45, colour = "grey15") +
    scale_fill_gradient2(low = fill_low, mid = "#F2E8D8",
                          high = fill_high, midpoint = 0,
                          guide = "none") +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = x_lab, y = NULL)
}
