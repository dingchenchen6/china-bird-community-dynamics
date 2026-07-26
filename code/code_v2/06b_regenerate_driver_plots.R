#!/usr/bin/env Rscript
## 06b_regenerate_driver_plots.R
##
## 把 stage-5 落到 derived_v2/ 的 brms driver fit RDS 重新出图：
##   - 单响应：雨林图（halfeye + jitter + 95% CRI）
##   - 多响应汇总：facet 雨林图
## 不再用之前简陋的 forest plot。

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(forcats); library(stringr)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_plots_advanced.R"))

P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_pilot_60sp_ar1")

message(sprintf("[stage-6b] regenerating driver raincloud plots for %s",
                RUN_LABEL))

# 找所有 brms_driver_fit_*_<RUN_LABEL>.rds
fit_files <- list.files(P$derived_v2,
  pattern = sprintf("^brms_driver_fit_.*_%s\\.rds$", RUN_LABEL),
  full.names = TRUE)
if (length(fit_files) == 0) {
  message("  no brms driver fits found for ", RUN_LABEL); quit(status = 0)
}
message(sprintf("  %d fit files found", length(fit_files)))

# 解析响应名（trend_xxx）
resp_of <- function(p) {
  bn <- basename(p)
  bn <- sub(paste0("_", RUN_LABEL, "\\.rds$"), "", bn)
  sub("^brms_driver_fit_", "", bn)
}

# 加载 + 抽 draws
fits <- lapply(fit_files, readRDS)
names(fits) <- vapply(fit_files, resp_of, character(1))

resp_pretty <- c(
  trend_corrected_richness = "Richness trend",
  trend_shannon            = "Shannon trend",
  trend_pd_prob            = "Faith's PD trend",
  trend_trait_volume       = "Trait-volume trend",
  trend_rao_q              = "Rao's Q trend",
  trend_mpd_prob           = "MPD trend"
)
title_for <- function(nm) {
  if (nm %in% names(resp_pretty)) resp_pretty[[nm]] else nm
}

## --- 1. 单响应雨林图（每个 brms fit 一张） ---------------------------------

for (nm in names(fits)) {
  draws <- brms_fixef_draws(fits[[nm]])
  draws$term <- gsub("^z_", "", draws$term)            # 去掉 z_ 前缀视觉简洁
  p <- raincloud_posterior(
    draws,
    title = sprintf("Drivers of %s — posterior distributions", title_for(nm)),
    subtitle = "brms (Gaussian) + cmdstanr | 4 chains | 95% CRI in solid bar",
    caption = "Halfeye = posterior density; bee swarm = thinned MCMC draws; bullet = posterior median.",
    x_lab = "Standardised coefficient"
  ) +
    theme_v2_pub(11) +
    theme(panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2))
  save_dual(p, sprintf("fig_driver_raincloud_%s_%s", nm, RUN_LABEL),
            width = 9.2, height = 5.6)
}
message(sprintf("  [ok] %d single-response raincloud figures", length(fits)))

## --- 2. 多响应汇总图 ------------------------------------------------------

draws_all <- purrr::imap_dfr(fits, function(fit, nm) {
  d <- brms_fixef_draws(fit)
  d$term <- gsub("^z_", "", d$term)
  d$response <- title_for(nm)
  d
})

med_within <- draws_all |>
  group_by(response, term) |>
  summarise(med = median(value), .groups = "drop")
draws_all <- draws_all |>
  left_join(med_within, by = c("response", "term")) |>
  mutate(term = forcats::fct_reorder(term, med))

multi_plot <- ggplot(draws_all, aes(value, term)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  ggdist::stat_halfeye(.width = c(0.5, 0.95), thickness = 0.55,
                       slab_alpha = 0.55, slab_colour = NA,
                       interval_colour = "grey25", point_size = 1.6,
                       fill = "#0E5A78") +
  facet_wrap(~ response, scales = "free", ncol = 2) +
  labs(
    title = "Drivers of community dynamic trends across responses",
    subtitle = sprintf("brms + cmdstanr | 4 chains | run=%s", RUN_LABEL),
    caption = "Halfeye = posterior density; thick bar = 50% CRI; thin bar = 95% CRI; point = median.",
    x = "Standardised coefficient", y = NULL
  ) +
  theme_v2_pub(11) +
  theme(panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2))
save_dual(multi_plot,
          sprintf("fig_driver_raincloud_multipanel_%s", RUN_LABEL),
          width = 12, height = 7.6)
message("  [ok] multi-response raincloud figure written")

## --- 3. 蜂群版（紧凑） ---------------------------------------------------

bee_plot <- ggplot(draws_all, aes(value, term, colour = response)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  ggbeeswarm::geom_quasirandom(alpha = 0.18, size = 0.4, width = 0.32,
                                groupOnX = FALSE, show.legend = FALSE) +
  geom_point(data = med_within |>
               left_join(draws_all |> distinct(response, term), by = c("response","term")) |>
               mutate(term = factor(term, levels = levels(draws_all$term))),
             aes(x = med, y = term, colour = response),
             inherit.aes = FALSE, size = 2.4) +
  scale_colour_manual(values = V2_PALETTES$qualitative[1:length(unique(draws_all$response))]) +
  labs(
    title = "Driver coefficients — bee-swarm overlay across responses",
    subtitle = sprintf("Each color = one response | run=%s", RUN_LABEL),
    caption = "Cloud = MCMC draws; large dot = posterior median.",
    x = "Standardised coefficient", y = NULL, colour = NULL
  ) +
  theme_v2_pub(11) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2))
save_dual(bee_plot,
          sprintf("fig_driver_beeswarm_overlay_%s", RUN_LABEL),
          width = 10.5, height = 6.8)
message("  [ok] bee-swarm overlay figure written")

message("[stage-6b] Done.")
