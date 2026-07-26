#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   补齐 v3 缺失的"标准多物种占域 MCMC 诊断图"（A 组），
#   让审稿人能直接看到链是否充分混合、ESS 是否够大、AR1 ρ 是否合理
#
# Objective / 分析目标:
#   - fig_v4_mcmc_rhat_histogram.png：所有参数 R-hat 直方 + 1.05 阈值线
#   - fig_v4_mcmc_ess_histogram.png：所有参数 ESS 直方 + 200 阈值线
#   - fig_v4_mcmc_trace_core_params.png：beta.comm / alpha.comm / rho / sigma.sq
#     每参数 4 链 trace 叠合
#   - fig_v4_mcmc_density_core_params.png：4 链 posterior density 叠合
#   - fig_v4_mcmc_ar1_rho_density_by_species.png：物种 ρ 后验密度
#
# Input data / 输入数据:
#   data/derived_v4/stMsPGOcc_fit_<run_label>_combined.qs
#   results_v4/table_convergence_diagnostics_<run_label>.csv
#
# Main workflow / 主要流程:
#   1. 读 fit + diag 表
#   2. R-hat 直方 / ESS 直方（ggplot）
#   3. 核心参数 trace + density（按 chain 着色）
#   4. AR1 ρ 物种密度
#
# Key assumptions / 关键假设:
#   多链合并后 phi/sigma.sq/rho/beta.comm samples 含 chain 维（3D）
#   单链时 fit$rhat 不可用，则只画 ESS + density（无 chain 叠合）
#
# Main packages / 主要包:
#   ggplot2, patchwork, dplyr, tidyr, qs
#
# Output directory / 输出路径:
#   figures_v4/fig_v4_mcmc_*_<run_label>.{png,pdf}
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(qs)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
P <- ensure_v4_dirs()

log_time("05b", "MCMC diagnostic plots")

is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

# ── 0. 通用主题 ─────────────────────────────────────────────────────
theme_v4_pub <- function(base_size = NATURE_PT) {
  theme_minimal(base_size = base_size, base_family = NATURE_FONT) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
      axis.line  = element_line(linewidth = NATURE_LINE_AXIS, colour = "grey20"),
      axis.ticks = element_line(linewidth = NATURE_LINE_TICK, colour = "grey20"),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
}

save_v4 <- function(p, stem, width_mm = NATURE_WIDTH_L, height_mm = 90) {
  ggsave(file.path(DIRS$figures, paste0(stem, ".png")),
         plot = p, width = width_mm, height = height_mm,
         units = "mm", dpi = NATURE_DPI)
  ggsave(file.path(DIRS$figures, paste0(stem, ".pdf")),
         plot = p, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
}

# ── 1. R-hat / ESS 直方 ─────────────────────────────────────────────
diag_path <- v4_file("results", paste0("table_convergence_diagnostics_", run_label))
diag_df <- read_csv_safe(diag_path)
if (!is.null(diag_df) && nrow(diag_df) > 0 && any(!is.na(diag_df$rhat))) {
  p_rhat <- ggplot(diag_df |> filter(!is.na(rhat)), aes(rhat)) +
    geom_histogram(bins = 40, fill = NATURE_ACCENT, alpha = 0.75) +
    geom_vline(xintercept = RHAT_THRESHOLD, linetype = 2, colour = "red", linewidth = 0.3) +
    labs(x = expression(hat(R)), y = "Count of parameters",
         title = sprintf("R-hat (threshold = %.2f)", RHAT_THRESHOLD)) +
    facet_wrap(~ group, scales = "free_y") +
    theme_v4_pub()
  save_v4(p_rhat, paste0("fig_v4_mcmc_rhat_histogram_", run_label))
} else {
  message("[05b] No valid R-hat values (likely single chain) — skip rhat histogram")
}

if (!is.null(diag_df) && any(!is.na(diag_df$ess))) {
  p_ess <- ggplot(diag_df |> filter(!is.na(ess)), aes(ess)) +
    geom_histogram(bins = 40, fill = NATURE_ACCENT, alpha = 0.75) +
    geom_vline(xintercept = ESS_THRESHOLD, linetype = 2, colour = "red", linewidth = 0.3) +
    labs(x = "Effective sample size",
         y = "Count of parameters",
         title = sprintf("ESS (threshold = %d)", ESS_THRESHOLD)) +
    facet_wrap(~ group, scales = "free") +
    theme_v4_pub()
  save_v4(p_ess, paste0("fig_v4_mcmc_ess_histogram_", run_label))
}

# ── 2. Trace + Density for core community params ────────────────────
fit_path <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, "_combined.qs"))
if (!file.exists(fit_path)) {
  fit_path <- file.path(DIRS$derived, paste0("stMsPGOcc_fit_", run_label, ".qs"))
}
fit <- safe_read(fit_path, quiet = TRUE)

if (!is.null(fit)) {
  # 把 samples slot 转为长表 [iter, chain, param, value]
  samples_to_long <- function(samples, param_prefix) {
    d <- dim(samples)
    if (is.null(d) || length(d) == 0) return(NULL)
    if (length(d) == 2) {
      # [param, sample] 单链
      sm <- as_tibble(t(samples))
      sm$iter <- seq_len(nrow(sm)); sm$chain <- 1L
      out <- sm |> pivot_longer(-c(iter, chain), names_to = "param", values_to = "value")
      out$param <- paste0(param_prefix, "_", out$param)
      return(out)
    }
    if (length(d) == 3) {
      # [param, sample, chain]
      out_list <- list()
      for (ch in seq_len(d[3])) {
        m <- samples[, , ch]
        sm <- as_tibble(t(m))
        sm$iter <- seq_len(nrow(sm)); sm$chain <- ch
        out_list[[ch]] <- sm |> pivot_longer(-c(iter, chain), names_to = "param", values_to = "value")
      }
      out <- bind_rows(out_list)
      out$param <- paste0(param_prefix, "_", out$param)
      return(out)
    }
    NULL
  }

  long_list <- list()
  for (slot in c("beta.comm.samples", "alpha.comm.samples",
                  "rho.samples", "sigma.sq.samples", "phi.samples")) {
    if (slot %in% names(fit) && !is.null(fit[[slot]])) {
      stem <- gsub("\\.samples$", "", slot)
      ll <- samples_to_long(fit[[slot]], stem)
      if (!is.null(ll)) long_list[[slot]] <- ll
    }
  }

  if (length(long_list) > 0) {
    long_df <- bind_rows(long_list)
    top_params <- long_df |>
      group_by(param) |>
      summarise(n = n(), .groups = "drop") |>
      slice_head(n = 12)
    long_df <- long_df |> filter(param %in% top_params$param)

    p_trace <- ggplot(long_df, aes(iter, value, colour = factor(chain))) +
      geom_line(alpha = 0.6, linewidth = 0.2) +
      facet_wrap(~ param, scales = "free_y", ncol = 3) +
      scale_colour_brewer(palette = "Set1", name = "Chain") +
      labs(x = "Iteration", y = "Value", title = "Trace (core community parameters)") +
      theme_v4_pub() +
      theme(strip.text = element_text(size = 6))
    save_v4(p_trace, paste0("fig_v4_mcmc_trace_core_params_", run_label),
            height_mm = 130)

    p_dens <- ggplot(long_df, aes(value, colour = factor(chain))) +
      geom_density(alpha = 0.6, linewidth = 0.3) +
      facet_wrap(~ param, scales = "free", ncol = 3) +
      scale_colour_brewer(palette = "Set1", name = "Chain") +
      labs(x = "Posterior value", y = "Density",
           title = "Posterior density by chain (core community parameters)") +
      theme_v4_pub() +
      theme(strip.text = element_text(size = 6))
    save_v4(p_dens, paste0("fig_v4_mcmc_density_core_params_", run_label),
            height_mm = 130)
  }

  # ── 3. AR1 ρ by species（如有 fit$rho.samples 且每物种独立） ────────
  # spOccupancy tMsPGOcc/stMsPGOcc 的 rho.samples 通常是 [sample, species]
  # 或 [sample, species, chain]
  if ("rho.samples" %in% names(fit)) {
    rho <- fit$rho.samples
    d <- dim(rho)
    if (!is.null(d) && length(d) >= 2) {
      if (length(d) == 2 && d[2] > 1) {
        # [sample, species]
        sp_names <- rownames(fit$y) %||% paste0("sp", seq_len(d[2]))
        rho_df <- as_tibble(rho)
        names(rho_df) <- sp_names
        rho_long <- rho_df |>
          pivot_longer(everything(), names_to = "species", values_to = "rho")

        # 高亮强 AR1 物种（mean |rho| > 0.3）
        strong_sp <- rho_long |>
          group_by(species) |>
          summarise(mean_abs = mean(abs(rho), na.rm = TRUE), .groups = "drop") |>
          filter(mean_abs > 0.3) |>
          pull(species)

        p_rho <- ggplot(rho_long, aes(rho, group = species)) +
          geom_density(colour = "grey60", alpha = 0.1, linewidth = 0.15) +
          geom_density(data = rho_long |> filter(species %in% strong_sp),
                       aes(colour = species), linewidth = 0.4) +
          geom_vline(xintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.2) +
          labs(x = expression(rho ~ "(AR1 temporal correlation)"),
               y = "Density",
               title = sprintf("AR1 ρ posterior by species (%d strong-AR1 highlighted)",
                                length(strong_sp))) +
          theme_v4_pub() +
          theme(legend.position = "none")
        save_v4(p_rho, paste0("fig_v4_mcmc_ar1_rho_density_by_species_", run_label))
      }
    }
  }
}

log_time("05b", "DONE")
