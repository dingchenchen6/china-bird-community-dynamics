#!/usr/bin/env Rscript
## 04b_recover_diagnostics.R
##
## 用途：stage-4 拟合 RDS 已存在但 PPC 阶段 OOM 中断时，
## 用本脚本接力跑：MCMC 诊断、trace 图、社区/物种系数表与图、psi.samples thin。
## 跳过 ppcOcc 以避免 OOM；PPC 留作 stage-5 后处理时用更轻量方式补做。
##
## 调用：V2_RUN_LABEL=v2_full_200sp_ar1 Rscript code_v2/04b_recover_diagnostics.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(coda); library(ggplot2); library(patchwork); library(forcats)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_diagnostics.R"))

P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")
PSI_MAX_DRAWS <- as.integer(Sys.getenv("V2_PSI_MAX_DRAWS", "200"))

fit_path <- v2_file("derived", paste0("tMsPGOcc_fit_", RUN_LABEL), "rds")
stopifnot(file.exists(fit_path))
message(sprintf("[recover] loading fit: %s (%.1f GB)",
                fit_path, file.info(fit_path)$size / 1024^3))
fit <- readRDS(fit_path)
message(sprintf("  n.chains=%d", fit$n.chains))

## --- 1. 诊断 ---------------------------------------------------------------

diag_tbl <- bind_rows(
  extract_mcmc_diag(fit, "beta.comm.samples"),
  extract_mcmc_diag(fit, "alpha.comm.samples"),
  extract_mcmc_diag(fit, "tau.sq.beta.samples"),
  extract_mcmc_diag(fit, "tau.sq.alpha.samples"),
  extract_mcmc_diag(fit, "sigma.sq.t.samples"),
  extract_mcmc_diag(fit, "rho.samples")
)
write_csv(diag_tbl, v2_file("results", paste0("mcmc_diagnostics_", RUN_LABEL)))

p_diag <- plot_rhat_ess(diag_tbl, RUN_LABEL) +
  patchwork::plot_annotation(
    title = sprintf("MCMC convergence diagnostics - %s", RUN_LABEL),
    theme = theme_v2_pub(11))
save_dual(p_diag, paste0("fig_mcmc_diag_", RUN_LABEL),
          width = 11, height = 4.6)
message(sprintf("  R-hat median=%.3f  90%%=%.3f  max=%.3f",
        median(diag_tbl$rhat, na.rm=TRUE),
        quantile(diag_tbl$rhat, 0.9, na.rm=TRUE),
        max(diag_tbl$rhat, na.rm=TRUE)))

## --- 2. 社区系数 ----------------------------------------------------------

beta_comm <- as.matrix(fit$beta.comm.samples)
alpha_comm <- as.matrix(fit$alpha.comm.samples)

summarise_post <- function(M, kind) tibble(
  parameter = colnames(M),
  mean   = apply(M, 2, mean),
  median = apply(M, 2, median),
  l95    = apply(M, 2, quantile, 0.025),
  u95    = apply(M, 2, quantile, 0.975),
  p_pos  = apply(M, 2, function(x) mean(x > 0)),
  kind   = kind
)
post_summary <- bind_rows(
  summarise_post(beta_comm, "occupancy"),
  summarise_post(alpha_comm, "detection")
) |> mutate(sig95 = (l95 > 0) | (u95 < 0))
write_csv(post_summary,
          v2_file("results", paste0("table_community_coefficients_", RUN_LABEL)))

cp <- post_summary |>
  mutate(parameter = forcats::fct_reorder(parameter, mean)) |>
  ggplot(aes(mean, parameter, xmin = l95, xmax = u95, colour = sig95)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#8B2E1E"),
                       guide = "none") +
  facet_wrap(~ kind, scales = "free", ncol = 2) +
  labs(title = sprintf("Community-level coefficients - %s", RUN_LABEL),
       subtitle = "Posterior mean +/- 95% CRI; red = CRI excludes 0",
       x = "Coefficient (logit scale)", y = NULL) +
  theme_v2_pub(11)
save_dual(cp, paste0("fig_community_caterpillar_", RUN_LABEL),
          width = 11, height = 5.6)

## --- 3. 物种系数 heatmap --------------------------------------------------

beta_sp <- as.matrix(fit$beta.samples)
sp_long <- as.data.frame(beta_sp) |>
  mutate(iter = row_number()) |>
  pivot_longer(-iter, names_to = "term", values_to = "value") |>
  separate(term, into = c("covariate", "species"), sep = "-",
           extra = "merge") |>
  group_by(covariate, species) |>
  summarise(mean = mean(value),
            l95 = quantile(value, 0.025),
            u95 = quantile(value, 0.975),
            .groups = "drop") |>
  mutate(sig95 = (l95 > 0) | (u95 < 0))
write_csv(sp_long,
          v2_file("results", paste0("table_species_coefficients_", RUN_LABEL)))

n_species_plot <- length(unique(sp_long$species))
heat <- sp_long |>
  ggplot(aes(covariate, species, fill = mean)) +
  geom_tile(colour = "white", linewidth = 0.05) +
  geom_text(aes(label = ifelse(sig95, "*", "")), size = 2.2,
            colour = "grey15") +
  scale_fill_v2_diverging(name = "Posterior mean beta", limits = c(-2, 2)) +
  labs(title = sprintf("Species x covariate occupancy effects - %s", RUN_LABEL),
       subtitle = "* = 95% CRI excludes 0",
       x = NULL, y = NULL) +
  theme_v2_pub(8.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = max(3.5, 8 - n_species_plot / 30)))
save_dual(heat, paste0("fig_species_coef_heatmap_", RUN_LABEL),
          width = 11, height = max(8, n_species_plot * 0.07))

## --- 4. psi.samples thin --------------------------------------------------

psi_samps <- fit$psi.samples
n_draws <- dim(psi_samps)[1]
keep <- if (n_draws > PSI_MAX_DRAWS)
  round(seq(1, n_draws, length.out = PSI_MAX_DRAWS)) else seq_len(n_draws)
psi_thin <- psi_samps[keep, , , , drop = FALSE]
saveRDS(list(psi_samples_thinned = psi_thin,
              keep_draws = keep, dims = dim(psi_samps)),
        v2_file("derived", paste0("psi_samples_thinned_", RUN_LABEL), "rds"))
message(sprintf("  psi.samples thinned: %d / %d draws kept (%.1f GB)",
                length(keep), n_draws,
                file.info(v2_file("derived",
                  paste0("psi_samples_thinned_", RUN_LABEL), "rds"))$size / 1024^3))

## --- 5. 写说明 ------------------------------------------------------------

writeLines(c(
  sprintf("# v2 stage-4 main fit (%s)", RUN_LABEL),
  sprintf("- 物种数: %d", dim(psi_samps)[2]),
  sprintf("- 站点数: %d", dim(psi_samps)[3]),
  sprintf("- 主期数: %d", dim(psi_samps)[4]),
  sprintf("- 总后验 draws: %d", n_draws),
  sprintf("- thinned 留下: %d", length(keep)),
  sprintf("- R-hat median: %.3f / 90%%: %.3f / max: %.3f",
          median(diag_tbl$rhat, na.rm=TRUE),
          quantile(diag_tbl$rhat, 0.9, na.rm=TRUE),
          max(diag_tbl$rhat, na.rm=TRUE)),
  "- PPC: 跳过（OOM 风险）；stage-5 用 site-stratified subset 重做。"
), v2_file("results",
           paste0("multispecies_dynamic_main_summary_", RUN_LABEL), "md"))

message("[recover] Done.")
