#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   物种水平占域趋势能否被生活史/形态/生态位性状解释（Q4）？
#   diet_specialization 与 habitat_breadth 是否是显著驱动？
#   系统发育随机效应是否必要？
#
# Objective / 分析目标:
#   v4 相对 v3：
#     (C1) loo::loo 替代 lo
#     (C9) 同时拟合含/不含系统发育随机效应两个模型，做 LOO 比较
#     正确报告 marginal R² vs conditional R²
#
# Input data / 输入数据:
#   results_v4/table_beta_species_<run_label>.csv（year_scaled 系数）
#   data/derived_v4/trait_extended_v4.rds（或 v3 fallback）
#   data/derived_v4/phylogeny_matched_v4.rds（或 v3 fallback）
#
# Main workflow / 主要流程:
#   1. 提取每物种 year_scaled 系数 → trend_i
#   2. 加载性状 + 系统发育 vcv
#   3. 拟合 brms 含 phylo random effect（M1）+ 不含（M0）
#   4. LOO 比较 + DHARMa 残差
#   5. 输出系数表 + R² 表 + LOO 表
#
# Key assumptions / 关键假设:
#   - brms + cmdstanr 可用
#   - 物种数 ≥ 30 才拟合 phylo 模型
#
# Main packages / 主要包:
#   brms, cmdstanr, loo, ape, DHARMa
#
# Output directory / 输出路径:
#   results_v4/table_trait_regression_coefs_<run_label>.csv
#   results_v4/table_trait_regression_R2_<run_label>.csv
#   results_v4/table_trait_regression_loo_comparison_<run_label>.csv
#   data/derived_v4/brms_trait_regression_<M0|M1>_<run_label>.qs
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(brms); library(cmdstanr); library(loo)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_seeds.R"))
source(file.path(CODE_V4, "utils_diagnostics.R"))
P <- ensure_v4_dirs()

set_seeds("14_trait_regression")
log_time("14", "Species-trait regression (Q4)")

is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

# ── 1. 加载种级 year_scaled 系数 ─────────────────────────────────────
beta_sp <- read_csv_safe(v4_file("results", paste0("table_beta_species_", run_label)))
if (is.null(beta_sp)) stop("[14] table_beta_species not found. Run 04b first.")

# year_scaled 系数取出（v4 中 param 可能叫 "year_scaled" 或类似）
trend_df <- beta_sp |>
  filter(grepl("year", param, ignore.case = TRUE) |
         grepl("year_scaled", covariate %||% param, ignore.case = TRUE)) |>
  group_by(species) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(species, trend_i = mean, trend_q025 = q025, trend_q975 = q975)

if (nrow(trend_df) == 0) {
  # fallback: 用 sp_trend 的 theil_sen 作为响应
  sp_trd <- read_csv_safe(v4_file("results", paste0("table_species_trend_", run_label)))
  if (!is.null(sp_trd)) {
    trend_df <- sp_trd |>
      filter(method == "theil_sen") |>
      select(species, trend_i = mean, trend_q025 = q025, trend_q975 = q975)
    message("[14] Fallback: using species_trend theil_sen as trend_i")
  } else {
    stop("[14] No species trend available")
  }
}

# ── 2. 加载性状 ──────────────────────────────────────────────────────
trait_ext <- safe_read(file.path(DIRS$derived, "trait_extended_v4.rds"), quiet = TRUE) %||%
             safe_read(file.path(DIRS$v3_derived, "trait_extended_v3.rds"), quiet = TRUE)
if (is.null(trait_ext)) stop("[14] trait_extended not found. Run 03b first.")

reg_df <- trend_df |>
  left_join(trait_ext, by = "species") |>
  drop_na(trend_i)

trait_cols_reg <- intersect(
  c("body_mass_g", "clutch_size", "avonet_hwi", "avonet_range_size",
    "diet_specialization", "habitat_breadth"),
  names(reg_df)
)

log_vars_14   <- intersect(TRAIT_VARS_LOG10, trait_cols_reg)
no_log_vars_14 <- intersect(TRAIT_VARS_NO_LOG, trait_cols_reg)

for (cc in log_vars_14) {
  reg_df[[paste0("z_", cc)]] <- as.numeric(scale(
    if_else(reg_df[[cc]] > 0, log10(reg_df[[cc]]), NA_real_)
  ))
}
for (cc in no_log_vars_14) {
  reg_df[[paste0("z_", cc)]] <- as.numeric(scale(reg_df[[cc]]))
}

z_trait_cols <- intersect(paste0("z_", trait_cols_reg), names(reg_df))
message(sprintf("[14] %d species × %d traits", nrow(reg_df), length(z_trait_cols)))

# ── 3. 加载系统发育 ─────────────────────────────────────────────────
phylo <- safe_read(file.path(DIRS$derived, "phylogeny_matched_v4.rds"), quiet = TRUE) %||%
         safe_read(file.path(DIRS$v3_derived, "phylogeny_matched_v3.rds"), quiet = TRUE) %||%
         safe_read(file.path(DIRS$v2_derived, "phylogeny_matched.rds"), quiet = TRUE)

A <- NULL
if (!is.null(phylo)) {
  matched <- intersect(reg_df$species, phylo$tip.label)
  if (length(matched) >= 30) {
    pruned <- ape::drop.tip(phylo, setdiff(phylo$tip.label, matched))
    A <- ape::vcv.phylo(pruned, corr = TRUE)
    reg_df <- reg_df |> filter(species %in% matched)
    message(sprintf("[14] Phylogeny matched: %d species", length(matched)))
  }
}

# ── 4. 拟合 M0（无 phylo）+ M1（含 phylo） ─────────────────────────
brms_prior_base <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "Intercept")
)

f_str_base <- paste0("trend_i ~ ", paste(z_trait_cols, collapse = " + "))

# M0: 无 phylo
f_m0 <- as.formula(f_str_base)
fit_m0 <- tryCatch(
  brm(formula = f_m0, data = reg_df, family = gaussian(),
      prior = brms_prior_base,
      iter = BRMS_ITER, warmup = BRMS_WARMUP,
      chains = BRMS_CHAINS, cores = max(1, parallel::detectCores() - 1),
      control = list(adapt_delta = BRMS_ADAPT_DELTA, max_treedepth = BRMS_MAX_TREED),
      seed = BRMS_SEED, backend = "cmdstanr", silent = 2, refresh = 0),
  error = function(e) { message("[14] M0 failed: ", e$message); NULL }
)

# M1: 含 phylo random effect
fit_m1 <- NULL
if (!is.null(A)) {
  brms_prior_phylo <- c(brms_prior_base,
                        prior(exponential(1), class = "sd"))
  f_m1 <- bf(as.formula(paste0(f_str_base, " + (1 | gr(species, cov = A))")))
  fit_m1 <- tryCatch(
    brm(formula = f_m1, data = reg_df, data2 = list(A = A),
        family = gaussian(), prior = brms_prior_phylo,
        iter = BRMS_ITER, warmup = BRMS_WARMUP,
        chains = BRMS_CHAINS, cores = max(1, parallel::detectCores() - 1),
        control = list(adapt_delta = BRMS_ADAPT_DELTA, max_treedepth = BRMS_MAX_TREED),
        seed = BRMS_SEED + 1, backend = "cmdstanr", silent = 2, refresh = 0),
    error = function(e) { message("[14] M1 failed: ", e$message); NULL }
  )
}

# ── 5. 系数表（用 M1，如有；否则 M0） ───────────────────────────────
fit_use <- fit_m1 %||% fit_m0
if (!is.null(fit_use)) {
  coef_df <- as_tibble(fixef(fit_use), rownames = "term") |>
    rename(estimate = Estimate, q025 = `Q2.5`, q975 = `Q97.5`) |>
    mutate(term = gsub("^b_", "", term),
            model = ifelse(identical(fit_use, fit_m1), "M1_phylo", "M0_basic"))
  write_csv(coef_df, v4_file("results", paste0("table_trait_regression_coefs_", run_label)))
  checkpoint_save(fit_use,
                   paste0("brms_trait_regression_",
                          ifelse(identical(fit_use, fit_m1), "M1", "M0"),
                          "_", run_label),
                   subdir = "derived")

  # R²
  r2_rows <- list()
  for (mod_name in c("M0_basic", "M1_phylo")) {
    fit <- if (mod_name == "M0_basic") fit_m0 else fit_m1
    if (is.null(fit)) next
    r2c <- tryCatch(bayes_R2(fit), error = function(e) NULL)
    r2m <- tryCatch(bayes_R2(fit, re_formula = NA), error = function(e) NULL)
    r2_rows[[mod_name]] <- tibble(
      model = mod_name,
      R2_conditional = r2c[1, "Estimate"] %||% NA_real_,
      R2_conditional_q025 = r2c[1, "Q2.5"] %||% NA_real_,
      R2_conditional_q975 = r2c[1, "Q97.5"] %||% NA_real_,
      R2_marginal = r2m[1, "Estimate"] %||% NA_real_,
      R2_marginal_q025 = r2m[1, "Q2.5"] %||% NA_real_,
      R2_marginal_q975 = r2m[1, "Q97.5"] %||% NA_real_
    )
  }
  if (length(r2_rows) > 0) {
    write_csv(bind_rows(r2_rows),
              v4_file("results", paste0("table_trait_regression_R2_", run_label)))
  }
}

# ── 6. LOO 比较（C1/C9 修复：loo::loo + loo_compare） ───────────────
loo_rows <- list()
loo_m0 <- if (!is.null(fit_m0)) tryCatch(loo::loo(fit_m0), error = function(e) NULL) else NULL
loo_m1 <- if (!is.null(fit_m1)) tryCatch(loo::loo(fit_m1), error = function(e) NULL) else NULL

if (!is.null(loo_m0)) {
  loo_rows[["M0"]] <- tibble(
    model    = "M0_basic",
    elpd_loo = loo_m0$estimates["elpd_loo", "Estimate"],
    elpd_se  = loo_m0$estimates["elpd_loo", "SE"],
    p_loo    = loo_m0$estimates["p_loo",    "Estimate"],
    looic    = loo_m0$estimates["looic",    "Estimate"]
  )
  checkpoint_save(loo_m0, paste0("loo_M0_trait_", run_label), subdir = "derived")
}
if (!is.null(loo_m1)) {
  loo_rows[["M1"]] <- tibble(
    model    = "M1_phylo",
    elpd_loo = loo_m1$estimates["elpd_loo", "Estimate"],
    elpd_se  = loo_m1$estimates["elpd_loo", "SE"],
    p_loo    = loo_m1$estimates["p_loo",    "Estimate"],
    looic    = loo_m1$estimates["looic",    "Estimate"]
  )
  checkpoint_save(loo_m1, paste0("loo_M1_trait_", run_label), subdir = "derived")
}

if (length(loo_rows) > 0) {
  loo_df <- bind_rows(loo_rows)

  # 比较
  if (!is.null(loo_m0) && !is.null(loo_m1)) {
    cmp <- loo::loo_compare(loo_m0, loo_m1)
    cmp_df <- as_tibble(cmp, rownames = "rank") |>
      mutate(model = rownames(cmp))
    write_csv(cmp_df,
              v4_file("results", paste0("table_trait_regression_loo_comparison_", run_label)))
    message(sprintf("[14] LOO compare: best model = %s, Δelpd = %.2f ± %.2f",
                    rownames(cmp)[1], cmp[2, "elpd_diff"], cmp[2, "se_diff"]))
  }
  write_csv(loo_df, v4_file("results", paste0("table_trait_regression_loo_", run_label)))
}

# ── 7. DHARMa 残差 ──────────────────────────────────────────────────
if (!is.null(fit_use)) {
  check_dharma_gate(fit_use, save_dir = DIRS$figures,
                     stem = paste0("dharma_trait_regression_", run_label))
}

log_time("14", "DONE")
