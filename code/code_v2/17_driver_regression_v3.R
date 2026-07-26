#!/usr/bin/env Rscript
## 17_driver_regression_v3.R
##
## 驱动回归 v3：
##   1) VIF 协变量筛选（car::vif 迭代删 >5）
##   2) brms + cmdstanr + Student + 6 个生态意义交互项 + gp(centroid_lon, centroid_lat, k=20)
##   3) ANOVA Type II + III（car::Anova）
##   4) 方差分解三轨：LMG (relaimpo) + permutation importance (ranger) + varpart (vegan)
##   5) DHARMa 残差诊断
##
## 响应变量：v2 已有 4 个 trend + v3 新指标 2-3 个（如 nri_trend, fric_trend）
##
## Output:
##   data/derived_v3/brms_driver_v3_<resp>.rds (per response)
##   results_v3/vif_v3.csv
##   results_v3/anova_v3_<resp>.csv (II + III)
##   results_v3/importance_v3_<resp>.csv (LMG + RF + varpart)
##   figures_v3/fig_v3_vif_heatmap.{png,pdf}
##   figures_v3/fig_v3_driver_importance_triple.{png,pdf}
##   figures_v3/fig_v3_interaction_marginal.{png,pdf}
##   figures_v3/dharma_v3_<resp>.png

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(forcats); library(ggplot2); library(patchwork)
  library(car); library(relaimpo); library(ranger)   # vegan loaded ad-hoc to avoid masking dplyr::select
  library(brms); library(DHARMa); library(ggcorrplot)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
source(file.path(CODE_V2, "utils_diagnostics.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")
N_BRMS_ITER <- as.integer(Sys.getenv("V2_BRMS_ITER", "2000"))

message(sprintf("[stage-17] v3 driver regression | %s", RUN_LABEL))

## --- 1. 读取趋势 + 环境 ------------------------------------------------

trends_v2 <- read_csv(v2_file("results",
                paste0("table_grid_trends_with_cri_", RUN_LABEL)),
                show_col_types = FALSE)
grid_env <- readRDS(file.path(P$derived_v2,
                "grid_environment_dynamic_occupancy.rds"))

# v3 新指标趋势：若 stage 16 已完成则加入
v3_metric_path <- v3_file("results",
                  paste0("table_community_metrics_v3_", RUN_LABEL))
if (file.exists(v3_metric_path)) {
  v3_metrics <- read_csv(v3_metric_path, show_col_types = FALSE)
  # 按 grid 算 5 期 trend
  period_idx <- seq_len(5) - 3
  v3_trends <- v3_metrics |>
    filter(metric %in% c("nri", "nti", "fric", "fdiv", "fdis")) |>
    group_by(grid_cell, metric) |>
    summarise(
      trend_mean = if (sum(is.finite(value_mean)) >= 3)
        sum(period_idx[is.finite(value_mean)] *
              (value_mean[is.finite(value_mean)] -
                 mean(value_mean[is.finite(value_mean)]))) /
          sum(period_idx[is.finite(value_mean)]^2) else NA_real_,
      .groups = "drop") |>
    mutate(trend_l95 = NA_real_, trend_u95 = NA_real_,
            p_pos = NA_real_, sig95 = NA)
  trends_all <- bind_rows(trends_v2, v3_trends)
} else {
  trends_all <- trends_v2
  message("  v3 metrics not ready; using v2 trends only")
}

## --- 2. 协变量集（含森林面积 + 6 生态交互） ------------------------------

main_vars <- c("bio1", "bio4", "bio12", "elev_sd", "npp_mean",
                "hfi_mean", "landcover_trees",          # NEW: forest
                "landcover_built", "habitat_diversity_shannon",
                "centroid_lat")
main_vars <- intersect(main_vars, names(grid_env))

interaction_pairs <- list(
  c("bio4", "elev_sd"),
  c("npp_mean", "landcover_trees"),
  c("hfi_mean", "landcover_built"),
  c("bio1", "centroid_lat"),
  c("landcover_trees", "habitat_diversity_shannon"),
  c("bio12", "npp_mean")
)
# 仅当两项都在 main_vars 中保留
interaction_pairs <- Filter(function(p) all(p %in% main_vars), interaction_pairs)

## --- 3. VIF 迭代 ------------------------------------------------------

# 用 trend_corrected_richness 做 VIF（其他响应共用同一组 X）
resp_focal <- "corrected_richness"
trend_focal <- trends_all |> filter(metric == resp_focal) |>
  select(grid_cell, trend_mean) |>
  rename(y = trend_mean)
df_vif <- trend_focal |>
  inner_join(grid_env |> select(grid_cell, all_of(main_vars)),
              by = "grid_cell") |>
  drop_na()
for (v in main_vars) df_vif[[paste0("z_", v)]] <-
  as.numeric(scale(df_vif[[v]]))

keep <- main_vars
vif_log <- tibble(step = character(), variable = character(), vif = numeric())
repeat {
  if (length(keep) <= 2) break
  rhs <- paste0("z_", keep, collapse = " + ")
  f <- as.formula(sprintf("y ~ %s", rhs))
  m <- lm(f, data = df_vif)
  v <- car::vif(m)
  vif_log <- bind_rows(vif_log,
    tibble(step = "iterative", variable = names(v), vif = unname(v)))
  if (all(v <= 5)) break
  drop_v <- names(v)[which.max(v)]
  drop_v <- sub("^z_", "", drop_v)
  keep <- setdiff(keep, drop_v)
  message("  drop ", drop_v, " (VIF ", round(max(v), 2), ")")
}
write_csv(vif_log, v3_file("results", "vif_v3"))
message("  VIF retained: ", paste(keep, collapse = ", "))

# 仅保留两项都在 keep 中的交互
interaction_pairs <- Filter(function(p) all(p %in% keep), interaction_pairs)
message("  Interactions: ",
        paste(sapply(interaction_pairs, paste, collapse=":"), collapse="; "))

## --- 4. VIF heatmap ---------------------------------------------------

cor_mat <- cor(grid_env[, main_vars], use = "pairwise.complete.obs",
                method = "spearman")
p_corr <- ggcorrplot::ggcorrplot(
  cor_mat, type = "lower", outline.col = "white",
  lab = TRUE, lab_size = 2.7,
  colors = c("#0E5A78", "white", "#8B2E1E"),
  ggtheme = ggplot2::theme_minimal(base_family = "Helvetica")) +
  labs(title = "Driver covariate correlation matrix",
       subtitle = "Spearman | lower triangle | full candidate set",
       caption = sprintf("VIF screening retained: %s", paste(keep, collapse=", "))) +
  theme_v2_pub(10.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_dual_v3(p_corr, "fig_v3_vif_heatmap", width = 9, height = 7.5)

## --- 5. 响应 × brms 双轨拟合（gp + interactions）--------------------------

resp_set <- intersect(
  c("corrected_richness", "shannon", "pd_prob", "trait_volume",
    "nri", "nti", "fric", "fdiv", "fdis"),
  unique(trends_all$metric))

build_df <- function(resp_metric) {
  trends_all |> filter(metric == resp_metric) |>
    select(grid_cell, trend_mean) |>
    rename(y = trend_mean) |>
    inner_join(grid_env |> select(grid_cell, all_of(keep),
                                    centroid_lon, centroid_lat),
                by = "grid_cell") |>
    drop_na() |>
    mutate(across(all_of(keep), ~ as.numeric(scale(.x)),
                  .names = "z_{.col}"))
}

interactions_chr <- sapply(interaction_pairs,
                            function(p) paste0("z_", p[1], ":z_", p[2]))
formula_str <- function() {
  paste0("y ~ ",
         paste(c(paste0("z_", keep), interactions_chr,
                 "gp(centroid_lon, centroid_lat, k = 20)"),
                collapse = " + "))
}

driver_results <- list()
for (resp_metric in resp_set) {
  message(sprintf("[stage-17] brms for %s", resp_metric))
  df <- build_df(resp_metric)
  if (nrow(df) < 100) { message("  skip n<100"); next }

  fit_path <- v3_file("derived",
                       paste0("brms_driver_v3_", resp_metric), "rds")
  if (file.exists(fit_path) &&
      identical(Sys.getenv("V2_REUSE_FIT", "false"), "true")) {
    fit <- readRDS(fit_path)
  } else {
    fit <- tryCatch(
      brms::brm(as.formula(formula_str()), data = df,
                family = student(),
                chains = 4, iter = N_BRMS_ITER, warmup = N_BRMS_ITER %/% 2,
                cores = 4, backend = "cmdstanr",
                refresh = 0, save_pars = brms::save_pars(all = TRUE),
                control = list(adapt_delta = 0.95, max_treedepth = 12)),
      error = function(e) { message("    brms err: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    saveRDS(fit, fit_path)
  }

  # ANOVA (II + III) on lm equivalent（去掉 gp 项）
  lm_form <- as.formula(paste0("y ~ ",
              paste(c(paste0("z_", keep), interactions_chr),
                     collapse = " + ")))
  lm_fit <- lm(lm_form, data = df)
  ano2 <- tryCatch(broom::tidy(car::Anova(lm_fit, type = 2)),
                    error = function(e) NULL)
  ano3 <- tryCatch(broom::tidy(car::Anova(lm_fit, type = 3)),
                    error = function(e) NULL)
  anova_out <- bind_rows(
    if (!is.null(ano2)) ano2 |> mutate(type = "II"),
    if (!is.null(ano3)) ano3 |> mutate(type = "III"))
  write_csv(anova_out,
            v3_file("results", paste0("anova_v3_", resp_metric)))

  # LMG (relaimpo): 不带交互
  lm_main <- lm(as.formula(paste0("y ~ ",
                paste(paste0("z_", keep), collapse = " + "))),
                data = df)
  lmg <- tryCatch(relaimpo::calc.relimp(lm_main, type = "lmg",
                                          rela = TRUE)@lmg,
                    error = function(e) {message("    lmg err"); NULL})

  # ranger permutation importance
  rf_df <- df |> select(y, all_of(paste0("z_", keep)))
  rf <- ranger::ranger(y ~ ., data = rf_df, num.trees = 2000,
                        importance = "permutation",
                        respect.unordered.factors = TRUE)
  rf_imp <- rf$variable.importance
  rf_imp <- rf_imp / sum(abs(rf_imp))

  # varpart 4 组
  groups <- list(
    climate = intersect(c("bio1","bio4","bio12"), keep),
    topo    = intersect(c("elev_sd"), keep),
    human   = intersect(c("hfi_mean","landcover_trees","landcover_built","habitat_diversity_shannon"), keep),
    space   = intersect(c("centroid_lat"), keep))
  groups <- Filter(function(x) length(x) > 0, groups)
  make_X <- function(vars) as.matrix(df[, paste0("z_", vars), drop = FALSE])
  vp_call <- tryCatch({
    if (length(groups) >= 2)
      do.call(vegan::varpart, c(list(df$y), lapply(groups, make_X)))
    else NULL
  }, error = function(e) NULL)
  vp_tab <- if (!is.null(vp_call))
    as_tibble(vp_call$part$indfract, rownames = "fraction") else NULL

  imp_tbl <- tibble(
    variable = paste0("z_", keep),
    lmg      = if (!is.null(lmg)) unname(lmg[paste0("z_", keep)]) else NA_real_,
    rf       = unname(rf_imp[paste0("z_", keep)])
  ) |> mutate(variable = sub("^z_", "", variable))
  if (!is.null(vp_tab)) {
    # 把 group-level 贡献附加到表（用 [a]/[b]/[c]/[d] 单组分）
    pure <- vp_tab |> filter(grepl("^\\[\\w\\]$", fraction))
    if (nrow(pure) > 0) {
      ord <- c("[a]", "[b]", "[c]", "[d]")
      imp_tbl$varpart_group <- NA_real_
      for (gi in seq_along(groups)) {
        gname <- names(groups)[gi]
        adj_col <- grep("^[Aa]dj", names(pure), value = TRUE)[1]
        val <- pure[[adj_col]][gi]
        for (v in groups[[gi]]) imp_tbl$varpart_group[imp_tbl$variable == v] <- val
      }
    }
  }
  write_csv(imp_tbl,
            v3_file("results", paste0("importance_v3_", resp_metric)))

  # DHARMa 残差
  dh <- tryCatch(dharma_from_brms(fit, n_sim = 300,
                                    run_label = paste0("v3_", resp_metric)),
                 error = function(e) {message("    DHARMa err: ", conditionMessage(e)); NULL})

  driver_results[[resp_metric]] <- list(fit = fit, imp = imp_tbl, anova = anova_out)
}

## --- 6. 全部响应的 importance 三轨对比图 -------------------------------

all_imp <- bind_rows(lapply(names(driver_results), function(nm) {
  driver_results[[nm]]$imp |> mutate(response = nm)
}))
imp_long <- all_imp |>
  pivot_longer(c(lmg, rf, varpart_group),
                names_to = "method", values_to = "importance") |>
  mutate(method = recode(method,
    lmg = "LMG (relaimpo)",
    rf = "RF permutation (ranger)",
    varpart_group = "varpart group (vegan)"))

p_imp <- ggplot(imp_long,
                 aes(forcats::fct_reorder(variable, importance,
                                            .fun = function(x) sum(x, na.rm = TRUE)),
                     importance, fill = method)) +
  geom_col(position = position_dodge(0.78), width = 0.7) +
  facet_wrap(~ response, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("LMG (relaimpo)" = "#0E5A78",
                                  "RF permutation (ranger)" = "#C9784A",
                                  "varpart group (vegan)" = "#3C8C5A"),
                     name = NULL) +
  coord_flip() +
  labs(title = "Driver importance across responses: LMG vs RF vs varpart",
       subtitle = sprintf("Standardised relative importance | run=%s", RUN_LABEL),
       x = NULL, y = "Relative importance (proportion of explained variation)") +
  theme_v2_pub(10.5) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))
save_dual_v3(p_imp, "fig_v3_driver_importance_triple",
              width = 12, height = max(6, 2 * length(driver_results)))

## --- 7. 交互项 conditional_effects 图 ---------------------------------

if (length(driver_results) > 0 && length(interaction_pairs) > 0) {
  fit1 <- driver_results[[1]]$fit
  cond_plots <- list()
  for (pair in interaction_pairs) {
    e <- paste0("z_", pair[1], ":z_", pair[2])
    ce <- tryCatch(brms::conditional_effects(fit1, effects = e,
                                                prob = 0.9),
                    error = function(err) NULL)
    if (!is.null(ce)) {
      p <- plot(ce, plot = FALSE)[[1]] +
        labs(title = sprintf("%s x %s", pair[1], pair[2])) +
        theme_v2_pub(9.5)
      cond_plots[[paste(pair, collapse=":")]] <- p
    }
  }
  if (length(cond_plots) > 0) {
    p_inter <- patchwork::wrap_plots(cond_plots, ncol = 2) +
      patchwork::plot_annotation(
        title = "Driver interaction marginal effects",
        subtitle = sprintf("brms conditional_effects on %s response | run=%s",
                            names(driver_results)[1], RUN_LABEL),
        theme = theme_v2_pub(11))
    save_dual_v3(p_inter, "fig_v3_interaction_marginal",
                  width = 12, height = 8.5)
  }
}

message(sprintf("[stage-17] done. %d responses fitted.", length(driver_results)))
