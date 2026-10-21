#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题: 驱动因子对校正物种数趋势的解释力，在空间外推时还剩多少？
#   How much explanatory power survives spatially blocked validation?
# Objective / 分析目标: 以空间分块交叉验证（blockCV）替代样本内随机森林重要性，报告样本外 R²
#   与按预测变量组汇总的置换重要性；同时报告协变量是否完整。（BLOCKERS Gate 4）
# Input / 输入: table_community_metrics_with_cri_<RUN_LABEL>_extended.csv；grid_environment_v3.rds（S0 重建）；
#   climate_change_v3.rds、hfi_change_v3.rds、landuse_change_v3.rds（若存在）
# Output / 输出: results_v3/table_rf_spatial_cv_<RUN_LABEL>.csv、table_rf_group_importance_<RUN_LABEL>.csv
# Key assumptions / 关键假设: 分块边长 500 km、5 折；组重要性 = 组内变量置换重要性之和。
# Packages / 主要包: ranger, blockCV, sf
# ============================================================
suppressPackageStartupMessages({ library(ranger); library(blockCV); library(sf) })
CODE <- Sys.getenv("V3_CODE_DIR", "code_v4_refit")
source(file.path(CODE, "00_config.R")); source(file.path(CODE, "utils_paths.R")); set.seed(RF_SEED)
cm <- read.csv(v3_file("results", paste0("table_community_metrics_with_cri_", RUN_LABEL, "_extended")))
cm <- cm[cm$metric == "corrected_richness", ]; cm$p <- as.integer(sub("P", "", cm$period))
y <- sapply(split(cm, cm$grid_cell), function(d) coef(lm(value_mean ~ p, d))[2])
d <- data.frame(grid_cell = as.integer(names(y)), trend = y)
rd <- function(b) { p <- v3_file("derived", b, "rds"); if (file.exists(p)) as.data.frame(readRDS(p)) else NULL }
env <- rd("grid_environment_v3"); d <- merge(d, env, by = "grid_cell")
for (b in c("climate_change_v3", "hfi_change_v3", "landuse_change_v3")) {
  x <- rd(b); if (!is.null(x)) d <- merge(d, x[, setdiff(names(x), c("centroid_lon", "centroid_lat")), drop = FALSE], by = "grid_cell", all.x = TRUE)
}
groups <- list(
  baseline = intersect(c("centroid_lon", "centroid_lat", "elev_mean", "bio11"), names(d)),
  climate_change = intersect(c("delta_t_mean", "delta_t_extreme"), names(d)),
  landuse_change = grep("^delta_(forest|impervious|natural)", names(d), value = TRUE),
  human_pressure_change = intersect("delta_hfi", names(d)))
X <- unlist(groups); d <- d[complete.cases(d[, c("trend", X)]), ]
message(sprintf("[09b] 完整样本 %d 格；变量 %d 个", nrow(d), length(X)))
pts <- st_as_sf(d, coords = c("centroid_lon", "centroid_lat"), crs = 4326, remove = FALSE)
cv <- cv_spatial(x = pts, size = 500000, k = 5, selection = "random", iteration = 50, progress = FALSE, plot = FALSE)
pred <- rep(NA_real_, nrow(d))
for (k in seq_along(cv$folds_list)) {
  tr <- cv$folds_list[[k]][[1]]; te <- cv$folds_list[[k]][[2]]
  m <- ranger(trend ~ ., d[tr, c("trend", X)], num.trees = 500, seed = RF_SEED)
  pred[te] <- predict(m, d[te, c("trend", X)])$predictions
}
r2_cv <- 1 - sum((d$trend - pred)^2) / sum((d$trend - mean(d$trend))^2)
full <- ranger(trend ~ ., d[, c("trend", X)], num.trees = 1000, importance = "permutation", seed = RF_SEED)
write.csv(data.frame(n_cells = nrow(d), r2_in_sample_oob = round(full$r.squared, 3), r2_spatial_block_cv = round(r2_cv, 3)),
          v3_file("results", paste0("table_rf_spatial_cv_", RUN_LABEL)), row.names = FALSE)
imp <- full$variable.importance
gi <- data.frame(group = names(groups), importance = sapply(groups, function(v) sum(imp[v])))
write.csv(gi, v3_file("results", paste0("table_rf_group_importance_", RUN_LABEL)), row.names = FALSE)
message(sprintf("[09b] OOB R² = %.3f；空间分块交叉验证 R² = %.3f", full$r.squared, r2_cv))
