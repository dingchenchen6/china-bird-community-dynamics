#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   占域模型的格点环境协变量是否正确对应到各自的格点？
#   Are grid-level environmental covariates attached to the right cells?
# Objective / 分析目标:
#   重建 1,247 个 100 km 格点的全部环境协变量，修复 v1 提取的行号错配缺陷。
#   Rebuild all grid covariates, fixing the v1 row-index mismatch.
# Root cause / 根因（2026-10-20 核实）:
#   v1 `code/prepare_dynamic_occupancy_environment.R` 中 extract_mean_tbl() 把
#   terra::extract() 返回的 ID（多边形行号 1..1739）改名为 grid_cell（栅格单元编号 54..3665）
#   再 left_join。结果：编号 > 1739 的格点（全部南方格点）缺失；编号 ≤ 1739 的格点拿到了
#   别的格点的值（与正确值相关 −0.58）。v2/v3 环境表原样继承了该错误。
# Input / 输入:
#   - 格网: data/derived_v2/china_grid_100km_v2.rds（与 v3 格点编号、中心坐标完全一致）
#   - 站点: data/derived_v3/survey_history_v3.rds 的 $sites（1,247 格）
#   - 栅格: WorldClim 2.1 10'（bio1–19、海拔）、EarthEnv 纹理、NPP、WorldCover 30"、HFI 2000–2020
# Workflow / 流程: 读格网 → 按 grid_cell 显式对接提取（exactextractr 面积加权）→ 派生生境指标
#   → 缺失率与旧值比对 → 输出新环境表与 QA 报告
# Output / 输出（写入 BIRD_OUTPUT_ROOT，默认 <project>/refit_v4）:
#   data/derived_v3/grid_environment_v3.rds, results_v3/table_grid_environment_v3.csv,
#   results_v3/QA_grid_environment_rebuild.csv
# Key assumptions / 关键假设: 均值按格点与栅格单元的重叠面积加权；海拔标准差为格内 10' 单元的
#   加权标准差（与 v1 口径一致，分辨率较粗）。
# Packages / 主要包: sf, terra, exactextractr, dplyr
# ============================================================
suppressPackageStartupMessages({ library(sf); library(terra); library(exactextractr); library(dplyr) })
sf_use_s2(FALSE)   # 个别格点多边形在球面引擎下自相交；用平面几何，与 v1 口径一致 / planar geometry as in v1

PROJ <- Sys.getenv("BIRD_PROJECT_ROOT", path.expand("~/Documents/New project/bird_dynamic_occupancy_analysis"))
OUT  <- Sys.getenv("BIRD_OUTPUT_ROOT", file.path(PROJ, "refit_v4"))
# 栅格根目录可由环境变量覆盖（服务器上路径不同）/ raster roots, overridable on the server
CLIM <- Sys.getenv("REFIT_WORLDCLIM_DIR", file.path(dirname(PROJ), "bird_full_community_analysis/data/external/climate/wc2.1_10m"))
EXT  <- Sys.getenv("REFIT_ENV_EXTERNAL_DIR", file.path(dirname(PROJ), "bird_grid_community_analysis/data/external"))
HFI  <- Sys.getenv("REFIT_HFI_DIR", file.path(PROJ, "data/external/hfi"))
dir.create(file.path(OUT, "data/derived_v3"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "results_v3"), recursive = TRUE, showWarnings = FALSE)

# ── 1. 格网与站点 / grid and sites ─────────────────────────────────────
grid <- readRDS(file.path(PROJ, "data/derived_v2/china_grid_100km_v2.rds"))
sites <- readRDS(file.path(PROJ, "data/derived_v3/survey_history_v3.rds"))$sites
grid <- grid[match(as.integer(sites), as.integer(grid$grid_cell)), ]
stopifnot(nrow(grid) == length(sites), !anyNA(grid$grid_cell))
message(sprintf("[S0] %d grid cells aligned to survey sites", nrow(grid)))

# 按 grid_cell 显式对接的提取函数：结果行序与 grid 一致，并带上 grid_cell 列
# Extraction keyed explicitly by grid_cell (row order follows `grid`)
ext <- function(r, fun = "mean", prefix = NULL) {
  g <- st_transform(grid, crs(r))
  v <- exact_extract(r, g, fun = fun, progress = FALSE)
  v <- as.data.frame(v)
  names(v) <- if (nlyr(r) == 1) names(r) else names(r)   # exactextractr 给 "mean.<layer>"；统一回层名
  if (!is.null(prefix)) names(v) <- paste0(prefix, names(v))
  data.frame(grid_cell = as.integer(grid$grid_cell), v, check.names = FALSE)
}
rx <- function(files, nm) { r <- rast(files); names(r) <- nm; r }

# ── 2. 各类栅格 / rasters ──────────────────────────────────────────────
bio  <- rx(file.path(CLIM, sprintf("wc2.1_10m_bio_%d.tif", 1:19)), paste0("bio", 1:19))
elev <- rx(file.path(CLIM, "wc2.1_10m_elev.tif"), "elev")
tex  <- rx(file.path(EXT, "earthenv_texture", c("shannon_01_05_5km_uint16.tif", "Entropy_01_05_5km_uint16.tif", "Contrast_01_05_5km_uint32.tif")),
           c("texture_shannon", "texture_entropy", "texture_contrast"))
npp  <- rx(file.path(EXT, "npp", "npp_mean_china.tif"), "npp_mean")
lcn  <- c("trees", "cropland", "built", "shrubs", "grassland", "water")
lc   <- rx(file.path(EXT, "landcover/landuse", sprintf("WorldCover_%s_30s.tif", lcn)), paste0("landcover_", lcn))
hfi_files <- sort(list.files(HFI, pattern = "^hfp20[0-9]{2}\\.tif$", full.names = TRUE))
hfi  <- rx(hfi_files, paste0("hfi_", sub("^hfp(20[0-9]{2})\\.tif$", "\\1", basename(hfi_files))))

env <- ext(bio) |>
  left_join(ext(elev) |> rename(elev_mean = elev), by = "grid_cell") |>
  left_join(ext(elev, fun = "stdev") |> rename(elev_sd = elev), by = "grid_cell") |>
  left_join(ext(tex), by = "grid_cell") |>
  left_join(ext(npp), by = "grid_cell") |>
  left_join(ext(lc), by = "grid_cell") |>
  left_join(ext(hfi), by = "grid_cell")

# ── 3. 派生指标（口径同 v1）/ derived metrics, same definitions as v1 ──────
lc_cols <- paste0("landcover_", lcn)
M <- as.matrix(env[, lc_cols]); M[!is.finite(M)] <- 0
P <- M / pmax(rowSums(M), 1e-12)
sh <- apply(P, 1, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else -sum(x * log(x)) })
ri <- rowSums(P > 0.01)
hcols <- grep("^hfi_20", names(env), value = TRUE)
env <- env |> mutate(
  habitat_diversity_shannon  = sh,
  habitat_diversity_richness = ri,
  habitat_diversity_evenness = ifelse(ri > 1, sh / log(ri), 0),
  habitat_dominance          = apply(P, 1, max),
  habitat_diversity_simpson  = 1 - rowSums(P^2),
  natural_landcover_fraction = rowSums(M[, c("landcover_trees", "landcover_shrubs", "landcover_grassland", "landcover_water")]),
  human_modified_fraction    = rowSums(M[, c("landcover_cropland", "landcover_built")]),
  hfi_mean = rowMeans(across(all_of(hcols)), na.rm = TRUE),
  hfi_sd   = apply(across(all_of(hcols)), 1, sd, na.rm = TRUE))
cen <- st_coordinates(suppressWarnings(st_centroid(st_transform(grid, 4326))))
env <- data.frame(grid_cell = env$grid_cell, grid_id = grid$grid_id, centroid_lon = cen[, 1], centroid_lat = cen[, 2],
                  area_km2 = as.numeric(st_area(st_transform(grid, "+proj=aea +lat_1=25 +lat_2=47 +lon_0=105 +datum=WGS84 +units=m"))) / 1e6, env[, -1], year_scaled = 0)

# ── 4. 质检：缺失率与旧值比对 / QA: missingness and agreement with the old table ──
old <- as.data.frame(readRDS(file.path(PROJ, "data/derived_v3/grid_environment_v3.rds")))
vars <- setdiff(names(env), c("grid_cell", "grid_id", "year_scaled"))
qa <- do.call(rbind, lapply(vars, function(v) {
  o <- if (v %in% names(old)) as.numeric(unlist(old[[v]]))[match(env$grid_cell, as.integer(unlist(old$grid_cell)))] else rep(NA_real_, nrow(env))
  n <- as.numeric(env[[v]]); k <- is.finite(o) & is.finite(n)
  data.frame(variable = v, n_missing_new = sum(!is.finite(n)), n_missing_old = sum(!is.finite(o)),
             n_both = sum(k), cor_old_new = if (sum(k) > 2) round(cor(o[k], n[k]), 3) else NA,
             share_identical = if (sum(k)) round(mean(abs(o[k] - n[k]) <= 1e-6 * pmax(1, abs(n[k]))), 3) else NA)
}))
write.csv(qa, file.path(OUT, "results_v3/QA_grid_environment_rebuild.csv"), row.names = FALSE)
print(qa[qa$variable %in% c("bio4", "bio7", "bio11", "bio13", "elev_mean", "elev_sd", "texture_shannon",
                           "habitat_diversity_shannon", "hfi_mean", "landcover_built", "landcover_cropland",
                           "centroid_lon", "centroid_lat"), ], row.names = FALSE)

model_vars <- c("bio4", "bio7", "bio11", "bio13", "elev_mean", "elev_sd", "texture_shannon",
                "habitat_diversity_shannon", "hfi_mean", "landcover_built", "landcover_cropland")
miss <- sapply(model_vars, function(v) mean(!is.finite(env[[v]])))
if (any(miss > 0.02)) stop("[S0] 模型协变量缺失率超过 2%：", paste(names(miss)[miss > 0.02], collapse = ", "))
saveRDS(env, file.path(OUT, "data/derived_v3/grid_environment_v3.rds"))
write.csv(env, file.path(OUT, "results_v3/table_grid_environment_v3.csv"), row.names = FALSE)
message(sprintf("[S0] 完成：%d 格 × %d 变量；模型协变量最大缺失率 %.1f%%", nrow(env), ncol(env), 100 * max(miss)))
