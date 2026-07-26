#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   在显式控制不完全探测、调查强度、空间相关性与时间相关性后，
#   2000–2024 年中国鸟类多物种动态占域如何演化？（Q1 核心模型）
#
# Objective / 分析目标:
#   v4 模型：stMsPGOcc + AR1 时间 + NNGP 空间。
#   v4 关键改进相对 v3：
#     (1) 保存前显式 fit$z.samples <- NULL；fit$w.samples <- NULL
#         减少 60%+ 文件体积（C7 修复）
#     (2) qs::qsave(preset="fast") 替代 saveRDS，压缩快 ~5×（C7）
#     (3) 立刻 thin 一份 psi_samples_thinned_v4_*.qs，下游不读完整 fit
#     (4) 检测协变量统一 ~ log_events + log_duration + has_duration
#         has_duration=0 时 log_duration=0（missingness indicator，D4）
#     (5) 候选物种 fallback 抽到 utils_core::limit_candidate_species
#     (6) 单独写出 model_io_audit.csv 让 04b 知道 fit 路径与 thinned 路径
#
# Input data / 输入数据:
#   data/derived_v4/survey_history{GRID_TAG}_v4.rds
#   data/derived_v4/grid_environment{GRID_TAG}_v4.rds
#   data/derived_v4/detection_covariates{GRID_TAG}_v4.rds
#   （以上若不存在，按优先级 fallback v3 → v2）
#
# Main workflow / 主要流程:
#   1. 加载 survey / grid_env / det_cov
#   2. 构建 y array (species × site × period × secondary)
#   3. 准备 occ.covs / det.covs / coords
#   4. 运行 stMsPGOcc（pilot 20 sp 或 full 200 sp，单链或链并行）
#   5. 清理 fit slots（z/w）+ qs 保存 + thinned 保存
#   6. 写出运行摘要 + IO 审计
#
# Key assumptions / 关键假设:
#   - 服务器 RAM ≥ 256 GB 用于 200 sp × 4 chain
#   - 链并行模式：V4_CHAIN_ID=1..N 单链 + 04b 合并
#
# Main packages / 主要包:
#   spOccupancy, qs, readr, dplyr, tidyr, geosphere
#
# Output directory / 输出路径:
#   data/derived_v4/stMsPGOcc_fit_<run_label>[_chainK].qs
#   data/derived_v4/psi_samples_thinned_<run_label>[_chainK].qs
#   results_v4/table_model_summary_<run_label>.csv
#   results_v4/table_model_io_audit_<run_label>.csv
# ============================================================

gc()

suppressPackageStartupMessages({
  library(spOccupancy); library(readr); library(dplyr); library(tidyr)
  library(qs)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
source(file.path(CODE_V4, "utils_seeds.R"))
P <- ensure_v4_dirs()

log_time("04", "Starting stMsPGOcc v4 model fit")

# ── 1. 加载数据（v4 → v3 → v2 fallback） ─────────────────────────────
load_with_fallback <- function(stems, dir_chain) {
  for (d in dir_chain) {
    for (s in stems) {
      p_qs <- file.path(d, paste0(s, ".qs"))
      p_rds <- file.path(d, paste0(s, ".rds"))
      for (p in c(p_qs, p_rds)) {
        if (file.exists(p)) {
          obj <- safe_read(p, quiet = TRUE)
          if (!is.null(obj)) {
            message(sprintf("[04] Loaded %s from %s", s, p))
            return(obj)
          }
        }
      }
    }
  }
  NULL
}

survey <- load_with_fallback(
  stems = c(paste0("survey_history", GRID_TAG, "_v4"),
            paste0("survey_history", GRID_TAG, "_v3"),
            "survey_history_v3"),
  dir_chain = c(DIRS$derived, DIRS$v3_derived)
)
if (is.null(survey)) stop("[04] survey_history not found. Run 02 first.")

grid_env <- load_with_fallback(
  stems = c(paste0("grid_environment", GRID_TAG, "_v4"),
            paste0("grid_environment", GRID_TAG, "_v3"),
            "grid_environment_v3",
            "grid_environment_dynamic_occupancy"),
  dir_chain = c(DIRS$derived, DIRS$v3_derived, DIRS$v2_derived)
)
if (is.null(grid_env)) stop("[04] grid_environment not found. Run 03 first.")

det_cov <- load_with_fallback(
  stems = c(paste0("detection_covariates", GRID_TAG, "_v4"),
            paste0("detection_covariates", GRID_TAG, "_v3"),
            "detection_covariates_v3"),
  dir_chain = c(DIRS$derived, DIRS$v3_derived)
)

visit_effort  <- safe_read(file.path(DIRS$v2_derived, "visit_effort_2000_2024.rds"), quiet = TRUE)
species_visit <- safe_read(file.path(DIRS$v2_derived, "species_visit_2000_2024.rds"), quiet = TRUE)

# ── 2. 站点、物种、主期 ──────────────────────────────────────────────
sites      <- survey$sites
n_sites    <- length(sites)
periods    <- survey$periods
n_periods  <- length(periods)

MAX_N_SP <- as.integer(Sys.getenv("V4_MAX_SPECIES", "200"))
candidate_species <- limit_candidate_species(survey, max_n = MAX_N_SP)
n_sp <- length(candidate_species)
message(sprintf("[04] %d species × %d sites × %d periods", n_sp, n_sites, n_periods))

# ── 3. 构建 y array ───────────────────────────────────────────────────
use_v2_y <- FALSE
if (!is.null(species_visit) && !is.null(visit_effort)) {
  sites_ve <- sort(unique(visit_effort$grid_cell))
  sites_ge <- sort(unique(grid_env$grid_cell))
  overlap_ve <- length(intersect(sites_ve, sites))
  overlap_ge <- length(intersect(sites_ge, sites))
  message(sprintf("[04] Site overlap: visit_effort=%d/%d, grid_env=%d/%d, survey=%d",
                  overlap_ve, length(sites_ve), overlap_ge, length(sites_ge), n_sites))
  if (overlap_ve > n_sites * 0.5 && overlap_ge > n_sites * 0.5) use_v2_y <- TRUE
}

if (use_v2_y) {
  message("[04] Building y from v2 species_visit + visit_effort (good overlap)")
  sites_ve <- sort(unique(visit_effort$grid_cell))
  sites_ge <- sort(unique(grid_env$grid_cell))
  sites    <- sort(intersect(sites_ve, sites_ge))
  n_sites  <- length(sites)
  n_primary   <- n_periods
  n_secondary <- max(visit_effort$year_in_block, na.rm = TRUE)

  eff <- visit_effort |>
    filter(grid_cell %in% sites) |>
    mutate(site_index = match(grid_cell, sites))
  det <- species_visit |>
    filter(species %in% candidate_species, grid_cell %in% sites) |>
    mutate(species_index = match(species, candidate_species),
           site_index    = match(grid_cell, sites))

  y <- array(NA_integer_, dim = c(n_sp, n_sites, n_primary, n_secondary),
             dimnames = list(candidate_species, as.character(sites),
                             paste0("P", seq_len(n_primary)),
                             paste0("rep", seq_len(n_secondary))))

  visited_mask     <- array(FALSE, dim = c(n_sites, n_primary, n_secondary))
  log_events_arr   <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
  has_duration_arr <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
  log_duration_arr <- array(0,        dim = c(n_sites, n_primary, n_secondary))

  for (ii in seq_len(nrow(eff))) {
    rr <- eff[ii, ]
    if (rr$block_id > n_primary || rr$year_in_block > n_secondary) next
    visited_mask[rr$site_index, rr$block_id, rr$year_in_block] <- TRUE
    log_events_arr[rr$site_index, rr$block_id, rr$year_in_block] <- rr$log_events
    has_dur <- as.integer(!is.na(rr$mean_duration_min) && rr$mean_duration_min > 0)
    has_duration_arr[rr$site_index, rr$block_id, rr$year_in_block] <- has_dur
    if (has_dur == 1) {
      log_duration_arr[rr$site_index, rr$block_id, rr$year_in_block] <- log10(rr$mean_duration_min)
    }
  }

  for (ii in seq_len(nrow(det))) {
    rr <- det[ii, ]
    if (rr$species_index > n_sp || rr$site_index > n_sites) next
    if (rr$block_id > n_primary || rr$year_in_block > n_secondary) next
    y[rr$species_index, rr$site_index, rr$block_id, rr$year_in_block] <- rr$detected
  }

  visited_idx <- which(visited_mask, arr.ind = TRUE)
  for (ii in seq_len(nrow(visited_idx))) {
    j <- visited_idx[ii, 1]; t <- visited_idx[ii, 2]; k <- visited_idx[ii, 3]
    miss <- is.na(y[, j, t, k])
    if (any(miss)) y[miss, j, t, k] <- 0L
  }

  fill_arr_na <- function(arr, mask) {
    obs <- is.finite(arr)
    fill_val <- if (any(obs)) median(arr[obs], na.rm = TRUE) else 0
    arr[mask & !obs] <- fill_val
    arr
  }
  log_events_arr   <- fill_arr_na(log_events_arr, visited_mask)
  has_duration_arr <- fill_arr_na(has_duration_arr, visited_mask)
  has_duration_arr[!is.finite(has_duration_arr)] <- 0

} else {
  if (Sys.getenv("ALLOW_SINGLE_REPEAT", "0") != "1") {
    stop(paste(
      "Repeated-detection data are unavailable or misaligned.",
      "Refusing the rep=1 fallback because occupancy and detection are not separately identifiable.",
      "Rebuild/sync visit_effort and species_visit, or set ALLOW_SINGLE_REPEAT=1 only for debugging."
    ))
  }
  message("[04] Building y from v3 survey_history (simplified rep=1)")
  n_secondary <- 1L
  y <- array(0L, dim = c(n_sp, n_sites, n_periods, n_secondary),
             dimnames = list(candidate_species, as.character(sites),
                             paste0("P", seq_len(n_periods)), "rep1"))

  y_mat <- survey$y
  sp_idx <- match(candidate_species, rownames(y_mat))
  if (any(is.na(sp_idx))) {
    missing_sp <- candidate_species[is.na(sp_idx)]
    message(sprintf("[04] Warning: %d species not in survey$y (e.g. %s)",
                    length(missing_sp), paste(head(missing_sp, 3), collapse = ", ")))
    candidate_species <- candidate_species[!is.na(sp_idx)]
    sp_idx <- sp_idx[!is.na(sp_idx)]
    n_sp <- length(sp_idx)
    y <- array(0L, dim = c(n_sp, n_sites, n_periods, n_secondary),
               dimnames = list(candidate_species, as.character(sites),
                               paste0("P", seq_len(n_periods)), "rep1"))
  }

  for (t in seq_len(n_periods)) {
    col_start <- (t - 1) * n_sites + 1
    col_end <- t * n_sites
    if (col_end <= ncol(y_mat)) {
      y[, , t, 1] <- as.integer(y_mat[sp_idx, col_start:col_end, drop = FALSE])
    }
  }

  log_events_arr   <- array(0, dim = c(n_sites, n_periods, 1))
  has_duration_arr <- array(1, dim = c(n_sites, n_periods, 1))
  log_duration_arr <- array(0, dim = c(n_sites, n_periods, 1))
  if (!is.null(det_cov)) {
    for (t in seq_len(n_periods)) {
      p_data <- det_cov |> filter(period == paste0("P", t))
      log_events_arr[, t, 1] <- p_data$log_events[match(sites, p_data$grid_cell)]
      has_duration_arr[, t, 1] <- p_data$has_duration[match(sites, p_data$grid_cell)]
      if ("log_duration" %in% names(p_data)) {
        log_duration_arr[, t, 1] <- p_data$log_duration[match(sites, p_data$grid_cell)]
      }
    }
    log_events_arr[!is.finite(log_events_arr)] <- 0
    has_duration_arr[!is.finite(has_duration_arr)] <- 0
    log_duration_arr[!is.finite(log_duration_arr)] <- 0
  }
}

# ── 4. 占有率协变量 ──────────────────────────────────────────────────
standardize <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(0, length(x)))
  x[!is.finite(x)] <- median(x[is.finite(x)], na.rm = TRUE)
  as.numeric(scale(x))
}

grid_env_aligned <- grid_env |>
  filter(grid_cell %in% sites) |>
  arrange(match(grid_cell, sites))

occ_vars <- c("bio4", "bio7", "bio11", "bio13",
              "elev_mean", "elev_sd", "texture_shannon",
              "habitat_diversity_shannon",
              "hfi_mean", "landcover_built", "landcover_cropland",
              "centroid_lon", "centroid_lat")
occ_vars <- intersect(occ_vars, names(grid_env_aligned))
message("[04] Occupancy covariates: ", paste(occ_vars, collapse = ", "))

occ_covs <- setNames(lapply(occ_vars, function(v) standardize(grid_env_aligned[[v]])), occ_vars)
occ_covs$year_scaled <- matrix(
  rep(scale(seq_len(n_periods)), each = n_sites),
  nrow = n_sites, ncol = n_periods
)

det_covs <- list(
  log_events   = log_events_arr,
  log_duration = log_duration_arr,
  has_duration = has_duration_arr
)

coords_sf <- sf::st_as_sf(
  grid_env_aligned[, c("grid_cell", "centroid_lon", "centroid_lat")],
  coords = c("centroid_lon", "centroid_lat"), crs = 4326
)
coords <- sf::st_coordinates(sf::st_transform(
  coords_sf,
  "+proj=aea +lat_1=25 +lat_2=47 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
))[, 1:2, drop = FALSE]
grid_index <- seq_len(n_sites)

# ── 5. 模型数据列表 ─────────────────────────────────────────────────
data_list <- list(
  y          = y,
  occ.covs   = occ_covs,
  det.covs   = det_covs,
  coords     = coords,
  grid.index = grid_index
)

occ_formula <- reformulate(c(occ_vars, "year_scaled"))
det_formula <- ~ log_events + log_duration + has_duration

# ── 6. Pilot / Full / Chain parallel ─────────────────────────────────
is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

n_batch  <- if (is_pilot) PILOT_N_BATCH  else FULL_N_BATCH
n_burn   <- if (is_pilot) PILOT_N_BURN   else FULL_N_BURN
n_thin   <- if (is_pilot) PILOT_N_THIN   else FULL_N_THIN
n_chains <- if (is_pilot) PILOT_N_CHAINS else FULL_N_CHAINS
n_omp <- as.integer(Sys.getenv("V4_N_OMP_THREADS", "1"))

chain_id_raw <- Sys.getenv("V4_CHAIN_ID", "")
chain_id <- NA_integer_
if (nzchar(chain_id_raw)) {
  chain_id <- as.integer(chain_id_raw)
  n_chains <- 1L
  message(sprintf("[04] CHAIN PARALLEL mode: chain %d of %d", chain_id, FULL_N_CHAINS))
} else {
  message("[04] SEQUENTIAL mode: all chains in one process")
}

n_sp_pilot <- n_sp
if (is_pilot) {
  pilot_n <- min(PILOT_MAX_SP, n_sp)
  set_seeds("04_pilot_subsample")
  pilot_sp_idx <- sort(sample(n_sp, pilot_n))
  data_list$y <- data_list$y[pilot_sp_idx, , , , drop = FALSE]
  n_sp_pilot <- pilot_n
  message(sprintf("[04] PILOT mode: %d species", pilot_n))
}

n_factors_use <- min(5, max(1, n_sp_pilot - 1))

# ── 7. 起始值与先验 ─────────────────────────────────────────────────
z_init <- apply(data_list$y, c(1, 2, 3), function(a) as.integer(any(a == 1, na.rm = TRUE)))
if (length(dim(z_init)) != 3 || dim(z_init)[1] != n_sp_pilot ||
    dim(z_init)[2] != n_sites || dim(z_init)[3] != n_periods) {
  z_init <- array(as.integer(z_init), dim = c(n_sp_pilot, n_sites, n_periods))
}

chain_seed_offset <- if (!is.na(chain_id)) chain_id * 1000L else 0L
set_seeds(paste0("04_stMsPGOcc_chain", ifelse(is.na(chain_id), "0", chain_id)))

dmat_km <- as.matrix(dist(coords))
pairwise_km <- dmat_km[upper.tri(dmat_km)]
nn_km <- apply(replace(dmat_km, row(dmat_km) == col(dmat_km), Inf), 1, min)
min_range_km <- max(100, stats::median(nn_km[is.finite(nn_km)]))
max_range_km <- max(pairwise_km)
phi_bounds <- sort(c(3 / max_range_km, 3 / min_range_km))
phi_init <- 3 / stats::median(pairwise_km)
phi_init <- min(max(phi_init, phi_bounds[1] * 1.01), phi_bounds[2] * 0.99)
message(sprintf("[04] projected coords (km); phi %.6f, prior [%.6f, %.6f]",
                phi_init, phi_bounds[1], phi_bounds[2]))
rm(dmat_km, pairwise_km, nn_km); gc()

inits <- list(
  alpha.comm = 0, beta.comm = 0, beta = 0, alpha = 0,
  tau.sq.beta = 1, tau.sq.alpha = 1,
  z = z_init, sigma.sq = 1, phi = phi_init
)
priors <- list(
  beta.comm.normal  = list(mean = 0, var = 2.72),
  alpha.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.beta.ig    = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig   = list(a = 0.1, b = 0.1),
  sigma.sq.ig       = list(a = 2, b = 1),
  phi.unif          = list(a = phi_bounds[1], b = phi_bounds[2])
)
tuning <- list(phi = 0.5, sigma.sq = 1, rho = 0.2)

# ── 8. 运行 stMsPGOcc ────────────────────────────────────────────────
message(sprintf("[04] Running stMsPGOcc: %s", run_label))
message(sprintf("[04] MCMC: %d batch × 25, burn=%d, thin=%d, chains=%d",
                n_batch, n_burn, n_thin, n_chains))
message(sprintf("[04] Spatial: NNGP %d neighbors, %s covariance",
                N_NEIGHBORS, COV_MODEL))

t_start <- Sys.time()
fit <- stMsPGOcc(
  occ.formula   = occ_formula,
  det.formula   = det_formula,
  data          = data_list,
  inits         = inits,
  priors        = priors,
  tuning        = tuning,
  cov.model     = COV_MODEL,
  NNGP          = TRUE,
  n.neighbors   = N_NEIGHBORS,
  n.factors     = n_factors_use,
  n.batch       = n_batch,
  batch.length  = 25,
  accept.rate   = 0.43,
  n.omp.threads = n_omp,
  verbose       = TRUE,
  ar1           = TRUE,
  n.report      = max(1, floor(n_batch / 10)),
  n.burn        = n_burn,
  n.thin        = n_thin,
  n.chains      = n_chains
)
t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
message(sprintf("[04] Model fit completed in %.1f minutes", as.numeric(t_elapsed)))

# ── 9. 内存清理（C7） ────────────────────────────────────────────────
# z.samples / w.samples 占用 60%+，下游不需要 → 丢
fit$z.samples <- NULL
fit$w.samples <- NULL
gc()

# ── 10. 保存：完整 fit（qs） + thinned psi.samples ───────────────────
fit_stem <- if (!is.na(chain_id)) {
  paste0("stMsPGOcc_fit_", run_label, "_chain", chain_id)
} else {
  paste0("stMsPGOcc_fit_", run_label)
}
fit_path <- file.path(DIRS$derived, paste0(fit_stem, ".qs"))
qs_save_safe(fit, fit_path, preset = "fast")

# Thinned psi.samples（下游 05/06 只读这个）
psi <- fit$psi.samples
psi_dim <- dim(psi)
n_draws_total <- psi_dim[1]
n_draws_use   <- min(PSI_MAX_DRAWS, n_draws_total)
draw_idx <- round(seq(1, n_draws_total, length.out = n_draws_use))

if (length(psi_dim) >= 4) {
  psi_thinned <- psi[draw_idx, , , , drop = FALSE]
} else {
  psi_thinned <- psi[draw_idx, , , drop = FALSE]
}

thin_stem <- if (!is.na(chain_id)) {
  paste0("psi_samples_thinned_", run_label, "_chain", chain_id)
} else {
  paste0("psi_samples_thinned_", run_label)
}
thin_path <- file.path(DIRS$derived, paste0(thin_stem, ".qs"))
qs_save_safe(list(
  psi_samples_thinned = psi_thinned,
  species   = rownames(fit$y),
  sites     = colnames(fit$y),
  draw_idx  = draw_idx,
  psi_dim   = dim(psi_thinned),
  n_periods = if (length(psi_dim) >= 4) psi_dim[4] else 1L,
  run_label = run_label,
  chain_id  = chain_id
), thin_path)

# ── 11. WAIC（保护性 try） ──────────────────────────────────────────
waic_val <- tryCatch(spOccupancy::waicOcc(fit)$waic, error = function(e) NA_real_)
if (is.na(waic_val)) waic_val <- tryCatch(waic(fit)$waic, error = function(e) NA_real_)
message(sprintf("[04] WAIC = %.2f", waic_val))

# ── 12. 摘要 + IO 审计 ──────────────────────────────────────────────
write_csv(tibble(
  run_label    = run_label,
  chain_id     = chain_id,
  n_species    = if (is_pilot) n_sp_pilot else n_sp,
  n_sites      = n_sites,
  n_periods    = n_periods,
  cov_model    = COV_MODEL,
  n_neighbors  = N_NEIGHBORS,
  n_batch      = n_batch,
  n_burn       = n_burn,
  n_thin       = n_thin,
  n_chains     = n_chains,
  elapsed_mins = as.numeric(t_elapsed),
  waic         = waic_val,
  is_pilot     = is_pilot
), v4_file("results", paste0("table_model_summary_", run_label)))

write_csv(tibble(
  run_label = run_label,
  chain_id  = chain_id,
  fit_path  = fit_path,
  thin_path = thin_path,
  fit_size_mb  = round(file.info(fit_path)$size / 1e6, 1),
  thin_size_mb = round(file.info(thin_path)$size / 1e6, 1)
), v4_file("results", paste0("table_model_io_audit_", run_label,
                              if (!is.na(chain_id)) paste0("_chain", chain_id) else "")))

log_time("04", sprintf("DONE: %s, WAIC=%.2f, %.0f min",
                       run_label, waic_val, as.numeric(t_elapsed)))
