#!/usr/bin/env Rscript
## 16_diversity_v3_phylo_func.R
##
## 在 200 thinned psi draws × 1308 grids × 5 periods × 200 species 上
## 计算扩展多样性指标：
##   - Phylogenetic (picante): PSV, PSR, MPD/SES.MPD/NRI, MNTD/SES.MNTD/NTI
##   - Functional (fundiversity + 自写): FRic, FEve, FDiv, FDis, Rao Q, FSpe, FOri
##
## 设计：
##   - 用 utils_diversity.R 的 div_phylo_extended() / div_fd_extended()
##   - V2_POST_DRAWS_PHYLO=200 用全 200 draws 算 PSV/MPD/MNTD（快）
##   - V2_POST_DRAWS_FD=50    用 50 draws 算 FRic（凸包慢）+ 自写 FSpe/FOri
##   - future.apply::future_lapply 并行 grid×period
##
## Output:
##   data/derived_v3/diversity_v3_per_draw_<RUN_LABEL>.rds
##   results_v3/table_community_metrics_v3_<RUN_LABEL>.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ape); library(picante); library(fundiversity); library(future.apply)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_diversity.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")
N_DRAWS_PHYLO <- as.integer(Sys.getenv("V2_POST_DRAWS_PHYLO", "200"))
N_DRAWS_FD    <- as.integer(Sys.getenv("V2_POST_DRAWS_FD", "50"))
N_WORKERS     <- as.integer(Sys.getenv("V2_WORKERS", "6"))

message(sprintf("[stage-16] v3 diversity expansion | label=%s | phylo=%d FD=%d workers=%d",
                RUN_LABEL, N_DRAWS_PHYLO, N_DRAWS_FD, N_WORKERS))

## --- 1. 加载 psi 后验 + 元数据 -------------------------------------------

psi_obj <- readRDS(v2_file("derived",
                paste0("psi_samples_thinned_", RUN_LABEL), "rds"))
psi_arr <- psi_obj$psi_samples_thinned
n_draws_avail <- dim(psi_arr)[1]
n_sp    <- dim(psi_arr)[2]
n_site  <- dim(psi_arr)[3]
n_per   <- dim(psi_arr)[4]

primary_blocks <- read_csv(v2_file("results", "table_primary_5year_blocks"),
                            show_col_types = FALSE)
visit_effort <- readRDS(file.path(P$derived_v2, "visit_effort_2000_2024.rds"))
sites <- sort(unique(visit_effort$grid_cell))
cand_all <- read_csv(v2_file("results",
                "table_dynamic_occupancy_candidate_species_all"),
                show_col_types = FALSE)
candidate_species <- head(cand_all$species, n_sp)

## --- 2. 性状矩阵 + PCoA + 距离 ------------------------------------------

trait_ext <- readRDS(file.path(P$derived_v2, "trait_extended.rds"))
trait_vars <- c("body_mass_g", "clutch_size", "longevity_y",
                "maturity_y", "avonet_hwi", "avonet_range_size",
                "Habitat.Density")
trait_vars <- intersect(trait_vars, names(trait_ext))

trait_aligned <- tibble(species = candidate_species) |>
  left_join(trait_ext |> select(species, all_of(trait_vars)),
             by = "species") |>
  mutate(across(where(is.numeric),
                ~ if_else(.x > 0, log10(.x), .x)))
trait_mat <- as.matrix(trait_aligned[, trait_vars])
rownames(trait_mat) <- candidate_species
for (j in seq_len(ncol(trait_mat))) {
  m <- !is.finite(trait_mat[, j])
  if (any(m)) trait_mat[m, j] <- median(trait_mat[, j], na.rm = TRUE)
}
trait_mat_z <- scale(trait_mat)
trait_dist <- as.matrix(stats::dist(trait_mat_z))
pco <- ape::pcoa(stats::as.dist(trait_dist))
trait_pcoa <- pco$vectors[, seq_len(min(5, ncol(pco$vectors)))]
rownames(trait_pcoa) <- candidate_species

## --- 3. clootl 系统发育 -------------------------------------------------

phy <- NULL; phy_dist <- NULL; species_match <- NULL
phy_obj <- tryCatch({
    suppressMessages(library(clootl))
    data("clootl_data", package = "clootl", envir = globalenv())
    clootl::extractTree(species = candidate_species,
                         label_type = "scientific",
                         taxonomy_year = 2025, version = "1.6",
                         force = TRUE)
  }, error = function(e) {
    message("  clootl failed: ", conditionMessage(e)); NULL })
if (!is.null(phy_obj)) {
  phy <- phy_obj
  phy_dist <- ape::cophenetic.phylo(phy)
  species_match <- match(candidate_species, phy$tip.label)
  message(sprintf("  phylogeny: %d / %d species matched",
                  sum(!is.na(species_match)), length(candidate_species)))
}

## --- 4. 选择 draws 索引 ------------------------------------------------

draws_phylo <- if (n_draws_avail > N_DRAWS_PHYLO)
  round(seq(1, n_draws_avail, length.out = N_DRAWS_PHYLO)) else
  seq_len(n_draws_avail)
draws_fd <- if (n_draws_avail > N_DRAWS_FD)
  round(seq(1, n_draws_avail, length.out = N_DRAWS_FD)) else
  seq_len(n_draws_avail)

## --- 5. 主计算循环（按 (grid×period) 并行） -----------------------------

options(future.globals.maxSize = 4 * 1024^3)
plan(multisession, workers = N_WORKERS)

worker_one <- function(s, t, draws_idx, kind = c("phylo","fd")) {
  kind <- match.arg(kind)
  out <- vector("list", length(draws_idx))
  for (k in seq_along(draws_idx)) {
    d <- draws_idx[k]
    psi_v <- psi_arr[d, , s, t]
    res <- if (kind == "phylo")
             div_phylo_extended(psi_v, phy, species_match, phy_dist)
           else
             div_fd_extended(psi_v, trait_pcoa, trait_dist)
    out[[k]] <- c(draw_id = d, res)
  }
  do.call(rbind, out)
}

message(sprintf("[stage-16] phylo (%d draws) over %d sites x %d periods ...",
                length(draws_phylo), n_site, n_per))
t0 <- Sys.time()
phylo_res <- future_lapply(seq_len(n_site), function(s) {
  per_period <- lapply(seq_len(n_per), function(t) {
    M <- worker_one(s, t, draws_phylo, "phylo")
    cbind(grid_cell = sites[s], block_id = t, M)
  })
  do.call(rbind, per_period)
}, future.seed = TRUE)
phylo_mat <- do.call(rbind, phylo_res)
message(sprintf("  phylo done in %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

message(sprintf("[stage-16] FD (%d draws) over %d sites x %d periods ...",
                length(draws_fd), n_site, n_per))
t1 <- Sys.time()
fd_res <- future_lapply(seq_len(n_site), function(s) {
  per_period <- lapply(seq_len(n_per), function(t) {
    M <- worker_one(s, t, draws_fd, "fd")
    cbind(grid_cell = sites[s], block_id = t, M)
  })
  do.call(rbind, per_period)
}, future.seed = TRUE)
fd_mat <- do.call(rbind, fd_res)
message(sprintf("  FD done in %.1f min",
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))

plan(sequential)

## --- 6. 整理 + 落盘 -----------------------------------------------------

phylo_df <- as_tibble(phylo_mat) |>
  mutate(grid_cell = as.integer(grid_cell),
         block_id = as.integer(block_id),
         draw_id = as.integer(draw_id),
         block_label = primary_blocks$block_label[block_id])

fd_df <- as_tibble(fd_mat) |>
  mutate(grid_cell = as.integer(grid_cell),
         block_id = as.integer(block_id),
         draw_id = as.integer(draw_id),
         block_label = primary_blocks$block_label[block_id])

saveRDS(list(phylo = phylo_df, fd = fd_df),
        v3_file("derived",
                paste0("diversity_v3_per_draw_", RUN_LABEL), "rds"))

summarise_metric <- function(df, metric_cols) {
  long <- df |>
    pivot_longer(all_of(metric_cols), names_to = "metric",
                  values_to = "value")
  long |>
    group_by(grid_cell, block_label, block_id, metric) |>
    summarise(value_mean = mean(value, na.rm = TRUE),
              value_l95  = quantile(value, 0.025, na.rm = TRUE),
              value_u95  = quantile(value, 0.975, na.rm = TRUE),
              value_sd   = sd(value, na.rm = TRUE),
              .groups = "drop")
}

phylo_summary <- summarise_metric(phylo_df,
  c("psv","psr","mpd_ses","mntd_ses","nri","nti"))
fd_summary <- summarise_metric(fd_df,
  c("fric","feve","fdiv","fdis","fraoq","fspe","fori"))

combined <- bind_rows(phylo_summary, fd_summary)
write_csv(combined,
          v3_file("results",
                  paste0("table_community_metrics_v3_", RUN_LABEL)))

message(sprintf("[stage-16] done; %d rows in summary table",
                nrow(combined)))
