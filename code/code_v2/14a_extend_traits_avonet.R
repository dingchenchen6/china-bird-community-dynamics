#!/usr/bin/env Rscript
## 14a_extend_traits_avonet.R
##
## 从 AVONET 1 (BirdLife taxonomy) 抽取 Habitat / Habitat.Density / Trophic.Level /
## Trophic.Niche / Migration / Primary.Lifestyle，与已插补 trait 表合并，并派生：
##   - habitat_breadth      : 基于 occupancy-weighted WorldCover 多样性（数据驱动 Levins B）
##   - diet_specialization  : 1 - 候选物种池中该 Trophic.Niche 的相对频率（罕见 niche => 更专化）
##   - migration_class      : Resident / Partially / Migratory / Nomadic（factor）
##
## 输出：data/derived_v2/trait_extended.rds + results_v2/table_traits_extended.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_full_200sp_ar1")

avonet_path <- Sys.getenv("V2_AVONET_PATH",
  "/Users/dingchenchen/lucc/AVONET1_BirdLife.csv")
trait_path  <- file.path(P$project_root, "bird_grid_community_analysis",
                          "results", "table_species_traits_imputed.csv")

avonet <- read_csv(avonet_path, show_col_types = FALSE) |>
  rename(species = Species1) |>
  select(species, Habitat, Habitat.Density, Migration, Trophic.Level,
          Trophic.Niche, Primary.Lifestyle)
trait_imp <- read_csv(trait_path, show_col_types = FALSE)

cand_all <- read_csv(v2_file("results",
                "table_dynamic_occupancy_candidate_species_all"),
                show_col_types = FALSE)
cand_species <- cand_all$species

# 物种学名匹配率
n_match <- sum(cand_species %in% avonet$species)
message(sprintf("[14a] AVONET match: %d / %d candidate species", n_match,
                length(cand_species)))

trait_ext <- tibble(species = cand_species) |>
  left_join(trait_imp, by = "species") |>
  left_join(avonet, by = "species")

# diet_specialization：1 - rel.frequency in candidate pool
niche_freq <- trait_ext |>
  filter(!is.na(Trophic.Niche)) |>
  count(Trophic.Niche) |>
  mutate(rel = n / sum(n))
trait_ext <- trait_ext |>
  left_join(niche_freq |> select(Trophic.Niche, niche_rel = rel),
             by = "Trophic.Niche") |>
  mutate(diet_specialization = 1 - niche_rel)

# habitat_breadth：用 psi posterior mean 加权 WorldCover 多样性（Shannon）
# 仅在 stage4 输出存在时计算；否则置 NA 并标 fallback
psi_path <- v2_file("derived",
                    paste0("psi_samples_thinned_", RUN_LABEL), "rds")
trait_ext$habitat_breadth_data <- NA_real_

if (file.exists(psi_path)) {
  psi_obj <- readRDS(psi_path)
  psi_arr <- psi_obj$psi_samples_thinned
  psi_mean_grid <- apply(psi_arr, c(2, 3), mean)   # sp × site avg over draws/periods
  visit_effort <- readRDS(file.path(P$derived_v2, "visit_effort_2000_2024.rds"))
  sites <- sort(unique(visit_effort$grid_cell))
  grid_env <- readRDS(file.path(P$derived_v2,
                                "grid_environment_dynamic_occupancy.rds")) |>
    filter(grid_cell %in% sites) |>
    arrange(match(grid_cell, sites))
  lc_cols <- c("landcover_trees", "landcover_cropland", "landcover_built",
               "landcover_shrubs", "landcover_grassland", "landcover_water")
  lc_mat <- as.matrix(grid_env[, lc_cols, drop = FALSE])
  lc_mat[!is.finite(lc_mat)] <- 0
  # 按 site psi 加权汇总每物种的 WorldCover 占比
  for (i in seq_len(nrow(psi_mean_grid))) {
    w <- psi_mean_grid[i, ]
    if (sum(w, na.rm = TRUE) <= 0) next
    p <- colSums(lc_mat * w, na.rm = TRUE)
    if (sum(p) <= 0) next
    p <- p / sum(p)
    p <- p[p > 0]
    trait_ext$habitat_breadth_data[i] <- -sum(p * log(p))
  }
}

# habitat_breadth_avonet：把 AVONET Habitat 转 ordinal "openness"（粗代理）
hab_to_open <- c(Forest = 1, Woodland = 1.5, Shrubland = 2, Grassland = 2.5,
                  Wetland = 3, Riverine = 3, Coastal = 3,
                  Desert = 3.5, Marine = 4, Rock = 4, `Human Modified` = 4)
trait_ext <- trait_ext |>
  mutate(habitat_openness_avonet = unname(hab_to_open[Habitat]),
         migration_class = factor(Migration,
                                    levels = c("1", "2", "3", "Resident",
                                                "Partial", "Migratory")))
# Migration in AVONET is 1=Resident, 2=Partial, 3=Migratory (numeric); coerce
trait_ext <- trait_ext |>
  mutate(migration_score = case_when(
    Migration == 1 | Migration == "1" | Migration == "Resident"   ~ 1L,
    Migration == 2 | Migration == "2" | Migration == "Partial"    ~ 2L,
    Migration == 3 | Migration == "3" | Migration == "Migratory"  ~ 3L,
    TRUE ~ NA_integer_
  ))

# 用列中位数 fill 数值列（保持下游兼容）
num_fill_cols <- c("diet_specialization", "habitat_breadth_data",
                    "habitat_openness_avonet", "migration_score",
                    "Habitat.Density")
for (cc in num_fill_cols) {
  if (cc %in% names(trait_ext)) {
    v <- as.numeric(trait_ext[[cc]])
    if (any(is.na(v))) v[is.na(v)] <- median(v, na.rm = TRUE)
    trait_ext[[cc]] <- v
  }
}

saveRDS(trait_ext, file.path(P$derived_v2, "trait_extended.rds"))
write_csv(trait_ext, v2_file("results", "table_traits_extended"))
message(sprintf("[14a] trait_extended: %d species x %d cols (data/derived_v2/trait_extended.rds)",
                nrow(trait_ext), ncol(trait_ext)))
