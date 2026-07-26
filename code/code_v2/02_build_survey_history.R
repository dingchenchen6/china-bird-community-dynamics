#!/usr/bin/env Rscript
## 02_build_survey_history.R
##
## 阶段 2：基于去重后事件构建动态占域所需的调查史。
## 输入：data/derived_v2/combined_events_merged_dedup_2000_2025.rds
##       (来自 01_merge_birdwatch_ebird.R)
## 输出：
##   data/derived_v2/visit_effort_2000_2024.rds
##   data/derived_v2/species_visit_2000_2024.rds
##   results_v2/table_primary_5year_blocks.csv
##   results_v2/table_dynamic_occupancy_candidate_species_all.csv
##   results_v2/table_species_detection_coverage.csv
##   results_v2/table_v1_vs_v2_coverage_comparison.csv
##   results_v2/table_survey_coverage_by_year.csv
##   results_v2/table_survey_coverage_by_block.csv
##   figures_v2/fig_audit_coverage_by_year_v2.{png,pdf}
##   figures_v2/fig_audit_coverage_by_block_v2.{png,pdf}
##   figures_v2/fig_audit_grid_coverage_map_v2.{png,pdf}

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(sf)
  library(ggplot2)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_spatial.R"))
source(file.path(CODE_V2, "utils_mapping.R"))

P <- ensure_v2_dirs()

ANALYSIS_YR_LO <- 2000L
ANALYSIS_YR_HI <- 2024L
BLOCK_SIZE     <- 5L
MIN_DETECTIONS <- as.integer(Sys.getenv("V2_MIN_DETECTIONS", "80"))
MIN_BLOCKS     <- as.integer(Sys.getenv("V2_MIN_BLOCKS", "2"))
MIN_GRIDS      <- as.integer(Sys.getenv("V2_MIN_GRIDS", "10"))

message("[stage-2] Loading dedup'd events")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                                "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]
message(sprintf("  events 2000-2024 after dedup: %s",
                format(nrow(events), big.mark=",")))

## --- 1. 5 年主期表 ---------------------------------------------------------

block_starts <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI - BLOCK_SIZE + 1, by = BLOCK_SIZE)
primary_blocks <- tibble(
  block_id    = seq_along(block_starts),
  block_start = block_starts,
  block_end   = block_starts + BLOCK_SIZE - 1L,
  block_label = sprintf("%d-%d", block_starts, block_starts + BLOCK_SIZE - 1L)
)
write_csv(primary_blocks, v2_file("results", "table_primary_5year_blocks"))

block_idx_for_year <- function(y) ((y - ANALYSIS_YR_LO) %/% BLOCK_SIZE) + 1L
events[, block_id := block_idx_for_year(year)]
events[, year_in_block := year - block_starts[block_id] + 1L]

## --- 2. 100km 网格 + 把事件分配到网格 -------------------------------------

grid_cache <- file.path(P$derived_v2, "china_grid_100km_v2.rds")
if (!file.exists(grid_cache)) {
  grid_sf <- load_or_build_v2_grid(P$china_boundary, resolution_m = 1e5,
                                   cache_path = grid_cache)
} else {
  grid_sf <- readRDS(grid_cache)
}
message(sprintf("  China 100 km grid: %d cells", nrow(grid_sf)))

assign_events_to_grid <- function(ev_dt, grid_sf) {
  pts <- sf::st_as_sf(
    ev_dt[, .(.I, longitude, latitude)],
    coords = c("longitude", "latitude"),
    crs = 4326, remove = FALSE
  )
  joined <- sf::st_join(pts, grid_sf[, "grid_cell"],
                        join = sf::st_intersects, left = TRUE)
  ev_dt[, grid_cell := joined$grid_cell]
  ev_dt[!is.na(grid_cell)]
}

message("  assigning events to grid cells (sf with planar geometry)")
old_s2 <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)
on.exit(sf::sf_use_s2(old_s2), add = TRUE)
grid_sf <- sf::st_make_valid(grid_sf)
events <- assign_events_to_grid(events, grid_sf)
message(sprintf("  events with grid_cell: %s",
                format(nrow(events), big.mark=",")))

## --- 3. visit_effort：(grid × year) 访问强度 ------------------------------

events[, event_visit_id := paste(grid_cell,
                                  sprintf("%04d-%02d-%02d", year, month, day),
                                  tolower(coalesce(username, "anon")), sep = "|")]
events[, has_dur := !is.na(duration_min)]

visit_effort <- events[, .(
  n_events           = uniqueN(event_visit_id),
  n_observers        = uniqueN(tolower(coalesce(username, "anon"))),
  n_unique_dates     = uniqueN(sprintf("%04d-%02d-%02d", year, month, day)),
  mean_duration_min  = if (any(has_dur)) mean(duration_min[has_dur]) else NA_real_
), by = .(grid_cell, year, block_id, year_in_block)]

visit_effort[, `:=`(
  log_events    = log1p(n_events),
  log_observers = log1p(n_observers)
)]

setorder(visit_effort, grid_cell, year)
saveRDS(as_tibble(visit_effort), file.path(P$derived_v2, "visit_effort_2000_2024.rds"))

## --- 4. species_visit：(grid × year × species) 检测 -----------------------

species_visit <- events[, .(
  detected           = 1L,
  n_detection_events = uniqueN(event_visit_id)
), by = .(grid_cell, year, block_id, year_in_block, species)]

setorder(species_visit, species, grid_cell, year)
saveRDS(as_tibble(species_visit), file.path(P$derived_v2,
                                            "species_visit_2000_2024.rds"))

## --- 5. 物种检测覆盖 + 候选筛选 -------------------------------------------

species_coverage <- species_visit[, .(
  n_grid_year_detections = .N,
  n_blocks_detected      = uniqueN(block_id),
  n_grids_detected       = uniqueN(grid_cell),
  first_year             = min(year),
  last_year              = max(year)
), by = species][order(-n_grid_year_detections)]

candidate_species_all <- species_coverage[
  n_grid_year_detections >= MIN_DETECTIONS &
  n_blocks_detected      >= MIN_BLOCKS &
  n_grids_detected       >= MIN_GRIDS
][, candidate_rank := .I]

write_csv(species_coverage,      v2_file("results", "table_species_detection_coverage"))
write_csv(candidate_species_all, v2_file("results",
                                  "table_dynamic_occupancy_candidate_species_all"))

message(sprintf("  species observed: %d", nrow(species_coverage)))
message(sprintf("  candidate (≥%d det, ≥%d blocks, ≥%d grids): %d",
                MIN_DETECTIONS, MIN_BLOCKS, MIN_GRIDS,
                nrow(candidate_species_all)))

## --- 6. v1 vs v2 覆盖对比 -------------------------------------------------

v1_cov_path <- file.path(P$results_v1, "table_species_detection_coverage.csv")
if (file.exists(v1_cov_path)) {
  v1_cov <- as.data.table(read_csv(v1_cov_path, show_col_types = FALSE))
  comparison <- merge(
    species_coverage[, .(species,
                          v2_n_detections = n_grid_year_detections,
                          v2_n_grids      = n_grids_detected,
                          v2_n_blocks     = n_blocks_detected)],
    v1_cov[, .(species,
                v1_n_detections = n_grid_year_detections,
                v1_n_grids      = n_grids_detected,
                v1_n_blocks     = n_blocks_detected)],
    by = "species", all = TRUE
  )
  comparison[, `:=`(
    delta_n_detections   = v2_n_detections - v1_n_detections,
    pct_change_detections = 100 * (v2_n_detections - v1_n_detections) /
                                  v1_n_detections
  )]
  setorder(comparison, -v2_n_detections)
  write_csv(as_tibble(comparison),
            v2_file("results", "table_v1_vs_v2_coverage_comparison"))
}

## --- 7. 调查覆盖审计：年 / 块 / 空间 --------------------------------------

cov_by_year <- visit_effort[, .(
  n_grids_visited = uniqueN(grid_cell),
  n_events_total  = sum(n_events),
  n_observers     = sum(n_observers),
  n_unique_dates  = sum(n_unique_dates)
), by = year][order(year)]
write_csv(as_tibble(cov_by_year),
          v2_file("results", "table_survey_coverage_by_year"))

cov_by_block <- visit_effort[, .(
  n_grids_visited   = uniqueN(grid_cell),
  n_events_total    = sum(n_events),
  n_observers       = sum(n_observers),
  mean_duration_min = mean(mean_duration_min, na.rm = TRUE)
), by = block_id][order(block_id)] |>
  merge(primary_blocks, by = "block_id")
write_csv(as_tibble(cov_by_block),
          v2_file("results", "table_survey_coverage_by_block"))

## --- 8. 审计图 -------------------------------------------------------------

scale_factor <- max(cov_by_year$n_events_total) / 1000 /
                max(cov_by_year$n_grids_visited)
plot_year_coverage <- ggplot(cov_by_year, aes(year)) +
  geom_col(aes(y = n_events_total / 1000), fill = "#0E5A78", alpha = 0.85) +
  geom_line(aes(y = n_grids_visited * scale_factor),
            colour = "#C9784A", linewidth = 0.9) +
  geom_point(aes(y = n_grids_visited * scale_factor),
             colour = "#C9784A", size = 1.7) +
  scale_y_continuous(
    name = "Visit events (thousands; bars)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Grids visited (line)")
  ) +
  scale_x_continuous(breaks = seq(2000, 2024, 4)) +
  labs(
    title    = "Survey coverage across years (merged + dedup)",
    subtitle = "Bars: # unique visit events; Line: # 100 km grids visited",
    x = NULL,
    caption  = sprintf("Total: %s events across %d grids; %s candidate species",
                       format(sum(cov_by_year$n_events_total), big.mark=","),
                       max(cov_by_year$n_grids_visited),
                       format(nrow(candidate_species_all), big.mark=","))
  ) +
  theme_v2_pub(11)
save_dual(plot_year_coverage, "fig_audit_coverage_by_year_v2",
          width = 10, height = 5.4)

cov_by_block_long <- as_tibble(cov_by_block) |>
  select(block_label, n_grids_visited, n_events_total, n_observers) |>
  pivot_longer(-block_label, names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
    n_grids_visited = "Grids visited",
    n_events_total  = "Events (total)",
    n_observers     = "Observer-visits"
  ))
plot_block_coverage <- ggplot(cov_by_block_long,
                              aes(block_label, value, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = V2_PALETTES$qualitative[1:3], guide = "none") +
  labs(
    title    = "Coverage by 5-year primary period",
    subtitle = "Merged + deduplicated event base",
    x = NULL, y = NULL
  ) +
  theme_v2_pub(11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
save_dual(plot_block_coverage, "fig_audit_coverage_by_block_v2",
          width = 12, height = 5)

grid_year_count <- visit_effort[, .(n_years_with_data = uniqueN(year)),
                                 by = grid_cell]
china_layers <- load_china_layers(P$china_boundary, P$province_line)
grid_cov_sf <- grid_sf |>
  left_join(as_tibble(grid_year_count), by = "grid_cell") |>
  mutate(n_years_with_data = ifelse(is.na(n_years_with_data), 0, n_years_with_data))

plot_grid_map <- ggplot(grid_cov_sf) +
  geom_sf(aes(fill = n_years_with_data),
          colour = scales::alpha("white", 0.10), linewidth = 0.05) +
  china_map_layers(china_layers) +
  scale_fill_v2_sequential(name = "Years with ≥1 visit (out of 25)",
                           palette = "lajolla", direction = 1,
                           limits = c(0, 25)) +
  v2_china_coord() +
  theme_v2_map(11) +
  labs(
    title    = "Per-grid temporal coverage (merged + deduplicated)",
    subtitle = "Grid cells coloured by # of years (2000-2024) with ≥1 visit event",
    caption  = sprintf("%d / %d grids have ≥1 visit during 2000-2024",
                       sum(grid_cov_sf$n_years_with_data > 0),
                       nrow(grid_cov_sf))
  )
save_dual(plot_grid_map, "fig_audit_grid_coverage_map_v2",
          width = 9.4, height = 6.6)

message("[stage-2] Survey history rebuilt. Outputs ready.")
