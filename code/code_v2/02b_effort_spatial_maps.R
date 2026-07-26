#!/usr/bin/env Rscript
## 02b_effort_spatial_maps.R
##
## 用合并去重后的事件数据，按 (grid × 5 年主期 × source) 输出调查努力时空图。
## 关键产出：
##   - 网格 × 主期：events / observers / unique-dates / mean duration  地图
##   - eBird 与中国观鸟记录平台分别 vs 合并后 三套地图
##   - 双源差值图（中国观鸟 - eBird）显示哪个源在哪里贡献更多
##   - 时间序列：年度 events / observers / sources 共存
##
## 输入：data/derived_v2/combined_events_merged_dedup_2000_2025.rds
## 输出：figures_v2/fig_effort_*.{png,pdf}
##         results_v2/table_effort_by_grid_period_source.csv

suppressPackageStartupMessages({
  library(data.table)
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(forcats)
  library(ggplot2); library(patchwork); library(sf)
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

message("[stage-2b] loading dedup'd events + assigning to grid")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                          "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]

block_starts <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI - BLOCK_SIZE + 1, by = BLOCK_SIZE)
events[, block_id := ((year - ANALYSIS_YR_LO) %/% BLOCK_SIZE) + 1L]
primary_blocks <- tibble(
  block_id = seq_along(block_starts),
  block_label = sprintf("%d-%d", block_starts, block_starts + BLOCK_SIZE - 1L)
)
events[, block_label := primary_blocks$block_label[block_id]]

# 把事件分到 100km 网格
grid_cache <- file.path(P$derived_v2, "china_grid_100km_v2.rds")
grid_sf <- if (file.exists(grid_cache)) readRDS(grid_cache) else
  load_or_build_v2_grid(P$china_boundary, 1e5, grid_cache)

# 如果事件已经有 grid_cell（来自 stage-2 缓存），不重做空间 join
if (!"grid_cell" %in% names(events) || any(is.na(events$grid_cell))) {
  message("  assigning events to grid (sf planar)")
  old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  pts <- sf::st_as_sf(events[, .(longitude, latitude)],
                      coords = c("longitude", "latitude"),
                      crs = 4326, remove = FALSE)
  joined <- sf::st_join(pts, grid_sf[, "grid_cell"],
                        join = sf::st_intersects, left = TRUE)
  events[, grid_cell := joined$grid_cell]
}
events <- events[!is.na(grid_cell)]

events[, event_visit_id := paste(grid_cell,
                                  sprintf("%04d-%02d-%02d", year, month, day),
                                  tolower(coalesce(username, "anon")), sep = "|")]
events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]

## --- 1. 按 (grid × period × source) 汇总 -----------------------------------

eff_grid_block <- events[, .(
  n_events    = uniqueN(event_visit_id),
  n_observers = uniqueN(tolower(coalesce(username, "anon"))),
  n_dates     = uniqueN(sprintf("%04d-%02d-%02d", year, month, day)),
  mean_duration_min = mean(duration_min, na.rm = TRUE)
), by = .(grid_cell, block_label, source_short)]

eff_grid_block_total <- events[, .(
  n_events    = uniqueN(event_visit_id),
  n_observers = uniqueN(tolower(coalesce(username, "anon"))),
  n_dates     = uniqueN(sprintf("%04d-%02d-%02d", year, month, day)),
  mean_duration_min = mean(duration_min, na.rm = TRUE)
), by = .(grid_cell, block_label)][, source_short := "Combined"]

eff_all <- rbindlist(list(eff_grid_block, eff_grid_block_total),
                      use.names = TRUE, fill = TRUE)
eff_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
eff_all[, block_label := factor(block_label, levels = primary_blocks$block_label)]
write_csv(as_tibble(eff_all), v2_file("results",
                            "table_effort_by_grid_period_source"))

## --- 2. 大图：events × 5期 × 3源（合并 + Birdwatch + eBird） -------------

china_layers <- load_china_layers(P$china_boundary, P$province_line)

events_log_max <- log10(quantile(eff_all$n_events, 0.99, na.rm = TRUE) + 1)

eff_sf <- grid_sf |> inner_join(as_tibble(eff_all), by = "grid_cell")

map_events <- ggplot(eff_sf) +
  geom_sf(aes(fill = log10(n_events + 1)),
          colour = scales::alpha("white", 0.06), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_v2_sequential(name = "log10(events+1)", palette = "lajolla",
                           direction = 1,
                           limits = c(0, events_log_max),
                           oob = scales::squish) +
  facet_grid(source_short ~ block_label) +
  v2_china_coord() +
  theme_v2_map(9.4) +
  labs(title = "Survey effort across 100 km grids — events per source × 5-year period",
       subtitle = "log10(events+1); rows = data source (China Birdwatch / eBird / merged)",
       caption = "After dedup. Tile colour clipped at 99% quantile of events count.")
save_dual(map_events, "fig_effort_events_by_source_period_v2",
          width = 14, height = 7.6)

map_observers <- ggplot(eff_sf) +
  geom_sf(aes(fill = log10(n_observers + 1)),
          colour = scales::alpha("white", 0.06), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_v2_sequential(name = "log10(observers+1)", palette = "lajolla",
                           direction = 1, oob = scales::squish) +
  facet_grid(source_short ~ block_label) +
  v2_china_coord() +
  theme_v2_map(9.4) +
  labs(title = "Observer effort across 100 km grids — unique observers per source × 5-year period",
       subtitle = "log10(unique observers+1)",
       caption = "Rows = source.")
save_dual(map_observers, "fig_effort_observers_by_source_period_v2",
          width = 14, height = 7.6)

map_dur <- ggplot(eff_sf |> filter(!is.na(mean_duration_min)),
                   aes(fill = pmin(mean_duration_min, 240))) +
  geom_sf(aes(fill = pmin(mean_duration_min, 240)),
          colour = scales::alpha("white", 0.06), linewidth = 0.04) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_v2_sequential(name = "Mean checklist duration (min, capped at 240)",
                           palette = "lajolla", direction = 1,
                           limits = c(0, 240), oob = scales::squish) +
  facet_grid(source_short ~ block_label) +
  v2_china_coord() +
  theme_v2_map(9.4) +
  labs(title = "Survey duration across 100 km grids",
       subtitle = "Mean checklist duration in minutes (capped at 240 for visibility)",
       caption = "Note: eBird-derived records often lack duration → mostly NA in eBird row.")
save_dual(map_dur, "fig_effort_duration_by_source_period_v2",
          width = 14, height = 7.6)

## --- 3. 双源差值图（中国观鸟 - eBird） ------------------------------------

diff_tbl <- eff_grid_block[, .(grid_cell, block_label, source_short, n_events)] |>
  as_tibble() |>
  pivot_wider(names_from = source_short, values_from = n_events,
              values_fill = 0L) |>
  mutate(diff_log = log10((China_Birdwatch + 1) / (eBird_GBIF + 1)))
diff_sf <- grid_sf |>
  inner_join(diff_tbl, by = "grid_cell") |>
  mutate(block_label = factor(block_label, levels = primary_blocks$block_label))
diff_lim <- quantile(abs(diff_sf$diff_log), 0.97, na.rm = TRUE)
diff_map <- ggplot(diff_sf) +
  geom_sf(aes(fill = diff_log),
          colour = scales::alpha("white", 0.08), linewidth = 0.05) +
  china_map_layers(china_layers, country_lwd = 0.22, province_lwd = 0.12) +
  scale_fill_v2_diverging(name = "log10(China-Birdwatch / eBird) per grid",
                          limits = c(-diff_lim, diff_lim),
                          oob = scales::squish) +
  facet_wrap(~ block_label, ncol = 5) +
  v2_china_coord() +
  theme_v2_map(10.5) +
  labs(title = "Source dominance per grid: China Birdwatch vs eBird/GBIF",
       subtitle = "Red = China-Birdwatch dominant; Blue = eBird/GBIF dominant; White = balanced.",
       caption = sprintf("Values clipped at +/- %.2f (97%% quantile of |log-ratio|)", diff_lim))
save_dual(diff_map, "fig_effort_source_dominance_v2",
          width = 16, height = 4.8)

## --- 4. 年度时序：events × source ----------------------------------------

annual <- events[, .(
  n_events    = uniqueN(event_visit_id),
  n_observers = uniqueN(tolower(coalesce(username, "anon"))),
  n_grids     = uniqueN(grid_cell)
), by = .(year, source_short)]
annual_total <- events[, .(
  n_events    = uniqueN(event_visit_id),
  n_observers = uniqueN(tolower(coalesce(username, "anon"))),
  n_grids     = uniqueN(grid_cell)
), by = year][, source_short := "Combined"]
annual <- rbindlist(list(annual, annual_total), use.names = TRUE)
annual[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]

annual_long <- as_tibble(annual) |>
  pivot_longer(c(n_events, n_observers, n_grids), names_to = "metric",
               values_to = "value") |>
  mutate(metric = recode(metric,
    n_events    = "Visit events",
    n_observers = "Unique observers",
    n_grids     = "Grids visited"))

ts_plot <- ggplot(annual_long, aes(year, value, colour = source_short,
                                     fill = source_short)) +
  geom_area(data = annual_long |> filter(source_short != "Combined"),
            aes(group = source_short), alpha = 0.35,
            position = "stack", colour = NA) +
  geom_line(data = annual_long |> filter(source_short == "Combined"),
            aes(group = 1), colour = "#1F1F1F", linewidth = 0.8) +
  geom_point(data = annual_long |> filter(source_short == "Combined"),
             aes(group = 1), colour = "#1F1F1F", size = 1.4) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = seq(2000, 2024, 4)) +
  scale_fill_manual(values = c(China_Birdwatch = "#0E5A78",
                                eBird_GBIF = "#C9784A",
                                Combined = "#1F1F1F"),
                     name = "Source (stacked area = per source; black line = merged total)") +
  scale_colour_manual(values = c(China_Birdwatch = "#0E5A78",
                                  eBird_GBIF = "#C9784A",
                                  Combined = "#1F1F1F"), guide = "none") +
  labs(title = "Survey effort time series — by source",
       subtitle = "Stacked area = per-source contribution; black line = merged total after dedup",
       x = NULL, y = NULL) +
  theme_v2_pub(11) +
  theme(legend.position = "top")
save_dual(ts_plot, "fig_effort_timeseries_by_source_v2",
          width = 11, height = 8)

## --- 5. 物种贡献分布（每源覆盖的物种数） ----------------------------------

sp_per_source <- events[, .(n_species = uniqueN(species)),
                         by = .(source_short, year)] |>
  as_tibble()
sp_plot <- ggplot(sp_per_source,
                   aes(year, n_species, colour = source_short,
                        fill = source_short)) +
  geom_area(aes(group = source_short), alpha = 0.30,
            position = "identity", colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c(China_Birdwatch = "#0E5A78",
                                  eBird_GBIF = "#C9784A"), name = "Source") +
  scale_fill_manual(values = c(China_Birdwatch = "#0E5A78",
                                eBird_GBIF = "#C9784A"), guide = "none") +
  scale_x_continuous(breaks = seq(2000, 2024, 4)) +
  labs(title = "Species detected per year by data source",
       subtitle = "Note: this counts unique species observed; not corrected for occupancy.",
       x = NULL, y = "Unique species detected") +
  theme_v2_pub(11)
save_dual(sp_plot, "fig_effort_species_per_source_v2",
          width = 10, height = 5.4)

message("[stage-2b] Effort + source-comparison figures done.")
