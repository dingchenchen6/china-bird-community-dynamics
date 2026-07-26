#!/usr/bin/env Rscript
## 02f_effort_by_province_species.R
##
## 按 年份 × 省份 × 数据源 汇总调查努力指标（网格级 + 省份级）
## 按 物种 × 年份 × 省份 × 数据源 汇总物种水平调查努力指标
##
## 产出：
##   results_v2/table_effort_year_province_source.csv         — 年份 × 省份 × 源
##   results_v2/table_effort_species_year_province_source.csv — 物种 × 年份 × 省份 × 源

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
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

## =============================================================================
## 0. 加载去重事件
## =============================================================================

message("[02f] Loading dedup'd events")
events <- as.data.table(readRDS(file.path(P$derived_v2,
                            "combined_events_merged_dedup_2000_2025.rds")))
events <- events[year >= ANALYSIS_YR_LO & year <= ANALYSIS_YR_HI]
message(sprintf("  events 2000-2024: %s", format(nrow(events), big.mark = ",")))

events[, source_short := ifelse(grepl("Birdwatch", as.character(source)),
                                 "China_Birdwatch", "eBird_GBIF")]
events[, day_key := sprintf("%04d-%02d-%02d", year, month, day)]
events[, observer_key := tolower(coalesce(username, "anon"))]
events[, visit_key := paste(observer_key, day_key, sep = "|")]
events[, has_dur := !is.na(duration_min)]

## =============================================================================
## 1. 省份空间 join：把事件分配到省份
## =============================================================================

message("[02f] Loading province boundaries and assigning events to provinces")

shp_dir <- file.path(P$task_root, "data", "中国shp")
china_prov <- suppressWarnings(sf::st_read(file.path(shp_dir, "省.shp"), quiet = TRUE)) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

# 去除鹰眼图子多边形
china_prov <- strip_nansha_inset(china_prov)

message(sprintf("  provinces: %d", nrow(china_prov)))

# 空间 join（使用 planar 避免 s2 慢）
old_s2 <- sf::sf_use_s2(FALSE)
on.exit(sf::sf_use_s2(old_s2), add = TRUE)

pts <- sf::st_as_sf(events[, .(longitude, latitude)],
                    coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
joined <- sf::st_join(pts, china_prov[, c("省名", "省代码")],
                      join = sf::st_intersects, left = TRUE)
events[, province := joined$省名]
events[, province_code := joined$省代码]

sf::sf_use_s2(old_s2)

events_prov <- events[!is.na(province)]
message(sprintf("  events with province: %s (%.1f%%)",
                format(nrow(events_prov), big.mark = ","),
                100 * nrow(events_prov) / nrow(events)))

## =============================================================================
## 2. 年份 × 省份 × 数据源 调查努力指标
## =============================================================================

message("[02f] Computing year x province x source effort metrics")

# 分源
yr_prov_src <- events_prov[, .(
  n_records         = .N,
  n_visits          = uniqueN(visit_key),
  n_observers       = uniqueN(observer_key),
  n_birding_days    = uniqueN(day_key),
  n_species         = uniqueN(species),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE),
  total_duration_min = sum(duration_min[has_dur], na.rm = TRUE)
), by = .(year, province, province_code, source_short)]

# 合并
yr_prov_comb <- events_prov[, .(
  n_records         = .N,
  n_visits          = uniqueN(visit_key),
  n_observers       = uniqueN(observer_key),
  n_birding_days    = uniqueN(day_key),
  n_species         = uniqueN(species),
  mean_duration_min = mean(duration_min[has_dur], na.rm = TRUE),
  total_duration_min = sum(duration_min[has_dur], na.rm = TRUE)
), by = .(year, province, province_code)][, source_short := "Combined"]

yr_prov_all <- rbindlist(list(yr_prov_src, yr_prov_comb), use.names = TRUE)
yr_prov_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
setorder(yr_prov_all, province, source_short, year)

write_csv(as_tibble(yr_prov_all),
          file.path(P$results_v2, "table_effort_year_province_source.csv"))

message(sprintf("  Year x Province table: %d rows, %d provinces",
                nrow(yr_prov_all), uniqueN(yr_prov_all$province)))

## =============================================================================
## 3. 物种 × 年份 × 省份 × 数据源 调查努力指标
## =============================================================================

message("[02f] Computing species x year x province x source effort metrics")

# 物种-省份-年份-源 汇总
# 注意：这里 n_visits / n_observers 等是"该物种被检测到的访问/观测者"数
# 即只计算检测到该物种的那些 visit，而非省份全部 visit
sp_yr_prov_src <- events_prov[, .(
  n_detection_records = .N,
  n_detection_visits  = uniqueN(visit_key),
  n_detection_observers = uniqueN(observer_key),
  n_detection_days    = uniqueN(day_key)
), by = .(species, year, province, province_code, source_short)]

sp_yr_prov_comb <- events_prov[, .(
  n_detection_records = .N,
  n_detection_visits  = uniqueN(visit_key),
  n_detection_observers = uniqueN(observer_key),
  n_detection_days    = uniqueN(day_key)
), by = .(species, year, province, province_code)][, source_short := "Combined"]

sp_yr_prov_all <- rbindlist(list(sp_yr_prov_src, sp_yr_prov_comb), use.names = TRUE)
sp_yr_prov_all[, source_short := factor(source_short,
  levels = c("China_Birdwatch", "eBird_GBIF", "Combined"))]
setorder(sp_yr_prov_all, species, province, source_short, year)

write_csv(as_tibble(sp_yr_prov_all),
          file.path(P$results_v2, "table_effort_species_year_province_source.csv"))

message(sprintf("  Species x Year x Province table: %d rows, %d species, %d provinces",
                nrow(sp_yr_prov_all),
                uniqueN(sp_yr_prov_all$species),
                uniqueN(sp_yr_prov_all$province)))

## =============================================================================
## 4. 省份级年度动态地图（补充）
## =============================================================================

message("[02f] Rendering province-level annual maps")

china_layers <- load_china_layers(
  file.path(shp_dir, "省.shp"),
  file.path(shp_dir, "省_境界线.shp"),
  file.path(shp_dir, "十段线.shp"),
  strip_inset = TRUE
)

anim_dir_prov <- file.path(P$figures_v2, "anim_province")
dir.create(anim_dir_prov, recursive = TRUE, showWarnings = FALSE)

# 合并版各指标的最大值（用于统一色标）
yr_prov_comb_only <- yr_prov_all[source_short == "Combined"]
prov_visits_max   <- quantile(yr_prov_comb_only$n_visits, 0.99, na.rm = TRUE)
prov_obs_max      <- quantile(yr_prov_comb_only$n_observers, 0.99, na.rm = TRUE)
prov_species_max  <- quantile(yr_prov_comb_only$n_species, 0.99, na.rm = TRUE)
prov_days_max     <- quantile(yr_prov_comb_only$n_birding_days, 0.99, na.rm = TRUE)

prov_metrics <- list(
  list(col = "n_visits",       label = "Visit events",    max = prov_visits_max,  stem = "visits"),
  list(col = "n_observers",    label = "Unique observers", max = prov_obs_max,     stem = "observers"),
  list(col = "n_species",      label = "Species detected", max = prov_species_max, stem = "species"),
  list(col = "n_birding_days", label = "Birding days",    max = prov_days_max,    stem = "days")
)

library(gifski)

for (m in prov_metrics) {
  message(sprintf("  Province Combined: %s", m$stem))

  frames <- lapply(ANALYSIS_YR_LO:ANALYSIS_YR_HI, function(yr) {
    sub <- yr_prov_comb_only[year == yr]
    if (nrow(sub) == 0) return(NULL)

    sf_yr <- china_prov |> left_join(as_tibble(sub), by = c("省名" = "province"))
    fv <- rlang::sym(m$col)

    p <- ggplot(sf_yr) +
      geom_sf(aes(fill = !!fv), colour = "#AAAAAA", linewidth = 0.15) +
      china_map_layers(china_layers,
                       country_colour = "#1F1F1F", country_lwd = 0.5,
                       province_colour = NA, province_lwd = 0,
                       with_ten_dash = TRUE, ten_dash_colour = "#1F1F1F",
                       ten_dash_lwd = 0.45) +
      scale_fill_v2_sequential(name = m$label, palette = "lajolla", direction = 1,
                               limits = c(0, m$max), oob = scales::squish,
                               na.value = "white") +
      v2_china_coord() + theme_v2_map(11) +
      labs(title = sprintf("Survey effort by province — %d", yr),
           subtitle = sprintf("%s (Combined sources)", m$label))

    out <- file.path(anim_dir_prov, sprintf("frame_%s_%04d.png", m$stem, yr))
    ggsave(out, p, width = 9.4, height = 6.6, dpi = 200, bg = "white")
    out
  })
  frames <- unlist(frames)

  gif_path <- file.path(P$figures_v2,
                        sprintf("fig_effort_annual_province_%s.gif", m$stem))
  gifski::gifski(frames, gif_path, width = 1880, height = 1320,
                 delay = 0.8, loop = TRUE)
}

message("[02f] === All done ===")
message("  Tables:")
message("    table_effort_year_province_source.csv")
message("    table_effort_species_year_province_source.csv")
message("  Province GIFs:")
for (m in prov_metrics)
  message(sprintf("    fig_effort_annual_province_%s.gif", m$stem))
