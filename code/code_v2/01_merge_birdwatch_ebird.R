#!/usr/bin/env Rscript
## 01_merge_birdwatch_ebird.R
##
## 阶段 1：合并中国观鸟记录平台 + GBIF/eBird China 事件，并按统一 key 去重。
## 输入（v1 已缓存的 RDS）：
##   data/derived/birdwatch_events_1980_2025.rds   (11.7M rows, 19 cols)
##   data/derived/gbif_ebird_events_2000_2025.rds  ( 2.9M rows, 19 cols)
## 输出：
##   data/derived_v2/combined_events_merged_dedup_2000_2025.rds
##   results_v2/table_dedup_audit_by_source_year.csv
##   results_v2/table_dedup_audit_summary.csv
##
## 去重策略：
##   key = species_canonical × event_date × round(lon,4) × round(lat,4) × username
##   优先保留 source=China_Birdwatch_Platform（元数据更全：duration_min, taxon_count_event 等）；
##   同源时保留信息完整度更高的一行（按 non-NA 字段计数 + duration_min 非缺失优先）。

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))

P <- ensure_v2_dirs()

YEAR_LO <- as.integer(Sys.getenv("V2_YEAR_LO", "2000"))
YEAR_HI <- as.integer(Sys.getenv("V2_YEAR_HI", "2025"))
COORD_DIGITS <- as.integer(Sys.getenv("V2_COORD_DIGITS", "4"))   # 4 ≈ 10 m

message(sprintf("[stage-1] Loading raw event caches"))
bw_path <- file.path(P$derived_v1, "birdwatch_events_1980_2025.rds")
gb_path <- file.path(P$derived_v1, "gbif_ebird_events_2000_2025.rds")
stopifnot(file.exists(bw_path), file.exists(gb_path))

bw <- as.data.table(readRDS(bw_path))
gb <- as.data.table(readRDS(gb_path))

message(sprintf("  birdwatch raw: %s rows", format(nrow(bw), big.mark = ",")))
message(sprintf("  gbif/ebird raw: %s rows", format(nrow(gb), big.mark = ",")))

## --- 1. schema 对齐 + 年份过滤 ---------------------------------------------

shared_cols <- intersect(names(bw), names(gb))
bw <- bw[, ..shared_cols]
gb <- gb[, ..shared_cols]

clean_username <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x <- str_trim(x)
  tolower(x)
}

clean_species <- function(x) str_trim(as.character(x))

filter_year <- function(dt) {
  dt[!is.na(year) & year >= YEAR_LO & year <= YEAR_HI &
       !is.na(longitude) & !is.na(latitude) &
       longitude >= 70 & longitude <= 140 &
       latitude >= 15 & latitude <= 55]
}

bw <- filter_year(bw)
gb <- filter_year(gb)

bw[, source := factor(source, levels = c("China_Birdwatch_Platform", "GBIF_eBird_China"))]
gb[, source := factor(source, levels = c("China_Birdwatch_Platform", "GBIF_eBird_China"))]

events <- rbindlist(list(bw, gb), use.names = TRUE, fill = TRUE)
rm(bw, gb); gc(verbose = FALSE)

message(sprintf("  after year/coord filter: %s rows", format(nrow(events), big.mark = ",")))

## --- 2. 构建 dedup key -----------------------------------------------------

events[, `:=`(
  species_key  = clean_species(species),
  username_key = clean_username(username),
  lon_key      = round(longitude, COORD_DIGITS),
  lat_key      = round(latitude,  COORD_DIGITS),
  date_key     = sprintf("%04d-%02d-%02d", year, month, day)
)]

# 缺 username 的记录用经纬度 + 时间 + 物种作为备选 key（避免 NA username 全部归一类）
events[is.na(username_key), username_key := paste0("anon_",
  sprintf("%.4f_%.4f_%s", lon_key, lat_key, date_key))]

events[, dedup_key := paste(species_key, date_key, lon_key, lat_key, username_key,
                             sep = "|")]

n_before <- nrow(events)
n_unique_keys <- uniqueN(events$dedup_key)
message(sprintf("  total rows: %s; unique dedup keys: %s; expected drops: %s",
                format(n_before, big.mark = ","),
                format(n_unique_keys, big.mark = ","),
                format(n_before - n_unique_keys, big.mark = ",")))

## --- 3. 同 key 内排序、保留首条 -------------------------------------------

# 优先级：source（Birdwatch 优先）→ duration_min 非 NA 优先 → 非 NA 字段更多优先
events[, n_nonna := rowSums(!is.na(.SD)),
       .SDcols = setdiff(names(events), c("species_key","username_key","lon_key",
                                          "lat_key","date_key","dedup_key"))]
events[, has_duration := !is.na(duration_min)]

setorder(events, dedup_key, source, -has_duration, -n_nonna)
events_dedup <- events[!duplicated(dedup_key)]

n_after <- nrow(events_dedup)
n_dropped <- n_before - n_after
message(sprintf("  dedup done: kept=%s, dropped=%s (%.2f%%)",
                format(n_after, big.mark = ","),
                format(n_dropped, big.mark = ","),
                100 * n_dropped / n_before))

## --- 4. 审计：按 source × year 的留存 / 删除统计 ----------------------------

events[, kept_after_dedup := !duplicated(dedup_key)]

audit_by_source_year <- events[, .(
  n_total   = .N,
  n_kept    = sum(kept_after_dedup),
  n_dropped = sum(!kept_after_dedup)
), by = .(source, year)][order(source, year)]
audit_by_source_year[, pct_dropped := 100 * n_dropped / n_total]

audit_summary <- events[, .(
  n_total          = .N,
  n_kept           = sum(kept_after_dedup),
  n_dropped        = sum(!kept_after_dedup),
  pct_dropped      = 100 * sum(!kept_after_dedup) / .N
), by = source]
audit_summary <- rbind(
  audit_summary,
  data.table(source = "_TOTAL_",
             n_total = n_before,
             n_kept = n_after,
             n_dropped = n_dropped,
             pct_dropped = 100 * n_dropped / n_before)
)

write_csv(audit_by_source_year, v2_file("results", "table_dedup_audit_by_source_year"))
write_csv(audit_summary,        v2_file("results", "table_dedup_audit_summary"))

## --- 5. 清理临时列 + 落盘 ---------------------------------------------------

events_dedup[, c("species_key","username_key","lon_key","lat_key","date_key",
                 "dedup_key","n_nonna","has_duration","kept_after_dedup") := NULL]

out_path <- file.path(P$derived_v2, "combined_events_merged_dedup_2000_2025.rds")
saveRDS(events_dedup, out_path, compress = FALSE)
message(sprintf("  wrote %s (%.1f MB)",
                out_path,
                file.info(out_path)$size / 1024^2))

## --- 6. 简短打印审计摘要 ---------------------------------------------------

cat("\n=== Dedup audit summary ===\n")
print(audit_summary)
cat("\n=== Top 5 source-year drop rates ===\n")
print(head(audit_by_source_year[order(-pct_dropped)], 5))

message("[stage-1] Merge + dedup done.")
