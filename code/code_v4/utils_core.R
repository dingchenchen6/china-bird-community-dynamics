#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   v4 基础 I/O 与日志工具；解决 v3 多处重复的 safe_read / log_time / candidate species fallback
#
# Objective / 分析目标:
#   - safe_read / read_csv_safe：兼容 NULL、缺失、压缩、qs/rds
#   - qs_save_safe：优先 qs::qsave 大对象；fallback saveRDS
#   - log_time：统一日志格式（stage + 时间戳）
#   - limit_candidate_species：单点定义候选物种限定逻辑
#   - checkpoint_save：brms/loo 等大对象统一落盘
#
# Input data / 输入数据:
#   00_config.R 与 utils_paths.R 已 source
#
# Main workflow / 主要流程:
#   纯函数定义
#
# Key assumptions / 关键假设:
#   readr / qs 包按需加载（qs 缺失时 fallback rds）
#
# Main packages / 主要包:
#   readr, qs（可选）
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

# ── safe_read：rds / qs 双格式（自动判别后缀） ───────────────────────
# tries .rds first, then .qs; returns NULL if both fail
safe_read <- function(path, quiet = FALSE) {
  if (is.null(path) || !is.character(path) || length(path) == 0) return(NULL)

  # 直接路径存在
  if (file.exists(path)) {
    ext <- tolower(tools::file_ext(path))
    return(tryCatch({
      if (ext == "qs" && requireNamespace("qs", quietly = TRUE)) {
        qs::qread(path)
      } else if (ext == "rds" || ext == "") {
        readRDS(path)
      } else if (ext == "csv") {
        if (requireNamespace("readr", quietly = TRUE)) {
          readr::read_csv(path, show_col_types = FALSE)
        } else {
          read.csv(path, stringsAsFactors = FALSE)
        }
      } else {
        readRDS(path)
      }
    }, error = function(e) {
      if (!quiet) message("[safe_read] FAIL: ", basename(path), " — ", e$message)
      NULL
    }))
  }

  # 路径未带后缀：尝试 .rds 或 .qs
  for (ext in c("rds", "qs")) {
    p2 <- paste0(path, ".", ext)
    if (file.exists(p2)) return(safe_read(p2, quiet = quiet))
  }

  if (!quiet) message("[safe_read] NOT FOUND: ", path)
  NULL
}

# ── read_csv_safe：自动加 .csv 后缀 ──────────────────────────────────
read_csv_safe <- function(path_or_stem, quiet = FALSE) {
  path <- if (grepl("\\.(csv|tsv|txt)$", path_or_stem, ignore.case = TRUE)) {
    path_or_stem
  } else {
    paste0(path_or_stem, ".csv")
  }
  if (!file.exists(path)) {
    if (!quiet) message("[read_csv_safe] NOT FOUND: ", path)
    return(NULL)
  }
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) {
      if (!quiet) message("[read_csv_safe] FAIL: ", path, " — ", e$message)
      NULL
    }
  )
}

# ── qs_save_safe：大对象优先 qs，fallback rds ────────────────────────
qs_save_safe <- function(obj, path, preset = "fast") {
  ext <- tolower(tools::file_ext(path))
  if (ext == "" || ext == "rds") path <- paste0(tools::file_path_sans_ext(path), ".qs")

  if (requireNamespace("qs", quietly = TRUE)) {
    tryCatch({
      qs::qsave(obj, path, preset = preset)
      sz_mb <- file.info(path)$size / 1e6
      message(sprintf("[qs_save] %s (%.1f MB, qs preset=%s)",
                      basename(path), sz_mb, preset))
      return(invisible(path))
    }, error = function(e) {
      message("[qs_save] qsave failed (", e$message, ") — fallback rds")
    })
  }

  rds_path <- paste0(tools::file_path_sans_ext(path), ".rds")
  saveRDS(obj, rds_path, compress = "xz")
  sz_mb <- file.info(rds_path)$size / 1e6
  message(sprintf("[qs_save] %s (%.1f MB, rds xz fallback)",
                  basename(rds_path), sz_mb))
  invisible(rds_path)
}

# ── log_time：统一时间戳日志 ─────────────────────────────────────────
log_time <- function(stage, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s | %s", ts, stage, msg))
}

# ── limit_candidate_species：单点定义候选物种限定 ────────────────────
# v3 中此逻辑在 04 / 15 重复实现，v4 统一
limit_candidate_species <- function(survey, max_n = 200L,
                                     v3_results_dir = DIRS$v3_results,
                                     v2_results_dir = DIRS$v2_results) {
  # 优先使用 v3 候选清单，否则 v2，否则 survey 内默认顺序
  paths <- c(
    file.path(v3_results_dir, "table_dynamic_occupancy_candidate_species_all.csv"),
    file.path(v2_results_dir, "table_dynamic_occupancy_candidate_species_all.csv")
  )
  for (p in paths) {
    if (file.exists(p)) {
      df <- read_csv_safe(p)
      if (!is.null(df) && "species" %in% names(df)) {
        out <- head(df$species, max_n)
        message(sprintf("[limit_species] %d species from %s",
                        length(out), basename(p)))
        return(out)
      }
    }
  }
  out <- head(survey$species, max_n)
  message(sprintf("[limit_species] %d species from survey fallback", length(out)))
  out
}

# ── checkpoint_save：brms/loo 等中间对象通用落盘 ──────────────────────
checkpoint_save <- function(obj, stem, subdir = "derived", ext = "qs") {
  dir <- DIRS[[subdir]] %||% DIRS$derived
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(stem, ".", ext))
  qs_save_safe(obj, path)
}

message("[utils_core_v4] loaded")
