#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   v4 单点定义项目根 + 文件路径函数，消除 v2/v3 中分散的 fallback
#
# Objective / 分析目标:
#   提供 find_project_root() / v4_file() / v3_file() / v2_file() /
#   project_file() / log_path() 等纯函数；其他 utils 一律 source 之
#
# Input data / 输入数据:
#   00_config.R 已 source（提供 DIRS / PROJECT_ROOT）
#
# Main workflow / 主要流程:
#   定义路径拼装函数；提供 %||% 操作符（v3 多处依赖但定义分散）
#
# Key assumptions / 关键假设:
#   00_config.R 已先 source；DIRS 列表存在
#
# Main packages / 主要包:
#   base R
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

# ── %||% 操作符（NULL 合并） ─────────────────────────────────────────
# 已有定义则不覆盖（兼容 rlang 等包）
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

# ── find_project_root（单点定义，其他 utils 一律调用此函数） ──────────
find_project_root <- function() {
  envroot <- Sys.getenv("BIRD_PROJECT_ROOT", "")
  if (nzchar(envroot) && dir.exists(envroot)) return(normalizePath(envroot))

  candidate <- file.path("~", "Documents", "New project",
                        "bird_dynamic_occupancy_analysis")
  if (dir.exists(candidate)) return(normalizePath(candidate))

  # fallback：检测当前 wd 向上回溯
  wd <- getwd()
  while (wd != "/" && wd != "") {
    if (file.exists(file.path(wd, "code_v3")) ||
        file.exists(file.path(wd, "code_v4"))) {
      return(normalizePath(wd))
    }
    wd <- dirname(wd)
  }
  normalizePath(getwd())
}

# ── 路径拼装：v4 / v3 / v2 / project / external ──────────────────────
v4_file <- function(category, stem, ext = "csv") {
  dir <- DIRS[[category]] %||% file.path(PROJECT_ROOT, paste0(category, "_v4"))
  file.path(dir, paste0(stem, ".", ext))
}

v3_file <- function(category, stem, ext = "csv") {
  dir <- switch(category,
                derived = DIRS$v3_derived,
                results = DIRS$v3_results,
                figures = file.path(PROJECT_ROOT, "figures_v3"),
                logs    = file.path(PROJECT_ROOT, "logs_v3"),
                DIRS[[category]])
  file.path(dir, paste0(stem, ".", ext))
}

v2_file <- function(category, stem, ext = "csv") {
  dir <- switch(category,
                derived = DIRS$v2_derived,
                results = DIRS$v2_results,
                figures = file.path(PROJECT_ROOT, "figures_v2"),
                logs    = file.path(PROJECT_ROOT, "logs_v2"),
                DIRS[[category]])
  file.path(dir, paste0(stem, ".", ext))
}

project_file <- function(...) file.path(PROJECT_ROOT, ...)

external_file <- function(...) file.path(DIRS$external, ...)

# ── 日志路径 ─────────────────────────────────────────────────────
log_path <- function(stem, stage = "") {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  fn <- if (nzchar(stage)) paste0(stage, "_", stem, "_", ts, ".log") else
                            paste0(stem, "_", ts, ".log")
  file.path(DIRS$logs, fn)
}

message("[utils_paths_v4] loaded")
