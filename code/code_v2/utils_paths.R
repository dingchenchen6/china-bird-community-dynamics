## utils_paths.R
## v2 共享路径解析工具：所有 v2 脚本统一从此 source。
## 不允许在业务脚本里再写绝对路径。

suppressPackageStartupMessages({
  library(fs)
})

## 全局工具：null-coalesce（必须在所有用到 %||% 之前 source 此文件）。
`%||%` <- function(a, b) if (is.null(a)) b else a

## v2_script_dir(): 不依赖 sys.frame —— 用 commandArgs 或环境变量回退。
v2_script_dir <- function() {
  # 1) 显式环境变量（CI / batch 推荐）
  env_dir <- Sys.getenv("V2_CODE_DIR", unset = "")
  if (nzchar(env_dir) && dir.exists(env_dir)) return(normalizePath(env_dir))
  # 2) Rscript --file=
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) {
    f <- sub("^--file=", "", hit[1])
    if (file.exists(f)) return(normalizePath(dirname(f)))
  }
  # 3) sys.frame ofile（source() 时可用）
  if (sys.nframe() > 0) {
    for (k in seq_len(sys.nframe())) {
      f <- sys.frame(k)$ofile
      if (!is.null(f) && file.exists(f)) return(normalizePath(dirname(f)))
    }
  }
  # 4) 最后回退到约定路径
  fallback <- "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2"
  if (dir.exists(fallback)) return(fallback)
  stop("Cannot resolve v2 code directory; set V2_CODE_DIR.")
}

resolve_task_root <- function() {
  candidates <- c(
    Sys.getenv("BIRD_DOA_ROOT", ""),
    "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
  )
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (length(candidates) == 0) stop("Cannot resolve task root; set BIRD_DOA_ROOT.")
  normalizePath(candidates[1])
}

v2_paths <- function() {
  task_root <- resolve_task_root()
  project_root <- normalizePath(dirname(task_root))
  list(
    task_root        = task_root,
    project_root     = project_root,
    code_v2          = file.path(task_root, "code_v2"),
    derived_v2       = file.path(task_root, "data", "derived_v2"),
    results_v2       = file.path(task_root, "results_v2"),
    figures_v2       = file.path(task_root, "figures_v2"),
    logs_v2          = file.path(task_root, "logs_v2"),
    derived_v3       = file.path(task_root, "data", "derived_v3"),
    results_v3       = file.path(task_root, "results_v3"),
    figures_v3       = file.path(task_root, "figures_v3"),
    derived_v1       = file.path(task_root, "data", "derived"),
    results_v1       = file.path(task_root, "results"),
    external_grid    = file.path(project_root, "bird_grid_community_analysis", "data", "external"),
    external_climate = file.path(project_root, "bird_full_community_analysis", "data", "external", "climate", "wc2.1_10m"),
    # 优先用 data/中国shp/（最新标准，CGCS2000，含十段线，无鹰眼图 inset）。
    # 回退顺序：data/中国shp/ -> data/china_boundary/ -> bird_grid_community_analysis/...
    china_boundary   = {
      candidates <- c(
        file.path(task_root, "data", "中国shp", "省.shp"),
        file.path(task_root, "data", "china_boundary", "省面.shp"),
        file.path(project_root, "bird_grid_community_analysis",
                  "data", "external", "china_boundary", "省面.shp"))
      Find(file.exists, candidates)
    },
    province_line    = {
      candidates <- c(
        file.path(task_root, "data", "中国shp", "省_境界线.shp"),
        file.path(task_root, "data", "china_boundary", "省界线.shp"),
        file.path(project_root, "bird_grid_community_analysis",
                  "data", "external", "china_boundary", "省界线.shp"))
      Find(file.exists, candidates)
    },
    ten_dash_line    = {
      candidates <- c(
        file.path(task_root, "data", "中国shp", "十段线.shp"),
        file.path(task_root, "data", "china_boundary", "十段线.shp"),
        file.path(project_root, "bird_grid_community_analysis",
                  "data", "external", "china_boundary", "十段线.shp"))
      Find(file.exists, candidates)
    }
  )
}

ensure_v2_dirs <- function() {
  P <- v2_paths()
  for (d in c(P$derived_v2, P$results_v2, P$figures_v2, P$logs_v2,
              P$derived_v3, P$results_v3, P$figures_v3)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(P)
}

## v3_file: 与 v2_file 同接口，但导向 derived_v3 / results_v3 / figures_v3
v3_file <- function(kind = c("results", "figure", "derived"), stem, ext = "csv") {
  kind <- match.arg(kind)
  P <- v2_paths()
  base <- switch(kind,
    results  = P$results_v3,
    figure   = P$figures_v3,
    derived  = P$derived_v3
  )
  file.path(base, paste0(stem, ".", ext))
}

v2_file <- function(kind = c("results", "figure", "derived", "log"), stem, ext = "csv") {
  kind <- match.arg(kind)
  P <- v2_paths()
  base <- switch(kind,
    results  = P$results_v2,
    figure   = P$figures_v2,
    derived  = P$derived_v2,
    log      = P$logs_v2
  )
  file.path(base, paste0(stem, ".", ext))
}
