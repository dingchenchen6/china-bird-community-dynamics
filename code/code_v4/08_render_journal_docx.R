#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把 results_v4/manuscript_v4_zh.md 渲染成 .docx（投稿格式）
#
# Objective / 分析目标:
#   - pandoc 一次直接生成 docx（如有 reference docx 模板则使用）
#   - 同时生成英文 abstract 单文件 docx
#
# Input data / 输入数据:
#   results_v4/manuscript_v4_zh.md
#   results_v4/manuscript_v4_en_abstract.md（如有）
#
# Main workflow / 主要流程:
#   call system pandoc
#
# Key assumptions / 关键假设:
#   pandoc 系统安装
#
# Main packages / 主要包:
#   无（system call）
#
# Output directory / 输出路径:
#   results_v4/manuscript_v4_zh.docx
#   results_v4/manuscript_v4_en_abstract.docx
# ============================================================

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
P <- ensure_v4_dirs()

log_time("08", "Rendering manuscript to docx")

if (Sys.which("pandoc") == "") {
  warning("[08] pandoc not found in PATH. Skipping docx render.")
  quit(status = 0)
}

renders <- list(
  list(md = "manuscript_v4_zh.md",            docx = "manuscript_v4_zh.docx"),
  list(md = "manuscript_v4_en_abstract.md",   docx = "manuscript_v4_en_abstract.docx"),
  list(md = "cover_letter_v4.md",             docx = "cover_letter_v4.docx"),
  list(md = "supplementary_outline_v4.md",    docx = "supplementary_outline_v4.docx")
)

for (r in renders) {
  in_md <- file.path(DIRS$results, r$md)
  out_docx <- file.path(DIRS$results, r$docx)
  if (!file.exists(in_md)) {
    message(sprintf("[08] Skip (input missing): %s", r$md))
    next
  }
  cmd <- sprintf("pandoc %s -o %s --from markdown --to docx",
                 shQuote(in_md), shQuote(out_docx))
  rv <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (rv == 0) {
    message(sprintf("[08] Rendered: %s → %s", r$md, r$docx))
  } else {
    warning(sprintf("[08] pandoc failed for %s", r$md))
  }
}

log_time("08", "DONE")
