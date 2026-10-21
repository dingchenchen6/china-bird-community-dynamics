#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题: 重拟合所用的输入是否完整、结构正确、可追溯？
#   Are the refit inputs complete, correctly shaped and traceable?
# Objective / 分析目标: 在隔离的输出根目录（BIRD_OUTPUT_ROOT）中准备输入，校验调查历史为四维并含
#   检测协变量，校验环境与变化表缺失率，冻结 200/500 种名单，生成 SHA-256 清单。
# Input / 输入: <PROJECT>/data/derived_v3/*（只读）；S0_01 已写入重建的 grid_environment_v3.rds
# Output / 输出: <OUT>/data/derived_v3/（输入副本）、<OUT>/results_v3/species_list_{200,500}.csv、
#   <OUT>/results_v3/INPUT_SHA256.csv、<OUT>/results_v3/QA_inputs.csv
# Key assumptions / 关键假设: 大型只读输入（去重事件表）以符号链接引用，其余输入复制，
#   避免下游脚本写回时覆盖 v3 原件。
# Packages / 主要包: digest
# ============================================================
suppressPackageStartupMessages(library(digest))
PROJ <- Sys.getenv("BIRD_PROJECT_ROOT", path.expand("~/Documents/New project/bird_dynamic_occupancy_analysis"))
OUT  <- Sys.getenv("BIRD_OUTPUT_ROOT", file.path(PROJ, "refit_v4"))
src  <- file.path(PROJ, "data/derived_v3"); dst <- file.path(OUT, "data/derived_v3")
dir.create(dst, recursive = TRUE, showWarnings = FALSE); dir.create(file.path(OUT, "results_v3"), showWarnings = FALSE)

# 1. 准备输入（不带入任何旧模型输出）/ stage inputs, excluding old model outputs
skip <- "^(stMsPGOcc_fit_|tMsPGOcc_fit_|psi_samples|metric_arrays|diversity_v3_per_draw|grid_environment)"
for (f in list.files(src, full.names = TRUE)) {
  b <- basename(f); to <- file.path(dst, b)
  if (grepl(skip, b) || file.exists(to)) next
  if (grepl("^combined_events_dedup", b)) file.symlink(f, to) else file.copy(f, to)   # 大型只读输入用链接 / link big read-only input
}

if (Sys.getenv("REFIT_STAGE_ONLY", "0") == "1") { message("[S0_02] 已准备输入（仅准备模式）"); quit(save = "no") }

# 2. 调查历史结构 / survey-history structure
s <- readRDS(file.path(dst, "survey_history_v3.rds"))
need_det <- c("log_events", "log_duration", "has_duration", "source_ebird_prop")
ok4d <- length(dim(s$y)) == 4 && !is.null(s$det_covs) && all(need_det %in% names(s$det_covs)) && !is.null(s$visited_mask)
if (!ok4d) stop("[S0_02] 调查历史不是四维或缺少检测协变量（当前维度：", paste(dim(s$y), collapse = " × "),
                "）。请先运行 02_build_survey_history.R 重建，再重跑本脚本。")
reps <- apply(s$visited_mask, c(1, 2), sum)
message(sprintf("[S0_02] y = %s；每格每期访问次数 1 次占比 %.1f%%（rep = 1 退化检查）",
                paste(dim(s$y), collapse = " × "), 100 * mean(reps[reps > 0] == 1)))

# 3. 环境与变化表缺失率 / missingness of environment and change tables
qa <- do.call(rbind, lapply(c("grid_environment_v3.rds", "climate_change_v3.rds", "hfi_change_v3.rds", "landuse_change_v3.rds"), function(b) {
  p <- file.path(dst, b); if (!file.exists(p)) return(data.frame(table = b, variable = NA, missing_share = NA, note = "缺文件"))
  x <- as.data.frame(readRDS(p)); x <- x[match(s$sites, x$grid_cell), , drop = FALSE]
  num <- names(x)[vapply(x, is.numeric, TRUE)]
  data.frame(table = b, variable = num, missing_share = round(vapply(num, function(v) mean(!is.finite(x[[v]])), 0), 4), note = "")
}))
write.csv(qa, file.path(OUT, "results_v3/QA_inputs.csv"), row.names = FALSE)
bad <- qa[!is.na(qa$missing_share) & qa$missing_share > 0.02, ]
if (nrow(bad)) { print(bad, row.names = FALSE); warning("[S0_02] 上述变量缺失率 > 2%，用于模型或驱动分析前须先修复") }

# 4. 冻结物种名单（与 04 的取法一致：按记录数排序取前 n）/ freeze species lists as 04 selects them
for (n in c(200L, 500L)) write.csv(data.frame(rank = seq_len(n), species = head(s$species, n)),
                                   file.path(OUT, sprintf("results_v3/species_list_%d.csv", n)), row.names = FALSE)

# 5. SHA-256 清单 / checksum manifest
fs <- list.files(dst, full.names = TRUE)
man <- data.frame(file = basename(fs), bytes = file.size(fs),
                  sha256 = vapply(fs, function(f) digest(file = f, algo = "sha256"), ""))
write.csv(man, file.path(OUT, "results_v3/INPUT_SHA256.csv"), row.names = FALSE)
message(sprintf("[S0_02] 完成：%d 个输入文件已登记校验和", nrow(man)))
