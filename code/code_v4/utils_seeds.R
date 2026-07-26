#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   v4 统一种子管理；解决 v2/v3 中 20260501 / 20260502 / 随机调用混杂
#
# Objective / 分析目标:
#   set_seeds(stage)：按 stage 名派生固定种子，让 missForest / ranger /
#   brms / spOccupancy / sample 等所有随机操作可复现。
#
# Input data / 输入数据:
#   stage 字符串（如 "01_merge", "04_stMsPGOcc", "05_postprocess"）
#
# Main workflow / 主要流程:
#   1. hash stage 名到稳定整数
#   2. 同时设置 base R / future / brms / parallel 种子
#
# Key assumptions / 关键假设:
#   - digest 包可用（如不可用 fallback 到 nchar-based hash）
#
# Main packages / 主要包:
#   digest（可选）
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

# 基础种子（同一 stage 在不同会话稳定）
.BASE_SEED <- 20260511L

set_seeds <- function(stage, base_seed = .BASE_SEED) {
  # 1. stage hash → 整数偏移
  offset <- if (requireNamespace("digest", quietly = TRUE)) {
    # 取 hash 前 8 个字符转为整数
    h <- digest::digest(stage, algo = "xxhash32")
    strtoi(substr(h, 1, 6), base = 16L) %% 100000L
  } else {
    # fallback：用 stage 字符的 ASCII 和
    sum(utf8ToInt(stage)) %% 100000L
  }

  seed <- base_seed + offset
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")

  # 2. 同步设置 RNGkind 让 parallel / future 也用同一 stream
  if (requireNamespace("parallel", quietly = TRUE)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  message(sprintf("[set_seeds] stage='%s' seed=%d", stage, seed))
  invisible(seed)
}

message("[utils_seeds_v4] loaded")
