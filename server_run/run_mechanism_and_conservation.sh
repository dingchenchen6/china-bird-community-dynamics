#!/usr/bin/env bash
# ============================================================
# run_mechanism_and_conservation.sh
#
# 一键运行"机制归因 + 保护成效/空缺/规划"分析链（脚本 30 → 31 → 32）。
# One-command driver for the mechanism + conservation analysis chain.
#
# 用法 / Usage:
#   bash run_mechanism_and_conservation.sh                # 全部跑
#   bash run_mechanism_and_conservation.sh 30             # 只跑脚本 30
#   bash run_mechanism_and_conservation.sh 31 32          # 跑 31 与 32
#
# 前置 / Prerequisites:
#   1) code_v3/ 全部脚本已同步到服务器
#   2) data/derived_v3/psi_samples_thinned_<RUN_LABEL>.rds 存在且为 4 维
#      （draws × species × sites × periods；若为 3 维需先重跑 05）
#   3) data/external/protected_areas/ 已同步（脚本 31/32 需要）
#   4) R 包：sf, dplyr, readr, matrixStats, MatchIt, prioritizr, highs
# ============================================================
set -uo pipefail

# ── 路径与运行标签（按需修改）────────────────────────────────────────
PROJ="${V3_PROJ_DIR:-$HOME/Documents/New project/bird_dynamic_occupancy_analysis}"
export V3_CODE_DIR="${V3_CODE_DIR:-$PROJ/code_v3}"
export V3_RUN_LABEL="${V3_RUN_LABEL:-v3_full_200sp_ar1_spatial}"
export V3_OMP_THREADS="${V3_OMP_THREADS:-8}"

LOG_DIR="$PROJ/logs_v3"; mkdir -p "$LOG_DIR"
STAMP="$(date +%Y%m%d_%H%M%S)"

echo "=============================================="
echo " 项目目录 : $PROJ"
echo " 代码目录 : $V3_CODE_DIR"
echo " 运行标签 : $V3_RUN_LABEL"
echo " 日志目录 : $LOG_DIR"
echo "=============================================="

# ── 前置检查：缺什么先说清楚，别跑到一半才失败 ──────────────────────
echo ""; echo "[检查] 环境与输入 ..."
command -v Rscript >/dev/null 2>&1 || { echo "  ✗ 未找到 Rscript"; exit 1; }
echo "  ✓ R: $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"

PSI="$PROJ/data/derived_v3/psi_samples_thinned_${V3_RUN_LABEL}.rds"
if [ -f "$PSI" ]; then
  DIMS=$(Rscript -e "x<-readRDS('$PSI'); cat(paste(dim(x\$psi_samples_thinned), collapse='x'))" 2>/dev/null)
  NDIM=$(awk -F'x' '{print NF}' <<< "$DIMS")
  if [ "${NDIM:-0}" -ge 4 ]; then echo "  ✓ psi 4维: $DIMS"
  else echo "  ✗ psi 维度为 $DIMS（需 4 维）。请先重跑 05_postprocess_diversity_extended.R"; exit 1; fi
else
  echo "  ✗ 未找到 $PSI —— 请先跑 04/05 或确认 RUN_LABEL"; exit 1
fi

PA_SHP="$PROJ/data/external/protected_areas/全国自然保护区名录+矢量边界/保护区.shp"
if [ -f "$PA_SHP" ]; then echo "  ✓ 保护区矢量已就位"
else echo "  ⚠ 未找到保护区矢量，脚本 31/32 会失败；请同步 data/external/protected_areas/"; fi

echo "  · 检查 R 包 ..."
Rscript -e '
pk <- c("sf","dplyr","readr","matrixStats","MatchIt","prioritizr","highs","tidyr")
miss <- pk[!sapply(pk, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  cat("  ⚠ 缺少 R 包:", paste(miss, collapse=", "), "\n")
  cat("    安装命令: install.packages(c(", paste(sprintf("\"%s\"", miss), collapse=","), "))\n")
} else cat("  ✓ R 包齐全\n")' 2>/dev/null

# ── 依次执行 ────────────────────────────────────────────────────────
declare -A SCRIPTS=(
  [30]="30_range_expansion_mechanism.R|物种数增长的机制归因（扩张vs填充、重心位移、CTI、类群归因）"
  [31]="31_protected_area_effectiveness.R|保护成效（匹配ATT+双重差分）与保护空缺"
  [32]="32_conservation_prioritization.R|系统保护规划（三方案对比与扩建缺口）"
)
TARGETS=("$@"); [ ${#TARGETS[@]} -eq 0 ] && TARGETS=(30 31 32)

FAILED=()
for id in "${TARGETS[@]}"; do
  entry="${SCRIPTS[$id]:-}"
  [ -z "$entry" ] && { echo "跳过未知脚本编号: $id"; continue; }
  f="${entry%%|*}"; desc="${entry##*|}"
  LOG="$LOG_DIR/${id}_${STAMP}.log"
  echo ""; echo "──────────────────────────────────────────────"
  echo "▶ [$id] $desc"
  echo "  脚本: $f"; echo "  日志: $LOG"
  START=$(date +%s)
  if Rscript "$V3_CODE_DIR/$f" > "$LOG" 2>&1; then
    echo "  ✓ 完成（$(( $(date +%s) - START )) 秒）"
    grep -E "^\[3[012]\] >>>" "$LOG" | sed 's/^/    /' || true
  else
    echo "  ✗ 失败 —— 末尾日志："
    tail -15 "$LOG" | sed 's/^/    /'
    FAILED+=("$id")
  fi
done

# ── 汇总 ────────────────────────────────────────────────────────────
echo ""; echo "=============================================="
if [ ${#FAILED[@]} -eq 0 ]; then echo " 全部完成"; else echo " 失败脚本: ${FAILED[*]}"; fi
echo " 新产出结果表："
ls -lt "$PROJ/results_v3"/table_{range,colonization,centroid,cti,guild,pa,conservation,priority}_*.csv 2>/dev/null \
  | head -15 | awk '{printf "   %8.1f KB  %s\n", $5/1024, $9}'
echo "=============================================="
