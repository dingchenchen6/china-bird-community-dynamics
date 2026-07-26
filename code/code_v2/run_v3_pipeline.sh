#!/usr/bin/env bash
## run_v3_pipeline.sh
## 在 stage 16 完成后自动串行运行 stage 17 → 18。
## 用法：bash code_v2/run_v3_pipeline.sh [v2_full_200sp_ar1]

set -e
cd "$(dirname "$0")/.."
mkdir -p logs_v2

LABEL="${1:-v2_full_200sp_ar1}"
TS="$(date +%Y%m%d_%H%M%S)"
PIPE_LOG="logs_v2/v3_pipeline_${TS}.log"
echo "[v3-pipeline] LABEL=${LABEL}" | tee "$PIPE_LOG"

# 等 stage 16 完成
echo "[v3-pipeline] waiting for stage 16 to finish ..." | tee -a "$PIPE_LOG"
while pgrep -f "16_diversity_v3_phylo_func.R" > /dev/null 2>&1; do
  sleep 60
done
echo "[v3-pipeline] stage 16 process ended" | tee -a "$PIPE_LOG"

# 验证 stage 16 输出存在
V3_METRIC="results_v3/table_community_metrics_v3_${LABEL}.csv"
if [ ! -f "$V3_METRIC" ]; then
  echo "[v3-pipeline] ERROR: $V3_METRIC missing. Aborting." | tee -a "$PIPE_LOG"
  exit 1
fi
echo "[v3-pipeline] stage 16 output verified: $(wc -l < $V3_METRIC) rows" | tee -a "$PIPE_LOG"

# Stage 17
echo "[v3-pipeline] running stage 17 ..." | tee -a "$PIPE_LOG"
S17_LOG="logs_v2/stage17_${LABEL}_${TS}.log"
R_MAX_VSIZE=160Gb V2_RUN_LABEL="$LABEL" \
  Rscript code_v2/17_driver_regression_v3.R > "$S17_LOG" 2>&1
S17_EXIT=$?
echo "[v3-pipeline] stage 17 exit=$S17_EXIT (log: $S17_LOG)" | tee -a "$PIPE_LOG"

# Stage 18
echo "[v3-pipeline] running stage 18 ..." | tee -a "$PIPE_LOG"
S18_LOG="logs_v2/stage18_${LABEL}_${TS}.log"
R_MAX_VSIZE=80Gb V2_RUN_LABEL="$LABEL" \
  Rscript code_v2/18_figures_v3_drivers_maps.R > "$S18_LOG" 2>&1
S18_EXIT=$?
echo "[v3-pipeline] stage 18 exit=$S18_EXIT (log: $S18_LOG)" | tee -a "$PIPE_LOG"

echo "[v3-pipeline] all done. summary:" | tee -a "$PIPE_LOG"
echo "  figures_v3/ count: $(ls figures_v3/ 2>/dev/null | wc -l)" | tee -a "$PIPE_LOG"
echo "  results_v3/ count: $(ls results_v3/ 2>/dev/null | wc -l)" | tee -a "$PIPE_LOG"
