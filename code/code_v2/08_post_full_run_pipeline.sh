#!/usr/bin/env bash
# 08_post_full_run_pipeline.sh
#
# 在 stage-4 full run 拟合完成（saved tMsPGOcc_fit_<LABEL>.rds + psi_samples_thinned_<LABEL>.rds）之后
# 一键跑：stage 5 → 6 → 6b → 7。
#
# 默认环境变量：V2_RUN_LABEL=v2_full_200sp_ar1
# 用法：bash code_v2/08_post_full_run_pipeline.sh
#       或 V2_RUN_LABEL=xxx bash code_v2/08_post_full_run_pipeline.sh

set -euo pipefail

cd "$(dirname "$0")/.."     # repo root
mkdir -p logs_v2

LABEL="${V2_RUN_LABEL:-v2_full_200sp_ar1}"
DRAWS="${V2_POST_DRAWS:-200}"          # full run 时一次性拿 200 draws
TS="$(date +%Y%m%d_%H%M%S)"

echo "[pipeline] running stages 5 -> 6 -> 6b -> 7 for ${LABEL}"

run_stage () {
  local script="$1"; local stem="$2"
  local LOG="logs_v2/${stem}_${LABEL}_${TS}.log"
  echo "[pipeline] ($stem) -> $LOG"
  R_MAX_VSIZE=160Gb V2_RUN_LABEL="${LABEL}" V2_POST_DRAWS="${DRAWS}" \
    Rscript "${script}" >"${LOG}" 2>&1
  echo "[pipeline] ($stem) done"
}

run_stage "code_v2/05_postprocess_diversity.R"     "stage5"
run_stage "code_v2/06_figures_publication.R"       "stage6"
run_stage "code_v2/06b_regenerate_driver_plots.R"  "stage6b"
run_stage "code_v2/07_render_manuscript.R"         "stage7"

echo "[pipeline] All stages done for ${LABEL}."
echo "  Outputs: results_v2/, figures_v2/"
ls "results_v2/" | grep "${LABEL}" | head -20
echo "  Figures (newest 10):"
ls -t figures_v2/*.png | head -10
