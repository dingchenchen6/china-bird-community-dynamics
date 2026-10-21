#!/usr/bin/env bash
# ============================================================
# 科学问题 / Question: 按修正后的数据与代码完整重跑群落动态分析，结果能否通过投稿前 Gate 0–5？
# 目标 / Objective: 分阶段、可断点续跑地执行全部重拟合与重分析；输出隔离在 BIRD_OUTPUT_ROOT，不覆盖 v3。
# 用法 / Usage:
#   export BIRD_PROJECT_ROOT=~/Documents/New\ project/bird_dynamic_occupancy_analysis
#   export BIRD_OUTPUT_ROOT=$BIRD_PROJECT_ROOT/refit_v4      # 新目录 / new directory
#   export V3_CODE_DIR=$BIRD_PROJECT_ROOT/code_v4_refit
#   bash code_v4_refit/run_refit.sh stage0      # 依次 stage0 → stage1 → stage2 / stage3 → stage4，或 all
# 栅格路径在服务器上不同时，设置 REFIT_WORLDCLIM_DIR、REFIT_ENV_EXTERNAL_DIR、REFIT_HFI_DIR（见 S0_01）。
# ============================================================
set -euo pipefail
STAGE="${1:-all}"
INPUT_ROOT="${BIRD_PROJECT_ROOT:?Set BIRD_PROJECT_ROOT}"
OUTPUT_ROOT="${BIRD_OUTPUT_ROOT:?Set BIRD_OUTPUT_ROOT (new directory, e.g. \$BIRD_PROJECT_ROOT/refit_v4)}"
CODE_DIR="${V3_CODE_DIR:?Set V3_CODE_DIR to code_v4_refit}"
N_CHAINS="${V3_N_CHAINS:-4}"
TOTAL_CORES="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 8)"
OMP_PER_CHAIN="${V3_N_OMP_PER_CHAIN:-$(( TOTAL_CORES / N_CHAINS ))}"; [[ "$OMP_PER_CHAIN" -lt 1 ]] && OMP_PER_CHAIN=1
export BIRD_PROJECT_ROOT="$INPUT_ROOT" BIRD_OUTPUT_ROOT="$OUTPUT_ROOT" V3_CODE_DIR="$CODE_DIR" V3_N_CHAINS="$N_CHAINS"
L200="refit_20261020_200sp_spatial_base"; L200S="refit_20261020_200sp_spatial_source"; L500="refit_20261020_500sp_temporal"
mkdir -p "$OUTPUT_ROOT"/{logs_v3,data/derived_v3,results_v3,figures_v3,status}

stamp() { date '+%Y-%m-%dT%H:%M:%S%z'; }
# 单个脚本：已完成则跳过（断点续跑）/ run one script, skip if already done
run_r() {
  local tag="$1" script="$2"; local done="$OUTPUT_ROOT/status/${tag}.done"
  [[ -f "$done" ]] && { echo "[$(stamp)] SKIP  $tag（已完成）"; return 0; }
  local log="$OUTPUT_ROOT/logs_v3/${tag}_$(date +%Y%m%d_%H%M%S).log"
  echo "[$(stamp)] START $tag" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
  if Rscript "$CODE_DIR/$script" >"$log" 2>&1; then touch "$done"; echo "[$(stamp)] DONE  $tag" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
  else echo "[$(stamp)] FAIL  $tag → $log" | tee -a "$OUTPUT_ROOT/status/pipeline.log"; return 1; fi
}
validate_chain() { [[ -s "$1" ]] && Rscript -e 'x<-readRDS(commandArgs(TRUE)[1]); d<-dim(x$psi.samples); stopifnot(length(d)==4L,all(d>0))' "$1" >/dev/null 2>&1; }
# 多链并行拟合，每链绑定一组 CPU，已有有效链则续跑 / parallel chains with resume
run_chains() {
  local tag="$1" prefix="$2" script="$3" label="$4" nsp="$5" source_flag="$6"
  export V3_RUN_LABEL="$label" V3_MAX_SPECIES="$nsp" V3_INCLUDE_SOURCE="$source_flag"
  local pids=()
  for chain in $(seq 1 "$N_CHAINS"); do
    local out="$OUTPUT_ROOT/data/derived_v3/${prefix}_fit_${label}_chain${chain}.rds"
    validate_chain "$out" && { echo "[$(stamp)] RESUME $tag chain $chain"; continue; }
    [[ -e "$out" ]] && mv "$out" "$out.invalid.$(date +%Y%m%d_%H%M%S)"
    local log="$OUTPUT_ROOT/logs_v3/${tag}_chain${chain}_$(date +%Y%m%d_%H%M%S).log"
    local c0=$(( (chain - 1) * OMP_PER_CHAIN )); local c1=$(( c0 + OMP_PER_CHAIN - 1 )); [[ "$c1" -ge "$TOTAL_CORES" ]] && c1=$(( TOTAL_CORES - 1 ))
    echo "[$(stamp)] START $tag chain $chain cores ${c0}-${c1}" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
    if command -v taskset >/dev/null 2>&1; then
      V3_CHAIN_ID="$chain" V3_N_OMP_THREADS="$OMP_PER_CHAIN" OMP_NUM_THREADS="$OMP_PER_CHAIN" taskset -c "${c0}-${c1}" Rscript "$CODE_DIR/$script" >"$log" 2>&1 &
    else
      V3_CHAIN_ID="$chain" V3_N_OMP_THREADS="$OMP_PER_CHAIN" OMP_NUM_THREADS="$OMP_PER_CHAIN" Rscript "$CODE_DIR/$script" >"$log" 2>&1 &
    fi
    pids+=("$!")
  done
  local failed=0; for p in "${pids[@]:-}"; do [[ -z "$p" ]] || wait "$p" || failed=$((failed + 1)); done
  [[ "$failed" -eq 0 ]] || { echo "$tag: $failed 条链失败" >&2; return 1; }
  for chain in $(seq 1 "$N_CHAINS"); do validate_chain "$OUTPUT_ROOT/data/derived_v3/${prefix}_fit_${label}_chain${chain}.rds" || { echo "$tag chain $chain 校验失败" >&2; return 1; }; done
  V3_MODEL_PREFIX="$prefix" run_r "${tag}_combine" "26_combine_validate_chains.R"
}
use() { export V3_RUN_LABEL="$1" V3_MAX_SPECIES="$2" V3_INCLUDE_SOURCE="${3:-0}"; }

stage0() {   # 数据冻结 Gate 0
  run_r S0_01_environment S0_01_rebuild_grid_environment.R
  REFIT_STAGE_ONLY=1 run_r S0_02_stage_inputs S0_02_check_inputs_and_freeze.R
  run_r S0_build_survey_4d 02_build_survey_history.R
  run_r S0_prepare_environment 03_prepare_environment.R
  run_r S0_02_check_freeze S0_02_check_inputs_and_freeze.R
  run_r S0_preflight 27_submission_preflight.R
}
stage1() {   # 主模型重拟合 Gate 1（最耗时）
  run_chains M200_base stMsPGOcc 04_run_stMsPGOcc_main.R "$L200" 200 0
  run_chains M500_temporal tMsPGOcc 04_run_tMsPGOcc_500sp.R "$L500" 500 0
}
stage2() {   # 观测过程检验 Gate 2
  run_chains M200_source stMsPGOcc 04_run_stMsPGOcc_main.R "$L200S" 200 1
  export V3_BASE_LABEL="$L200" V3_SOURCE_LABEL="$L200S"; run_r S2_compare_source 29_compare_source_detection.R
  use "$L200" 200; run_r S2_ppc 05c_ppc_bayesian_pvalue.R
  # 三年窗口：在独立目录中按 3 年期重建调查历史并重拟合 / 3-year window in its own directory
  ( export BIRD_OUTPUT_ROOT="${OUTPUT_ROOT}_3yr" V3_PERIOD_LENGTH=3; OUTPUT_ROOT="$BIRD_OUTPUT_ROOT"
    mkdir -p "$OUTPUT_ROOT"/{logs_v3,data/derived_v3,results_v3,figures_v3,status}
    cp -n "${BIRD_OUTPUT_ROOT%_3yr}/data/derived_v3/grid_environment_v3.rds" "$OUTPUT_ROOT/data/derived_v3/" 2>/dev/null || true
    REFIT_STAGE_ONLY=1 run_r S2_3yr_stage S0_02_check_inputs_and_freeze.R
    run_r S2_3yr_survey 02_build_survey_history.R; run_r S2_3yr_env 03_prepare_environment.R
    run_chains M200_3yr stMsPGOcc 04_run_stMsPGOcc_main.R refit_20261020_200sp_spatial_3yr 200 0
    use refit_20261020_200sp_spatial_3yr 200; run_r S2_3yr_metrics 05_postprocess_diversity_extended.R )
}
stage3() {   # 下游重算 Gate 3–4（依赖 stage1；可与 stage2 并行）
  use "$L200" 200
  run_r S3_metrics200 05_postprocess_diversity_extended.R
  run_r S3_trend_draws200 05e_export_trend_draws.R
  run_r S3_functional_pdecline200 05f_functional_trend_pdecline.R
  run_r S3_phylo_ses200 05g_phylo_ses.R
  run_r S3_homogenization200 12_homogenization_spatiotemporal_maps.R
  run_r S3_drivers200 09_extended_analyses_extended.R
  run_r S3_rf_spatialcv200 09b_rf_spatial_block_cv.R
  run_r S3_traits200 14_species_trait_regression.R
  run_r S3_mechanism200 30_range_expansion_mechanism.R
  run_r S3_pa_effect200 31_protected_area_effectiveness.R
  run_r S3_prioritization200 32_conservation_prioritization.R
  run_r S3_pa_contrast200 33_pa_spatiotemporal_contrast.R
  run_r S3_effort_gain200 S2_04_effort_gain_test.R
  use "$L500" 500
  run_r S3_metrics500 05_postprocess_diversity_extended.R
  run_r S3_functional_pdecline500 05f_functional_trend_pdecline.R
  run_r S3_traits500 14_species_trait_regression.R
}
stage4() {   # 图件与收尾
  use "$L200" 200
  for s in 24b_main_figures_standard_map_20261011.R 24c_figure1_layout_20261011.R 24d_figures4to6_20261011.R \
           21b_supporting_figures_native_pptx_20261011.R 34b_pa_figures_english_standard_map_20261011.R 34c_pa_balance_and_maps_20261011.R; do
    run_r "S4_${s%%_*}" "$s"; done
  python3 "$CODE_DIR/35_qa_figure_package_20261011.py" > "$OUTPUT_ROOT/logs_v3/S4_figure_qa.log" 2>&1 || echo "图件质检有问题，见 S4_figure_qa.log"
  run_r S4_postflight 28_submission_postflight.R
  Rscript -e 'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("BIRD_OUTPUT_ROOT"), "results_v3", "sessionInfo.txt"))'
  touch "$OUTPUT_ROOT/status/REFIT_COMPLETE"
}
case "$STAGE" in
  stage0|stage1|stage2|stage3|stage4) "$STAGE" ;;
  all) stage0; stage1; stage2; stage3; stage4 ;;
  *) echo "用法: bash run_refit.sh stage0|stage1|stage2|stage3|stage4|all" >&2; exit 2 ;;
esac
