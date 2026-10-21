#!/usr/bin/env bash
# Offline, resumable submission rerun. Upload all inputs before launching.
set -euo pipefail

INPUT_ROOT="${BIRD_PROJECT_ROOT:?Set BIRD_PROJECT_ROOT to the immutable input project}"
OUTPUT_ROOT="${BIRD_OUTPUT_ROOT:?Set BIRD_OUTPUT_ROOT to the versioned rerun directory}"
CODE_DIR="${V3_CODE_DIR:?Set V3_CODE_DIR to the uploaded code_v3 directory}"
N_CHAINS="${V3_N_CHAINS:-4}"
TOTAL_CORES="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 8)"
OMP_PER_CHAIN="${V3_N_OMP_PER_CHAIN:-$(( TOTAL_CORES / N_CHAINS ))}"
[[ "$OMP_PER_CHAIN" -lt 1 ]] && OMP_PER_CHAIN=1
mkdir -p "$OUTPUT_ROOT/logs_v3" "$OUTPUT_ROOT/data/derived_v3" \
  "$OUTPUT_ROOT/results_v3" "$OUTPUT_ROOT/figures_v3" "$OUTPUT_ROOT/status"
export BIRD_PROJECT_ROOT="$INPUT_ROOT" BIRD_OUTPUT_ROOT="$OUTPUT_ROOT" V3_CODE_DIR="$CODE_DIR"
export V3_N_CHAINS="$N_CHAINS"

stamp() { date '+%Y-%m-%dT%H:%M:%S%z'; }
run_r() {
  local tag="$1" script="$2"
  local log="$OUTPUT_ROOT/logs_v3/${tag}_$(date +%Y%m%d_%H%M%S).log"
  echo "[$(stamp)] START $tag" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
  Rscript "$CODE_DIR/$script" >"$log" 2>&1
  echo "[$(stamp)] DONE  $tag" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
  touch "$OUTPUT_ROOT/status/${tag}.done"
}
validate_chain() {
  local f="$1"
  [[ -s "$f" ]] && Rscript -e 'x<-readRDS(commandArgs(TRUE)[1]); d<-dim(x$psi.samples); stopifnot(length(d)==4L,all(d>0))' "$f" >/dev/null 2>&1
}
run_chains() {
  local tag="$1" prefix="$2" script="$3" label="$4" nsp="$5" source_flag="$6"
  export V3_RUN_LABEL="$label" V3_MAX_SPECIES="$nsp" V3_INCLUDE_SOURCE="$source_flag"
  local pids=() chains=()
  for chain in $(seq 1 "$N_CHAINS"); do
    local out="$OUTPUT_ROOT/data/derived_v3/${prefix}_fit_${label}_chain${chain}.rds"
    if validate_chain "$out"; then
      echo "[$(stamp)] RESUME $tag chain $chain already valid" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
      continue
    fi
    if [[ -e "$out" ]]; then mv "$out" "$out.invalid.$(date +%Y%m%d_%H%M%S)"; fi
    local log="$OUTPUT_ROOT/logs_v3/${tag}_chain${chain}_$(date +%Y%m%d_%H%M%S).log"
    local cpu_start=$(( (chain - 1) * OMP_PER_CHAIN ))
    local cpu_end=$(( cpu_start + OMP_PER_CHAIN - 1 ))
    [[ "$cpu_end" -ge "$TOTAL_CORES" ]] && cpu_end=$(( TOTAL_CORES - 1 ))
    echo "[$(stamp)] START $tag chain $chain cores ${cpu_start}-${cpu_end}" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
    if command -v taskset >/dev/null 2>&1; then
      V3_CHAIN_ID="$chain" V3_N_OMP_THREADS="$OMP_PER_CHAIN" OMP_NUM_THREADS="$OMP_PER_CHAIN" \
        taskset -c "${cpu_start}-${cpu_end}" Rscript "$CODE_DIR/$script" >"$log" 2>&1 &
    else
      V3_CHAIN_ID="$chain" V3_N_OMP_THREADS="$OMP_PER_CHAIN" OMP_NUM_THREADS="$OMP_PER_CHAIN" \
        Rscript "$CODE_DIR/$script" >"$log" 2>&1 &
    fi
    pids+=("$!"); chains+=("$chain")
  done
  local failed=0
  for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then failed=$((failed + 1)); fi
  done
  [[ "$failed" -eq 0 ]] || { echo "$tag: $failed chain processes failed" >&2; return 1; }
  for chain in $(seq 1 "$N_CHAINS"); do
    validate_chain "$OUTPUT_ROOT/data/derived_v3/${prefix}_fit_${label}_chain${chain}.rds" || {
      echo "$tag chain $chain failed checkpoint validation" >&2; return 1; }
  done
  V3_MODEL_PREFIX="$prefix" run_r "${tag}_combine" "26_combine_validate_chains.R"
}

echo "[$(stamp)] Submission rerun begins" | tee "$OUTPUT_ROOT/status/pipeline.log"
run_r preflight_initial 27_submission_preflight.R
run_r build_survey_4d 02_build_survey_history.R
run_r prepare_environment 03_prepare_environment.R
run_r preflight_canonical 27_submission_preflight.R

# Reuse the frozen offline trait table when available; 03b otherwise uses local fallbacks.
if [[ -f "$INPUT_ROOT/data/derived_v3/trait_extended_v3.rds" ]]; then
  cp -p "$INPUT_ROOT/data/derived_v3/trait_extended_v3.rds" "$OUTPUT_ROOT/data/derived_v3/"
else
  unset IUCN_API_KEY
  run_r prepare_traits_offline 03b_extend_traits.R
fi

run_chains model200_base stMsPGOcc 04_run_stMsPGOcc_main.R \
  submission_20260722_200sp_spatial_base 200 0
run_chains model200_source stMsPGOcc 04_run_stMsPGOcc_main.R \
  submission_20260722_200sp_spatial_source 200 1
export V3_BASE_LABEL=submission_20260722_200sp_spatial_base
export V3_SOURCE_LABEL=submission_20260722_200sp_spatial_source
run_r compare_source_detection 29_compare_source_detection.R
run_chains model500_temporal tMsPGOcc 04_run_tMsPGOcc_500sp.R \
  submission_20260722_500sp_temporal 500 0

export V3_RUN_LABEL=submission_20260722_500sp_temporal V3_MAX_SPECIES=500 V3_INCLUDE_SOURCE=0
run_r metrics500_extended 05_postprocess_diversity_extended.R
run_r trend_probability500 05f_functional_trend_pdecline.R
run_r homogenization500 12_homogenization_spatiotemporal_maps.R
run_r trait_regression500 14_species_trait_regression.R
run_r figures500_all_formats 22_make_500sp_all_analysis_figures.R
run_r submission_postflight 28_submission_postflight.R
touch "$OUTPUT_ROOT/status/SUBMISSION_RERUN_COMPLETE"
echo "[$(stamp)] Submission rerun complete" | tee -a "$OUTPUT_ROOT/status/pipeline.log"
