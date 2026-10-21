#!/usr/bin/env bash
# ============================================================
# 无人值守连续重跑驱动 / unattended pipeline driver
#   waiting for the running stage1 workload → stage1 收尾 → stage2 → stage3 → stage4
#
# 依赖 run_refit.sh 的 status/*.done 断点标记；每阶段失败自动重试 1 次。
# 阶段依赖：stage2/stage3 只依赖 stage1（可先后执行）；stage4 依赖 stage3。
# 输出隔离在 BIRD_OUTPUT_ROOT，不触碰 v3。
# 用法：
#   export BIRD_PROJECT_ROOT=<项目根> BIRD_OUTPUT_ROOT=<新输出根> V3_CODE_DIR=<code_v4_refit>
#   nohup setsid bash $V3_CODE_DIR/run_refit_auto.sh >> $BIRD_OUTPUT_ROOT/logs_v3/auto_driver.bootstrap.log 2>&1 &
# ============================================================
set -uo pipefail

: "${BIRD_PROJECT_ROOT:?Set BIRD_PROJECT_ROOT}"
: "${BIRD_OUTPUT_ROOT:?Set BIRD_OUTPUT_ROOT}"
: "${V3_CODE_DIR:?Set V3_CODE_DIR to code_v4_refit}"

OUT="$BIRD_OUTPUT_ROOT"
CODE="$V3_CODE_DIR"
mkdir -p "$OUT/status" "$OUT/logs_v3"
DRIVER_LOG="$OUT/logs_v3/auto_driver.log"
export V3_N_CHAINS="${V3_N_CHAINS:-4}"
export V3_N_OMP_PER_CHAIN="${V3_N_OMP_PER_CHAIN:-32}"

stamp() { date '+%Y-%m-%dT%H:%M:%S%z'; }
say()   { echo "[$(stamp)] $*" | tee -a "$DRIVER_LOG"; }

busy() {
  pgrep -f "run_refit\.sh stage[0-9]"      >/dev/null 2>&1 && return 0
  pgrep -f "04_run_stMsPGOcc_main\.R"      >/dev/null 2>&1 && return 0
  pgrep -f "04_run_tMsPGOcc_500sp\.R"      >/dev/null 2>&1 && return 0
  pgrep -f "26_combine_validate_chains\.R" >/dev/null 2>&1 && return 0
  return 1
}

run_stage() {  # stage, retries
  local stage="$1" retries="${2:-1}" n=0 rc=1
  while (( n <= retries )); do
    say "START $stage (attempt $((n+1))/$((retries+1)))"
    bash "$CODE/run_refit.sh" "$stage" >> "$OUT/logs_v3/${stage}_driver.log" 2>&1
    rc=$?
    if [[ $rc -eq 0 ]]; then say "DONE  $stage"; return 0; fi
    say "FAIL  $stage (exit $rc) — see logs_v3/${stage}_driver.log"
    n=$((n+1)); (( n <= retries )) && sleep 120
  done
  say "GIVEUP $stage"
  return 1
}

say "AUTO driver started (pid $$); project=$BIRD_PROJECT_ROOT output=$OUT code=$CODE"

# ── 1. 等待已在运行的 stage1 结束 / wait for the running stage1 workload ──
n=0
while busy; do
  n=$((n+1))
  if (( n % 12 == 1 )); then
    say "waiting: $(pgrep -fc '04_run_stMsPGOcc_main\.R' 2>/dev/null || echo 0) stMsPGOcc chain proc(s), $(pgrep -fc '04_run_tMsPGOcc_500sp\.R' 2>/dev/null || echo 0) tMsPGOcc proc(s)"
  fi
  sleep 300
done
say "stage1 workload idle — supervising the rest of the pipeline"

# ── 2. stage1 收尾（补齐缺失链 + 合并；已完成步骤自动跳过） ──────────────
run_stage stage1 1 || say "stage1 reported failure — validating chain files before deciding"

# ── 3. stage1 产物完整性闸门 / gate on the 200sp chain files ────────────
ok=1
for c in 1 2 3 4; do
  f="$OUT/data/derived_v3/stMsPGOcc_fit_refit_20261020_200sp_spatial_base_chain${c}.rds"
  if [[ ! -s "$f" ]]; then ok=0; say "GATE missing: $(basename "$f")"; continue; fi
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if (( sz < 10485760 )); then ok=0; say "GATE too small ($sz B): $(basename "$f")"; fi
done
if [[ $ok -eq 0 ]]; then
  say "GATE FAILED: M200 chains incomplete — stopping before downstream (a later rerun resumes)"
  echo "$(stamp)" > "$OUT/status/AUTO_PIPELINE_INCOMPLETE"
  exit 1
fi
say "GATE passed: 4 × M200 chain files present and >10 MB"

# ── 3b. 收敛闸门 / convergence gate ─────────────────────────────────────
# v4 2026-10-21: 闸门语义升级为披露制（26_combine_validate_chains.R）：
#   gate_pass = 严格层(R̂≤1.05,ESS≥200) 全过 且 观察名单(beta.comm/tau.sq.beta,
#   R̂≤1.4,ESS≥5) 兜底全过。优先读 gate_pass；旧表无该列时退回 rhat/ess 两列。
GATES=$(ls "$OUT"/results_v3/table_convergence_gate_*.csv 2>/dev/null || true)
if [[ -n "$GATES" ]]; then
  if ! Rscript -e 'f <- commandArgs(TRUE); ok <- all(vapply(f, function(p) { g <- read.csv(p); if ("gate_pass" %in% names(g)) isTRUE(all(g$gate_pass)) else isTRUE(all(g$rhat_pass)) && isTRUE(all(g$ess_pass)) }, TRUE)); quit(status = if (ok) 0L else 1L)' $GATES; then
    say "CONVERGENCE GATE FAILED — not converging models must be refitted before downstream analyses"
    echo "$(stamp)" > "$OUT/status/AUTO_PIPELINE_PARTIAL"
    exit 1
  fi
  say "CONVERGENCE GATE passed: $(echo "$GATES" | wc -l | tr -d ' ') model(s)"
else
  say "WARNING: no convergence gate table found — proceeding (combine step may have been skipped)"
fi

# ── 4. 下游阶段 / downstream stages ────────────────────────────────────
STATUS=0
run_stage stage2 1 || say "stage2 failed twice — continuing (stage3 does not depend on stage2)"
run_stage stage3 1 || { say "stage3 failed twice"; STATUS=1; }
if [[ $STATUS -eq 0 ]]; then
  run_stage stage4 1 || { say "stage4 failed twice"; STATUS=1; }
else
  say "stage4 skipped: stage3 did not complete"
fi

if [[ $STATUS -eq 0 ]]; then
  echo "$(stamp)" > "$OUT/status/AUTO_PIPELINE_COMPLETE"
  say "AUTO PIPELINE COMPLETE"
else
  echo "$(stamp)" > "$OUT/status/AUTO_PIPELINE_PARTIAL"
  say "AUTO PIPELINE ENDED WITH FAILURES — inspect stage logs, rerun is resumable"
fi
exit $STATUS
