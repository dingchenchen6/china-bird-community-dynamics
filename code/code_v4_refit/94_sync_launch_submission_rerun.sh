#!/usr/bin/env bash
# Sync code plus immutable local inputs, then launch a network-independent server job.
set -euo pipefail
LOCAL_ROOT="${1:-$HOME/Documents/New project/bird_dynamic_occupancy_analysis}"
REMOTE_HOST="${2:-server23}"
REMOTE_INPUT="${3:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3}"
REMOTE_OUTPUT="${4:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_rerun_20260722}"
CODE_SOURCE="${5:-$LOCAL_ROOT/code_v3}"
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=15 -o ServerAliveInterval=15 -o ServerAliveCountMax=3)

ssh "${SSH_OPTS[@]}" "$REMOTE_HOST" "mkdir -p '$REMOTE_INPUT/data/derived_v3' '$REMOTE_INPUT/data/derived_v2' '$REMOTE_OUTPUT/code_v3' '$REMOTE_OUTPUT/logs_v3' '$REMOTE_OUTPUT/status'"
rsync -a --checksum -e "ssh ${SSH_OPTS[*]}" "$CODE_SOURCE/" "$REMOTE_HOST:$REMOTE_OUTPUT/code_v3/"
for rel in \
  data/derived_v3/combined_events_dedup_v3.rds \
  data/derived_v3/trait_extended_v3.rds \
  data/derived_v3/climate_change_v3.rds \
  data/derived_v3/hfi_change_v3.rds \
  data/derived_v2/china_grid_100km_v2.rds \
  data/derived_v2/grid_environment_dynamic_occupancy.rds; do
  rsync -a --checksum --partial -e "ssh ${SSH_OPTS[*]}" "$LOCAL_ROOT/$rel" "$REMOTE_HOST:$REMOTE_INPUT/$(dirname "$rel")/"
done

ssh "${SSH_OPTS[@]}" "$REMOTE_HOST" \
  "BIRD_PROJECT_ROOT='$REMOTE_INPUT' BIRD_OUTPUT_ROOT='$REMOTE_OUTPUT' V3_CODE_DIR='$REMOTE_OUTPUT/code_v3' \
   nohup setsid bash '$REMOTE_OUTPUT/code_v3/93_server_submission_rerun.sh' \
   >'$REMOTE_OUTPUT/logs_v3/submission_rerun_launcher.log' 2>&1 < /dev/null & echo \$! >'$REMOTE_OUTPUT/status/launcher.pid'"
sleep 3
ssh "${SSH_OPTS[@]}" "$REMOTE_HOST" \
  "cat '$REMOTE_OUTPUT/status/launcher.pid'; tail -n 20 '$REMOTE_OUTPUT/logs_v3/submission_rerun_launcher.log'"
