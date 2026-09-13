#!/bin/bash
cd /N/slate/tg11/ebmcal_v2/dev
L=/N/slate/tg11/ebmcal_v2/control/takeover_20260910.jsonl
try() { # $1 mode $2 nchunks $3 datasets $4 array $5 throttle
  local jid
  jid=$(MODE=$1 NCHUNKS=$2 DATASETS=$3 sbatch --parsable --array=$4%$5 standardized_consensus_dev.sbatch 2>/dev/null) || return 1
  echo "{\"utc\":\"$(date -u +%FT%TZ)\",\"run\":\"dev_stdcons_$1\",\"job_id\":\"$jid\",\"n\":$2,\"note\":\"composition-standardised consensus (pooled/min ref) on shared and per-group mixture; $1; R=$3; dev seeds 31800000\"}" >> $L
  echo "$(date -u +%FT%TZ) submitted $1 as $jid"; return 0
}
until try mechanism 20 200 0-19 20; do sleep 600; done
until try calibration 30 300 0-29 16; do sleep 600; done
