#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
py="$HOME/ebmcal-env/bin/python"
proposal="$root/control/campaign_versions/confirmation_20260908_1455"
for snapshot in "$root"/snapshots/confirm_{mechanism_b,core_a,stage_a,power_global_a,pair_complete_ref_a,pair_power_k27_a,pair_partial_e2_topup_a,pair_partial_e4_topup_a,pair_partial_e33_a,power_global_k14_k55_a}; do
  (cd "$snapshot" && sha256sum --check --status SHA256SUMS)
done
"$py" "$root/toolkit/hpc/v2/reconcile_ledger.py" --ledger "$root/control/submission_ledger.json"
"$py" "$root/toolkit/hpc/v2/collect_queue.py" --registry "$root/control/job_registry.json" --ledger "$root/control/submission_ledger.json" --output "$proposal/queue_dryrun.json"
"$py" "$root/toolkit/hpc/v2/manage_campaign.py" --campaign "$proposal/campaign.json" --queue "$proposal/queue_dryrun.json" --ledger "$root/control/submission_ledger.json" --weekly-used-percent 30 --output "$proposal/plan_dryrun.json" --markdown "$proposal/plan_dryrun.md"
"$py" "$root/toolkit/hpc/v2/dispatch_campaign.py" --plan "$proposal/plan_dryrun.json" --queue "$proposal/queue_dryrun.json" --ledger "$root/control/submission_ledger.json" --project-root "$root" --weekly-used-percent 30 --output "$proposal/dispatch_dryrun.json"
cat "$proposal/plan_dryrun.md" "$proposal/dispatch_dryrun.json"
