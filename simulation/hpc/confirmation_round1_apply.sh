#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
py="$HOME/ebmcal-env/bin/python"
proposal="$root/control/campaign_versions/confirmation_20260908_1455"
"$py" "$root/toolkit/hpc/v2/reconcile_ledger.py" --ledger "$root/control/submission_ledger.json"
"$py" "$root/toolkit/hpc/v2/collect_queue.py" --registry "$root/control/job_registry.json" --ledger "$root/control/submission_ledger.json" --output "$proposal/queue_apply.json"
"$py" "$root/toolkit/hpc/v2/manage_campaign.py" --campaign "$proposal/campaign.json" --queue "$proposal/queue_apply.json" --ledger "$root/control/submission_ledger.json" --weekly-used-percent 30 --output "$proposal/plan_apply.json" --markdown "$proposal/plan_apply.md"
"$py" "$root/toolkit/hpc/v2/dispatch_campaign.py" --plan "$proposal/plan_apply.json" --queue "$proposal/queue_apply.json" --ledger "$root/control/submission_ledger.json" --project-root "$root" --weekly-used-percent 30 --output "$proposal/dispatch_apply.json" --apply
cat "$proposal/dispatch_apply.json"
