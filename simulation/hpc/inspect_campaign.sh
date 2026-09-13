#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
"$HOME/ebmcal-env/bin/python" "$root/toolkit/hpc/v2/collect_queue.py" --registry "$root/control/job_registry.json" --ledger "$root/control/submission_ledger.json" --output "$root/control/queue_snapshot.json"
"$HOME/ebmcal-env/bin/python" "$root/toolkit/hpc/v2/manage_campaign.py" --campaign "$root/control/campaign.json" --queue "$root/control/queue_snapshot.json" --ledger "$root/control/submission_ledger.json" --weekly-used-percent 21 --output "$root/control/campaign_plan.json" --markdown "$root/control/campaign_plan.md"
