#!/bin/bash
set -euo pipefail
sinfo -p general -o '%D %C %l'
sacctmgr -n show assoc user="$USER" format=Account,Partition,GrpCPUs,GrpTRESMins,MaxJobs,MaxSubmitJobs
scontrol show config | /usr/bin/grep -E 'MaxArraySize|MaxJobCount'
df -h /N/slate/tg11/ebmcal_v2
du -sh /N/slate/tg11/ebmcal_v2/runs
