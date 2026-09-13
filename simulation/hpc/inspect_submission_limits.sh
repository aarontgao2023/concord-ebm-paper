#!/bin/bash
set -euo pipefail
sacctmgr -n -P show assoc where user=tg11 account=r00579 format=User,Account,Partition,QOS,DefaultQOS,MaxJobs,MaxSubmitJobs,MaxTRES
sacctmgr -n -P show qos format=Name,MaxJobsPU,MaxSubmitPU,MaxTRESPU,MaxWall
squeue -r -u tg11 -h -o '%i|%T|%k'
sacct -X -n -P -u tg11 -S 2026-09-08T10:54:00 -o JobID%80,State%40,Comment%100,Submit,Start,End
