#!/bin/bash
set -euo pipefail
scontrol update JobId=10255861 TimeLimit=01:00:00
scontrol update JobId=10258286 TimeLimit=01:00:00
scontrol show job 10255861_0
