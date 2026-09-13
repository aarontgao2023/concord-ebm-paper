#!/bin/bash
set -euo pipefail
scontrol update JobId=10254376 ArrayTaskThrottle=8
scontrol update JobId=10255861 ArrayTaskThrottle=4
scontrol update JobId=10258286 ArrayTaskThrottle=4
squeue -u "$USER" -h -o '%.18i %.16j %.9T %.10M %.6C %R'
