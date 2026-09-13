#!/bin/bash
set -euo pipefail
scontrol update JobId=10255861 ArrayTaskThrottle=4
scontrol update JobId=10258286 ArrayTaskThrottle=4
scontrol show job 10255861_0
scontrol show job 10258286_0
