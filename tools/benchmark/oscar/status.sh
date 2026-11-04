#!/bin/bash
set -euo pipefail
root=${1:?campaign root required}
[[ -f "$root/campaign.json" ]]
ids=$(paste -sd, "$root/job-ids.txt")
squeue -u "$(id -un)" -h -o '%i %P %j %T %M %l %R' | awk -v ids="$ids" 'BEGIN {n=split(ids,a,","); for(i=1;i<=n;i++) wanted[a[i]]=1} $1 in wanted'
sacct -j "$ids" --format=JobID,State,ExitCode,Elapsed,AllocCPUS,MaxRSS,NodeList --parsable2
python3 - "$root" <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1])
for p in sorted(root.glob('logs/job-*-state.json'))+sorted(root.glob('results/*/status.json')):
 print(str(p.relative_to(root)),json.dumps(json.loads(p.read_text()),sort_keys=True))
PY
