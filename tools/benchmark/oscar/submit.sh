#!/bin/bash
# Run only after the actual checkquota/sinfo/sacctmgr preflight has been reviewed.
set -euo pipefail
root=${1:?new scratch campaign root required}
budget=${2:?selected campaign disk budget in GiB required}
harness=$(cd -- "$(dirname -- "$0")" && pwd)
scratch=$(realpath "${3:?explicit temporary-storage base required}")
root=$(realpath -m "$root")
[[ "$root" == "$scratch/"* ]]
[[ ! -e "$root" ]]
[[ -O "$scratch" && -w "$scratch" ]]
[[ "$budget" =~ ^[1-9][0-9]*$ ]]
command -v samtools
command -v sbatch
mkdir -p "$root/logs"
checkquota > "$root/logs/quota-before.txt"
sinfo -p batch -o '%P %a %l %D %c %m %f' > "$root/logs/partition-before.txt"
python3 - "$root" "$budget" "$harness" <<'PY'
import hashlib,json,sys
from pathlib import Path
root,budget,harness=sys.argv[1:];source=Path(harness)
(Path(root)/'campaign.json').write_text(json.dumps({'root':root,'disk_budget_gib':int(budget),
 'phase_walltime':'24:00:00','harness':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir() if p.is_file()}},indent=2)+'\n')
PY
# Each STAR/NuSTAR pair and its repetitions stay on one allocated node.
job1=$(sbatch --parsable --chdir="$root" --output="$root/logs/main-%j.log" "$harness/campaign.sbatch" "$root" "$budget" main "$harness")
job1=${job1%%;*}
printf '%s\n' "$job1" > "$root/job-ids.txt"
job2=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$job1" --chdir="$root" --output="$root/logs/full-bulk-%j.log" "$harness/campaign.sbatch" "$root" "$budget" full-bulk "$harness")
job2=${job2%%;*}
printf '%s\n' "$job2" >> "$root/job-ids.txt"
job3=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$job2" --chdir="$root" --output="$root/logs/full-pbmc-%j.log" "$harness/campaign.sbatch" "$root" "$budget" full-pbmc "$harness")
job3=${job3%%;*}
printf '%s\n' "$job3" >> "$root/job-ids.txt"
printf 'Campaign: %s\nJobs: %s %s %s\n' "$root" "$job1" "$job2" "$job3"
