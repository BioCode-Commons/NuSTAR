#!/bin/bash
#SBATCH --job-name=nustar-bench
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --signal=B:TERM@120
set -euo pipefail
set -x

root=${1:?campaign root required}
budget=${2:?disk budget in GiB required}
phase=${3:?phase required}
# Slurm executes a spool copy of this file; the original harness path is explicit.
harness=${4:?original harness directory required}
[[ -f "$root/campaign.json" ]]
[[ -n ${SLURM_JOB_ID:-} ]]
mkdir -p "$root/logs" "$root/tmp"
export TMPDIR="$root/tmp"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_DYNAMIC=FALSE
export LC_ALL=C
export PYTHONUNBUFFERED=1
step=starting
finish() {
    code=$?
    python3 - "$root" "$phase" "$step" "$code" <<'PY'
import json, os, sys
from pathlib import Path
root, phase, step, code = sys.argv[1:]
path=Path(root)/'logs'/f'job-{os.environ["SLURM_JOB_ID"]}-state.json'
path.write_text(json.dumps({'phase':phase,'last_step':step,'exit_code':int(code),
                            'state':'complete' if code=='0' else 'failed'},indent=2)+'\n')
PY
}
trap finish EXIT
trap 'exit 143' TERM
prep=(python3 "$harness/prepare_inputs.py" --sources "$harness/sources.json" --root "$root" --disk-budget-gib "$budget")
runner=(python3 "$harness/run_suite.py" --manifest "$harness/suite.json" --root "$root" --disk-budget-gib "$budget" --prune-validated)
star="$root/source/STAR/source/STAR"
nustar="$root/build/nustar/nustar"
star_index="$root/index/star"
nustar_index="$root/index/nustar"
samtools=$(command -v samtools)
runner+=(--star "$star" --nustar "$nustar" --star-index "$star_index" --nustar-index "$nustar_index" --samtools "$samtools")

case "$phase" in
  main)
    step=build
    bash "$harness/setup.sh" "$root" "$harness"
    step=reference
    "${prep[@]}" --dataset reference
    step=indexes
    python3 "$harness/build_indexes.py" --root "$root" --threads "$SLURM_CPUS_PER_TASK" --disk-budget-gib "$budget"
    step=bulk-subsets
    "${prep[@]}" --dataset bulk
    step=pbmc-subsets
    "${prep[@]}" --dataset pbmc
    step=pilot
    srun --cpu-bind=cores "${runner[@]}" --stage pilot --threads 1 32 --repeats 1
    step=main
    srun --cpu-bind=cores "${runner[@]}" --stage main --threads 1 8 32 --repeats 3
    ;;
  full-bulk|full-pbmc)
    dataset=${phase#full-}
    step=full-input
    "${prep[@]}" --dataset "$dataset" --full
    step=full-benchmark
    srun --cpu-bind=cores "${runner[@]}" --stage full --dataset "$dataset" --threads 32 --repeats 3
    step=release-full-input
    python3 "$harness/release_full_inputs.py" --root "$root" --dataset "$dataset" --job-id "$SLURM_JOB_ID"
    ;;
  *) printf 'Unknown campaign phase: %s\n' "$phase" >&2; exit 2 ;;
esac
step=finished
