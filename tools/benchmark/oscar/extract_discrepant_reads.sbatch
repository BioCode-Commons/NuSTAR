#!/bin/bash
#SBATCH --job-name=nustar-read-diff
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -euo pipefail
root=${1:?campaign root required}
harness=${2:?harness directory required}
module purge
module load gcc/15.1.0-hhh7 cmake/3.31.8-i66z python/3.11.11-5e66 samtools/1.19.2-fm2x
export TMPDIR="$root/tmp"
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export LC_ALL=C
export PYTHONUNBUFFERED=1
srun --cpu-bind=cores python3 "$harness/extract_discrepant_reads.py" "$root"
