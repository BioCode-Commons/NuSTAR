#!/bin/bash
set -euo pipefail
root=${1:?campaign root required}
harness=${2:?harness directory required}
[[ -f "$root/campaign.json" ]]
mkdir -p "$root/source" "$root/build" "$root/logs"
for tool in git cmake gcc g++ make python3 samtools; do command -v "$tool"; done
python3 -c 'import sys; assert sys.version_info >= (3, 10), "Python 3.10+ required"'
# These are the exact sources being benchmarked. The campaign harness is separate.
clone_exact() {
    url=$1; destination=$2; revision=$3
    if [[ ! -d "$destination/.git" ]]; then git clone "$url" "$destination"; fi
    [[ -z $(git -C "$destination" diff --name-only) ]]
    git -C "$destination" checkout --detach "$revision"
    [[ $(git -C "$destination" rev-parse HEAD) == "$revision" ]]
}
clone_exact https://github.com/BioCode-Commons/NuSTAR.git "$root/source/NuSTAR" 90de97ddce76b1f115f47e37d15add62784bafed
clone_exact https://github.com/alexdobin/STAR.git "$root/source/STAR" b1edc1208d91a53bf40ebae8669f71d50b994851
cmake -S "$root/source/NuSTAR" -B "$root/build/nustar" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DNUSTAR_ENABLE_RELEASE_IPO=ON > "$root/logs/nustar-configure.log" 2>&1
cmake --build "$root/build/nustar" --target nustar -j "$SLURM_CPUS_PER_TASK" > "$root/logs/nustar-build.log" 2>&1
make -C "$root/source/STAR/source" -j "$SLURM_CPUS_PER_TASK" STAR CXX=g++ > "$root/logs/star-build.log" 2>&1
"$root/build/nustar/nustar" --version > "$root/logs/nustar-version.txt"
"$root/source/STAR/source/STAR" --version > "$root/logs/star-version.txt"
gcc --version > "$root/logs/compiler.txt"
samtools --version > "$root/logs/samtools-version.txt"
(module list) > "$root/logs/modules.txt" 2>&1 || true
cp "$root/build/nustar/CMakeCache.txt" "$root/logs/CMakeCache.txt"
