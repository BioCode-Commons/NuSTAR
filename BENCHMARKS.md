# Benchmarks

[← NuSTAR README](readme.md#benchmarks)

| Test setup | Details |
| :--- | :--- |
| Datasets | SRR21781765: 535,736,296 single-end reads (R1) / paired-end read pairs (R1 + R2), full dataset, average 148 bases per mate; 10x Genomics PBMC 1k v3: 66,601,887 RNA/barcode pairs, full dataset (R2 RNA, R1 barcode/UMI). |
| Reference | GRCh38 primary assembly (`GRCh38.primary_assembly.genome.fa`); GENCODE v50 (`gencode.v50.primary_assembly.annotation.gtf`). |
| Hardware |  AMD EPYC 9634 84-Core Processor  |
| Software and run settings | NuSTAR `; STAR 2.7.11b; 32 threads; one run per program/case; serial NuSTAR then STAR; no CPU pinning or cache eviction; `--sjdbOverhang 150`; overlap: `--peOverlapNbasesMin 20 --peOverlapMMp 0.01`; two-pass: `--twopassMode Basic --twopass1readsN -1 --outSAMorder Paired`; Solo: Gene, CB 16 / UMI 12, `CellRanger2.2 3000 0.99 10`, `--outSAMtype None`. |

Single-end alignment, paired-end alignment, and paired overlap refinement use the [20261113 profile results](build/full-star-profile-unified-20261113/summary.csv).

| Workload | STAR time / peak RSS | NuSTAR time / peak RSS | Output comparison |
| :--- | :--- | :--- | :--- |
| Index construction (libsais40) | 1095.416 s / 70.016 GiB | 374.413 s / 47.799 GiB | passed |
| Single-end alignment | 658.994 s / 33.481 GiB | 503.396 s / 29.918 GiB | passed |
| Paired-end alignment | 1610.968 s / 33.488 GiB | 1217.181 s / 29.960 GiB | passed |
| Paired overlap refinement | 1904.047 s / 34.852 GiB | 1488.024 s / 30.058 GiB | passed |
| Two-pass single-end mapping | 1644.442 s / 49.492 GiB | 1212.230 s / 40.929 GiB | passed |
| Two-pass paired-end mapping | 4241.699 s / 50.198 GiB | 2926.378 s / 41.421 GiB | passed |
| Coordinate-sorted BAM, single-end | 926.828 s / 100.798 GiB | 751.109 s / 75.024 GiB | passed |
| Coordinate-sorted BAM, paired-end | 3408.854 s / 97.904 GiB | 1774.513 s / 61.716 GiB | passed |
| Solo counting | 118.887 s / 35.939 GiB | 81.207 s / 31.780 GiB | passed |

Output checks will report whether alignment records and tags, junctions, mapping statistics, and Solo matrices match. File-byte equality will be reported separately from decoded-record equality.

