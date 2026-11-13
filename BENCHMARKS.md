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
| Index construction | 1094.262 s / 70.013 GiB | 375.877 s / 47.798 GiB | passed |
| Single-end alignment | 658.994 s / 33.481 GiB | 503.396 s / 29.918 GiB | passed |
| Paired-end alignment | 1610.968 s / 33.488 GiB | 1217.181 s / 29.960 GiB | passed |
| Paired overlap refinement | 1904.047 s / 34.852 GiB | 1488.024 s / 30.058 GiB | passed |
| Two-pass mapping | 1644.326 s / 49.496 GiB | 1562.513 s / 40.732 GiB | passed |
| Solo counting | 121.487 s / 35.938 GiB | 84.823 s / 32.100 GiB | passed |

Output checks will report whether alignment records and tags, junctions, mapping statistics, and Solo matrices match. File-byte equality will be reported separately from decoded-record equality.

