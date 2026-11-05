# Oscar benchmark: partial results and paired-end mismatch

**Status: the campaign stopped on an output mismatch. It has not passed the full benchmark.** Jobs 8361933 and 8361934 were cancelled because their successful dependency was not met. All failed-pair outputs remain on Oscar. This report records the completed measurements; it does not treat the paired-end timing as an equivalent-output speedup.

## Inputs and execution

- Bulk: ENCODE K562 polyA-plus experiment ENCSR000AEM, files ENCFF001RED/ENCFF001RDZ, paired 101 bp. The complete source contains 92,172,367 read pairs. Formal subset: 10,000,000 pairs sampled uniformly without replacement, seed 20261104; the single-end case uses R1 alone.
- Single-cell: 10x 5k_pbmc_v3, four lanes, 383,941,607 read pairs. Its formal subset was prepared, but no Solo measurement has run.
- Reference: GRCh38-2020-A FASTA/GTF, separate native indexes, sjdbOverhang 100. See the [frozen sources](../../../tools/benchmark/oscar/sources.json) and [plan](../../OSCAR_BENCHMARK_PLAN.md).
- NuSTAR source: `90de97ddce76b1f115f47e37d15add62784bafed`. STAR 2.7.11b source: `b1edc1208d91a53bf40ebae8669f71d50b994851`.
- Oscar batch job 8361932, node1848: Intel Xeon Platinum 8268 @ 2.90 GHz; 32 allocated physical CPUs, 192 GiB requested RAM, 24-hour limit. The node is shared. GPFS jobtmp storage; campaign disk budget 350 GiB.
- Index pages are warmed before each timed process. Input cache state is not forced. Timing covers the process, including native-index loading; preparation, warming, and validation are outside the interval. Linux peak RSS is recorded in KiB.
- Pilot was skipped at the user's request. The [execution policy](evidence/execution-policy.json) and [change record](evidence/logs/harness-revisions/skip-pilot/change.json) retain that decision.

## Completed single-end measurements

All nine STAR/NuSTAR pairs passed exact junction-file comparison and equality of non-time mapping statistics. No SAM/BAM was requested in this case, so these checks do not prove read-level record or tag identity. Each row reports the median of three runs on the same node.

| Threads | STAR seconds | NuSTAR seconds | STAR / NuSTAR | STAR peak RSS, GiB | NuSTAR peak RSS, GiB |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 354.33 | 319.38 | 1.11× | 27.78 | 27.84 |
| 8 | 59.31 | 49.56 | 1.20× | 28.47 | 27.88 |
| 32 | 31.17 | 21.49 | 1.45× | 30.78 | 28.05 |

Individual commands, measurements, and comparison results are in `evidence/results/main-job-8361932/bulk-se-*/validated.json`. Hardware, input hashes, and index manifests are captured in [metadata.json](evidence/results/main-job-8361932/metadata.json).

## Paired-end mismatch

The first paired-end run used the same ten-million-pair subset, one thread, one pass, and no SAM/BAM output. Both processes exited successfully, but their outputs differ. STAR took 870.856 seconds with peak RSS 27.777 GiB; NuSTAR took 764.465 seconds with peak RSS 27.960 GiB. These timings describe different outputs and are not a passing compatibility result.

Both junction files have 164,877 rows and the same genomic coordinates. Two rows have different support counts:

| Junction, 1-based inclusive | STAR unique support | NuSTAR unique support | STAR multimapping support | NuSTAR multimapping support |
| --- | ---: | ---: | ---: | ---: |
| chr12:6668116–6669030 | 39 | 26 | 0 | 5 |
| chr12:6669190–6670759 | 41 | 25 | 0 | 5 |

Their strand, motif, annotation status, and maximum overhang fields match. Mapping statistics also differ: STAR reports 8,213,543 uniquely mapped fragments and 1,305,608 multimapping fragments; NuSTAR reports 8,213,538 and 1,305,613 respectively. Total splice counts are 6,049,989 versus 6,049,960. Both processed exactly 10,000,000 input pairs.

A diagnostic replay has now identified the affected reads. Job 8381167 ran both unchanged engines with one thread on node1827 and enabled SAM output. Each engine exactly reproduced its original junction table and non-time mapping statistics. This diagnostic run is separate from the performance measurements on node1848.

There are **24 read pairs with different SAM records**: five change from `NH:i:1` in STAR to `NH:i:2` in NuSTAR; the other 19 retain `NH:i:1` but have different alignment records. All 24 belong to the 56 pairs that span either affected junction in at least one engine. The handoff includes those 56 pairs, including 32 additional junction-support pairs whose records match. These observations do not establish the implementation bug's cause.

Replaying only the 56 extracted pairs reproduced the same 24 differing read IDs. For each engine, all SAM records for all 56 pairs were also identical to that engine's corresponding records in the full ten-million-pair replay, after sorting optional tags and record order. The original failed outputs and full diagnostic outputs remain on Oscar.

Download the [developer handoff archive](read-discrepancy/paired-read-handoff.tar.gz), or read its [instructions](read-discrepancy/paired-read-handoff/README.md). It contains paired FASTQs, both full-replay and subset SAMs, exact differences, read IDs, executable hashes, reference/index provenance, and a reproduction script. The root cause is still unresolved; neither aligner was changed.

Evidence: [STAR junctions](evidence/results/main-job-8361932/bulk-pe-1-r1/star/output.SJ.out.tab), [NuSTAR junctions](evidence/results/main-job-8361932/bulk-pe-1-r1/nustar/output.SJ.out.tab), [STAR statistics](evidence/results/main-job-8361932/bulk-pe-1-r1/star/output.Log.final.out), [NuSTAR statistics](evidence/results/main-job-8361932/bulk-pe-1-r1/nustar/output.Log.final.out), and [failed suite status](evidence/results/main-job-8361932/status.json).

## Index construction

One measured build per program, 32 threads: STAR 1,385.415 seconds / 34.839 GiB peak RSS; NuSTAR 802.266 seconds / 85.271 GiB. NuSTAR was 1.73× faster in this build and used about 2.45× the memory. Index construction is a separate workload from alignment; the indexes use different native formats. Raw records are in [STAR timing](evidence/logs/index-star/timing.json) and [NuSTAR timing](evidence/logs/index-nustar/timing.json).

## Remaining coverage and limits

Paired-end repetitions and multithreaded runs, two-pass, sorted BAM, Solo, full bulk, and full PBMC have not completed. There is no current two-pass performance or compatibility conclusion. Historical SRR21781765 figures are excluded: that accession is WXS, and its relationship to the old local subset remains unresolved. These new measurements use documented RNA-seq inputs. Output agreement with STAR is not a biological accuracy benchmark.
