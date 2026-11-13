# NuSTAR and STAR: runtime and correctness

**Latest full-human result:** after the [junction-padding fix](JUNCTION_PADDING_FIX.md), NuSTAR completes both passes on all 535.7 million reads. Both passes' junction tables and non-time mapping statistics match STAR. Whole-process time is 1906.23 seconds for NuSTAR versus the saved STAR baseline of 1598.68 seconds: NuSTAR is 19.2% slower. Peak RSS is 40.73 versus 49.49 GiB. The [earlier full-data report](DATA_FEATURE_VALIDATION.md) retains the failure that led to this fix and the real PBMC Solo baseline. These are single runs, and STAR's full-human result was reused rather than rerun alongside the repair.

The earlier saved human-reference runs show NuSTAR finishing the smaller tested two-pass, Solo, and BAM workloads sooner than STAR 2.7.11b. The compared outputs match in the cases listed below. This does not establish that NuSTAR is faster for every feature, or that every STAR option behaves identically.

The tables below summarize the earlier measurements; the linked full-data report contains the newer runs. The measurements cover several NuSTAR builds from the feature migration, Solo refactor, and BAM implementation. They are not all measurements of the current executable. Each run directory records its commands and binary hashes.

## What “correctness” means here

There are three different questions:

- **Agreement with STAR:** Do both programs produce the same alignments, junctions, counts, and tags? Most of the existing tests answer this question.
- **Mapping rate:** What fraction of reads maps uniquely or to multiple loci? This is available in `Log.final.out`, but mapping a read does not prove that its location is correct.
- **Accuracy against known truth:** Are the reported locations, genes, and cells actually correct? We do not yet have a systematic truth-based accuracy benchmark for NuSTAR.

“Exact match” below means agreement with STAR within the stated comparison. It does not mean 100% biological accuracy. Test pass rates are also not read-level accuracy estimates.

## Timed workloads

The FASTQ workloads use one million single-end reads from the human SRR21781765 sample against the complete GRCh38 reference. They use the local `stratified-small.fastq` sample, not the full sequencing run. Solo adds synthetic cell barcodes and UMIs to those bulk RNA reads. This is useful for testing the implementation, but it is not a real single-cell experiment.

Both programs ran on the same machine with 1 or 32 requested threads. The benchmark warmed the indexes and ran the programs serially under a benchmark lock. Each table entry is one process run. There are no error bars or repeated-run averages here.

Times are wall-clock seconds for the complete process: loading the index, reading input, mapping, counting where enabled, and writing output. Initial index construction is excluded. Two-pass index extension is included. The BAM runs use compression level 1 and a 2 GiB sort budget. All timed Solo runs use `--soloCellFilter None`.

Speedup is `STAR time / NuSTAR time`; a value above 1 means NuSTAR finished sooner.

| Workload | Threads | STAR (s) | NuSTAR (s) | Speedup | Evidence group |
| --- | ---: | ---: | ---: | ---: | --- |
| Two-pass, no alignment file | 1 | 139.27 | 126.28 | 1.10x | A |
| Two-pass, no alignment file | 32 | 89.76 | 78.84 | 1.14x | A |
| Solo gene counting, no alignment file | 1 | 44.75 | 31.54 | 1.42x | B |
| Solo gene counting, no alignment file | 32 | 20.74 | 7.79 | 2.66x | B |
| Solo + coordinate-sorted BAM with corrected CB/UB | 1 | 50.09 | 40.08 | 1.25x | C |
| Solo + coordinate-sorted BAM with corrected CB/UB | 32 | 24.30 | 9.03 | 2.69x | C |
| Solo + unsorted BAM | 1 | 48.57 | 36.66 | 1.32x | D |
| Solo + unsorted BAM | 32 | 22.58 | 7.40 | 3.05x | D |
| BAM input, Solo counting, no alignment file | 1 | 53.97 | 41.06 | 1.31x | E |
| BAM input, Solo counting, no alignment file | 32 | 24.05 | 9.02 | 2.67x | E |

The last workload remaps the 1,215,901 alignment records in a previous output BAM, so it processes more input records than the FASTQ workloads. NuSTAR reads BAM directly. STAR receives SAM text from the benchmark's BAM-to-SAM converter through `readFilesCommand`. That is a comparison of the usable input paths, not just the two decoders.

These are whole-run comparisons. The sorted-BAM rows do not show that NuSTAR's sorting algorithm alone is 2.69 times faster. Likewise, the Solo rows include reference loading and mapping; they do not measure UMI deduplication in isolation.

### What the times suggest

Two-pass has a much smaller advantage than Solo in these runs: about 9–12% less elapsed time. Going from 1 to 32 threads gives about 1.60x speedup for NuSTAR's two-pass run and 1.55x for STAR. There is still substantial work that does not benefit from the extra mapping threads. The saved whole-process times do not isolate which step accounts for it.

For the sorted-BAM single-thread run, user CPU time is nearly the same: 34.42 seconds for STAR and 34.52 seconds for NuSTAR. System CPU time is 15.61 versus 5.53 seconds. That supports looking at file access, loading, and output handling when explaining the wall-time gap; it does not support claiming that NuSTAR's alignment algorithm became faster by the full measured speedup. CPU times and wall time are different measurements, especially with multiple threads.

The no-BAM Solo result in group B is from the refactored Solo build. The older `20.61 / 7.67` second result at 32 threads belongs to the earlier feature-migration build; it is not an additional measurement of group B.

## Correctness by feature

| Feature | Comparison performed | Result and limit |
| --- | --- | --- |
| Ordinary single-end and paired-end mapping | CLI regression against STAR; records, junction tables, defaults, and one/two-thread execution | Passed. This is a small generated fixture, not a large paired-end accuracy study. |
| Two-pass mapping | 12 generated cases with 227 reads or pairs per case; annotated/unannotated indexes, SE/PE, 1/4 threads, all reads or a 40-read first pass | SAM record multisets, junctions, non-time statistics, and inserted junction tables match. The human timed runs also match the junction tables and mapping statistics; they used no SAM/BAM output, so those runs do not provide a full alignment-record comparison. |
| Solo core counting | 22 cases on 2,400 generated RNA reads, at 1/4 threads | All passed. Includes gene/full-gene counting, UMI modes, multimappers, strand settings, velocity, cell statistics, and filtering combinations. This covers selected combinations, not every possible parameter combination. |
| Solo input adapters and standalone filtering | 27 additional cases | All passed. Includes paired and inline barcodes, complex layouts, no whitelist, SmartSeq, read limits, two-pass integration, SJ counting, and standalone TopCells/CellRanger2.2/EmptyDrops_CR. |
| Human Solo gene counting | Groups B–E compare output file sets, raw matrices, barcode/feature statistics, gene summaries, junctions, and non-time mapping statistics | Compared files match exactly at 1 and 32 threads. Because barcodes are synthetic and cell filtering is disabled, these runs do not measure real cell recovery or EmptyDrops accuracy. |
| Coordinate-sorted BAM | Decode every human output record and compare in order | All 1,215,901 records match at both 1 and 32 threads, including tags and coordinate-tie ordering. This is decoded-record equality, not equality of compressed BAM bytes or program headers. |
| Unsorted BAM | Decode every human output record and compare multisets | All 1,215,901 records match at both thread counts, including multiplicity. Record order is not required to match in this comparison. |
| SAM/BAM input and BAM output modes | 37 integration cases, plus a focused paired-flag case | All passed. Cases cover SE/PE, raw/corrected tags, compression, two-pass output, tag retention, and a rejected mismatched-mate input. The separate flag case matches 10,336 decoded records. |
| Index construction and individual mapping stages | Unit tests and dedicated STAR comparison tools exist | No current NuSTAR-specific human timing and accuracy summary is established here. Older results copied from KrinonRefactored are not treated as new NuSTAR measurements. |
| Features not yet implemented | No comparison possible | No performance or correctness claim. See the feature table in [readme.md](../readme.md). |

For the human sorted-BAM comparison, the decoded-record agreement is 100% of the 1,215,901 checked records. For unsorted BAM, the complete record multisets agree. Neither number estimates the fraction of reads aligned to their true biological origin.

### Mapping rates in the human runs

These values come from group A's `Log.final.out` files, including its accompanying one-pass Solo runs. STAR and NuSTAR report the same percentages at both thread counts.

| Mode | STAR uniquely mapped | NuSTAR uniquely mapped | STAR mapped to multiple loci | NuSTAR mapped to multiple loci |
| --- | ---: | ---: | ---: | ---: |
| One-pass + Solo | 88.22% | 88.22% | 8.71% | 8.71% |
| Two-pass | 88.19% | 88.19% | 8.75% | 8.75% |

The shift between one-pass and two-pass is not evidence of an accuracy gain or loss. Adding junctions can change which candidate alignments are available. A truth set would be needed to decide whether those changed assignments are better.

## Comparisons still missing

| Area | What still needs to be measured |
| --- | --- |
| Index construction | A fresh NuSTAR/STAR build-time comparison from the same FASTA/GTF, with exact index-domain validation. |
| Seed search, windows, stitching, selection | Stage timings from matching current builds, alongside intermediate-output comparisons. Whole-run timing cannot identify which stage is faster. |
| Ordinary mapping without Solo | The newer full human first pass matches junctions and statistics, but complete per-read comparison and current-build SE/PE scaling measurements are still missing. |
| Solo substeps | Separate timings for barcode correction, annotation, UMI deduplication, multimapper redistribution, and cell filtering. |
| Real single-cell data | The newer full PBMC v3 test matches raw/filtered matrices and summaries. Real-data corrected BAM tags, more chemistries, and truth-based accuracy tests remain to be checked. |
| Less common modes | More coverage for Transcript3p, parameter combinations, long CIGARs, and unusual inputs. The full human two-pass boundary failure also needs to be fixed and retested. |
| Reproducibility | Repeated runs of the same current binaries, with hardware, cache state, and disk conditions recorded. |

## Where the numbers come from

These links point to local generated artifacts under `NuSTAR/build/`. They may be absent in a fresh checkout because `build/` is ignored by Git. The tables above retain the main results; keep the artifacts when reviewing or reproducing a claim.

| Group | Build/workload | Timing and commands | Output comparisons |
| --- | --- | --- | --- |
| A | Verified two-pass migration run | [results.json](../build/feature-migration/performance-verified/results.json), [metadata.json](../build/feature-migration/performance-verified/metadata.json) | [validation.json](../build/feature-migration/performance-verified/validation.json) |
| B | Solo after module refactoring, no alignment file | [results.json](../build/solo-refactor/performance/results.json), [metadata.json](../build/solo-refactor/performance/metadata.json) | [validation.json](../build/solo-refactor/performance/validation.json) |
| C | Final sorted-BAM run | [results.json](../build/bam-performance-final-sorted/results.json), [metadata.json](../build/bam-performance-final-sorted/metadata.json) | [feature outputs](../build/bam-performance-final-sorted/validation.json), [records, 1 thread](../build/bam-performance-final-sorted/bam-record-validation-1.json), [records, 32 threads](../build/bam-performance-final-sorted/bam-record-validation-32.json) |
| D | Final unsorted-BAM run | [results.json](../build/bam-performance-final-unsorted/results.json), [metadata.json](../build/bam-performance-final-unsorted/metadata.json) | [feature outputs](../build/bam-performance-final-unsorted/validation.json), [records](../build/bam-performance-final-unsorted/bam-record-validation.json) |
| E | Final BAM-input run | [results.json](../build/bam-performance-final-input/results.json), [metadata.json](../build/bam-performance-final-input/metadata.json) | [validation.json](../build/bam-performance-final-input/validation.json) |

Small-case reports: [two-pass](../build/feature-migration/two-pass/validation.json), [Solo core](../build/bam-solo-regression/validation.json), [Solo adapters](../build/bam-solo-regression/adapter-validation.json), [BAM](../build/bam-validation-final/validation.json), and [paired flags](../build/bam-validation-final/pair-flags-validation.txt).

The test entry points are in [VALIDATION.md](VALIDATION.md). [benchmark_features.py](../tools/benchmark/benchmark_features.py) records complete commands, wall time, CPU time, peak RSS, and output comparisons. Use its `--help` output when preparing a new run.
