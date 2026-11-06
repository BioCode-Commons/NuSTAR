# Benchmark template

Current implementation reviewed for documentation: `bdf6bfa`. **Measurements are pending.** Empty fields below are placeholders, not zeroes or passing checks. No benchmark was run as part of this documentation update.

## Inputs and machine

| Field | Value |
| --- | --- |
| Dataset name, public accession, and download link | To be supplied |
| Read layout, lengths, and number of reads/pairs | To be supplied |
| Full library or subset; selection method and seed | To be supplied |
| Reference FASTA and GTF source/version | To be supplied |
| CPU model, sockets, cores, and allocated threads | To be supplied |
| RAM, operating system, and storage | To be supplied |
| Compiler, build flags, STAR and NuSTAR revisions | To be supplied |
| Index-build and alignment commands | To be supplied |
| Cache conditions, run order, and repetitions | To be supplied |

## Runtime and memory

Measure both engines with the same inputs, explicit settings, thread count, and output mode. Report repeated wall times and variation, plus peak resident memory. Record index construction separately from alignment; include junction insertion in two-pass timing. Separate node and dataset results instead of pooling them.

| Workload | Dataset / threads / repeats | STAR time / RSS | NuSTAR time / RSS | Output comparison |
| --- | --- | --- | --- | --- |
| Index construction | Pending | Pending | Pending | Pending |
| Single-end, one pass | Pending | Pending | Pending | Pending |
| Paired-end, overlap disabled | Pending | Pending | Pending | Pending |
| Paired-end, overlap enabled | Pending | Pending | Pending | Pending |
| Two-pass | Pending | Pending | Pending | Pending |
| Coordinate-sorted BAM | Pending | Pending | Pending | Pending |
| Solo counts and corrected BAM tags | Pending | Pending | Pending | Pending |

## Are the outputs identical?

**Not established for the updated implementation by this document.** Added regression tests are available in the source, but their existence is not a result of a fresh benchmark run.

For each measured workload, report decoded alignment records and tag values, the treatment of record order and headers, splice-junction tables, and non-time mapping statistics. Compare both passes and inserted junction tables for two-pass runs. Compare matrices, cell filtering, and corrected tags for Solo. Preserve and describe any differences; equal mapping totals alone are insufficient. Agreement with STAR is a compatibility check, not biological ground-truth accuracy.

## Historical records

- [Earlier STAR comparison](STAR_COMPARISON.md): measurements from older builds, including figures whose input provenance is unresolved.
- [November 4 Oscar campaign](benchmarks/oscar-20261104/REPORT.md): partial RNA-seq results and a genuine paired-end mismatch from NuSTAR `90de97d`.
- [Discrepant read fixtures](benchmarks/oscar-20261104/read-discrepancy/discrepant-24-pairs/README.md): debugging evidence from that older run.

These records are retained as evidence. They have not been replaced with fabricated current results, and their conclusions should not be transferred to the updated binaries without rerunning the comparisons.
