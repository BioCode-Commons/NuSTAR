# NuSTAR / STAR benchmark on Oscar

The campaign is prepared locally. It has **not** run on Oscar. No result from the local integrity checks is a performance measurement.

## Frozen inputs

`sources.json` records the actual public URLs, byte counts, published ENCODE MD5 values, PBMC archive member offsets and ETag, and a pinned whitelist SHA256.

- Bulk: ENCSR000AEM replicate 1, ENCFF001RED/ENCFF001RDZ, PE101 K562 polyA-plus RNA-seq.
- Single-cell: 5k_pbmc_v3, all four lanes, with R2 as RNA and R1 as barcode/UMI. Only R1/R2 gzip members are downloaded through verified HTTP byte ranges. I1 files and a duplicate full tar archive are not stored. GEX compressed input totals 26,625,946,771 bytes.
- Reference: 10x `refdata-gex-GRCh38-2020-A`; extract its FASTA, GTF, and `reference.json` only. Discard the packaged STAR index and rebuild both programs' indexes with `sjdbOverhang=100`, the longest RNA read length in this campaign minus one. Both datasets use the same reference inputs. This compares STAR with NuSTAR; it is not a Cell Ranger accuracy comparison.
- Whitelist: 3M-february-2018, pinned to the Teichlab source revision in the manifest. The verified gzip contains 6,794,880 16-base barcodes. It is not a whitelist reconstructed from this sample's observed barcodes.
- Engines: STAR 2.7.11b (`b1edc12`) and NuSTAR `90de97d`. Those are immutable source revisions; new benchmark results must name the revisions actually built.

The old SRR21781765 measurements are excluded from this campaign: the cited accession is WXS and its relationship to the historical FASTQ subset is unresolved.

## Matrix

| Stage | Inputs | Cases | Threads | Repetitions per engine |
| --- | --- | --- | --- | ---: |
| Pilot | Uniform one-million-pair subsets | Bulk SE, PE, PE two-pass, PE sorted BAM; Solo gene counts and sorted BAM | 1, 32 | 1 |
| Main | Uniform ten-million-pair subsets | Same six cases | 1, 8, 32 | 3 |
| Full bulk | Entire selected K562 replicate | Bulk SE, PE, PE two-pass, PE sorted BAM | 32 | 3 |
| Full PBMC | Entire GEX library | Solo Gene + EmptyDrops_CR, expected cells 5000, no alignment file | 32 | 3 |

Subsets use seed 20261104, sampling without replacement after validating/counting the complete compressed source. The pilot is nested within the main subset. Both mates and all lanes stay synchronized, and input order is retained. Sampling is outside timing. Subset Solo uses `soloCellFilter None`; full-library filtering is reported separately. Full-library PBMC BAM tags are **not** covered; corrected tags are compared on the ten-million-pair main input.

Cases and arguments are in `suite.json`. Within every pair, both engines receive the same inputs, thread count, and biological options; only executable, native index, and output prefix differ. Alternate which engine runs first. Index pages are read before every timed process. Input-file cache state is not forced; do not label these cold-cache measurements.

## Slurm and storage

`submit.sh` submits three dependent `batch` jobs, each requesting **24 hours, one node, 32 CPUs, no SMT, and 192 GiB RAM**. Main includes preparation, builds, indexes, pilot, and main runs. Full bulk depends on main passing; full PBMC depends on full bulk passing. A failure cancels invalid downstream dependencies. Resources must first pass the actual account/partition checks; command-line submission options can be adjusted to the verified limits.

All downloads, source/build trees, native indexes, uncompressed inputs, sorting runs, logs, and results belong under one scratch campaign directory. The harness itself is also transferred to scratch. `TMPDIR` is inside the campaign. Home holds no datasets or indexes.

Set `--disk-budget-gib` from the real remaining **user quota**, not filesystem free space alone. An illustrative 350 GiB campaign budget is not assumed to be available. The scripts check allocated campaign bytes and retain a 32 GiB filesystem free-space reserve. Download/preparation checks run during writes; an external-writer watchdog samples every two seconds and terminates a process that exceeds the selected campaign budget. This is a monitored budget, not a filesystem quota reservation.

- Do not materialize the vendor reference's prebuilt index.
- Prepare full uncompressed libraries one at a time. Their exact required FASTQ bytes are measured in the initial scan and checked before expansion.
- Keep at most one STAR/NuSTAR output pair awaiting validation. Decode sorted BAM directly through two `samtools view -h` streams, without writing expanded SAM copies.
- Prune only this pair's generated BAMs **after** exact record validation; retain the record count and ordered-record SHA256. Preserve failed output pairs for diagnosis.
- Release full-library FASTQs only when their suite has completed and their hashes still match. Keep their manifest, source archive receipts, and subset inputs.
- Preserve logs, commands, binary/input/index hashes, raw measurements, Solo matrices, junctions, first differences, and summaries. Transfer the compact evidence bundle locally before scratch retention expires.

## One login, then launch

Do not trigger Duo until local preparation is finished. For the actual run:

1. Authenticate once with a dedicated SSH master, for example `ControlMaster=auto`, `ControlPersist=7d`, and a task-specific socket. No password or Duo code is stored in the harness. Scheduler jobs continue independently of the SSH connection.
2. Through that connection, inspect `checkquota`, `sinfo -p batch`, `scontrol show partition batch`, account/QOS limits, and available GCC/CMake/Python/samtools modules. Check that scratch resolves to its intended mount and that each 24-hour resource request is valid. These account-specific facts cannot be established offline.
3. Load the selected toolchain, capture module versions, transfer this harness to scratch, and run `sbatch --test-only` with the final resources. Choose a campaign budget within the remaining user quota. Do not silently change the dataset or test matrix to fit missing resources.
4. Run `bash submit.sh CAMPAIGN_ROOT DISK_BUDGET_GIB`, supplying the actual scratch path and chosen budget. The root must be new and below the resolved `~/scratch` directory. Submission IDs are saved immediately so partial submission cannot lose track of an accepted job.
5. Monitor with `bash status.sh CAMPAIGN_ROOT`. Use `sacct` for completed jobs and peak memory; inspect `results/*/status.json` and `logs/*` for progress/failures. Arrange a heartbeat against these concrete job IDs after submission, reporting meaningful transitions or action needed rather than every unchanged poll.

No job is automatically resubmitted after a failed output comparison. Inspect the retained evidence, resolve the cause, and resume deliberately through the same authenticated connection while it remains valid.

## Checks and retained evidence

For each timed process, retain wall time, user/system CPU time, Linux `ru_maxrss` in KiB, exact command, return code, and storage observations. Metadata captures CPU model/topology and affinity, RAM, allocated node, Slurm parameters/sharing, compiler, samtools, build options, and source/data identities. Warm-up, preparation, checksums, and validation are outside alignment timings. Index construction is measured separately.

Validate exact junction files, both passes' junction/insertion tables, non-time mapping statistics, Solo file sets/matrix contents, and ordered decoded BAM records including tags. Compare `@SQ` headers; exclude other headers and compressed BAM bytes explicitly. Require the input count to match the prepared dataset's count. Retain a first mismatch and fail the case rather than summarize it as compatible.

`summary.json` reports median/min/max wall time and STAR/NuSTAR ratios, with every raw run retained. Different phase jobs may use different hardware; compare engines and thread scaling within their allocation, and report the hardware for each phase rather than pooling incompatible nodes.

Local checks:

```sh
python3 -m unittest discover -s tools/benchmark/oscar -p 'test_*.py' -v
bash -n tools/benchmark/oscar/campaign.sbatch tools/benchmark/oscar/setup.sh tools/benchmark/oscar/submit.sh tools/benchmark/oscar/status.sh
```

These validate file/mate integrity, mismatches, capacity checks, and argument symmetry. Actual compilation, Slurm allocation, aligner behavior, and the pilot are verified on Oscar after authentication. Local checks do not simulate either aligner.
