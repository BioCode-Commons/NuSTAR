# Oscar benchmark plan

Status: preparation only. No Oscar job has been submitted. The [campaign harness](../tools/benchmark/oscar/README.md), [frozen sources](../tools/benchmark/oscar/sources.json), and [case matrix](../tools/benchmark/oscar/suite.json) are now available. Complete local checks before requesting Brown/Duo authentication.

## Runtime and storage

Request the `batch` partition with an initial **24-hour walltime**, subject to its actual account and partition limits. This is a scheduler walltime request, not an SSH timeout. Once `sbatch` accepts a job, it does not require the SSH connection to remain open.

Use an isolated directory under the account's actual Oscar scratch location for downloads, extracted FASTQs, both programs' indexes, build files, outputs, and temporary sorting runs. Discover that location and its quota after login; do not assume a path or write large files into home.

Before downloading or extracting a dataset:

1. Check the scratch quota, used space, filesystem free space, and existing shared inputs.
2. Record the files already available and their checksums; reuse matching data rather than duplicate it.
3. Include compressed inputs, extracted inputs, both indexes, output files, validation intermediates, and sort spill files in the capacity assessment.
4. Process datasets and large output workloads sequentially. Complete the record comparison before deleting that run's generated BAMs, decoded files, or temporary sort data. Keep manifests, exact commands, raw timing JSON, logs, validation summaries, and discrepancy evidence.
5. Delete only files created by this benchmark. Preserve failed outputs until the mismatch has been diagnosed.

The exact CPU, memory, and storage allocation will be recorded from the assigned compute node. Build, indexing, validation, and metadata collection stay outside alignment timing intervals. Index construction can be measured separately with its own memory request.

## Dataset provenance correction

The historical report cites **SRR21781765** as RNA-seq. ENA metadata checked on November 4, 2026 reports:

- Organism: `Homo sapiens`
- Library strategy: **WXS**
- Layout: **PAIRED**
- Reported read count: `535736296`
- FASTQ download bytes: `18762515622;20829356726`

Source: [ENA run record](https://www.ebi.ac.uk/ena/browser/view/SRR21781765), queried through the [ENA portal API](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR21781765&result=read_run&fields=run_accession,scientific_name,library_strategy,library_layout,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp&format=json).

This accession cannot establish the identity of the historical `stratified-small.fastq` without the original file or checksum. Keep the old numbers labeled as historical measurements with unresolved provenance. Do not use this accession as the main biological RNA-seq benchmark.

## Main bulk RNA-seq dataset

Use **ENCODE ENCSR000AEM, K562 polyA-plus RNA-seq**, biological replicate 1, with the mutually paired released FASTQs below:

| Mate | File accession | Read length | Compressed bytes | Published MD5 |
| --- | --- | ---: | ---: | --- |
| R1 | ENCFF001RED | 101 bp | 7,668,466,665 | `0e88b38c8d48806ee98929a5bad06fc4` |
| R2 | ENCFF001RDZ | 101 bp | 7,848,217,965 | `3a25713488522b0ea8a4d9635b5afa41` |

Primary sources: [experiment](https://www.encodeproject.org/experiments/ENCSR000AEM/), [R1 metadata](https://www.encodeproject.org/files/ENCFF001RED/?format=json), and [R2 metadata](https://www.encodeproject.org/files/ENCFF001RDZ/?format=json). Metadata was checked on November 4, 2026.

Run a one-million-pair pilot to check the environment and comparison machinery. Main measurements use a reproducible ten-million-pair subset, then the full selected replicate if the confirmed scratch quota and allocation support it. Keep mates synchronized and save the extraction method, source checksums, selected read count, and subset checksums. Do not call a newly sampled subset identical to a historical subset without matching its identity.

Cover paired-end one-pass alignment, `Basic` two-pass alignment, and coordinate-sorted BAM. Single-end results may use R1 alone, explicitly labeled as a single-end derivation of the paired library.

## Real single-cell dataset

Use **10x Genomics 5k PBMCs from a healthy donor, v3 chemistry**, dataset **`5k_pbmc_v3`**, released May 29, 2019 and analyzed by Cell Ranger 3.0.2.

The official page lists R1 as 28 bp (16-base cell barcode plus 12-base UMI) and R2 as 91 bp transcript sequence. Record the exact GEX FASTQ files, chemistry-specific whitelist, checksums, and any lane-merging procedure before running either aligner.

Source: [10x dataset page](https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-1-standard-3-0-2).

Use a synchronized subset for the pipeline pilot, then the complete selected library for cell filtering and full raw/filtered matrix comparison if scratch capacity permits. Label subset cell-calling results separately, since downsampling changes the observed library. Cover Solo gene counting without alignment output and corrected sorted-BAM tags on the ten-million-pair subset. Run full-library gene counting and cell filtering without a full-library BAM to control scratch use. Full-library corrected BAM tags are not covered by this campaign. Synthetic barcodes may supplement implementation checks; they are not a substitute for this real single-cell dataset.

## Reference and comparison

Both programs build their own indexes from the FASTA and GTF in the pinned 10x GRCh38-2020-A reference package. Extract only FASTA, GTF, and reference metadata, then discard its bundled STAR index. Record input hashes and use `sjdbOverhang=100` for both programs and both datasets. Reference sources are frozen in the campaign manifest; do not silently substitute another release.

Use STAR **2.7.11b** against the recorded NuSTAR commit. Run both executables serially on the same allocated node with matching thread counts and output settings. Use **1, 8, and 32 threads**, at least **three repetitions**, and alternate which program runs first. Record warm-cache behavior explicitly; do not claim a cold-cache test without controlling and measuring it.

Save each process's wall time, user/system CPU time, peak RSS, exit code, exact command, and binary hash. Record node CPU model, core topology, process CPU affinity, RAM, Slurm allocation, OS, compiler/build options, and storage mounts. Record any node sharing; do not describe a shared allocation as isolated.

Compare junction files, non-time mapping statistics, Solo file sets and matrix contents, and decoded BAM records including tags and coordinate-tie order. Exclude program-specific headers and compressed-file bytes only when that exclusion is explicitly recorded. Keep first differences and nonzero exit codes visible; a timing result without a passing output comparison does not establish compatibility.

Report the individual runs, median, range, and STAR/NuSTAR ratio. Separate input preparation and indexing measurements from end-to-end alignment. Do not fill the README hardware placeholders using hardware from the old runs or the login node.

## Execution package

See the [campaign guide](../tools/benchmark/oscar/README.md) for the exact matrix, three dependent 24-hour jobs, capacity enforcement, cleanup rules, source revisions, and the one-login launch/monitor procedure. Local checks are not a substitute for the Oscar pilot.
