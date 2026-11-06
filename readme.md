<div align="center">

# NuSTAR

**A refactored implementation of STAR for bulk and single-cell RNA-seq.**

[Quick start](#quick-start) · [Benchmarks](#benchmarks) · [Implementation](#implementation-status) · [Roadmap](#roadmap) · [Usage](docs/USAGE.md) · [Contributing](CONTRIBUTING.md)

</div>

---

NuSTAR reimplements STAR's core mapping stages and integrates refactored STARsolo code. It supports single-end and paired-end RNA-seq alignment, one-pass and two-pass workflows, and single-cell counting. The code is organized into modules for indexing, mapping, input/output, and Solo analysis.

For supported workflows, NuSTAR is designed to preserve STAR's command-line interface and behavior. It does not yet implement every STAR feature. See [implementation status](#implementation-status) for the current scope and [AI use and human review](#ai-use-and-human-review) for how the code was developed.

> [!IMPORTANT]
> Build a **NuSTAR index** before aligning reads; STAR binary indexes cannot be loaded directly. FASTQ input must be **uncompressed**. Chimeric/fusion alignment is not implemented.

## Quick start

### 1. Build

Requirements: **CMake 3.20+**, a **C++20 compiler**, a **C99 compiler**, **zlib development files**, and **Python 3** for the default test-enabled build. OpenMP is optional; HTSlib and libsais are bundled.

```sh
git clone https://github.com/BioCode-Commons/NuSTAR.git
cd NuSTAR

cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 8
ctest --test-dir build/release --output-on-failure
```

The executable is `build/release/nustar`. Check it with `build/release/nustar --version` or `build/release/nustar --help`.

<details>
<summary><strong>Build configuration</strong></summary>

Release builds enable interprocedural optimization when supported by the compiler. Add `-DNUSTAR_ENABLE_RELEASE_IPO=OFF` to the configure command to disable it.

To build without tests, configure with `-DBUILD_TESTING=OFF`; Python is then not required by the test configuration. Skip the `ctest` command in that case.

The default CTest run does not enable every STAR differential suite. Additional comparisons require a STAR baseline and, for some suites, reference data or instrumented binaries. See [contribution guide](CONTRIBUTING.md#tests).

</details>

### 2. Build a genome index

Use your own reference FASTA and matching GTF. The filenames below illustrate a human reference; these data files are not bundled.

```sh
build/release/nustar --runMode genomeGenerate \
  --genomeDir index/human \
  --genomeFastaFiles reference/GRCh38.fa \
  --sjdbGTFfile reference/genes.gtf \
  --sjdbOverhang 100 \
  --runThreadN 8
```

Choose `--sjdbOverhang` for your reads; read length minus one is the usual starting point. Omit `--sjdbGTFfile` to build without annotation. **Solo gene counting requires the annotation tables produced from a GTF.**

NuSTAR stores its native index in `reference.bin` and `index.sa`. For a small test genome, reduce `--genomeSAindexNbases`; the human-reference default is not appropriate for every genome.

### 3. Align paired-end reads

Replace the FASTQ filenames with your uncompressed reads.

```sh
mkdir -p results

build/release/nustar --genomeDir index/human \
  --readFilesIn mate1.fastq mate2.fastq \
  --runThreadN 8 \
  --outFileNamePrefix results/sample. \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --limitBAMsortRAM 2147483648
```

This writes a coordinate-sorted BAM, a splice-junction table, and mapping statistics. `--runMode alignReads` is the default and can be omitted.

> [!NOTE]
> `--outFileNamePrefix` is a literal prefix: `results/sample.` produces `results/sample.Aligned.sortedByCoord.out.bam`. Use `results/sample/` to place outputs inside a sample directory. The BAM sort budget is not a total-process memory limit.


More examples—including single-end, two-pass, 10x, SmartSeq, and SAM/BAM input—and the full parameter reference are in the [usage guide](docs/USAGE.md).

## Benchmarks

**Pending a new benchmark of the updated implementation.** Runtime, peak memory, and output-equivalence results will be added after the current binaries are measured together. No current speedup or full compatibility claim is made here.

| Item | Status |
| :--- | :--- |
| Datasets, accessions, and subset selection | To be supplied |
| Reference genome and annotation versions | To be supplied |
| CPU, RAM, operating system, and storage | To be supplied |
| STAR/NuSTAR versions, build settings, and commands | To be supplied |
| Runtime and peak RSS, with repeated runs | To be measured |
| Alignment records, junctions, statistics, and Solo output comparisons | To be compared |

See the [benchmark template](docs/BENCHMARKS.md) for the workloads and fields to fill in. Earlier runs remain in the [historical reports](docs/BENCHMARKS.md#historical-records); they describe older builds and are not results for this revision.

## Implementation status

The two figures below compare STAR 2.7.11b with NuSTAR. The first shows the processing flow; the second expands implementation changes found in the current source. See the [source comparison](docs/IMPLEMENTATION_COMPARISON.md) for the compared revisions, code links, and scope.

### Updates in the current implementation

Documentation follows the implementation at [`bdf6bfa`](https://github.com/BioCode-Commons/NuSTAR/tree/bdf6bfa7d465c965edac091c6f13d292a81be1ef). The latest changes include:

- **Overlapping paired-end reads:** `--peOverlapNbasesMin` enables mate-overlap merge/remapping. NuSTAR maps the merged sequence with a separate single-read workspace, projects the result back to the mates, and recalculates paired scores. The option defaults to `0` (disabled); `--peOverlapMMp` defaults to `0.01`.
- **Junction strand handling:** annotation preparation now follows STAR's ordering and winner-selection rules when proposals share coordinates but differ in strand, motif, or priority. New regression fixtures cover these cases; the older human-data mismatch has not been rerun here.
- **Two-pass index extension:** the suffix-array merge is split into independent output-rank chunks. Each chunk finds its starting insertion by binary search and writes a separate packed range; OpenMP schedules the chunks when available.
- **Index-build memory lifetime:** temporary suffix positions are released while packing on Linux. Prefix tables are built after that temporary array is freed, and construction objects are released before index reload verification.
- **Clipped-read statistics:** input-base counts use clipped read lengths, following STAR's statistics convention.

These are implementation changes, not measured performance improvements. See the [detailed source comparison](docs/IMPLEMENTATION_COMPARISON.md#changes-in-the-updated-implementation).

### Processing flow

[![STAR and NuSTAR flows: index setup, read loops, mapping stages, ordered output, and finalization](docs/assets/implementation-flow.svg)](docs/assets/implementation-flow.svg)

Both programs already use worker-owned chunks. NuSTAR exposes the stage results as `SeedBatch`, `WindowBatch`, `CandidateSet`, and `FragmentEffects`, with reusable workspaces and a separate ordered-commit runtime. The figure expands the normal single-end path; paired-end mapping uses corresponding typed stages.

### Implementation differences

[![Source-level differences in suffix sorting, window reset, read-loop state, and BAM sorting](docs/assets/implementation-differences.gif)](docs/assets/implementation-differences.svg)

[Static version](docs/assets/implementation-differences.svg) · [Source comparison and code links](docs/IMPLEMENTATION_COMPARISON.md) · [Figure source](tools/documentation/render_implementation_figures.py)

- **Suffix sorting:** STAR sorts prefix-range chunks with `qsort`; NuSTAR starts with `libsais64`, then normalizes STAR-equivalent suffix groups.
- **Window reset:** STAR clears both full bin tables for each read; NuSTAR logs and resets the bins that were touched.
- **Pipeline state:** STAR methods share `ReadAlign` fields; NuSTAR passes explicit results between stages. The recursive stitching kernel is retained.
- **BAM sorting:** STAR sorts coordinate bins and concatenates them; NuSTAR sorts buffered runs and merges them with a priority queue.

The figures describe source-level mechanisms. Current runtime, memory use, and output agreement remain to be measured.

| STAR feature or workflow | NuSTAR implementation | Status and limits |
| :--- | :--- | :--- |
| Reference indexing and annotation | `reference/`, `index/`, `annotation/` | Implemented with a native index format. STAR binary indexes and shared-memory load/remove modes are not supported. |
| Seed search, windows, stitching, and selection | `seed_search/`, `window/`, `stitching/`, `selection/` | Implemented for single-end and paired-end spliced alignment, including optional overlap merge/remapping. Most configurable STAR scoring and filtering options are not exposed. |
| Basic two-pass alignment | `runtime/two_pass.cpp` and junction insertion in `index/` | Implemented: discover junctions, extend the index, and remap the reads. |
| STARsolo | `solo/` | STAR-derived code, refactored and adapted to NuSTAR. Supports barcode correction, UMI deduplication, gene/SJ/velocity counting, cell filtering, and SmartSeq. `Transcript3p` has less validation. |
| SAM/BAM output, junctions, and statistics | `sam/`, `bam/`, `output/`, `splice_junction/`, `statistics/` | Implemented, including sorted BAM and corrected Solo tags. See the [output limits](docs/USAGE.md#compatibility). |
| Input handling | `reads/`, `bam/source.cpp` | Plain FASTQ and direct SAM/BAM remapping; SmartSeq manifests. Compressed FASTQ and `readFilesCommand` are not supported. |
| Chimeric/fusion alignment and WASP/VCF workflows | — | Not implemented. |
| Transcriptome BAM, bulk `GeneCounts`, and signal tracks | — | Not implemented. |

The [usage guide](docs/USAGE.md#compatibility) lists the remaining input, output, and parameter limits. The [code guide](CODE_STRUCTURE.md) explains module boundaries and the Solo adaptation; [solo/UPSTREAM.json](solo/UPSTREAM.json) records its upstream source.

## Roadmap

The proposed next steps below follow the gaps in implementation and current validation. They have no release dates yet; feature requests and contributions can help set priorities.

| Area | Next step | How to check the result |
| :--- | :--- | :--- |
| Reproducible benchmarks | Fill in the benchmark template using paired runs of the current binaries | Record commands, binary hashes, hardware, cache conditions, repeated timings, and output comparisons |
| Bulk RNA-seq validation | Add large single-end and paired-end record comparisons and thread-scaling runs | Compare full alignment records, junctions, and non-time statistics with STAR |
| Two-pass performance | Measure the updated insertion and merge stages | Repeat identical inputs with both engines; compare both passes and final outputs |
| Single-cell coverage | Test more real datasets and chemistries, corrected BAM tags, and less-tested mode combinations | Compare raw/filtered matrices, summaries, and corrected tags; add truth-based evaluations where data allow |
| STAR feature coverage | Scope compressed FASTQ, additional output modes, and advanced alignment workflows | Define supported options and comparison fixtures before implementing each feature |

The unimplemented features above are candidates for future work, not a promise that full STAR compatibility will arrive in a particular release.

## AI use and human review

NuSTAR's maintainers report AI-assisted development and comprehensive human code review in the project documentation. This documentation refresh does not independently verify a comprehensive audit of every change in the updated revision.

No separate audit report listing reviewers, reviewed revisions, and scope is included here. The source tests and [historical comparison report](docs/STAR_COMPARISON.md) describe output checks and limitations; neither establishes that every STAR behavior has been reproduced. Current benchmark and output-comparison results remain pending.

## Questions and contributions

For questions about NuSTAR—or STAR behavior we can help explain—open an issue in [this repository](https://github.com/BioCode-Commons/NuSTAR/issues). Include the command, NuSTAR version or commit, relevant logs, and what you expected to happen. For a STAR comparison, include the STAR version and both commands.

Contributions are welcome: bug reports, small reproducible datasets, documentation fixes, tests, and implementation work. The [contribution guide](CONTRIBUTING.md) explains how to build, run the comparison suites, and prepare a pull request. For a large change, open an issue first to discuss scope and expected behavior.

## Acknowledgments and references

NuSTAR builds on the algorithms and methods developed by **Alexander Dobin and the STAR contributors**. We thank them for STAR and their contributions to the RNA-seq community. Please cite the original paper when using NuSTAR:

**Dobin A, Davis CA, Schlesinger F, et al.** STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013;29(1):15–21. [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635).

- **Original implementation:** [alexdobin/STAR](https://github.com/alexdobin/STAR)
- **Behavior baseline and Solo provenance:** [STAR 2.7.11b](https://github.com/alexdobin/STAR/tree/2.7.11b) and [solo/UPSTREAM.json](solo/UPSTREAM.json)
- **STAR-derived Solo code:** [retained STAR license](solo/STAR.LICENSE)
- **Bundled dependencies:** HTSlib, libsais, and Simple Good–Turing code; see [third_party/](third_party/) for source and license information.
