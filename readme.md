<div align="center">

# NuSTAR

**NuSTAR is a reimplementation of STAR, a popular RNA-seq aligner.**

[Quick start](#quick-start) · [Codebase](#codebase) · [Benchmarks](#benchmarks) · [Roadmap](#roadmap) · [Usage](docs/USAGE.md) · [Community](https://biocode-commons.github.io/)

</div>

---

NuSTAR supports single-end and paired-end reads, overlap refinement, two-pass mapping, and STARsolo.

## Quick start

### 1. Build

Requires CMake 3.20+, C++20 and C99 compilers, zlib development files, and Python 3. OpenMP is optional; HTSlib and libsais are bundled.

```sh
git clone https://github.com/BioCode-Commons/NuSTAR.git
cd NuSTAR

cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 8
ctest --test-dir build/release --output-on-failure
```

<details>
<summary><strong>Build configuration</strong></summary>

Disable release IPO with `-DNUSTAR_ENABLE_RELEASE_IPO=OFF`. To build without tests, use `-DBUILD_TESTING=OFF` and skip `ctest`; Python is then optional.

See [build options and tests](CONTRIBUTING.md#tests).

</details>

### 2. Build a genome index

Supply a reference FASTA and matching GTF.

```sh
build/release/nustar --runMode genomeGenerate \
  --genomeDir index/human \
  --genomeFastaFiles reference/GRCh38.fa \
  --sjdbGTFfile reference/genes.gtf \
  --sjdbOverhang 100 \
  --runThreadN 8
```

Set `--sjdbOverhang` to read length minus one. Solo gene counting needs GTF annotation tables; small genomes need a smaller `--genomeSAindexNbases`.

### 3. Align paired-end reads

Use a native NuSTAR index and uncompressed FASTQ reads.

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

Produces sorted BAM, splice junctions, and mapping statistics. See [usage and parameters](docs/USAGE.md) for single-end, two-pass, overlap, and single-cell workflows.

## Codebase

| Module | Responsibility |
| --- | --- |
| [reference/](reference/), [index/](index/), [annotation/](annotation/) | Native indexes, suffix arrays, and gene/junction annotations |
| [seed_search/](seed_search/), [window/](window/), [stitching/](stitching/) | Seed lookup, window construction, and candidate alignments |
| [pipeline/](pipeline/), [selection/](selection/) | Single/paired mapping, overlap refinement, selection, and result projection |
| [runtime/](runtime/) | Input batches, worker workspaces, ordered output, and two-pass coordination |
| [reads/](reads/), [sam/](sam/), [bam/](bam/), [output/](output/) | Read input, alignment records, compression, and coordinate sorting |
| [solo/](solo/) | Barcode correction, UMI deduplication, gene counting, and cell filtering |
| [tests/](tests/), [tools/](tools/) | Regression tests, tracing, and benchmarking |

Read the [code guide](CODE_STRUCTURE.md) for execution flow and entry points; see [supported modes and limits](docs/USAGE.md#compatibility).

## Benchmarks

Updated measurements are pending.

| Item | Status |
| --- | --- |
| Datasets and reference versions | To be supplied |
| Hardware, builds, and commands | To be supplied |
| Runtime and peak memory | To be measured |
| Alignment and counting outputs | To be checked |

[Benchmark template](docs/BENCHMARKS.md).

## Roadmap

- Publish repeatable runtime, memory, and output checks.
- Expand bulk and single-cell validation across datasets and modes.
- Profile two-pass indexing and extend input/output support.

## AI use and review

NuSTAR was developed with AI assistance. Maintainers report comprehensive human code review.

## Questions and contributions

Visit [BioCode Commons](https://biocode-commons.github.io/) for projects, contribution guidance, and contact information.

## Acknowledgments and references

NuSTAR builds on [STAR](https://github.com/alexdobin/STAR) and STARsolo by Alexander Dobin and collaborators. Please cite the original work using [references.bib](references.bib).

Upstream attribution is recorded in [solo/UPSTREAM.json](solo/UPSTREAM.json); licenses for STAR-derived code, HTSlib, libsais, and Simple Good–Turing are retained in [solo/STAR.LICENSE](solo/STAR.LICENSE) and [third_party/](third_party/).
