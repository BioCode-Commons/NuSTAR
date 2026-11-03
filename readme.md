# NuSTAR

**NuSTAR is a refactored implementation of STAR, a widely used RNA-seq aligner.** Using AI-assisted development, we have refactored and improved STAR’s implementation to support its continued maintenance. We thank Dr. Alexander Dobin and the STAR contributors for their work and contributions to the community over the past decade. Although AI was used during development, we have **conducted comprehensive human code reviews and testing** to verify that NuSTAR produces the same outputs as the original STAR for the workflows tested. For supported workflows, NuSTAR is designed to preserve STAR’s user-facing interface and behavior.

NuSTAR does not yet support all STAR features. Currently supported workflows include **bulk RNA-seq** alignment in both one-pass and two-pass modes, as well as single-cell RNA-seq analysis with **STARsolo**. Advanced analysis, such as XXX and XXX Chimeric are Not supported yet. This is still a partial implementation of STAR. Use the options below rather than assuming every STAR command will work. The [code guide](CODE_STRUCTURE.md) explains where the work happens. [Validation notes](docs/VALIDATION.md) explain how we compare results. [Runtime and correctness](docs/STAR_COMPARISON.md) compares the measured workloads with STAR and lists the remaining gaps. The [full-data tests](docs/DATA_FEATURE_VALIDATION.md) include real PBMC Solo results. The [junction-padding fix](docs/JUNCTION_PADDING_FIX.md) covers the completed full human two-pass run and the remaining runtime gap.

**NuSTAR builds on the algorithms and methods introduced in STAR. We therefore recommend that all users cite the original STAR publication:**

Dobin, Alexander, et al. *STAR: ultrafast universal RNA-seq aligner.* Bioinformatics 29.1 (2013): 15–21.

The original STAR repository is available at [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR).

Questions about STAR are also welcome in this repository. Having studied STAR’s source code in detail, we will do our best to answer them.


## Build

You need CMake 3.20 or newer, a C++20 compiler, a C99 compiler, and the zlib development package. OpenMP is optional. HTSlib and libsais are included in `third_party/`.

Run these commands inside `NuSTAR/`:

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 8
ctest --test-dir build/release --output-on-failure
```

The executable is `build/release/nustar`. Release builds use interprocedural optimization when the compiler supports it. Set `-DNUSTAR_ENABLE_RELEASE_IPO=OFF` when configuring to disable it.

## Build an index

```sh
build/release/nustar --runMode genomeGenerate \
  --genomeDir index/human \
  --genomeFastaFiles reference/GRCh38.fa \
  --sjdbGTFfile reference/genes.gtf \
  --sjdbOverhang 100 \
  --runThreadN 32
```

Choose `sjdbOverhang` for your reads; read length minus one is the usual starting point. Omit `--sjdbGTFfile` to build without annotation. Solo gene counting needs the annotation tables produced from a GTF.

NuSTAR uses its own `reference.bin` and `index.sa`. A STAR `Genome`/`SA`/`SAindex` directory cannot be used directly. For a small test genome, reduce `--genomeSAindexNbases`, for example to `3`.

## Align reads

Single-end:

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn reads.fastq \
  --runThreadN 32 --outFileNamePrefix results/sample.
```

Paired-end, with coordinate-sorted BAM:

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn mate1.fastq mate2.fastq \
  --runThreadN 32 --outFileNamePrefix results/paired. \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within --limitBAMsortRAM 2147483648
```

The default run mode is `alignReads`. FASTQ files must be uncompressed. `--outFileNamePrefix` is a literal prefix: `results/sample.` produces `results/sample.Aligned.out.sam`; `results/sample/` puts files in that directory.

Common outputs are `Aligned.out.sam`, `Aligned.out.bam`, or `Aligned.sortedByCoord.out.bam`, plus `SJ.out.tab` and `Log.final.out`. Use `--outSAMtype None` to skip alignment records while keeping junctions, statistics, and any requested Solo counts. Request both BAM orders with `--outSAMtype BAM Unsorted SortedByCoordinate`.

### Two-pass alignment

Add these options to an alignment command:

```sh
--twopassMode Basic --twopass1readsN -1 --sjdbInsertSave Basic
```

The first pass finds junctions. NuSTAR inserts the accepted junction sequences into the index, then rereads the input for the second pass. Flanks may extend into chromosome padding, as in STAR; see the [boundary fix and tests](docs/JUNCTION_PADDING_FIX.md). `-1` uses all reads in pass one; a positive value limits that pass. `--sjdbInsertSave All` also saves the updated native index. Input must be available for both passes.

### SAM and BAM input

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn input.bam --readFilesType SAM SE \
  --readFilesSAMattrKeep All \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results/remapped.
```

Use `SAM PE` for paired input. Despite the option name, NuSTAR detects and reads both SAM and BAM directly. It remaps the sequences; it does not just convert the file format. Paired input must contain adjacent records with the same read name and complementary mate flags. A coordinate-sorted paired BAM usually needs name grouping first. Each record needs sequence and quality data; input records are not implicitly discarded just because they are secondary or supplementary.

## Single-cell reads

For a typical 10x 3' v3 layout, R2 contains RNA and R1 contains a 16-base cell barcode followed by a 12-base UMI:

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn RNA_R2.fastq barcode_R1.fastq \
  --soloType CB_UMI_Simple --soloCBwhitelist whitelist.txt \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 --soloFeatures Gene \
  --soloCellFilter EmptyDrops_CR \
  --runThreadN 32 --outFileNamePrefix results/cells. \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CR CY UR UY CB UB GX GN
```

The barcode file goes last. `soloBarcodeReadLength 0` disables the exact barcode-read length check; the selected CB and UMI positions must still exist. Use the whitelist and layout for the actual chemistry. For v2, the usual UMI length is 10 rather than 12.

Raw and filtered matrices are written below `<prefix>Solo.out/`, grouped by feature. Corrected `CB`/`UB` tags for counting modes need `BAM SortedByCoordinate`: correction finishes after all reads have been counted. `CR`/`CY` and `UR`/`UY` hold the raw barcode/UMI and their qualities.

To read barcodes already stored in SAM/BAM, use `--readFilesType SAM SE` (or `PE`) and add:

```sh
--soloType CB_UMI_Simple --soloCBwhitelist whitelist.txt \
--soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY
```

For SmartSeq, provide a tab-separated manifest with three columns: mate 1 path, mate 2 path (or `-` for single-end), and cell/read-group information. Do not mix single-end and paired-end rows.

```sh
build/release/nustar --genomeDir index/human \
  --readFilesManifest cells.tsv --soloType SmartSeq \
  --soloUMIdedup Exact --soloCellFilter None \
  --outSAMtype None --outFileNamePrefix results/smartseq.
```

To filter an existing raw matrix without mapping reads:

```sh
build/release/nustar --runMode soloCellFiltering \
  results/cells.Solo.out/Gene/raw results/refiltered/ \
  --soloCellFilter TopCells 3000
```

## Options

Names and argument forms follow STAR where an equivalent exists. The old `--threads` alias and positional `alignSingleSam`, `alignPairedSam`, `genomeGenerate`, and `build-index` commands have been removed. Use `--runThreadN`, `--runMode`, `--readFilesIn`, and `--outFileNamePrefix` instead.

### Shared and index options

| Option | Default | Use |
| --- | --- | --- |
| `--help`, `--version` | — | Print help or the NuSTAR version; use either on its own. |
| `--runMode` | `alignReads` | `alignReads`, `genomeGenerate`, or `soloCellFiltering RAW_DIR OUTPUT_PREFIX`. |
| `--runThreadN N` | `1` | Worker count, from 1 to 256 for mapping and index building. |
| `--genomeDir DIR` | `./GenomeDir/` | Native index directory. |
| `--genomeFastaFiles FILE...` | required for index building | One or more reference FASTA files. |
| `--sjdbGTFfile FILE` | `-` | GTF used during index building; `-` means no GTF. |
| `--sjdbOverhang N` | `100` | Junction flank length, used when building or extending the index. |
| `--genomeChrBinNbits N` | `18` | Chromosome padding/bin bits, 8–30. |
| `--genomeSAindexNbases N` | `14` | SA prefix length, 1–15. Reduce for small genomes. |
| `--limitSjdbInsertNsj N` | `1000000` | Total junction count after insertion, including existing junctions. The full human test requires at least 1,237,619; both programs were tested with 2,000,000. |
| `--sjdbGTFchrPrefix TEXT` | `-` | Prefix for GTF chromosome names; `-` adds nothing. |
| `--sjdbGTFfeatureExon NAME` | `exon` | GTF feature used as an exon. |
| `--sjdbGTFtagExonParentTranscript NAME` | `transcript_id` | Transcript attribute. |
| `--sjdbGTFtagExonParentGene NAME` | `gene_id` | Gene ID attribute. |
| `--sjdbGTFtagExonParentGeneName NAME...` | `gene_name` | Gene name attributes. |
| `--sjdbGTFtagExonParentGeneType NAME...` | `gene_type gene_biotype` | Gene type attributes. |

### Alignment and output options

| Option | Default | Use |
| --- | --- | --- |
| `--readFilesIn FILE...` | required unless using a manifest | One or two RNA FASTQ files, followed by a separate barcode file when needed; or one SAM/BAM file. |
| `--readFilesManifest FILE` | none | SmartSeq input groups. |
| `--readFilesType Fastx` or `SAM SE/PE` | `Fastx` | FASTQ, or SAM/BAM with the specified mate layout. FASTA reads are not implemented. |
| `--readFilesSAMattrKeep TAG...` | `All` | BAM input tags to retain: `All`, `None`, or tag names. The SAM text reader currently retains all input tags. |
| `--readMapNumber N` | `-1` | Maximum reads/fragments to map; `-1` means all. |
| `--clip5pNbases N...`, `--clip3pNbases N...` | `0` | Fixed clipping at each end. Supply both mate values for paired reads. |
| `--outFileNamePrefix PREFIX` | `./` | Prefix for output paths. |
| `--outSAMtype` | `SAM` | `SAM`, `None`, or `BAM` followed by `Unsorted`, `SortedByCoordinate`, or both. |
| `--outSAMattributes TAG...` | `Standard` | `Standard` expands to `NH HI AS nM`; otherwise list tags explicitly. |
| `--outSAMunmapped MODE` | `None` | `None` or `Within`; `KeepPairs` is not implemented. |
| `--outSAMorder MODE` | `Paired` | `Paired` or `PairedKeepInputOrder`. NuSTAR commits unsorted batches in input order in either mode. |
| `--outBAMcompression N` | `1` | BGZF compression level, -1 through 9. |
| `--limitBAMsortRAM BYTES` | `268435456` | Sort-record budget, not a total-process memory limit. Zero also selects 256 MiB in NuSTAR. |
| `--twopassMode MODE` | `None` | `None` or `Basic`. |
| `--twopass1readsN N` | `-1` | First-pass read limit; use `-1` or a positive count. |
| `--sjdbInsertSave MODE` | `Basic` | `Basic` or `All`; `All` saves the updated index. |

Ordinary output tags are `NH HI AS nM RG`. Single-end mapping also supports `NM MD jM jI XS`. Solo adds `CR CY UR UY CB UB GX GN sM sS sQ sF`. `GX`/`GN` currently require Solo. `outSAMattributes None` and `All` are not implemented.

BAM sorting spills runs to disk. During mapping, the record cache is capped at the smaller of the requested sort budget and 8 MiB per thread; other buffers and the genome index consume additional memory. The `limitBAMsortRAM 0` behavior differs from STAR's automatic budget choice.

### Solo options

These options apply to an enabled Solo mode. Defaults are in [solo/configuration/defaults.hpp](solo/configuration/defaults.hpp).

| Option | Default | Use |
| --- | --- | --- |
| `--soloType MODE` | `None` | `None`, `CB_UMI_Simple`, `CB_UMI_Complex`, `CB_samTagOut`, or `SmartSeq`. `CB_samTagOut` writes barcode tags without matrix counting. |
| `--soloCBtype TYPE` | `Sequence` | `Sequence` or `String`; string barcodes require `--soloCBwhitelist None`. |
| `--soloCBwhitelist FILE...` | `-` | Barcode whitelist(s); `None` disables matching for simple layouts. The default `-` means no file was supplied and is not valid for simple-layout counting. |
| `--soloCBstart N`, `--soloCBlen N` | `1`, `16` | Simple-layout CB start (1-based) and length. |
| `--soloUMIstart N`, `--soloUMIlen N` | `17`, `10` | Simple-layout UMI start (1-based) and length. |
| `--soloBarcodeReadLength N` | `1` | `1`: require CB length + UMI length; `0`: skip exact length check; larger values: required length. |
| `--soloBarcodeMate N` | `0` | `0`: separate barcode file; `1` or `2`: barcode is in an RNA mate. Set fixed clipping for that mate. |
| `--soloCBposition DESCRIPTOR...` | `-` | Complex-layout CB pieces in STAR's `anchor_offset_anchor_offset` form. |
| `--soloUMIposition DESCRIPTOR` | `-` | Complex-layout UMI position. |
| `--soloAdapterSequence SEQUENCE` | `-` | Adapter anchor for complex layouts. |
| `--soloAdapterMismatchesNmax N` | `1` | Allowed adapter mismatches. |
| `--soloCBmatchWLtype MODE` | `1MM_multi` | `Exact`, `1MM`, `1MM_multi`, `1MM_multi_pseudocounts`, `1MM_multi_Nbase_pseudocounts`, or complex-layout `EditDist_2`. |
| `--soloInputSAMattrBarcodeSeq TAG...` | `-` | Input sequence tags, commonly `CR UR`. |
| `--soloInputSAMattrBarcodeQual TAG...` | `-` | Input quality tags, commonly `CY UY`. |
| `--soloStrand MODE` | `Forward` | `Forward`, `Reverse`, or `Unstranded`. |
| `--soloFeatures FEATURE...` | `Gene` | `Gene`, `GeneFull`, `GeneFull_Ex50pAS`, `GeneFull_ExonOverIntron`, `SJ`, `Velocyto`, `VelocytoSimple`, or `Transcript3p`. |
| `--soloMultiMappers MODE...` | `Unique` | `Unique`, `Uniform`, `Rescue`, `PropUnique`, or `EM`. |
| `--soloUMIdedup MODE...` | `1MM_All` | `Exact`, `1MM_All`, `1MM_Directional`, `1MM_Directional_UMItools`, `1MM_CR`, or `NoDedup`. |
| `--soloUMIfiltering MODE` | `-` | `-`, `MultiGeneUMI`, `MultiGeneUMI_All`, or `MultiGeneUMI_CR`. The last requires `1MM_CR` deduplication. |
| `--soloCellFilter MODE ...` | `CellRanger2.2 3000 0.99 10` | `None`, `TopCells N`, `CellRanger2.2`, or `EmptyDrops_CR`. Numeric settings are described below. |
| `--soloOutFileNames DIR FEATURES BARCODES MATRIX` | `Solo.out/ features.tsv barcodes.tsv matrix.mtx` | Output directory suffix and matrix component names. |
| `--soloOutFormatFeaturesGeneField3 TEXT` | `Gene Expression` | Third field of gene feature rows; quote text containing spaces. |
| `--soloCellReadStats MODE` | `None` | `None` or `Standard`. |
| `--soloClusterCBfile FILE` | `-` | Barcode clusters for `Transcript3p` counting. This path has less validation than gene counting. |

Custom `--runRNGseed` values are not supported.

Not every combination is valid. For example, `CB_UMI_Complex` uses `Exact`, `1MM`, or `EditDist_2` matching; `CB_samTagOut` uses `Exact` or `1MM`. SmartSeq does not support `Velocyto` or multimapper redistribution. Check the error message when choosing a non-default combination.

`CellRanger2.2` accepts three numbers: expected cells, maximum percentile, and maximum/minimum ratio. Supplying the mode alone uses `3000 0.99 10`.

`EmptyDrops_CR` accepts ten numbers in this order: expected cells, maximum percentile, maximum/minimum ratio, lower ambient rank, upper ambient rank, minimum UMI count, minimum fraction of median UMIs, maximum candidate count, FDR, and simulation count. Its defaults are `3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000`.

### NuSTAR-only controls

| Option or command | Use |
| --- | --- |
| `--io-backend mmap` or `buffered` | FASTQ transport; default `mmap`. |
| `--profile-stages` | Detailed mapping-stage timing, currently single-end only. |
| `--profile-runtime` | Runtime timing for single-end or paired-end runs. |
| `--profileGenomeGenerate` | Index-building timing. |
| `upgradeIndex SOURCE_DIR NEW_DIR` | Upgrade an older native index without rebuilding the suffix array. |

These controls have no direct STAR equivalent. Benchmark and trace scripts also have their own tool-specific arguments; use their `--help` output.

## What works, and what is missing

| Area | Implemented | Missing or limited |
| --- | --- | --- |
| Index | FASTA/GTF loading, packed suffix array, annotation tables, incremental junction insertion | Loading STAR's binary index; STAR shared-memory load/remove modes. |
| Mapping | Single-end and paired-end spliced alignment, seed/window/stitch/select stages, fixed clipping, `Basic` two-pass | Chimeric/fusion alignment, WASP/VCF workflows, adapter clipping, most configurable STAR scoring and filtering options. |
| Input | Plain FASTQ; direct SAM/BAM; SmartSeq manifests | Compressed FASTQ, FASTA reads, comma-separated file groups, `readFilesCommand`, general manifest input outside SmartSeq. |
| Output | SAM, unsorted/sorted BAM, junction tables, mapping statistics | CRAM as a supported interface, BAI/CSI generation, transcriptome BAM, bulk `GeneCounts`, signal tracks, all STAR header/tag customization. BAM read names are limited to 254 bytes and CIGARs to 65,535 operations. |
| Solo | Barcode extraction/correction, UMI deduplication, gene/SJ/velocity counting, cell filtering, matrix output, SmartSeq, corrected BAM tags | Full validation of every mode/parameter combination on real single-cell datasets; `Transcript3p` remains less tested. |
| Parallelism and profiling | Shared read-only index, worker-local workspaces, ordered output, BAM compression/sorting, stage timing | Detailed paired-end stage timing. |

Matching one dataset does not prove that every STAR behavior is reproduced. The tests compare records, junctions, statistics, matrices, and selected internal stages. Performance depends on input and output mode; the stored measurements in [docs/VALIDATION.md](docs/VALIDATION.md) are examples, not a promise for every workload.

The Solo algorithms retain STAR-derived code and its [license](solo/STAR.LICENSE); [solo/UPSTREAM.json](solo/UPSTREAM.json) records the source revision. Third-party libraries keep their own licenses under `third_party/`.
