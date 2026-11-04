# NuSTAR usage and parameter reference

[Back to the README](../readme.md) · [Build and quick start](../readme.md#quick-start) · [Contributing](../CONTRIBUTING.md)

Run the examples from the repository root after building `build/release/nustar` and creating a NuSTAR index. Supply your own input files and create the output parent directory, such as `results/`, before running a command.

## Workflows

### Single-end alignment

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn reads.fastq \
  --runThreadN 8 \
  --outFileNamePrefix results/single.
```

The default alignment output is SAM. Use `--outSAMtype None` to retain junctions, statistics, and requested Solo counts without writing alignment records.

### Two-pass alignment

Add the following options to a single-end or paired-end alignment command:

```sh
--twopassMode Basic --twopass1readsN -1 --sjdbInsertSave Basic
```

Pass one discovers junctions. NuSTAR inserts accepted junction sequences into the index, then rereads the input for pass two.

- `--twopass1readsN -1` uses all reads in pass one; a positive value limits that pass.
- `--sjdbInsertSave All` also saves the updated native index.
- Input must remain available for both passes.

### Single-cell: 10x 3′ v3

For this layout, **R2 contains RNA** and **R1 contains a 16-base cell barcode followed by a 12-base UMI**. Supply the barcode FASTQ last and use the whitelist for your chemistry.

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn RNA_R2.fastq barcode_R1.fastq \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist whitelist.txt \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene \
  --soloCellFilter EmptyDrops_CR \
  --runThreadN 8 \
  --outFileNamePrefix results/cells. \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CR CY UR UY CB UB GX GN
```

Raw and filtered matrices appear under `results/cells.Solo.out/`, grouped by feature. `--soloBarcodeReadLength 0` disables the exact barcode-read length check; the selected barcode and UMI positions must still exist. For a typical v2 layout, the UMI length is 10 rather than 12.

> [!NOTE]
> Corrected `CB`/`UB` tags in counting modes require `BAM SortedByCoordinate`, because correction finishes after counting all reads. `CR`/`CY` and `UR`/`UY` contain the raw barcode/UMI sequences and qualities.

<details>
<summary><strong>Remap SAM/BAM input</strong></summary>

```sh
build/release/nustar --genomeDir index/human \
  --readFilesIn input.bam \
  --readFilesType SAM SE \
  --readFilesSAMattrKeep All \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix results/remapped.
```

Despite the option name, NuSTAR detects and reads both SAM and BAM directly. It **remaps the sequences**, rather than only converting the format.

Use `SAM PE` for paired input. Paired records must be adjacent, share the same read name, and have complementary mate flags; coordinate-sorted paired BAM usually needs name grouping first. Each record needs sequence and quality data. Secondary and supplementary records are not implicitly discarded.

For barcodes already stored in SAM/BAM tags, add:

```sh
--soloType CB_UMI_Simple --soloCBwhitelist whitelist.txt \
--soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY
```

</details>

<details>
<summary><strong>Analyze SmartSeq reads</strong></summary>

Provide a tab-separated manifest with three columns: mate 1 path, mate 2 path (or `-` for single-end), and cell/read-group information. Do not mix single-end and paired-end rows.

```sh
build/release/nustar --genomeDir index/human \
  --readFilesManifest cells.tsv \
  --soloType SmartSeq \
  --soloUMIdedup Exact \
  --soloCellFilter None \
  --outSAMtype None \
  --outFileNamePrefix results/smartseq.
```

</details>

<details>
<summary><strong>Filter an existing raw matrix</strong></summary>

Cell filtering can run independently of read alignment:

```sh
build/release/nustar --runMode soloCellFiltering \
  results/cells.Solo.out/Gene/raw results/refiltered/ \
  --soloCellFilter TopCells 3000
```

</details>

## Outputs

Paths below are appended to `--outFileNamePrefix`. Alignment files depend on the selected output mode; Solo matrices require a counting workflow.

| Output | Contents |
| :--- | :--- |
| `Aligned.out.sam` | SAM alignment records; default output |
| `Aligned.out.bam` | Unsorted BAM records |
| `Aligned.sortedByCoord.out.bam` | Coordinate-sorted BAM records |
| `SJ.out.tab` | Splice-junction observations |
| `Log.final.out` | Mapping statistics |
| `Solo.out/<feature>/raw/` | Raw feature-by-barcode matrix files |
| `Solo.out/<feature>/filtered/` | Filtered matrix files when cell filtering is enabled |

Request both BAM orders with `--outSAMtype BAM Unsorted SortedByCoordinate`. BAM index files (`BAI`/`CSI`) are not generated.

## Compatibility

NuSTAR follows STAR option names and argument forms where implemented. It is a partial implementation; a successful comparison on one dataset does not establish compatibility with all STAR workflows.

| Area | Implemented | Missing or limited |
| --- | --- | --- |
| Index | FASTA/GTF loading, packed suffix array, annotation tables, incremental junction insertion | Loading STAR's binary index; STAR shared-memory load/remove modes. |
| Mapping | Single-end and paired-end spliced alignment, seed/window/stitch/select stages, fixed clipping, `Basic` two-pass | Chimeric/fusion alignment, WASP/VCF workflows, adapter clipping, most configurable STAR scoring and filtering options. |
| Input | Plain FASTQ; direct SAM/BAM; SmartSeq manifests | Compressed FASTQ, FASTA reads, comma-separated file groups, `readFilesCommand`, general manifest input outside SmartSeq. |
| Output | SAM, unsorted/sorted BAM, junction tables, mapping statistics | CRAM as a supported interface, BAI/CSI generation, transcriptome BAM, bulk `GeneCounts`, signal tracks, all STAR header/tag customization. BAM read names are limited to 254 bytes and CIGARs to 65,535 operations. |
| Solo | Barcode extraction/correction, UMI deduplication, gene/SJ/velocity counting, cell filtering, matrix output, SmartSeq, corrected BAM tags | Full validation of every mode/parameter combination on real single-cell datasets; `Transcript3p` remains less tested. |
| Parallelism and profiling | Shared read-only index, worker-local workspaces, ordered output, BAM compression/sorting, stage timing | Detailed paired-end stage timing. |

## Parameter reference

Use `--runThreadN`, `--runMode`, `--readFilesIn`, and `--outFileNamePrefix`. The legacy `--threads` alias and positional `alignSingleSam`, `alignPairedSam`, `genomeGenerate`, and `build-index` commands have been removed.

Expand a group for defaults, accepted values, and constraints.

<details>
<summary><strong>Shared and index options</strong></summary>

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

</details>

<details>
<summary><strong>Alignment and output options</strong></summary>

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

</details>

<details>
<summary><strong>Solo options</strong></summary>

These options apply to an enabled Solo mode. Defaults are in [solo/configuration/defaults.hpp](../solo/configuration/defaults.hpp).

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

</details>

<details>
<summary><strong>NuSTAR-only controls</strong></summary>

| Option or command | Use |
| --- | --- |
| `--io-backend mmap` or `buffered` | FASTQ transport; default `mmap`. |
| `--profile-stages` | Detailed mapping-stage timing, currently single-end only. |
| `--profile-runtime` | Runtime timing for single-end or paired-end runs. |
| `--profileGenomeGenerate` | Index-building timing. |
| `upgradeIndex SOURCE_DIR NEW_DIR` | Upgrade an older native index without rebuilding the suffix array. |

These controls have no direct STAR equivalent. Benchmark and trace scripts also have their own tool-specific arguments; use their `--help` output.

</details>
