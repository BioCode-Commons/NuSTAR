# Reading the NuSTAR code

Start with [cli/align_reads.cpp](cli/align_reads.cpp). It turns command-line arguments into a run configuration and chooses the single-end or paired-end runner. The alignment algorithms live below that layer; they do not parse arguments or open output files.

Most alignment code uses the namespace `nustar::star2711b`. The suffix names the STAR version whose behavior it follows. Files such as `star_search.cpp` and `star_window.cpp` are NuSTAR's implementations of those stages, not wrappers that call a STAR executable.

## From a read to an alignment

The main path is:

```text
CLI -> runtime setup -> input batches
    -> prepare sequence -> seed search -> windows -> stitching
    -> select/filter alignments -> SAM records + junctions + statistics
    -> ordered output
    -> merge Solo counts -> correct tags -> finish coordinate-sorted BAM
```

The Solo steps run only when requested. A run without alignment output still computes the alignments needed for junctions, statistics, and counting.

| Step | Read these files | What the code does |
| --- | --- | --- |
| Parse the command | [cli/main.cpp](cli/main.cpp), [cli/align_reads.cpp](cli/align_reads.cpp) | Select a run mode, parse STAR-style arguments, and choose the mate layout and output format. |
| Set up a run | [runtime/single_alignment.cpp](runtime/single_alignment.cpp), [runtime/paired_alignment.cpp](runtime/paired_alignment.cpp) | Open reference/index resources, connect the input and output, create Solo state if needed, and run the worker loop. |
| Read fragments | [reads/fastq.cpp](reads/fastq.cpp), [reads/dispatch.cpp](reads/dispatch.cpp), [bam/source.cpp](bam/source.cpp) | Read FASTQ or SAM/BAM and form one fragment containing one or two RNA reads. SAM/BAM input is converted back to read orientation before remapping. |
| Prepare bases | [reads/fastq_prepare.cpp](reads/fastq_prepare.cpp) | Apply fixed clipping and create the encoded sequence views used by mapping. |
| Find seeds | [seed_search/star_search.cpp](seed_search/star_search.cpp) | Search the suffix array for maximal mappable pieces, retaining read positions, lengths, and match intervals. |
| Build windows | [window/star_window.cpp](window/star_window.cpp) | Place seed hits into genomic bins, form candidate windows, and arrange their seeds and anchors for stitching. |
| Stitch candidates | [stitching/single_stitch.cpp](stitching/single_stitch.cpp), [stitching/paired_stitch.cpp](stitching/paired_stitch.cpp), [stitching/stitch_window.cpp](stitching/stitch_window.cpp) | Explore compatible seed chains within each window. Paired stitching also handles the relationship between mates. |
| Extend a chain | [stitching/stitch_segments.cpp](stitching/stitch_segments.cpp), [stitching/extend_alignment.cpp](stitching/extend_alignment.cpp) | Extend matches and score mismatches, indels, and splice junctions. |
| Choose alignments | [selection/star_select.cpp](selection/star_select.cpp) | Apply candidate selection and mapping filters, including multimapping decisions. |
| Describe results | [sam/star_single.cpp](sam/star_single.cpp), [sam/star_paired.cpp](sam/star_paired.cpp) | Project selected alignments into SAM fields and tags. These files do not compress BAM. |
| Gather side results | [splice_junction/star_junction.cpp](splice_junction/star_junction.cpp), [statistics/star_statistics.cpp](statistics/star_statistics.cpp) | Collect junction observations and mapping counters from the same mapping result. |
| Write records | [sam/encode.cpp](sam/encode.cpp), [output/file_sink.cpp](output/file_sink.cpp), [bam/codec.cpp](bam/codec.cpp), [bam/sink.cpp](bam/sink.cpp) | Encode text or BAM records, write chunks, compress BGZF blocks, and sort BAM when requested. |

The stage order is wired in [pipeline/single_end/mapping_composition.cpp](pipeline/single_end/mapping_composition.cpp) and [pipeline/paired_end/mapping_composition.cpp](pipeline/paired_end/mapping_composition.cpp). Their neighboring `pipeline.cpp` files add selection, record projection, junctions, and statistics. This is the place to review how one stage's output becomes the next stage's input.

## Overlapping mates

[pipeline/paired_end/star_overlap.cpp](pipeline/paired_end/star_overlap.cpp) implements `peOverlapNbasesMin` and `peOverlapMMp`. The paired pipeline calls it after ordinary paired candidate mapping and before selection. It searches the two mate offsets using STAR's match/mismatch rule, maps a merged read with a separate single-read workspace, projects blocks back into paired coordinates, and recalculates the paired scores. The original sequences remain available for SAM/BAM output.

[runtime/paired_alignment.cpp](runtime/paired_alignment.cpp) creates this helper only when overlap is enabled. It shares the seed finder, window builder and reference with the normal mapper. The helper does not invoke a complete single-end output pipeline or consume random selection state.

## Index and annotation

[cli/genome_generate.cpp](cli/genome_generate.cpp) coordinates index building. [reference/genome.cpp](reference/genome.cpp) owns the reference representation. [annotation/gtf_loader.cpp](annotation/gtf_loader.cpp) reads genes, transcripts, and exons; [annotation/junction_database.cpp](annotation/junction_database.cpp) builds the annotated junction data.

The suffix-array backend is selected with `--genomeSAalgorithm`: [index/wide_suffix_array.cpp](index/wide_suffix_array.cpp) supplies the default `libsais64` path, and [index/libsais40.cpp](index/libsais40.cpp) adapts the kernels in [third_party/libsais40/](third_party/libsais40/) to the owning five-byte workspace. [index/packed_suffix_array.cpp](index/packed_suffix_array.cpp) handles packed entries. [index/star_prefix.cpp](index/star_prefix.cpp) builds and queries the prefix search data. [index/star_rank.cpp](index/star_rank.cpp) preserves the required suffix ordering.

During construction, `PackedSuffixArray::packAndRelease` consumes temporary positions. The 64-bit path streams into a separate packed allocation while releasing consumed source pages on Linux. The five-byte path packs in place with a fixed 32 KiB buffer, transfers ownership of the original mapping, and unmaps unused tail pages on Unix. Prefix ranges are built after conversion, and the CLI releases construction objects before reloading the saved index. Suffix ranks and the v4 file format are unchanged.

[index/suffix_array_index.cpp](index/suffix_array_index.cpp) loads and saves the native index. [index/backend.hpp](index/backend.hpp) is the interface used by mapping code. The index representation stays behind that interface so seed search does not need to know how files were opened.

The reference and index are shared between workers. Their owners keep the backing storage alive; per-read stages use views into it. Per-read scratch buffers belong to workers and are reused rather than recreated for each read.

## Two-pass mapping

[runtime/two_pass.cpp](runtime/two_pass.cpp) runs pass one with junction collection, filters the junctions, extends the index, and starts pass two from the original input. [index/junction_insert.cpp](index/junction_insert.cpp) and [index/extension.cpp](index/extension.cpp) do the index update. Pass two uses the same mapping pipeline as an ordinary run.

This split matters when changing the code: junction discovery belongs to the mapping run, while index insertion belongs to `index/`. The two-pass coordinator should not implement another seed search or stitcher.

## Solo

Solo receives NuSTAR alignments. It does not run a second aligner.

| Directory or file | Responsibility |
| --- | --- |
| [solo/config.hpp](solo/config.hpp), [solo/session.hpp](solo/session.hpp) | Public configuration and run-level interface. The runtime uses these rather than the counting internals. |
| [solo/session.cpp](solo/session.cpp), [solo/worker.cpp](solo/worker.cpp) | Create shared and worker-local Solo state, adapt mapping results, finish counting, and make corrected tags available to BAM output. |
| [solo/configuration/](solo/configuration/) | Parse options, apply defaults, and check combinations. |
| [solo/barcode/](solo/barcode/) | Extract CB/UMI sequences, locate complex layouts, and match barcodes against whitelists. |
| [solo/annotation/](solo/annotation/) | Match alignment blocks to genes, exons, introns, and transcript features. |
| [solo/counting/](solo/counting/) | Record read assignments, merge workers, count features, and handle multimappers and SmartSeq. |
| [solo/umi/](solo/umi/) | Deduplicate UMIs using the selected rule. |
| [solo/filtering/](solo/filtering/) | Select cells, including EmptyDrops-style filtering. |
| [solo/matrix/](solo/matrix/) | Load and write matrices and summaries. |
| [solo/filter.cpp](solo/filter.cpp) | Run cell filtering on an existing matrix without mapping reads. |
| [solo/support/](solo/support/) | Small helpers shared by the Solo implementation. |

For example, [solo/annotation/alignment.cpp](solo/annotation/alignment.cpp) contains `Alignment::extractSpliceJunctions()`. It walks adjacent exon blocks, skips indels and mate gaps, and records intron starts and lengths. `chrStartLengthExtended()` computes an extended interval used by annotation. These functions describe an existing alignment; they do not choose how reads align.

Corrected `CB`/`UB` tags can depend on counts from later reads. Coordinate BAM output therefore retains the read identity and obtains final tags after Solo finishes. Raw barcode tags can be produced during mapping.

Solo's counting algorithms come from STAR. The source has been split into modules and adapted to NuSTAR's pipeline. [solo/UPSTREAM.json](solo/UPSTREAM.json) and [solo/STAR.LICENSE](solo/STAR.LICENSE) keep that origin visible.

## Threads and output order

[runtime/ordered_alignment.cpp](runtime/ordered_alignment.cpp) owns the mapping worker loop. Input access is serialized when a worker takes a batch. The expensive mapping work runs outside that input lock.

Each worker has its own workspace, random state, and temporary output buffers. The reference and suffix array are shared read-only. Solo also has per-worker recording state, followed by a merge after mapping.

Batches get sequence numbers. Workers may finish in a different order, but output commits and junction/statistic reductions follow input order. File sinks can spool chunks concurrently and concatenate them later. At most twice the worker count of batches can be in flight, which bounds the reorder backlog when one batch is slow.

Coordinate sorting is a separate step in [bam/sink.cpp](bam/sink.cpp). It writes bounded runs and merges them, rather than keeping the whole BAM in memory. The runtime releases mapping resources it no longer needs before the remaining counting and sorting work; shared owners can keep resources alive longer, including during two-pass runs.

## Other directories

| Directory | Responsibility |
| --- | --- |
| [common/](common/) | Shared result types, workspace contracts, and file identity helpers. |
| [profiling/](profiling/) | Stage and progress timing. Timings are observations, not inputs to mapping decisions. |
| [tools/trace/](tools/trace/) | Dump or replay stage inputs and outputs for inspection. |
| [tools/oracle/](tools/oracle/) | Build instrumented STAR binaries for stage comparisons. These are test tools, not runtime dependencies. |
| [tools/benchmark/](tools/benchmark/) | Compare runtimes and memory, with inputs and commands saved alongside results. |
| [tests/](tests/) | Unit tests, module-boundary checks, and comparisons against STAR. |
| [third_party/](third_party/) | libsais, HTSlib, and other imported code with their licenses. |

For a change to seed search, start with its stage test and then compare full alignments. For a change to counting, compare matrices as well as BAM tags. Equal final record counts alone can hide different alignments or gene assignments. [docs/STAR_COMPARISON.md](docs/STAR_COMPARISON.md#running-the-comparisons) lists the test entry points.
