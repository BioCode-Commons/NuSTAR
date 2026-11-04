# STAR and NuSTAR: flow and implementation differences

[Back to the README](../readme.md#implementation-status) · [Code guide](../CODE_STRUCTURE.md) · [Output comparisons](STAR_COMPARISON.md)

This comparison reads STAR **2.7.11b**, revision [`b1edc12`](https://github.com/alexdobin/STAR/tree/b1edc1208d91a53bf40ebae8669f71d50b994851), against NuSTAR revision [`fa16481`](https://github.com/BioCode-Commons/NuSTAR/tree/fa1648142582107665f1ff562d62666c6d8b61cb). The figures show specific mechanisms found in those sources. The abbreviated code below describes control flow; it is not a patch that can be applied between the repositories.

## End-to-end flows

![STAR and NuSTAR flows, including index setup, per-read and batch loops, stage boundaries, and finalization](assets/implementation-flow.svg)

Both implementations already process reads in worker-owned chunks and perform seed search, window construction, stitching, selection, and output. NuSTAR makes intermediate results and workspace ownership explicit. Splitting these interfaces does not by itself change the alignment algorithm.

The NuSTAR lane expands the normal single-end path through `runOrderedAlignment`; it is not a diagram of the separate profiling implementation. Its paired-end pipeline uses corresponding typed mapping stages and handles mate overlap as well. Error exits and STAR workflows outside NuSTAR's supported scope are omitted.

- **STAR control flow:** [chunk processing](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlignChunk_processChunks.cpp#L11), [per-read processing](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_oneRead.cpp#L78-L111), and [seed search followed by window/stitch processing](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_mapOneRead.cpp#L34-L118).
- **NuSTAR control flow:** [normal single-end runtime](../runtime/single_alignment.cpp), [ordered batch loop](../runtime/ordered_alignment.cpp), [candidate mapping composition](../pipeline/single_end/mapping_composition.cpp), and [selection, projection, and effects](../pipeline/single_end/pipeline.cpp).
- **Finalization:** STAR finishes Solo counting before [coordinate BAM sorting](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/STAR.cpp#L255-L272); NuSTAR calls `solo_session->finish()` before `sink->finalize()` in its single-end and paired-end runtimes.
- **Two-pass:** both discover junctions, extend their respective indexes, and remap the input. See STAR's [pass-one coordinator](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/twoPassRunPass1.cpp) and NuSTAR's [two-pass coordinator](../runtime/two_pass.cpp).

## Implementation differences

![Four implementation differences: suffix sorting, window reset, intermediate state, and BAM sorting](assets/implementation-differences.svg)

### 1. Suffix-array construction

**STAR:** count prefix ranges, divide them into chunks, collect suffix offsets for each chunk, call `qsort` with `funCompareSuffixes`, and write the sorted chunk to a temporary `SA_i` file. Read those chunks back in order and pack the suffix array. The comparator includes STAR's finite comparison length and chromosome-padding behavior.

**NuSTAR:** the default production path calls `libsais64` or its OpenMP variant to construct the initial suffix order. It removes suffix starts that are not ordinary bases, normalizes STAR-equivalent groups, builds prefix ranges, and packs the native index. Normalization still uses `std::sort` within groups that STAR treats as equal, ordered by the translated STAR coordinate. Replacing initial suffix construction therefore does not remove all comparison sorts from the builder.

```text
STAR:    prefix chunks → qsort(custom suffix comparator) → disk chunks → pack
NuSTAR:  libsais64 → filter starts → normalize STAR ranks → prefix ranges → pack
```

The distinction matters for compatibility: a plain lexicographic suffix order is not sufficient to reproduce every STAR rank choice. NuSTAR also retains a separate 32-bit builder interface; `index/sais.cpp` is not the default production path shown here. The two programs' binary index formats remain different.

Sources: [STAR chunk construction and qsort](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/Genome_genomeGenerate.cpp#L213-L325), [NuSTAR's default builder route](../index/suffix_array_index.cpp), [libsais64 call](../index/wide_suffix_array.cpp), and [STAR rank normalization](../index/star_rank.cpp).

### 2. Per-read window reset

**STAR:** `stitchPieces()` begins with two `memset` calls covering all `winBinN` entries of the strand-specific bin-to-window tables.

**NuSTAR:** the tables are allocated when the worker workspace is created. `setBin()` logs a bin when it first changes from `no_window` to an assigned window. `reset()` walks these `touched_bins`, restores them to `no_window`, and clears the log. Updating an already assigned bin does not add another log entry.

```text
STAR:    reset every bin on both strands → build this read's windows
NuSTAR:  reset logged bins → clear the log → log this read's first bin mutations
```

If the table has B entries and the previous read touched T of them, the reset visits the full table in STAR and the touched entries in NuSTAR. This changes the amount of clearing work, not the window-selection rules or the allocation size of the bin table. No isolated reset-time measurement is provided here.

Sources: [STAR full-table reset](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_stitchPieces.cpp#L12-L16), [NuSTAR workspace initialization and touched-bin reset](../window/star_window.cpp), and [existing workspace-reuse regression cases](../tests/window/workspace_reuse.cpp).

### 3. Read loop and intermediate state

**STAR:** a worker's `ReadAlignChunk` loops over reads. `ReadAlign::oneRead()` calls mapping, selection, filtering, and output methods on its `ReadAlign` object. `mapOneRead()` fills seed-related fields; `stitchPieces()` reads those fields to construct windows and stitch candidates. Arrays and members such as `PC`, `WC`, `WA`, and `trAll` connect the steps.

**NuSTAR:** the worker creates reusable stage workspaces, then calls the composed pipeline for each fragment. The candidate mapper passes a `SeedBatch` into the window builder, receives a `WindowBatch`, and passes it into the stitcher to produce a `CandidateSet`. The outer pipeline performs selection, filtering, record projection, and junction/statistic derivation, returning `FragmentEffects`. Input access, encoding, and commit order are handled by the runtime and sink interfaces.

```text
STAR:    oneRead → mapOneRead → stitchPieces → selection/filter/output
         intermediate state lives in ReadAlign fields

NuSTAR:  find → SeedBatch → build → WindowBatch → stitch → CandidateSet
         outer pipeline → FragmentEffects → sink / ordered commit
```

The read loop still exists. The recursive `stitchWindowAligns()` kernel is also retained and adapted; its source explicitly preserves branch order, score arithmetic, and candidate tie rules. The change shown here is the boundary and ownership of each stage, not a replacement of recursive stitching with a new algorithm.

Sources: [STAR per-read calls](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_oneRead.cpp#L78-L111), [STAR seed/window connection](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_mapOneRead.cpp#L34-L118), [NuSTAR candidate mapper](../pipeline/single_end/mapping_composition.cpp), [outer pipeline](../pipeline/single_end/pipeline.cpp), [ordered commit](../runtime/ordered_alignment.cpp), and [retained stitching kernel](../stitching/stitch_window.cpp).

### 4. Coordinate-sorted BAM

**STAR:** mapping writes temporary coordinate-bin files per worker. For a mapped bin, `BAMbinSortByCoordinate()` loads the bin's data, extracts coordinate/read-order keys and byte offsets into `startPos`, and uses `qsort`. It writes each sorted bin, then `bam_cat()` concatenates the bins in coordinate order. Unmapped records use a separate sorting path. The sorting coordinator checks that the largest mapped bin fits its configured memory budget.

**NuSTAR:** `bam/sink.cpp` accumulates records and calls `std::sort` before spilling an ordered run to disk. A priority queue holds the current head of each run for a k-way merge. Merge fan-in is limited to 64; more runs are combined in additional passes. If no spill was needed, finalization sorts the remaining records in memory and writes them directly.

```text
STAR:    coordinate bins → load a bin → qsort → sorted bin files → concatenate
NuSTAR:  record buffer → std::sort → sorted runs → priority-queue merge
```

NuSTAR's sort key contains reference ID, position, input-fragment ordinal, and within-fragment record order. This makes its coordinate-tie ordering explicit. The record-buffer threshold is the smaller of the requested sort budget and 8 MiB per mapping thread; this is not a total-process memory limit. The change is the sorting and temporary-I/O strategy, not just a rename from `qsort` to `std::sort`.

Sources: [STAR per-bin sorting](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/BAMbinSortByCoordinate.cpp#L12-L52), [STAR memory check and concatenation](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/bamSortByCoordinate.cpp#L14-L94), and [NuSTAR record keys, run generation, merging, and finalization](../bam/sink.cpp).

## How to interpret the figures

These are source comparisons, not new performance measurements. End-to-end benchmarks cannot isolate the speed of suffix construction, bin clearing, stitching, or BAM sorting. Use the [runtime and correctness report](STAR_COMPARISON.md) for measured workloads and the scope of the output comparisons. Solo retains STAR-derived algorithms; splitting its code into modules does not mean all of those algorithms were independently redesigned.
