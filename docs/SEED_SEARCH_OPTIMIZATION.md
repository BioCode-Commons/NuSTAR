# Seed search optimization experiment

This change did **not** reach the 1.5x target. On the measured ARM64 machine,
complete seed search was about 1.07x faster on real yeast reads and 1.21x faster
on the synthetic workload. Native AVX2 execution and server performance remain
unverified.

## Measurements

Baseline: `3da465796b14f6f1e3f7ffb427d5b1a53fae8163`.
Platform: macOS ARM64, 10 logical CPUs, AppleClang 21, Release
`-O3 -DNDEBUG -flto=thin`. Both benchmark binaries used the same benchmark
source. The baseline used its original headers and seed library. Building the
complete baseline product on macOS also required removing an unused `omp.h`
include; this does not change runtime behavior.

| Workload | Baseline median seed time | Candidate median seed time | Speedup | Median whole-process peak RSS |
|---|---:|---:|---:|---:|
| Real yeast SE, 50,000 reads | 0.115379 s | 0.107809 s | 1.070x | 130.656 to 130.625 MiB |
| Real yeast PE, 50,000 pairs | 0.316094 s | 0.294672 s | 1.073x | 152.078 to 152.078 MiB |
| Synthetic mixed PE, 100,000 pairs | 0.731636 s | 0.602761 s | 1.214x | 203.750 to 203.828 MiB |

Each sweep visits every listed read or pair. Real-data measurements used seven
processes per binary and ten timed sweeps per process; synthetic measurements
used five processes per binary and five sweeps per process. Process order
alternated A/B and B/A. Each process warmed up before timing. Input parsing,
encoding, index loading, and full-result SHA-256 validation were outside the
timed region. The timed loop reused the production seed-result storage.

The benchmark records wall time, process CPU time, commands, input and binary
hashes, and all individual measurements. RSS includes the index, preloaded
reads, and validation. Small RSS differences should not be interpreted as
memory savings. This is a single-threaded seed benchmark, not a measurement of
complete aligner memory, 32-thread throughput, or human-genome performance.
CPU affinity and frequency were not fixed.

## Implementation

- Specialize query traversal direction at compile time.
- Compare the first two bases directly, then use bounded ARM NEON 16-byte or
  x86 AVX2 32-byte comparisons. Short tails use the existing word/scalar path.
  AVX2 is detected at runtime; no global `-mavx2` flag is required.
- Pack four canonical bases at a time when constructing SAi prefix codes.
  Noncanonical bases retain the original shift/OR behavior.
- Reuse identical uncapped queries within one split and direction when the
  split is at most 256 bases. The fixed opposite split boundary makes the
  query anchor identify its complete input. Reset validity bits for each
  direction and preserve every original `store()` call and counter update.

The query memo uses about 8 KiB of temporary stack storage per worker, or about
257 KiB of payload across 32 workers. It adds no genome-sized table, index
initialization pass, or on-disk format change. Longer and capped queries use
the exact search path. SA interval logic, tie handling, N semantics, thresholds,
and the native v4 index format remain unchanged.

## Correctness evidence

- CTest: 24 executed tests passed; one AVX2 hardware test was skipped because
  the local execution environment does not expose AVX2.
- 61 successful seed differential cases covering 337,553 fragments matched
  complete ordered result digests. These include positions, lengths,
  directions, multiplicities, closed SA intervals, splits, every counter, and
  the map marker.
- One invalid-parameter case produced the same rejection in both versions;
  it was recorded separately from successful searches.
- Cache-boundary cases cover lengths 255/256/257, longer reads, N bases,
  multiple starts, multimap filtering, and capped searches.
- Four product comparisons used 5,000 yeast reads or pairs: SE/PE at one and
  four threads. Every SAM field, tag, and record multiplicity matched, along
  with junction files and statistics excluding time-dependent fields.
- Independent byte-oracle and guard-page tests passed for NEON, prefix
  packing, and the comparison wrapper. ASan/UBSan runs of the changed search
  and benchmark, including ten cache-boundary configurations, reported no
  errors. Linked non-seed libraries were not rebuilt with instrumentation.
- x86_64 builds and Rosetta scalar-fallback comparisons passed. Actual AVX2
  kernel execution and performance still require an AVX2 host.

These are regression results for the covered inputs, not a formal proof for
all possible inputs.

## Reproduction

The Python utilities require Python 3.9 or newer. Build from the repository
root:

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 8
ctest --test-dir build/release --output-on-failure
```

On macOS, if the architecture header test cannot find standard library headers,
set `SDKROOT` to the output of `xcrun --show-sdk-path` before running CTest.

Download the bounded real dataset and construct a native NuSTAR index:

```sh
python3 tools/benchmark/download_yeast.py --out work/yeast --reads 50000
build/release/nustar --runMode genomeGenerate \
  --genomeDir work/yeast/index \
  --genomeFastaFiles work/yeast/reference.fa \
  --sjdbGTFfile work/yeast/annotation.gtf \
  --sjdbOverhang 100 --genomeSAindexNbases 10 \
  --genomeChrBinNbits 14 --runThreadN 1
```

The reads are the first 50,000 pairs from [SRR453566](https://www.ncbi.nlm.nih.gov/sra/SRR453566),
not a random sample. They are 101 bases per mate. The reference and annotation
are Ensembl release 112, *S. cerevisiae* R64-1-1, with 12,157,105 bases across
17 contigs. The CENPK sample is aligned to the standard S288C reference, as in
RNA-seq teaching exercises. The downloader saves exact URLs and SHA-256 hashes.

Build the baseline with the same compiler and flags. Copy this revision's
`tools/benchmark/seed_benchmark.cpp` into the baseline at the same path, and add
these lines to the baseline's `tools/CMakeLists.txt`:

```cmake
add_executable(seed-benchmark tools/benchmark/seed_benchmark.cpp)
target_link_libraries(seed-benchmark PRIVATE nustar_seed nustar_utils)
```

Then compare the two benchmark binaries:

```sh
python3 tools/benchmark/benchmark_seed_search.py \
  --baseline BEFORE/seed-benchmark --candidate build/release/seed-benchmark \
  --index work/yeast/index \
  --mate1 work/yeast/reads_1.fastq --mate2 work/yeast/reads_2.fastq \
  --reads 50000 --repeats 7 --sweeps 10 --output comparison.json

python3 tests/integration/seed_regression.py \
  --baseline BEFORE/seed-benchmark --candidate build/release/seed-benchmark \
  --index work/yeast/index \
  --mate1 work/yeast/reads_1.fastq --mate2 work/yeast/reads_2.fastq \
  --reads 1000 --output regression.json
```

Generate the synthetic workload with
`python3 tools/benchmark/make_seed_dataset.py work/synthetic`, then build its
index with prefix depth 10, chromosome-bin bits 18, and junction overhang 149.
Its manifest records the deterministic recipe and input hashes.

On the intended AVX2 server, first require `base_compare_avx2_test` to execute
and pass rather than skip. Run the differential suite before benchmarking, and
measure the real 32-thread mapping workload separately. The local results do
not establish a 1.5x speedup on that server.

## Other experiments

Word-only comparison, scalar-head lengths, inlining, endpoint/child prefetching,
interval narrowing, reuse of existing pieces, and interleaving independent
queries were evaluated. Interleaved variants measured roughly 1.07-1.13x on
real reads and 1.24-1.27x on some synthetic runs, with greater maintenance
complexity. The retained implementation favors the simpler query-local memo.
The proposed genome-wide two-bit auxiliary table was removed.

In one pre-change diagnostic of 5,000 pairs, about 78% of interval probes
matched fewer than two additional bases. This points toward short-probe and
memory-access overhead as areas to investigate. No hardware-counter profile
was collected to attribute the runtime precisely.
