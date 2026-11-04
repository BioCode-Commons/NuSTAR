# Contributing to NuSTAR

Bug reports, documentation fixes, tests, and code contributions are welcome. Start with the [README](readme.md) for supported workflows and the [code guide](CODE_STRUCTURE.md) for module responsibilities.

## Ask a question or report a bug

Open a [GitHub issue](https://github.com/BioCode-Commons/NuSTAR/issues). For a bug report, include:

- The NuSTAR version or commit, operating system, and build command.
- The complete command that failed, the relevant log or error, and the result you expected.
- A small input that reproduces the problem, or enough detail to construct one.
- For differences from STAR: the STAR version, both commands, the reference/annotation versions, and the exact records, tags, counts, or statistics that differ.
- For performance reports: CPU model and core count, RAM, storage, thread count, compiler/build flags, dataset size, and whether the index was already cached.

Share only data you have permission to publish. A generated example is useful when the original input cannot be shared.

For a new feature or a large refactor, open an issue before starting the implementation. Describe the intended behavior, affected workflows, and how the result will be compared with STAR.

## Prepare a change

1. Fork the repository and create a branch for your change.
2. Keep the change focused. Explain behavior changes separately from code cleanup.
3. Update the relevant documentation. Add a regression case when fixing a behavior bug.
4. Run the tests that cover the change and record the commands and results in your pull request.
5. Open a pull request against `main`, link the issue if there is one, and explain the problem and the resulting behavior.

Preserve upstream attribution and license notices when changing STAR-derived or bundled code. If AI tools helped produce the change, describe what they were used for and what you reviewed and tested yourself.

## Tests

The default build requires CMake 3.20+, C++20 and C99 compilers, zlib development files, and Python 3. OpenMP is optional. From the repository root:

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 8
ctest --test-dir build/release --output-on-failure
```

This runs the tests enabled by the default configuration. It does **not** enable every comparison against STAR.

### Compare two-pass, Solo, and BAM workflows with STAR

Build or obtain STAR **2.7.11b** first. Set `NUSTAR_FEATURE_STAR_BASELINE` to its executable path; replace `/absolute/path/to/STAR` below with your actual path.

```sh
cmake -S . -B build/star-comparison \
  -DCMAKE_BUILD_TYPE=Release \
  -DNUSTAR_FEATURE_STAR_BASELINE=/absolute/path/to/STAR
cmake --build build/star-comparison -j 8
ctest --test-dir build/star-comparison --output-on-failure
```

This configuration adds the junction-padding, two-pass, Solo, Solo-mode, BAM, and paired-flag comparison suites. They build their own small indexes. CTest manages their shared fixtures and ordering.

Other index and mapping-stage comparisons use separate settings, including `NUSTAR_STAR_BASELINE`, index paths, or instrumented STAR binaries. Read [tests/CMakeLists.txt](tests/CMakeLists.txt) and the relevant scripts in [tests/integration/](tests/integration/) before enabling them; some defaults refer to data outside this repository. The [oracle tools](tools/oracle/) build instrumented STAR binaries for stage comparisons.

### Choose the right comparison

| Change | Checks to include |
| :--- | :--- |
| Index construction or seed search | Index/stage comparison plus final alignment records, junctions, and statistics |
| Windows, stitching, or selection | Relevant stage comparison and single-end/paired-end output comparisons |
| Two-pass mapping | Both passes' junctions and statistics, inserted junctions, and final records |
| Solo counting or correction | Matrices, summaries, barcode/UMI tags, and relevant filtering modes |
| SAM/BAM input or output | Decoded records, tags, paired flags, and ordering where required |
| Runtime or performance | Output comparisons as well as wall time, CPU time, and peak resident memory |

Equal record counts alone do not establish equal alignments. State whether your comparison checks file bytes, decoded records in order, record multisets, or selected statistics. Explain any fields intentionally excluded.

## Report a benchmark

Use the [benchmark tools](tools/benchmark/) and consult each script's `--help` output. Save the input identifiers, reference/annotation versions, exact commands, binary hashes, machine and build details, cache conditions, and output checks alongside the measurements. Use repeated runs when making a performance claim, and report variation as well as the summary.

The [existing comparison report](docs/STAR_COMPARISON.md) explains the scope of the saved results and the missing evidence. Keep performance claims specific to the dataset and output mode you tested.
