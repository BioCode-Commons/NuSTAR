# Published benchmarks

Each project maintainer updates `benchmarks/published.json` and the corresponding benchmark documentation in their own repository, then GitHub Actions validates the data and publishes only the benchmark updates to the BioCode Commons website after the changes merge into `main`.

- Update this file when publishing a reviewed benchmark result. Keep experimental runs out of it until maintainers select them for publication.
- Record measured baseline and candidate wall times and peak RSS, with explicit units, resource counts, memory scope and timing scope. Do not infer missing measurements or select the fastest repeat.
- Link `evidence_url` to benchmark evidence at a full commit SHA. The evidence must identify the measured code, inputs, machine, repetitions and output checks. A documentation commit is not necessarily the measured binary revision.
- Keep captions and notes accurate. Report index construction separately when excluded from sample runtime and RSS. Do not silently change comparison conditions or claim validation beyond the tests performed.
- Preserve the JSON structure. Numbers are positive numeric values, seconds for time and the declared MiB/GiB unit for RSS. Each dataset/resource pair is unique. The website calculates speedup from the two times.
- Update the corresponding benchmark documentation in the same PR. Run the Published benchmark workflow; fix validation errors before merging to `main`.
- Do not edit the website manually to publish these results. Only benchmark tables and performance summaries sync; project descriptions and validation prose remain under website maintainers' control.

See `benchmarks/README.md` for the publication workflow. Do not add heuristic placeholders or data-specific exceptions to make invalid benchmark records pass.

## Maintainer workflow

Read root `CONTRIBUTING.md` before publishing benchmark changes. Keep numbers repeated in the root README consistent with `benchmarks/published.json` and its evidence. Prose edits alone do not update the website. Do not configure a personal token per maintainer or change another project's results. The organization administrator owns the publishing credential; project reviewers own measurement quality.
