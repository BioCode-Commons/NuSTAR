# Publishing benchmarks

Edit `published.json` and the supporting benchmark report in one PR. Keep the evidence URL pinned to the commit containing that report; if the report is new, commit it first and use that commit SHA. Merging the data file into `main` publishes the selected results to the BioCode Commons website.

`rows` contains the displayed comparisons: dataset, resources, baseline/candidate seconds and baseline/candidate RSS. `rss_unit`, `rss_scope` and `time_scope` define the measurement. `caption` and `note` are arrays of plain-text segments with optional HTTPS GitHub links. `display` controls decimal places and the speedup column. Speedups are calculated by the website, never entered manually.

The action validates against the website's shared validator and notifies the website after a main-branch data update. `BENCHMARK_PUBLISH_TOKEN` is configured in this repository's Actions Secrets and maintained by the organization administrator. Project maintainers do not need to create, copy or rotate their own token. The website also checks hourly as a fallback; GitHub may delay scheduled jobs.

The website reads all four current main branches, validates the complete snapshot and deploys only after every record passes. A failed fetch, validation or build leaves the previous site live. No benchmark program is executed by this publication workflow.
