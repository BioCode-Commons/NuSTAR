# Contributing

## Benchmark publication

Each project maintainer owns this repository's published measurements. The website publisher renders the data; it does not choose results, run benchmarks or certify output equivalence.

1. Retain commands, measured revisions, machine details, inputs, repetitions and output checks in the benchmark report. Separate index construction from sample measurements when that is the stated comparison.
2. Update `benchmarks/published.json` with the selected results. Declare units and timing scope. Update the report and any benchmark numbers repeated in the root README so they agree.
3. Pin `evidence_url` to the commit containing the report. Commit new evidence first, then reference that SHA in the data commit. Preserve the evidence commit when merging; if squashing removes it from branch history, repin the data to a retained commit before publication.
4. Open a PR in this repository. The project reviewer checks the measurements, comparison conditions and output-validation claims. The **Published benchmark** action checks the data format; a green action alone does not establish scientific correctness.
5. Merge after review and a passing check. The action notifies the website, which fetches the data and deploys the benchmark tables. Check **Publish benchmark tables** on the [website Actions page](https://github.com/BioCode-Commons/BioCode-Commons.github.io/actions/workflows/benchmarks.yml) for deployment status.

Only benchmark data is synchronized. Editing a README alone does not publish new numbers. Project descriptions, validation prose and other website content are maintained separately.

### Failed publication or corrections

- Validation failed: fix the data here and rerun the check. Do not relax the shared validator to admit incomplete measurements.
- Notification failed: contact the organization administrator. The hourly website check remains a fallback; GitHub may delay it.
- Website publication failed: the live site keeps its previous results. The workflow log identifies the source or field needing correction. Fix that source and rerun publication.
- Incorrect published result: revert the relevant `published.json` change through review. The next successful publication uses the restored data.

### Responsibility and handover

Project maintainers and reviewers own benchmark selection and evidence in their repository. The organization administrator maintains the shared publisher and rotates the organization GitHub App private key. Credentials belong in Actions Secrets, never in source files, reports, issues or pull requests. A new project maintainer needs neither a personal publishing token nor access to the website repository. The workflow contains no individual maintainer names.
