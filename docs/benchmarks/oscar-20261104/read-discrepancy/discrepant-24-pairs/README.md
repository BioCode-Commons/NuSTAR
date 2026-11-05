# Individual read pairs for debugging

This package contains only the 24 read pairs whose SAM records differ between
STAR and NuSTAR. It excludes the 32 additional matching junction-support pairs
included in the earlier 56-pair handoff. FASTQ sequences, quality strings, and
read names are copied exactly from the original paired inputs.

`per-read/pair-001/` through `per-read/pair-024/` each contain one R1 record,
one R2 record, both engines' original SAM records (with headers), exact JSON
differences, and a short comparison table. `pairs.tsv` maps directories to
original read IDs. The top-level FASTQs/SAMs contain all 24 pairs together.
Keep R1 and R2 together when debugging: the discrepancy was found in paired-end
alignment. A single mate alone exercises a different workload.

Five pairs change from NH:i:1 to NH:i:2; 19 retain NH:i:1 but their records differ.
The earlier 56-pair input reproduced all 24 differences and the same records
as the full 10-million-pair replay. Neither this 24-pair-only input nor each
individual pair has yet been independently replayed; batch dependence remains
possible. The supplied SAM files are observed results, not predictions of a
one-pair replay. No implementation root cause is claimed.

## Run one pair

Use the original engine versions and full reference indexes recorded in
`metadata.json` and `context/*-index.json`. STAR: 2.7.11b commit
`b1edc1208d91a53bf40ebae8669f71d50b994851`; NuSTAR:
`90de97ddce76b1f115f47e37d15add62784bafed`. Reference: 10x GRCh38-2020-A,
including the original GTF junctions. Index data and executables are not bundled.

```sh
python3 reproduce.py --pair pair-001 \
  --star /path/to/STAR --star-index /path/to/star-index \
  --nustar /path/to/nustar --nustar-index /path/to/nustar-index \
  --output /path/to/new-output-directory
```

Add `--engine nustar --dry-run` to print a shell-quoted command without running
the aligner or creating output directories. Paste the printed command after
`gdb --args` (Linux) or `lldb --` to step through that engine. For source-level
debugging, build the same source revision with debug symbols; the recorded
binary hashes identify the original benchmark builds. Without `--pair`, the
script uses all 24 pairs. All alignment flags remain those of the diagnostic
replay; input/output and executable/index paths are replaced explicitly.

`checksums.json` covers all package files except itself. This is a debugging
fixture and is not a passing benchmark or a minimized causal reproducer.
