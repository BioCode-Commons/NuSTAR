# Paired-end alignment discrepancy

This package comes from a diagnostic replay of the failed 10-million-pair K562 comparison. Both engines ran with one thread and the original indexes. The replay reproduced each engine's original splice-junction table and non-timing mapping statistics exactly.

24 read pairs have different SAM records; 56 pairs span the junctions with different support counts. The union contains 56 pairs in R1.fastq/R2.fastq. read-ids.tsv identifies both groups. These are observed differences, not a claim about their cause.

star.sam and nustar.sam contain every alignment for the exported pairs from the full replay. Optional SAM tags were sorted for comparison; all mandatory fields and tag values were retained. record-differences.jsonl contains the exact differing record groups.

The extracted-pair replay found 24 differing pairs. It reproduced the same differing read IDs: True. See subset-*.sam for the actual subset results; a reduced input is not assumed to preserve batch-dependent behavior.

metadata.json records executable SHA-256 values and full-input replay commands. subset-*.timing.json records commands for the extracted FASTQs. Change only the executable, genomeDir, input and output paths for another machine. Use the original GRCh38-2020-A reference and engine-specific indexes built with sjdbOverhang=100. The parent campaign report records source commits, reference provenance and index-build commands. Neither implementation was changed for this diagnostic run.

## Reproduce on another machine

Use STAR 2.7.11b at commit `b1edc1208d91a53bf40ebae8669f71d50b994851` and NuSTAR at `90de97ddce76b1f115f47e37d15add62784bafed`. The reference is the FASTA and GTF from 10x GRCh38-2020-A. Frozen download details are in `context/sources.json`; `context/*-index.json` records exact index-build commands, source hashes and generated file hashes. Build separate native indexes for each engine. The indexes themselves are not included.

```sh
python3 reproduce.py --star /path/to/STAR --nustar /path/to/nustar \
  --star-index /path/to/star-index --nustar-index /path/to/nustar-index \
  --output /path/to/new-output-directory
```

Five pairs change from `NH:i:1` to `NH:i:2`; 19 retain `NH:i:1` but have different records. The remaining 32 exported pairs support the affected junctions and have matching records. Each engine's subset SAM records also match its full-input records for all 56 pairs, after tag and record sorting. No root cause has been established.

The package file hashes were checked after download. `remote-checksums.json` preserves the original remote manifest; `checksums.json` covers the final handoff including added provenance. The full campaign report is included for context; its repository-relative links refer to the NuSTAR repository.
