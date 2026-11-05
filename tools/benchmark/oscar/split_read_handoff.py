#!/usr/bin/env python3
"""Split verified SAM-discrepant pairs into exact per-pair debugging inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def name(text):
    return text.split()[0].removesuffix('/1').removesuffix('/2')


def fastq(path):
    records = {}
    with path.open() as source:
        while header := source.readline():
            record = [header, source.readline(), source.readline(), source.readline()]
            if (not header.startswith('@') or not record[2].startswith('+') or not record[3]
                    or len(record[1].rstrip()) != len(record[3].rstrip())):
                raise ValueError(f'invalid FASTQ: {path}')
            key = name(header[1:])
            if key in records:
                raise ValueError(f'duplicate read name: {key}')
            records[key] = ''.join(record)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    source, destination = args.source.resolve(), args.destination.resolve()
    discrepancies = [json.loads(line) for line in (source / 'record-differences.jsonl').read_text().splitlines()]
    by_name = {row['read']: row for row in discrepancies}
    if len(by_name) != len(discrepancies):
        raise ValueError('duplicate discrepancy IDs')
    mates = {mate: fastq(source / f'R{mate}.fastq') for mate in (1, 2)}
    if mates[1].keys() != mates[2].keys() or not by_name.keys() <= mates[1].keys():
        raise ValueError('discrepancy names absent from paired FASTQ')
    headers = {engine: ''.join(line for line in (source / f'{engine}.sam').read_text().splitlines(True)
                              if line.startswith('@')) for engine in ('star', 'nustar')}
    destination.mkdir(parents=True, exist_ok=False)
    manifest = []
    for index, (key, row) in enumerate(sorted(by_name.items()), 1):
        pair = f'pair-{index:03}'
        folder = destination / 'per-read' / pair
        folder.mkdir(parents=True)
        for mate in (1, 2):
            (folder / f'R{mate}.fastq').write_text(mates[mate][key])
        summary = ['# ' + key, '', 'Observed alignments from the full-input diagnostic replay. '
                   'A replay of these 56 junction-support pairs reproduced the same records; '
                   'this individual pair has not been replayed alone.', '',
                   '| Engine | FLAG | Reference | Position | MAPQ | CIGAR | RNEXT | PNEXT | TLEN | Tags |',
                   '| --- | ---: | --- | ---: | ---: | --- | --- | ---: | ---: | --- |']
        nh = {}
        for engine in ('star', 'nustar'):
            records = row[engine]
            for line in records:
                if name(line.split('\t')[0]) != key:
                    raise ValueError('SAM record name differs from discrepancy ID')
                fields = line.rstrip().split('\t')
                summary.append('| ' + ' | '.join([engine] + fields[1:9] + [', '.join(fields[11:])]) + ' |')
            nh[engine] = sorted({tag for line in records for tag in line.rstrip().split('\t')[11:]
                                 if tag.startswith('NH:')})
            (folder / f'{engine}.sam').write_text(headers[engine] + ''.join(records))
        (folder / 'difference.json').write_text(json.dumps(row, indent=2) + '\n')
        (folder / 'README.md').write_text('\n'.join(summary) + '\n')
        manifest.append({'pair': pair, 'read_id': key, 'NH': nh})
    for mate in (1, 2):
        (destination / f'R{mate}.fastq').write_text(''.join(mates[mate][key] for key in sorted(by_name)))
    for engine in ('star', 'nustar'):
        (destination / f'{engine}.sam').write_text(headers[engine] + ''.join(
            ''.join(by_name[key][engine]) for key in sorted(by_name)))
    (destination / 'pairs.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (destination / 'pairs.tsv').write_text('pair\tread_id\tSTAR_NH\tNuSTAR_NH\n' + ''.join(
        f'{x["pair"]}\t{x["read_id"]}\t{",".join(x["NH"]["star"])}\t{",".join(x["NH"]["nustar"])}\n' for x in manifest))
    shutil.copy2(source / 'metadata.json', destination / 'metadata.json')
    shutil.copytree(source / 'context', destination / 'context')
    shutil.copy2(Path(__file__).with_name('replay_one_pair.py'), destination / 'reproduce.py')
    (destination / 'README.md').write_text('''# Individual read pairs for debugging

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
python3 reproduce.py --pair pair-001 \\
  --star /path/to/STAR --star-index /path/to/star-index \\
  --nustar /path/to/nustar --nustar-index /path/to/nustar-index \\
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
''')
    (destination / 'checksums.json').write_text(json.dumps({str(p.relative_to(destination)):
        hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(destination.rglob('*')) if p.is_file()}, indent=2) + '\n')
    with tarfile.open(destination.with_suffix('.tar.gz'), 'w:gz') as archive:
        archive.add(destination, arcname=destination.name)
    print(f'Exported {len(manifest)} exact discrepant pairs to {destination}')


if __name__ == '__main__':
    main()
