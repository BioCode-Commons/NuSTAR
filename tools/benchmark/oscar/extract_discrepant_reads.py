#!/usr/bin/env python3
"""Replay the failed paired-end comparison and export read-level evidence."""
import argparse
from collections import Counter
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

from common import Capacity, mapping_statistics, sha256, write_json
from run_suite import measure


def read_name(value):
    name = value.split()[0]
    return re.sub(r'/[12]$', '', name)


def junctions(fields):
    position = int(fields[3])
    for length, operation in re.findall(r'(\d+)([MIDNSHP=X])', fields[5]):
        length = int(length)
        if operation == 'N':
            yield fields[2], str(position), str(position + length - 1)
        if operation in 'MDN=X':
            position += length


def normalize_sam(source, destination, changed_junctions, capacity):
    selected = set()
    headers = []
    with source.open() as incoming, destination.open('w') as outgoing:
        for number, line in enumerate(incoming):
            if line.startswith('@'):
                headers.append(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                raise ValueError('incomplete SAM record')
            fields[0] = read_name(fields[0])
            if changed_junctions.intersection(junctions(fields)):
                selected.add(fields[0])
            outgoing.write('\t'.join(fields[:11] + sorted(fields[11:])) + '\n')
            if number % 100000 == 0:
                capacity.check()
    return selected, headers


def groups(path):
    with path.open() as source:
        for name, records in itertools.groupby(source, lambda line: line.split('\t', 1)[0]):
            yield name, tuple(records)


def differences(left, right):
    a, b = iter(groups(left)), iter(groups(right))
    x, y = next(a, None), next(b, None)
    while x is not None or y is not None:
        name = min(item[0] for item in (x, y) if item is not None)
        lx = x[1] if x is not None and x[0] == name else ()
        ry = y[1] if y is not None and y[0] == name else ()
        if lx != ry:
            yield name, lx, ry
        if x is not None and x[0] == name:
            x = next(a, None)
        if y is not None and y[0] == name:
            y = next(b, None)


def fastq_record(source):
    first = source.readline()
    if not first:
        return None
    lines = [first, source.readline(), source.readline(), source.readline()]
    if (not first.startswith('@') or not lines[2].startswith('+')
            or not lines[3] or len(lines[1].rstrip()) != len(lines[3].rstrip())):
        raise ValueError('invalid FASTQ record')
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or os.uname().sysname != 'Linux':
        parser.error('run inside a Linux Slurm allocation')
    root = args.root.resolve()
    work = root / 'diagnostics' / ('paired-reads-job-' + os.environ['SLURM_JOB_ID'])
    work.mkdir(parents=True, exist_ok=False)
    capacity = Capacity(root, 350 << 30, 32 << 30)
    original = root / 'results/main-job-8361932/bulk-pe-1-r1'
    package = work / 'handoff'
    package.mkdir()
    state = {'state': 'running', 'step': 'replay', 'job_id': os.environ['SLURM_JOB_ID']}
    write_json(work / 'status.json', state)
    try:
        metadata = {'original_comparison': str(original), 'job_id': os.environ['SLURM_JOB_ID'],
                    'hostname': os.uname().nodename, 'commands': {}, 'files': {}}
        changed = set()
        original_sj = {}
        for engine in ('star', 'nustar'):
            original_sj[engine] = {tuple(row.split('\t')[:6]): row for row in
                                     (original / engine / 'output.SJ.out.tab').read_text().splitlines()}
        for key in original_sj['star'].keys() | original_sj['nustar'].keys():
            if original_sj['star'].get(key) != original_sj['nustar'].get(key):
                changed.add(key[:3])
        metadata['changed_junctions'] = sorted(changed)
        headers, selected = {}, set()
        for engine in ('star', 'nustar'):
            destination = work / engine
            destination.mkdir()
            original_command = json.loads((original / engine / 'timing.json').read_text())['command']
            command = list(original_command)
            command[command.index('--outSAMtype') + 1] = 'SAM'
            command[command.index('--outFileNamePrefix') + 1] = str(destination / 'output.')
            command += ['--outSAMunmapped', 'Within', '--outSAMattributes', 'NH', 'HI', 'AS', 'nM']
            metadata['commands'][engine] = command
            metadata['files'][command[0]] = sha256(command[0])
            measure(command, destination, capacity)
            if sha256(destination / 'output.SJ.out.tab') != sha256(original / engine / 'output.SJ.out.tab'):
                raise AssertionError(f'{engine}: diagnostic SJ output differs from original failure')
            if mapping_statistics(destination / 'output.Log.final.out') != mapping_statistics(original / engine / 'output.Log.final.out'):
                raise AssertionError(f'{engine}: diagnostic mapping statistics differ from original failure')
            hits, headers[engine] = normalize_sam(destination / 'output.Aligned.out.sam',
                                                 destination / 'canonical.sam', changed, capacity)
            selected.update(hits)
            sort_work = destination / 'sort'
            sort_work.mkdir()
            measure(['sort', '--parallel=1', '-S', '2G', '-T', str(sort_work),
                     '-o', str(destination / 'canonical.sorted.sam'), str(destination / 'canonical.sam')],
                    sort_work, capacity)
        if [x for x in headers['star'] if x.startswith('@SQ\t')] != [x for x in headers['nustar'] if x.startswith('@SQ\t')]:
            raise AssertionError('reference dictionaries differ')
        state['step'] = 'extract'
        write_json(work / 'status.json', state)
        differing = set()
        with (package / 'record-differences.jsonl').open('w') as output:
            for name, left, right in differences(work / 'star/canonical.sorted.sam', work / 'nustar/canonical.sorted.sam'):
                differing.add(name)
                output.write(json.dumps({'read': name, 'star': left, 'nustar': right}) + '\n')
        junction_reads = selected.copy()
        selected.update(differing)
        if not differing:
            raise AssertionError('no SAM differences despite reproduced summary discrepancy')
        for engine in ('star', 'nustar'):
            with (package / f'{engine}.sam').open('w') as output:
                output.writelines(headers[engine])
                for name, records in groups(work / engine / 'canonical.sorted.sam'):
                    if name in selected:
                        output.writelines(records)
            for filename in ('output.SJ.out.tab', 'output.Log.final.out', 'timing.json'):
                shutil.copy2(work / engine / filename, package / f'{engine}.{filename}')
        found = Counter()
        inputs = [root / 'prepared/bulk/main' / f'R{mate}.fastq' for mate in (1, 2)]
        with inputs[0].open() as r1, inputs[1].open() as r2, (package / 'R1.fastq').open('w') as o1, (package / 'R2.fastq').open('w') as o2:
            count = 0
            while True:
                left, right = fastq_record(r1), fastq_record(r2)
                if left is None and right is None:
                    break
                if left is None or right is None or read_name(left[0][1:]) != read_name(right[0][1:]):
                    raise ValueError('FASTQ mate mismatch')
                name = read_name(left[0][1:])
                if name in selected:
                    found[name] += 1
                    o1.writelines(left)
                    o2.writelines(right)
                count += 1
                if count % 100000 == 0:
                    capacity.check()
        if set(found) != selected or any(value != 1 for value in found.values()):
            raise AssertionError('selected names missing or duplicated in original FASTQ')
        metadata.update({'input_pairs_scanned': count, 'differing_read_pairs': len(differing),
                         'junction_support_read_pairs': len(junction_reads), 'exported_read_pairs': len(selected),
                         'original_summary_reproduced': True})
        with (package / 'read-ids.tsv').open('w') as output:
            output.write('read_id\tSAM_records_differ\tspans_changed_junction\n')
            for name in sorted(selected):
                output.write(f'{name}\t{int(name in differing)}\t{int(name in junction_reads)}\n')
        # Test whether the extracted pairs alone reproduce the full-input record differences.
        state['step'] = 'subset-replay'
        write_json(work / 'status.json', state)
        for engine in ('star', 'nustar'):
            destination = work / ('subset-' + engine)
            destination.mkdir()
            command = list(metadata['commands'][engine])
            offset = command.index('--readFilesIn')
            command[offset + 1:offset + 3] = [str(package / 'R1.fastq'), str(package / 'R2.fastq')]
            command[command.index('--outFileNamePrefix') + 1] = str(destination / 'output.')
            measure(command, destination, capacity)
            normalize_sam(destination / 'output.Aligned.out.sam', destination / 'canonical.sam', changed, capacity)
            subprocess.run(['sort', '-S', '2G', '-T', str(destination), '-o', str(destination / 'canonical.sorted.sam'),
                            str(destination / 'canonical.sam')], check=True)
            shutil.copy2(destination / 'output.Aligned.out.sam', package / f'subset-{engine}.sam')
            shutil.copy2(destination / 'timing.json', package / f'subset-{engine}.timing.json')
        subset_differences = list(differences(work / 'subset-star/canonical.sorted.sam', work / 'subset-nustar/canonical.sorted.sam'))
        metadata['subset_differing_read_pairs'] = len(subset_differences)
        metadata['subset_reproduces_same_differing_read_ids'] = {x[0] for x in subset_differences} == differing
        write_json(package / 'metadata.json', metadata)
        (package / 'reproduce.py').write_text('''#!/usr/bin/env python3
"""Replay the extracted pairs with existing engine-specific reference indexes."""
import argparse, json, pathlib, subprocess
p = argparse.ArgumentParser(description=__doc__)
for option in ('star', 'nustar', 'star-index', 'nustar-index', 'output'):
    p.add_argument('--' + option, type=pathlib.Path, required=True)
a = p.parse_args()
here = pathlib.Path(__file__).resolve().parent
metadata = json.loads((here / 'metadata.json').read_text())
for engine in ('star', 'nustar'):
    destination = a.output.resolve() / engine
    destination.mkdir(parents=True, exist_ok=False)
    command = list(metadata['commands'][engine])
    command[0] = str(getattr(a, engine).resolve())
    command[command.index('--genomeDir') + 1] = str(getattr(a, engine + '_index').resolve())
    i = command.index('--readFilesIn')
    command[i + 1:i + 3] = [str(here / 'R1.fastq'), str(here / 'R2.fastq')]
    command[command.index('--outFileNamePrefix') + 1] = str(destination / 'output.')
    print(command, flush=True)
    subprocess.run(command, check=True)
''')
        (package / 'README.md').write_text(
            '# Paired-end alignment discrepancy\n\n'
            'This package comes from a diagnostic replay of the failed 10-million-pair K562 comparison. '
            'Both engines ran with one thread and the original indexes. The replay reproduced each engine\'s '
            'original splice-junction table and non-timing mapping statistics exactly.\n\n'
            f'{len(differing)} read pairs have different SAM records; {len(junction_reads)} pairs span the '
            f'junctions with different support counts. The union contains {len(selected)} pairs in R1.fastq/R2.fastq. '
            'read-ids.tsv identifies both groups. These are observed differences, not a claim about their cause.\n\n'
            'star.sam and nustar.sam contain every alignment for the exported pairs from the full replay. '
            'Optional SAM tags were sorted for comparison; all mandatory fields and tag values were retained. '
            'record-differences.jsonl contains the exact differing record groups.\n\n'
            f'The extracted-pair replay found {len(subset_differences)} differing pairs. '
            f'It reproduced the same differing read IDs: {metadata["subset_reproduces_same_differing_read_ids"]}. '
            'See subset-*.sam for the actual subset results; a reduced input is not assumed to preserve batch-dependent behavior.\n\n'
            'metadata.json records executable SHA-256 values and full-input replay commands. '
            'subset-*.timing.json records commands for the extracted FASTQs. Change only the executable, '
            'genomeDir, input and output paths for another machine. Use the original GRCh38-2020-A reference '
            'and engine-specific indexes built with sjdbOverhang=100. '
            'The parent campaign report records source commits, reference provenance and index-build commands. '
            'Neither implementation was changed for this diagnostic run.\n')
        write_json(package / 'checksums.json', {p.name: sha256(p) for p in package.iterdir() if p.is_file()})
        archive = work / 'paired-read-handoff.tar.gz'
        with tarfile.open(archive, 'w:gz') as output:
            output.add(package, arcname='paired-read-handoff')
        state.update({'state': 'complete', 'step': 'finished', 'archive': str(archive), 'summary': metadata})
        write_json(work / 'status.json', state)
    except BaseException as error:
        state.update({'state': 'failed', 'error': repr(error)})
        write_json(work / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
