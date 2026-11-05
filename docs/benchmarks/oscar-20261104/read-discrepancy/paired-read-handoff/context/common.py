"""File integrity, capacity, and output checks shared by the Oscar campaign."""
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(8 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def allocated_bytes(root):
    total = 0
    for parent, dirs, files in os.walk(root, followlinks=False):
        for name in files:
            path = Path(parent) / name
            try:
                total += path.lstat().st_blocks * 512
            except FileNotFoundError:
                pass  # Sort writers can remove their own completed temporary runs.
    return total


class Capacity:
    def __init__(self, root, budget, reserve):
        self.root = Path(root).resolve()
        self.budget = int(budget)
        self.reserve = int(reserve)
        if self.budget <= 0 or self.reserve < 0:
            raise ValueError('capacity budget must be positive; reserve must be nonnegative')
        self.peak = 0

    def check(self):
        used = allocated_bytes(self.root)
        free = shutil.disk_usage(self.root).free
        self.peak = max(self.peak, used)
        if used > self.budget:
            raise RuntimeError(f'campaign disk budget exceeded: {used} > {self.budget}')
        if free < self.reserve:
            raise RuntimeError(f'filesystem free-space reserve reached: {free} < {self.reserve}')
        return {'campaign_allocated_bytes': used, 'filesystem_free_bytes': free,
                'campaign_peak_observed_bytes': self.peak}


def same_file(left, right):
    with Path(left).open('rb') as a, Path(right).open('rb') as b:
        digest = hashlib.sha256()
        blocks = 0
        for blocks, (x, y) in enumerate(itertools.zip_longest(
                iter(lambda: a.read(8 << 20), b''), iter(lambda: b.read(8 << 20), b'')), 1):
            if x != y:
                raise AssertionError(f'file differs at block {blocks}: {left}, {right}')
            digest.update(x)
    return {'sha256': digest.hexdigest(), 'bytes': Path(left).stat().st_size}


def mapping_statistics(path):
    changing = {'Started job on', 'Started mapping on', 'Finished on',
                'Mapping speed, Million of reads per hour'}
    rows = {}
    for line in Path(path).read_text().splitlines():
        if '|' in line:
            key, value = line.split('|', 1)
            key = key.strip()
            if key not in changing:
                rows[key] = value.strip()
    if not rows or 'Number of input reads' not in rows:
        raise ValueError(f'no mapping statistics: {path}')
    return rows


def compare_sam_streams(left, right):
    """Every record and @SQ line, in order; all other headers are excluded."""
    def records(stream, header):
        for line in stream:
            if line.startswith(b'@'):
                if line.startswith(b'@SQ\t'):
                    header.append(line.rstrip(b'\r\n'))
                continue
            if line.count(b'\t') < 10:
                raise ValueError('decoded SAM record has fewer than 11 fields')
            yield line.rstrip(b'\r\n')
    headers = ([], [])
    digest = hashlib.sha256()
    count = 0
    for count, (a, b) in enumerate(itertools.zip_longest(
            records(left, headers[0]), records(right, headers[1])), 1):
        if a != b:
            raise AssertionError(f'decoded BAM differs at ordered record {count}: {a!r} / {b!r}')
        digest.update(a + b'\n')
    if headers[0] != headers[1]:
        raise AssertionError('BAM sequence dictionaries differ')
    return {'records': count, 'ordered_sam_sha256': digest.hexdigest(),
            'sequence_dictionary_lines': len(headers[0]),
            'excluded': 'headers other than @SQ; compressed BAM bytes'}


def compare_bam(left, right, samtools, log_dir):
    processes = []
    logs = []
    try:
        for i, source in enumerate((left, right)):
            log = (Path(log_dir) / f'bam-decode-{i}.stderr').open('wb')
            logs.append(log)
            processes.append(subprocess.Popen([str(samtools), 'view', '-h', str(source)],
                                               stdout=subprocess.PIPE, stderr=log))
        result = compare_sam_streams(processes[0].stdout, processes[1].stdout)
        for process in processes:
            process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError('samtools decoding failed; see decoder stderr')
        return result
    finally:
        for process in processes:
            if process.stdout:
                process.stdout.close()
            if process.poll() is None:
                process.terminate()
            process.wait()
        for log in logs:
            log.close()


def validate_pair(star, nustar, case, samtools, log_dir):
    result = {'files': {}}
    suffixes = ['SJ.out.tab']
    if case['two_pass']:
        suffixes += ['_STARpass1/SJ.out.tab', '_STARgenome/sjdbInfo.txt',
                     '_STARgenome/sjdbList.out.tab']
    if case['solo']:
        folders = (Path(star + 'Solo.out'), Path(nustar + 'Solo.out'))
        sets = [{p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file()}
                for folder in folders]
        if not sets[0] or sets[0] != sets[1]:
            raise AssertionError('Solo output file sets differ or are empty')
        suffixes += ['Solo.out/' + name for name in sorted(sets[0])]
    for suffix in suffixes:
        result['files'][suffix] = same_file(star + suffix, nustar + suffix)
    stats = [mapping_statistics(prefix + 'Log.final.out') for prefix in (star, nustar)]
    if stats[0] != stats[1]:
        raise AssertionError(f'non-time mapping statistics differ: {stats}')
    result['statistics'] = stats[0]
    if case['bam']:
        result['bam'] = compare_bam(star + 'Aligned.sortedByCoord.out.bam',
                                    nustar + 'Aligned.sortedByCoord.out.bam', samtools, log_dir)
    return result
