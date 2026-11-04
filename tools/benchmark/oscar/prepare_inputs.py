#!/usr/bin/env python3
"""Fetch pinned inputs into the campaign root and create synchronized FASTQs."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import random
import tarfile
from urllib.request import Request, urlopen

from common import Capacity, sha256, write_json


def download(resource, folder, capacity):
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / resource['filename']
    receipt = target.with_suffix(target.suffix + '.json')
    identity = hashlib.sha256(json.dumps(resource, sort_keys=True).encode()).hexdigest()
    if target.exists():
        if not receipt.exists():
            raise RuntimeError(f'untracked existing input: {target}')
        saved = json.loads(receipt.read_text())
        if saved['resource_identity'] != identity or sha256(target) != saved['sha256']:
            raise RuntimeError(f'existing input failed identity check: {target}')
        return target
    part = target.with_suffix(target.suffix + '.partial')
    partial_receipt = part.with_suffix(part.suffix + '.json')
    if part.exists():
        if not partial_receipt.exists() or json.loads(partial_receipt.read_text())['resource_identity'] != identity:
            raise RuntimeError('partial download belongs to an unknown or different resource')
    else:
        write_json(partial_receipt, {'resource_identity': identity})
    start = part.stat().st_size if part.exists() else 0
    if start > resource['bytes']:
        raise RuntimeError(f'partial download exceeds declared size: {part}')
    base = resource.get('byte_range', [0, resource['bytes'] - 1])[0]
    ranged = 'byte_range' in resource or start > 0
    headers = {}
    if ranged:
        headers['Range'] = f'bytes={base + start}-{base + resource["bytes"] - 1}'
    if resource.get('etag'):
        headers['If-Match'] = resource['etag']
    if start < resource['bytes']:
        with urlopen(Request(resource['url'], headers=headers), timeout=120) as remote:
            expected = f'bytes {base + start}-{base + resource["bytes"] - 1}/'
            if ranged and (remote.status != 206 or not remote.headers.get('Content-Range', '').startswith(expected)):
                raise RuntimeError('download server did not honor the exact byte range')
            if not ranged and remote.status != 200:
                raise RuntimeError('whole-object download did not return HTTP 200')
            if resource.get('etag') and remote.headers.get('ETag') != resource['etag']:
                raise RuntimeError('source archive identity changed')
            with part.open('ab') as output:
                for block in iter(lambda: remote.read(8 << 20), b''):
                    capacity.check()
                    if output.tell() + len(block) > resource['bytes']:
                        raise RuntimeError('download exceeds the declared object size')
                    output.write(block)
    if part.stat().st_size != resource['bytes']:
        raise RuntimeError(f'incomplete download: {part}')
    digest = sha256(part)
    if resource.get('sha256') and digest != resource['sha256']:
        raise RuntimeError(f'SHA256 mismatch: {part}')
    if resource.get('md5'):
        md5 = hashlib.md5()
        with part.open('rb') as source:
            for block in iter(lambda: source.read(8 << 20), b''):
                md5.update(block)
        if md5.hexdigest() != resource['md5']:
            raise RuntimeError(f'published MD5 mismatch: {part}')
    part.replace(target)
    partial_receipt.unlink()
    write_json(receipt, {'resource_identity': identity, 'sha256': digest,
                         'bytes': target.stat().st_size, 'source': resource})
    return target


def atomic_copy(source, target, capacity):
    part = target.with_suffix(target.suffix + '.partial')
    digest = hashlib.sha256()
    with part.open('wb') as output:
        for block in iter(lambda: source.read(8 << 20), b''):
            capacity.check()
            output.write(block)
            digest.update(block)
    part.replace(target)
    return {'sha256': digest.hexdigest(), 'bytes': target.stat().st_size}


def reference(resources, root, capacity):
    folder = root / 'reference'
    folder.mkdir(exist_ok=True)
    done = folder / 'prepared.json'
    identity = hashlib.sha256(json.dumps({key: resources[key] for key in ['reference', 'whitelist']},
                                        sort_keys=True).encode()).hexdigest()
    if done.exists():
        saved = json.loads(done.read_text())
        if saved['resource_identity'] != identity:
            raise RuntimeError('prepared reference uses different sources')
        for name, info in saved['files'].items():
            if sha256(folder / name) != info['sha256']:
                raise RuntimeError('prepared reference file changed: ' + name)
        return
    archive = download(resources['reference'], root / 'downloads/reference', capacity)
    wanted = resources['reference']['members']
    found = {}
    # Only these exact regular files are extracted. The packaged STAR index is discarded.
    with tarfile.open(archive, 'r|gz') as tar:
        for member in tar:
            if member.name in wanted:
                if not member.isfile():
                    raise RuntimeError('reference member is not a regular file')
                with tar.extractfile(member) as stream:
                    found[wanted[member.name]] = atomic_copy(stream, folder / wanted[member.name], capacity)
    if set(found) != set(wanted.values()):
        raise RuntimeError(f'reference members missing: {set(wanted.values()) - set(found)}')
    whitelist = download(resources['whitelist'], root / 'downloads/reference', capacity)
    with gzip.open(whitelist, 'rb') as source:
        found['3M-february-2018.txt'] = atomic_copy(source, folder / '3M-february-2018.txt', capacity)
    count = 0
    with (folder / '3M-february-2018.txt').open('rb') as source:
        for count, line in enumerate(source, 1):
            value = line.rstrip(b'\r\n')
            if len(value) != 16 or not set(value) <= set(b'ACGT'):
                raise RuntimeError('invalid v3 barcode whitelist')
    if count != resources['whitelist']['barcode_count']:
        raise RuntimeError('barcode whitelist count changed')
    write_json(done, {'resource_identity': identity, 'files': found,
                      'reference_metadata': json.loads((folder / 'reference.json').read_text())})
    # Archive hash and source are retained in its receipt; never materialize its prebuilt index.
    archive.unlink()


def fastq_record(source):
    header = source.readline()
    if not header:
        return None
    sequence, plus, quality = source.readline(), source.readline(), source.readline()
    if (not header.startswith(b'@') or not plus.startswith(b'+') or not quality.endswith(b'\n')
            or len(sequence.rstrip(b'\r\n')) != len(quality.rstrip(b'\r\n'))):
        raise ValueError('truncated or invalid four-line FASTQ record')
    return header, sequence, plus, quality


def read_name(record):
    token = record[0].split()[0]
    return token.removesuffix(b'/1').removesuffix(b'/2')


def pairs(files):
    for r1, r2 in files:
        with gzip.open(r1, 'rb') as left, gzip.open(r2, 'rb') as right:
            while True:
                a, b = fastq_record(left), fastq_record(right)
                if a is None and b is None:
                    break
                if a is None or b is None or read_name(a) != read_name(b):
                    raise ValueError(f'FASTQ mates are not synchronized: {r1}, {r2}')
                yield a, b


def prepare_dataset(name, resources, root, capacity, full):
    selected = {key: value for key, value in resources.items() if key.startswith(name + '_')}
    raw = {key: download(value, root / 'downloads' / name, capacity) for key, value in selected.items()}
    lanes = sorted({value.get('lane', '') for value in selected.values()})
    files = []
    for lane in lanes:
        mate_paths = {value['mate']: raw[key] for key, value in selected.items() if value.get('lane', '') == lane}
        if set(mate_paths) != {'R1', 'R2'}:
            raise ValueError('incomplete lane mate set')
        files.append((mate_paths['R1'], mate_paths['R2']))
    folder = root / 'prepared' / name
    folder.mkdir(parents=True, exist_ok=True)
    identity = {path.name: sha256(path) for pair in files for path in pair}
    scan_path = folder / 'scan.json'
    if scan_path.exists():
        scan = json.loads(scan_path.read_text())
        if scan['compressed_sha256'] != identity:
            raise RuntimeError('source FASTQs changed since their scan')
    else:
        count, sizes = 0, [0, 0]
        for count, record_pair in enumerate(pairs(files), 1):
            for mate, record in enumerate(record_pair):
                sizes[mate] += sum(map(len, record))
        scan = {'pairs': count, 'plaintext_bytes': sizes, 'compressed_sha256': identity,
                'lane_order': [[str(x) for x in pair] for pair in files]}
        if count < 10_000_000:
            raise RuntimeError('selected library has fewer than ten million pairs')
        write_json(scan_path, scan)
    stage_names = ['full'] if full else ['pilot', 'main']
    if all((folder / stage / 'prepared.json').exists() for stage in stage_names):
        for stage in stage_names:
            saved = json.loads((folder / stage / 'prepared.json').read_text())
            if saved['compressed_sha256'] != identity:
                raise RuntimeError('prepared subset provenance differs')
            for mate, info in saved['files'].items():
                if sha256(folder / stage / mate) != info['sha256']:
                    raise RuntimeError('prepared FASTQ checksum changed')
        return
    selections = {}
    if full:
        required = sum(scan['plaintext_bytes'])
        used = capacity.check()
        if required > min(capacity.budget - used['campaign_allocated_bytes'],
                          used['filesystem_free_bytes'] - capacity.reserve):
            raise RuntimeError('insufficient declared capacity for the full uncompressed library')
    else:
        sample = random.Random(20261104).sample(range(scan['pairs']), 10_000_000)
        selections = {'pilot': set(sample[:1_000_000]), 'main': set(sample)}
    outputs, digests, counts = {}, {}, {}
    try:
        for stage in stage_names:
            destination = folder / stage
            destination.mkdir(exist_ok=True)
            outputs[stage] = [(destination / f'R{mate}.fastq.partial').open('wb') for mate in [1, 2]]
            digests[stage] = [hashlib.sha256(), hashlib.sha256()]
            counts[stage] = 0
        for ordinal, record_pair in enumerate(pairs(files)):
            for stage in stage_names:
                if full or ordinal in selections[stage]:
                    for mate, record in enumerate(record_pair):
                        data = b''.join(record)
                        outputs[stage][mate].write(data)
                        digests[stage][mate].update(data)
                    counts[stage] += 1
            if ordinal % 100_000 == 0:
                capacity.check()
    finally:
        for streams in outputs.values():
            for stream in streams:
                stream.close()
    for stage in stage_names:
        expected = scan['pairs'] if full else len(selections[stage])
        if counts[stage] != expected:
            raise RuntimeError('prepared subset count differs from selection')
        files_meta = {}
        for mate in [1, 2]:
            target = folder / stage / f'R{mate}.fastq'
            target.with_suffix('.fastq.partial').replace(target)
            files_meta[target.name] = {'sha256': digests[stage][mate - 1].hexdigest(), 'bytes': target.stat().st_size}
        write_json(folder / stage / 'prepared.json', {'pairs': counts[stage], 'files': files_meta,
                   'compressed_sha256': identity, 'seed': None if full else 20261104,
                   'sampling': 'all records' if full else 'uniform without replacement; pilot nested in main; source order retained'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--dataset', choices=['reference', 'bulk', 'pbmc'], required=True)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--disk-budget-gib', type=int, required=True)
    parser.add_argument('--reserve-gib', type=int, default=32)
    args = parser.parse_args()
    root = args.root.resolve()
    if not (root / 'campaign.json').is_file():
        parser.error('root must be an initialized campaign directory')
    capacity = Capacity(root, args.disk_budget_gib << 30, args.reserve_gib << 30)
    resources = json.loads(args.sources.read_text())['resources']
    if args.dataset == 'reference':
        reference(resources, root, capacity)
    else:
        prepare_dataset(args.dataset, resources, root, capacity, args.full)


if __name__ == '__main__':
    main()
