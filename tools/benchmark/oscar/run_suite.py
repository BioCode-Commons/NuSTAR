#!/usr/bin/env python3
"""Run and validate a serial STAR/NuSTAR matrix inside one Slurm allocation."""
import argparse
import itertools
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import threading
import time

from common import Capacity, sha256, validate_pair, write_json


def measure(command, work, capacity):
    stopped = threading.Event()
    failures = []
    with (work / 'console.log').open('wb') as log:
        capacity.check()
        start = time.perf_counter()
        process = subprocess.Popen(command, cwd=work, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        def watch():
            while not stopped.wait(2):
                try:
                    capacity.check()
                except Exception as error:
                    failures.append(str(error))
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    if not stopped.wait(15):
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    return
        watcher = threading.Thread(target=watch, daemon=True)
        watcher.start()
        try:
            _, status, usage = os.wait4(process.pid, 0)
            seconds = time.perf_counter() - start
            process.returncode = os.waitstatus_to_exitcode(status)
        finally:
            stopped.set()
            watcher.join()
    # ru_maxrss is KiB on the Linux compute node; this runner rejects other hosts.
    try:
        storage = capacity.check()
    except Exception as error:
        failures.append(str(error))
        storage = {'error': str(error), 'campaign_peak_observed_bytes': capacity.peak}
    result = {'command': command, 'seconds': seconds, 'rss_kib': usage.ru_maxrss,
              'user_seconds': usage.ru_utime, 'system_seconds': usage.ru_stime,
              'returncode': process.returncode, 'capacity_failures': failures,
              'storage': storage}
    write_json(work / 'timing.json', result)
    if process.returncode or failures:
        raise RuntimeError(f'run failed in {work}; see timing.json and console.log')
    return result


def warm_index(directory, label):
    names = {'star': ['Genome', 'SA', 'SAindex'], 'nustar': ['reference.bin', 'index.sa']}
    for name in names[label]:
        with (directory / name).open('rb') as source:
            for _ in iter(lambda: source.read(8 << 20), b''):
                pass


def command_for(binary, index, inputs, case, threads, prefix, whitelist):
    command = [str(binary), '--genomeDir', str(index), '--readFilesIn', *map(str, inputs),
               '--runThreadN', str(threads), '--outFileNamePrefix', prefix,
               '--sjdbOverhang', '100', '--limitSjdbInsertNsj', '2000000',
               '--outSAMtype', *(['BAM', 'SortedByCoordinate'] if case['bam'] else ['None'])]
    if case['two_pass']:
        command += ['--twopassMode', 'Basic', '--twopass1readsN', '-1', '--sjdbInsertSave', 'Basic']
    if case['solo']:
        command += ['--soloType', 'CB_UMI_Simple', '--soloCBwhitelist', str(whitelist),
                    '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17',
                    '--soloUMIlen', '12', '--soloBarcodeReadLength', '0', '--soloFeatures', 'Gene',
                    '--soloStrand', 'Forward', '--soloCellFilter', *case['cell_filter']]
    if case['bam']:
        command += ['--outSAMunmapped', 'Within', '--outBAMcompression', '1',
                    '--limitBAMsortRAM', str(2 << 30), '--outSAMattributes', 'NH', 'HI', 'AS', 'nM']
        if case['solo']:
            command += ['CR', 'CY', 'UR', 'UY', 'CB', 'UB', 'GX', 'GN']
    return command


def execution_policy(root):
    """Read an explicit stage selection, including for already-spooled Slurm jobs."""
    path = root / 'execution-policy.json'
    if not path.exists():
        return None
    policy = json.loads(path.read_text())
    stages = policy.get('enabled_stages')
    if (policy.get('schema') != 1 or not isinstance(stages, list) or not stages
            or any(stage not in ['pilot', 'main', 'full'] for stage in stages)
            or len(stages) != len(set(stages)) or not policy.get('reason')):
        raise ValueError('invalid campaign execution policy')
    return {'policy': policy, 'sha256': sha256(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['root', 'manifest', 'star', 'nustar', 'star-index', 'nustar-index', 'samtools']:
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--stage', choices=['pilot', 'main', 'full'], required=True)
    parser.add_argument('--dataset', choices=['bulk', 'pbmc'], nargs='+', default=['bulk', 'pbmc'])
    parser.add_argument('--threads', type=int, nargs='+', required=True)
    parser.add_argument('--repeats', type=int, required=True)
    parser.add_argument('--disk-budget-gib', type=int, required=True)
    parser.add_argument('--reserve-gib', type=int, default=32)
    parser.add_argument('--prune-validated', action='store_true')
    args = parser.parse_args()
    if platform.system() != 'Linux' or not os.environ.get('SLURM_JOB_ID'):
        parser.error('run_suite must execute on a Linux Slurm compute allocation')
    if args.repeats < 1 or not args.threads or min(args.threads) < 1:
        parser.error('repeats and thread counts must be positive')
    affinity = os.sched_getaffinity(0)
    if max(args.threads) > len(affinity):
        parser.error('requested thread count exceeds process CPU affinity')
    root = args.root.resolve()
    policy = execution_policy(root)
    if policy is not None and args.stage not in policy['policy']['enabled_stages']:
        record = {'state': 'skipped', 'stage': args.stage, 'execution_policy': policy,
                  'job_id': os.environ['SLURM_JOB_ID']}
        write_json(root / 'logs' / f'skipped-{args.stage}-job-{os.environ["SLURM_JOB_ID"]}.json', record)
        print(json.dumps(record), flush=True)
        return
    manifest = json.loads(args.manifest.read_text())
    suite = root / 'results' / f'{args.stage}-job-{os.environ["SLURM_JOB_ID"]}'
    suite.mkdir(parents=True, exist_ok=False)
    capacity = Capacity(root, args.disk_budget_gib << 30, args.reserve_gib << 30)
    hardware = {'execution_policy': policy, 'platform': platform.platform(), 'node': platform.node(),
                'cpu_affinity': sorted(affinity), 'slurm': {k: v for k, v in os.environ.items() if k.startswith('SLURM_')},
                'proc_cpuinfo': Path('/proc/cpuinfo').read_text(),
                'proc_meminfo': Path('/proc/meminfo').read_text(),
                'binaries': {k: {'path': str(getattr(args, k).resolve()), 'sha256': sha256(getattr(args, k))}
                             for k in ['star', 'nustar', 'samtools']},
                'manifest_sha256': sha256(args.manifest), 'warm_index': True,
                'inputs': {dataset: json.loads((root / 'prepared' / dataset / args.stage / 'prepared.json').read_text())
                           for dataset in args.dataset},
                'cpu_sharing': 'See captured scontrol allocation; no isolation claim is inferred.'}
    for name, command in {'lscpu': ['lscpu'], 'allocation': ['scontrol', 'show', 'job', os.environ['SLURM_JOB_ID'], '-dd'],
                          'mounts': ['findmnt', '-T', str(root)], 'compiler': ['c++', '--version'],
                          'samtools_version': [str(args.samtools), '--version']}.items():
        completed = subprocess.run(command, capture_output=True, text=True)
        hardware[name] = {'returncode': completed.returncode, 'stdout': completed.stdout, 'stderr': completed.stderr}
    hardware['runner_arguments'] = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    hardware['index_manifests'] = {}
    for label, index in [('star', args.star_index), ('nustar', args.nustar_index)]:
        saved = json.loads((index / 'validated-build.json').read_text())
        for name, info in saved['files'].items():
            if sha256(index / name) != info['sha256']:
                raise RuntimeError('index checksum changed: ' + name)
        hardware['index_manifests'][label] = saved
    for dataset, prepared in hardware['inputs'].items():
        for name, info in prepared['files'].items():
            if sha256(root / 'prepared' / dataset / args.stage / name) != info['sha256']:
                raise RuntimeError('prepared FASTQ checksum changed')
    write_json(suite / 'metadata.json', hardware)
    programs = {'star': args.star.resolve(), 'nustar': args.nustar.resolve()}
    indices = {'star': args.star_index.resolve(), 'nustar': args.nustar_index.resolve()}
    results = []
    try:
        cases = [c for c in manifest['cases'] if args.stage in c['stages'] and c['dataset'] in args.dataset]
        if not cases:
            raise ValueError('no cases match the selected stage and dataset')
        for case, threads, repeat in itertools.product(cases, args.threads, range(args.repeats)):
            inputs = [root / path.format(stage=args.stage) for path in case['inputs']]
            for path in inputs:
                if not path.is_file():
                    raise FileNotFoundError(path)
            work = suite / f'{case["name"]}-{threads}-r{repeat + 1}'
            work.mkdir()
            write_json(suite / 'status.json', {'state': 'running', 'case': case['name'], 'threads': threads, 'repeat': repeat + 1})
            entry = {'case': case['name'], 'threads': threads, 'repeat': repeat + 1, 'engines': {}}
            order = ('star', 'nustar')[::1 if repeat % 2 == 0 else -1]
            entry['run_order'] = order
            for label in order:
                engine = work / label
                engine.mkdir()
                prefix = str(engine / 'output.')
                warm_index(indices[label], label)
                command = command_for(programs[label], indices[label], inputs, case, threads, prefix,
                                      root / 'reference/3M-february-2018.txt')
                entry['engines'][label] = measure(command, engine, capacity)
            entry['validation'] = validate_pair(str(work / 'star/output.'), str(work / 'nustar/output.'),
                                                case, args.samtools, work)
            expected_pairs = hardware['inputs'][case['dataset']]['pairs']
            if int(entry['validation']['statistics']['Number of input reads']) != expected_pairs:
                raise AssertionError('mapped input count differs from the prepared dataset')
            write_json(work / 'validated.json', entry)
            if args.prune_validated:
                for label in programs:
                    for path in (work / label).glob('output.Aligned*.bam'):
                        path.unlink()
            results.append(entry)
            write_json(suite / 'results.json', results)
        summary = []
        for name, threads in sorted({(r['case'], r['threads']) for r in results}):
            group = [r for r in results if r['case'] == name and r['threads'] == threads]
            row = {'case': name, 'threads': threads, 'repeats': len(group)}
            for label in programs:
                times = [r['engines'][label]['seconds'] for r in group]
                row[label] = {'median_seconds': statistics.median(times), 'min_seconds': min(times),
                              'max_seconds': max(times), 'max_rss_kib': max(r['engines'][label]['rss_kib'] for r in group)}
            row['star_over_nustar'] = row['star']['median_seconds'] / row['nustar']['median_seconds']
            summary.append(row)
        write_json(suite / 'summary.json', summary)
        write_json(suite / 'status.json', {'state': 'complete', 'validated_pairs': len(results)})
    except Exception as error:
        write_json(suite / 'status.json', {'state': 'failed', 'error': repr(error), 'validated_pairs': len(results)})
        raise


if __name__ == '__main__':
    main()
