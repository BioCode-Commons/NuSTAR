#!/usr/bin/env python3
"""Measure one STAR index build per reference, verify products, then remove them."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--star', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=32)
    parser.add_argument('--ram-bytes', type=int, default=120000000000)
    args = parser.parse_args()
    if args.threads < 1 or args.ram_bytes < 1:
        parser.error('threads and ram-bytes must be positive')
    root, work = args.root.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    binary = work / 'STAR'
    shutil.copy2(args.star, binary)
    shutil.copy2(__file__, work / Path(__file__).name)
    state = {'status': 'waiting-for-benchmark-lock', 'rounds': 1, 'threads': args.threads,
             'ram_limit_bytes': args.ram_bytes, 'program_sha256': sha(binary),
             'version': subprocess.check_output([binary, '--version'], text=True).strip(),
             'cpu_affinity': sorted(os.sched_getaffinity(0)), 'runs': [], 'checks': [], 'removed': []}
    cases = [
        ('test', root/'data/test/test.fa', root/'data/test/test.gtf', 10,
         root/'NuSTAR/build/defensive-check-cleanup/test-star-index'),
        ('human', root/'data/references/human/GRCh38.primary_assembly.genome.fa',
         root/'data/references/human/gencode.v50.primary_assembly.annotation.gtf', 14,
         root/'GitHubAction/runs/human-profile-fix/work/star-index')]
    state['inputs'] = {str(p): {'bytes': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
                       for _, fasta, gtf, _, _ in cases for p in [fasta, gtf]}

    def save():
        tmp = work / 'measurement.tmp.json'
        tmp.write_text(json.dumps(state, indent=2) + '\n')
        tmp.replace(work / 'measurement.json')
        lines = ['# STAR index construction measurement', '',
                 f"Status: {state['status']}; current: {state.get('current', 'waiting')}.", '',
                 f"STAR {state['version']}, {args.threads} threads, one round per reference; limitGenomeGenerateRAM={args.ram_bytes} bytes.",
                 'Full reference and GTF; overhang 150, chromosome bin bits 18, sparse SA distance 1, junction limit 2000000.',
                 'Wall time and peak RSS cover the STAR process only. Product checks and cleanup are excluded.', '',
                 '| Reference | Seconds | Peak RSS GiB |', '| --- | ---: | ---: |']
        for r in state['runs']:
            lines.append(f"| {r['case']} | {r['seconds']:.3f} | {r['rss_kib']/2**20:.3f} |")
        lines += ['', 'Products are compared in full against retained, previously verified STAR indices before cleanup.',
                  'Detailed commands, input metadata, checksums and cleanup records are in `measurement.json`.']
        if state.get('error'):
            lines += ['', 'Error: ' + state['error']]
        (work / 'report.md').write_text('\n'.join(lines) + '\n')

    save()
    try:
        with (root/'GitHubAction/benchmark.lock').open('r') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            for case, fasta, gtf, prefix, baseline in cases:
                directory = work / (case + '-index')
                directory.mkdir()
                command = [str(binary), '--runMode', 'genomeGenerate', '--genomeDir', str(directory),
                           '--genomeFastaFiles', str(fasta), '--sjdbGTFfile', str(gtf),
                           '--sjdbOverhang', '150', '--genomeChrBinNbits', '18',
                           '--genomeSAindexNbases', str(prefix), '--genomeSAsparseD', '1',
                           '--limitSjdbInsertNsj', '2000000', '--runThreadN', str(args.threads),
                           '--limitGenomeGenerateRAM', str(args.ram_bytes),
                           '--outFileNamePrefix', str(work/(case + '.'))]
                state.update(status='running', current=case+'-build')
                save()
                start = time.monotonic()
                with (work/(case+'.console')).open('wb') as output:
                    process = subprocess.Popen(command, cwd=work, stdout=output, stderr=subprocess.STDOUT)
                    state['pid'] = process.pid
                    save()
                    _, status, usage = os.wait4(process.pid, 0)
                    process.returncode = os.waitstatus_to_exitcode(status)
                record = {'case': case, 'command': command, 'exit': process.returncode,
                          'seconds': time.monotonic()-start, 'rss_kib': usage.ru_maxrss,
                          'user_seconds': usage.ru_utime, 'system_seconds': usage.ru_stime}
                state['runs'].append(record)
                state.pop('pid', None)
                save()
                if process.returncode:
                    raise RuntimeError('STAR build failed: '+case)
                state['current'] = case+'-verify'
                save()
                products = ['Genome', 'SA', 'SAindex', 'chrName.txt', 'chrLength.txt', 'chrStart.txt',
                            'chrNameLength.txt', 'geneInfo.tab', 'exonGeTrInfo.tab', 'transcriptInfo.tab',
                            'exonInfo.tab', 'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'sjdbInfo.txt']
                evidence = {}
                for name in products:
                    left, right = directory/name, baseline/name
                    evidence[name] = {'bytes': left.stat().st_size, 'sha256': sha(left),
                                      'baseline_bytes': right.stat().st_size, 'baseline_sha256': sha(right)}
                    item = evidence[name]
                    if item['bytes'] != item['baseline_bytes'] or item['sha256'] != item['baseline_sha256']:
                        raise RuntimeError('STAR product differs: '+case+'/'+name)
                state['checks'].append({'case': case, 'equal': True, 'baseline': str(baseline), 'products': evidence})
                shutil.copy2(directory/'genomeParameters.txt', work/(case+'.genomeParameters.txt'))
                if (directory/'Log.out').exists():
                    shutil.copy2(directory/'Log.out', work/(case+'.Log.out'))
                save()
                entry = {'path': str(directory), 'reason': 'All index products verified', 'deleted': False}
                state['removed'].append(entry)
                save()
                shutil.rmtree(directory)
                entry['deleted'] = True
                save()
            state.update(status='complete', current='finished')
            save()
    except BaseException as error:
        state.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    main()
