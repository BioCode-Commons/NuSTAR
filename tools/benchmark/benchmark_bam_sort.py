#!/usr/bin/env python3
"""Serial full-input BAM acceptance: frozen old/new NuSTAR and STAR, rotated runs.

Experimental evidence only. Never updates published benchmark tables.
"""
import argparse
import csv
import datetime
import fcntl
import json
import os
from pathlib import Path
import shutil
import statistics as stats
import subprocess
import time

from benchmark_alignment import metadata, sha, statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--star', type=Path, required=True)
    parser.add_argument('--compare', type=Path, required=True)
    parser.add_argument('--rounds', type=int, default=1)
    parser.add_argument('--cases', nargs='+', choices=['se', 'pe'], default=['se', 'pe'])
    parser.add_argument('--reuse-se', type=Path,
                        help='Reuse completed SE candidate/STAR pilot and its explicitly preserved STAR BAM')
    parser.add_argument('--reuse-se-star', type=Path,
                        help='Reuse only the measured STAR SE baseline after changing the candidate')
    parser.add_argument('--reuse-reference', type=Path,
                        help='Reuse a completed SE/PE run: its candidate becomes the baseline; reuse STAR too')
    parser.add_argument('--reference-bams', type=Path,
                        help='Directory of preserved hardlinks named <run-directory>-<original-BAM-name>')
    parser.add_argument('--wait-lock', action='store_true', help='Wait for an active benchmark to finish')
    parser.add_argument('--keep-bams', action='store_true', help='Keep validated outputs for subsequent comparisons')
    parser.add_argument('--fail-fast-se', action='store_true',
                        help='After validating a slower SE candidate, stop before old NuSTAR and PE runs')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    work = args.output.resolve()
    work.mkdir(parents=True, exist_ok=False)
    lock = (root / 'GitHubAction/benchmark.lock').open('a')
    sources = {'candidate': args.candidate.resolve(), 'star': args.star.resolve()}
    if args.baseline: sources['baseline'] = args.baseline.resolve()
    programs = {}
    for name, source in sources.items():
        programs[name] = work / name
        shutil.copy2(source, programs[name])
    comparator = work / 'compare-bam'
    shutil.copy2(args.compare.resolve(), comparator)
    indexes = {
        'candidate': root / 'NuSTAR/build/defensive-check-cleanup/full-sam/human-after-index',
        'star': root / 'GitHubAction/runs/human-profile-fix/work/star-index',
    }
    indexes['baseline'] = indexes['candidate']
    reads = [root / f'data/fastq/SRR21781765_{i}.fastq' for i in (1, 2)]
    state = {'status': 'waiting', 'queued_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'rounds': args.rounds, 'cases': args.cases, 'runs': [], 'comparisons': [], 'cleanup': [],
             'programs': {k: {'source': str(sources[k]), 'sha256': sha(p)} for k, p in programs.items()},
             'comparator_sha256': sha(comparator), 'inputs': {str(p): metadata(p) for p in reads},
             'indexes': {k: str(p) for k, p in indexes.items()},
             'method': 'Serial rotated runs; 32 threads, 96 GiB sort budget, compression 1; '
                       'no pinning, prewarming or cache eviction. Wall/RSS cover loading through output. '
                       'Validation and cleanup excluded. No automatic retry.',
             'git_revision': subprocess.check_output(['git', '-C', root/'NuSTAR', 'rev-parse', 'HEAD'], text=True).strip()}
    (work/'source.diff').write_bytes(subprocess.check_output(['git', '-C', root/'NuSTAR', 'diff', 'HEAD']))
    state['sort_sources'] = {}
    for source in sorted((root/'NuSTAR/bam').glob('*')):
        if source.is_file():
            target = work/'source'/source.name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(source, target)
            state['sort_sources'][source.name] = sha(source)
    shutil.copy2(Path(__file__).resolve(), work/'benchmark_bam_sort.py')

    def save():
        tmp = work/'state.tmp.json'
        tmp.write_text(json.dumps(state, indent=2)+'\n'); tmp.replace(work/'state.json')
        with (work/'summary.csv').open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['case', 'round', 'program', 'seconds', 'peak_RSS_GiB', 'exit'])
            for r in state['runs']:
                writer.writerow([r['case'], r['round'], r['program'], r['seconds'], r['rss_kib']/2**20, r['exit']])
        lines = ['# Experimental BAM sort acceptance', '', f"Status: {state['status']}; current: {state.get('current', '')}.", '', state['method'], '',
                 '| Case | Round | Program | Seconds | Peak RSS GiB | Exit |', '| --- | ---: | --- | ---: | ---: | ---: |']
        for r in state['runs']:
            lines.append(f"| {r['case']} | {r['round']} | {r['program']} | {r['seconds']:.3f} | {r['rss_kib']/2**20:.3f} | {r['exit']} |")
        if state.get('validation_progress'):
            lines += ['', 'Validation input bytes read (not a correctness result): `'+json.dumps(state['validation_progress'])+'`']
        lines += ['', 'Correctness requires exact decoded record order, references, coordinate monotonicity, '
                  'BGZF EOF, all primary templates, SJ, and all 28 non-time statistics. '
                  'Only validated BAM outputs are deleted. Failed artifacts are retained.', '',
                  'Acceptance requires both SE and PE, all correctness checks passing, '
                  'and candidate time below both STAR and the frozen old NuSTAR. '
                  'One round was requested; single measurements do not establish statistical significance.']
        if 'acceptance' in state: lines += ['', 'Acceptance: `'+json.dumps(state['acceptance'])+'`']
        if 'reused_se' in state:
            lines += ['', 'The measured STAR SE run is reused from `'+state['reused_se']+
                      '` with matching binary/input identities and a preserved BAM; STAR SE is not rerun.']
        if 'reused_reference' in state:
            lines += ['', 'STAR and the preceding candidate (now baseline) are reused from `'+state['reused_reference']+
                      '`. Binary identities, inputs, indexes, commands, and prior correctness are checked. '
                      'The new candidate runs once per case and is compared directly with both preserved BAMs. '
                      'These are sequential historical measurements, not interleaved contemporary repeats.']
        if 'error' in state: lines += ['', 'Error: '+state['error']]
        (work/'report.md').write_text('\n'.join(lines)+'\n')

    save()
    try:
        # Freeze the executable and matching sources before waiting, so later
        # workspace edits cannot be mistaken for the measured implementation.
        fcntl.flock(lock, fcntl.LOCK_EX | (0 if args.wait_lock else fcntl.LOCK_NB))
        state.update(status='running', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()); save()
        audited = json.loads((root/'NuSTAR/build/full-star-profile-unified-20261113/state.json').read_text())
        for p in reads:
            if metadata(p) != audited['inputs'][str(p)]:
                raise RuntimeError('Full input metadata differs from audited input: '+str(p))
        if subprocess.check_output([programs['star'], '--version'], text=True).strip() != '2.7.11b':
            raise RuntimeError('Expected unmodified STAR 2.7.11b')
        env = os.environ.copy(); env.update(OMP_NUM_THREADS='32', LC_ALL='C')
        for k in list(env):
            if k.startswith(('STAR_STAGE_TIMING_', 'NUSTAR_STAGE_TIMING_')): del env[k]
        labels = list(programs)
        reused = None
        reused_labels = set()
        reference = None
        if args.reuse_reference:
            if args.reuse_se or args.reuse_se_star or not args.baseline or args.rounds != 1:
                raise RuntimeError('Reference reuse requires a baseline, one round, and no SE reuse flags')
            reference = args.reuse_reference.resolve()
            previous = json.loads((reference/'state.json').read_text())
            if previous['status'] not in ('complete', 'not-accepted') or previous['rounds'] != 1:
                raise RuntimeError('Reference run has not finished successfully')
            if not set(args.cases).issubset(previous['cases']):
                raise RuntimeError('Reference cases missing')
            if previous['inputs'] != state['inputs']:
                raise RuntimeError('Reference inputs differ')
            if previous['comparator_sha256'] != state['comparator_sha256']:
                raise RuntimeError('Reference comparator differs')
            for old, new in [('candidate', 'baseline'), ('star', 'star')]:
                if previous['programs'][old]['sha256'] != state['programs'][new]['sha256']:
                    raise RuntimeError('Reference binary differs: '+new)
                if previous['indexes'][old] != state['indexes'][new]:
                    raise RuntimeError('Reference index differs: '+new)
            expected = {'--runThreadN': '32', '--sjdbOverhang': '150', '--limitSjdbInsertNsj': '2000000',
                        '--outBAMcompression': '1', '--limitBAMsortRAM': str(96 << 30),
                        '--outSAMunmapped': 'Within', '--outSAMorder': 'Paired'}
            for case in args.cases:
                checks = [c for c in previous['comparisons'] if c['case']==case and c['program']=='candidate']
                if len(checks)!=1 or not all(checks[0].get(k) for k in ('equal','complete','sj_equal','statistics_equal')):
                    raise RuntimeError('Reference candidate correctness missing: '+case)
                for old, new in [('candidate', 'baseline'), ('star', 'star')]:
                    runs = [r for r in previous['runs'] if r['case']==case and r['program']==old]
                    if len(runs)!=1 or runs[0]['exit']:
                        raise RuntimeError('Reference measurement missing: '+case+' '+old)
                    run = runs[0]; command = run['command']
                    read_start = command.index('--readFilesIn')+1
                    if command[read_start:read_start+(1 if case=='se' else 2)] != list(map(str, reads[:1 if case=='se' else 2])):
                        raise RuntimeError('Reference reads differ')
                    if command[command.index('--genomeDir')+1] != str(indexes[new]):
                        raise RuntimeError('Reference command index differs')
                    for flag, value in expected.items():
                        if flag not in command or command[command.index(flag)+1] != value:
                            raise RuntimeError('Reference command differs: '+flag)
                    if command[command.index('--outSAMtype')+1:command.index('--outSAMtype')+3] != ['BAM','SortedByCoordinate']:
                        raise RuntimeError('Reference output mode differs')
                    state['runs'].append(dict(run, program=new, reused_from=str(reference), source_program=old))
            state['reused_reference'] = str(reference); save()

        def reference_output(case, label):
            run = next(r for r in state['runs'] if r['case']==case and r['program']==label)
            command = run['command']
            prefix = Path(command[command.index('--outFileNamePrefix')+1])
            bam = Path(str(prefix)+'Aligned.sortedByCoord.out.bam')
            candidates = [bam, Path(str(prefix)+'preserved.bam')]
            if args.reference_bams:
                candidates += [args.reference_bams/(prefix.parent.name+'-'+p.name) for p in candidates[:]]
            found = next((p for p in candidates if p.is_file()), None)
            if found is None: raise RuntimeError('Missing preserved reference BAM: '+str(bam))
            return prefix, found
        if args.reuse_se and args.reuse_se_star: raise RuntimeError('Choose one SE reuse mode')
        if args.reuse_se or args.reuse_se_star:
            reused = (args.reuse_se or args.reuse_se_star).resolve()
            reused_labels = {'candidate', 'star'} if args.reuse_se else {'star'}
            previous = json.loads((reused/'state.json').read_text())
            if args.rounds != 1 or previous['rounds'] != 1 or previous['cases'] != ['se'] or previous['status'] not in ('complete', 'not-accepted'):
                raise RuntimeError('Only a completed one-round SE pilot can be reused')
            for label in reused_labels:
                if previous['programs'][label]['sha256'] != state['programs'][label]['sha256']:
                    raise RuntimeError('Pilot binary differs: '+label)
            if 'candidate' in reused_labels and previous['sort_sources'] != state['sort_sources']:
                raise RuntimeError('Candidate BAM sources changed since pilot; old timings cannot validate the new implementation')
            if any(previous['inputs'][str(p)] != metadata(p) for p in reads):
                raise RuntimeError('Pilot inputs differ')
            if not (reused/'se-1-star.preserved.bam').is_file():
                raise RuntimeError('Missing explicitly preserved pilot STAR BAM')
            checks = previous['comparisons']
            if len(checks) != 1 or not all(checks[0].get(k) for k in ('equal','complete','sj_equal','statistics_equal')):
                raise RuntimeError('Pilot correctness did not pass')
            pilot_times = {r['program']: r['seconds'] for r in previous['runs']}
            if 'candidate' in reused_labels and pilot_times['candidate'] >= pilot_times['star']:
                raise RuntimeError('SE pilot did not beat STAR; stop before running the remaining acceptance cases')
            state['runs'].extend(dict(r, reused_from=str(reused)) for r in previous['runs'] if r['program'] in reused_labels)
            if 'candidate' in reused_labels:
                state['comparisons'].extend(dict(c, reused_from=str(reused)) for c in checks)
            state['reused_se'] = str(reused); save()
        for repeat in range(args.rounds):
            order = labels[repeat % len(labels):] + labels[:repeat % len(labels)]
            for case in args.cases:
                prefixes = {}
                bam_paths = {}
                slow_se = False
                for label in order:
                    if slow_se and label == 'baseline': continue
                    if reference and label in ('star', 'baseline'):
                        prefixes[label], bam_paths[label] = reference_output(case, label)
                        continue
                    if reused and case == 'se' and label in reused_labels:
                        prefixes[label] = reused/f'se-1-{label}.'
                        bam_paths[label] = reused/f'se-1-{label}.preserved.bam'
                        continue
                    if shutil.disk_usage(work).free < 700 * 2**30:
                        raise RuntimeError('Less than 700 GiB free for full BAM run; no automatic cleanup')
                    prefix = work / f'{case}-{repeat+1}-{label}.'
                    prefixes[label] = prefix
                    bam_paths[label] = Path(str(prefix)+'Aligned.sortedByCoord.out.bam')
                    command = [str(programs[label]), '--genomeDir', str(indexes[label]),
                               '--readFilesIn', *map(str, reads[:1 if case == 'se' else 2]),
                               '--runThreadN', '32', '--sjdbOverhang', '150', '--limitSjdbInsertNsj', '2000000',
                               '--outFileNamePrefix', str(prefix), '--outSAMtype', 'BAM', 'SortedByCoordinate',
                               '--outSAMunmapped', 'Within', '--outSAMorder', 'Paired', '--outBAMcompression', '1',
                               '--limitBAMsortRAM', str(96 << 30)]
                    command += ['--outBAMsortingThreadN', '32'] if label == 'star' else ['--profile-runtime', '--io-backend', 'pread']
                    state['current'] = f'{case} round {repeat+1} {label}'; save()
                    with Path(str(prefix)+'console.log').open('w') as log:
                        start = time.monotonic()
                        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env)
                        _, status, usage = os.wait4(child.pid, 0)
                        elapsed = time.monotonic()-start
                        child.returncode = os.waitstatus_to_exitcode(status)
                    result = {'case': case, 'round': repeat+1, 'program': label, 'command': command,
                              'seconds': elapsed, 'rss_kib': usage.ru_maxrss, 'user_seconds': usage.ru_utime,
                              'system_seconds': usage.ru_stime, 'exit': child.returncode}
                    state['runs'].append(result); save()
                    if child.returncode: raise RuntimeError('Aligner failed: '+state['current'])
                    if args.fail_fast_se and reused and case == 'se' and label == 'candidate':
                        star_time = next(r['seconds'] for r in state['runs'] if r['case']=='se' and r['program']=='star')
                        slow_se = elapsed >= star_time
                comparisons = [('candidate','star'),('candidate','baseline')] if reference else [(label,'star') for label in labels if label!='star']
                for label, reference_label in comparisons:
                    if label not in prefixes: continue
                    if reused and case == 'se' and label == 'candidate' and label in reused_labels: continue
                    state['current'] = f'validate {case} round {repeat+1} {label} vs {reference_label}'; save()
                    left, right = prefixes[label], prefixes[reference_label]
                    suffix = f'-vs-{reference_label}' if reference else ''
                    diff = work/f'{case}-{repeat+1}-{label}{suffix}.differences.txt'
                    comparison_path = work/f'{case}-{repeat+1}-{label}{suffix}.comparison.json'
                    stderr_path = work/f'{case}-{repeat+1}-{label}{suffix}.comparison.stderr'
                    with comparison_path.open('w') as output, stderr_path.open('w') as errors:
                        checked = subprocess.Popen([comparator, bam_paths[label], bam_paths[reference_label], diff], stdout=output, stderr=errors)
                        paths = {str(bam_paths[k].resolve()): k for k in (label, reference_label)}
                        while checked.poll() is None:
                            time.sleep(5)
                            if checked.poll() is not None: break
                            progress = {}
                            # Inspect only this comparator child's known BAM descriptors.
                            # This observes sequential compressed input reads; it does
                            # not instrument or change the frozen comparator binary.
                            for fd in (Path('/proc')/str(checked.pid)/'fd').glob('*'):
                                try:
                                    path = str(fd.resolve())
                                    if path not in paths: continue
                                    info = (fd.parent.parent/'fdinfo'/fd.name).read_text()
                                    offset = int(next(line.split()[1] for line in info.splitlines() if line.startswith('pos:')))
                                    progress[paths[path]] = round(100*offset/Path(path).stat().st_size, 1)
                                except (OSError, StopIteration, ValueError): pass
                            state['validation_progress'] = progress; save()
                        checked.wait()
                    check = json.loads(comparison_path.read_text()) if comparison_path.stat().st_size else {'error': stderr_path.read_text()}
                    check.update(case=case, round=repeat+1, program=label, reference=reference_label,
                                 sj_equal=Path(str(left)+'SJ.out.tab').read_bytes() == Path(str(right)+'SJ.out.tab').read_bytes(),
                                 statistics_equal=statistics(Path(str(left)+'Log.final.out')) == statistics(Path(str(right)+'Log.final.out')))
                    check['complete'] = check.get('left_primary_templates') == check.get('right_primary_templates') == 535736296
                    state['comparisons'].append(check); state.pop('validation_progress', None); save()
                    if checked.returncode or not all(check.get(k) for k in ('equal', 'sj_equal', 'statistics_equal', 'complete')):
                        raise RuntimeError('BAM validation failed: '+state['current'])
                if slow_se:
                    state.update(status='not-accepted', current='stopped after SE',
                                 acceptance={'passed': False, 'qualified_protocol': False,
                                             'reason': 'SE candidate did not beat STAR; remaining cases not run'})
                    save(); return
                for label, path in bam_paths.items():
                    if args.keep_bams or (reference and label in ('star','baseline')): continue
                    if reused and case == 'se' and label == 'candidate' and label in reused_labels: continue
                    state['cleanup'].append({'path': str(path), 'bytes': path.stat().st_size, 'validated': True}); save()
                    path.unlink()
        medians = {case: {label: stats.median(r['seconds'] for r in state['runs'] if r['case']==case and r['program']==label)
                          for label in labels} for case in args.cases}
        qualified = args.rounds >= 1 and set(args.cases)=={'se', 'pe'} and 'baseline' in programs
        faster = all(m['candidate'] < m['star'] and m['candidate'] < m.get('baseline', float('inf')) for m in medians.values())
        state['acceptance'] = {'qualified_protocol': qualified, 'medians_seconds': medians, 'faster': faster,
                               'passed': qualified and faster}
        state['status'] = 'complete' if qualified and faster else 'not-accepted'
        state['current'] = 'finished'; save()
    except BaseException as error:
        state.update(status='failed', error=repr(error)); save(); raise


if __name__ == '__main__': main()
