"""Serial STAR/NuSTAR feature benchmark; wait4 reports child peak RSS and CPU.

A provided barcode FASTQ is biological input. --synthetic-barcodes instead adds
artificial CB/UMIs to bulk RNA reads, for implementation measurements only.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import time


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('nustar', 'star', 'nustar-index', 'star-index', 'reads', 'work', 'lock'):
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--barcodes', type=Path)
    p.add_argument('--whitelist', type=Path)
    p.add_argument('--synthetic-barcodes', action='store_true')
    p.add_argument('--runThreadN', dest='threads', nargs='+', type=int, default=[1, 32])
    p.add_argument('--modes', nargs='+', choices=['twopass', 'solo'], default=['twopass', 'solo'])
    p.add_argument('--input-bam', action='store_true')
    p.add_argument('--bam-reader', type=Path, help='BAM to SAM converter required by the upstream STAR reader')
    p.add_argument('--bam-output', choices=['Unsorted', 'SortedByCoordinate'])
    p.add_argument('--warm-index', action='store_true', help='Read index files before each timed process')
    return p.parse_args()


def prepare_barcodes(reads, work):
    rng = random.Random(911)
    cbs = [''.join(rng.choice('ACGT') for _ in range(16)) for _ in range(128)]
    whitelist, barcodes = work / 'whitelist.txt', work / 'barcodes.fq'
    whitelist.write_text('\n'.join(cbs) + '\n')
    count = 0
    with reads.open() as source, barcodes.open('w') as output:
        while header := source.readline():
            for _ in range(3):
                source.readline()
            seq = cbs[count % len(cbs)] + ''.join(rng.choice('ACGT') for _ in range(12))
            output.write(header + seq + '\n+\n' + 'I' * 28 + '\n')
            count += 1
    return barcodes, whitelist, count


def warm_index(directory, label):
    names = ['Genome', 'SA', 'SAindex'] if label == 'star' else ['reference.bin', 'index.sa']
    for name in names:
        with (directory / name).open('rb') as f:
            while f.read(32 << 20):
                pass


def compare(work, threads, mode):
    a, b = (str(work / f'{mode}-{label}-{threads}') + '.' for label in ('star', 'nustar'))
    files = ['SJ.out.tab']
    if mode == 'twopass':
        files += ['_STARpass1/SJ.out.tab', '_STARgenome/sjdbInfo.txt', '_STARgenome/sjdbList.out.tab']
    else:
        solo_a, solo_b = Path(a + 'Solo.out'), Path(b + 'Solo.out')
        paths = lambda d: {str(p.relative_to(d)) for p in d.rglob('*') if p.is_file()}
        assert paths(solo_a) == paths(solo_b), 'Solo output file lists differ'
        files += ['Solo.out/' + p for p in sorted(paths(solo_a))]
    for suffix in files:
        assert Path(a + suffix).read_bytes() == Path(b + suffix).read_bytes(), suffix
    # Exclude timestamps and measured throughput, retain every biological statistic.
    def stats(prefix):
        return {line.split('|')[0].strip(): line.split('|')[1].strip()
                for line in Path(prefix + 'Log.final.out').read_text().splitlines()
                if '|' in line and not any(word in line.lower() for word in ('time', 'speed', 'started', 'finished'))}
    assert stats(a) == stats(b), 'Log.final.out biological statistics differ'
    return files


def main():
    args = arguments()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    programs = {'star': args.star.resolve(), 'nustar': args.nustar.resolve()}
    indices = {'star': args.star_index.resolve(), 'nustar': args.nustar_index.resolve()}
    reads = args.reads.resolve()
    metadata = {'RNA_source': str(reads), 'warm_index': args.warm_index,
                'bam_output': args.bam_output, 'bam_input': args.input_bam,
                'binary_sha256': {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in programs.items()}}
    if args.synthetic_barcodes:
        args.barcodes, args.whitelist, metadata['reads'] = prepare_barcodes(reads, work)
        metadata['barcodes'] = 'synthetic CB/UMI, seed 911; not biological single-cell data'
    elif args.barcodes:
        metadata['barcodes'] = str(args.barcodes.resolve())
    if 'solo' in args.modes and ((not args.barcodes and not args.input_bam) or not args.whitelist):
        raise ValueError('Solo requires --barcodes/--whitelist or --synthetic-barcodes')
    (work / 'metadata.json').write_text(json.dumps(metadata, indent=2))
    results, checks = [], []
    with args.lock.open('r') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for mode in args.modes:
            for threads in args.threads:
                for label in ('star', 'nustar'):
                    name = f'{mode}-{label}-{threads}'
                    if args.warm_index:
                        warm_index(indices[label], label)
                    command = [str(programs[label]), '--genomeDir', str(indices[label]), '--readFilesIn', str(reads)]
                    if mode == 'solo' and not args.input_bam:
                        command += [str(args.barcodes.resolve())]
                    if args.input_bam:
                        command += ['--readFilesType', 'SAM', 'SE']
                        if label == 'star':
                            if not args.bam_reader: raise ValueError('--input-bam needs --bam-reader for STAR')
                            command += ['--readFilesCommand', str(args.bam_reader.resolve())]
                    command += ['--runThreadN', str(threads), '--outFileNamePrefix', str(work/name) + '.', '--outSAMtype', *(['BAM', args.bam_output] if args.bam_output else ['None'])]
                    if args.bam_output:
                        command += ['--outSAMunmapped', 'Within', '--outBAMcompression', '1', '--limitBAMsortRAM', str(2 << 30)]
                    if mode == 'twopass':
                        command += ['--twopassMode', 'Basic']
                    else:
                        command += ['--soloType', 'CB_UMI_Simple', '--soloCBwhitelist', str(args.whitelist.resolve()),
                                    '--soloUMIlen', '12', '--soloCellFilter', 'None']
                        if args.input_bam:
                            command += ['--soloInputSAMattrBarcodeSeq', 'CR', 'UR', '--soloInputSAMattrBarcodeQual', 'CY', 'UY']
                        if args.bam_output:
                            command += ['--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'CR', 'CY', 'UR', 'UY']
                            if args.bam_output == 'SortedByCoordinate':
                                command += ['CB', 'UB', 'GX', 'GN']
                    with (work / (name + '.console')).open('w') as log:
                        started = time.monotonic()
                        process = subprocess.Popen(command, cwd=work, stdout=log, stderr=subprocess.STDOUT)
                        _, status, usage = os.wait4(process.pid, 0)
                        process.returncode = os.waitstatus_to_exitcode(status)
                        result = {'name': name, 'command': command, 'seconds': time.monotonic()-started,
                                  'rss_kib': usage.ru_maxrss, 'returncode': process.returncode,
                                  'user_seconds': usage.ru_utime, 'system_seconds': usage.ru_stime}
                    (work / (name + '.json')).write_text(json.dumps(result, indent=2))
                    results.append(result)
                    print(json.dumps(result), flush=True)
                    if process.returncode:
                        raise RuntimeError(name)
                checks.append({'mode': mode, 'threads': threads, 'matching_files': compare(work, threads, mode)})
    (work / 'results.json').write_text(json.dumps(results, indent=2))
    (work / 'validation.json').write_text(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
