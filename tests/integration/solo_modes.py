"""Additional native Solo adapters: paired reads, inline CB, SmartSeq and filters.
Run solo.py first to construct the common annotated reference/reads.
"""
from pathlib import Path
import json, subprocess, sys
binary, star, work = (Path(x).resolve() for x in sys.argv[1:4])

def run(exe, args, name):
    with (work / (name + '.log')).open('w') as log:
        result = subprocess.run([str(exe), *map(str, args)], cwd=work, stdout=log, stderr=subprocess.STDOUT)
    assert result.returncode == 0, (name, result.returncode, (work / (name + '.log')).read_text()[-4000:])

def compare_dirs(a, b, name):
    left = {str(p.relative_to(a)): p.read_bytes() for p in a.rglob('*') if p.is_file()}
    right = {str(p.relative_to(b)): p.read_bytes() for p in b.rglob('*') if p.is_file()}
    assert left.keys() == right.keys(), (name, 'file list', left.keys(), right.keys())
    for filename in left:
        assert left[filename] == right[filename], (name, filename)

def fastq(path):
    lines = path.read_text().splitlines()
    return [lines[i:i+4] for i in range(0, len(lines), 4)]

rna = fastq(work / 'reads.fq')
barcodes = fastq(work / 'barcodes.fq')
(work/'short-barcodes.fq').write_text(''.join('\n'.join(b)+'\n' for b in barcodes[:40]))
inline = [[r[0], b[1] + r[1], '+', b[3] + r[3]] for r, b in zip(rna, barcodes)]
(work / 'inline.fq').write_text(''.join('\n'.join(r) + '\n' for r in inline))
(work / 'manifest-se.tsv').write_text(''.join(f'{work / "reads.fq"}\t-\tcell{i}\n' for i in range(3)))
(work / 'manifest.tsv').write_text(''.join(f'{work / "reads.fq"}\t{work / "mates.fq"}\tcell{i}\n' for i in range(3)))
base = ['--soloType', 'CB_UMI_Simple', '--soloCBstart', 1, '--soloCBlen', 8,
        '--soloUMIstart', 9, '--soloUMIlen', 6, '--soloCBwhitelist', work / 'whitelist.txt']
cases = [
    ('read-limit', ['--readFilesIn', work/'reads.fq', work/'short-barcodes.fq', '--readMapNumber', 40], base),
    ('paired', ['--readFilesIn', work/'reads.fq', work/'mates.fq', work/'barcodes.fq'], base),
    ('no-whitelist', ['--readFilesIn', work/'reads.fq', work/'barcodes.fq'], base[:-2] + ['--soloCBwhitelist', 'None']),
    ('inline', ['--readFilesIn', work/'inline.fq'], base + ['--soloBarcodeMate', 1, '--clip5pNbases', 14, '--soloBarcodeReadLength', 0]),
    ('inline-paired', ['--readFilesIn', work/'inline.fq', work/'mates.fq'], base + ['--soloBarcodeMate', 1, '--clip5pNbases', 14, 0, '--soloBarcodeReadLength', 0]),
    ('complex', ['--readFilesIn', work/'reads.fq', work/'barcodes.fq'],
        ['--soloType', 'CB_UMI_Complex', '--soloCBposition', '0_0_0_7', '--soloUMIposition', '0_8_0_13',
         '--soloCBwhitelist', work/'whitelist.txt', '--soloCBmatchWLtype', '1MM']),
    ('smartseq-se', ['--readFilesManifest', work/'manifest-se.tsv'], ['--soloType', 'SmartSeq', '--soloUMIdedup', 'Exact', 'NoDedup']),
    ('profile', ['--readFilesIn', work/'reads.fq', work/'barcodes.fq'], base),
    ('quoted-feature', ['--readFilesIn', work/'reads.fq', work/'barcodes.fq'], base + ['--soloOutFormatFeaturesGeneField3', 'RNA Expression']),
    ('smartseq', ['--readFilesManifest', work/'manifest.tsv'], ['--soloType', 'SmartSeq', '--soloUMIdedup', 'Exact', 'NoDedup']),
    ('two-pass', ['--readFilesIn', work/'reads.fq', work/'barcodes.fq', '--twopassMode', 'Basic'], base),
    ('sj', ['--readFilesIn', work/'reads.fq', work/'barcodes.fq'], base + ['--soloFeatures', 'Gene', 'SJ']),
]
checks = []
for case, inputs, solo in cases:
    for threads in (1, 4):
        name = f'adapter-{case}-{threads}'
        for exe, label in ((binary, 'native'), (star, 'star')):
            run(exe, ['--genomeDir', work/f'{label}-index', *inputs, '--runThreadN', threads,
                '--outFileNamePrefix', str(work/f'{label}-{name}') + '.', '--outSAMunmapped', 'Within',
                '--soloCellFilter', 'None', *solo, *(['--profile-stages'] if case=='profile' and label=='native' else [])], f'run-{label}-{name}')
        compare_dirs(work/f'native-{name}.Solo.out', work/f'star-{name}.Solo.out', name)
        # Headers contain executable paths; compare the actual SAM records.
        def sam(label):
            return sorted(line for line in (work/f'{label}-{name}.Aligned.out.sam').read_text().splitlines() if not line.startswith('@'))
        assert sam('native') == sam('star'), (name, 'SAM')
        checks.append(name)
        print('PASS', name, flush=True)
for method in [('TopCells', '3'), ('CellRanger2.2', '3', '0.99', '10')]:
    name = 'standalone-' + method[0]
    for exe, label in ((binary, 'native'), (star, 'star')):
        output = work / (label + '-' + name)
        output.mkdir(exist_ok=True)
        run(exe, ['--runMode', 'soloCellFiltering', work/'star-default-1.Solo.out/Gene/raw', str(output)+'/',
                  '--outFileNamePrefix', str(work / (label + '-' + name)) + '.', '--soloCellFilter', *method], label+'-'+name)
    compare_dirs(work/('native-'+name), work/('star-'+name), name)
    checks.append(name)
    print('PASS', name, flush=True)
# Exercise EmptyDrops' ambient model and simulation, rather than its early exit.
import random
rng = random.Random(719)
raw = work / 'ambient-raw'
raw.mkdir(exist_ok=True)
entries = []
for cb in range(200):
    counts = [0] * 50
    for _ in range(1000 if cb < 3 else 100 if cb < 30 else 20):
        counts[rng.randrange(5 if cb < 3 else 50)] += 1
    entries += [(g+1, cb+1, n) for g, n in enumerate(counts) if n]
(raw/'matrix.mtx').write_text('%%MatrixMarket matrix coordinate integer general\n' +
    f'50 200 {len(entries)}\n' + ''.join(f'{g} {cb} {n}\n' for g, cb, n in entries))
(raw/'features.tsv').write_text(''.join(f'g{i}\tgene{i}\tGene Expression\n' for i in range(50)))
(raw/'barcodes.tsv').write_text(''.join(f'cell{i}\n' for i in range(200)))
name = 'standalone-EmptyDrops_CR'
for exe, label in ((binary, 'native'), (star, 'star')):
    output = work / (label + '-' + name)
    output.mkdir(exist_ok=True)
    run(exe, ['--runMode', 'soloCellFiltering', raw, str(output)+'/',
              '--outFileNamePrefix', str(work/(label+'-'+name))+'.', '--soloCellFilter',
              'EmptyDrops_CR', 3, 0.99, 2, 50, 150, 50, 0.01, 20, 0.01, 100], label+'-'+name)
compare_dirs(work/('native-'+name), work/('star-'+name), name)
assert 'finished simulations' in (work/('native-'+name+'.Log.out')).read_text()
checks.append(name)
print('PASS', name, flush=True)
(work/'adapter-validation.json').write_text(json.dumps(checks, indent=2))
