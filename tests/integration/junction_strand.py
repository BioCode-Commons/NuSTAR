"""Compare conflicting GTF junction strands with STAR, including paired mapping.
Usage: junction_strand.py NUSTAR STAR WORK_DIRECTORY
"""
from collections import Counter
from pathlib import Path
import itertools
import json
import random
import subprocess
import sys

binary, star, work = (Path(x).resolve() for x in sys.argv[1:4])
work.mkdir(parents=True, exist_ok=True)
rng = random.Random(61105)
sequence = list(''.join(rng.choice('ACGT') for _ in range(160000)))
annotations, pairs, cases = [], [], []
strands = [('+', '-'), ('-', '+'), ('.', '+', '-'), ('-', '.', '+'),
           ('+', '.'), ('-', '.'), ('.',), ('+',), ('-',)]
for number, (motif, choices) in enumerate(itertools.product(
        ['GTAG', 'CTAC', 'GCAG', 'CTGC', 'ATAC', 'GTAT', 'AAAA'], strands)):
    start = 1000 + number * 2400
    a, b, c, d = start + 180, start + 500, start + 680, start + 1000
    sequence[a:a+2], sequence[b-2:b] = motif[:2], motif[2:]
    sequence[a-1] = 'C' if motif[-1] == 'A' else 'A'
    sequence[b] = 'A' if motif[0] != 'A' else 'C'
    reverse = motif in ('CTAC', 'CTGC', 'GTAT')
    second_motif = 'CTAC' if reverse else 'GTAG'
    sequence[c:c+2], sequence[d-2:d] = second_motif[:2], second_motif[2:]
    gene = f'g{number}'
    for t, strand in enumerate(choices):
        for left, right in [(start, a), (b, c)]:
            annotations.append(f'chr1\ttest\texon\t{left+1}\t{right}\t.\t{strand}\t.\t'
                               f'gene_id "{gene}"; transcript_id "{gene}_{t}";\n')
    for left, right in [(b, c), (d, d+180)]:
        annotations.append(f'chr1\ttest\texon\t{left+1}\t{right}\t.\t{"-" if reverse else "+"}\t.\t'
                           f'gene_id "{gene}"; transcript_id "{gene}_second";\n')
    first = ''.join(sequence[a-45:a] + sequence[b:b+56])
    second = ''.join(sequence[c-45:c] + sequence[d:d+56])
    second = second.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    pairs.append((gene, first, second))
    cases.append({'name': gene, 'motif': motif, 'strands': choices})

(work / 'genome.fa').write_text('>chr1\n' + ''.join(sequence) + '\n')
(work / 'genes.gtf').write_text(''.join(annotations))
for mate in (1, 2):
    (work / f'R{mate}.fastq').write_text(''.join(
        f'@{p[0]}\n{p[mate]}\n+\n'+ 'I'*len(p[mate])+'\n' for p in pairs))

commands = []
def run(command):
    command = list(map(str, command))
    with (work / f'command-{len(commands)}.log').open('w') as log:
        result = subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT)
    commands.append({'command': command, 'exit': result.returncode})
    (work / 'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    if result.returncode:
        raise RuntimeError(f'command {len(commands)-1} failed')

for name, exe in [('star', star), ('nustar', binary)]:
    index = work / (name + '-index')
    index.mkdir(exist_ok=True)
    run([exe, '--runMode', 'genomeGenerate', '--genomeDir', index,
         '--genomeFastaFiles', work/'genome.fa', '--sjdbGTFfile', work/'genes.gtf',
         '--sjdbOverhang', 100, '--genomeChrBinNbits', 18,
         '--genomeSAindexNbases', 6, '--runThreadN', 1])
    run([exe, '--genomeDir', index, '--readFilesIn', work/'R1.fastq', work/'R2.fastq',
         '--runThreadN', 1, '--outSAMunmapped', 'Within',
         '--outFileNamePrefix', str(work / (name + '.'))])

checks = {}
for filename in ('sjdbInfo.txt', 'sjdbList.out.tab'):
    checks[filename] = (work/'star-index'/filename).read_bytes() == (work/'nustar-index'/filename).read_bytes()

# STAR also collapses adjacent equal coordinates at a strand-group boundary
# in its first pass. Keep this behavior for a reference with one duplicate SJ.
(work/'isolated.gtf').write_text(''.join(row for row in annotations
    if 'gene_id "g9";' in row and 'g9_second' not in row))
for name, exe in [('star', star), ('nustar', binary)]:
    index = work / (name + '-isolated')
    index.mkdir(exist_ok=True)
    run([exe, '--runMode', 'genomeGenerate', '--genomeDir', index,
         '--genomeFastaFiles', work/'genome.fa', '--sjdbGTFfile', work/'isolated.gtf',
         '--sjdbOverhang', 100, '--genomeChrBinNbits', 18,
         '--genomeSAindexNbases', 6, '--runThreadN', 1])
for filename in ('sjdbInfo.txt', 'sjdbList.out.tab'):
    checks['isolated-' + filename] = (work/'star-isolated'/filename).read_bytes() == (work/'nustar-isolated'/filename).read_bytes()

def records(name):
    return Counter(tuple(f[:11]+sorted(f[11:]))
                   for line in (work/(name+'.Aligned.out.sam')).read_text().splitlines()
                   if not line.startswith('@') for f in [line.split('\t')])
checks['SAM'] = records('star') == records('nustar')
checks['SJ.out.tab'] = (work/'star.SJ.out.tab').read_bytes() == (work/'nustar.SJ.out.tab').read_bytes()
def stats(name):
    return [line for line in (work/(name+'.Log.final.out')).read_text().splitlines()
            if '|' in line and not any(s in line.lower() for s in ('started', 'finished', 'speed'))]
checks['statistics'] = stats('star') == stats('nustar')
(work/'validation.json').write_text(json.dumps({'cases': cases, 'checks': checks}, indent=2)+'\n')
if not all(checks.values()):
    raise RuntimeError(f'STAR junction strand comparison failed: {checks}')
print(f'Passed {len(cases)} motif/strand combinations: junction metadata, paired SAM, SJ and statistics')
