"""Compare ordinary PE and peOverlap mapping with STAR on adversarial mate geometry.
Usage: paired_overlap.py NUSTAR STAR WORK_DIRECTORY
"""
from pathlib import Path
from collections import Counter
import json, random, subprocess, sys

binary, star, work = (Path(x).resolve() for x in sys.argv[1:4])
work.mkdir(parents=True, exist_ok=True)
programs = {'star': star, 'nustar': binary}

def run(command, name):
    with (work / (name + '.log')).open('w') as log:
        result = subprocess.run(list(map(str, command)), cwd=work, stdout=log, stderr=subprocess.STDOUT)
    assert result.returncode == 0, (name, (work / (name + '.log')).read_text()[-4000:])

def records(path):
    return Counter(line for line in path.read_text().splitlines() if not line.startswith('@'))

def stats(path):
    return [line for line in path.read_text().splitlines() if '|' in line and
            not any(word in line.lower() for word in ('time', 'speed', 'started', 'finished'))]

rng = random.Random(110527)
genome = list(''.join(rng.choice('ACGT') for _ in range(90000)))
for start in range(4000, 72001, 4000):
    for a, b in ((start+300, start+700), (start+1000, start+1500)):
        genome[a:a+2] = 'GT'; genome[b-2:b] = 'AG'
genome = ''.join(genome)
# Two identical loci exercise multimapping and primary alignment choices.
genome = genome[:80000] + genome[4000:8000] + genome[84000:]
(work/'genome.fa').write_text('>chr1\n'+genome+'\n>chr2\n'+genome[40000:48000]+'\n')
gtf = []
for number, start in enumerate(range(4000, 72001, 8000)):
    for left, right in ((start, start+300), (start+700, start+1000), (start+1500, start+1800)):
        gtf.append(f'chr1\ttest\texon\t{left+1}\t{right}\t.\t+\t.\tgene_id "g{number}"; transcript_id "t{number}";\n')
(work/'genes.gtf').write_text(''.join(gtf))
rc = lambda s: s.translate(str.maketrans('ACGT','TGCA'))[::-1]
geometries = [(0,150,60,150), (0,100,80,100), (0,100,99,100), (0,100,100,100),
              (0,100,180,100), (0,150,20,75), (30,75,0,150), (40,150,0,100),
              (0,100,0,100), (0,75,50,150), (0,150,140,100)]
pairs=[]; categories=Counter()
for start in range(4000, 72001, 4000):
    for kind, sequence in [('contiguous',genome[start:start+1800]),
            ('spliced',genome[start:start+300]+genome[start+700:start+1000]+genome[start+1500:start+1800])]:
        for geometry, (x, n0, y, n1) in enumerate(geometries):
            for change in ('exact','mismatch1','mismatch4','N','insertion','deletion','one-unmapped'):
                # Offset 240 makes overlap windows cross an annotated or novel splice.
                a=sequence[240+x:240+x+n0]; b=sequence[240+y:240+y+n1]
                if change=='mismatch1': b=b[:30]+('C' if b[30]!='C' else 'A')+b[31:]
                if change=='mismatch4':
                    for i in (5,20,35,50): b=b[:i]+('T' if b[i]!='T' else 'G')+b[i+1:]
                if change=='N': b=b[:25]+'NN'+b[27:]
                if change=='insertion': a=a[:45]+'A'+a[45:-1]
                if change=='deletion': a=a[:45]+a[46:]+sequence[240+x+n0]
                if change=='one-unmapped': b='N'*len(b)
                for reverse in (False,True):
                    pair=(a,rc(b)) if not reverse else (rc(b),a)
                    name=f'{kind}_{start}_g{geometry}_{change}_r{int(reverse)}'
                    pairs.append((name,*pair));categories[kind+'-'+change]+=1
pairs.extend([('all-N','N'*100,'N'*100), ('discordant',genome[1000:1100],rc(genome[70000:70100]))])
for mate in (1,2):
    (work/f'mate{mate}.fq').write_text(''.join(f'@{p[0]}\n{p[mate]}\n+\n'+('I'*len(p[mate]))+'\n' for p in pairs))
for label, exe in programs.items():
    run([exe,'--runMode','genomeGenerate','--genomeDir',work/(label+'-index'),
         '--genomeFastaFiles',work/'genome.fa','--sjdbGTFfile',work/'genes.gtf',
         '--sjdbOverhang',150,'--genomeChrBinNbits',18,'--genomeSAindexNbases',5,
         '--runThreadN',2], 'build-'+label)
checks=[]; baseline=None; changed=0
settings=[(0,.01,[], 'off'),(10,.01,[], 'on10'),(20,.05,[], 'on20-mm05'),
          (50,0,[], 'on50-exact'),(100,.1,[], 'on100-mm10'),
          (0,.01,['--clip5pNbases',3,7,'--clip3pNbases',5,2],'clipped-off'),
          (20,.05,['--clip5pNbases',3,7,'--clip3pNbases',5,2],'clipped'),
          (20,.05,['--twopassMode','Basic','--sjdbOverhang',150],'two-pass')]
for threads in (1,4):
    for minimum, mismatch, extra, case in settings:
        name=f'{case}-t{threads}'
        for label, exe in programs.items():
            run([exe,'--genomeDir',work/(label+'-index'),'--readFilesIn',work/'mate1.fq',work/'mate2.fq',
                 '--peOverlapNbasesMin',minimum,'--peOverlapMMp',mismatch,'--runThreadN',threads,
                 '--outSAMunmapped','Within','--outFileNamePrefix',str(work/(label+'-'+name))+'.',*extra],label+'-'+name)
        a,b=(str(work/(label+'-'+name))+'.' for label in programs)
        left,right=records(Path(a+'Aligned.out.sam')),records(Path(b+'Aligned.out.sam'))
        if left!=right:
            diff={'case':name,'star_only':list((left-right).elements())[:40], 'nustar_only':list((right-left).elements())[:40]}
            (work/'difference.json').write_text(json.dumps(diff,indent=2)+'\n')
            raise AssertionError((name,'SAM differs; see difference.json'))
        for suffix in ['SJ.out.tab','Log.final.out']+(['_STARpass1/SJ.out.tab','_STARpass1/Log.final.out'] if case=='two-pass' else []):
            read=stats if suffix.endswith('Log.final.out') else Path.read_bytes
            assert read(Path(a+suffix))==read(Path(b+suffix)),(name,suffix)
        if name=='off-t1': baseline=left
        if threads==1 and minimum and not extra: changed+=sum((left-baseline).values())
        checks.append(name)
        (work/'validation.json').write_text(json.dumps({'passed':checks,'pairs':len(pairs),'categories':categories,
            'changed_STAR_records_with_overlap':changed},indent=2)+'\n')
        print('PASS',name,flush=True)
assert changed>0, 'fixture did not exercise a change from peOverlap'
print('Exact PE/overlap SAM, junctions and non-time statistics:',len(checks),'cases',flush=True)
