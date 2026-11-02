"""Compare native two-pass mapping with STAR on novel/known junctions and mates."""
from pathlib import Path
import random, subprocess, sys, json
binary,star,work=map(lambda s:Path(s).resolve(),sys.argv[1:4])
work.mkdir(parents=True,exist_ok=True)
def run(exe,args,name):
    with (work/(name+'.log')).open('w') as log:
        p=subprocess.run([str(exe),*map(str,args)],stdout=log,stderr=subprocess.STDOUT,cwd=work)
    assert p.returncode==0,(name,(work/(name+'.log')).read_text()[-3000:])
def body(p):return sorted(s for s in p.read_text().splitlines() if not s.startswith('@'))
def stats(p):return [s for s in p.read_text().splitlines() if '|' in s and not any(t in s.lower() for t in ('time','speed','started','finished'))]
rng=random.Random(730)
genome=list(''.join(rng.choice('ACGT') for _ in range(16000)))
junctions=[(500,800),(2000,2600),(4000,4600),(7000,7600),(12000,12500)]
for start,end in junctions:genome[start:start+2]=list('GT');genome[end-2:end]=list('AG')
genome=''.join(genome)
(work/'genome.fa').write_text('>chr1\n'+genome+'\n')
(work/'genes.gtf').write_text('chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\nchr1\ttest\texon\t801\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
reads=[];mates=[]
for start,end in junctions:
    for anchor in (8,12,20,35,50,70,80,88,92):
        for repeat in range(5):
            s=genome[start-anchor:start]+genome[end:end+100-anchor]
            if repeat%2:s=rc(s)
            reads.append(s);mates.append(rc(genome[end+110:end+210]))
reads.extend([genome[10000:10100],'N'*100]);mates.extend([rc(genome[10200:10300]),'N'*100])
for filename,seqs in [('reads.fq',reads),('mates.fq',mates)]:
    (work/filename).write_text(''.join(f'@r{i}\n{s}\n+\n'+('I'*len(s))+'\n' for i,s in enumerate(seqs)))
checks=[]
for annotated in (False,True):
    for exe,label in ((binary,'native'),(star,'star')):
        index=work/f'{label}-idx-{annotated}'
        run(exe,['--runMode','genomeGenerate','--genomeDir',index,'--genomeFastaFiles',work/'genome.fa',
            '--genomeSAindexNbases',3,'--genomeChrBinNbits',8,'--runThreadN',2,
            *(['--sjdbGTFfile',work/'genes.gtf','--sjdbOverhang',99] if annotated else [])],f'index-{label}-{annotated}')
    for paired in (False,True):
      for threads,limit in ((1,-1),(1,40),(4,-1)):
        suffix=f'{annotated}-{paired}-{threads}-{limit}'
        for exe,label in ((binary,'native'),(star,'star')):
            prefix=str(work/f'{label}-{suffix}.')
            run(exe,['--genomeDir',work/f'{label}-idx-{annotated}','--readFilesIn',work/'reads.fq',
                *([work/'mates.fq'] if paired else []),'--runThreadN',threads,'--twopassMode','Basic',
                '--twopass1readsN',limit,'--sjdbOverhang',99,'--outSAMunmapped','Within',
                '--outSAMorder','Paired','--outFileNamePrefix',prefix],f'run-{label}-{suffix}')
        n=str(work/f'native-{suffix}.');s=str(work/f'star-{suffix}.')
        for name in ['Aligned.out.sam','SJ.out.tab','Log.final.out','_STARpass1/SJ.out.tab','_STARpass1/Log.final.out','_STARgenome/sjdbInfo.txt','_STARgenome/sjdbList.out.tab']:
            a,b=Path(n+name),Path(s+name)
            reader=body if name.endswith('.sam') else stats if name.endswith('Log.final.out') else lambda p:p.read_text().splitlines()
            assert reader(a)==reader(b),(suffix,name)
        assert not Path(n+'_STARpass1/Aligned.out.sam').exists()
        checks.append(suffix)
(work/'validation.json').write_text(json.dumps({'passed':checks,'reads':len(reads)},indent=2))
print('Exact SAM, junctions, non-time statistics and inserted junction tables:',len(checks),'two-pass cases')
# Persisted pass-2 indices must reload with the new reference identity.
for exe,label in ((binary,'native'),(star,'star')):
    prefix=str(work/f'{label}-saved.')
    run(exe,['--genomeDir',work/f'{label}-idx-True','--readFilesIn',work/'reads.fq',
             '--twopassMode','Basic','--sjdbInsertSave','All','--outSAMunmapped','Within',
             '--outFileNamePrefix',prefix],f'save-{label}')
    run(exe,['--genomeDir',prefix+'_STARgenome/','--readFilesIn',work/'reads.fq',
             '--outSAMunmapped','Within','--outFileNamePrefix',str(work/f'{label}-reload.')],f'reload-{label}')
    assert body(Path(prefix+'Aligned.out.sam'))==body(work/f'{label}-reload.Aligned.out.sam')
assert body(work/'native-reload.Aligned.out.sam')==body(work/'star-reload.Aligned.out.sam')
print('Saved pass-2 indices reload with identical SAM')
