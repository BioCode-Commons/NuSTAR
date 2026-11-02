"""Differential STARsolo tests: exact feature matrices, barcodes and statistics."""
from pathlib import Path
import random, subprocess, sys, json
binary,star,work=map(lambda s:Path(s).resolve(),sys.argv[1:4])
work.mkdir(parents=True,exist_ok=True)
def run(exe,args,name):
    with (work/(name+'.log')).open('w') as log:
        p=subprocess.run([str(exe),*map(str,args)],stdout=log,stderr=subprocess.STDOUT,cwd=work)
    assert p.returncode==0,(name,p.returncode,(work/(name+'.log')).read_text()[-3000:])
rng=random.Random(741)
genome=''.join(rng.choice('ACGT') for _ in range(30000))
# Two genes share their complete sequence, so multi-gene UMI redistribution is exercised.
genome=genome[:20500]+genome[18500:20500]+genome[22500:]
(work/'genome.fa').write_text('>chr1\n'+genome+'\n')
gtf=[]
for g in range(12):
    start=500+g*2000
    for left,right in [(start,start+500),(start+800,start+1400)]:
        gtf.append(f'chr1\ttest\texon\t{left+1}\t{right}\t.\t+\t.\tgene_id "g{g}"; transcript_id "t{g}"; gene_name "gene{g}"; gene_type "test";\n')
(work/'genes.gtf').write_text(''.join(gtf))
cbs=['AACCGGTT','CCGGAATT','GGCCTTAA','TTCCAAGG','AGCTTCGA','CGATAGCT','AACCGGTA']
(work/'whitelist.txt').write_text('\n'.join(cbs)+'\n')
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
reads=[];barcodes=[];mates=[]
for i in range(2400):
    g=rng.randrange(12); start=500+g*2000+rng.randrange(400)
    seq=genome[start:start+100]
    if i%17==0:seq=genome[500+g*2000+550:500+g*2000+650] # intronic
    if i%19==0:
        end=500+g*2000+500; second=500+g*2000+800
        seq=genome[end-50:end]+genome[second:second+50]
    if i%23==0:seq='N'*100
    if i%29==0:seq=rc(seq)
    cb=cbs[rng.randrange(len(cbs))]
    umi=['ACGTAC','ACGTAT','ACGTTC','CATGCA','GTACGT','AAAAAA','CCCCCC','GGGGGG','TTTTTT','NNACGT'][rng.randrange(10)]
    if i%13==0:cb='T'+cb[1:]
    if i%31==0:cb='N'+cb[1:]
    if i%37==0:cb=cb[:-1]+'N'
    if i%47==0:cb='NN'+cb[2:]
    reads.append(seq);barcodes.append(cb+umi);mates.append(rc(genome[start+160:start+260]))
for filename,seqs in [('reads.fq',reads),('barcodes.fq',barcodes),('mates.fq',mates)]:
    (work/filename).write_text(''.join(f'@r{i}\n{s}\n+\n'+('I'*len(s))+'\n' for i,s in enumerate(seqs)))
for exe,label in ((binary,'native'),(star,'star')):
    run(exe,['--runMode','genomeGenerate','--genomeDir',work/f'{label}-index','--genomeFastaFiles',work/'genome.fa',
        '--genomeSAindexNbases',3,'--genomeChrBinNbits',8,'--sjdbOverhang',99,'--runThreadN',2,
        '--sjdbGTFfile',work/'genes.gtf'],f'index-{label}')
cases=[('default',[]),('gene-full',['--soloFeatures','Gene','GeneFull','GeneFull_Ex50pAS','GeneFull_ExonOverIntron']),
 ('dedup',['--soloUMIdedup','Exact','1MM_All','1MM_Directional','1MM_CR']),
 ('exact',['--soloCBmatchWLtype','Exact']),
 ('multi',['--soloCBmatchWLtype','1MM_multi_Nbase_pseudocounts','--soloMultiMappers','Uniform','Rescue','PropUnique','EM']),
 ('filtered',['--soloCellFilter','TopCells','3']),
 ('velocity',['--soloFeatures','Gene','Velocyto']),
 ('reverse',['--soloStrand','Reverse']),('unstranded',['--soloStrand','Unstranded']),
 ('cell-read-stats',['--soloCellReadStats','Standard']),
 ('cr-filter',['--soloUMIdedup','1MM_CR','--soloUMIfiltering','MultiGeneUMI_CR'])]
checks=[]
for name,extra in cases:
 for threads in [1,4]:
    suffix=f'{name}-{threads}'
    for exe,label in ((binary,'native'),(star,'star')):
        prefix=str(work/f'{label}-{suffix}.')
        run(exe,['--genomeDir',work/f'{label}-index','--readFilesIn',work/'reads.fq',work/'barcodes.fq',
            '--runThreadN',threads,'--outFileNamePrefix',prefix,'--outSAMtype','None',
            '--soloType','CB_UMI_Simple','--soloCBstart',1,'--soloCBlen',8,'--soloUMIstart',9,'--soloUMIlen',6,
            '--soloCBwhitelist',work/'whitelist.txt',*(['--soloCellFilter','None'] if '--soloCellFilter' not in extra else []),*extra],f'run-{label}-{suffix}')
    n=work/f'native-{suffix}.Solo.out';s=work/f'star-{suffix}.Solo.out'
    a={str(p.relative_to(n)):p.read_bytes() for p in n.rglob('*') if p.is_file()}
    b={str(p.relative_to(s)):p.read_bytes() for p in s.rglob('*') if p.is_file()}
    assert a.keys()==b.keys(),(suffix,'files',a.keys(),b.keys())
    for f in a:
        assert a[f]==b[f],(suffix,f)
    checks.append(suffix)
(work/'validation.json').write_text(json.dumps({'passed':checks,'reads':len(reads)},indent=2))
print('Exact Solo outputs:',len(checks),'cases')
