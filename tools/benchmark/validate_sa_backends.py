#!/usr/bin/env python3
"""Index-only SA-backend correctness, construction time and memory validation.

Each timed command runs alone under benchmark.lock. Each backend/configuration
is built once using the full test or human reference. Human uses 32 threads
exclusively. No FASTQ input, read mapping or SAM output is used.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while block := f.read(8 << 20): h.update(block)
    return h.hexdigest()


def index_header(f):
    def number(fmt):
        b=f.read(struct.calcsize(fmt))
        if len(b)!=struct.calcsize(fmt): raise RuntimeError('truncated index')
        return struct.unpack(fmt,b)[0]
    if f.read(4)!=b'KSAI': raise RuntimeError('invalid index magic')
    fields={k:number(fmt) for k,fmt in [('version','<I'),('prefix','<I'),('limit','<I'),
                                      ('physical','<Q'),('count','<Q'),('prefix_count','<Q')]}
    for k in ['reference_semantic','options_hash','algorithm']:
        n=number('<I')
        if n>4096: raise RuntimeError('invalid index string')
        fields[k]=f.read(n).decode()
    fields['width']=number('<I')
    if fields['version']!=4: raise RuntimeError('expected v4 index')
    return fields


def index_evidence(directory):
    with (directory/'index.sa').open('rb') as f:
        header=index_header(f);h=hashlib.sha256()
        while b:=f.read(8<<20):h.update(b)
    content={k:v for k,v in header.items() if k not in ['algorithm','options_hash']}
    return {'header':header,'semantic_header':content,'payload_sha256':h.hexdigest(),
            'reference_sha256':sha(directory/'reference.bin'),
            'sj':{f.name:sha(f) for f in directory.glob('sjdb*.tab')} |
                 {'sjdbInfo.txt':sha(directory/'sjdbInfo.txt')}}


def same_index(a,b):
    return all(a[k]==b[k] for k in ['semantic_header','payload_sha256','reference_sha256','sj'])


def exact_payload(left,right):
    with (left/'index.sa').open('rb') as a,(right/'index.sa').open('rb') as b:
        index_header(a);index_header(b)
        while True:
            x=a.read(8<<20);y=b.read(8<<20)
            if x!=y: raise RuntimeError('exact packed SA/SAi comparison failed')
            if not x:return


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nustar',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--human-backends',nargs='+',choices=['libsais64','libsais40'],
                        default=['libsais64','libsais40'])
    args=parser.parse_args();root=args.root.resolve();work=args.work.resolve()
    if 'libsais64' not in args.human_backends or len(set(args.human_backends))!=len(args.human_backends):
        parser.error('human-backends must include libsais64 and contain no duplicates')
    work.mkdir(parents=True,exist_ok=True)
    if (work/'validation.json').exists():raise SystemExit('Refusing to overwrite existing evidence')
    program=work/'nustar';shutil.copy2(args.nustar,program)
    script=work/'validate_sa_backends.py'
    if Path(__file__).resolve()!=script:shutil.copy2(__file__,script)
    state={'status':'waiting for cleanup/benchmark lock','runs':[],'checks':[],
           'mode':'index-only','program_sha256':sha(program),'script_sha256':sha(script),
           'rounds':1,'first_round_passed':False,'removed':[]}
    input_paths=[root/'data/test'/name for name in ['test.fa','test.gtf']]
    input_paths += [root/'data/references/human'/name for name in ['GRCh38.primary_assembly.genome.fa','gencode.v50.primary_assembly.annotation.gtf']]
    state['inputs']={str(p):{'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in input_paths}
    def save():
        tmp=work/'validation.tmp.json';tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(work/'validation.json')
        lines=['# SA backend validation','',f"Status: {state['status']}; current: {state.get('current','waiting for benchmark lock')}.",'',
               'Index construction only, using the full test and human references. No read mapping or SAM output. Human builds use 32 threads exclusively.',
               'Each backend/configuration is built once. Index comparison and cleanup are excluded from timings.',
               'The default backend remains libsais64. Memory acceptance and timing results use this single construction round.','',
               '| Run | Seconds | Peak RSS GiB |','| --- | ---: | ---: |']
        for r in state['runs']:lines.append(f"| {r['name']} | {r['seconds']:.3f} | {r['rss_kib']/2**20:.3f} |")
        if 'summary' in state:
            lines+=['','| Configuration | Algorithm | Seconds | Peak RSS GiB |','| --- | --- | ---: | ---: |']
            for r in state['summary']:lines.append(f"| {r['case']} | {r['algorithm']} | {r['seconds']:.3f} | {r['rss_gib']:.3f} |")
        if 'error' in state:lines+=['','Stopped because: '+state['error']]
        lines+=['','See `validation.json` for individual checks and cleanup records. Artifacts with mismatches are retained.']
        (work/'report.md').write_text('\n'.join(lines)+'\n')
    def run(name,command):
        state.update(status='running',current=name);save()
        started=time.monotonic()
        with (work/(name+'.console')).open('wb') as f:
            p=subprocess.Popen(list(map(str,command)),stdout=f,stderr=subprocess.STDOUT)
            _,status,usage=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(status)
        record={'name':name,'command':list(map(str,command)),'exit':p.returncode,
                'seconds':time.monotonic()-started,'rss_kib':usage.ru_maxrss,
                'user_seconds':usage.ru_utime,'system_seconds':usage.ru_stime}
        state['runs'].append(record);save()
        if p.returncode:raise RuntimeError('command failed: '+name)
        return record
    def remove(path,reason):
        if path.parent!=work or path.is_symlink():raise RuntimeError('unexpected cleanup path')
        entry={'path':str(path),'reason':reason,'deleted':False};state['removed'].append(entry);save()
        if path.is_dir():shutil.rmtree(path)
        else:path.unlink()
        entry['deleted']=True;save()
    def build(case,dataset,threads,algorithm,rnd):
        name=f'r{rnd}-{case}-{algorithm}';directory=work/name;directory.mkdir()
        fasta=root/('data/test/test.fa' if dataset=='test' else 'data/references/human/GRCh38.primary_assembly.genome.fa')
        gtf=root/('data/test/test.gtf' if dataset=='test' else 'data/references/human/gencode.v50.primary_assembly.annotation.gtf')
        r=run('build-'+name,[program,'--runMode','genomeGenerate','--genomeFastaFiles',fasta,
              '--sjdbGTFfile',gtf,'--sjdbOverhang',150,'--limitSjdbInsertNsj',2000000,
              '--genomeChrBinNbits',18,'--genomeSAindexNbases',10 if dataset=='test' else 14,
              '--genomeDir',directory,'--runThreadN',threads,'--genomeSAalgorithm',algorithm,'--profileGenomeGenerate'])
        state['current']='verify-index-'+name;save()
        r['index']=index_evidence(directory);save()
        return directory,r
    save()
    try:
        with (root/'GitHubAction/benchmark.lock').open('r') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            audit=root/'NuSTAR/build/defensive-check-cleanup/full-sam/cleanup-audit.json'
            if audit.exists() and json.loads(audit.read_text())['status']!='complete':
                raise RuntimeError('previous SAM cleanup has not completed successfully')
            cases=[('test-t1','test',1,['libsais64','libsais40']),
                   ('test-t32','test',32,['libsais64','libsais40']),
                   ('human-t32','human',32,['libsais64']+[a for a in args.human_backends if a!='libsais64'])]
            first={}
            for case,data,threads,algorithms in cases:
                built={a:build(case,data,threads,a,1) for a in algorithms}
                baseline=built['libsais64']
                for a in algorithms[1:]:
                    state['current']='compare-index-'+case+'-'+a;save()
                    if not same_index(baseline[1]['index'],built[a][1]['index']):raise RuntimeError('index mismatch: '+case+' '+a)
                    exact_payload(baseline[0],built[a][0]);state['checks'].append({'type':'exact-index','case':case,'algorithm':a,'equal':True});save()
                first[case]={a:r['index'] for a,(_,r) in built.items()}
                if case=='test-t32':
                    for a in algorithms:
                        if not same_index(first['test-t1'][a],first[case][a]):
                            raise RuntimeError('index differs between 1 and 32 threads: '+a)
                        state['checks'].append({'type':'index-thread-invariance','algorithm':a,
                                                'threads':[1,32],'equal':True});save()
                for directory,_ in built.values():remove(directory,'complete index verified; comparison evidence retained')
            state['first_round_passed']=True;save()
            summary=[]
            for case,_,_,algorithms in cases:
                for a in algorithms:
                    r=next(r for r in state['runs'] if r['name']==f'build-r1-{case}-{a}')
                    summary.append({'case':case,'algorithm':a,'seconds':r['seconds'],
                                    'rss_gib':r['rss_kib']/2**20})
            state['summary']=summary
            human={r['algorithm']:r for r in summary if r['case']=='human-t32'}
            for algorithm in human.keys()-{'libsais64'}:
                if human[algorithm]['rss_gib']>=human['libsais64']['rss_gib']:
                    raise RuntimeError('human peak memory did not improve: '+algorithm)
            state.update(status='complete',current='finished');save()
    except BaseException as e:
        state['status']='failed';state['error']=repr(e);save();raise

if __name__=='__main__':main()
