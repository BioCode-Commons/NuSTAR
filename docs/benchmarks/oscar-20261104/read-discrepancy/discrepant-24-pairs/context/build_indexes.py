#!/usr/bin/env python3
"""Build both native indexes from identical verified reference inputs."""
import argparse
import json
import os
from pathlib import Path
from common import Capacity, sha256, write_json
from run_suite import measure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--threads',type=int,required=True)
    parser.add_argument('--disk-budget-gib',type=int,required=True)
    args=parser.parse_args()
    root=args.root.resolve()
    if not (root/'campaign.json').is_file() or not os.environ.get('SLURM_JOB_ID'):
        parser.error('initialized Slurm campaign required')
    capacity=Capacity(root,args.disk_budget_gib<<30,32<<30)
    reference=root/'reference'
    identity={'fasta_sha256':sha256(reference/'genome.fa'),'gtf_sha256':sha256(reference/'genes.gtf'),
              'sjdbOverhang':100,'genomeSAindexNbases':14,'genomeChrBinNbits':18,'threads':args.threads}
    for label,binary in [('star',root/'source/STAR/source/STAR'),('nustar',root/'build/nustar/nustar')]:
        index=root/'index'/label
        index.mkdir(parents=True,exist_ok=True)
        marker=index/'validated-build.json'
        expected={**identity,'binary_sha256':sha256(binary)}
        if marker.exists():
            if json.loads(marker.read_text())['identity']!=expected:
                raise RuntimeError('existing index was built from different inputs')
            continue
        if any(index.iterdir()):
            raise RuntimeError(f'incomplete index retained at {index}; diagnose before removing it')
        log=root/'logs'/f'index-{label}'
        log.mkdir(exist_ok=True)
        command=[str(binary),'--runMode','genomeGenerate','--genomeDir',str(index),
                 '--genomeFastaFiles',str(reference/'genome.fa'),'--sjdbGTFfile',str(reference/'genes.gtf'),
                 '--sjdbOverhang','100','--genomeSAindexNbases','14','--genomeChrBinNbits','18',
                 '--runThreadN',str(args.threads)]
        result=measure(command,log,capacity)
        write_json(marker,{'identity':expected,'measurement':result,
                           'files':{p.name:{'bytes':p.stat().st_size,'sha256':sha256(p)} for p in index.iterdir() if p.is_file()}})


if __name__=='__main__':main()
