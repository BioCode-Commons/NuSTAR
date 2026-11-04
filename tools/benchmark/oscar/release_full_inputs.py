#!/usr/bin/env python3
"""Release only this campaign's full FASTQs after every requested pair passed."""
import argparse,json
from pathlib import Path
from common import sha256,write_json

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
p.add_argument('--dataset',choices=['bulk','pbmc'],required=True)
p.add_argument('--job-id',required=True)
a=p.parse_args();root=a.root.resolve()
assert (root/'campaign.json').is_file(),'campaign marker missing'
suite=root/'results'/f'full-job-{a.job_id}'
status=json.loads((suite/'status.json').read_text())
if status['state']!='complete':raise RuntimeError('full suite has not passed')
metadata=json.loads((suite/'metadata.json').read_text())
record=metadata['inputs'][a.dataset]
folder=root/'prepared'/a.dataset/'full'
if json.loads((folder/'prepared.json').read_text())!=record:raise RuntimeError('input manifest changed')
for name,info in record['files'].items():
    path=folder/name
    if path.parent!=folder or path.name not in ['R1.fastq','R2.fastq']:raise RuntimeError('unowned cleanup path')
    if sha256(path)!=info['sha256']:raise RuntimeError('full input changed')
for name in record['files']:(folder/name).unlink()
write_json(folder/'released.json',{'validated_suite':str(suite),'inputs':record})
(folder/'prepared.json').unlink()
