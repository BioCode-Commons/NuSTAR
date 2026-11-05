"""Local integrity checks using real files; no aligner is simulated or timed."""
import gzip
import io
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest

from common import Capacity, compare_sam_streams, same_file, validate_pair, sha256
from prepare_inputs import pairs
from run_suite import command_for, measure, execution_policy


SAM=b'@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\tNH:i:1\n'

class Integrity(unittest.TestCase):
    def test_explicit_stage_policy(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            self.assertIsNone(execution_policy(root))
            path=root/'execution-policy.json'
            policy={'schema':1,'enabled_stages':['main','full'],'reason':'User requested skipping pilot'}
            path.write_text(json.dumps(policy))
            result=execution_policy(root)
            self.assertEqual(result['policy'],policy)
            self.assertEqual(result['sha256'],sha256(path))
            self.assertNotIn('pilot',result['policy']['enabled_stages'])
            for stages in [[],['main','main'],['unknown'],'main']:
                path.write_text(json.dumps({**policy,'enabled_stages':stages}))
                with self.assertRaises(ValueError):execution_policy(root)

    def test_ordered_records_headers_and_tags(self):
        result=compare_sam_streams(io.BytesIO(SAM),io.BytesIO(SAM.replace(b'VN:1.6',b'VN:1.5')))
        self.assertEqual(result['records'],1)
        for changed in [SAM.replace(b'NH:i:1',b'NH:i:2'),SAM.replace(b'LN:1000',b'LN:1001'),SAM+SAM.splitlines(True)[-1]]:
            with self.assertRaises(AssertionError):compare_sam_streams(io.BytesIO(SAM),io.BytesIO(changed))

    def test_truncated_record_rejected(self):
        with self.assertRaises(ValueError):compare_sam_streams(io.BytesIO(b'bad\n'),io.BytesIO(b'bad\n'))

    def test_fastq_mates_and_truncation(self):
        with tempfile.TemporaryDirectory() as temp:
            left,right=Path(temp)/'left.gz',Path(temp)/'right.gz'
            def put(path,data):
                with gzip.open(path,'wb') as out:out.write(data)
            put(left,b'@r1/1\nACGT\n+\nIIII\n');put(right,b'@r1/2\nTGCA\n+\nIIII\n')
            self.assertEqual(len(list(pairs([(left,right)]))),1)
            put(right,b'@r2/2\nTGCA\n+\nIIII\n')
            with self.assertRaises(ValueError):list(pairs([(left,right)]))
            put(right,b'@r1/2\nTGCA\n+\nII')
            with self.assertRaises(ValueError):list(pairs([(left,right)]))

    def test_exact_files_and_capacity(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);a=root/'a';b=root/'b';a.write_bytes(b'ABC');b.write_bytes(b'ABC')
            self.assertEqual(same_file(a,b)['bytes'],3)
            b.write_bytes(b'ABD')
            with self.assertRaises(AssertionError):same_file(a,b)
            with self.assertRaises(RuntimeError):Capacity(root,1,0).check()

    def test_validation_keeps_real_difference(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for label in ['star','nustar']:
                (root/f'{label}.SJ.out.tab').write_text('chr1\t1\t2\n')
                (root/f'{label}.Log.final.out').write_text('Number of input reads | 1\nStarted job on | clock\n')
            result=validate_pair(str(root/'star.'),str(root/'nustar.'),{'two_pass':False,'solo':False,'bam':False},None,root)
            self.assertEqual(result['statistics'],{'Number of input reads':'1'})
            (root/'nustar.SJ.out.tab').write_text('chr1\t1\t3\n')
            with self.assertRaises(AssertionError):validate_pair(str(root/'star.'),str(root/'nustar.'),{'two_pass':False,'solo':False,'bam':False},None,root)
            self.assertTrue((root/'nustar.SJ.out.tab').is_file())

    def test_process_failure_preserves_measurement(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with self.assertRaises(RuntimeError):
                measure([sys.executable,'-c','raise SystemExit(7)'],root,Capacity(root,1<<20,0))
            self.assertEqual(json.loads((root/'timing.json').read_text())['returncode'],7)

    def test_full_cleanup_requires_completed_validation(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'campaign.json').write_text('{}')
            folder=root/'prepared/bulk/full';folder.mkdir(parents=True)
            files={}
            for name in ['R1.fastq','R2.fastq']:
                path=folder/name;path.write_bytes(b'@r1\nACGT\n+\nIIII\n')
                files[name]={'bytes':path.stat().st_size,'sha256':sha256(path)}
            record={'pairs':1,'files':files};(folder/'prepared.json').write_text(json.dumps(record))
            (folder/'keep.txt').write_text('retained')
            suite=root/'results/full-job-check';suite.mkdir(parents=True)
            (suite/'metadata.json').write_text(json.dumps({'inputs':{'bulk':record}}))
            (suite/'status.json').write_text(json.dumps({'state':'failed'}))
            command=[sys.executable,str(Path(__file__).parent/'release_full_inputs.py'),'--root',str(root),'--dataset','bulk','--job-id','check']
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)
            self.assertTrue((folder/'R1.fastq').exists())
            (suite/'status.json').write_text(json.dumps({'state':'complete'}))
            subprocess.run(command,check=True,capture_output=True)
            self.assertFalse((folder/'R1.fastq').exists())
            self.assertTrue((folder/'keep.txt').exists())
            self.assertTrue((folder/'released.json').exists())

    def test_same_biological_arguments(self):
        manifest=json.loads((Path(__file__).parent/'suite.json').read_text())
        for case in manifest['cases']:
            inputs=[Path(p.format(stage='main')) for p in case['inputs']]
            a=command_for(Path('STAR'),Path('star-index'),inputs,case,8,'star.',Path('white'))
            b=command_for(Path('nustar'),Path('nustar-index'),inputs,case,8,'nustar.',Path('white'))
            for args in [a,b]:
                args[0]='binary';args[args.index('--genomeDir')+1]='index';args[args.index('--outFileNamePrefix')+1]='prefix'
            self.assertEqual(a,b)
            if case['solo']:self.assertTrue(str(inputs[0]).endswith('R2.fastq'))

if __name__=='__main__':unittest.main()
