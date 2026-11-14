#!/usr/bin/env python3
"""Download a bounded real yeast RNA-seq subset and its Ensembl reference.

SRR453566 is the yeast RNA-seq accession used in the Babraham RNA-seq course.
The first N paired records are a deterministic prefix subset, not a random sample.
Only the gzip prefix necessary for N records is transferred from each ENA file.
"""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess


def fetch(url):
    return subprocess.check_output(['curl', '--fail', '--location', '--silent',
                                    '--show-error', '--connect-timeout', '30',
                                    '--max-time', '180', url])


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fastq_prefix(url, path, count):
    proc = subprocess.Popen(['curl', '--fail', '--location', '--silent',
                             '--show-error', '--connect-timeout', '30',
                             '--max-time', '300', url], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL)
    lengths = {}
    names = []
    try:
        with gzip.GzipFile(fileobj=proc.stdout) as src, path.open('wb') as dst:
            for i in range(count):
                record = [src.readline() for _ in range(4)]
                if not record[0].startswith(b'@') or not record[2].startswith(b'+'):
                    raise RuntimeError('incomplete FASTQ at record %d: %s' % (i, url))
                length = len(record[1].rstrip())
                if length != len(record[3].rstrip()):
                    raise RuntimeError('sequence/quality mismatch')
                dst.writelines(record)
                names.append(record[0].split()[0].removesuffix(b'/1').removesuffix(b'/2'))
                lengths[length] = lengths.get(length, 0) + 1
    finally:
        if proc.poll() is None:
            proc.terminate()
        proc.wait()
        proc.stdout.close()
    return {'url': url, 'records': count, 'length_histogram': lengths,
            'bytes': path.stat().st_size, 'sha256': digest(path)}, names


def reference(url, path):
    data = gzip.decompress(fetch(url))
    path.write_bytes(data)
    return {'url': url, 'bytes': len(data), 'sha256': digest(path)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent / 'yeast')
    parser.add_argument('--reads', type=int, default=50000)
    args = parser.parse_args()
    if args.reads < 1 or args.reads > 1000000:
        parser.error('--reads must be 1..1000000')
    args.out.mkdir(parents=True, exist_ok=True)
    api = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR453566'
           '&result=read_run&fields=run_accession,scientific_name,library_strategy,'
           'library_layout,fastq_ftp,fastq_bytes,fastq_md5&format=tsv')
    metadata = fetch(api).decode()
    (args.out / 'ena_metadata.tsv').write_text(metadata)
    run = next(csv.DictReader(io.StringIO(metadata), delimiter='\t'))
    urls = ['https://' + item for item in run['fastq_ftp'].split(';')]
    mates = [[url for url in urls if url.endswith('_%d.fastq.gz' % i)][0] for i in (1, 2)]
    base = 'https://ftp.ensembl.org/pub/release-112/'
    fasta = base + 'fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz'
    gtf = base + 'gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz'
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        tasks = [pool.submit(fastq_prefix, url, args.out / ('reads_%d.fastq' % i), args.reads)
                 for i, url in enumerate(mates, 1)]
        fa_task = pool.submit(reference, fasta, args.out / 'reference.fa')
        gtf_task = pool.submit(reference, gtf, args.out / 'annotation.gtf')
        first, names1 = tasks[0].result()
        second, names2 = tasks[1].result()
        if names1 != names2:
            raise RuntimeError('paired read identifiers differ')
        fa, annotation = fa_task.result(), gtf_task.result()
    contigs = {}
    name = None
    for line in (args.out / 'reference.fa').read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            contigs[name] = 0
        elif name:
            contigs[name] += len(line)
    gtf_contigs = {line.split('\t')[0] for line in (args.out / 'annotation.gtf').read_text().splitlines()
                   if line and not line.startswith('#')}
    if not gtf_contigs <= contigs.keys():
        raise RuntimeError('reference and annotation contigs differ')
    fa.update({'total_bases': sum(contigs.values()), 'contigs': contigs})
    result = {'accession': 'SRR453566', 'subset': 'first %d paired records' % args.reads,
              'course': 'https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf',
              'study': 'https://www.ncbi.nlm.nih.gov/sra/SRR453566',
              'metadata_api': api, 'reads_1': first, 'reads_2': second,
              'reference': fa, 'annotation': annotation,
              'note': 'CENPK RNA-seq aligned to standard S288C R64-1-1 reference, as a teaching/benchmark dataset; native NuSTAR index must be built locally.'}
    (args.out / 'provenance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
