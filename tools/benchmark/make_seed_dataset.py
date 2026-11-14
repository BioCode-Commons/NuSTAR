#!/usr/bin/env python3
"""Generate deterministic synthetic inputs for seed timing and differential tests.

This is a stress/coverage dataset, not a substitute for real RNA-seq benchmarks.
The default writes 8 Mb of sequence and 100,000 paired reads under work/seed-data.
Both mates can also be used independently for single-end tests. The JSON manifest
records the recipe and SHA-256 of every input; no expected alignment is asserted.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random

BASES = b"ACGT"
RANDOM_BASES = bytes(BASES[i % 4] for i in range(256))
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")
# A complete cycle covers common long matches and difficult interval searches.
CASES = (['exact'] * 35 + ['mismatch'] * 15 + ['insertion'] * 5 +
         ['deletion'] * 5 + ['splice'] * 15 + ['repeat'] * 10 +
         ['low_complexity'] * 5 + ['ambiguous'] * 5 + ['short'] * 3 +
         ['unmapped'] * 2)


def random_sequence(rng, length):
    return rng.randbytes(length).translate(RANDOM_BASES)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def positive(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return result


def write_fastq(handle, name, sequence):
    handle.write(b'@' + name + b'\n' + sequence + b'\n+\n' +
                 b'I' * len(sequence) + b'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out', type=Path, nargs='?', default=Path('work/seed-data'))
    parser.add_argument('--genome-size', type=positive, default=8_000_000,
                        help='total biological bases, before index padding')
    parser.add_argument('--chromosomes', type=positive, default=2)
    parser.add_argument('--reads', type=positive, default=100_000,
                        help='read pairs; each FASTQ contains this many records')
    parser.add_argument('--read-length', type=positive, default=150)
    parser.add_argument('--seed', type=int, default=2711)
    args = parser.parse_args()
    n = args.read_length
    if n < 30 or n > 1000:
        parser.error('--read-length must be between 30 and 1000')
    if args.genome_size // args.chromosomes < max(100_000, n * 100):
        parser.error('each chromosome must have at least max(100000, read-length*100) bases')
    rng = random.Random(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    fragment_length = 2 * n + 100
    repeat_length = max(1024, fragment_length + 32)
    families = [random_sequence(rng, repeat_length) for _ in range(16)]
    genomes, repeat_regions, low_regions, transcripts, gtf = [], [], [], [], []
    for chrom in range(args.chromosomes):
        size = args.genome_size // args.chromosomes + (chrom < args.genome_size % args.chromosomes)
        sequence = bytearray(random_sequence(rng, size))
        # The first fifth is partitioned into repeated and low-complexity blocks.
        # Genes are placed later, so motif edits cannot destroy these blocks.
        repeat_limit = size // 6
        for copy, start in enumerate(range(1000, repeat_limit - repeat_length, repeat_length)):
            family = copy % len(families)
            sequence[start:start + repeat_length] = families[family]
            repeat_regions.append((chrom, start, repeat_length))
        low_start, low_end = size // 6, size // 5
        for copy, start in enumerate(range(low_start, low_end - repeat_length, repeat_length)):
            motif = (b'A', b'AC', b'ACGT', b'GATTACA')[copy % 4]
            sequence[start:start + repeat_length] = (motif * (repeat_length // len(motif) + 1))[:repeat_length]
            low_regions.append((chrom, start, repeat_length))
        # Isolated N runs cover chromosome-adjacent and interior prefix fallbacks.
        sequence[:31] = b'N' * 31
        sequence[-31:] = b'N' * 31
        for start in range(size // 5 + 200, size - 1000, 9973):
            sequence[start:start + 23] = b'N' * 23
        exon_length = max(450, 2 * n)
        gene_span = 3 * exon_length + 2 * 300
        spacing = gene_span + 2500
        for gene, start in enumerate(range(size // 4, size - gene_span - 100, spacing)):
            exons = [(start + i * (exon_length + 300),
                      start + i * (exon_length + 300) + exon_length) for i in range(3)]
            minus = gene % 2 == 1
            # Clean transcript sequence and set canonical motifs for each strand.
            for a, b in exons:
                sequence[a:b] = random_sequence(rng, b - a)
            for (_, donor), (acceptor, _) in zip(exons, exons[1:]):
                sequence[donor:donor + 2] = b'CT' if minus else b'GT'
                sequence[acceptor - 2:acceptor] = b'AC' if minus else b'AG'
            annotated = gene % 4 != 3
            transcripts.append((chrom, exons, minus, annotated))
            if annotated:
                for a, b in exons:
                    gtf.append(f'chr{chrom + 1}\tsynthetic\texon\t{a + 1}\t{b}\t.\t'
                               f'{"-" if minus else "+"}\t.\tgene_id "g{chrom}_{gene}"; '
                               f'transcript_id "t{chrom}_{gene}";\n')
        genomes.append(sequence)
    if not repeat_regions or not low_regions or not transcripts:
        parser.error('genome too small for the selected read length')
    fasta = args.out / 'genome.fa'
    with fasta.open('wb') as handle:
        for chrom, sequence in enumerate(genomes):
            handle.write(f'>chr{chrom + 1}\n'.encode())
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + b'\n')
    (args.out / 'genes.gtf').write_text(''.join(gtf))
    counts = Counter()
    with (args.out / 'mate1.fastq').open('wb') as first, (args.out / 'mate2.fastq').open('wb') as second:
        for ordinal in range(args.reads):
            case = CASES[ordinal % len(CASES)]
            counts[case] += 1
            if case in ('repeat', 'low_complexity'):
                chrom, begin, size = rng.choice(repeat_regions if case == 'repeat' else low_regions)
                start = begin + rng.randrange(size - fragment_length + 1)
                fragment = bytes(genomes[chrom][start:start + fragment_length])
            elif case == 'splice':
                chrom, exons, minus, annotated = rng.choice(transcripts)
                transcript = b''.join(genomes[chrom][a:b] for a, b in exons)
                if minus:
                    transcript = reverse_complement(transcript)
                # One mate crosses a junction, with both long and short anchors.
                junction = exons[0][1] - exons[0][0]
                start = junction - rng.randint(5, n - 5)
                fragment = transcript[start:start + fragment_length]
                counts['splice_annotated' if annotated else 'splice_unannotated'] += 1
            elif case == 'unmapped':
                fragment = random_sequence(rng, fragment_length)
            else:
                chrom = rng.randrange(len(genomes))
                start = rng.randrange(len(genomes[chrom]) - fragment_length + 1)
                fragment = bytes(genomes[chrom][start:start + fragment_length])
            if rng.randrange(2):
                fragment = reverse_complement(fragment)
                counts['fragment_reverse'] += 1
            mates = [bytearray(fragment[:n]), bytearray(reverse_complement(fragment[-n:]))]
            for mate in mates:
                if case == 'mismatch':
                    for position in rng.sample(range(n), rng.randint(1, 4)):
                        base = mate[position]
                        mate[position] = rng.choice([value for value in BASES if value != base])
                elif case == 'insertion':
                    position = rng.randint(5, n - 5)
                    mate[position:position] = random_sequence(rng, rng.randint(1, 3))
                elif case == 'deletion':
                    position = rng.randint(5, n - 5)
                    del mate[position:position + rng.randint(1, 3)]
                elif case == 'ambiguous':
                    if ordinal % 1000 < 100:
                        mate[:] = b'N' * len(mate)
                    else:
                        for position in rng.sample(range(n), rng.randint(1, 6)):
                            mate[position] = ord('N')
                elif case == 'short':
                    del mate[rng.randint(5, 24):]
            name = f'r{ordinal:09d}_{case}'.encode()
            write_fastq(first, name, bytes(mates[0]))
            write_fastq(second, name, bytes(mates[1]))
    files = {}
    for name in ('genome.fa', 'genes.gtf', 'mate1.fastq', 'mate2.fastq'):
        path = args.out / name
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        files[name] = {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}
    manifest = {
        'kind': 'synthetic seed-search stress data; not representative real RNA-seq',
        'recipe_version': 1,
        'seed': args.seed,
        'genome_size': args.genome_size,
        'chromosomes': args.chromosomes,
        'read_pairs': args.reads,
        'nominal_read_length': n,
        'suggested_sjdb_overhang': n - 1,
        'transcripts': len(transcripts),
        'annotated_transcripts': sum(annotated for _, _, _, annotated in transcripts),
        'counts': dict(sorted(counts.items())),
        'files': files,
    }
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'output': str(args.out.resolve()), 'read_pairs': args.reads,
                      'genome_size': args.genome_size, 'counts': manifest['counts']}, indent=2))


if __name__ == '__main__':
    main()
