"""Write a synthetic genome, GTF and spliced single-end reads for two-pass timing.

The genome is random sequence, so it has none of a real genome's repeats. It is
meant to time the junction-insertion step at a realistic suffix-array size and
junction count, not to measure alignment accuracy. Half of the generated genes
go into the GTF; reads come from all genes, so pass 1 also finds unannotated
junctions that pass 2 must insert into the index.

usage: make_twopass_dataset.py OUT_DIR [--chromosomes 4] [--chromosome-mb 50]
                               [--reads 1000000] [--seed 1]
"""
import argparse
import os
import random
from pathlib import Path

BASES = bytes((b'ACGT'[i % 4] for i in range(256)))


def random_sequence(rng, length):
    return bytearray(rng.randbytes(length).translate(BASES))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('out', type=Path)
    parser.add_argument('--chromosomes', type=int, default=4)
    parser.add_argument('--chromosome-mb', type=int, default=50)
    parser.add_argument('--reads', type=int, default=1_000_000)
    parser.add_argument('--read-length', type=int, default=100)
    parser.add_argument('--seed', type=int, default=1)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    transcripts = []  # (chromosome index, [(exon start, exon end)], annotated), 0-based half-open
    genome = []
    gtf = []
    for c in range(args.chromosomes):
        seq = random_sequence(rng, args.chromosome_mb * 1_000_000)
        pos, gene = 20_000, 0
        while True:
            exons = []
            start = pos
            for _ in range(rng.randint(3, 12)):
                length = rng.randint(80, 300)
                exons.append((start, start + length))
                start += length + rng.randint(300, 6000)
            if exons[-1][1] + 20_000 > len(seq):
                break
            for (_, end), (next_start, _) in zip(exons, exons[1:]):
                seq[end:end + 2] = b'GT'
                seq[next_start - 2:next_start] = b'AG'
            annotated = gene % 2 == 0
            transcripts.append((c, exons, annotated))
            if annotated:
                for a, b in exons:
                    gtf.append(f'chr{c + 1}\tsynthetic\texon\t{a + 1}\t{b}\t.\t+\t.\t'
                               f'gene_id "g{c + 1}_{gene}"; transcript_id "t{c + 1}_{gene}";\n')
            pos = exons[-1][1] + rng.randint(2_000, 20_000)
            gene += 1
        genome.append(seq)
    with (args.out / 'genome.fa').open('wb') as fa:
        for c, seq in enumerate(genome):
            fa.write(f'>chr{c + 1}\n'.encode())
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i + 60] + b'\n')
    (args.out / 'genes.gtf').write_text(''.join(gtf))

    complement = bytes.maketrans(b'ACGT', b'TGCA')
    n = args.read_length
    with (args.out / 'reads.fq').open('wb') as fq:
        quality = b'I' * n
        for r in range(args.reads):
            c, exons, _ = rng.choice(transcripts)
            mrna = b''.join(genome[c][a:b] for a, b in exons)
            offset = rng.randrange(len(mrna) - n)
            read = bytearray(mrna[offset:offset + n])
            for _ in range(rng.choice((0, 0, 1, 2))):
                read[rng.randrange(n)] = rng.choice(b'ACGT')
            if rng.random() < 0.5:
                read = read.translate(complement)[::-1]
            fq.write(b'@r%d\n%s\n+\n%s\n' % (r, read, quality))
    introns = sum(len(e) - 1 for _, e, _ in transcripts)
    annotated = sum(len(e) - 1 for _, e, a in transcripts if a)
    print(f'{len(transcripts)} genes, {introns} introns ({annotated} annotated), '
          f'{args.chromosomes * args.chromosome_mb} Mb, {args.reads} reads')


if __name__ == '__main__':
    main()
