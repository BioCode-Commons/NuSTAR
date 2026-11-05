#!/usr/bin/env python3
"""Replay a selected discrepant pair, or print its command for a debugger."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess


def main():
    root = Path(__file__).resolve().parent
    pairs = {row['pair']: row for row in json.loads((root / 'pairs.json').read_text())}
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pair', choices=sorted(pairs))
    parser.add_argument('--engine', choices=('star', 'nustar', 'both'), default='both')
    for option in ('star', 'nustar', 'star-index', 'nustar-index'):
        parser.add_argument('--' + option, type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    engines = ('star', 'nustar') if args.engine == 'both' else (args.engine,)
    inputs = root / 'per-read' / args.pair if args.pair else root
    metadata = json.loads((root / 'metadata.json').read_text())
    for engine in engines:
        binary, index = getattr(args, engine), getattr(args, engine + '_index')
        if binary is None or index is None:
            parser.error(f'--{engine} and --{engine}-index are required for this engine')
        command = list(metadata['commands'][engine])
        command[0] = str(binary.resolve())
        command[command.index('--genomeDir') + 1] = str(index.resolve())
        offset = command.index('--readFilesIn')
        command[offset + 1:offset + 3] = [str(inputs / 'R1.fastq'), str(inputs / 'R2.fastq')]
        destination = args.output.resolve() / engine
        command[command.index('--outFileNamePrefix') + 1] = str(destination / 'output.')
        print(shlex.join(command), flush=True)
        if not args.dry_run:
            destination.mkdir(parents=True, exist_ok=False)
            subprocess.run(command, check=True)


if __name__ == '__main__':
    main()
