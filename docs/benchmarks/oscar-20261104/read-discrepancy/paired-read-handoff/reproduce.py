#!/usr/bin/env python3
"""Replay the extracted pairs with existing engine-specific reference indexes."""
import argparse, json, pathlib, subprocess
p = argparse.ArgumentParser(description=__doc__)
for option in ('star', 'nustar', 'star-index', 'nustar-index', 'output'):
    p.add_argument('--' + option, type=pathlib.Path, required=True)
a = p.parse_args()
here = pathlib.Path(__file__).resolve().parent
metadata = json.loads((here / 'metadata.json').read_text())
for engine in ('star', 'nustar'):
    destination = a.output.resolve() / engine
    destination.mkdir(parents=True, exist_ok=False)
    command = list(metadata['commands'][engine])
    command[0] = str(getattr(a, engine).resolve())
    command[command.index('--genomeDir') + 1] = str(getattr(a, engine + '_index').resolve())
    i = command.index('--readFilesIn')
    command[i + 1:i + 3] = [str(here / 'R1.fastq'), str(here / 'R2.fastq')]
    command[command.index('--outFileNamePrefix') + 1] = str(destination / 'output.')
    print(command, flush=True)
    subprocess.run(command, check=True)
