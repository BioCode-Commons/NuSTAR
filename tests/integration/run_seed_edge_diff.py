#!/usr/bin/env python3
"""Build tiny annotated STAR/NuSTAR indexes and compare edge-case seed traces."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--nustar", type=Path, required=True)
    parser.add_argument("--star-oracle", type=Path, required=True)
    parser.add_argument("--nustar-trace", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    for name in ("star", "nustar", "star_oracle", "nustar_trace",
                 "fasta", "gtf", "reads"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    here = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="nustar-seed-edge-") as temporary:
        root = Path(temporary)
        star_index = root/"star-index"
        nustar_index = root/"nustar-index"
        star_index.mkdir()
        subprocess.run([
            str(args.star), "--runMode", "genomeGenerate",
            "--genomeDir", str(star_index),
            "--genomeFastaFiles", str(args.fasta),
            "--sjdbGTFfile", str(args.gtf), "--sjdbOverhang", "16",
            "--genomeChrBinNbits", "8", "--genomeSAindexNbases", "4",
            "--limitGenomeGenerateRAM", "4000000000",
        ], cwd=root, stdout=subprocess.DEVNULL, check=True)
        subprocess.run([
            str(args.nustar), "genomeGenerate",
            "--genomeDir", str(nustar_index),
            "--genomeFastaFiles", str(args.fasta),
            "--sjdbGTFfile", str(args.gtf), "--sjdbOverhang", "16",
            "--genomeChrBinNbits", "8", "--genomeSAindexNbases", "4",
        ], cwd=root, stdout=subprocess.DEVNULL, check=True)
        comparison = [
            "python3", str(here/"run_seed_diff.py"),
            "--star-oracle", str(args.star_oracle),
            "--nustar-trace", str(args.nustar_trace),
            "--star-index", str(star_index),
            "--nustar-index", str(nustar_index),
            "--mate1", str(args.reads), "--sample", "511",
        ]
        if args.json:
            comparison.extend(["--json", str(args.json.resolve())])
        subprocess.run(comparison, cwd=root, check=True)


if __name__ == "__main__":
    main()
