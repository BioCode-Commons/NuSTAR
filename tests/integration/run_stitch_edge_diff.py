#!/usr/bin/env python3
"""Build tiny annotated indexes and compare SE stitching on edge-case reads."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--nustar", type=Path, required=True)
    parser.add_argument("--compare-index", type=Path, required=True)
    parser.add_argument("--star-oracle", type=Path, required=True)
    parser.add_argument("--nustar-trace", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    for name in ("star", "nustar", "compare_index", "star_oracle",
                 "nustar_trace", "fasta", "gtf", "reads"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    here = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="nustar-stitch-edge-") as temporary:
        root = Path(temporary)
        index = root / "index"
        subprocess.run([
            "python3", str(here / "run_genome_generate_diff.py"),
            "--star", str(args.star), "--nustar", str(args.nustar),
            "--compare-index", str(args.compare_index),
            "--fasta", str(args.fasta), "--gtf", str(args.gtf),
            "--overhang", "16", "--chr-bin-bits", "8", "--prefix-bases", "4",
            "--work-dir", str(index),
        ], cwd=root, check=True)
        comparison = [
            "python3", str(here / "run_stitch_diff.py"),
            "--star-oracle", str(args.star_oracle),
            "--nustar-trace", str(args.nustar_trace),
            "--star-index", str(index / "star"),
            "--nustar-index", str(index / "nustar"),
            "--reads", str(args.reads), "--sample", "511",
        ]
        if args.json:
            comparison.extend(["--json", str(args.json.resolve())])
        subprocess.run(comparison, cwd=root, check=True)


if __name__ == "__main__":
    main()
