#!/usr/bin/env python3
"""Build the same reference in STAR and NuSTAR, then compare sequence, tables and index."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def run(command: list[str], cwd: Path, label: str) -> None:
    log = cwd / f"{label}.command.log"
    with log.open("w") as output:
        completed = subprocess.run(command, cwd=cwd, stdout=output,
                                   stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-20:])
        raise RuntimeError(f"{label} failed ({completed.returncode}):\n{tail}")


def compare(root: Path, args: argparse.Namespace) -> None:
    star_dir, nustar_dir = root / "star", root / "nustar"
    star_dir.mkdir()
    nustar_dir.mkdir()
    common = ["--genomeFastaFiles", str(args.fasta),
              "--genomeChrBinNbits", str(args.chr_bin_bits),
              "--genomeSAindexNbases", str(args.prefix_bases)]
    if args.gtf:
        common += ["--sjdbGTFfile", str(args.gtf),
                   "--sjdbOverhang", str(args.overhang)]
    run([str(args.star), "--runMode", "genomeGenerate", "--genomeDir", str(star_dir),
         *common, "--runThreadN", str(args.star_threads),
         "--limitGenomeGenerateRAM", "4000000000"], star_dir, "star")
    run([str(args.nustar), "--runMode", "genomeGenerate", "--genomeDir", str(nustar_dir),
         *common], nustar_dir, "nustar")
    here = Path(__file__).resolve().parent
    subprocess.run(["python3", str(here / "compare_genome_generate.py"),
                    "--star", str(star_dir), "--nustar", str(nustar_dir)], check=True)
    subprocess.run([str(args.compare_index), str(star_dir), str(nustar_dir)], check=True)
    print(f"Differential comparison passed. STAR: {star_dir}; NuSTAR: {nustar_dir}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--nustar", type=Path, required=True)
    parser.add_argument("--compare-index", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--overhang", type=int, default=100)
    parser.add_argument("--chr-bin-bits", type=int, default=18)
    parser.add_argument("--prefix-bases", type=int, default=10)
    parser.add_argument("--star-threads", type=int, default=1)
    parser.add_argument("--work-dir", type=Path,
                        help="Keep artifacts here; directory must not already exist")
    args = parser.parse_args()
    args.star = args.star.resolve(strict=True)
    args.nustar = args.nustar.resolve(strict=True)
    args.compare_index = args.compare_index.resolve(strict=True)
    args.fasta = args.fasta.resolve(strict=True)
    if args.gtf:
        args.gtf = args.gtf.resolve(strict=True)
    if args.work_dir:
        root = args.work_dir.resolve()
        root.mkdir(parents=True, exist_ok=False)
        compare(root, args)
    else:
        with tempfile.TemporaryDirectory(prefix="nustar-genome-diff-") as temporary:
            compare(Path(temporary), args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
