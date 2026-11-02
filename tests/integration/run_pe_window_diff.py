#!/usr/bin/env python3
"""Compare paired-read seed -> Window rows against the project STAR oracle."""

import argparse
import itertools
import os
from pathlib import Path
import subprocess
import tempfile
from compare_paired_stitch import compare as compare_stitch


KINDS = ("STITCH_READ\t", "STITCH_WIN\t", "STITCH_WA\t")


def sample_fastq(source: Path, target: Path, count: int) -> None:
    with source.open() as inp, target.open("w") as out:
        for _ in range(count):
            record = [inp.readline() for _ in range(4)]
            if any(not line for line in record) or not record[0].startswith("@"):
                raise ValueError(f"FASTQ has fewer than {count} reads: {source}")
            out.writelines(record)


def star_rows(path: Path):
    with path.open() as source:
        for row in source:
            if row.startswith(KINDS):
                yield row


def compare(star_trace: Path, nustar_trace: Path) -> int:
    with nustar_trace.open() as actual:
        count = 0
        for count, (wanted, got) in enumerate(
            itertools.zip_longest(star_rows(star_trace), actual, fillvalue=""), 1
        ):
            if wanted != got:
                raise AssertionError(f"PE Window row {count} differs:\n"
                                     f"STAR:   {wanted.rstrip()}\n"
                                     f"NuSTAR: {got.rstrip()}")
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-oracle", type=Path, required=True)
    parser.add_argument("--nustar-trace", type=Path, required=True)
    parser.add_argument("--stitch-trace", type=Path)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--nustar-index", type=Path, required=True)
    parser.add_argument("--mate1", type=Path, required=True)
    parser.add_argument("--mate2", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=1000)
    parser.add_argument("--work-dir", type=Path)
    args = parser.parse_args()
    if args.sample < 1:
        parser.error("--sample must be positive")

    def run(root: Path) -> None:
        root.mkdir(parents=True, exist_ok=True)
        first, second = root / "mate1.fastq", root / "mate2.fastq"
        sample_fastq(args.mate1, first, args.sample)
        sample_fastq(args.mate2, second, args.sample)
        star_stderr = root / "star.stderr"
        with (root / "star.stdout").open("w") as stdout, star_stderr.open("w") as stderr:
            subprocess.run([
                str(args.star_oracle.resolve()), "--runMode", "alignReads",
                "--runThreadN", "1", "--genomeDir", str(args.star_index.resolve()),
                "--readFilesIn", str(first), str(second),
                "--outFileNamePrefix", str(root / "star."), "--outSAMtype", "None",
            ], stdout=stdout, stderr=stderr, check=True,
               env=dict(os.environ, STAR_STITCH_TRACE="1"))
        nustar_trace = root / "nustar.windows.tsv"
        subprocess.run([
            str(args.nustar_trace.resolve()), str(args.nustar_index.resolve()),
            str(first), str(second), str(nustar_trace),
        ], check=True, stdout=subprocess.DEVNULL)
        count = compare(star_stderr, nustar_trace)
        print(f"PE Window exact: {args.sample} pairs, {count} rows")
        if args.stitch_trace:
            stitch_output = root / "nustar.stitch.tsv"
            subprocess.run([
                str(args.stitch_trace.resolve()),
                str(args.nustar_index.resolve()), str(first), str(second),
                str(star_stderr), str(stitch_output),
            ], check=True, stdout=subprocess.DEVNULL)
            stitch_rows = compare_stitch(star_stderr, stitch_output)
            print(f"PE stitching exact: {args.sample} pairs, {stitch_rows} rows "
                  "(time and uninitialized mate-gap shifts excluded)")

    if args.work_dir:
        run(args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="nustar-pe-window-") as temp:
            run(Path(temp))


if __name__ == "__main__":
    main()
