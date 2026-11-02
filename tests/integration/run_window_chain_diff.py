#!/usr/bin/env python3
"""Compare the new seed -> window -> stitch chain with saved STAR runtime rows.

The STAR rows are produced by build_star_stitch_oracle.py, from this project's
STAR source copy. This script keeps all replay code in validation, not modules.
"""

import argparse
import itertools
import subprocess
import tempfile
from pathlib import Path


def compare_rows(expected: Path, actual: Path, ignore_window_time: bool) -> int:
    count = 0
    with expected.open() as left, actual.open() as right:
        for count, (wanted, got) in enumerate(
            itertools.zip_longest(left, right, fillvalue=""), start=1
        ):
            if ignore_window_time and wanted.startswith("STITCH_OUT\t"):
                wanted = "\t".join(wanted.rstrip("\n").split("\t")[:-1])
                got = "\t".join(got.rstrip("\n").split("\t")[:-1])
            if wanted != got:
                raise AssertionError(
                    f"{expected.name} vs {actual.name}, line {count}:\n"
                    f"STAR:   {wanted.rstrip()}\nNuSTAR: {got.rstrip()}"
                )
    return count


def make_stitch_replay(windows: Path, target: Path) -> None:
    current_read = None
    current_window = None
    with windows.open() as source, target.open("w") as output:
        def close_window() -> None:
            if current_window is not None:
                output.write(f"STITCH_OUT\t{current_read}\t{current_window}\t0\t0\n")

        def close_read() -> None:
            if current_read is not None:
                output.write(f"STITCH_BEST\t{current_read}\t0\t0\t0\t0\t0\t0\n")

        for row in source:
            kind, name, *_ = row.split("\t")
            if kind == "STITCH_READ":
                close_window()
                close_read()
                current_window = None
                current_read = name
            elif kind == "STITCH_WIN":
                close_window()
                current_window = row.split("\t")[2]
            output.write(row)
        close_window()
        close_read()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--star-windows", type=Path, required=True)
    parser.add_argument("--star-stitch", type=Path, required=True)
    parser.add_argument("--window-trace", type=Path, required=True)
    parser.add_argument("--stitch-trace", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path)
    args = parser.parse_args()

    def run(directory: Path) -> None:
        windows = directory / "nustar.windows.tsv"
        playback = directory / "nustar.window_replay.tsv"
        stitch = directory / "nustar.stitch.tsv"
        subprocess.run([args.window_trace, args.index, args.reads, windows], check=True)
        window_rows = compare_rows(args.star_windows, windows, False)
        make_stitch_replay(windows, playback)
        subprocess.run([args.stitch_trace, args.index, args.reads, playback, stitch],
                       check=True)
        stitch_rows = compare_rows(args.star_stitch, stitch, True)
        print(f"exact windows: {window_rows} rows; exact stitching: {stitch_rows} rows")

    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=True)
        run(args.work_dir)
    else:
        with tempfile.TemporaryDirectory(prefix="nustar-window-diff-") as directory:
            run(Path(directory))


if __name__ == "__main__":
    main()
