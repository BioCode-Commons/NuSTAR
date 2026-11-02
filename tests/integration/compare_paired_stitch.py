#!/usr/bin/env python3
"""Compare paired stitching traces, excluding time and STAR uninitialized bytes."""

import argparse
import itertools
from pathlib import Path

from run_stitch_diff import RECORD_TYPES


def rows(path: Path):
    with path.open() as source:
        for row in source:
            if row.startswith(RECORD_TYPES):
                yield row.rstrip("\n")


def normalized(row: str) -> str:
    fields = row.split("\t")
    if fields[0] == "STITCH_OUT":
        fields[-1] = "<elapsed-ns>"
    if fields[0] == "STITCH_JUNCTION" and fields[5] == "-3":
        # STAR leaves shiftSJ uninitialized for a mate gap; those bytes are
        # never a splice-junction output field and vary between processes.
        fields[-2:] = ["<mate-gap>", "<mate-gap>"]
    return "\t".join(fields)


def compare(star_path: Path, nustar_path: Path) -> int:
    count = 0
    for count, (star, nustar) in enumerate(
        itertools.zip_longest(rows(star_path), rows(nustar_path), fillvalue=""), 1
    ):
        if normalized(star) != normalized(nustar):
            raise AssertionError(f"paired stitch row {count} differs:\n"
                                 f"STAR:   {star}\nNuSTAR: {nustar}")
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("star", type=Path)
    parser.add_argument("nustar", type=Path)
    args = parser.parse_args()
    count = compare(args.star, args.nustar)
    print(f"paired stitching exact: {count} rows (excluding time/mate-gap garbage)")


if __name__ == "__main__":
    main()
