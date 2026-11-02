#!/usr/bin/env python3
"""Strict SAM differential checker for STAR and NuSTAR.

All body columns, optional-field order, multiplicity, and record order are
compared by default. Use --unordered only when the run's cross-read order is
allowed to vary; it still compares every complete SAM line as a multiset.
"""

import argparse
import collections
import itertools
import sys
from pathlib import Path


FIELDS = ("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT",
          "PNEXT", "TLEN", "SEQ", "QUAL")


def read_sam(path: Path) -> tuple[list[str], list[str]]:
    header: list[str] = []
    records: list[str] = []
    with path.open() as source:
        for line_number, line in enumerate(source, 1):
            line = line.rstrip("\n")
            if line.startswith("@") and not records:
                header.append(line)
                continue
            if len(line.split("\t")) < 11:
                raise ValueError(f"{path}:{line_number}: fewer than 11 SAM fields")
            records.append(line)
    return header, records


def compare_ordered_sam(star_path: Path, nustar_path: Path,
                        header_mode: str = "sq") -> int:
    """Compare ordered SAM records while retaining only one line per file."""
    def read_header(source):
        header = []
        while line := source.readline():
            if line.startswith("@"):
                header.append(line.rstrip("\n"))
            else:
                return header, line.rstrip("\n")
        return header, None

    with star_path.open() as star_file, nustar_path.open() as nustar_file:
        star_header, star_first = read_header(star_file)
        nustar_header, nustar_first = read_header(nustar_file)
        if header_mode == "sq":
            star_header = [line for line in star_header if line.startswith("@SQ\t")]
            nustar_header = [line for line in nustar_header if line.startswith("@SQ\t")]
        if header_mode != "none" and star_header != nustar_header:
            first = next((number for number, (left, right) in enumerate(
                itertools.zip_longest(star_header, nustar_header), 1)
                if left != right), None)
            raise AssertionError(f"SAM header mismatch at line {first}")
        star_records = itertools.chain(
            () if star_first is None else (star_first,),
            (line.rstrip("\n") for line in star_file))
        nustar_records = itertools.chain(
            () if nustar_first is None else (nustar_first,),
            (line.rstrip("\n") for line in nustar_file))
        count = 0
        for count, (star, nustar) in enumerate(itertools.zip_longest(
                star_records, nustar_records), 1):
            if star is None or nustar is None or star != nustar:
                raise AssertionError(
                    f"SAM body differs at record {count}: "
                    f"{describe_difference(star or '', nustar or '')}")
            if star.count("\t") < 10:
                raise ValueError(f"SAM body has fewer than 11 fields at record {count}")
        return count


def describe_difference(left: str, right: str) -> str:
    a = left.split("\t")
    b = right.split("\t")
    for column, (x, y) in enumerate(itertools.zip_longest(a, b), 1):
        if x != y:
            label = FIELDS[column - 1] if column <= 11 else f"TAG[{column - 12}]"
            return f"{label}: STAR={x!r}, NuSTAR={y!r}"
    return "record differs"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("star", type=Path)
    parser.add_argument("nustar", type=Path)
    parser.add_argument("--unordered", action="store_true")
    parser.add_argument("--header", choices=("none", "sq", "all"), default="sq")
    args = parser.parse_args()
    if not args.unordered:
        try:
            count = compare_ordered_sam(args.star, args.nustar, args.header)
        except (AssertionError, ValueError) as error:
            print(error, file=sys.stderr)
            return 1
        print(f"SAM exact: {count} records (ordered); header={args.header}")
        return 0
    star_header, star = read_sam(args.star)
    nustar_header, nustar = read_sam(args.nustar)
    if args.header == "sq":
        star_header = [line for line in star_header if line.startswith("@SQ\t")]
        nustar_header = [line for line in nustar_header if line.startswith("@SQ\t")]
    if args.header != "none" and star_header != nustar_header:
        print("SAM header mismatch", file=sys.stderr)
        for number, (left, right) in enumerate(itertools.zip_longest(
                star_header, nustar_header), 1):
            if left != right:
                print(f"first header difference at line {number}:\n"
                      f"STAR: {left}\nNuSTAR: {right}", file=sys.stderr)
                break
        return 1
    if args.unordered:
        missing = collections.Counter(star) - collections.Counter(nustar)
        extra = collections.Counter(nustar) - collections.Counter(star)
        if missing or extra:
            print(f"SAM differs: STAR={len(star)} NuSTAR={len(nustar)} "
                  f"missing={sum(missing.values())} extra={sum(extra.values())}",
                  file=sys.stderr)
            for title, counter in (("missing", missing), ("extra", extra)):
                for line, count in counter.most_common(3):
                    print(f"{title} {count}x: {line}", file=sys.stderr)
            return 1
    print(f"SAM exact: {len(star)} records "
          f"({'unordered' if args.unordered else 'ordered'}); "
          f"header={args.header}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
