#!/usr/bin/env python3
"""Compare a NuSTAR SAIS build with STAR's GenomeGenerate reference products."""

import argparse
import struct
from pathlib import Path


def nustar_genome(path: Path) -> bytes:
    data = path.read_bytes()
    magic, version, _, chromosomes, forward_size, junctions, synthetic_size, _ = (
        struct.unpack_from("<4sIIIQQQQ", data)
    )
    if magic != b"KREF" or version != 1:
        raise ValueError("unsupported NuSTAR reference format")
    offset = 48
    for _ in range(chromosomes):
        name_size = struct.unpack_from("<I", data, offset)[0]
        offset += 4 + name_size + 16
    forward = data[offset : offset + forward_size]
    offset += forward_size + junctions * 60
    synthetic = data[offset : offset + synthetic_size]
    if len(forward) != forward_size or len(synthetic) != synthetic_size:
        raise ValueError("truncated NuSTAR reference")
    return forward + synthetic


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--nustar", type=Path, required=True)
    args = parser.parse_args()
    failures = []
    products = [
        "chrName.txt", "chrLength.txt", "chrStart.txt", "chrNameLength.txt",
        "geneInfo.tab", "exonGeTrInfo.tab", "transcriptInfo.tab", "exonInfo.tab",
        "sjdbList.fromGTF.out.tab", "sjdbList.out.tab", "sjdbInfo.txt",
    ]
    for name in products:
        star_file, nustar_file = args.star / name, args.nustar / name
        if not star_file.exists() and not nustar_file.exists():
            continue
        if not star_file.exists() or not nustar_file.exists() or (
            star_file.read_bytes() != nustar_file.read_bytes()
        ):
            failures.append(name)
    if nustar_genome(args.nustar / "reference.bin") != (args.star / "Genome").read_bytes():
        failures.append("Genome sequence")
    if failures:
        print("DIFFERENT: " + ", ".join(failures))
        return 1
    print("MATCH: Genome sequence and all available chromosome/GTF/SJ tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
