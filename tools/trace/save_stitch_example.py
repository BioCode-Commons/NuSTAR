#!/usr/bin/env python3
"""Save a readable, complete 1000-read STAR/NuSTAR SE stitching example."""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import tempfile

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "tests/integration"))

from run_stitch_diff import compare, run, sample_fastq


WINDOW_TYPES = ("STITCH_READ\t", "STITCH_WIN\t", "STITCH_WA\t")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_oracle_trace(source: Path, windows: Path, candidates: Path) -> None:
    with source.open() as inp, windows.open("w") as win, candidates.open("w") as out:
        for line in inp:
            if line.startswith(WINDOW_TYPES):
                win.write(line)
            elif line.startswith(("STITCH_OUT\t", "STITCH_TR\t", "STITCH_META\t",
                                  "STITCH_EXON\t", "STITCH_JUNCTION\t", "STITCH_BEST\t")):
                out.write(line)


def save_summary(reads: Path, trace: Path, star_index: Path,
                 destination: Path, chromosomes: Path) -> None:
    names = (star_index / "chrName.txt").read_text().splitlines()
    starts = [int(x) for x in (star_index / "chrStart.txt").read_text().splitlines()]
    lengths = [int(x) for x in (star_index / "chrLength.txt").read_text().splitlines()]
    if len(names) != len(lengths) or len(starts) < len(names):
        raise ValueError("STAR chromosome tables have incompatible lengths")
    with chromosomes.open("w") as out:
        out.write("chromosome_id\tname\tglobal_start_0based\tlength\n")
        for i, (name, start, length) in enumerate(zip(names, starts, lengths)):
            out.write(f"{i}\t{name}\t{start}\t{length}\n")
    best = {}
    candidate_counts = Counter()
    window_counts = Counter()
    with trace.open() as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "STITCH_OUT":
                window_counts[fields[1]] += 1
            elif fields[0] == "STITCH_TR":
                candidate_counts[fields[1]] += 1
            elif fields[0] == "STITCH_BEST":
                best[fields[1]] = fields[2:]
    with reads.open() as inp, destination.open("w") as out:
        out.write("read\tbest_score\tchromosome\tposition_1based\tstrand"
                  "\tglobal_start_0based\tgenome_span\taccepted_windows"
                  "\tvisited_windows\tcandidates\n")
        for header in inp:
            name = header.rstrip("\n").split()[0]
            for _ in range(3):
                if not inp.readline():
                    raise ValueError("truncated FASTQ")
            if name not in best:
                raise ValueError(f"missing STITCH_BEST for {name}")
            score, global_start, genome_span, chromosome_id, strand, accepted = best[name]
            chromosome = int(chromosome_id)
            if accepted == "0":
                chromosome_name, position, strand_name = "*", "0", "*"
            else:
                chromosome_name = names[chromosome]
                position = str(int(global_start) - starts[chromosome] + 1)
                strand_name = "+" if strand == "0" else "-"
            out.write("\t".join([name, score, chromosome_name, position,
                                  strand_name, global_start, genome_span,
                                  accepted, str(window_counts[name]),
                                  str(candidate_counts[name])]) + "\n")


def verify_sam_reads(reads: Path, sam: Path) -> int:
    with reads.open() as inp:
        read_names = {line.split()[0][1:] for index, line in enumerate(inp)
                      if index % 4 == 0}
    sam_names = set()
    records = 0
    with sam.open() as inp:
        for line in inp:
            if not line.startswith("@"):
                sam_names.add(line.split("\t", 1)[0])
                records += 1
    if read_names != sam_names:
        raise ValueError(f"STAR SAM read names differ: missing={len(read_names - sam_names)}, "
                         f"extra={len(sam_names - read_names)}")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True,
                        help="Unmodified STAR executable in this project")
    parser.add_argument("--star-oracle", type=Path, required=True)
    parser.add_argument("--nustar-trace", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--nustar-index", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--count", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.count < 1:
        parser.error("--count must be positive")
    for name in ("star", "star_oracle", "nustar_trace", "star_index", "nustar_index", "reads"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sampled = output / "reads.fastq"
    sample_fastq(args.reads, sampled, args.count)

    with tempfile.TemporaryDirectory(prefix="nustar-stitch-example-") as temporary:
        raw = Path(temporary) / "oracle.stderr"
        oracle_wall = run([
            str(args.star_oracle), "--runMode", "alignReads", "--runThreadN", "1",
            "--genomeDir", str(args.star_index), "--readFilesIn", str(sampled),
            "--outFileNamePrefix", str(output / "trace."), "--outSAMtype", "None",
        ], output, output / "trace.stdout", raw,
            dict(os.environ, STAR_STITCH_TRACE="1"))
        nustar = output / "nustar.stitch.tsv"
        nustar_wall = run([
            str(args.nustar_trace), str(args.nustar_index), str(sampled),
            str(raw), str(nustar),
        ], output, output / "nustar.stdout", output / "nustar.stderr")
        windows, candidates, exons, star_seconds, nustar_seconds = compare(raw, nustar)
        split_oracle_trace(raw, output / "star.windows.tsv", output / "star.stitch.tsv")

    baseline_wall = run([
        str(args.star), "--runMode", "alignReads", "--runThreadN", "1",
        "--genomeDir", str(args.star_index), "--readFilesIn", str(sampled),
        "--outFileNamePrefix", str(output / "star."), "--outSAMtype", "SAM",
        "--outSAMunmapped", "Within",
    ], output, output / "star.stdout", output / "star.stderr")
    sam = output / "star.Aligned.out.sam"
    if not sam.is_file():
        raise RuntimeError("STAR did not write star.Aligned.out.sam")
    sam_records = verify_sam_reads(sampled, sam)
    save_summary(sampled, output / "nustar.stitch.tsv", args.star_index,
                 output / "summary.tsv", output / "chromosomes.tsv")
    (output / "FORMAT.md").write_text(
        (PROJECT_ROOT / "docs/stitch_trace_format.md").read_text())

    files = ("reads.fastq", "star.windows.tsv", "star.stitch.tsv",
             "nustar.stitch.tsv", "star.Aligned.out.sam", "summary.tsv",
             "chromosomes.tsv", "FORMAT.md")
    manifest = {
        "sample_reads": args.count,
        "star_source_executable": str(args.star),
        "star_oracle_executable": str(args.star_oracle),
        "nustar_trace_executable": str(args.nustar_trace),
        "source_fastq": str(args.reads),
        "star_index": str(args.star_index),
        "nustar_index": str(args.nustar_index),
        "windows": windows,
        "candidates": candidates,
        "exons": exons,
        "star_sam_records": sam_records,
        "exact_stitch_match_excluding_elapsed_ns": True,
        "star_stitch_seconds": star_seconds,
        "nustar_stitch_seconds": nustar_seconds,
        "oracle_wall_seconds": oracle_wall,
        "nustar_trace_wall_seconds": nustar_wall,
        "star_baseline_wall_seconds": baseline_wall,
        "sha256": {name: sha256(output / name) for name in files},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "README.md").write_text(
        f"# Stitching trace for {args.count} single-end reads\n\n"
        "`reads.fastq` contains the input reads. `star.windows.tsv` contains the ordered "
        "windows and seeds captured from STAR, which are also fed to the NuSTAR trace. "
        "`star.stitch.tsv` and `nustar.stitch.tsv` contain every window candidate, metric, "
        "exon, junction and best candidate in the same order. All fields match except "
        "the elapsed nanoseconds in the final column of `STITCH_OUT`.\n\n"
        "`summary.tsv` lists each read's best score, chromosome, 1-based chromosome position, "
        "strand, and window/candidate counts. `chromosomes.tsv` maps chromosome IDs to "
        "names and global offsets. IDs follow the STAR index order; genome_start is a "
        "0-based offset in the concatenated reference. See `FORMAT.md` for column definitions. "
        "`manifest.json` stores counts, paths and file checksums.\n\n"
        f"`star.Aligned.out.sam` contains {sam_records:,} records from the unmodified STAR "
        f"executable for all {args.count:,} input reads. `--outSAMunmapped Within` keeps "
        "unmapped reads. Other `star.*` and `trace.*` files contain logs and junction output. "
        "This tool compares the stitching stage; it does not run the full NuSTAR alignment "
        "pipeline or establish final SAM equivalence.\n\n"
        "Generated by `NuSTAR/tools/trace/save_stitch_example.py`. Both trace programs "
        "receive the same STAR windows. `run_stitch_diff.py` performs the comparison.\n"
    )
    print(json.dumps({key: manifest[key] for key in
                      ("sample_reads", "windows", "candidates", "exons",
                       "exact_stitch_match_excluding_elapsed_ns")}, indent=2))
    print(output)


if __name__ == "__main__":
    main()
