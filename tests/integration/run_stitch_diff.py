#!/usr/bin/env python3
"""Compare STAR and NuSTAR single-end stitching from STAR's live WA arrays."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time


RECORD_TYPES = ("STITCH_OUT\t", "STITCH_TR\t", "STITCH_META\t", "STITCH_EXON\t",
                "STITCH_JUNCTION\t", "STITCH_BEST\t")


def sample_fastq(source: Path, target: Path, count: int) -> None:
    with source.open() as inp, target.open("w") as out:
        for number in range(count):
            lines = [inp.readline() for _ in range(4)]
            if any(not line for line in lines) or not lines[0].startswith("@"):
                raise RuntimeError(f"FASTQ contains fewer than {count} valid reads: {source}")
            out.writelines(lines)


def run(command: list[str], cwd: Path, stdout: Path, stderr: Path,
        env: dict[str, str] | None = None) -> float:
    started = time.perf_counter()
    with stdout.open("w") as out, stderr.open("w") as err:
        completed = subprocess.run(command, cwd=cwd, stdout=out, stderr=err,
                                   env=env, check=False)
    if completed.returncode:
        tail = "\n".join(stderr.read_text(errors="replace").splitlines()[-15:])
        raise RuntimeError(f"command exited {completed.returncode}: {' '.join(command)}\n{tail}")
    return time.perf_counter() - started


def star_records(path: Path):
    with path.open() as stream:
        for line in stream:
            if line.startswith(RECORD_TYPES):
                yield line.rstrip("\n")


def compare(star_path: Path, nustar_path: Path) -> tuple[int, int, int, float, float]:
    expected = star_records(star_path)
    windows = candidates = exons = 0
    star_seconds = nustar_seconds = 0.0
    with nustar_path.open() as stream:
        for line_number, actual in enumerate(stream, 1):
            actual = actual.rstrip("\n")
            star = next(expected, None)
            if star is None:
                raise RuntimeError(f"NuSTAR has extra stitch row {line_number}: {actual}")
            if actual.startswith("STITCH_OUT\t"):
                if not star.startswith("STITCH_OUT\t"):
                    raise RuntimeError(f"stitch row {line_number} differs:\nSTAR: {star}\nNuSTAR: {actual}")
                star_fields = star.split("\t")
                actual_fields = actual.split("\t")
                star_seconds += int(star_fields[-1]) / 1e9
                nustar_seconds += int(actual_fields[-1]) / 1e9
                star = "\t".join(star_fields[:-1])
                actual = "\t".join(actual_fields[:-1])
                windows += 1
            elif actual.startswith("STITCH_TR\t"):
                candidates += 1
            elif actual.startswith("STITCH_EXON\t"):
                exons += 1
            if star != actual:
                raise RuntimeError(f"stitch row {line_number} differs:\nSTAR: {star}\nNuSTAR: {actual}")
    extra = next(expected, None)
    if extra is not None:
        raise RuntimeError(f"STAR has extra stitch row: {extra}")
    return windows, candidates, exons, star_seconds, nustar_seconds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-oracle", type=Path, required=True)
    parser.add_argument("--nustar-trace", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--nustar-index", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=1000)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.sample < 1:
        parser.error("--sample must be positive")
    for name in ("star_oracle", "nustar_trace", "star_index", "nustar_index", "reads"):
        setattr(args, name, getattr(args, name).resolve(strict=True))

    def execute(root: Path) -> dict:
        reads = root / "reads.fastq"
        sample_fastq(args.reads, reads, args.sample)
        star_stdout, star_stderr = root / "star.stdout", root / "star.stderr"
        nustar_stdout, nustar_stderr = root / "nustar.stdout", root / "nustar.stderr"
        nustar_trace = root / "nustar.tsv"
        star_wall = run([str(args.star_oracle), "--runMode", "alignReads",
                         "--runThreadN", "1", "--genomeDir", str(args.star_index),
                         "--readFilesIn", str(reads),
                         "--outFileNamePrefix", str(root / "star."),
                         "--outSAMtype", "None"], root, star_stdout, star_stderr,
                        dict(os.environ, STAR_STITCH_TRACE="1"))
        nustar_wall = run([str(args.nustar_trace), str(args.nustar_index),
                           str(reads), str(star_stderr), str(nustar_trace)],
                          root, nustar_stdout, nustar_stderr)
        windows, candidates, exons, star_seconds, nustar_seconds = compare(
            star_stderr, nustar_trace)
        return {
            "sample_reads": args.sample,
            "windows": windows,
            "candidates": candidates,
            "exons": exons,
            "exact_stitch_trace_match": True,
            "star_stitch_seconds": star_seconds,
            "nustar_stitch_seconds": nustar_seconds,
            "stitch_speedup": star_seconds / nustar_seconds if nustar_seconds else None,
            "star_process_wall_seconds": star_wall,
            "nustar_process_wall_seconds": nustar_wall,
            "star_index": str(args.star_index),
            "nustar_index": str(args.nustar_index),
            "reads": str(args.reads),
        }

    if args.work_dir:
        root = args.work_dir.resolve()
        root.mkdir(parents=True, exist_ok=False)
        result = execute(root)
    else:
        with tempfile.TemporaryDirectory(prefix="nustar-stitch-diff-") as temporary:
            result = execute(Path(temporary))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
