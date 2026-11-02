#!/usr/bin/env python3
"""Compare STAR's seed-stage PC trace against NuSTAR on FASTQ reads.

The STAR oracle must be built by build_star_seed_oracle.py from a separate
source copy. Trace generation and comparison live outside production code.
"""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time

TRACE_PREFIXES = ("SEED_SUM\t", "SEED_SPLIT\t", "SEED_PC\t")


def sample_fastq(source: Path, target: Path, count: int) -> int:
    copied = 0
    with source.open() as inp, target.open("w") as out:
        while copied < count:
            lines = [inp.readline() for _ in range(4)]
            if not lines[0]:
                break
            if any(not line for line in lines) or not lines[0].startswith("@"):
                raise RuntimeError(f"invalid FASTQ record in {source}")
            out.writelines(lines)
            copied += 1
    return copied


def run(command: list[str], cwd: Path, stdout: Path, stderr: Path,
        env: dict[str, str] | None = None) -> float:
    started = time.perf_counter()
    with stdout.open("w") as out, stderr.open("w") as err:
        completed = subprocess.run(command, cwd=cwd, stdout=out, stderr=err,
                                   env=env, check=False)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): "
                           f"{' '.join(command)}\n" +
                           "\n".join(stderr.read_text(errors="replace").splitlines()[-15:]))
    return time.perf_counter()-started


def compare(star_trace: Path, nustar_trace: Path) -> tuple[int, int, float]:
    def records():
        with star_trace.open() as source:
            for line in source:
                if line.startswith(TRACE_PREFIXES):
                    yield line
    star_records = records()
    reads = pieces = lines = 0
    with nustar_trace.open() as source:
        for actual in source:
            expected = next(star_records, None)
            lines += 1
            if expected != actual:
                raise RuntimeError(f"seed trace differs at line {lines}:\n"
                                   f"STAR: {expected!r}\nNuSTAR: {actual!r}")
            reads += actual.startswith("SEED_SUM\t")
            pieces += actual.startswith("SEED_PC\t")
    extra = next(star_records, None)
    if extra is not None:
        raise RuntimeError(f"STAR has an extra seed record: {extra!r}")
    star_ns = 0
    timed_reads = 0
    with star_trace.open() as source:
        for line in source:
            if line.startswith("SEED_TIME\t"):
                star_ns += int(line.rstrip().split("\t")[2])
                timed_reads += 1
    if timed_reads != reads:
        raise RuntimeError(f"STAR timed {timed_reads} reads, compared {reads}")
    return reads, pieces, star_ns/1e9


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-oracle", type=Path, required=True)
    parser.add_argument("--nustar-trace", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--nustar-index", type=Path, required=True)
    parser.add_argument("--mate1", type=Path, required=True)
    parser.add_argument("--mate2", type=Path)
    parser.add_argument("--sample", type=int, default=10000)
    parser.add_argument("--seed-search-lmax", type=int)
    parser.add_argument("--seed-start-lmax", type=int)
    parser.add_argument("--seed-multimap-nmax", type=int)
    parser.add_argument("--seed-split-min", type=int)
    parser.add_argument("--seed-map-min", type=int)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.sample < 1:
        parser.error("--sample must be positive")
    for field in ("star_oracle", "nustar_trace", "star_index", "nustar_index",
                  "mate1", "mate2"):
        value = getattr(args, field)
        if value is not None:
            setattr(args, field, value.resolve(strict=True))

    def execute(root: Path) -> dict:
        first = root/"mate1.fastq"
        second = root/"mate2.fastq"
        count = sample_fastq(args.mate1, first, args.sample)
        if count != args.sample:
            raise RuntimeError(f"only {count} reads available")
        if args.mate2 and sample_fastq(args.mate2, second, args.sample) != count:
            raise RuntimeError("paired FASTQ record counts differ")
        star_stderr = root/"star.stderr"
        star_stdout = root/"star.stdout"
        nustar_stderr = root/"nustar.stderr"
        nustar_stdout = root/"nustar.stdout"
        nustar_trace = root/"nustar.tsv"
        reads = [str(first)] + ([str(second)] if args.mate2 else [])
        seed_options = []
        for value, star_name in (
            (args.seed_search_lmax, "seedSearchLmax"),
            (args.seed_start_lmax, "seedSearchStartLmax"),
            (args.seed_multimap_nmax, "seedMultimapNmax"),
            (args.seed_split_min, "seedSplitMin"),
            (args.seed_map_min, "seedMapMin"),
        ):
            if value is not None:
                seed_options.extend(["--"+star_name, str(value)])
        env = dict(os.environ, STAR_SEED_TRACE="1")
        star_wall = run([str(args.star_oracle), "--runMode", "alignReads",
                         "--runThreadN", "1", "--genomeDir", str(args.star_index),
                         "--readFilesIn", *reads,
                         "--outFileNamePrefix", str(root/"star."),
                         "--outSAMtype", "None", *seed_options], root,
                        star_stdout, star_stderr, env)
        nustar_wall = run([str(args.nustar_trace), str(args.nustar_index),
                           str(first), *([str(second)] if args.mate2 else []),
                           str(nustar_trace), *seed_options], root,
                          nustar_stdout, nustar_stderr)
        compared_reads, pieces, star_seed_seconds = compare(star_stderr, nustar_trace)
        timing = re.search(r"search_seconds=([0-9.]+)", nustar_stderr.read_text())
        if not timing:
            raise RuntimeError("NuSTAR seed timing is absent")
        nustar_seed_seconds = float(timing.group(1))
        return {
            "sample_reads": compared_reads,
            "seed_pc_records": pieces,
            "exact_trace_match": True,
            "star_seed_seconds": star_seed_seconds,
            "nustar_seed_seconds": nustar_seed_seconds,
            "seed_speedup": star_seed_seconds/nustar_seed_seconds,
            "star_process_wall_seconds": star_wall,
            "nustar_process_wall_seconds": nustar_wall,
            "star_index": str(args.star_index),
            "nustar_index": str(args.nustar_index),
            "mate1": str(args.mate1),
            "mate2": str(args.mate2) if args.mate2 else None,
            "seed_options": seed_options,
        }

    if args.work_dir:
        root = args.work_dir.resolve()
        root.mkdir(parents=True, exist_ok=False)
        result = execute(root)
    else:
        with tempfile.TemporaryDirectory(prefix="nustar-seed-diff-") as temporary:
            result = execute(Path(temporary))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
