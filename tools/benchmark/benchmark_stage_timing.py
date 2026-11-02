#!/usr/bin/env python3
"""Measure opt-in STAR and NuSTAR seed timing on identical reads."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "tests/integration"))

from run_stitch_diff import sample_fastq


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_timing(path: Path, expected_reads: int) -> None:
    if not path.is_file():
        raise RuntimeError(f"stage timing file missing: {path}")
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or rows[-1]["status"] != "final" or \
            int(rows[-1]["reads"]) != expected_reads:
        raise RuntimeError(f"incorrect final timing count: {path}")
    counts = [int(row["reads"]) for row in rows[:-1]]
    if any(row["status"] != "progress" for row in rows[:-1]) or \
            any(left >= right for left, right in zip(counts, counts[1:])):
        raise RuntimeError(f"progress counts not increasing: {path}")
    if expected_reads < 10000 and counts:
        raise RuntimeError(f"progress row before 10000 reads: {path}")


def run_star(binary: Path, index: Path, reads: Path, prefix: Path,
             threads: int, every: int | None, expected_reads: int) -> float:
    command = [str(binary), "--runMode", "alignReads", "--runThreadN", str(threads),
               "--genomeDir", str(index), "--readFilesIn", str(reads),
               "--outFileNamePrefix", str(prefix), "--outSAMtype", "None"]
    env = os.environ.copy()
    env.pop("STAR_STAGE_TIMING_EVERY", None)
    env.pop("STAR_STAGE_TIMING_FILE", None)
    if every is not None:
        env["STAR_STAGE_TIMING_EVERY"] = str(every)
    started = time.perf_counter()
    with (Path(str(prefix) + "stdout")).open("w") as out, \
         (Path(str(prefix) + "stderr")).open("w") as err:
        completed = subprocess.run(command, stdout=out, stderr=err, env=env,
                                   check=False)
    elapsed = time.perf_counter()-started
    if completed.returncode:
        raise RuntimeError(f"STAR failed ({completed.returncode}): {prefix}")
    timing = Path(str(prefix) + "Log.stage.progress.tsv")
    if every is not None:
        verify_timing(timing, expected_reads)
    elif timing.exists():
        raise RuntimeError(f"disabled run wrote stage timing file: {timing}")
    return elapsed


def run_nustar_seed(binary: Path, index: Path, reads: Path, prefix: Path,
                    every: int | None, expected_reads: int) -> float:
    trace = Path(str(prefix) + "seed.trace.tsv")
    command = [str(binary), str(index), str(reads), str(trace)]
    env = os.environ.copy()
    env.pop("NUSTAR_STAGE_TIMING_EVERY", None)
    env.pop("NUSTAR_STAGE_TIMING_FILE", None)
    if every is not None:
        env["NUSTAR_STAGE_TIMING_EVERY"] = str(every)
    started = time.perf_counter()
    with Path(str(prefix) + "stdout").open("w") as out, \
         Path(str(prefix) + "stderr").open("w") as err:
        completed = subprocess.run(command, stdout=out, stderr=err, env=env,
                                   check=False)
    elapsed = time.perf_counter()-started
    if completed.returncode:
        raise RuntimeError(f"NuSTAR seed trace failed ({completed.returncode}): {prefix}")
    timing = Path(str(trace) + ".timing.tsv")
    if every is not None:
        verify_timing(timing, expected_reads)
    elif timing.exists():
        raise RuntimeError(f"disabled run wrote stage timing file: {timing}")
    return elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-before", type=Path,
                        help="Optional unmodified executable for the same STAR source")
    parser.add_argument("--star-after", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--nustar-seed-trace", type=Path,
                        help="Optional NuSTAR seed-trace executable")
    parser.add_argument("--nustar-index", type=Path,
                        help="Index required with --nustar-seed-trace")
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=100000)
    parser.add_argument("--every", type=int, default=10000)
    parser.add_argument("--runThreadN", dest="threads", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if min(args.sample, args.threads, args.repeats) < 1 or args.every < 10000:
        parser.error("sample, threads and repeats must be positive; every >= 10000")
    if bool(args.nustar_seed_trace) != bool(args.nustar_index):
        parser.error("--nustar-seed-trace and --nustar-index must be supplied together")
    for name in ("star_after", "star_index", "reads"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    if args.star_before:
        args.star_before = args.star_before.resolve(strict=True)
    if args.nustar_seed_trace:
        args.nustar_seed_trace = args.nustar_seed_trace.resolve(strict=True)
        args.nustar_index = args.nustar_index.resolve(strict=True)
    variants = {"after_disabled": (args.star_after, None),
                "after_enabled": (args.star_after, args.every)}
    if args.star_before:
        variants["before"] = (args.star_before, None)
    times = {name: [] for name in variants}
    if args.nustar_seed_trace:
        times["nustar_seed_disabled"] = []
        times["nustar_seed_enabled"] = []
    with tempfile.TemporaryDirectory(prefix="star-stage-benchmark-") as temporary:
        root = Path(temporary)
        sampled = root / "reads.fastq"
        sample_fastq(args.reads, sampled, args.sample)
        names = list(variants)
        for repeat in range(args.repeats):
            for name in (names if repeat % 2 == 0 else list(reversed(names))):
                binary, every = variants[name]
                prefix = root / f"{repeat}.{name}."
                times[name].append(run_star(binary, args.star_index, sampled,
                                            prefix, args.threads, every,
                                            args.sample))
            if args.nustar_seed_trace:
                seed_names = ["nustar_seed_disabled", "nustar_seed_enabled"]
                if repeat % 2: seed_names.reverse()
                for name in seed_names:
                    prefix = root / f"{repeat}.{name}."
                    interval = args.every if name.endswith("enabled") else None
                    times[name].append(run_nustar_seed(
                        args.nustar_seed_trace, args.nustar_index, sampled,
                        prefix, interval, args.sample))
    medians = {name: statistics.median(values) for name, values in times.items()}
    result = {
        "sample_reads": args.sample,
        "threads": args.threads,
        "interval_reads": args.every,
        "repeats": args.repeats,
        "star_before": str(args.star_before) if args.star_before else None,
        "star_before_sha256": sha256(args.star_before) if args.star_before else None,
        "star_after": str(args.star_after),
        "star_after_sha256": sha256(args.star_after),
        "star_index": str(args.star_index),
        "nustar_seed_trace": str(args.nustar_seed_trace) if args.nustar_seed_trace else None,
        "nustar_seed_trace_sha256": (
            sha256(args.nustar_seed_trace) if args.nustar_seed_trace else None),
        "nustar_index": str(args.nustar_index) if args.nustar_index else None,
        "source_reads": str(args.reads),
        "wall_seconds": times,
        "median_wall_seconds": medians,
        "enabled_over_disabled_ratio": medians["after_enabled"] / medians["after_disabled"],
    }
    if args.star_before:
        result["disabled_over_before_ratio"] = (
            medians["after_disabled"] / medians["before"])
    if args.nustar_seed_trace:
        result["nustar_seed_enabled_over_disabled_ratio"] = (
            medians["nustar_seed_enabled"] / medians["nustar_seed_disabled"])
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
