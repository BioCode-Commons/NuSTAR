#!/usr/bin/env python3
"""Compare two seed-benchmark binaries with alternating process order.

Example:
  python3 tools/benchmark/benchmark_seed_search.py \
      --baseline build/before/seed-benchmark \
      --candidate build/after/seed-benchmark --index index \
      --mate1 reads_1.fastq --mate2 reads_2.fastq \
      --reads 100000 --repeats 5 --sweeps 5 --output seed-comparison.json

Both binaries must be built from the same seed_benchmark.cpp, using the same
compiler and flags. Each process performs its own untimed warmup, --sweeps bulk
timed search sweeps, and one untimed full-result SHA-256 validation sweep. Every
loaded fragment participates in validation. Peak RSS includes the index and
preloaded reads, and is reported in bytes by seed-benchmark. Run without other
CPU-intensive jobs; this is a single-threaded seed microbenchmark, not a claim
about end-to-end throughput or other hardware/data.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
from typing import Any


SEED_OPTIONS = (
    "seedSearchLmax", "seedSearchStartLmax", "seedMultimapNmax",
    "seedSplitMin", "seedMapMin",
)


def positive(value: str) -> int:
    parsed = nonnegative(value)
    if parsed == 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def nonnegative(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if parsed < 0 or parsed > 2**64-1:
        raise argparse.ArgumentTypeError("must fit an unsigned 64-bit integer")
    return parsed


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def input_prefix(path: Path, reads: int) -> dict[str, Any]:
    """Hash exactly the FASTQ records that participated, without hashing the tail."""
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        for ordinal in range(reads):
            for line in range(4):
                value = source.readline()
                if not value:
                    raise RuntimeError(f"{path}: truncated selected FASTQ record {ordinal+1}")
                if line == 0 and not value.startswith(b"@"):
                    raise RuntimeError(f"{path}: invalid selected FASTQ record {ordinal+1}")
                digest.update(value)
                size += len(value)
    return {**metadata(path), "selected_records": reads,
            "selected_bytes": size, "selected_prefix_sha256": digest.hexdigest()}


def summarize(values: list[float | int]) -> dict[str, float | int]:
    return {"median": statistics.median(values), "min": min(values), "max": max(values)}


def run(command: list[str]) -> dict[str, Any]:
    started = time.perf_counter()
    process = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, check=False)
    wall_seconds = time.perf_counter()-started
    if process.returncode:
        raise RuntimeError(
            f"benchmark failed ({process.returncode}): {command!r}\n{process.stderr[-8000:]}")
    try:
        result = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"benchmark did not return JSON: {command!r}\n"
                           f"{process.stdout[-4000:]}\n{process.stderr[-4000:]}") from error
    if result.get("schema") != "nustar-seed-benchmark-v1":
        raise RuntimeError(f"unsupported benchmark schema: {result.get('schema')!r}")
    rounds = result.get("rounds", [])
    if len(rounds) != 1:
        raise RuntimeError("expected exactly one timed round from each process")
    measured = rounds[0]
    seconds = measured.get("seconds")
    if not isinstance(seconds, (int, float)) or not math.isfinite(seconds) or seconds <= 0:
        raise RuntimeError(f"invalid search time: {seconds!r}")
    cpu_seconds = measured.get("cpu_seconds")
    if not isinstance(cpu_seconds, (int, float)) or not math.isfinite(cpu_seconds) or cpu_seconds <= 0:
        raise RuntimeError(f"invalid process CPU time: {cpu_seconds!r}")
    total_seconds = measured.get("total_seconds")
    sweeps = measured.get("timed_sweeps")
    if (not isinstance(sweeps, int) or sweeps < 1 or result.get("timed_sweeps") != sweeps
            or not isinstance(total_seconds, (int, float)) or not math.isfinite(total_seconds)
            or total_seconds <= 0
            or not math.isclose(seconds*sweeps, total_seconds, rel_tol=1e-9)):
        raise RuntimeError("invalid or inconsistent timed sweep count/total time")
    digest = result.get("digest", "")
    if (len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest)
            or measured.get("digest") != digest):
        raise RuntimeError("missing or inconsistent SHA-256 output digest")
    if measured.get("reads") != result.get("reads") or result["reads"] < 1:
        raise RuntimeError("invalid or inconsistent fragment count")
    rss = result.get("peak_rss_bytes")
    if not isinstance(rss, int) or rss <= 0:
        raise RuntimeError("peak RSS is unavailable or invalid on this platform")
    return {"command": command, "process_wall_seconds": wall_seconds,
            "seconds": seconds, "total_seconds": total_seconds, "timed_sweeps": sweeps,
            "cpu_seconds": cpu_seconds,
            "peak_rss_bytes": rss, "result": result,
            "stderr": process.stderr}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--mate1", type=Path, required=True)
    parser.add_argument("--mate2", type=Path)
    parser.add_argument("--reads", type=positive, default=100000,
                        help="maximum fragments preloaded in each process (default: 100000)")
    parser.add_argument("--repeats", type=positive, default=5,
                        help="process runs per binary, alternating A/B and B/A (default: 5)")
    parser.add_argument("--sweeps", type=positive, default=5,
                        help="timed sweeps per process; seconds is per sweep (default: 5)")
    parser.add_argument("--output", type=Path, required=True)
    for option in SEED_OPTIONS:
        parser.add_argument("--"+option, type=nonnegative)
    args = parser.parse_args()
    binaries = {name: getattr(args, name).resolve(strict=True)
                for name in ("baseline", "candidate")}
    for name, path in binaries.items():
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f"{name} must be an executable file: {path}")
    index = args.index.resolve(strict=True)
    mate1 = args.mate1.resolve(strict=True)
    mate2 = args.mate2.resolve(strict=True) if args.mate2 else None
    index_files = [index/"index.sa", index/"reference.bin"]
    input_files = [mate1] + ([mate2] if mate2 else [])
    file_state = {path: metadata(path) for path in [*index_files, *input_files]}
    binary_hashes = {name: sha256(path) for name, path in binaries.items()}
    policies: list[str] = []
    for option in SEED_OPTIONS:
        value = getattr(args, option)
        if value is not None:
            policies += ["--"+option, str(value)]
    commands = {
        name: [str(path), str(index), str(mate1), str(mate2) if mate2 else "--",
               "1", "--reads", str(args.reads), "--sweeps", str(args.sweeps), *policies]
        for name, path in binaries.items()
    }
    observed: list[dict[str, Any]] = []
    expected: dict[str, Any] | None = None
    checked_fields = ("reads", "input_bases", "paired", "pieces_per_sweep",
                      "digest_algorithm", "digest", "policy")
    for repetition in range(args.repeats):
        order = ("baseline", "candidate") if repetition % 2 == 0 else ("candidate", "baseline")
        for name in order:
            measurement = run(commands[name])
            current = {key: measurement["result"][key] for key in checked_fields}
            if expected is None:
                expected = current
            if current != expected:
                difference = {key: {"expected": expected[key], "actual": current[key]}
                              for key in checked_fields if expected[key] != current[key]}
                raise RuntimeError(f"baseline/candidate output differs in repetition "
                                   f"{repetition+1} ({name}): {json.dumps(difference)}")
            if current["reads"] > args.reads or measurement["timed_sweeps"] != args.sweeps:
                raise RuntimeError("benchmark exceeded requested input limit or used wrong sweep count")
            observed.append({"variant": name, "repeat": repetition+1, **measurement})
            print(f"repeat {repetition+1}/{args.repeats} {name}: "
                  f"{measurement['seconds']:.6f} seed seconds, "
                  f"{measurement['peak_rss_bytes']} peak RSS bytes, digest matched",
                  file=sys.stderr, flush=True)

    # Reject concurrent recompilation or a changed input/index. Full selected
    # input hashes are recorded below; stat checks avoid scanning an entire
    # multi-gigabyte suffix array during a lightweight seed benchmark.
    for name, path in binaries.items():
        if sha256(path) != binary_hashes[name]:
            raise RuntimeError(f"{name} binary changed during benchmarking")
    for path, before in file_state.items():
        if metadata(path) != before:
            raise RuntimeError(f"input/index changed during benchmarking: {path}")
    assert expected is not None
    summaries = {}
    for name in binaries:
        rows = [row for row in observed if row["variant"] == name]
        summaries[name] = {
            "seed_seconds": summarize([row["seconds"] for row in rows]),
            "seed_cpu_seconds": summarize([row["cpu_seconds"] for row in rows]),
            "timed_total_seconds": summarize([row["total_seconds"] for row in rows]),
            "peak_rss_bytes": summarize([row["peak_rss_bytes"] for row in rows]),
            "process_wall_seconds": summarize([row["process_wall_seconds"] for row in rows]),
        }
    baseline_time = summaries["baseline"]["seed_seconds"]["median"]
    candidate_time = summaries["candidate"]["seed_seconds"]["median"]
    baseline_rss = summaries["baseline"]["peak_rss_bytes"]["median"]
    candidate_rss = summaries["candidate"]["peak_rss_bytes"]["median"]
    source = Path(__file__).resolve()
    report = {
        "schema": "nustar-seed-comparison-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "platform": {"system": platform.system(), "release": platform.release(),
                     "machine": platform.machine(), "processor": platform.processor(),
                     "python": sys.version, "cpu_count": os.cpu_count()},
        "provenance": {
            "runner": {"path": str(source), "sha256": sha256(source)},
            "binaries": {name: {**metadata(path), "sha256": binary_hashes[name]}
                         for name, path in binaries.items()},
            "commands": commands,
            "index": {"directory": str(index),
                      "files": [file_state[path] for path in index_files],
                      "identity_note": "index files identified by size/mtime; no full index hash"},
            "inputs": [input_prefix(path, expected["reads"]) for path in input_files],
        },
        "read_limit": args.reads, "repeats_per_binary": args.repeats,
        "timed_sweeps_per_process": args.sweeps,
        "order": [row["variant"] for row in observed],
        "correctness": {"all_runs_same_full_result_digest": True,
                        "validation": "every loaded fragment in separate untimed sweeps",
                        **expected},
        "summary": summaries,
        "seed_speedup_baseline_over_candidate": baseline_time/candidate_time,
        "seed_cpu_speedup_baseline_over_candidate": (
            summaries["baseline"]["seed_cpu_seconds"]["median"] /
            summaries["candidate"]["seed_cpu_seconds"]["median"]),
        "median_peak_rss_delta_bytes": candidate_rss-baseline_rss,
        "median_peak_rss_ratio_candidate_over_baseline": candidate_rss/baseline_rss,
        "timing_scope": "mean seconds per single-threaded search sweep with reused output; excludes preload and digest",
        "memory_scope": "whole-process peak RSS including preloaded reads, index and validation",
        "runs": observed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"output": str(args.output.resolve()),
                      "seed_speedup": report["seed_speedup_baseline_over_candidate"],
                      "median_peak_rss_delta_bytes": report["median_peak_rss_delta_bytes"],
                      "reads": expected["reads"], "digest": expected["digest"]}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, KeyError) as error:
        print(f"benchmark_seed_search: {error}", file=sys.stderr)
        sys.exit(1)
