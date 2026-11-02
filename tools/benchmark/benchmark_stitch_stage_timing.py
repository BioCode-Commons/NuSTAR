#!/usr/bin/env python3
"""Measure the optional progress timer in NuSTAR's STAR-window stitch replay."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

from benchmark_stage_timing import sha256, verify_timing


def merge_trace(windows: Path, candidates: Path, target: Path) -> tuple[int, int]:
    best_rows = {}
    out_rows = {}
    with candidates.open() as stream:
        for line in stream:
            if line.startswith("STITCH_BEST\t"):
                name = line.split("\t", 2)[1]
                if name in best_rows:
                    raise RuntimeError(f"duplicate STITCH_BEST: {name}")
                best_rows[name] = line
            elif line.startswith("STITCH_OUT\t"):
                parts = line.split("\t", 3)
                key = (parts[1], parts[2])
                if key in out_rows:
                    raise RuntimeError(f"duplicate STITCH_OUT: {key}")
                out_rows[key] = line
    count = 0
    window_count = len(out_rows)
    previous = None
    pending_window = None
    with windows.open() as source, target.open("w") as output:
        for line in source:
            if line.startswith("STITCH_READ\t"):
                if pending_window is not None:
                    output.write(out_rows.pop(pending_window))
                    pending_window = None
                if previous is not None:
                    output.write(best_rows.pop(previous))
                previous = line.split("\t", 2)[1]
                count += 1
            elif line.startswith("STITCH_WIN\t"):
                if pending_window is not None:
                    output.write(out_rows.pop(pending_window))
                parts = line.split("\t", 4)
                pending_window = (parts[1], parts[2])
            output.write(line)
        if pending_window is not None:
            output.write(out_rows.pop(pending_window))
        if previous is not None:
            output.write(best_rows.pop(previous))
    if count == 0 or window_count == 0 or best_rows or out_rows:
        raise RuntimeError("STAR window and candidate rows differ")
    return count, window_count


def run(binary: Path, index: Path, reads: Path, windows: Path,
        output: Path, every: int | None, expected_reads: int,
        expected_windows: int) -> tuple[float, str]:
    env = os.environ.copy()
    env.pop("NUSTAR_STAGE_TIMING_EVERY", None)
    env.pop("NUSTAR_STAGE_TIMING_FILE", None)
    if every is not None:
        env["NUSTAR_STAGE_TIMING_EVERY"] = str(every)
    command = [str(binary), str(index), str(reads), str(windows), str(output)]
    started = time.perf_counter()
    with Path(str(output) + ".stdout").open("w") as out, \
         Path(str(output) + ".stderr").open("w") as err:
        completed = subprocess.run(command, stdout=out, stderr=err, env=env,
                                   check=False)
    elapsed = time.perf_counter()-started
    if completed.returncode:
        raise RuntimeError(f"stitch trace failed ({completed.returncode}): {output}")
    timing = Path(str(output) + ".timing.tsv")
    if every is not None:
        verify_timing(timing, expected_reads)
    elif timing.exists():
        raise RuntimeError(f"disabled run wrote timing file: {timing}")
    digest = hashlib.sha256()
    observed_windows = 0
    with output.open("rb") as stream:
        for line in stream:
            if line.startswith(b"STITCH_OUT\t"):
                observed_windows += 1
                line = line.rsplit(b"\t", 1)[0] + b"\n"
            digest.update(line)
    if observed_windows != expected_windows:
        raise RuntimeError(f"stitch window count differs: {observed_windows}")
    return elapsed, digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nustar-stitch-trace", type=Path, required=True)
    parser.add_argument("--nustar-index", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--star-window-trace", type=Path, required=True)
    parser.add_argument("--star-candidates", type=Path, required=True,
                        help="STAR candidate trace supplying STITCH_BEST boundaries")
    parser.add_argument("--every", type=int, default=10000)
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.every < 10000 or args.repeats < 1:
        parser.error("every must be at least 10000 and repeats must be positive")
    for name in ("nustar_stitch_trace", "nustar_index", "reads",
                 "star_window_trace", "star_candidates"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    times = {"disabled": [], "enabled": []}
    expected_digest = None
    with tempfile.TemporaryDirectory(prefix="stitch-stage-benchmark-") as temporary:
        root = Path(temporary)
        merged = root / "star.replay.tsv"
        read_count, window_count = merge_trace(args.star_window_trace,
                                               args.star_candidates, merged)
        for repeat in range(args.repeats):
            names = ("disabled", "enabled") if repeat % 2 == 0 else ("enabled", "disabled")
            for name in names:
                output = root / f"{repeat}.{name}.tsv"
                interval = args.every if name == "enabled" else None
                elapsed, digest = run(args.nustar_stitch_trace,
                    args.nustar_index, args.reads, merged,
                    output, interval, read_count, window_count)
                if expected_digest is not None and digest != expected_digest:
                    raise RuntimeError("stitch traces differ after ignoring elapsed_ns")
                expected_digest = digest
                times[name].append(elapsed)
    medians = {name: statistics.median(values) for name, values in times.items()}
    result = {
        "nustar_stitch_trace": str(args.nustar_stitch_trace),
        "nustar_stitch_trace_sha256": sha256(args.nustar_stitch_trace),
        "nustar_index": str(args.nustar_index),
        "reads_sha256": sha256(args.reads),
        "star_window_trace_sha256": sha256(args.star_window_trace),
        "star_candidates_sha256": sha256(args.star_candidates),
        "reads": read_count,
        "windows": window_count,
        "interval_reads": args.every,
        "repeats": args.repeats,
        "wall_seconds": times,
        "median_wall_seconds": medians,
        "enabled_over_disabled_ratio": medians["enabled"] / medians["disabled"],
    }
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
