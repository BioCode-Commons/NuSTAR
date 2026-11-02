#!/usr/bin/env python3
"""Compare seed-stage time using legacy raw and packed NuSTAR SA files.

Both indexes must contain the same reference and suffix ranks. Run the STAR
stage differential test separately to prove the PC trace is unchanged.
"""

import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "tests/integration"))

from run_seed_diff import sample_fastq


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-trace", type=Path, required=True)
    parser.add_argument("--legacy-index", type=Path, required=True)
    parser.add_argument("--packed-index", type=Path, required=True)
    parser.add_argument("--mate1", type=Path, required=True)
    parser.add_argument("--mate2", type=Path)
    parser.add_argument("--sample", type=int, default=100000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.sample < 1 or args.repeats < 1:
        parser.error("--sample and --repeats must be positive")
    for name in ("seed_trace", "legacy_index", "packed_index", "mate1", "mate2"):
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.resolve(strict=True))

    times = {"legacy": [], "packed": []}
    counts = None
    with tempfile.TemporaryDirectory(prefix="nustar-packed-sa-bench-") as temporary:
        root = Path(temporary)
        first = root/"mate1.fastq"
        second = root/"mate2.fastq"
        if sample_fastq(args.mate1, first, args.sample) != args.sample:
            raise RuntimeError("mate1 has fewer reads than requested")
        if args.mate2 and sample_fastq(args.mate2, second, args.sample) != args.sample:
            raise RuntimeError("mate2 has fewer reads than requested")
        reads = [str(first)] + ([str(second)] if args.mate2 else [])
        for repeat in range(args.repeats):
            order = ("legacy", "packed") if repeat % 2 == 0 else ("packed", "legacy")
            for name in order:
                index = args.legacy_index if name == "legacy" else args.packed_index
                run = subprocess.run([
                    str(args.seed_trace), str(index), *reads, os.devnull],
                    text=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                    check=False)
                if run.returncode:
                    raise RuntimeError(f"{name} seed trace failed: {run.stderr}")
                match = re.search(
                    r"reads=(\d+) pieces=(\d+) search_seconds=([0-9.]+)",
                    run.stderr)
                if not match:
                    raise RuntimeError(f"{name} seed timing missing: {run.stderr}")
                current_counts = (int(match[1]), int(match[2]))
                if counts is not None and current_counts != counts:
                    raise RuntimeError("raw and packed SA produced different counts")
                counts = current_counts
                times[name].append(float(match[3]))

    legacy_bytes = (args.legacy_index/"index.sa").stat().st_size
    packed_bytes = (args.packed_index/"index.sa").stat().st_size
    result = {
        "sample_reads": counts[0],
        "seed_pc_records": counts[1],
        "repeats": args.repeats,
        "legacy_index_bytes": legacy_bytes,
        "packed_index_bytes": packed_bytes,
        "file_reduction_percent": 100*(legacy_bytes-packed_bytes)/legacy_bytes,
        "legacy_seed_seconds": times["legacy"],
        "packed_seed_seconds": times["packed"],
        "legacy_median_seconds": statistics.median(times["legacy"]),
        "packed_median_seconds": statistics.median(times["packed"]),
        "packed_over_legacy_time_ratio":
            statistics.median(times["packed"])/statistics.median(times["legacy"]),
    }
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
