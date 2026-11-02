#!/usr/bin/env python3
"""Alternating, exact-SAM benchmark of project STAR and native mapping.

STAR runs one thread for stable record order; NuSTAR can use multiple workers.
STAR output headers contain its command line;
@SQ rows and every ordered SAM body row are compared exactly.
"""

import argparse
import itertools
import json
import statistics
import subprocess
import time
from pathlib import Path

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "tests/integration"))

from compare_sam import compare_ordered_sam


def numeric_statistics(path: Path) -> dict[str, str]:
    changing = {"Started job on", "Started mapping on", "Finished on",
                "Mapping speed, Million of reads per hour"}
    rows = {}
    for line in path.read_text().splitlines():
        if "|" not in line:
            continue
        label, value = line.split("|", 1)
        label = label.strip()
        if label not in changing:
            rows[label] = value.strip()
    return rows


def run(command: list[str], stdout: Path, stderr: Path) -> float:
    start = time.perf_counter()
    with stdout.open("w") as out, stderr.open("w") as err:
        subprocess.run(command, stdout=out, stderr=err, check=True)
    return time.perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--nustar", type=Path, required=True)
    parser.add_argument("--nustar-mode", choices=("trace", "product", "paired"),
                        default="trace")
    parser.add_argument("--nustar-threads", type=int, default=1)
    parser.add_argument("--io-backend", choices=("mmap", "buffered"), default="mmap")
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--nustar-index", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--mate2", type=Path)
    parser.add_argument("--compare-sj", action="store_true")
    parser.add_argument("--compare-stats", action="store_true")
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if not 1 <= args.nustar_threads <= 256:
        parser.error("--nustar-threads must be between 1 and 256")
    if args.nustar_mode == "trace" and args.nustar_threads != 1:
        parser.error("multi-threaded benchmarking requires product or paired mode")
    if (args.nustar_mode == "paired") != (args.mate2 is not None):
        parser.error("--mate2 is required exactly when --nustar-mode paired")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for repetition in range(args.repeats):
        star_prefix = args.work_dir / f"star-{repetition}."
        nustar_prefix = str(args.work_dir / f"nustar-{repetition}.")
        nustar_sam = Path(nustar_prefix + "Aligned.out.sam")
        commands = {
            "star": [str(args.star.resolve()), "--runMode", "alignReads",
                     "--genomeDir", str(args.star_index.resolve()),
                     "--readFilesIn", str(args.reads.resolve())] +
                    ([str(args.mate2.resolve())] if args.mate2 else []) + [
                     "--runThreadN", "1", "--outSAMtype", "SAM",
                     "--outSAMunmapped", "Within",
                     "--outFileNamePrefix", str(star_prefix), "--outStd", "Log"],
            "nustar": ([str(args.nustar.resolve()), str(args.nustar_index.resolve()),
                        str(args.reads.resolve()), str(nustar_sam)] if args.nustar_mode == "trace" else
                       [str(args.nustar.resolve()), "--runMode", "alignReads",
                        "--genomeDir", str(args.nustar_index.resolve()),
                        "--readFilesIn", str(args.reads.resolve())] +
                       ([str(args.mate2.resolve())] if args.mate2 else []) +
                       ["--outFileNamePrefix", nustar_prefix, "--outSAMtype", "SAM",
                        "--outSAMunmapped", "Within", "--outSAMorder", "PairedKeepInputOrder",
                        "--runThreadN", str(args.nustar_threads), "--io-backend", args.io_backend]),
        }
        times = {}
        for name in ("star", "nustar") if repetition % 2 == 0 else ("nustar", "star"):
            times[name] = run(commands[name],
                args.work_dir / f"{name}-{repetition}.stdout",
                args.work_dir / f"{name}-{repetition}.stderr")
        star_sam = Path(str(star_prefix) + "Aligned.out.sam")
        sam_records = compare_ordered_sam(star_sam, nustar_sam)
        sj_rows = None
        if args.compare_sj:
            star_sj = Path(str(star_prefix) + "SJ.out.tab")
            nustar_sj = Path((str(nustar_sam) + "." if args.nustar_mode == "trace" else nustar_prefix) + "SJ.out.tab")
            wanted = star_sj.read_bytes()
            got = nustar_sj.read_bytes()
            if wanted != got:
                star_lines = wanted.decode().splitlines()
                nustar_lines = got.decode().splitlines()
                first = next((i for i, (a, b) in enumerate(
                    itertools.zip_longest(star_lines, nustar_lines), 1)
                    if a != b), None)
                raise AssertionError(f"SJ.out.tab mismatch in repetition {repetition}: "
                                     f"STAR={len(star_lines)}, NuSTAR={len(nustar_lines)}, "
                                     f"first_diff={first}")
            sj_rows = len(wanted.splitlines())
        stats_fields = None
        if args.compare_stats:
            star_stats = numeric_statistics(Path(str(star_prefix) + "Log.final.out"))
            nustar_stats = numeric_statistics(
                Path((str(nustar_sam) + "." if args.nustar_mode == "trace" else nustar_prefix) + "Log.final.out"))
            if star_stats != nustar_stats:
                first = next(((key, star_stats.get(key), nustar_stats.get(key))
                              for key in sorted(star_stats.keys() | nustar_stats.keys())
                              if star_stats.get(key) != nustar_stats.get(key)), None)
                raise AssertionError(f"final statistics mismatch in repetition "
                                     f"{repetition}: {first}")
            stats_fields = len(star_stats)
        results.append({"repeat": repetition, "star_seconds": times["star"],
                        "nustar_seconds": times["nustar"],
                        "exact_sam_records": sam_records,
                        "exact_sj_rows": sj_rows,
                        "exact_statistics_fields": stats_fields})
        print(results[-1], flush=True)
    star_median = statistics.median(row["star_seconds"] for row in results)
    nustar_median = statistics.median(row["nustar_seconds"] for row in results)
    report = {
        "star": str(args.star.resolve()), "nustar": str(args.nustar.resolve()),
        "nustar_mode": args.nustar_mode,
        "reads": str(args.reads.resolve()), "star_threads": 1,
        "nustar_threads": args.nustar_threads, "io_backend": args.io_backend,
        "mate2": str(args.mate2.resolve()) if args.mate2 else None,
        "sam_header_comparison": "excluded (STAR records its own command line)",
        "sam_body_comparison": "ordered, all fields, tags and multiplicity exact",
        "sj_table_comparison": "byte exact" if args.compare_sj else "not compared",
        "statistics_comparison": "numeric fields exact; timestamps and mapping speed excluded"
            if args.compare_stats else "not compared",
        "runs": results, "star_median_seconds": star_median,
        "nustar_median_seconds": nustar_median,
        "star_over_nustar_speed_ratio": star_median / nustar_median,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"median STAR={star_median:.3f}s NuSTAR={nustar_median:.3f}s "
          f"ratio={star_median / nustar_median:.3f}x")


if __name__ == "__main__":
    main()
