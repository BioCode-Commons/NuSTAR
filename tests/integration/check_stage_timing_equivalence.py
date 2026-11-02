#!/usr/bin/env python3
"""Check that stage timing leaves STAR SAM and NuSTAR seed traces unchanged."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "tools/benchmark"))

from benchmark_stage_timing import verify_timing


def environment(variable: str, every: int | None) -> dict[str, str]:
    env = os.environ.copy()
    env.pop(variable, None)
    env.pop(variable.replace("EVERY", "FILE"), None)
    if every is not None:
        env[variable] = str(every)
    return env


def execute(command: list[str], env: dict[str, str], prefix: Path) -> None:
    with Path(str(prefix) + ".stdout").open("w") as out, \
         Path(str(prefix) + ".stderr").open("w") as err:
        completed = subprocess.run(command, env=env, stdout=out, stderr=err,
                                   check=False)
    if completed.returncode:
        detail = Path(str(prefix) + ".stderr").read_text(errors="replace")[-2000:]
        raise RuntimeError(f"command failed ({completed.returncode}): {detail}")


def sam_records(path: Path) -> list[bytes]:
    return [line for line in path.read_bytes().splitlines(keepends=True)
            if not line.startswith(b"@")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--star-index", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--nustar-seed-trace", type=Path)
    parser.add_argument("--nustar-index", type=Path)
    parser.add_argument("--every", type=int, default=10000)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.every < 10000 or bool(args.nustar_seed_trace) != bool(args.nustar_index):
        parser.error("--every >= 10000 and paired NuSTAR executable/index are required")
    for name in ("star", "star_index", "reads"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    if args.nustar_seed_trace:
        args.nustar_seed_trace = args.nustar_seed_trace.resolve(strict=True)
        args.nustar_index = args.nustar_index.resolve(strict=True)
    with args.reads.open() as source:
        line_count = sum(1 for _ in source)
    if line_count % 4 or line_count == 0:
        parser.error("reads must be a nonempty four-line FASTQ")
    read_count = line_count // 4
    with tempfile.TemporaryDirectory(prefix="stage-timing-equivalence-") as temporary:
        root = Path(temporary)
        sam_outputs = []
        for label, every in (("disabled", None), ("enabled", args.every)):
            prefix = root / f"star.{label}."
            command = [str(args.star), "--runMode", "alignReads", "--runThreadN", "1",
                       "--genomeDir", str(args.star_index), "--readFilesIn", str(args.reads),
                       "--outFileNamePrefix", str(prefix), "--outSAMtype", "SAM",
                       "--outSAMunmapped", "Within"]
            execute(command, environment("STAR_STAGE_TIMING_EVERY", every), prefix)
            timing = Path(str(prefix) + "Log.stage.progress.tsv")
            if every is None and timing.exists():
                raise RuntimeError("STAR wrote timing file while disabled")
            if every is not None:
                verify_timing(timing, read_count)
            sam_outputs.append(sam_records(Path(str(prefix) + "Aligned.out.sam")))
        if sam_outputs[0] != sam_outputs[1]:
            raise RuntimeError("STAR SAM records changed with timing enabled")
        result = {"reads": read_count, "star_sam_records": len(sam_outputs[0]),
                  "star_sam_identical": True}
        if args.nustar_seed_trace:
            traces = []
            for label, every in (("disabled", None), ("enabled", args.every)):
                prefix = root / f"seed.{label}"
                trace = Path(str(prefix) + ".tsv")
                execute([str(args.nustar_seed_trace), str(args.nustar_index),
                         str(args.reads), str(trace)],
                        environment("NUSTAR_STAGE_TIMING_EVERY", every), prefix)
                timing = Path(str(trace) + ".timing.tsv")
                if every is None and timing.exists():
                    raise RuntimeError("NuSTAR wrote timing file while disabled")
                if every is not None:
                    verify_timing(timing, read_count)
                traces.append(trace.read_bytes())
            if traces[0] != traces[1]:
                raise RuntimeError("NuSTAR seed trace changed with timing enabled")
            result["nustar_seed_trace_identical"] = True
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
