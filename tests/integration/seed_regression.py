#!/usr/bin/env python3
"""Compare complete seed-result digests across a seed-policy regression matrix.

Example:
  python3 tests/integration/seed_regression.py \
      --baseline build/before/seed-benchmark \
      --candidate build/after/seed-benchmark --index index \
      --mate1 reads_1.fastq --mate2 reads_2.fastq --reads 1000 \
      --output seed-regression.json

The primary sample runs single-end and, when --mate2 is supplied, paired-end.
Optionally supply --fixture-index and --fixture-fastq for a second annotated
edge-case dataset. It runs both single-end and paired-end with the same FASTQ
as both mates, exercising the mate spacer without requiring another fixture.
Build that fixture index from tests/fixtures/star_rank_probe.fa and its .gtf
with --sjdbOverhang 16 --genomeChrBinNbits 8 --genomeSAindexNbases 4; use
tests/fixtures/seed_edge_reads.fastq as its FASTQ.

Both executables must use the same seed_benchmark.cpp harness. The harness
hashes every ordered output field, including map_marker and
no_good_piece_length. Each invocation performs one timed sweep plus separate
warmup/validation sweeps. Timing is diagnostic only; use benchmark_seed_search.py
for performance measurements. No STAR installation is required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any


FIELDS = ("reads", "input_bases", "paired", "pieces_per_sweep",
          "digest_algorithm", "digest", "policy")


def positive(value: str) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if not 0 < number <= 2**64-1:
        raise argparse.ArgumentTypeError("must be a positive unsigned 64-bit integer")
    return number


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def policies() -> list[tuple[str, list[str]]]:
    result = [("default", [])]
    result += [(f"lmax_{value}", ["--seedSearchLmax", str(value)])
               for value in (1, 10, 25, 100)]
    result += [(f"start_{value}", ["--seedSearchStartLmax", str(value)])
               for value in (0, 1, 37)]
    result += [(f"multimap_{value}", ["--seedMultimapNmax", str(value)])
               for value in (1, 10)]
    result += [("split_map_min_1", ["--seedSplitMin", "1", "--seedMapMin", "1"])]
    result += [("split_map_min_1_start_0", ["--seedSplitMin", "1", "--seedMapMin", "1",
                                           "--seedSearchStartLmax", "0"])]
    result += [("split_map_min_1_start_1", ["--seedSplitMin", "1", "--seedMapMin", "1",
                                           "--seedSearchStartLmax", "1"])]
    return result


def run(command: list[str]) -> dict[str, Any]:
    try:
        process = subprocess.run(command, capture_output=True, text=True,
                                 timeout=300, check=False)
    except subprocess.TimeoutExpired:
        return {"command": command, "error": "benchmark exceeded 300 seconds", "returncode": None}
    result: dict[str, Any] = {"command": command, "returncode": process.returncode,
                              "stdout": process.stdout, "stderr": process.stderr}
    if process.returncode:
        result["error"] = "benchmark returned a failure"
        return result
    try:
        value = json.loads(process.stdout)
        if value["schema"] != "nustar-seed-benchmark-v1":
            raise ValueError("unsupported seed benchmark schema")
        digest = value["digest"]
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("invalid full-result SHA-256 digest")
        if value["digest_algorithm"] != "sha256-full-ordered-result-v1":
            raise ValueError("unsupported output digest algorithm")
        if (value["reads"] < 1 or len(value["rounds"]) != 1 or value["timed_sweeps"] != 1
                or value["rounds"][0]["digest"] != digest
                or value["rounds"][0]["reads"] != value["reads"]):
            raise ValueError("inconsistent sweep count, digest or fragment count")
        result["comparison"] = {field: value[field] for field in FIELDS}
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
        result["error"] = str(error)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path, required=True,
                        help="seed-benchmark linked against the original seed search")
    parser.add_argument("--candidate", type=Path, required=True,
                        help="seed-benchmark linked against the modified seed search")
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--mate1", type=Path, required=True)
    parser.add_argument("--mate2", type=Path)
    parser.add_argument("--reads", type=positive, default=1000,
                        help="maximum fragments per invocation (default: 1000)")
    parser.add_argument("--fixture-index", type=Path)
    parser.add_argument("--fixture-fastq", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if bool(args.fixture_index) != bool(args.fixture_fastq):
        parser.error("--fixture-index and --fixture-fastq must be supplied together")
    binaries = {name: getattr(args, name).resolve(strict=True)
                for name in ("baseline", "candidate")}
    for name, path in binaries.items():
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f"{name} must be an executable file: {path}")
    datasets = [("sample", args.index.resolve(strict=True), args.mate1.resolve(strict=True),
                 args.mate2.resolve(strict=True) if args.mate2 else None)]
    if args.fixture_index:
        fixture = args.fixture_fastq.resolve(strict=True)
        datasets.append(("fixture", args.fixture_index.resolve(strict=True), fixture, fixture))
    for _, index, _, _ in datasets:
        for file in ("index.sa", "reference.bin"):
            if not (index/file).is_file():
                parser.error(f"index is missing {file}: {index}")
    hashes = {name: sha256(path) for name, path in binaries.items()}
    report: dict[str, Any] = {
        "schema": "nustar-seed-regression-v1", "status": "running",
        "read_limit_per_case": args.reads,
        "binaries": {name: {"path": str(path), "sha256": hashes[name]}
                     for name, path in binaries.items()},
        "datasets": [{"name": name, "index": str(index), "mate1": str(first),
                      "mate2": str(second) if second else None}
                     for name, index, first, second in datasets],
        "validation": "complete ordered result SHA-256, all loaded fragments",
        "cases": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save() -> None:
        args.output.write_text(json.dumps(report, indent=2)+"\n")
    save()
    for dataset, index, first, second in datasets:
        for paired in ((False, True) if second else (False,)):
            for policy, options in policies():
                case: dict[str, Any] = {
                    "name": f"{dataset}_{'pe' if paired else 'se'}_{policy}", "runs": {}}
                for name, path in binaries.items():
                    command = [str(path), str(index), str(first), str(second) if paired else "--",
                               "1", "--sweeps", "1", "--reads", str(args.reads), *options]
                    case["runs"][name] = run(command)
                before = case["runs"]["baseline"].get("comparison")
                after = case["runs"]["candidate"].get("comparison")
                case["equal_full_result_digest"] = bool(before and after and before == after
                                                        and before["reads"] <= args.reads)
                a, b = case["runs"]["baseline"], case["runs"]["candidate"]
                case["matching_rejection"] = (
                    isinstance(a["returncode"], int) and a["returncode"] > 0
                    and a["returncode"] == b["returncode"] and bool(a.get("stderr"))
                    and a.get("stderr") == b.get("stderr") and a.get("stdout") == b.get("stdout"))
                case["passed"] = case["equal_full_result_digest"] or case["matching_rejection"]
                if before and after:
                    case["differences"] = {
                        field: {"baseline": before[field], "candidate": after[field]}
                        for field in FIELDS if before[field] != after[field]}
                report["cases"].append(case)
                save()
                outcome = ("PASS" if case["equal_full_result_digest"] else
                           "MATCHING_REJECTION" if case["matching_rejection"] else "FAIL")
                print(f"{case['name']}: {outcome}", flush=True)
    changed = [name for name, path in binaries.items() if sha256(path) != hashes[name]]
    report["binaries_changed_during_test"] = changed
    report["status"] = "passed" if not changed and all(case["passed"] for case in report["cases"]) else "failed"
    report["case_count"] = len(report["cases"])
    report["passed_case_count"] = sum(case["passed"] for case in report["cases"])
    report["successful_search_case_count"] = sum(case["equal_full_result_digest"] for case in report["cases"])
    report["matching_rejection_case_count"] = sum(case["matching_rejection"] for case in report["cases"])
    report["compared_fragments_across_cases"] = sum(
        case["runs"]["baseline"]["comparison"]["reads"] for case in report["cases"]
        if case["equal_full_result_digest"])
    save()
    print(f"{report['status']}: {report['passed_case_count']}/{report['case_count']} cases; "
          f"report={args.output.resolve()}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except OSError as error:
        print(f"seed_regression: {error}", file=sys.stderr)
        sys.exit(1)
