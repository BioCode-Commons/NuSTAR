#!/usr/bin/env python3
"""Benchmark STAR and NuSTAR GenomeGenerate on the same FASTA/GTF inputs."""

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time


def run_one(label: str, command: list[str], directory: Path) -> dict:
    directory.mkdir()
    log = directory / "run.log"
    with log.open("w") as output:
        started = time.perf_counter()
        with subprocess.Popen(command, cwd=directory, stdout=output,
                              stderr=subprocess.STDOUT) as process:
            _, status, usage = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
        elapsed = time.perf_counter() - started
    if process.returncode:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-20:])
        raise RuntimeError(f"{label} failed ({process.returncode}):\n{tail}")
    return {
        "label": label,
        "command": command,
        "wall_seconds": elapsed,
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
        "max_rss_kib": usage.ru_maxrss,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--nustar", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--overhang", type=int, default=100)
    parser.add_argument("--chr-bin-bits", type=int, default=18)
    parser.add_argument("--prefix-bases", type=int, default=10)
    parser.add_argument("--runThreadN", "--star-threads", dest="threads", type=int, default=1,
                        help="threads for both programs (--star-threads is a legacy alias)")
    parser.add_argument("--limitGenomeGenerateRAM", type=int, default=40000000000,
                        help="STAR construction RAM setting in bytes")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.repeat < 1 or args.threads < 1 or args.limitGenomeGenerateRAM < 1:
        parser.error("repeat, runThreadN and limitGenomeGenerateRAM must be positive")
    star, nustar, fasta = (p.resolve(strict=True) for p in
                           (args.star, args.nustar, args.fasta))
    gtf = args.gtf.resolve(strict=True) if args.gtf else None
    records = []
    with tempfile.TemporaryDirectory(prefix="nustar-genome-benchmark-") as temp:
        root = Path(temp)
        for iteration in range(args.repeat):
            labels = ("STAR", "NuSTAR") if iteration % 2 == 0 else ("NuSTAR", "STAR")
            for label in labels:
                directory = root / f"{iteration}-{label.lower()}"
                if label == "STAR":
                    command = [str(star), "--runMode", "genomeGenerate",
                               "--genomeFastaFiles", str(fasta), "--genomeDir", str(directory),
                               "--sjdbOverhang", str(args.overhang)] if gtf else [
                               str(star), "--runMode", "genomeGenerate",
                               "--genomeFastaFiles", str(fasta), "--genomeDir", str(directory)]
                    command += ["--genomeChrBinNbits", str(args.chr_bin_bits),
                                "--genomeSAindexNbases", str(args.prefix_bases),
                                "--runThreadN", str(args.threads),
                                "--limitGenomeGenerateRAM", str(args.limitGenomeGenerateRAM)]
                else:
                    command = [str(nustar), "--runMode", "genomeGenerate", "--genomeFastaFiles",
                               str(fasta), "--genomeDir", str(directory),
                               "--genomeChrBinNbits", str(args.chr_bin_bits),
                               "--genomeSAindexNbases", str(args.prefix_bases),
                               "--runThreadN", str(args.threads)]
                    if gtf:
                        command += ["--sjdbOverhang", str(args.overhang)]
                if gtf:
                    command += ["--sjdbGTFfile", str(gtf)]
                result = run_one(label, command, directory)
                result["iteration"] = iteration + 1
                result["threads"] = args.threads
                records.append(result)
                print(f"{label} {iteration+1}: {result['wall_seconds']:.2f}s, "
                      f"{result['max_rss_kib']/1024:.0f} MiB peak RSS", flush=True)
    for label in ("STAR", "NuSTAR"):
        values = [x["wall_seconds"] for x in records if x["label"] == label]
        rss = [x["max_rss_kib"] for x in records if x["label"] == label]
        print(f"{label} median: {statistics.median(values):.2f}s, "
              f"max RSS: {max(rss)/1024:.0f} MiB")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps({
            "inputs": {"fasta": str(fasta), "gtf": str(gtf) if gtf else None,
                       "overhang": args.overhang, "chr_bin_bits": args.chr_bin_bits,
                       "prefix_bases": args.prefix_bases,
                       "runThreadN": args.threads,
                       "limitGenomeGenerateRAM": args.limitGenomeGenerateRAM},
            "records": records,
        }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
