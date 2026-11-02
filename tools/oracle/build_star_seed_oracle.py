#!/usr/bin/env python3
"""Build a seed-tracing STAR in a separate source copy.

This script never edits the baseline STAR tree. The copied object files allow a
single instrumented translation unit to be rebuilt before the normal link.
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess


TRACE = r'''
    if (std::getenv("STAR_SEED_TRACE") != nullptr) {
        struct timespec seed_time_end;
        clock_gettime(CLOCK_MONOTONIC, &seed_time_end);
        unsigned long long seed_ns =
            seed_quality_ns +
            (unsigned long long)(seed_time_end.tv_sec-seed_time_start.tv_sec)*1000000000ULL +
            (unsigned long long)(seed_time_end.tv_nsec-seed_time_start.tv_nsec);
        std::fprintf(stderr, "SEED_TIME\t%s\t%llu\n", readName, seed_ns);
        std::fprintf(stderr,
            "SEED_SUM\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
            readName, (unsigned long long)Lread, (unsigned long long)Nsplit,
            (unsigned long long)nP, (unsigned long long)nA,
            (unsigned long long)nUM[0], (unsigned long long)nUM[1],
            (unsigned long long)multNmin, (unsigned long long)multNminL,
            (unsigned long long)uniqLmax);
        for (uint it = 0; it < Nsplit; ++it)
            std::fprintf(stderr, "SEED_SPLIT\t%s\t%llu\t%llu\t%llu\t%llu\n",
                readName, (unsigned long long)it,
                (unsigned long long)splitR[0][it],
                (unsigned long long)splitR[1][it],
                (unsigned long long)splitR[2][it]);
        for (uint it = 0; it < nP; ++it)
            std::fprintf(stderr,
                "SEED_PC\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                readName, (unsigned long long)it,
                (unsigned long long)PC[it][PC_rStart],
                (unsigned long long)PC[it][PC_Length],
                (unsigned long long)PC[it][PC_Dir],
                (unsigned long long)PC[it][PC_Nrep],
                (unsigned long long)PC[it][PC_SAstart],
                (unsigned long long)PC[it][PC_SAend],
                (unsigned long long)PC[it][PC_iFrag]);
    }
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    source = args.star_source.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"output already exists: {output}")
    shutil.copytree(source, output)
    path = output / "ReadAlign_mapOneRead.cpp"
    body = path.read_text()
    start_anchor = "int ReadAlign::mapOneRead() {"
    if body.count(start_anchor) != 1:
        raise SystemExit("STAR seed timing anchor changed")
    body = body.replace(start_anchor,
        start_anchor + '\n    struct timespec seed_time_start, seed_quality_end;\n'
        '    unsigned long long seed_quality_ns = 0;\n'
        '    clock_gettime(CLOCK_MONOTONIC, &seed_time_start);', 1)
    quality_anchor = "    resetN(); //reset aligns counters to 0"
    if body.count(quality_anchor) != 1:
        raise SystemExit("STAR quality split timing anchor changed")
    quality_timing = (
        '    clock_gettime(CLOCK_MONOTONIC, &seed_quality_end);\n'
        '    seed_quality_ns =\n'
        '        (unsigned long long)(seed_quality_end.tv_sec-seed_time_start.tv_sec)*1000000000ULL +\n'
        '        (unsigned long long)(seed_quality_end.tv_nsec-seed_time_start.tv_nsec);\n'
    )
    body = body.replace(quality_anchor, quality_timing + quality_anchor, 1)
    search_anchor = "    trBest=trInit;"
    if body.count(search_anchor) != 1:
        raise SystemExit("STAR seed search timing anchor changed")
    body = body.replace(search_anchor,
        search_anchor + '\n    clock_gettime(CLOCK_MONOTONIC, &seed_time_start);', 1)
    anchor = "    #ifdef OFF_AFTER_SEEDING"
    if body.count(anchor) != 1:
        raise SystemExit("STAR seed tracing anchor changed")
    body = body.replace(anchor, TRACE + "\n" + anchor, 1)
    # This checkout already contains STAR_DEBUG_PIECES diagnostics. Disable
    # their getenv/printing branch in the oracle while retaining every
    # maxMappableLength2strands call and its PC side effects.
    body = body.replace('getenv("STAR_DEBUG_PIECES")', 'nullptr')
    path.write_text(body)
    make_options = []
    if not Path("/usr/include/zlib.h").exists():
        # Only the already compiled BGZF objects use zlib. This translation unit
        # needs the opaque pointer declaration from bgzf.h and a linker soname.
        compatibility = output / "seed_oracle_compat"
        compatibility.mkdir()
        (compatibility / "zlib.h").write_text(
            "struct z_stream_s; typedef struct z_stream_s z_stream;\n")
        libz = Path(subprocess.check_output(
            ["g++", "-print-file-name=libz.so.1"], text=True).strip())
        if not libz.is_file():
            raise SystemExit("cannot find the existing zlib runtime")
        (compatibility / "libz.so").symlink_to(libz)
        dependencies = output / "Depend.list"
        if dependencies.exists():
            dependencies.write_text(re.sub(
                r"/[^\s\\]*?/zlib\.h", str(compatibility / "zlib.h"),
                dependencies.read_text()))
        make_options = [f"CPPFLAGS=-I{compatibility}",
                        f"CXXFLAGSextra=-I{compatibility}",
                        f"LDFLAGSextra=-L{compatibility}"]
    subprocess.run(["make", "-C", str(output), f"-j{args.jobs}",
                    "-o", "Depend.list", "-o", "parametersDefault.xxd",
                    *make_options, "STAR"],
                   check=True)
    print(output / "STAR")


if __name__ == "__main__":
    main()
