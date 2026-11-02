#!/usr/bin/env python3
"""Build a trace-only STAR copy for single-end window stitching comparison."""

import argparse
from pathlib import Path
import subprocess
import sys


BEFORE = r'''
            const bool traceStitch = std::getenv("STAR_STITCH_TRACE") != nullptr;
            if (traceStitch) {
                std::fprintf(stderr, "STITCH_WIN\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                    readName, (unsigned long long)iW,
                    (unsigned long long)WC[iW][WC_Chr],
                    (unsigned long long)WC[iW][WC_Str],
                    (unsigned long long)trA.roStr,
                    (unsigned long long)Lread,
                    (unsigned long long)outFilterMismatchNmaxTotal,
                    (unsigned long long)nWA[iW],
                    (unsigned long long)WlastAnchor[iW]);
                for (uint traceSeed = 0; traceSeed < nWA[iW]; ++traceSeed) {
                    std::fprintf(stderr,
                        "STITCH_WA\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                        readName, (unsigned long long)iW,
                        (unsigned long long)traceSeed,
                        (unsigned long long)WA[iW][traceSeed][WA_rStart],
                        (unsigned long long)WA[iW][traceSeed][WA_gStart],
                        (unsigned long long)WA[iW][traceSeed][WA_Length],
                        (unsigned long long)WA[iW][traceSeed][WA_Nrep],
                        (unsigned long long)WA[iW][traceSeed][WA_Anchor],
                        (unsigned long long)WA[iW][traceSeed][WA_iFrag],
                        (unsigned long long)WA[iW][traceSeed][WA_sjA]);
                }
            }
            const auto traceStitchStart = std::chrono::steady_clock::now();
'''

READ = r'''
    if (std::getenv("STAR_STITCH_TRACE") != nullptr) {
        std::fprintf(stderr, "STITCH_READ\t%s\t%llu\t%llu\n",
            readName, (unsigned long long)Lread,
            (unsigned long long)outFilterMismatchNmaxTotal);
    }
'''

AFTER = r'''
            if (traceStitch) {
                const auto traceNanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - traceStitchStart).count();
                std::fprintf(stderr, "STITCH_OUT\t%s\t%llu\t%llu\t%lld\n",
                    readName, (unsigned long long)iW,
                    (unsigned long long)nWinTr[iW1], (long long)traceNanoseconds);
                for (uint traceCandidate = 0; traceCandidate < nWinTr[iW1]; ++traceCandidate) {
                    const Transcript& traceTr = *trAll[iW1][traceCandidate];
                    std::fprintf(stderr,
                        "STITCH_TR\t%s\t%llu\t%llu\t%d\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                        readName, (unsigned long long)iW,
                        (unsigned long long)traceCandidate, traceTr.maxScore,
                        (unsigned long long)traceTr.rStart,
                        (unsigned long long)traceTr.roStart,
                        (unsigned long long)traceTr.rLength,
                        (unsigned long long)traceTr.gStart,
                        (unsigned long long)traceTr.gLength,
                        (unsigned long long)traceTr.nExons,
                        (unsigned long long)traceTr.nMatch,
                        (unsigned long long)traceTr.nMM,
                        (unsigned long long)traceTr.mappedLength,
                        (unsigned long long)traceTr.nGap,
                        (unsigned long long)traceTr.lGap,
                        (unsigned long long)traceTr.nDel,
                        (unsigned long long)traceTr.lDel,
                        (unsigned long long)traceTr.nIns,
                        (unsigned long long)traceTr.lIns,
                        (unsigned long long)traceTr.nAnchor);
                    std::fprintf(stderr,
                        "STITCH_META\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%u\t%u\t%d\n",
                        readName, (unsigned long long)iW,
                        (unsigned long long)traceCandidate,
                        (unsigned long long)traceTr.nUnique,
                        (unsigned long long)traceTr.extendL,
                        (unsigned long long)traceTr.intronMotifs[0],
                        (unsigned long long)traceTr.intronMotifs[1],
                        (unsigned long long)traceTr.intronMotifs[2],
                        (unsigned)traceTr.sjMotifStrand,
                        (unsigned)traceTr.sjYes, traceTr.iFrag);
                    for (uint traceExon = 0; traceExon < traceTr.nExons; ++traceExon) {
                        std::fprintf(stderr,
                            "STITCH_EXON\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                            readName, (unsigned long long)iW,
                            (unsigned long long)traceCandidate,
                            (unsigned long long)traceExon,
                            (unsigned long long)traceTr.exons[traceExon][EX_R],
                            (unsigned long long)traceTr.exons[traceExon][EX_G],
                            (unsigned long long)traceTr.exons[traceExon][EX_L],
                            (unsigned long long)traceTr.exons[traceExon][EX_iFrag],
                            (unsigned long long)traceTr.exons[traceExon][EX_sjA]);
                        if (traceExon + 1 < traceTr.nExons) {
                            std::fprintf(stderr,
                                "STITCH_JUNCTION\t%s\t%llu\t%llu\t%llu\t%d\t%u\t%u\t%llu\t%llu\n",
                                readName, (unsigned long long)iW,
                                (unsigned long long)traceCandidate,
                                (unsigned long long)traceExon,
                                traceTr.canonSJ[traceExon],
                                (unsigned)traceTr.sjAnnot[traceExon],
                                (unsigned)traceTr.sjStr[traceExon],
                                (unsigned long long)(traceTr.canonSJ[traceExon]!=-2 ? traceTr.shiftSJ[traceExon][0] : 0),
                                (unsigned long long)(traceTr.canonSJ[traceExon]!=-2 ? traceTr.shiftSJ[traceExon][1] : 0));
                        }
                    }
                }
            }
'''

BEST = r'''
    if (std::getenv("STAR_STITCH_TRACE") != nullptr) {
        std::fprintf(stderr, "STITCH_BEST\t%s\t%d\t%llu\t%llu\t%llu\t%llu\t%llu\n",
            readName, trBest->maxScore,
            (unsigned long long)trBest->gStart,
            (unsigned long long)trBest->gLength,
            (unsigned long long)trBest->Chr,
            (unsigned long long)trBest->Str,
            (unsigned long long)nW);
    }
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    builder = Path(__file__).with_name("build_star_seed_oracle.py")
    subprocess.run([sys.executable, str(builder), "--star-source",
                    str(args.star_source), "--output", str(args.output),
                    "--jobs", str(args.jobs)], check=True)
    output = args.output.resolve(strict=True)
    source = output / "ReadAlign_stitchPieces.cpp"
    body = source.read_text()
    read_anchor = "void ReadAlign::stitchPieces(char **R, uint Lread) {"
    if body.count(read_anchor) != 1:
        raise SystemExit("STAR stitch read anchor changed")
    body = body.replace(read_anchor, read_anchor + "\n" + READ, 1)
    anchor = "            stitchWindowAligns(0, nWA[iW], 0, WAincl, 0, 0, trA, Lread, WA[iW], R[trA.roStr==0 ? 0:2], mapGen, P, trAll[iW1], nWinTr+iW1, this);"
    if body.count(anchor) != 1:
        raise SystemExit("STAR stitch call anchor changed")
    body = body.replace(anchor, BEFORE + anchor + AFTER, 1)
    best_anchor = "    nW=iW1;//only count windows that had alignments"
    if body.count(best_anchor) != 1:
        raise SystemExit("STAR stitch best anchor changed")
    body = body.replace(best_anchor, best_anchor + "\n" + BEST, 1)
    source.write_text(body)
    make_options = []
    compatibility = output / "seed_oracle_compat"
    if compatibility.is_dir():
        make_options = [f"CPPFLAGS=-I{compatibility}",
                        f"CXXFLAGSextra=-I{compatibility}",
                        f"LDFLAGSextra=-L{compatibility}"]
    subprocess.run(["make", "-C", str(output), f"-j{args.jobs}",
                    "-o", "Depend.list", "-o", "parametersDefault.xxd",
                    *make_options, "STAR"], check=True)
    print(output / "STAR")


if __name__ == "__main__":
    main()
