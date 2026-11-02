#include "index.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"

namespace nustar::star2711b::solo::detail {

AnnotationIndex::AnnotationIndex(Context& context_in) : context_(context_in) {

    if (!context_.quant.yes)
        return;

    if (!context_.pGe.transform.outQuant) { // standard
        trInfoDir = context_.pGe.sjdbGTFfile == "-"
                        ? context_.pGe.gDir
                        : context_.sjdbInsert.outDir; // if GTF file is given at the mapping stage,
                                                      // it's always used for transcript info
    } else { // transformed genome
        trInfoDir = context_.pGeOut.gDir;
    };

    std::ifstream& geStream = ifstrOpen(trInfoDir + "/geneInfo.tab", sourceLocation(),
                                        "SOLUTION: utilize --sjdbGTFfile /path/to/annotations.gtf "
                                        "option at the genome generation step or mapping step",
                                        context_);
    geStream >> nGe;
    geID.resize(nGe);
    geName.resize(nGe);
    geBiotype.resize(nGe);
    geStream.ignore(999, '\n');
    for (UnsignedWord ii = 0; ii < nGe; ii++) {
        std::string line1;
        std::getline(geStream, line1);
        std::istringstream stream1(line1);
        stream1 >> geID[ii] >> geName[ii] >> geBiotype[ii];
    };
    geStream.close();

    if (context_.quant.trSAM.yes || context_.quant.gene.yes ||
        context_.quant.geneFull_Ex50pAS.yes) { // load exon-transcript structures
        // load tr and ex info
        std::ifstream& trinfo =
            ifstrOpen(trInfoDir + "/transcriptInfo.tab", sourceLocation(),
                      "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the "
                      "genome generation step or mapping step",
                      context_);
        trinfo >> nTr;
        trS = std::make_unique_for_overwrite<UnsignedWord[]>(nTr);
        trE = std::make_unique_for_overwrite<UnsignedWord[]>(nTr);
        trEmax = std::make_unique_for_overwrite<UnsignedWord[]>(nTr);
        trExI = std::make_unique_for_overwrite<std::uint32_t[]>(nTr);
        trExN = std::make_unique_for_overwrite<std::uint16_t[]>(nTr);
        trStr = std::make_unique_for_overwrite<std::uint8_t[]>(nTr);
        trID.resize(nTr);
        trGene = std::make_unique_for_overwrite<std::uint32_t[]>(nTr);
        trLen = std::make_unique_for_overwrite<std::uint32_t[]>(nTr);

        for (std::uint32_t itr = 0; itr < nTr; itr++) {
            std::uint16_t str1;
            trinfo >> trID[itr] >> trS[itr] >> trE[itr] >> trEmax[itr] >> str1 >> trExN[itr] >>
                trExI[itr] >> trGene[itr];
            trStr[itr] = str1;

            if (!trinfo.good()) {
                std::ostringstream errOut;
                errOut << "EXITING because of FATAL GENOME INDEX FILE error: transcriptInfo.tab is "
                          "corrupt, or is incompatible with the current STAR version\n";
                errOut << "SOLUTION: re-generate genome index";
                exitWithError(errOut.str(), std::cerr, context_.inOut->logMain,
                              EXIT_CODE_GENOME_FILES, context_);
            };
        };
        context_.inOut->logMain << "Loaded transcript database, nTr=" << nTr << std::endl;
        trinfo.close();

        std::ifstream& exinfo =
            ifstrOpen(trInfoDir + "/exonInfo.tab", sourceLocation(),
                      "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the "
                      "genome generation step or mapping step",
                      context_);
        exinfo >> nEx;
        exSE = std::make_unique_for_overwrite<std::uint32_t[]>(2 * nEx);
        exLenCum = std::make_unique_for_overwrite<std::uint32_t[]>(nEx);
        for (std::uint32_t iex = 0; iex < nEx; iex++) {
            exinfo >> exSE[2 * iex] >> exSE[2 * iex + 1] >>
                exLenCum[iex]; // reading all elements one after another
        };
        for (std::uint32_t ii = 0; ii < nTr; ii++) {
            std::uint32_t iex1 = trExI[ii] + trExN[ii] - 1; // last exon of the transcript
            trLen[ii] = exLenCum[iex1] + exSE[2 * iex1 + 1] - exSE[2 * iex1] + 1;
        };

        context_.inOut->logMain << "Loaded exon database, nEx=" << nEx << std::endl;
        exinfo.close();
    };
    // load exon-gene structures
    if (context_.quant.geCount.yes) {
        std::ifstream& exinfo =
            ifstrOpen(trInfoDir + "/exonGeTrInfo.tab", sourceLocation(),
                      "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the "
                      "genome generation step or mapping step",
                      context_);
        exinfo >> exG.nEx;
        exG.s = std::make_unique_for_overwrite<UnsignedWord[]>(exG.nEx);
        exG.e = std::make_unique_for_overwrite<UnsignedWord[]>(exG.nEx);
        exG.eMax = std::make_unique_for_overwrite<UnsignedWord[]>(exG.nEx);
        exG.str = std::make_unique_for_overwrite<std::uint8_t[]>(exG.nEx);
        exG.g = std::make_unique_for_overwrite<std::uint32_t[]>(exG.nEx);
        exG.t = std::make_unique_for_overwrite<std::uint32_t[]>(exG.nEx);
        for (UnsignedWord ii = 0; ii < exG.nEx; ii++) {
            int str1;
            exinfo >> exG.s[ii] >> exG.e[ii] >> str1 >> exG.g[ii] >> exG.t[ii];
            exG.str[ii] = (std::uint8_t)str1;
        };
        exinfo.close();
        // calculate eMax
        exG.eMax[0] = exG.e[0];
        for (UnsignedWord iex = 1; iex < exG.nEx; iex++) {
            exG.eMax[iex] = std::max(exG.eMax[iex - 1], exG.e[iex]);
        };
    };

    if (context_.quant.geneFull.yes || context_.quant.geneFull_ExonOverIntron.yes) {
        std::ifstream& exinfo =
            ifstrOpen(trInfoDir + "/exonGeTrInfo.tab", sourceLocation(),
                      "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the "
                      "genome generation step or mapping step",
                      context_);
        exinfo >> exG.nEx;

        geneFull.s = std::make_unique_for_overwrite<UnsignedWord[]>(nGe);
        geneFull.e = std::make_unique_for_overwrite<UnsignedWord[]>(nGe);
        geneFull.eMax = std::make_unique_for_overwrite<UnsignedWord[]>(nGe);
        geneFull.g = std::make_unique_for_overwrite<std::uint32_t[]>(nGe);
        geneFull.str = std::make_unique_for_overwrite<std::uint8_t[]>(nGe);

        for (UnsignedWord ig = 0; ig < nGe; ig++) {
            geneFull.s[ig] = -1; // largest uint64
            geneFull.e[ig] = 0;
        };

        for (UnsignedWord ii = 0; ii < exG.nEx; ii++) {
            UnsignedWord s1, e1, str1, g1, t1;
            exinfo >> s1 >> e1 >> str1 >> g1 >> t1;
            geneFull.s[g1] = std::min(geneFull.s[g1], s1);
            geneFull.e[g1] = std::max(geneFull.e[g1], e1);
            geneFull.str[g1] = (std::uint8_t)str1;
        };
        exinfo.close();

        auto gF = std::make_unique_for_overwrite<UnsignedWord[]>(4 * nGe);
        for (UnsignedWord ii = 0; ii < nGe; ii++) {
            gF[4 * ii] = geneFull.s[ii];
            gF[4 * ii + 1] = geneFull.e[ii];
            gF[4 * ii + 2] = geneFull.str[ii];
            gF[4 * ii + 3] = ii;
        };

        qsort(gF.get(), nGe, sizeof(UnsignedWord) * 4, funCompareArrays<UnsignedWord, 2>);

        for (UnsignedWord ii = 0; ii < nGe; ii++) {
            geneFull.s[ii] = gF[4 * ii];
            geneFull.e[ii] = gF[4 * ii + 1];
            geneFull.str[ii] = gF[4 * ii + 2];
            geneFull.g[ii] = gF[4 * ii + 3];
        };

        // calculate eMax
        geneFull.eMax[0] = geneFull.e[0];
        for (UnsignedWord iex = 1; iex < nGe; iex++) {
            geneFull.eMax[iex] = std::max(geneFull.eMax[iex - 1], geneFull.e[iex]);
        };
    };
};

} // namespace nustar::star2711b::solo::detail
