#include "../annotation/feature_types.hpp"
#include "../barcode/sequence.hpp"
#include "../counting/feature_counter.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"
#include <unistd.h> // for get_current_dir

namespace nustar::star2711b::solo::detail {

void FeatureCounter::outputResults(bool cellFilterYes, std::string outputPrefixMat) {
    // make directories
    createDirectory(outputPrefixMat, context_.runDirPerm, "Solo output directory", context_);
    /* old way, does not work when parent directores are needed
    if (mkdir(outputPrefixMat.c_str(),P.runDirPerm)!=0 && errno!=EEXIST) {//create directory
        exitWithError("EXITING because of fatal OUTPUT FILE error: could not create Solo output
    directory" + outputPrefixMat +
                      "\nSOLUTION: check the path and permisssions\n",
                       std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    */
    // write features.tsv
    switch (featureType) {
    case FeatureTypes::Gene:
    case FeatureTypes::GeneFull:
    case FeatureTypes::GeneFull_Ex50pAS:
    case FeatureTypes::GeneFull_ExonOverIntron:
    case FeatureTypes::Velocyto:
    case FeatureTypes::VelocytoSimple: {
        std::ofstream& geneStr =
            ofstrOpen(outputPrefixMat + options.outFileNames[1], sourceLocation(), context_);
        for (std::uint32_t ii = 0; ii < annotation_.nGe; ii++) {
            geneStr << annotation_.geID[ii] << "\t"
                    << (annotation_.geName[ii].empty() ? annotation_.geID[ii]
                                                       : annotation_.geName[ii]);
            if (options.outFormat.featuresGeneField3 != "-") {
                geneStr << '\t' << options.outFormat.featuresGeneField3;
            };
            geneStr << '\n';
        };
        geneStr.close();
        break;
    };
    case FeatureTypes::SJ:
        std::string sjout; // full path to SJ.out.tab with user-defined prefix
        if (context_.outFileNamePrefix[0] == '/') {
            sjout = context_.outFileNamePrefix + "SJ.out.tab";
        } else {
            char cwd[4096];
            getcwd(cwd, 4096);
            sjout.append(cwd);
            sjout += '/' + context_.outFileNamePrefix + "SJ.out.tab";
        };

        remove((outputPrefixMat + options.outFileNames[1])
                   .c_str()); // remove symlink before creating it
        if (symlink(sjout.c_str(), (outputPrefixMat + options.outFileNames[1]).c_str()) != 0)
            exitWithError("EXITING because of fatal OUTPUT FILE error: could not sym-link " +
                              sjout + " into " + outputPrefixMat + options.outFileNames[1] +
                              "\nSOLUTION: check the path and permisssions\n",
                          std::cerr, context_.inOut->logMain, EXIT_CODE_PARAMETER, context_);
    };
    // write barcodes.tsv
    std::ofstream& cbStr =
        ofstrOpen(outputPrefixMat + options.outFileNames[2], sourceLocation(), context_);
    UnsignedWord nCellGeneEntries =
        0; // total number of non-zero cell/gene combinations (entries in the output matrix)
    if (cellFilterYes) { // filtered cells
        for (std::uint32_t icb = 0; icb < nCB; icb++) {
            if (filteredCells.filtVecBool[icb]) {
                cbStr << options.cbWLstr[indCB[icb]] << '\n';
                nCellGeneEntries += nGenePerCB[icb];
            };
        };
    } else { // unfiltered cells
        for (UnsignedWord ii = 0; ii < options.cbWLsize; ii++) {
            cbStr << options.cbWLstr[ii] << '\n';
        };
        for (std::uint32_t icb = 0; icb < nCB; icb++) {
            nCellGeneEntries += nGenePerCB[icb];
        };
    };
    cbStr.flush();
    // output counting matrix

    for (std::uint32_t iCol = 1; iCol < countMatStride; iCol++) {

        std::string matrixFileName = outputPrefixMat;
        if (featureType == FeatureTypes::Velocyto) {
            const std::array<std::string, 3> velonames = {"spliced.mtx", "unspliced.mtx",
                                                          "ambiguous.mtx"};
            matrixFileName += velonames[iCol - 1];

        } else if (iCol > 1 && cellFilterYes) {
            break; // if cellFilterYes, only output first iCol, since filtering is only done for it

        } else if (options.umiDedup.types.size() > 1) {
            matrixFileName +=
                "umiDedup-" + options.umiDedup.typeNames[options.umiDedup.types[iCol - 1]] + ".mtx";

        } else {
            matrixFileName += options.outFileNames[3];
        };
        std::ofstream& countMatrixStream = ofstrOpen(matrixFileName, sourceLocation(), context_);

        // header
        countMatrixStream << "%%MatrixMarket matrix coordinate integer general\n";
        countMatrixStream << "%\n";
        countMatrixStream << featuresNumber << ' '
                          << (cellFilterYes ? filteredCells.nCells : options.cbWLsize) << ' '
                          << nCellGeneEntries << '\n';

        std::uint32_t cbInd1 = 0;
        for (std::uint32_t icb = 0; icb < nCB; icb++) {
            if (cellFilterYes) {
                if (filteredCells.filtVecBool[icb]) {
                    ++cbInd1;
                } else {
                    continue;
                };
            } else {
                cbInd1 = indCB[icb] + 1;
            };
            for (std::uint32_t ig = 0; ig < nGenePerCB[icb]; ig++) {

                std::uint32_t indG1 = countCellGeneUMIindex[icb] + ig * countMatStride;

                // feature index, CB index, count
                countMatrixStream << countCellGeneUMI[indG1] + 1 << ' ' << cbInd1 << ' '
                                  << countCellGeneUMI[indG1 + iCol] << '\n';
            };
        };

        countMatrixStream.close();
    };
    if (options.multiMap.yes.multi && !cellFilterYes &&
        (featureType == FeatureTypes::Gene || featureType == FeatureTypes::GeneFull ||
         featureType == FeatureTypes::GeneFull_ExonOverIntron ||
         featureType == FeatureTypes::GeneFull_Ex50pAS)) {
        // skipping unique
        nUMIperCBmulti.resize(nCB);
        nGenePerCBmulti.resize(nCB);
        bool nGeneUMIperCBmultiFill = true; // if true, need to fill the arrays

        for (const auto& iMult : options.multiMap.types) {
            for (std::uint32_t iDed = 0; iDed < options.umiDedup.yes.N; iDed++) {
                std::string matrixFileName =
                    outputPrefixMat + "UniqueAndMult-" + options.multiMap.typeNames[iMult];
                if (options.umiDedup.types.size() > 1) {
                    matrixFileName +=
                        "_umiDedup-" + options.umiDedup.typeNames[options.umiDedup.types[iDed]];
                };
                matrixFileName += ".mtx";

                std::uint32_t mIndex = options.multiMap.countInd.I[iMult] + iDed;

                // string matOutString;
                // matOutString.reserve(100000000);
                std::ostringstream matOutStringStream; //(matOutString);

                nCellGeneEntries = 0;

                std::uint32_t cbInd1 = 0;
                for (std::uint32_t icb = 0; icb < nCB; icb++) {
                    cbInd1 = indCB[icb] + 1;

                    /* just unique
                    for (uint32 igm=countCellGeneUMIindex[icb]; igm<countCellGeneUMIindex[icb+1];
                    igm+=countMatStride) {

                        //feature index, CB index, count
                        matOutString += to_string(countCellGeneUMI[igm]+1) +' '+ to_string(cbInd1)
                    +' '+ to_string(countCellGeneUMI[igm+iDed]) + '\n';

                        ++nCellGeneEntries;
                    };
                    */

                    /*//just multiple
                    for (uint32 igm=countMatMult.i[icb]; igm<countMatMult.i[icb+1];
                    igm+=countMatMult.s) { matOutStringStream << int(countMatMult.m[igm])+1 <<' '<<
                    cbInd1 <<' '<< countMatMult.m[igm+mIndex] <<'\n';
                        ++nCellGeneEntries;
                    };
                    */
                    { // sum unique+multiple: go over sorted lists. TODO: use a map to combine
                      // (check efficiency)
                        auto igm1 = countCellGeneUMIindex[icb];
                        auto igm2 = countMatMult.i[icb];
                        while (igm1 < countCellGeneUMIindex[icb + 1] ||
                               igm2 < countMatMult.i[icb + 1]) {
                            std::uint32_t g1, c1 = 0, g2;
                            double c2 = 0;

                            if (igm1 < countCellGeneUMIindex[icb + 1]) {
                                g1 = countCellGeneUMI[igm1];
                                c1 = countCellGeneUMI[igm1 + 1 + iDed];
                            } else {
                                g1 = (std::uint32_t)-1;
                            };

                            if (igm2 < countMatMult.i[icb + 1]) {
                                g2 = countMatMult.m[igm2];
                                c2 = countMatMult.m[igm2 + mIndex];
                            } else {
                                g2 = (std::uint32_t)-1;
                            };

                            if (g1 < g2) { // only unique counts for this gene
                                matOutStringStream << g1 + 1 << ' ' << cbInd1 << ' ' << c1 << '\n';
                                igm1 += countMatStride;
                            } else if (g1 > g2) { // only multiple counts for this gene
                                matOutStringStream << g2 + 1 << ' ' << cbInd1 << ' ' << c2 << '\n';
                                igm2 += countMatMult.s;
                                if (nGeneUMIperCBmultiFill) {
                                    nUMIperCBmulti[icb] += c2;
                                    nGenePerCBmulti[icb]++;
                                };
                            } else { // unique and multiple counts for this gene
                                matOutStringStream << g1 + 1 << ' ' << cbInd1 << ' ' << c1 + c2
                                                   << '\n';
                                igm1 += countMatStride;
                                igm2 += countMatMult.s;
                                if (nGeneUMIperCBmultiFill) {
                                    nUMIperCBmulti[icb] += c2;
                                };
                            };

                            ++nCellGeneEntries;
                        };
                    };
                };
                // matOutStringStream.flush();
                nGeneUMIperCBmultiFill = false; // arrays were filled - only do it once

                std::ofstream& countMatrixStream =
                    ofstrOpen(matrixFileName, sourceLocation(), context_);

                // header
                countMatrixStream << "%%MatrixMarket matrix coordinate real general\n";
                countMatrixStream << "%\n";
                countMatrixStream << featuresNumber << ' ' << options.cbWLsize << ' '
                                  << nCellGeneEntries << '\n';
                countMatrixStream << matOutStringStream.str();

                countMatrixStream.close();
            };
        };
    };
};

} // namespace nustar::star2711b::solo::detail
