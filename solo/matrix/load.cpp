#include "../annotation/feature_types.hpp"
#include "../counting/feature_counter.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"

namespace nustar::star2711b::solo::detail {

void FeatureCounter::loadRawMatrix() {
    // make directories
    if (context_.runModeIn.size() < 3) {
        std::string errOut =
            "Exiting because of fatal PARAMETER error: --runMode soloCellFiltering should contain "
            "paths to count matrix input directorry and output prefix.";
        errOut += "\nSOLUTION: re-run with --runMode soloCellFiltering </path/to/raw/count/dir/> "
                  "</path/to/output/prefix>\n";
        exitWithError(errOut, std::cerr, context_.inOut->logMain, EXIT_CODE_PARAMETER, context_);
    };

    std::string inputPrefix = context_.runModeIn[1] + '/';
    outputPrefix = context_.runModeIn[2];
    outputPrefixFiltered = outputPrefix;
    // load counting matrix
    std::string matrixFileName = inputPrefix + options.outFileNames[3];
    std::ifstream& matStream = ifstrOpen(
        matrixFileName, sourceLocation(),
        "SOLUTION: check path and permission for the matrix file" + matrixFileName, context_);

    // header
    while (matStream.peek() == '%') {
        std::string dummy;
        matStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    };

    std::uint32_t nCB1; // number of features read from file
    UnsignedWord nTot;  // total number of entries
    matStream >> featuresNumber >> nCB1 >> nTot;

    if (nTot == 0) {
        exitWithError("Exiting because of fatal INPUT FILE error: no counts detected in " +
                          matrixFileName + "\nSOLUTION: check the formatting of the matrix file.\n",
                      std::cerr, context_.inOut->logMain, EXIT_CODE_PARAMETER, context_);
    };

    /* do not need it - read the information from features file
    if (nfeat1 != featuresNumber) {
        ostringstream errOut;
        errOut <<"Exiting because of fatal INPUT FILE error: number of features in "<<
    matrixFileName <<": "<< nfeat1; errOut <<"\n   is not equal to the number of features from
    annotations from genome index: " << featuresNumber << "\n"; errOut <<"SOLUTION: use the same
    genome index that was use in the STARsolo run that generated the matrix\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    */

    // struct {uint32 gene; uint32 cb; uint32 count} matEntry;
    /*vector<array<uint32,3>> matSparse(nTot);
    for (uint64 ii=0; ii<nTot; ii++) {
        for (uint32 jj=0; jj<3; jj++)
            matStream >> matSparse[ii][jj];
    };

    //sort by cell index
    std::sort(matSparse.begin(), matSparse.end(), [](const array<uint32,3> &m1, const
    array<uint32,3> &m2) { return m1[1] < m2[1] ;
                                                    });
    */

    countMatStride = 3; // gene, cell, count. Recording cell at shift=1 is temporary: later will
                        // replace cell with count
    countCellGeneUMI.resize(nTot * countMatStride, 0);

    // countMatMult.s = 3;
    // countMatMult.m.resize(nTot*countMatStride,0.0);

    for (UnsignedWord ii = 0; ii < nTot; ii++) {
        matStream >> countCellGeneUMI[ii * countMatStride + 0]; // gene
        --countCellGeneUMI[ii * countMatStride + 0];            // 0-based gene

        matStream >> countCellGeneUMI[ii * countMatStride + 1]; // cell
        --countCellGeneUMI[ii * countMatStride + 1];            // 0-based cell

        double count1;
        matStream >> count1;
        countCellGeneUMI[ii * countMatStride + 2] = std::round(count1);

        /* to keep fractional part
        matStream >> countMatMult.m[ii*countMatStride+2];//count: for now, only allow one value per
        cell/gene countCellGeneUMI[ii*countMatStride+2] =
        std::round(countMatMult.m[ii*countMatStride+2]); //round to integer
        countMatMult.m[ii*countMatStride+2] -= countCellGeneUMI[ii*countMatStride+2]; //just the
        fractional part countMatMult.m[ii*countMatStride+0] = countCellGeneUMI[ii*countMatStride+0];
        countMatMult.m[ii*countMatStride+1] = countCellGeneUMI[ii*countMatStride+1];
        */
    };

    qsort((void*)countCellGeneUMI.data(), nTot, countMatStride * sizeof(countCellGeneUMI[0]),
          funCompareTypeSecondFirst<std::uint32_t>);

    // count number of detected cell
    nCB = 0;
    std::uint32_t ciprev = (std::uint32_t)-1;
    for (std::uint32_t ii = 0; ii < nTot; ii++) {
        std::uint32_t ci1 = countCellGeneUMI[ii * countMatStride + 1];
        if (ci1 != ciprev) { // new cell
            ciprev = ci1;
            nCB++;
        };
    };

    indCB.resize(nCB);
    countCellGeneUMIindex.resize(nCB);
    nUMIperCB.resize(nCB, 0);
    nGenePerCB.resize(nCB, 0);
    nReadPerCB.resize(nCB, 0);

    nCB = (std::uint32_t)-1;
    ciprev = (std::uint32_t)-1;
    for (std::uint32_t ii = 0; ii < nTot; ii++) {
        std::uint32_t ci1 = countCellGeneUMI[ii * countMatStride + 1];
        if (ci1 != ciprev) { // new cell
            ciprev = ci1;
            nCB++;
            indCB[nCB] = ci1;
            countCellGeneUMIindex[nCB] = ii * countMatStride;
        };
        nGenePerCB[nCB]++;
        nUMIperCB[nCB] += countCellGeneUMI[ii * countMatStride + 2];
        countCellGeneUMI[ii * countMatStride + 1] =
            countCellGeneUMI[ii * countMatStride + 2]; // replace cell with count to keep standard
                                                       // convention about countCellGeneUMI
    };

    { // load barcodes
        std::ifstream& wlstream =
            ifstrOpen(inputPrefix + options.outFileNames[2], sourceLocation(),
                      "SOLUTION: check the path and permissions of the barcodes file", context_);
        options.cbWLstr.resize(nCB1);
        for (auto& cb : options.cbWLstr)
            std::getline(wlstream, cb);
    };

    { // copy features
        std::ifstream& infeat =
            ifstrOpen(inputPrefix + options.outFileNames[1], sourceLocation(),
                      "SOLUTION: check the path and permissions of the features file", context_);
        createDirectory(outputPrefixFiltered, context_.runDirPerm, "Solo output directory",
                        context_);
        std::ofstream& outfeat =
            ofstrOpen(outputPrefixFiltered + options.outFileNames[1], sourceLocation(), context_);
        outfeat << infeat.rdbuf();
        outfeat.close();
    };

    return;
};

} // namespace nustar::star2711b::solo::detail
