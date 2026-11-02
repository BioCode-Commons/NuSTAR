#include "feature_counter.hpp"
#include "../support/streams.hpp"
#include <numeric>

namespace nustar::star2711b::solo::detail {
////
void FeatureCounter::redistributeReadsByCB() { // redistribute reads in files by CB - each file with
                                               // the approximately the same number of reads, each
                                               // CB is on one file only

    /* SoloFeature vars that have to be setup:
     * nCB
     * readFeatSum->cbReadCount[]
     */

    // find boundaries for cells
    UnsignedWord nReadRec =
        std::accumulate(readFeatSum->cbReadCount.begin(), readFeatSum->cbReadCount.end(), 0LLU);
    //    nReadRec += cbrc;

    UnsignedWord nReadRecBin = nReadRec / options.redistrReadsNfiles;

    context_.inOut->logMain << "     Redistributing reads into " << options.redistrReadsNfiles
                            << "files; nReadRec=" << nReadRec << ";   nReadRecBin=" << nReadRecBin
                            << std::endl;

    redistrFilesCBfirst.push_back(0);
    redistrFilesCBindex.resize(nCB);
    UnsignedWord nreads = 0;
    std::uint32_t ind = 0;
    for (std::uint32_t icb = 0; icb < nCB; icb++) {
        redistrFilesCBindex[icb] = ind;
        nreads += readFeatSum->cbReadCount[indCB[icb]];
        if (nreads >= nReadRecBin) {
            ind++;
            redistrFilesCBfirst.push_back(icb + 1);
            redistrFilesNreads.push_back(nreads);
            nreads = 0;
        };
    };
    if (nreads > 0) {
        redistrFilesCBfirst.push_back(nCB);
        redistrFilesNreads.push_back(nreads);
    };

    // open output files
    redistrFilesStreams.resize(redistrFilesNreads.size());
    for (std::uint32_t ii = 0; ii < redistrFilesNreads.size(); ii++) {
        // open file with flagDelete=true
        redistrFilesStreams[ii] =
            &fstrOpen(context_.outFileTmp + "solo" + FeatureTypes::Names[featureType] +
                          "_redistr_" + std::to_string(ii),
                      sourceLocation(), context_, true);
    };

    // main cycle
    for (int ii = 0; ii < context_.runThreadN; ii++) {
        readFeatAll[ii]->streamReads->clear(); // this is needed if eof was reached before
        readFeatAll[ii]->streamReads->seekg(0, std::ios::beg);

        while (true) {
            std::string line1;
            std::getline(*readFeatAll[ii]->streamReads, line1);
            if (line1.empty()) {
                break;
            };

            std::istringstream line1stream(line1);
            UnsignedWord cb1, umi;
            line1stream >> umi >> cb1 >> cb1;
            if (featureType == FeatureTypes::SJ)
                line1stream >> cb1;
            line1stream >> cb1;

            *redistrFilesStreams[redistrFilesCBindex[indCBwl[cb1]]] << line1 << '\n';
        };
        // TODO: delete streamReads files one by one to save disk space
    };

    // close files
    //     redistrFilesStreams[ii]->flush();
};

} // namespace nustar::star2711b::solo::detail
