#include "../barcode/sequence.hpp"
#include "feature_counter.hpp"
#include "input_record.hpp"
#include "../support/algorithms.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"

namespace nustar::star2711b::solo::detail {

void FeatureCounter::countSmartSeq() {
    time_t rawTime;

    // nCB=pSolo.cbWLsize; //all cells are recorded
    // nCB, indCB are recorded in sumThreads

    for (int ii = 0; ii < context_.runThreadN; ii++) {
        readFeatSum->addStats(*readFeatAll[ii]);
    };

    redistributeReadsByCB();

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Finished redistribution of reads from Solo read files"
                            << std::endl;
    nReadPerCB.resize(nCB);
    cbFeatureUMImap.resize(nCB);

    typedef struct {
        std::uint32_t feature;
        std::uint32_t count[2]; //[0]=NoDedup, [1]=Exact
    } typeFeatureCount;

    std::vector<std::vector<typeFeatureCount>> vCellFeatureCount(nCB);

    typedef struct {
        std::uint32_t feature;
        UnsignedWord umi;
    } typeFeatureUMI;
    std::vector<std::vector<typeFeatureUMI>::iterator> cbFeatUMI(nCB + 1);

#pragma omp parallel for num_threads(context_.runThreadN)
    for (std::uint32_t ired = 0; ired < redistrFilesCBfirst.size() - 1;
         ired++) { // TODO: parallelize, each ired is independent here!
        std::vector<typeFeatureUMI> vFeatureUMI(redistrFilesNreads[ired]);

        std::uint32_t iCB1 = redistrFilesCBfirst[ired];
        std::uint32_t iCB2 = redistrFilesCBfirst[ired + 1];

        // allocate arrays
        cbFeatUMI[iCB1] = vFeatureUMI.begin();
        for (std::uint32_t icb = iCB1 + 1; icb < iCB2; icb++) {
            cbFeatUMI[icb] = cbFeatUMI[icb - 1] + readFeatSum->cbReadCount[indCB[icb - 1]];
        };

        // input records
        redistrFilesStreams[ired]->flush();
        redistrFilesStreams[ired]->seekg(0, std::ios::beg);

        std::uint32_t feature;
        UnsignedWord umi, iread;
        std::int32_t cbmatch;
        SignedWord cb;
        std::vector<std::uint32_t> trIdDist; // not used
        bool readInfoYes = false;
        while (soloInputFeatureUMI(redistrFilesStreams[ired], featureType, readInfoYes,
                                   context_.sjAll, iread, cbmatch, feature, umi, trIdDist,
                                   readFlagCounts)) { // cycle over file records

            *redistrFilesStreams[ired] >> cb;

            if (feature + 1 == 0)
                continue;

            std::uint32_t icb = indCBwl[cb];
            *(cbFeatUMI[icb] + nReadPerCB[icb]) = {feature, umi};
            nReadPerCB[icb]++;
        };

        // collapse UMI, simple
        for (std::uint32_t icb = iCB1; icb < iCB2; icb++) {
            if (nReadPerCB[icb] == 0) {
                continue;
            };

            std::sort(cbFeatUMI[icb], cbFeatUMI[icb] + nReadPerCB[icb],
                      [](const typeFeatureUMI& a, const typeFeatureUMI& b) {
                          return (a.feature < b.feature) ||
                                 (a.feature == b.feature && a.umi < b.umi);
                      });

            vCellFeatureCount[icb].reserve(8192);
            vCellFeatureCount[icb].push_back({cbFeatUMI[icb]->feature, {1, 1}}); // first read
            for (auto fu = cbFeatUMI[icb] + 1; fu != cbFeatUMI[icb] + nReadPerCB[icb];
                 fu++) {                                // cycle over all reads for this icb
                if (fu->feature != (fu - 1)->feature) { // compare to previous feature
                    vCellFeatureCount[icb].push_back(
                        {fu->feature, {1, 1}});               // create next feature entry
                } else {                                      // same feature
                    vCellFeatureCount[icb].back().count[0]++; // non-collapsed UMI count

                    if (fu->umi != (fu - 1)->umi) {               // same feature, new umi
                        vCellFeatureCount[icb].back().count[1]++; // collapsed UMI count
                    };
                };
            };
        };
    };

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading / collapsing"
                            << std::endl;

    // convert to countCellGeneUMI. TODO - this is not necessary, recode output to use
    // vCellFeatureCount
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);

    countMatStride = options.umiDedup.yes.N + 1;

    UnsignedWord ccgN = 0; // total number of entries in the countCellGeneUMI
    for (auto& cbf : vCellFeatureCount) {
        ccgN += cbf.size();
    };
    countCellGeneUMI.resize(ccgN * countMatStride);
    countCellGeneUMIindex.resize(nCB + 1);
    countCellGeneUMIindex[0] = 0;

    for (std::uint32_t icb = 0; icb < nCB; icb++) { // copy vCellFeatureCount

        nGenePerCB[icb] = vCellFeatureCount[icb].size();
        countCellGeneUMIindex[icb + 1] =
            countCellGeneUMIindex[icb] + nGenePerCB[icb] * countMatStride;

        for (UnsignedWord ic = countCellGeneUMIindex[icb], ig = 0;
             ic < countCellGeneUMIindex[icb + 1]; ic += countMatStride, ++ig) { // loop over genes

            countCellGeneUMI[ic + 0] = vCellFeatureCount[icb][ig].feature;

            if (options.umiDedup.yes.NoDedup)
                countCellGeneUMI[ic + options.umiDedup.countInd.NoDedup] =
                    vCellFeatureCount[icb][ig].count[0];
            if (options.umiDedup.yes.Exact)
                countCellGeneUMI[ic + options.umiDedup.countInd.Exact] =
                    vCellFeatureCount[icb][ig].count[1]; // collapsed UMI count

            nUMIperCB[icb] += countCellGeneUMI[ic + options.umiDedup.countInd.main];
        };
    };

    // sum stats
    nReadPerCBtotal = nReadPerCB;
    nReadPerCBunique = nReadPerCB;

    std::uint32_t nReadPerCBmax = 0;
    for (std::uint32_t icb = 0; icb < nCB; icb++) {

        nReadPerCBmax = std::max(nReadPerCBmax, nReadPerCB[icb]);
        readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += nReadPerCBtotal[icb];
        readFeatSum->stats.V[readFeatSum->stats.yessubWLmatch_UniqueFeature] +=
            nReadPerCBunique[icb];

        readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[icb];
        if (nGenePerCB[icb] > 0)
            ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];
    };

    readFeatSum->stats.V[readFeatSum->stats.yessubWLmatchExact] =
        readFeatSum->stats.V[readFeatSum->stats.yesWLmatch];

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished SmartSeq counting"
                            << std::endl;
};

} // namespace nustar::star2711b::solo::detail
