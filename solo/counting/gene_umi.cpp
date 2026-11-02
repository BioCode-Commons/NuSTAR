#include "../barcode/sequence.hpp"
#include "feature_counter.hpp"
#include "../support/memory.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"

namespace nustar::star2711b::solo::detail {

void FeatureCounter::countCBgeneUMI() {
    time_t rawTime;

    rguStride = 2;
    if (options.readIndexYes[featureType])
        rguStride = 3; // to keep readI column

    if (options.readInfoYes[featureType]) {
        readInfo.resize(nReadsInput, {(UnsignedWord)-1, (std::uint32_t)-1});
        time(&rawTime);
        context_.inOut->logMain << timeMonthDayTime(rawTime)
                                << " ... Allocated and initialized readInfo array, nReadsInput = "
                                << nReadsInput << std::endl;
    };

    rGeneUMI =
        new std::uint32_t[rguStride *
                          nReadsMapped]; // big array for all CBs - each element is gene and UMI
    rCBp = new std::uint32_t*[nCB + 1];
    std::uint32_t** rCBpa = new std::uint32_t*[options.cbWLsize + 1];

    rCBp[0] = rGeneUMI;
    rCBpa[0] = rGeneUMI;
    nCB = 0; // will count it again below
    for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii] > 0) {
            rCBp[nCB + 1] = rCBp[nCB] + rguStride * readFeatSum->cbReadCount[ii];
            ++nCB;
        };
        rCBpa[ii + 1] = rCBp[nCB];
    };

    // read and store the CB/gene/UMI from files
    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Finished allocating arrays for Solo "
                            << nReadsMapped * rguStride * 4.0 / 1024 / 1024 / 1024 << " GiB"
                            << std::endl;
    readFlagCounts.flagCounts.reserve(nCB * 3 / 2);
    readFlagCounts.flagCountsNoCB = {};
    std::vector<std::uint32_t> nReadPerCBunique1(options.cbWLsize),
        nReadPerCBmulti1(
            options.cbWLsize); // temp arrays to record # of reads for all cells in the WL
    for (int ii = 0; ii < context_.runThreadN; ii++) { // TODO: this can be parallelized
        readFeatAll[ii]->inputRecords(rCBpa, rguStride, barcode_sum_->cbReadCountExact, readInfo,
                                      readFlagCounts, nReadPerCBunique1, nReadPerCBmulti1);
        readFeatSum->addStats(*readFeatAll[ii]); // sum stats: has to be done after inputRecords,
                                                 // since the stats values are updated there
    };
    readFlagCounts.countsAddNoCBarray(
        readFeatSum->readFlag
            .flagCountsNoCB); // add no-CB counts calculated in SoloReadFeature_record.cpp and not
                              // recorded to temp Solo files

    nReadPerCBtotal.resize(nCB);
    nReadPerCBunique.resize(nCB);
    for (std::uint32_t icb = 0; icb < nCB; icb++) {
        nReadPerCBunique[icb] = nReadPerCBunique1[indCB[icb]];
        nReadPerCBtotal[icb] = nReadPerCBunique[icb] + nReadPerCBmulti1[indCB[icb]];
    };

    // debug
    /*{
        uint64 n1=0,n2=0;
        for (uint32 icb=0; icb<nCB; icb++) {
            n1 += nReadPerCBtotal[icb];
        };

        for (auto &c: readFlagCounts.flagCounts) {
            n2+= c.second[readFlagCounts.countedU] + c.second[readFlagCounts.countedM];
        };

        cout << "n1,2=" << n1<<" "<<n2<<endl;
    };*/

    nReadPerCB.resize(nCB);
    nReadPerCBmax = 0;
    for (std::uint32_t iCB = 0; iCB < nCB; iCB++) {
        nReadPerCB[iCB] = (rCBpa[indCB[iCB]] - rCBp[iCB]) /
                          rguStride; // number of reads that were matched to WL, rCBpa accumulated
                                     // reference to the last element+1 for multimappers this is the
                                     // number of all alignments > number of reads
        nReadPerCBmax = std::max(nReadPerCBmax, nReadPerCB[iCB]);
        // readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += nReadPerCB[iCB];
    };

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Finished reading reads from Solo files nCB=" << nCB
                            << ", nReadPerCBmax=" << nReadPerCBmax;
    context_.inOut->logMain << ", yesWLmatch="
                            << readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] << std::endl;
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);

    // dedup options        //gene ID
    countMatStride = options.umiDedup.yes.N + 1;
    countCellGeneUMI.resize(nReadsMapped * countMatStride / 5 +
                            16); // 5 is heuristic, will be resized if needed
    countCellGeneUMIindex.resize(nCB + 1, 0);

    if (options.multiMap.yes.multi) {
        // gene   //uniform  //rescue
        countMatMult.s = 1 + options.multiMap.yes.N * options.umiDedup.yes.N;
        countMatMult.m.resize(nReadsMapped * countMatMult.s / 5 + 16);
        countMatMult.i.resize(nCB + 1, 0);
    };

    deduplicateAllCells();

    context_.inOut->logMain << "RAM for solo feature " << FeatureTypes::Names[featureType] << "\n"
                            << linuxProcMemory() << std::flush;
    delete[] rGeneUMI;
    delete[] rCBp;
    delete[] rCBpa;

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs"
                            << std::endl;
};

} // namespace nustar::star2711b::solo::detail
